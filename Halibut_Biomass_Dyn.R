#************************************************************************************************
#Project Name: HALIBUT PSC BIOECONOMIC MODEL  
#Creator: Curry J. Cunningham, NOAA/NMFS
#Date: 12.19.16
#
#PURPOSE: Simulate Pacific Halibut population dynamics with an age and sex-structured model,
#           with three types of fishing mortaltity
#
#  MODEL SPECIFICATION:
#    Length-at-age: Von-bert
#    Weight-at-age: Exponential fxn. of length
#    Natural Mortality: Sex-specific, but NOT varying across ages or time.
#    Fishing Mortality: Input values by year and gear type
#    Selectivity: Joint probability of p(capture)*p(retain)
#       p(capture) = 
#*************************************************************************************************
#NOTES:
#  1) Total fishing mortality includes mortality for discards.
#  2) Dimensions for demographic parameters halibut$theta... are sex first, Female,Male
#  3) Successfully setup with GitHub for cross-platform use
#  4) Population viability analysis?
#       Only comfortable 
#*************************************************************************************************
require(ggplot2)
require(gridExtra)
require(cowplot)
require(dplyr)
require(snowfall)
require(parallel)
require(reshape2)
require(xlsx)
require(corrplot)

#Working Directory
if(.Platform$OS.type=='unix') { #Laptop
  wd <- '/Users/curryc2/Documents/2016/Reimer Halibut/Halibut_BioEcon'
}else { #Desktop
  wd <- "//nmfs.local/AKC-ABL/Users/curry.cunningham/My Documents/Projects/Halibut MSE/Halibut_BioEcon"  
}
setwd(wd)

#Sources
source('R/Halibut_Plot_Fxns.R')


#PRELIMINARY PLOTS
# plot.landings(dpi=600)



#READ IN DATA from excel input file
load('data/halibut.rda')
source('R/fisheryFootprint_plus.R') #This is an updated version of Steve's functions

#=============================================================
##### CONTROL SECTION #####
do.init.plots <- FALSE


#=============================================================
#Adjust halibut object values based on inputs from spreadsheet
#Growth
in.growth <- read.xlsx('Halibut Model Inputs.xlsx', sheetName='Growth')
halibut$theta$linf <- as.numeric(in.growth[in.growth$Par=='linf',(2:3)])
halibut$theta$vbk <- as.numeric(in.growth[in.growth$Par=='vbk',(2:3)])
halibut$theta$to <- as.numeric(in.growth[in.growth$Par=='to',(2:3)])
halibut$theta$a <- as.numeric(in.growth[in.growth$Par=='a',(2:3)])
halibut$theta$b <- as.numeric(in.growth[in.growth$Par=='b',(2:3)])
#Maturity
in.maturity <- read.xlsx('Halibut Model Inputs.xlsx', sheetName='Maturity')
halibut$theta$ahat <- as.numeric(in.maturity[in.maturity$Par=='ahat',(2:3)])
halibut$theta$ghat <- as.numeric(in.maturity[in.maturity$Par=='ghat',(2:3)])
#Mortality
in.mortality <- read.xlsx('Halibut Model Inputs.xlsx', sheetName='Mortality')
halibut$theta$m <- as.numeric(in.mortality[in.mortality$Par=='m',(2:3)])
halibut$theta$A <- max(in.mortality[in.mortality$Par=='A',2:3])
#Fishery Selectivity
in.selex <- read.xlsx2('Halibut Model Inputs.xlsx', sheetName='FisherySelectivity', stringsAsFactors=FALSE, header=TRUE)
halibut$MP$slx1 <- as.numeric(in.selex[in.selex$Par=='mu',c(2:5)])
halibut$MP$slx2 <- as.numeric(in.selex[in.selex$Par=='sigma',c(2:5)])
halibut$MP$slx3 <- as.numeric(in.selex[in.selex$Par=='gamma',c(2:5)])
halibut$MP$slx4 <- as.numeric(in.selex[in.selex$Par=='plus.age',c(2:5)])  #Plus group age for selectivity
halibut$MP$slim <- as.numeric(in.selex[in.selex$Par=='slim',c(2:5)])  #Minimum size for retention
halibut$MP$ulim <- as.numeric(in.selex[in.selex$Par=='ulim',c(2:5)])  #Maximum size for retention
halibut$MP$dmr <- as.numeric(in.selex[in.selex$Par=='dmr',c(2:5)])  #Discard Mortality Rate
# halibut$MP$pscLimit <- as.numeric(in.selex[in.selex$Par=='pscLimit',c(2:5)])  #Prohibited species catch limit
halibut$MP$pYPR <- as.numeric(in.selex[in.selex$Par=='pYPR',c(2:5)])
halibut$MP$pYPR <- as.numeric(in.selex[in.selex$Par=='pYPR',c(2:5)])
halibut$MP$pYPR <- as.numeric(in.selex[in.selex$Par=='pYPR',c(2:5)])


#INPUT FISHING MORTALITY RATES
fmort <- read.xlsx('Halibut Model Inputs.xlsx', sheetName='Fmort')

#Input control parameters
in.control <- read.xlsx('Halibut Model Inputs.xlsx', sheetName='Control')
n.year <- in.control$Value[in.control$Par=='n.yrs'] #Number of years to simulate
Bstart <- in.control$Value[in.control$Par=='Bstart'] #Starting Biomass

#=========================================
halibut <- getSelectivities(halibut)
#Determine Age Schedules
# ageSchedules <- getAgeSchedules(halibut)
#Plot Age Schedule
if(do.init.plots==TRUE) { plot.growth_allometry(ageSchedules=halibut, dpi=500) }
#

#=========================================
#Determine Selectivities
# selectivity <- getSelectivities(halibut)

if(do.init.plots==TRUE) { plot.selectivity(selectivity=halibut, dpi=500, pt.blk=FALSE) }

#=========================================
#Extract variables
n.age  <- halibut$theta$A
n.gear <- dim(halibut$MP)[1]
n.sex  <- halibut$theta$H
va    <- as.array(halibut$selex) #Overall selectivity
wa    <- halibut$ageSc$wa #Weight at Age
fa    <- halibut$ageSc$fa #Fecundity at age

#Age Schedule stuff
mx <- halibut$ageSc$mx

sexes <- c('Female','Male')


lz  <- matrix(1/n.sex,nrow=n.sex,ncol=n.age)
za  <- matrix(0,nrow=n.sex,ncol=n.age)
qa  <- array(0,dim=c(n.sex,n.age,n.gear))
pa  <- array(0,dim=c(n.sex,n.age,n.gear))
ra  <- array(0,dim=c(n.sex,n.age,n.gear))
dlz <- array(0,dim=c(n.sex,n.age,n.gear))

#=========================================
#Define Data Structures
N <- array(dim=c(n.sex, n.year, n.age), dimnames=list(sexes, c(1:n.year), halibut$theta$age))
surv <- array(dim=c(n.sex, n.year, n.age), dimnames=list(sexes, c(1:n.year), halibut$theta$age))
rec <- array(dim=c(n.sex, n.year), dimnames=list(sexes, c(1:n.year)))

#Define initial population structure based on equilibirum conditions
Bstart*1e6





































