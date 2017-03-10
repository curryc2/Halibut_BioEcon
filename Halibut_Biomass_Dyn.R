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
#*************************************************************************************************
require(ggplot2)
require(gridExtra)
require(cowplot)
require(dplyr)
require(snowfall)
require(parallel)
require(reshape2)
require(xlsx)

#Working Directory
if(.Platform$OS.type=='unix') { #Laptop
  wd <- '/Users/curryc2/Documents/2016/Reimer Halibut/Halibut_BioEcon'
}else { #Desktop
  wd <- "//nmfs.local/AKC-ABL/Users/curry.cunningham/Desktop/GitHub/Halibut_BioEcon"  
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
halibut$theta$m <- as.numeric(in.maturity[in.maturity$Par=='m',(2:3)])
halibut$theta



out <- getSelectivities(halibut)







