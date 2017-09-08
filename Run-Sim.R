#************************************************************************************************
#Project Name: HALIBUT PSC BIOECONOMIC MODEL - Run A Simple Simulation Example
#Creator: Curry J. Cunningham, NOAA/NMFS
#Date: 8.19.17
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
#  5) From 2016 RARA: Age-0 recruitment Coastwide ~55 million, Areas as fleets ~ 25 million
#  6) From 2013 Presentation SSB0= 
#     a) Coastwide = 709 million lbs.
#     b) Areas 2B/2C/3A (SE/BC) = 489 million lbs.
#  7) Finding ro:
#     a) Set Fmort to zero
#     b) Find ro that results in equilibrium SSB ~ 709 million
#     c) This happens to be ro = 2.5e6
#====================================================
# If Recruitment in lbs:
#  8a) In unfished condition w/ ro = 2.5e6
#     Starting biomass at 1500 (million lbs) leads to stability.
#
# If Recruitment in fish:
#  8b) In unfished condition w/ ro = 3e7
#     Starting biomass at 1600 (million lbs) leads to stability.
#
# Fishing Mortality rates that bring population into 
#====================================================
#Ricker
# R <- S*exp(1.25*log(5*h)*(1-(S/Seq)))
#Beverton-Holt
# R <- S/(1-((5*h-1)/(4*h))*(1-(S/(Seq))))

#Fmort Starts: 0.075, 0.05, 0.01, 0.01

#*************************************************************************************************
require(ggplot2)
require(gridExtra)
require(cowplot)
require(tidyverse)
require(snowfall)
require(parallel)
require(reshape2)
require(xlsx)
require(corrplot)
require(R.utils)
require(ggthemes)

#Working Directory

#Sources
source('R/Halibut-Plot-Fxns.R')
source('R/ricker-recruit.R')
source('R/beverton-holt-recruit.R')
source('R/get-fished-survivorship.R')

source('R/fisheryFootprint-plus.R') #This is an updated version of Steve's functions
source('R/read-update-params.R')
source('R/extract-params.R')
source('R/create-sim-objects.R')
source('R/HCR/HCR-linear.R')

source('R/calc-init-age-prop.R')

#=============================================================
##### CONTROL SECTION #####
in.control <- read.xlsx('Halibut Model Inputs.xlsx', sheetName='Control')
n.year <- in.control$Value[in.control$Par=='n.yrs'] #Number of years to simulate
years <- 1:n.year
Bstart <- 700#in.control$Value[in.control$Par=='Bstart'] #Starting Biomass

SSB0 <- 709e6 #709 million lbs from 2013 RARA

#For Harvest Control Rules
floors <- rep(0,4)
ceilings <- rep(0.4,4)#c(0.1,0.05,0.01,0.01)
ascent.range <- c(0.2,0.4)
#=============================================================
#1) Adjust halibut object values based on inputs from spreadsheet
#2) Calculate survival/growth by age
#3) Calculate selectivity by age
halibut <- read_update_params()

#=============================================================
#Extract variables
params <- extract_params(halibut)
attach(params)

#=============================================================
#Create data objects for simulation
objs <- create_sim_objects()
attach(objs)

#=============================================================
#Initialize Population (year 1)
init.prop <- calc_init_age_prop()
B[,1,] <- Bstart*1e6 * (init.prop)
N[,1,] <- B[,1,] / wa

#=============================================================
#Forward Simulation

y <- 2
for(y in 2:n.year) {
  
  #Initial Recruitment
  ssb[,,y-1] <- fa*N[,y-1,]
  
  
  #HARVEST CONTROL RULE
  temp.ssb <- sum(ssb[,,y-1])
  #Fishing Sectors
  temp.fmort <- vector(length=n.gear)
  g <- 1
  for(g in 1:n.gear) {
    temp.fmort[g] <- HCR_linear(curr.SSB=temp.ssb, SSB0=SSB0, floor.F=floors[g], ceiling.F=ceilings[g], 
                                  ascent.range=ascent.range, plot=FALSE)
  }#next g
  
  #Ricker
  # rec[,y-1] <- 0.5 * ricker_recruit(ssb[y-1], steep, bo)
  #Beverton-Holt
  rec[,y-1] <-  0.5 * beverton_holt_recruit(sum(ssb[,,y-1]), steep, bo=ro) #* exp(rnorm(1,0,sigma_rec) - ((sigma_rec^2)/2))
  
  for(a in 1:n.age) {
    #Update Numbers and Biomass Matrix
    if(a==1) { #Age-1
      N[,y,a] <- rec[,y-1]
      B[,y,a] <- rec[,y-1]*wa[,a]
      # N[,y,a] <- rec[,y-1]/wa[,a]
      # B[,y,a] <- rec[,y-1]
    }else {
      h <- 1
      for(h in 1:n.sex) {
        #Instantaneous Version
        F.a[h,y-1,a-1] <- sum(temp.fmort*va[h,a-1,])
        Z.a[h,y-1,a-1] <- F.a[h,y-1,a-1] + mx[h,a-1]  #Natural mortality is NOT time-varying
        
        #Continuous
        surv[h,y-1,a-1] <- exp(-Z.a[h,y-1,a-1])
        mort[h,y-1,a-1] <- 1-surv[h,y-1,a-1]
        
        #Update
        
        N[h,y,a] <- N[h,y-1,a-1]*surv[h,y-1,a-1]
        # B[h,y,a] <- B[h,y-1,a-1]*surv[h,y-1,a-1]
        B[h,y,a] <- N[h,y,a]*wa[h,a]
        #Total Catch
        C.n[h,y-1,a-1] <- N[h,y-1,a-1] * (F.a[h,y-1,a-1]/Z.a[h,y-1,a-1]) * (1-exp(-1*Z.a[h,y-1,a-1])) #Catch in number of halibut
        C.b[h,y-1,a-1] <- C.n[h,y-1,a-1] * wa[h,a-1]
        
        g <- 1
        for(g in 1:n.gear) {
          temp.F <- temp.fmort[g]*va[h,a-1,g]
          # temp.Z <- temp.F + mx[h,a-1]
          temp.Z <- sum(temp.fmort*va[h,a-1,]) + mx[h,a-1]
          
          # harvest.n[h,y-1,a-1,g] <- N[h,y-1,a-1] * (F.a[h,y-1,a-1]/Z.a[h,y-1,a-1]) * (1-exp(-1*Z.a[h,y-1,a-1]))
          harvest.n[h,y-1,a-1,g] <- N[h,y-1,a-1] * (temp.F/temp.Z) * (1-exp(-1*temp.Z))
          # harvest.n[h,y-1,a-1,g] <- N[h,y-1,a-1] * (va[h,a-1,g] * (1-exp(-1*(mx[h,a-1]*Z.a[h,y-1,a-1])))) /(Z.a[h,y-1,a-1])
          
          harvest.b[h,y-1,a-1,g] <- harvest.n[h,y-1,a-1,g] * wa[h,a-1]
        }#next gear
      }#next sex
    }
    
    if(a==plus.age) {
      h <- 1
      for(h in 1:n.sex) {
        #Fish in Plus Group
        F.a[h,y-1,a] <- sum(temp.fmort*va[h,a,])
        Z.a[h,y-1,a] <- F.a[h,y-1,a] + mx[h,a]  #Natural mortality is NOT time-varying        
        
        #Continuous
        surv[h,y-1,a] <- exp(-Z.a[h,y-1,a])
        mort[h,y-1,a] <- 1-surv[h,y-1,a]
        
        #Update
        N[h,y,a] <- N[h,y,a] + N[h,y-1,a]*surv[h,y-1,a] #New Entrants (calculated above), plus existing plus group occupants.
        # B[h,y,a] <- B[h,y,a] + B[h,y-1,a]*surv[h,y-1,a] 
        B[h,y,a] <- N[h,y,a] * wa[h,a]
        #Total Catch
        C.n[h,y-1,a] <- N[h,y-1,a] * (F.a[h,y-1,a]/Z.a[h,y-1,a]) * (1-exp(-1*Z.a[h,y-1,a])) #Catch in number of halibut
        C.b[h,y-1,a] <- C.n[h,y-1,a] * wa[h,a]
        
        g <- 1
        for(g in 1:n.gear) {
          temp.F <- temp.fmort[g]*va[h,a,g]
          # temp.Z <- temp.F + mx[h,a]
          temp.Z <- sum(temp.fmort*va[h,a,]) + mx[h,a]
          # 
          harvest.n[h,y-1,a,g] <- N[h,y-1,a] * (temp.F/temp.Z) * (1-exp(-1*temp.Z))
          harvest.b[h,y-1,a,g] <- harvest.n[h,y-1,a,g] * wa[h,a]
        }#next gear
      }#next sex
    }# If plus age group
  }#next age  
}#next y


list.B <- melt(B) 
names(list.B) <- c('Sex','Year','Age','value')

total.B <- list.B %>% group_by(Year, Age) %>% summarize(total=sum(value))

gt <- ggplot(total.B, aes(x=Year, y=total/1e6, fill=Age, group=Age)) + 
  theme_gray() + 
  geom_area(alpha=0.75) +
  scale_fill_gradient2(midpoint=plus.age/2, low='darkblue', mid='green', high='red') +
  labs(y='Total Biomass (millions of lbs)') 
gt

#Plotting ssb
list.ssb <- melt(ssb)
names(list.ssb) <- c('Sex', 'Age', 'Year', 'value')


SSB.eq <- 709
SSB.2017 <- 212

g <- ggplot(list.ssb[list.ssb$Sex=='Female',], aes(x=Year, y=value/1e6, fill=Age, group=Age)) +
  theme_gray() +
  geom_area(alpha=0.75) + 
  scale_fill_gradient2(midpoint=plus.age/2, low='darkblue', mid='green', high='red') +
  labs(y='Spawning Stock Biomass (millions of lbs)') +
  geom_hline(yintercept=SSB.eq, lty=2, col='red', show.legend=TRUE) +
  geom_hline(yintercept=SSB.2017, lty=2, col='black', show.legend=TRUE)
# scale_fill_brewer(palette="RdYlGn")
g


#Detatching section
detach(objs)
detach(params)
# detach(ageSc)
# detach(halibut)



