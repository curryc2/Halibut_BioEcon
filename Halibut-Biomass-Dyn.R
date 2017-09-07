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

#=============================================================
##### CONTROL SECTION #####
do.init.plots <- FALSE

#=============================================================
#1) Adjust halibut object values based on inputs from spreadsheet
#2) Calculate survival/growth by age
#3) Calculate selectivity by age
halibut <- read_update_params()

#=========================================
if(do.init.plots==TRUE) {
  #PRELIMINARY PLOTS
  plot.landings(dpi=600)
  plot.growth_allometry(ageSchedules=halibut, dpi=500)
  plot.selectivity(selectivity=halibut, dpi=500, pt.blk=FALSE)
}

#=========================================
#Extract variables


#INPUT FISHING MORTALITY RATES
fmort <- read.xlsx('Halibut Model Inputs.xlsx', sheetName='Fmort')[,-1]

#Input control parameters
in.control <- read.xlsx('Halibut Model Inputs.xlsx', sheetName='Control')
n.year <- in.control$Value[in.control$Par=='n.yrs'] #Number of years to simulate
years <- 1:n.year
Bstart <- in.control$Value[in.control$Par=='Bstart'] #Starting Biomass


n.age  <- halibut$theta$A
n.gear <- dim(halibut$MP)[1]
n.sex  <- halibut$theta$H
va <- as.array(halibut$selex) #Overall selectivity

gears <- as.vector(halibut$MP$sector)


probCap <- as.array(halibut$probCap) #Probability of capture @ age
probRetain <- as.array(halibut$probRetain)

wa <- halibut$ageSc$wa #Weight @ age
fa <- halibut$ageSc$fa #Female spawning biomass @ age

#Age Schedule stuff
mx <- halibut$ageSc$mx  #Natural mortality @ age
la <- halibut$ageSc$la  #Length @ age
wa <- halibut$ageSc$wa  #Weight @ age
ma <- halibut$ageSc$ma  #Maturity @ age
fa <- halibut$ageSc$fa  #Fecundity @ age
lx <- halibut$ageSc$lx  #Survivorship to age

#Age Information
ages <- halibut$theta$age
plus.age <- halibut$theta$A
sexes <- c('Female','Male')

#Recruitment
steep <- halibut$rec$steep
sigma_rec <- halibut$rec$sigma_rec
ro <- halibut$rec$ro
bo <- halibut$theta$bo*1e6

#=========================================

lz  <- matrix(1/n.sex,nrow=n.sex,ncol=n.age)
za  <- matrix(0,nrow=n.sex,ncol=n.age)
qa  <- array(0,dim=c(n.sex,n.age,n.gear))
pa  <- array(0,dim=c(n.sex,n.age,n.gear))
ra  <- array(0,dim=c(n.sex,n.age,n.gear))
dlz <- array(0,dim=c(n.sex,n.age,n.gear))

#========================================================
#Define Data Structures
B <- array(dim=c(n.sex, n.year, n.age), dimnames=list(sexes, years, ages)) #Biomass (pounds)
N <- array(dim=c(n.sex, n.year, n.age), dimnames=list(sexes, years, ages)) #Numbers
C.b <- array(dim=c(n.sex, n.year, n.age), dimnames=list(sexes, years, ages)) #Catch (lbs)
C.n <- array(dim=c(n.sex, n.year, n.age), dimnames=list(sexes, years, ages)) #Catch (number)
harvest.b <- array(dim=c(n.sex, n.year, n.age, n.gear), dimnames=list(sexes, years, ages, gears))  #Harvest (lbs) by gear type
harvest.n <- array(dim=c(n.sex, n.year, n.age, n.gear), dimnames=list(sexes, years, ages, gears))  #Harvest (number) by gear type


#Total Instantaneous mortality
Z.a <- array(dim=c(n.sex, n.year, n.age), dimnames=list(sexes, years, ages)) 
F.a <- array(dim=c(n.sex, n.year, n.age), dimnames=list(sexes, years, ages)) #Fishing mortality

#Continuous
surv <- array(dim=c(n.sex, n.year, n.age), dimnames=list(sexes, years, ages))
mort <- array(dim=c(n.sex, n.year, n.age), dimnames=list(sexes, years, ages))

#Recruitment
ssb <- array(dim=c(n.sex, n.age, n.year), dimnames=list(sexes, ages, years)) #Female spawning-stock biomass
rec <- array(dim=c(n.sex, n.year), dimnames=list(sexes, years))

#Define initial population structure based on equilibirum conditions

#Calculate expected biomass proportions by weight
# init.prop.B <- (lx*wa)^-1
# #Standardize
# init.prop.B <- init.prop.B/sum(init.prop.B)
# 
# bo/wa[1,1]

#This should be updated
init.prop <- matrix(nrow=n.sex, ncol=n.age, dimnames=list(sexes,ages))
for(a in 1:n.age) {
  if(a==1) {
    init.prop[,a] <- bo#*mx[,a]
  }else {
    init.prop[,a] <- bo*exp(-(a-1)*mx[,a])
  }
  if(a==plus.age) {
    init.prop[,a] <- init.prop[,a]/(1-exp(-mx[,a]))
  }
}
#Multiply by spawning sbpr @ age
# init.prop <- init.prop
init.prop <- init.prop * wa
# for(h in 1:n.sex) {
#   init.prop[h,] <- init.prop[h,]/sum(init.prop[h,])
# }

init.prop <- init.prop/sum(init.prop)

# lst <- melt(init.prop)
# names(lst) <- c('Sex','Age','Prop')
# ggplot(lst, aes(x=Age, y=Prop, color=Sex))+ geom_line() + geom_point(pch=21, fill='black') 

#Initial Biomass
B[,1,] <- Bstart*1e6 * (init.prop)

#Initial Numbers of Individuals
N[,1,] <- B[,1,] / wa


#Calculate equilibrium recruitment
# sp <- seq(from=1, to=5e8, length.out=1e3)
# bh <- beverton_holt_recruit(sp, steep, 1e7)
# # rk <- ricker_recruit(sp, steep, 1e7)
# 
# plot(bh~sp, type='l', col='blue')
# lines(rk~sp, col='red')

##################################################
#BEGIN SIMULATION

#Calculate age-specific fishing mortality rates
dim(F.a)
dim(va)

y <- 2
for(y in 2:n.year) {
  
  #Initial Recruitment
  ssb[,,y-1] <- fa*N[,y-1,]
  
  #Ricker
  # rec[,y-1] <- 0.5 * ricker_recruit(ssb[y-1], steep, bo)
  #Beverton-Holt
  rec[,y-1] <-  0.5 * beverton_holt_recruit(sum(ssb[,,y-1]), steep, ro)#bo) #* exp(rnorm(1,0,sigma_rec) - ((sigma_rec^2)/2))
  
  for(a in 1:n.age) {
    #Update Numbers and Biomass Matrix
    if(a==1) { #Age-1
      N[,y,a] <- rec[,y-1]
      B[,y,a] <- rec[,y-1]*wa[,a]
    }else {
      h <- 1
      for(h in 1:n.sex) {
        #Instantaneous Version
        F.a[h,y-1,a-1] <- sum(fmort[y-1,]*va[h,a-1,])
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
          temp.F <- fmort[y-1,g]*va[h,a-1,g]
          # temp.Z <- temp.F + mx[h,a-1]
          temp.Z <- sum(fmort[y-1,]*va[h,a-1,]) + mx[h,a-1]
          
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
        F.a[h,y-1,a] <- sum(fmort[y-1,]*va[h,a,])
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
          temp.F <- fmort[y-1,g]*va[h,a,g]
          # temp.Z <- temp.F + mx[h,a]
          temp.Z <- sum(fmort[y-1,]*va[h,a,]) + mx[h,a]
          # 
          harvest.n[h,y-1,a,g] <- N[h,y-1,a] * (temp.F/temp.Z) * (1-exp(-1*temp.Z))
          harvest.b[h,y-1,a,g] <- harvest.n[h,y-1,a,g] * wa[h,a]
        }#next gear
      }#next sex
    }# If plus age group
  }#next age  
}#next y


#=============================
#Some Exploratory Plotting
dim(C.n)
dim(C.b)

#Total Biomass plot
#Unfished Biomass is coastwide 709 million lbs
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

#


# g <- list.ssb[list.ssb$Sex=='Female',] %>% group_by(Sex, Year) %>% summarize(total=sum(value, na.rm=TRUE)) %>% ggplot() + aes(x=Year, y=total, fill=Sex, group=Sex) + geom_line()

# list.ssb[list.ssb$Sex=='Female',] %>% group_by(Sex, Year) %>% 
#     summarize(total=sum(value, na.rm=TRUE), prop=(value/sum(value)))

#Plotting abundance
dim(N)

list.N <- melt(N)
names(list.N) <- c('Sex', 'Year', 'Age', 'value')
# list.N$Age <- as.factor(list.N$Age)

g2 <- ggplot(list.N, aes(x=Year, y=value/1e6, fill=Age, group=Age)) +
      theme_gray() +
      geom_area() + 
      facet_wrap(~Sex, ncol=1) +
      scale_fill_gradient2(midpoint=plus.age/2, low='darkblue', mid='green', high='red') +
      labs(y='Abundance (millions)')   
      # scale_fill_brewer(palette="RdYlGn")

g2


#Plotting abundance proportions
list.N.prop <- data.frame(list.N %>% group_by(Sex, Year) %>% mutate(prop=value/sum(value)))

g3 <- ggplot(list.N.prop, aes(x=Year, y=prop, fill=Age, group=Age)) +
        theme_gray() +
        geom_area(alpha=0.75) +
        facet_wrap(~Sex, ncol=1) +
        scale_fill_gradient2(midpoint=plus.age/2, low='darkblue', mid='green', high='red') +
        labs(y='Proportion of Total Abundance (N)') 
g3

list.b <- melt(B)
names(list.b) <- c('Sex','Year','Age','value')

list.B.prop <- data.frame(list.b %>% group_by(Sex, Year) %>% mutate(prop=value/sum(value)))

g4 <- ggplot(list.B.prop, aes(x=Year, y=prop, fill=Age, group=Age)) +
  theme_gray() +
  geom_area(alpha=0.75) +
  facet_wrap(~Sex, ncol=1) +
  scale_fill_gradient2(midpoint=plus.age/2, low='darkblue', mid='green', high='red') +
  labs(y='Proportion of Total Biomass (B)') 
g4

g5 <- ggplot(list.B.prop, aes(x=Year, y=value/1e6, fill=Age, group=Age)) +
  theme_gray() +
  geom_area(alpha=0.75) +
  facet_wrap(~Sex, ncol=1) +
  scale_fill_gradient2(midpoint=plus.age/2, low='darkblue', mid='green', high='red') +
  labs(y='Total Biomass (millions of lbs)') 
g5 




# g.rec <- ggplot()

#2016
# Commercial fishery landings= 25 milion lbs
#Catch
C.b.list <- melt(C.b)
names(C.b.list) <- c('Sex','Year','Age','value')

g.6 <- ggplot(C.b.list, aes(x=Year, y=value/1e6, fill=Age, group=Age)) +
  theme_gray() +
  geom_area(alpha=0.75) +
  facet_wrap(~Sex, ncol=1) +
  scale_fill_gradient2(midpoint=plus.age/2, low='darkblue', mid='green', high='red') +
  labs(y='Total Catch (millions of lbs)') 
g.6

#Harvest by Sector
harv.b.list <- melt(harvest.b)
names(harv.b.list) <- c('Sex','Year','Age','Sector','value')

g.7 <- ggplot(harv.b.list, aes(x=Year, y=value/1e6, fill=Age, group=Age)) +
  theme_gray() +
  geom_area(alpha=0.75) +
  facet_grid(Sector~Sex) +
  scale_fill_gradient2(midpoint=plus.age/2, low='darkblue', mid='green', high='red') +
  labs(y='Total Catch (millions of lbs)') 

g.7


