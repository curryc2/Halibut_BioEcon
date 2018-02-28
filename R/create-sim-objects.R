#' Function to create data objects for simulation
#'
#' @param halibut List object containing demographic and fishery parameters for simulation.
#'
#' @return A named list with all of the data structures to be attached after function call
#' @export
#'
create_sim_objects <- function() {
  lz  <- matrix(1/n.sex,nrow=n.sex,ncol=n.age)
  za  <- matrix(0,nrow=n.sex,ncol=n.age)
  qa  <- array(0,dim=c(n.sex,n.age,n.gear))
  pa  <- array(0,dim=c(n.sex,n.age,n.gear))
  ra  <- array(0,dim=c(n.sex,n.age,n.gear))
  dlz <- array(0,dim=c(n.sex,n.age,n.gear))
  
  #========================================================
  #Define Data Structures
  
  sims <- paste0('sim',c(1:n.sims))
  
  B <- array(dim=c(n.sex, n.year, n.age, n.sims), dimnames=list(sexes, years, ages, sims)) #Biomass (pounds)
  N <- array(dim=c(n.sex, n.year, n.age, n.sims), dimnames=list(sexes, years, ages, sims)) #Numbers
  C.b <- array(dim=c(n.sex, n.year, n.age, n.sims), dimnames=list(sexes, years, ages, sims)) #Catch (lbs)
  C.n <- array(dim=c(n.sex, n.year, n.age, n.sims), dimnames=list(sexes, years, ages, sims)) #Catch (number)
  harvest.b <- array(dim=c(n.sex, n.year, n.age, n.gear, n.sims), dimnames=list(sexes, years, ages, gears, sims))  #Harvest (lbs) by gear type
  harvest.n <- array(dim=c(n.sex, n.year, n.age, n.gear, n.sims), dimnames=list(sexes, years, ages, gears, sims))  #Harvest (number) by gear type
  
  
  #Total Instantaneous mortality
  Z.a <- array(dim=c(n.sex, n.year, n.age, n.sims), dimnames=list(sexes, years, ages, sims)) 
  F.a <- array(dim=c(n.sex, n.year, n.age, n.sims), dimnames=list(sexes, years, ages, sims)) #Age-specific Fishing mortality
  Fmort <- array(dim=c(n.sex, n.year, n.sims), dimnames=list(sexes, years, sims)) #Annual Fishing mortality
  
  #Continuous
  surv <- array(dim=c(n.sex, n.year, n.age, n.sims), dimnames=list(sexes, years, ages, sims))
  mort <- array(dim=c(n.sex, n.year, n.age, n.sims), dimnames=list(sexes, years, ages, sims))
  
  #Recruitment
  ssb <- array(dim=c(n.sex, n.age, n.year, n.sims), dimnames=list(sexes, ages, years, sims)) #Female spawning-stock biomass
  rec <- array(dim=c(n.sex, n.year, n.sims), dimnames=list(sexes, years, sims))
  
  #Return section
  out <- NULL
  out$lz <- lz
  out$za <- za
  out$qa <- qa
  out$pa <- pa
  out$ra <- ra
  out$dlz <- dlz
  out$B <- B
  out$N <- N
  out$C.b <- C.b
  out$C.n <- C.n
  out$harvest.b <- harvest.b
  out$harvest.n <- harvest.n
  out$Z.a <- Z.a
  out$F.a <- F.a
  out$Fmort <- Fmort
  out$surv <- surv
  out$mort <- mort
  out$ssb <- ssb
  out$rec <- rec
  return(out)
  
}