#' Function to create data objects for simulation
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
  out$surv <- surv
  out$mort <- mort
  out$ssb <- ssb
  out$rec <- rec
  return(out)
}