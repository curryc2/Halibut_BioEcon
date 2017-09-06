

#' Convert catch C to fishing mortality rate F, by sector, given specified information on selectivity and natural mortality.
#'   F values are found through simple minimization.
#'
#' @param C.input Matrix of catches for each of the four fishing sectors
#' @param N.input Matrix of abundances by sex and age class
#' @param Catch.biom Boolean describing whether C.input is in units of biomass or numbers
#' @param halibut List object containing all of the input Pacific Halibut life history and fishery selectivity information 
#'
#' @return Named vector of estimated fishing mortality rates by fishing sector
#' 
#' @export
#'
C_to_F <- function(C.input, N.input, Catch.biom=TRUE, ...) {
  
  ### TESTING ###
  # C.input <- c(1e7,1e7,1e6,1e6)
  # N.input <- N[,5,]
  # Catch.biom <- TRUE
  
  ##############
  
  require(bbmle)
  # Parameters to be estimated #
  # F.input <- rep(0.2,4)#vector(length=n.gear)
  
  
  # #Extract Halibut Parameters
  # wa <- halibut$ageSc$wa
  # mx <- halibut$ageSc$mx
  # 
  # n.age  <- halibut$theta$A
  # n.gear <- dim(halibut$MP)[1]
  # n.sex  <- halibut$theta$H
  # va <- as.array(halibut$selex) #Overall selectivity
  
  
  
  diff_fun <- function(F.IFQ, F.PSC, F.SPT, F.PER, C.input, Catch.biom=Catch.biom, halibut=halibut) {
    #Get Parameters
    F.input <- c(F.IFQ, F.PSC, F.SPT, F.PER)
    
    # with(halibut, {
    pred.harvest.n <- array(dim=c(n.sex,n.age,n.gear))
    pred.harvest.b <- array(dim=c(n.sex,n.age,n.gear))
      
    h <- 1
    for(h in 1:n.sex) {
      g <- 1
      for(g in 1:n.gear) {
        temp.F <- F.input[g]*va[h,,g]
        temp.Z <- temp.F + mx[h,]
        pred.harvest.n[h,,g] <- N.input[h,] * (temp.F/temp.Z) * (1-exp(-1*temp.Z))
        pred.harvest.b[h,,g] <- pred.harvest.n[h,,g] * wa[h,]
      }#next g
    }#next h
    
    #Convert to prediction for Minimization
    pred.total.b <- apply(pred.harvest.b, 3, sum)
    pred.total.n <- apply(pred.harvest.n, 3, sum)
    
    #LIKELIHOOD COMPONENT
    if(Catch.biom==TRUE) {
      SSQ <- sum((C.input-pred.total.b)^2)
    }else {
      SSQ <- sum((C.input-pred.total.n)^2)
    }
    # })
    return(SSQ)
  }
  
  #Minimization
  # diff_fun(c(0.2,0.01,0.01,0.01), C.input, Catch.biom)
  
  # diff_fun(F.IFQ=0.2, F.PSC=0.2, F.SPT=0.2, F.PER=0.2, C.input, Catch.biom)
  
  fit <- mle2(diff_fun, start=list(F.IFQ=0.2, F.PSC=0.2, F.SPT=0.2, F.PER=0.2),
              data=list(C.input=C.input, Catch.biom=Catch.biom, halibut=halibut),
              method='Nelder-Mead', 
              control=list(trace=FALSE, maxit=1e4))

  #Shoot Again
  fit <- mle2(diff_fun, start=list(F.IFQ=coef(fit)[1], F.PSC=coef(fit)[2], F.SPT=coef(fit)[3], F.PER=coef(fit)[4]),
              data=list(C.input=C.input, Catch.biom=Catch.biom, halibut=halibut),
              method='Nelder-Mead', 
              control=list(trace=FALSE, maxit=1e4))
  
  return(coef(fit))
  
} 

# C_to_F(C.input=c(1e7,1e7,1e6,1e6), N.input=N[,5,], Catch.biom=TRUE)
