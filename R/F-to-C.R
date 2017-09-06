
#' Convert fishing mortality rate F to catch C, by sector
#'
#' @param F.input Vector of fishing mortality rates for each of the four fishing sectors
#' @param N.input Matrix of abundances by sex and age class
#' @param Catch.biom Boolean describing whether catch will be returned as biomass or numbers
#' @param halibut List object containing all of the input Pacific Halibut life history and fishery selectivity information 
#'
#' @return A vector of catches in biomass or numbers, for each of four fishing sectors
#' @export
#'
F_to_C <- function(F.input, N.input,  Catch.biom=TRUE, halibut=halibut) {
  # ### TESTING ###
  # F.input <- c(0.1,0.05,0.01,0.01)
  # N.input <- N[,5,]
  ##############
  
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
  if(Catch.biom==TRUE) {
    output <- pred.total.b
  }else {
    output <- pred.total.n
  }
  return(output)
}

# F_to_C(F.input=c(0.2,0.1,0.05,0), N.input=N[,5,], Catch.biom=TRUE, halibut=halibut)
