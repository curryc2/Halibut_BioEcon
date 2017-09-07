calc_SPR <- function() {
  
  ##### NOT YET COMPLETE, PSEUDOCODE ONLY #####
  
  #Create objects
  za  <- matrix(0,nrow=n.sex,ncol=n.age)
  qa  <- array(0,dim=c(n.sex,n.age,n.gear))
  pa  <- array(0,dim=c(n.sex,n.age,n.gear))
  ra  <- array(0,dim=c(n.sex,n.age,n.gear))
  #Determine
  
  g <- 1
  for(g in 1:n.gear) {
    fe <- Fmort[g]
    
    h <- 1
    for(h in 1:n.sex) {
      fage <- va[h,,g]*fe  
      za[h,] <- mx[h,] + fage
    }#next h
    #Survival or mortality 
    sa <- exp(-za)
    oa <- 1.0 - sa
    #
    
    
    pa[,,g] <- va[,,g] * oa / za
    qa[,,g] <- va[,,g] * wa * oa / za
    ra[,,g] <- va[,,g] * fa * oa / za
  }#next g
  
  
  
}

for(h in 1:nsex)
{
  # print(fe)
  if(dim(va)[3] > 1){
    fage   <- rowSums(fe*va[h,,])
  }
  else if(dim(va)[3] == 1){
    fage   <- fe * va[h,,]
  }
  za[h,] <- ageSc$mx[h,] + fage
}
sa <- exp(-za)
oa <- 1.0 - sa

# per recruit yield & numbers for gear k
for(k in 1:ngear)
{
  pa[,,k] <- va[,,k] * oa / za
  qa[,,k] <- va[,,k] * wa * oa / za
  ra[,,k] <- va[,,k] * fa * oa / za
}
  output <- NULL
  return(output)
}