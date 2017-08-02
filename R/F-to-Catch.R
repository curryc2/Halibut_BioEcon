F_to_Cb <- funciton(F.input, halibut, N., ...) {
  ### TESTING ###
  F.input <- c(0.1,0.05,0.01,0.01)
  N.input <- N[,5,]
  ##############
  
  h <- 1
  for(h in 1:n.sex) {
    g <- 1
    for(g in 1:n.gear) {
      temp.F <- F.input[g]*va[h,,g]
      # temp.Z <- temp.F + mx[h,a-1]
      temp.Z <- apply(F.input*va[h,,], 1, sum) + mx[h,]
      
      # harvest.n[h,y-1,a-1,g] <- N[h,y-1,a-1] * (F.a[h,y-1,a-1]/Z.a[h,y-1,a-1]) * (1-exp(-1*Z.a[h,y-1,a-1]))
      harvest.n[h,y-1,a-1,g] <- N[h,y-1,a-1] * (temp.F/temp.Z) * (1-exp(-1*temp.Z))
      # harvest.n[h,y-1,a-1,g] <- N[h,y-1,a-1] * (va[h,a-1,g] * (1-exp(-1*(mx[h,a-1]*Z.a[h,y-1,a-1])))) /(Z.a[h,y-1,a-1])
      
      harvest.b[h,y-1,a-1,g] <- harvest.n[h,y-1,a-1,g] * wa[h,a-1]
    }#next g
  }#next h

  
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