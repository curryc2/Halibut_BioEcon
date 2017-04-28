get_fished_survivorship <- function(halibut, F_eq) {
  # F_eq <- rep(0.1,4)

  #Define for easier use down the line
  age <- halibut$theta$age
  H <- halibut$theta$H
  A <- halibut$theta$A  
  mx <- halibut$ageSc$mx
  

  fx <- matrix(0, halibut$theta$H, halibut$theta$A)
  
  s <- 1
  for(s in 1:H) {
    a <- 1
    for(a in age) {
      if(a == min(age)) {
        fx[s,a] <- 1.0/H
      } else {
        fx[s,a] <- fx[s,a-1] * exp(-mx[s,a-1] - sum(F_eq*halibut$selex[s,a-1,]) )
      }
      if(a == A) {
        fx[s,a] <- fx[s,a] / (1-exp(-mx[s,a] - sum(F_eq*halibut$selex[s,a-1,]) ))
      }
    }#next age
  }#next sex
  
  return(fx)
}