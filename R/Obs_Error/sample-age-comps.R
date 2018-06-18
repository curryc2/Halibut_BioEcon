#' Sample age composition data with multinomial or dirichlet error.
#'
#' @param true.comps Matrix of age composition: rows=number of samples, cols=number of age classes
#' @param Nsamp Vector of sample sizes of length equal to number of samples
#' @param type Type of error, of type: 'multinomial' or 'dirichlet'
#' @param cpar 
#'
#' @return
#' @export
#'
#' @examples
sample_age_comps <- function(true.props, Nsamp, cpar=NULL) {
  require(gtools)
  
  #Initial Checking
  
  
  # cpar - A value of 1 indicates the same standard deviation as a multinomial of the given Nsamp, 2 indicates twice, etc. Values greater than one indicate overdispersion, and less underdispersion.
  
  cpar <- 2
  #True data
  probs <- as.numeric(newcomp[-(1:6)]/sum(newcomp[-(1:6)]))
  ## If cpar is NA this signifies to use the multinomial
  if(is.null(cpar)){ #MULTINOMIAL
    newcomp[-(1:6)] <-
      rmultinom(1, size=newcomp$Nsamp, prob=probs)#/newcomp$Nsamp
  } else { # use Dirichlet
    lambda <- newcomp$Nsamp/cpar[i]^2 - 1
    if(lambda<0)
      stop(paste("Invalid Dirichlet parameter: Lambda=", lambda))
    newcomp[-(1:6)] <- gtools::rdirichlet(1,probs * lambda)
    ## Use the effective sample size when using Dirichlet
    effectiveN <- newcomp$Nsamp/cpar[i]^2
    newcomp$Nsamp <- effectiveN
  }
  return(obs.comp)
}