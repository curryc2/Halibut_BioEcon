#' Calculates initial sex-specific proportions of biomass at age, based on unfished equilibrium
#'
#' @return init.prop - matrix of proportions by sex and age
#' @export
#'
calc_init_age_prop <- function() {
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
  
  init.prop <- init.prop/sum(init.prop)
  return(init.prop)
}