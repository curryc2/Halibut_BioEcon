#' Calculate Ricker Recruitment
#'
#' @param ssb value for spawning stock biomass
#' @param h steepness
#' @param ssb_eq equilibrium biomass
#'
#' @return recruitment biomass
#' @export
#'
ricker_recruit <- function(ssb, h, ssb_eq) {
  return(ssb*exp(1.25*log(5*h)*(1-(ssb/ssb_eq))))
}