#' Calculate Ricker Recruitment
#'
#' @param ssb value for spawning stock biomass
#' @param h steepness
#' @param bo unfished biomass
#'
#' @return recruitment biomass
#' @export
#'
ricker_recruit <- function(ssb, h, bo) {
  return(ssb*exp(1.25*log(5*h)*(1-(ssb/bo))))
}