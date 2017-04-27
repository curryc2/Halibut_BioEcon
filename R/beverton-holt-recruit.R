#' Calculate Beverton-Holt Recruitment
#'
#' @param ssb value for spawning stock biomass
#' @param h steepness
#' @param ssb_eq equilibrium biomass
#'
#' @return
#' @export
#'
beverton_holt_recruit <- function(ssb, h, ssb_eq) {
  return(ssb/(1-((5*h-1)/(4*h))*(1-(ssb/(ssb_eq)))))
}