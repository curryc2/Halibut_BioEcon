#' Calculate Beverton-Holt Recruitment
#'
#' @param ssb value for spawning stock biomass
#' @param h steepness
#' @param unfished unfished biomass
#'
#' @return
#' @export
#'
beverton_holt_recruit <- function(ssb, h, bo) {
  return(ssb/(1-((5*h-1)/(4*h))*(1-(ssb/(bo)))))
}