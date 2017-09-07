#' Harvest Control Rule: Basic linear increase
#'
#' @param curr.SSB Current spawning stock biomass 
#' @param SSB0 Equilibrium spawning stock biomass
#' @param floor.F Fishing mortality rate floor for HCR
#' @param ceiling.F Fishing mortality rate ceiling for HCR
#' @param ascent.range Vector of length 2, describing the 1. lower and 2. upper range of the ascending region of HCR
#' @param plot Boolean for whether to plot HCR example
#'
#' @return
#' @export
#'
HCR_linear <- function(curr.SSB, SSB0, floor.F, ceiling.F, ascent.range=c(0.2,0.4), plot=FALSE) {
  
  ### TESTING
  
  # curr.SSB <- 3e5
  # SSB0 <- 1e6
  # floor <- 0
  # ceiling <- 0.8
  # ascent.range=c(0.2,0.4)
  # plot=TRUE
  
  ###
  
  #Calculate SSB ration
  ratio.SSB <- curr.SSB/SSB0
  slope <- (ceiling.F-floor.F)/(ascent.range[2]-ascent.range[1])
  
  if(ratio.SSB>ascent.range[2]) {
    curr.F <- ceiling.F
  }else {
    if(ratio.SSB<ascent.range[1]) {
      curr.F <- floor.F
    }else {
      #ACTION ZONE
      curr.F <- floor.F + slope * (ratio.SSB - ascent.range[1])
    }
  }
  
  if(plot==TRUE) {
    sim.ratio <- seq(0,max(1,ratio.SSB), length.out=100)
    n.sim <- length(sim.ratio)
    sim.F <- vector(length=n.sim)
    
    s <- 1
    for(s in 1:n.sim) {
      if(sim.ratio[s]>ascent.range[2]) {
        sim.F[s] <- ceiling.F
      }else {
        if(sim.ratio[s]<ascent.range[1]) {
          sim.F[s] <- floor.F
        }else {
          #ACTION ZONE
          sim.F[s] <- floor.F + slope * (sim.ratio[s] - ascent.range[1])
        }
      }
    }#next s
    
    #Plotting Section
    plot(x=NULL, y=NULL, xlim=c(min(sim.ratio), max(sim.ratio)), ylim=c(0, 1.1*max(sim.F)),
           xlab='SSB Ratio', ylab='Fishing Mortality Rate (F)', main='Linear HCR')
    polygon(x=c(0, min(ascent.range), min(ascent.range), 0), y=c(0,0,1.5*max(sim.F),1.5*max(sim.F)), 
              col='red', border=FALSE)
    polygon(x=c(min(ascent.range), max(ascent.range), max(ascent.range), min(ascent.range)), y=c(0,0,1.5*max(sim.F),1.5*max(sim.F)), 
              col='yellow', border=FALSE)
    polygon(x=c(max(ascent.range), 1, 1, max(ascent.range)), y=c(0,0,1.5*max(sim.F),1.5*max(sim.F)), 
            col='green', border=FALSE)
    #Harvest Control Rule
    lines(x=sim.ratio, y=sim.F, lwd=5, col='blue')
    #Current
    segments(x0=ratio.SSB, y0=-10, x1=ratio.SSB, y1=curr.F, lty=3, lwd=2, col='black')
    segments(x0=-10, y0=curr.F, x1=ratio.SSB, y1=curr.F, lty=3, lwd=2, col='black')
    points(x=ratio.SSB, y=curr.F, pch=21, bg='gray', cex=1.25)
  }
  return(curr.F)  
}