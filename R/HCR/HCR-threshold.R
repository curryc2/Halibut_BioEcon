HCR_threshold <- function(curr.SSB, SSB0, floor.F, ceiling.F, threshold=0.4, plot=FALSE) {
  #Calculate SSB ration
  ratio.SSB <- curr.SSB/SSB0

  if(ratio.SSB>threshold) {
    curr.F <- ceiling.F
  }else {
    curr.F <- floor.F
  }
  
  if(plot==TRUE) {
    sim.ratio <- seq(0,max(1,ratio.SSB), length.out=100)
    n.sim <- length(sim.ratio)
    sim.F <- vector(length=n.sim)
    
    s <- 1
    for(s in 1:n.sim) {
      if(sim.ratio[s]>threshold) {
        sim.F[s] <- ceiling.F
      }else {
          sim.F[s] <- floor.F
      }
    }#next s
    
    #Plotting Section
    plot(x=NULL, y=NULL, xlim=c(min(sim.ratio), max(sim.ratio)), ylim=c(0, 1.1*max(sim.F)),
         xlab='SSB Ratio', ylab='Fishing Mortality Rate (F)', main='Threshold HCR')
    polygon(x=c(0, threshold, threshold, 0), y=c(0,0,1.5*max(sim.F),1.5*max(sim.F)), 
            col='red', border=FALSE)
    polygon(x=c(threshold, 1, 1, threshold), y=c(0,0,1.5*max(sim.F),1.5*max(sim.F)), 
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