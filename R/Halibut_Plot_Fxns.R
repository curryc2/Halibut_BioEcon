#************************************************************************************************
#Project Name: HALIBUT PSC BIOECONOMIC MODEL - Plotting Functions
#Creator: Curry J. Cunningham, NOAA/NMFS
#Date: 12.19.16
#
#PURPOSE: Plot data and results from Halibut Biomass dynamics model
#
#*************************************************************************************************
#NOTES:
#
#*************************************************************************************************

#PRELIMINARY PLOTS
#=================================================================
#' Plot landings data for halibut fishery by sector and year
#'
#'  @param dpi Resolution for plots
#'
#'  @export
plot.landings <- function(dpi=300) {
  require(ggplot2)
  require(gridExtra)
  require(cowplot)
  require(reshape2)
  
  #  Explore Landings and Survey and Assessment data
  assess.dat <- data.frame(read.table('data/HalibutAssessment.dat', header=TRUE))
  lands.dat <- data.frame(read.table('data/HalibutLandings.dat', header=TRUE))
  
  #Stacked area of landings over time
  lands.list <- melt(lands.dat, id.vars='Year')
  names(lands.list) <- c('Year','Sector', 'Landings')
  
  g <- ggplot(lands.list, aes(x=Year, y=Landings, fill=Sector)) +
    theme_gray() +
    geom_area(alpha=0.75)
  g
  ggsave('plots/Landings.pdf', height=4, width=6, dpi=dpi)
  
  #As proportion
  lands.prop <- lands.dat
  for(i in 1:nrow(lands.prop)) {
    lands.prop[i,2:5] <- lands.prop[i,2:5]/sum(lands.prop[i,2:5])
  }#next i
  prop.list <- melt(lands.prop, id.vars='Year')
  names(prop.list) <- c('Year','Sector', 'Proportion')
  
  g2 <- ggplot(prop.list, aes(x=Year, y=Proportion, fill=Sector)) +
    theme_gray() +
    geom_area(alpha=0.75)
  # g2
  
  plot_grid(g,g2, nrow=2, ncol=1)
  ggsave('plots/Landings2.pdf', height=6, width=6, dpi=dpi)
  
  #Plot Indices of Abundance
}



