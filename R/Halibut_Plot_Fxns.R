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



#' Plot growth and allometry characteristics 
#'
#' @param ageSchedules Object containing age schedules calculated by getAgeSchedules function
#' @param dpi Dots per inch specification for .png figures
#'
#' @export
plot.growth_allometry <- function(ageSchedules, dpi=500) {
  require(ggplot2)
  require(reshape2)
  

  #Create Directory in plots
  plot.dir <- 'plots/Growth and Allometry'
  dir.create(plot.dir)
  
  ages <- ageSchedules$theta$age
  
  ####################
  #Plot Length@Age
  lengths <- ageSchedules$ageSc$la
  dimnames(lengths) <- list(c('Female','Male'), ages)
  la.df <- data.frame(ages, t(lengths))
  la.list <- melt(la.df, id.vars=list('ages'))
  names(la.list) <- c('Age','Sex','value')
  #Plot
  plt.la <- ggplot(la.list, aes(x=Age, y=value, color=Sex)) +
              theme_gray() +
              geom_line(lwd=1.5) +
              geom_point(pch=21, fill='black', size=1.5) +
              ylab('Length (cm)')
  # plt.la
  ggsave(paste0(plot.dir,'/Length_at_Age.png'), plot=plt.la, height=4, width=6, units='in', dpi=dpi)
  
  ####################
  #Plot Weight@Age
  weights <- ageSchedules$ageSc$wa
  dimnames(weights) <- list(c('Female','Male'), ages)
  wa.df <- data.frame(ages, t(weights))
  wa.list <- melt(wa.df, id.vars=list('ages'))
  names(wa.list) <- c('Age','Sex','value')
  #Plot
  plt.wa <- ggplot(wa.list, aes(x=Age, y=value, color=Sex)) +
              theme_gray() +
              geom_line(lwd=1.5) +
              geom_point(pch=21, fill='black', size=1.5) +
              ylab('Weight (pounds)')
  
  ggsave(paste0(plot.dir,'/Weight_at_Age.png'), plot=plt.wa, height=4, width=6, units='in', dpi=dpi)
  
  ####################
  #Plot Length and Weight
  la.list.2 <- cbind(la.list,'Length (cm)')
  wa.list.2 <- cbind(wa.list,'Weight (pounds)')
  names(la.list.2)[4] <- 'variable'
  names(wa.list.2)[4] <- 'variable'
  
  la.wa.list <- rbind(la.list.2, wa.list.2)
  
  plt.la.wa <- ggplot(la.wa.list, aes(x=Age, y=value, color=Sex)) +
                 theme_gray() +
                 geom_line(lwd=1.5) +
                 geom_point(pch=21, fill='black', size=1.5) +
                 facet_wrap(~variable, ncol=1, scales='free')
  ggsave(paste0(plot.dir,'/Length_Weight_at_Age.png'), plot=plt.la.wa, 
           height=6, width=6, units='in', dpi=dpi)
                 
  ####################
  #Plot Matruity@Age
  Maturity <- ageSchedules$ageSc$ma[1,]
  Fecundity <- ageSchedules$ageSc$fa[1,]
  reprod.df <- data.frame(ages, Maturity, Fecundity)
  reprod.list <- melt(reprod.df, id.vars=list('ages'))
  names(reprod.list)[1] <- 'Age'
  #Plot
  plt.reprod <- ggplot(reprod.list, aes(x=Age, y=value)) +
    theme_gray() +
    # geom_area(aes(y=value),  fill='red', alpha=0.25)
    geom_line(lwd=1.5, colour='blue') +
    geom_point(pch=21, fill='red', size=1.5) +
    facet_wrap(~variable, scales='free', ncol=1)
  ggsave(paste0(plot.dir,'/Maturity.png'), plot=plt.reprod, height=4, width=6, units='in', dpi=dpi)
  
}



#' Plot selectivity components.
#' OVerall selectivity is the joint probability of probability of capturing and
#' probability of retaining, and individual at age... size
#'
#' @param selectivity Object containing selectivity information calculated by getSelectivities function
#' @param dpi Dots per inch specification for .png figures
#' @param pt.blk Boolean for whether interior of points will be black or colored
#'
#' @return
#' @export
#'
#' @examples
plot.selectivity <- function(selectivity, dpi=500, pt.blk=FALSE) {
  require(ggplot2)
  require(reshape2)

  #Create Directory in plots
  plot.dir <- 'plots/Selectivity'
  dir.create(plot.dir)
  
  ages <- selectivity$theta$age
  sectors <- selectivity$MP$sector
  
  ####################
  #Probability of Capture
  probCap <- selectivity$probCap
  dimnames(probCap) <- list(c('Female','Male'), ages, sectors)
  probCap.list <- melt(probCap)
  names(probCap.list) <- c('Sex','Age','Sector','value')
  #Plot
  plt.probCap <- ggplot(probCap.list, aes(x=Age, y=value, color=Sector)) +
    theme_gray() +
    geom_line(lwd=1) +
    facet_wrap(~Sex, ncol=1) +
    ylab('Probability of Capture')
    if(pt.blk==TRUE) { 
      plt.probCap <- plt.probCap + geom_point(pch=21, fill='black', size=1.5)
    }else {
      plt.probCap <- plt.probCap + geom_point(pch=21, colour='black', aes(fill=Sector), size=1.5)
    }

  # plt.probCap
  ggsave(paste0(plot.dir,'/Probability of Capture.png'), plot=plt.probCap, height=4, width=5, units='in', dpi=dpi)
  
  plt.probCap.2 <- ggplot(probCap.list, aes(x=Age, y=value, color=Sex)) +
    theme_gray() +
    geom_line(lwd=1) +
    facet_wrap(~Sector, ncol=1) +
    ylab('Probability of Capture')
    if(pt.blk==TRUE) { 
      plt.probCap.2 <- plt.probCap.2 + geom_point(pch=21, fill='black', size=1.5)
    }else {
      plt.probCap.2 <- plt.probCap.2 + geom_point(pch=21, colour='black', aes(fill=Sector), size=1.5)
    }
  # plt.probCap.2
  ggsave(paste0(plot.dir,'/Probability of Capture 2.png'), plot=plt.probCap.2, height=6, width=5, units='in', dpi=dpi)
  
  ####################
  #Probability of Retain
  probRet <- selectivity$probRetain
  dimnames(probRet) <- list(c('Female','Male'), ages, sectors)
  probRet.list <- melt(probRet)
  names(probRet.list) <- c('Sex','Age','Sector','value')
  #Plot
  plt.probRet <- ggplot(probRet.list, aes(x=Age, y=value, color=Sector)) +
    theme_gray() +
    geom_line(lwd=1) +
    facet_wrap(~Sex, ncol=1) +
    ylab('Probability of Retaining an Individual at Age')
    if(pt.blk==TRUE) { 
      plt.probRet <- plt.probRet + geom_point(pch=21, fill='black', size=1.5)
    }else {
      plt.probRet <- plt.probRet + geom_point(pch=21, colour='black', aes(fill=Sector), size=1.5)
    }
  # plt.probRet
  ggsave(paste0(plot.dir,'/Probability of Retain.png'), plot=plt.probRet, height=4, width=5, units='in', dpi=dpi)
  
  plt.probRet.2 <- ggplot(probRet.list, aes(x=Age, y=value, color=Sex)) +
    theme_gray() +
    geom_line(lwd=1) +
    facet_wrap(~Sector, ncol=1) +
    ylab('Probability of Retaining an Individual at Age')
    if(pt.blk==TRUE) { 
      plt.probRet.2 <- plt.probRet.2 + geom_point(pch=21, fill='black', size=1.5)
    }else {
      plt.probRet.2 <- plt.probRet.2 + geom_point(pch=21, colour='black', aes(fill=Sector), size=1.5)
    }
  # plt.probRet.2
  ggsave(paste0(plot.dir,'/Probability of Retain 2.png'), plot=plt.probRet.2, height=6, width=5, units='in', dpi=dpi)
  
  ####################
  #Selectivity@Age
  sel <- selectivity$selex
  dimnames(sel) <- list(c('Female','Male'), ages, sectors)
  sel.list <- melt(sel)
  names(sel.list) <- c('Sex','Age','Sector','value')
  #Plot
  plt.sel <- ggplot(sel.list, aes(x=Age, y=value, color=Sector)) +
    theme_gray() +
    geom_line(lwd=1) +
    facet_wrap(~Sex, ncol=1) +
    ylab('Fishery Selectivity at Age')
    if(pt.blk==TRUE) { 
      plt.sel <- plt.sel + geom_point(pch=21, fill='black', size=1.5)
    }else {
      plt.sel <- plt.sel + geom_point(pch=21, colour='black', aes(fill=Sector), size=1.5)
    }
  # plt.sel
  ggsave(paste0(plot.dir,'/Selectivity.png'), plot=plt.sel, height=4, width=5, units='in', dpi=dpi)
  
  plt.sel.2 <- ggplot(sel.list, aes(x=Age, y=value, color=Sex)) +
    theme_gray() +
    geom_line(lwd=1) +
    facet_wrap(~Sector, ncol=1) +
    ylab('Fishery Selectivity at Age')
    if(pt.blk==TRUE) { 
      plt.sel.2 <- plt.sel.2 + geom_point(pch=21, fill='black', size=1.5)
    }else {
      plt.sel.2 <- plt.sel.2 + geom_point(pch=21, colour='black', aes(fill=Sector), size=1.5)
    }
  # plt.probRet.2
  ggsave(paste0(plot.dir,'/Selectivity 2.png'), plot=plt.sel.2, height=6, width=5, units='in', dpi=dpi)
  
  
  
  
  ####################
  #Combined Plots
  
  #type_sex_sector
  cap <- cbind(probCap.list,'Capture')
  ret <- cbind(probRet.list, 'Retain')
  sel <- cbind(sel.list, 'Selectivity')
  names(cap)[5] <- 'name'
  names(ret)[5] <- 'name'
  names(sel)[5] <- 'name'
  all.list <- rbind(cap,ret,sel)
  
  plt.1 <- ggplot(all.list, aes(x=Age, y=value, color=Sector)) +
    theme_gray() +
    geom_line(lwd=1) +
    facet_grid(name~Sex) +
    ylab('Probability')
  if(pt.blk==TRUE) { 
    plt.1 <- plt.1 + geom_point(pch=21, fill='black', size=1.5)
  }else {
    plt.1 <- plt.1 + geom_point(pch=21, colour='black', aes(fill=Sector), size=1.5)
  }
  # plt.1
  ggsave(paste0(plot.dir,'/type_sex_sector.png'), plot=plt.1, height=5, width=7, units='in', dpi=dpi)
 
  #sex_type_sector
  plt.2 <- ggplot(all.list, aes(x=Age, y=value, color=Sector)) +
    theme_gray() +
    geom_line(lwd=1) +
    facet_grid(Sex~name) +
    ylab('Probability')
  if(pt.blk==TRUE) { 
    plt.2 <- plt.2 + geom_point(pch=21, fill='black', size=1.5)
  }else {
    plt.2 <- plt.2 + geom_point(pch=21, colour='black', aes(fill=Sector), size=1.5)
  }
  # plt.2
  ggsave(paste0(plot.dir,'/sex_type_sector.png'), plot=plt.2, height=5, width=7, units='in', dpi=dpi)
   
  #type_sector_sex
  plt.3 <- ggplot(all.list, aes(x=Age, y=value, color=Sex)) +
    theme_gray() +
    geom_line(lwd=1) +
    facet_grid(name~Sector) +
    ylab('Probability')
  if(pt.blk==TRUE) { 
    plt.3 <- plt.3 + geom_point(pch=21, fill='black', size=1.5)
  }else {
    plt.3 <- plt.3 + geom_point(pch=21, colour='black', aes(fill=Sex), size=1.5)
  }
  # plt.3
  ggsave(paste0(plot.dir,'/type_sector_sex.png'), plot=plt.3, height=6, width=8, units='in', dpi=dpi)
  
  #sector_type_sex
  plt.4 <- ggplot(all.list, aes(x=Age, y=value, color=Sex)) +
    theme_gray() +
    geom_line(lwd=1) +
    facet_grid(Sector~name) +
    ylab('Probability')
  if(pt.blk==TRUE) { 
    plt.4 <- plt.4 + geom_point(pch=21, fill='black', size=1.5)
  }else {
    plt.4 <- plt.4 + geom_point(pch=21, colour='black', aes(fill=Sex), size=1.5)
  }
  # plt.4
  ggsave(paste0(plot.dir,'/sector_type_sex.png'), plot=plt.3, height=5, width=8, units='in', dpi=dpi)
  
  
}




