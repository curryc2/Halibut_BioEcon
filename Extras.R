#Extra Pieces for Code Exploration

#=====================================================
#Plot Example HCRs
source('R/HCR/HCR-linear.R')
source('R/HCR/HCR-threshold.R')

png('plots/HCR Examples/Ex1.png', units='in', height=5, width=6, res=300)
HCR_linear(curr.SSB=200e6, SSB0=709e6, floor.F=0, ceiling.F=0.2, ascent.range=c(0.2,0.4), plot=TRUE)
dev.off()

png('plots/HCR Examples/Ex2.png', units='in', height=5, width=6, res=300)
HCR_linear(curr.SSB=200e6, SSB0=709e6, floor.F=0.05, ceiling.F=0.2, ascent.range=c(0.2,0.4), plot=TRUE)
dev.off()

png('plots/HCR Examples/Threshold.png', units='in', height=5, width=6, res=300)
HCR_threshold(curr.SSB=200e6, SSB0=709e6, floor.F=0.05, ceiling.F=0.2, threshold=0.4, plot=TRUE)
dev.off()