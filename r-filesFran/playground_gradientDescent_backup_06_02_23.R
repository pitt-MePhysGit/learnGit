library(pracma)
library(lattice)
library(ggplot2)
#library(corrplot)

# 1) run gradient descent for one subject
# 2) set up a pipeline to run all subjects automatically



rda <- function(){ return('/Volumes/EEGCyno/RES/ampOddClick/') }

# 1.1 Load data from one subject ====
tn <- load_ampOddClick_Human(Wave='P2', Subject='2005',  Type_signal='EEG', Hemisphere='L',Parcel='A1')
#table(tn$xpp$ISI_break, tn$xpp$clickIntensity)
#tapply(tn$xpp$DV, list(tn$xpp$ISI_break,tn$xpp$p2p<5000), mean)

# run gradient descent
tst <- callSTPSP(tn, DVstr='DV', gd=T, show=TRUE, showCh=1,FUN=ampFUN,maxISI=10, minISI=0.2, maxp2p=5000, DVsdscale=5, maxEOG=Inf, trialType=NA, fitType=NA, manStartPar = NULL, robust = FALSE)



parm <- c(1.946613e+00,  1.385492e-02,  3.302899e+00,  1.259871e+01,  1.000000e+01)
resMMN   <- ampFUN( parm, tn, rep(T,dim(tn$xpp)[1]), DVstr='DV',returnStr='Model', show=TRUE, chNr=1)




