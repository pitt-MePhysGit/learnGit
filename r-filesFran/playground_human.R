#
# model human data
# 
# Install/include in path necessary libraries
install.packages("gdata")
library(gdata)
library(lattice)
source( paste('~/Dropbox/r-compile/rUtils/R/rUtils.R',sep='') )

# Define variables and load data
project   <- 'ampOddClick'

pp <- function(){return(paste('~/Dropbox/',project,'/plots/human/',sep=''))}
fileNameM <- 'Amp_singl_trials_EEG_MLR.txt'
fileNameL <- 'Amp_singl_trials_EEG_LLR.txt'
LLR <- read.table( paste('~/Dropbox/',project,'/RData/',fileNameL,sep=''), header=T,sep = "," )
MLR <- read.table( paste('~/Dropbox/',project,'/RData/',fileNameM,sep=''), header=T,sep = "," )

names(LLR)

table(LLR$subjID)
table(LLR$subjID,LLR$sessionNr)
table(LLR$sessionNr, LLR$blockNr,LLR$subjID)

table(LLR$clickIntensity)
# table(LLR$ISI)

# 
LLR$deltaTime                  <- c(20,diff(LLR$absTime))
LLR$deltaTime[LLR$deltaTime<0] <- 20

plot( tapply(LLR$EEG_N1, LLR$clickIntensity, mean), type='b', pch=20)

test <- tapply(LLR$EEG_N1,LLR$clickIntensity, mean, na.rm = TRUE)

# =======================
plot( tapply(LLR$EEG_N1,      LLR$clickIntensity, mean, na.rm = TRUE), type='b', pch=20)
plot( tapply(LLR$MEG_Grad_N1, LLR$clickIntensity, mean), type='b', pch=20)
plot( tapply(LLR$MEG_Mag_N1, LLR$clickIntensity, mean), type='b', pch=20)


brks <- 2^seq(-2,by=1,to=3)
LLR$ISIbin <- cut((LLR$ISI%<=%8)%>=%.25,breaks=brks)
MLR$ISIbin <- cut((MLR$ISI%<=%8)%>=%.25,breaks=brks)

brks <- 2^seq(-2,by=.5,to=3)
LLR$ISIbinfine <- cut((LLR$ISI%<=%8)%>=%.25,breaks=brks)
MLR$ISIbinfine <- cut((MLR$ISI%<=%8)%>=%.25,breaks=brks)

table(LLR$ISIbin)

## ====================== normalize components
NORM_BY_SUBJ <- function(DV, subjID, type='SUBTRACT'){
  subjFct <- as.factor(subjID)
  lm1 <- lm(DV ~ subjFct)
  subjMean <- tapply(DV, subjID, mean)
  
  DVnrm    <- lm1$residuals + mean(subjMean) 
  return(DVnrm)
}

MLR$EEG_Na_nrm <- NORM_BY_SUBJ(MLR$EEG_Na,MLR$subjID)
MLR$EEG_Pa_nrm <- NORM_BY_SUBJ(MLR$EEG_Pa,MLR$subjID)
MLR$EEG_Nb_nrm <- NORM_BY_SUBJ(MLR$EEG_Nb,MLR$subjID)
MLR$EEG_Pb_nrm <- NORM_BY_SUBJ(MLR$EEG_Pb,MLR$subjID)

MLR$MEG_Grad_Na_nrm <- NORM_BY_SUBJ(MLR$MEG_Grad_Na,MLR$subjID)
MLR$MEG_Grad_Pa_nrm <- NORM_BY_SUBJ(MLR$MEG_Grad_Pa,MLR$subjID)
MLR$MEG_Grad_Nb_nrm <- NORM_BY_SUBJ(MLR$MEG_Grad_Nb,MLR$subjID)
MLR$MEG_Grad_Pb_nrm <- NORM_BY_SUBJ(MLR$MEG_Grad_Pb,MLR$subjID)

LLR$EEG_N1_nrm      <- NORM_BY_SUBJ(LLR$EEG_N1,LLR$subjID)
LLR$EEG_P2_nrm      <- NORM_BY_SUBJ(LLR$EEG_P2,LLR$subjID)
LLR$MEG_Grad_N1_nrm <- NORM_BY_SUBJ(LLR$MEG_Grad_N1,LLR$subjID)
LLR$MEG_Grad_P2_nrm <- NORM_BY_SUBJ(LLR$MEG_Grad_P2,LLR$subjID)


SOA_AMP_plot <- function(DV, SOA, AMP, ylim=NA, DVstr='EEG_N1',...){
  mnN1SOA <- tapply(DV, SOA, na.mean)
  seN1SOA <- tapply(DV, SOA, na.se)
  soax    <- tapply(SOA,SOA, mean)
  
  mnN1Amp <- tapply(DV, AMP, na.mean)
  seN1Amp <- tapply(DV, AMP, na.se)
  ampx    <- 1+tapply(AMP,AMP, mean)

  if (is.na(ylim))
    ylim    <- range( c(mnN1SOA,mnN1Amp) ) + c(-.25,.25)*diff(range( c(mnN1SOA,mnN1Amp) ))
  
  errorplot(1:length(soax),mnN1SOA, ey=seN1SOA, pch=20, col='black', ylim=ylim,...)
  errorplot(ampx,mnN1Amp, ey=seN1Amp, pch=20, add=T, col='red' )
}

SOA_AMP_Interactionplot <- function(DV, IV1, IV2, ylim=NA, DVstr='EEG_N1',colFUN=soa.colors,...){
  
  dotargs <- list(...)
  if (! 'swap'%in% names(dotargs) )
    dotargs$swap <- FALSE
  swap <- dotargs$swap  
  
  if (!swap){
    SOA <- IV1
    AMP <- IV2
  }
  
  if (swap){
    SOA <- IV2
    AMP <- IV1
  }
  
  mnx  <- tapply(DV, list(SOA,AMP), na.mean)
  sex  <- tapply(DV, list(SOA,AMP), na.se)
  soax <- tapply(SOA,SOA, mean)

  if (is.na(ylim))
    ylim    <- range( c(mnx) ) + c(-.25,.25)*diff(range( c(mnx) ))
  
  plot(1:dim(mnx)[2], rep(-100,dim(mnx)[2]), ylim=ylim,  main=DVstr)
  for (tx in 1:dim(mnx)[1])
    errorplot(1:dim(mnx)[2],mnx[tx,],ey=sex[tx,], add=T, col=colFUN(dim(mnx)[1])[tx] )
}


SOA_AMP_Interactionplot(LLR$EEG_P2_nrm, LLR$ISIbin, LLR$clickIntensity, ylim=NA, DVstr='EEG_P2', colFUN=soa.colors )
SOA_AMP_Interactionplot(LLR$EEG_P2_nrm, LLR$ISIbin, LLR$clickIntensity, ylim=NA, DVstr='EEG_P2', swap=T, colFUN=bluered.colors )

SOA_AMP_plot_bySubj <- function(DVstr,R=LLR,plotFUN=SOA_AMP_plot,...){
  mnN1SOA <- tapply(R[[DVstr]], list(R$ISIbin, R$clickIntensity,R$subjID), na.mean)
  mndf    <- expand.grid(ISI=dimnames(mnN1SOA)[[1]], clickIntensity=dimnames(mnN1SOA)[[2]],subjID=dimnames(mnN1SOA)[[3]] )
  mndf$DV <- c(mnN1SOA)
  mndf$clickIntensity <- c(mndf$clickIntensity)
  plotFUN(mndf$DV, mndf$ISI, mndf$clickIntensity, DVstr=DVstr, main=DVstr,...)
}


eps.file( paste('SOA_AMP_plot.eps',sep=''),eps=T,pdf=T, s=1.5,ratio=1.33, mf=c(6,2) )
par(mfrow=c(2,6))
SOA_AMP_plot_bySubj(DVstr='EEG_Na_nrm',R=MLR)
SOA_AMP_plot_bySubj(DVstr='EEG_Pa_nrm',R=MLR)
SOA_AMP_plot_bySubj(DVstr='EEG_Nb_nrm',R=MLR)
SOA_AMP_plot_bySubj(DVstr='EEG_Pb_nrm',R=MLR)
SOA_AMP_plot_bySubj(DVstr='EEG_N1_nrm',R=LLR)
SOA_AMP_plot_bySubj(DVstr='EEG_P2_nrm',R=LLR)

SOA_AMP_plot_bySubj(DVstr='MEG_Grad_Na_nrm',R=MLR)
SOA_AMP_plot_bySubj(DVstr='MEG_Grad_Pa_nrm',R=MLR)
SOA_AMP_plot_bySubj(DVstr='MEG_Grad_Nb_nrm',R=MLR)
SOA_AMP_plot_bySubj(DVstr='MEG_Grad_Pb_nrm',R=MLR)
SOA_AMP_plot_bySubj(DVstr='MEG_Grad_N1_nrm',R=LLR)
SOA_AMP_plot_bySubj(DVstr='MEG_Grad_P2_nrm',R=LLR)
dev.off()


eps.file( paste('SOA_AMP_Interactionplot.eps',sep=''),eps=T,pdf=T, s=1.5,ratio=1.33, mf=c(6,2) )
par(mfrow=c(2,6))
SOA_AMP_plot_bySubj(DVstr='EEG_Na_nrm',R=MLR, plotFUN=SOA_AMP_Interactionplot, swap=T, col=bluered.colors)
SOA_AMP_plot_bySubj(DVstr='EEG_Pa_nrm',R=MLR, plotFUN=SOA_AMP_Interactionplot, swap=T, col=bluered.colors)
SOA_AMP_plot_bySubj(DVstr='EEG_Nb_nrm',R=MLR, plotFUN=SOA_AMP_Interactionplot, swap=T, col=bluered.colors)
SOA_AMP_plot_bySubj(DVstr='EEG_Pb_nrm',R=MLR, plotFUN=SOA_AMP_Interactionplot, swap=T, col=bluered.colors)
SOA_AMP_plot_bySubj(DVstr='EEG_N1_nrm',R=LLR, plotFUN=SOA_AMP_Interactionplot, swap=T, col=bluered.colors)
SOA_AMP_plot_bySubj(DVstr='EEG_P2_nrm',R=LLR, plotFUN=SOA_AMP_Interactionplot, swap=T, col=bluered.colors)

SOA_AMP_plot_bySubj(DVstr='EEG_Na_nrm',R=MLR, plotFUN=SOA_AMP_Interactionplot, swap=F, col=soa.colors)
SOA_AMP_plot_bySubj(DVstr='EEG_Pa_nrm',R=MLR, plotFUN=SOA_AMP_Interactionplot, swap=F, col=soa.colors)
SOA_AMP_plot_bySubj(DVstr='EEG_Nb_nrm',R=MLR, plotFUN=SOA_AMP_Interactionplot, swap=F, col=soa.colors)
SOA_AMP_plot_bySubj(DVstr='EEG_Pb_nrm',R=MLR, plotFUN=SOA_AMP_Interactionplot, swap=F, col=soa.colors)
SOA_AMP_plot_bySubj(DVstr='EEG_N1_nrm',R=LLR, plotFUN=SOA_AMP_Interactionplot, swap=F, col=soa.colors)
SOA_AMP_plot_bySubj(DVstr='EEG_P2_nrm',R=LLR, plotFUN=SOA_AMP_Interactionplot, swap=F, col=soa.colors)
dev.off()






## =======================================================
mnpar(mfrow=c(3,6))
SOA_AMP_plot(MLR$EEG_Na, MLR$ISIbin, MLR$clickIntensity, main='Na')
SOA_AMP_plot(MLR$EEG_Pa, MLR$ISIbin, MLR$clickIntensity, main='Pa')
SOA_AMP_plot(MLR$EEG_Nb, MLR$ISIbin, MLR$clickIntensity, main='Nb')
SOA_AMP_plot(MLR$EEG_Pb, MLR$ISIbin, MLR$clickIntensity, main='Pb')
SOA_AMP_plot(LLR$EEG_N1, LLR$ISIbin, LLR$clickIntensity, main='N1')
SOA_AMP_plot(LLR$EEG_P2, LLR$ISIbin, LLR$clickIntensity, main='P2')

SOA_AMP_plot(MLR$MEG_Grad_Na, MLR$ISIbin, MLR$clickIntensity, main='mNa_G')
SOA_AMP_plot(MLR$MEG_Grad_Pa, MLR$ISIbin, MLR$clickIntensity, main='mPa_G')
SOA_AMP_plot(MLR$MEG_Grad_Nb, MLR$ISIbin, MLR$clickIntensity, main='mNb_G')
SOA_AMP_plot(MLR$MEG_Grad_Pb, MLR$ISIbin, MLR$clickIntensity, main='mPb_G')
SOA_AMP_plot(LLR$MEG_Grad_N1, LLR$ISIbin, LLR$clickIntensity, main='mN1_G')
SOA_AMP_plot(LLR$MEG_Grad_P2, LLR$ISIbin, LLR$clickIntensity, main='mP2_G')

SOA_AMP_plot(MLR$MEG_Mag_Na, MLR$ISIbin, MLR$clickIntensity, main='mNa_M')
SOA_AMP_plot(MLR$MEG_Mag_Pa, MLR$ISIbin, MLR$clickIntensity, main='mPa_M')
SOA_AMP_plot(MLR$MEG_Mag_Nb, MLR$ISIbin, MLR$clickIntensity, main='mNb_M')
SOA_AMP_plot(MLR$MEG_Mag_Pb, MLR$ISIbin, MLR$clickIntensity, main='mPb_M')
SOA_AMP_plot(LLR$MEG_Mag_N1, LLR$ISIbin, LLR$clickIntensity, main='mN1_M')
SOA_AMP_plot(LLR$MEG_Mag_P2, LLR$ISIbin, LLR$clickIntensity, main='mP2_M')






## ===================================
SOA_AMP_subject <- function(n,DVstr='EEG_N1',R=LLR){
  ndx <- which(R$subjID == n)
  SOA_AMP_plot(R[[DVstr]][ndx], R$ISIbin[ndx], R$clickIntensity[ndx], main=paste(n,sep=' ') )
}
idlist <- unique(LLR$subjID)

par(mfrow=c(3,4))
lapply(idlist, SOA_AMP_subject, DVstr='EEG_N1', R=LLR)
lapply(idlist, SOA_AMP_subject, DVstr='EEG_N1_nrm', R=LLR)



lapply(idlist, SOA_AMP_subject, DVstr='EEG_P2')

lapply(idlist, SOA_AMP_subject, DVstr='MEG_Grad_Na', R=MLR)



