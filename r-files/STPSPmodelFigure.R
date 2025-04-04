## Playground Pre-synaptic plasticity
baseDir <- '~/Dropbox/'
source('~/Dropbox/r-compile/rUtils/R/rUtils.R')
pp <- function(){paste(baseDir, 'ampOdd_click/plots/',sep='')}
source(paste(baseDir, 'RSkernel/r-files/functions/source_all_functions.R',sep=''))

library(lattice)
# three-parameter model of presynaptic plasticiy
# R  : currently available transmitter
# U  : fraction of available trasmitter being released per tone
# tau: time-constant of re-uptake
# tres : temporal resolution of simulation 

U     <- 0.5
tau  <-  0.0054

tres  <- 0.01
Tto   <- 25


runSimulation <- function(tone.at,amp,U,baseName='psp_example',extn='',eps=F,tres=0.01,Tto=25){
  taxis <- seq(0,to=Tto, by=tres)
  
  R     <- rep(NA,length(taxis))
  dR    <- rep(NA,length(taxis))
  R[1]  <- 1
  dR[1] <- 0
  
  tone.ind <- as.numeric(taxis %in% tone.at)
  tone.ind[tone.ind>.1] <- amp
  
  for (tnd in 1:length(taxis)){
    dR[tnd] <- (1-R[tnd])*tau - U*tone.ind[tnd]*R[tnd]
    if(tnd<length(taxis))
      R[tnd+1] <- (R[tnd]+dR[tnd])%>=%0
  }
  
  fileName <- paste(baseName,'tau',tau,'U',U,extn,sep='_')
  
  if(eps)
    eps.file( paste(fileName,'.eps',sep=''), pdf=T,ratio=0.66, s=5 )
  
  plot(    taxis, 1-(dR/max(abs(dR))), ylim=c(0,2), col='red', type='l',lwd=3 )
  points(  taxis, R, lwd=3, type='l' )
  points(  tone.at, rep(1,length(tone.at)), pch=20, col='blue',cex=2 )
  
  if(eps)
    dev.off()
}

disc <- runSimulation( seq(1,by=1,to=24), c(rep(0.5,10),rep(0.8,1),rep(0.5,6),rep(0.2,1),rep(0.5,6)), 1.0, eps=F)
disc <- runSimulation( seq(1,by=1,to=24), c(rep(0.5,10),rep(0.8,1),rep(0.5,6),rep(0.2,1),rep(0.5,6)), 0.50, eps=T)

disc <- runSimulation( seq(1,by=1,to=21), c(rep(0.2,10),rep(0.5,1),rep(0.2,10)), 0.50, eps=F)
disc <- runSimulation( seq(1,by=1,to=21), c(rep(0.8,10),rep(0.5,1),rep(0.8,10)), 0.50, eps=F)


disc <- runSimulation( c(0,4,6,7, 7.5, 7.75, 8, 8.5, 9.5, 11.5, 15.5),rep(1,11), 0.10, eps=F)
disc <- runSimulation( c(0,4,6,7, 7.5, 7.75, 8, 8.5, 9.5, 11.5, 15.5), 0.95, eps=F)

