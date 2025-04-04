## Preliminaries ====
projectName <- 'ampOddClick'
linkFileName <- 'ampOddClick_link.txt'


baseDir <- '~/Dropbox/'
source( paste(baseDir,'RSkernel/r-files/functions/call_fastICA.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/save_fastICA.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/lamiarLFPtheory/showICAcomponents.R',sep='') )
source( paste(baseDir,'RSkernel/r-files//call_fastPCA.R',sep='') )

library(animation)
library(fastICA)


baseDir <- '~/Dropbox/'
source('~/Dropbox/r-compile/SignalProcessoR/R/signalProcessoR.R')
df           <- prepDFsignalProcessor(projectName, linkFileName)


## Select good sessions ====
pp <- function(){'~/Dropbox/ampOddClick/plots_MePhys/'}
tdf <- df[which(df$animal=='Soleil'  & df$ElectrodeConfig=='ebbbbbbby'), ]
ndx <- dim(tdf)[1]


## load the grand averages ====
mxChan      <- (8*128) +1
badChannels <- c(722,724)
allChannels <- 1:mxChan
goodChannels <- allChannels[ !allChannels%in%badChannels ]
chStrTmp    <- paste('ch',goodChannels,sep='')



## try adding countours using the measured electrode locations
ref       <- 'skull'
chPosM    <- getMePhyspos(chStrTmp, type='exact', returnAll=F, layout='5rows', measured=T, ref=ref)
chPosP    <- getMePhyspos(chStrTmp, type='exact', returnAll=F, layout='5rows', measured=F, ref=ref)
GMC_M     <- getMePhysAnatomy( chPosM, ref=ref) 


amp                 <- loadSignalProcessor( tdf[ndx,], chStr=chPosM$chStr, epoch='raw', analysis='ampInd', blrng=c(-15,0))
amp$channelPosition <- chPosM
amp$GMC             <- getMePhysAnatomy( chPosM, ref=ref) 

# re-reference to the average intracranial potential
amp  <- reRefSignalProcessor( amp,  reRefChan = c(33:dim(amp$lfp)[3])  )


# remove mean for all channels (jointly across all trials) 
tLFP    <- three2twoDimSignalProcessor( amp$lfp)
mnLFP   <- array( rep( apply(tLFP,2,mean), each=dim(tLFP)[1]), dim(tLFP)) 
tLFP    <- tLFP - mnLFP
amp$lfp <- two2threeDimSignalProcessor(tLFP, dm=dim(amp$lfp))


## 

disp <- amp
disp$channelPosition[1:33, ] <- chPosP[1:33,]
disp$lfp[,,1:33] <- 2.5*disp$lfp[,,1:33]
disp$GMC <- GMC_M

eps.file('Sal_multiPlot_amp5_13slices.eps',pdf=T,s=10, ratio=2/3)
showSignalProcessor(disp, frml=~ampInd, trng=c(-25,250), axlab = T, scl='global', sclfctr=1, addyaxis=T,lwd=.5,tcex=.5)#, chIndx = 1:33)#,  vlineat=eval.at)
dev.off()

disp$channelPosition[1:33, ] <- chPosM[1:33,]
eval.at <- c(9, 17, 38)#, 65, 85, 130 ) # 30 
eps.file('Sal_map_amp5_13slices.eps',s=6,ratio=1/2)
showSignalProcessorMap(disp, eval.at=eval.at, dataStr = 'lfp', taxisStr = 'taxis', trialIndx=5, style='dots', dotscl=2, mfrow=c(1,3), zlim=c(-400,400), addContour = T, showTraces = T, trng=c(-10,150), colList=function(n)(rep('gray',n)))
dev.off()


eval.at <-c( seq(0,by=1,to=40), seq(42,by=2,to=80), seq(84,by=4,to=250) )
showSignalProcessorGIF(disp, gifName='LFP_13slices_long',       eval.at=eval.at,              zlim=c(-400,400), trng=c(-10,200), trialIndx=5, setmfrow=T ) 

disp <- amp
disp$GMC <- GMC_M

showSignalProcessorGIF(disp, gifName='LFP_13slices_long_zoom',       eval.at=1:40,              zlim=c(-200,200), trng=c(-10,200), trialIndx=5, setmfrow=T ) 
showSignalProcessorGIF(disp, gifName='LFP_13slices',       eval.at=seq( 0, by=1, to=80), zlim=c(-400,400), trng=c(-10,90), trialIndx=5, setmfrow=T ) 




showSignalProcessor(amp, frml=~ampInd, trng=c(-5,80), chIndx = which(amp$channelPosition$slice==0) )




## topo-plots ====
amp5sal <- subs.data(amp, amp$xpp$ampInd==5)
showSignalProcessor(amp5sal, frml=~sessionID, trng=c(-5,80), zlim=c(-250,250))


showSignalProcessorMap(amp5sal,eval.at=c(17,37,70), chIndx = which(ket$channelPosition$slice==7), style='map', zlim=c(-250,250), isolines=T,isostp=50 )
showSignalProcessorMap(amp5sal,eval.at=c(17,37,70), chIndx = which(ket$channelPosition$slice%in%c(6,7)), style='dots', zlim=c(-250,250), showTraces = T, trng=c(-25,150) )


doOneTimePoint <- function(tndx,...){   showSignalProcessorMap(amp5sal,eval.at=tndx,... )}
doOneTimePoint(6, style='dots', zlim=c(-60,60), showTraces = T, trng=c(-5,12), mfrow=c(1,1),setmfrow=F, dotscl=2, dot.colors=potential.colors, colList=function(n)(rep('gray',n)),addContour=T)

saveGIF({
  disc <- sapply( seq(0,by=.1,to=15), doOneTimePoint, style='dots', dotscl=2, zlim=c(-60,60), 
                  showTraces = T, trng=c(-5,15), dot.colors=potential.colors,addContour=T, 
                  mfrow=c(1,1),setmfrow=F,colList=function(n)(rep('gray',n)) )
}, movie.name = paste(pp(),"test.gif",sep=''), interval = 0.1, nmax = 30, 
ani.width = 600, ani.height = 600, loop=1)

