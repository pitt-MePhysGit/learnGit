## this script generates figures for the time-resolved
## ANOVA analyses

efpar <- par()

baseDir <- '~/Dropbox/'
source(paste(baseDir, 'RSkernel/r-files/functions/source_all_functions.R',sep=''))
source( paste(baseDir,'ampOdd_click/r-files/functions/prepDFVprobe.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/expand2chdf.R',sep='') )
source(paste(baseDir, 'ampOdd_click/r-files/source_all_functions.R',sep=''))



colClasses <- c('character',rep('numeric',2),'character',rep('numeric',5)  ,rep('character',3),rep('numeric',3) )
df         <- read.table(paste(baseDir,'ampOdd_click/ampOddclickVProbe_link.txt',sep=''),header=T,colClasses=colClasses)
df         <- prepDFVprobe(df)
df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))
#df$exp     <- df$session

## preliminaries
#rda  <- function(){'~/Dropbox/RSkernel/Rdata/'}
rda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
trda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
spd  <- function(){'/Volumes/Drobo5D3/EEG/EEGLab/ampOddClick/'}
options(warn=0)

basepp     <- function(){paste(baseDir, '/ampOdd_click/word/VProbe/plots/',sep='')}
pp         <- function(){paste(baseDir, '/ampOdd_click/word/VProbe/plots/',sep='')}

dataSet    <- 'g' 
whoList    <- c('all')

#badDays <- c(9,11,12)+1  # what about 7 and 10?
df$use  <- c(1)


animals   <- c('Sam')#,'Sam')#,'Walter')
animalStr <- '_all'
if (length(animals)==1)
  animalStr <- paste('_',animals,sep='')

cleanOnly <- FALSE
cleanStr  <- ''
if (cleanOnly){
  df$use[ which(df$toneArtefact+df$onsetArtefact>0) ] <- 0
  cleanStr <- '_clean'
}


useInd <- which(df$drug=='N' & !df$Layer <= 0 &  df$dataSet==dataSet &  df$animal%in%animals & df$use==1)

#uIdir  <- 'gclean'
uIdir <- paste('g', animalStr, cleanStr,sep='')

if (!file.exists( paste(pp(), uIdir,sep='') ))
  system( paste('mkdir ', pp(), uIdir,sep='') )


df[useInd,]



### ======================
bfh      <- T
Nepochs  <- 1                  ## total number of epochs (61 tests each)
sortP    <- FALSE ## sort according to earliest significant p-value
stat     <- 'Car'
extn     <- '_mua_all'
width    <- 0.050
alpha    <- .01
Nmin     <- 5                   ## minimum number of signficant tests per set in order to be considered significant
eps      <- F
show     <- F
glmdir   <- 'glm_work'
experiment <- 'ampOddClick'

chdf          <- expand2chdf(df[useInd,])
#chdf$depth    <- (chdf$Layer - chdf$chDepthInd) * chdf$spacing
chdf$channel  <- chdf$chIndx


Intercept      <- groupAnalysis('Intercept',extn=extn, stat=stat,lst=chdf,srt=F,alpha=alpha,bfh=bfh,eps=eps,pdf=T,experiment=experiment,
                                Nmin=15,dir=glmdir,show=show,width=width,showPar=F)
show.tres.glm(Intercept,showPar=F,bfh=F)

ampInd      <- groupAnalysis('ampInd',extn=extn, stat=stat,lst=chdf,srt=Intercept$index,alpha=alpha,bfh=bfh,eps=eps,pdf=T,experiment=experiment,
                                Nmin=Nmin,dir=glmdir,show=show,width=width,showPar=F)
show.tres.glm(ampInd,showPar=F,bfh=F)

isiLog      <- groupAnalysis('isiLog',extn=extn, stat=stat,lst=chdf,srt=F,alpha=alpha,bfh=bfh,eps=eps,pdf=T,experiment=experiment,
                             Nmin=Nmin,dir=glmdir,show=show,width=width,showPar=F)

trialNr      <- groupAnalysis('trialNr',extn=extn, stat=stat,lst=chdf,srt=sortP,alpha=alpha,bfh=bfh,eps=eps,pdf=T,experiment=experiment,
                              Nmin=Nmin,dir=glmdir,show=show,width=width,showPar=F)

absDA       <- groupAnalysis('absDeltaAmpInd',extn=extn, stat=stat,lst=chdf,srt=sortP,alpha=alpha,bfh=bfh,eps=eps,pdf=T,experiment=experiment,
                             Nmin=Nmin,dir=glmdir,show=show,width=width,showPar=F)

DA          <- groupAnalysis('deltaAmpInd',extn=extn, stat=stat,lst=chdf,srt=sortP,alpha=alpha,bfh=bfh,eps=eps,pdf=T,experiment=experiment,
                             Nmin=Nmin,dir=glmdir,show=show,width=width,showPar=F)


## ================== result of varND
extn     <- '_mua_all_ND'
ND       <- groupAnalysis('NDfct',extn=extn, stat=stat,lst=chdf,srt=sortP,alpha=alpha,bfh=bfh,eps=eps,pdf=T,experiment=experiment,
                          Nmin=Nmin,dir=glmdir,show=show,width=width,showPar=F)



picFile <- paste('trANOVA',extn,'.eps',sep='')
eps.file(picFile, pdf=T, s=3, ratio=1)
par(mar=c(5.1, 4.1, 4.1, 2.1))
xlim <- c(-20,200)
plot(  Intercept$eval.at,  apply(Intercept$bfhBin,1,na.mean),type='l', col='black', lwd=2,ylim=c(0,1),xlim=xlim,main=extn)
points(ampInd$eval.at,     apply(ampInd$bfhBin,1,na.mean),type='l', col='red',   lwd=2)
points(isiLog$eval.at,     apply(isiLog$bfhBin,   1,na.mean),type='l', col='blue',  lwd=2)
#points(trialNr$eval.at,    apply(trialNr$bfhBin,  1,na.mean),type='l', col='gray',  lwd=2)
#points(absDA$eval.at,      apply(absDA$bfhBin,  1,na.mean),type='l',   col='orange',lwd=2)
points(DA$eval.at,         apply(DA$bfhBin,  1,na.mean),type='l',   col='orange',lwd=2)
abline(h=0);abline(v=0)
dev.off()

picFile <- paste('trANOVA',extn,'_zoom.eps',sep='')
eps.file(picFile, pdf=T, s=3, ratio=1)
par(mar=c(5.1, 4.1, 4.1, 2.1))
xlim <- c(-10,50)
plot(  Intercept$eval.at,  apply(Intercept$bfhBin,1,na.mean),type='l', col='black', lwd=2,ylim=c(0,1),xlim=xlim,main=extn)
points(ampInd$eval.at,  apply(ampInd$bfhBin,1,na.mean),type='l', col='red',   lwd=2)
points(isiLog$eval.at,     apply(isiLog$bfhBin,   1,na.mean),type='l', col='blue',  lwd=2)
#points(trialNr$eval.at,    apply(trialNr$bfhBin,  1,na.mean),type='l', col='gray',  lwd=2)
#points(absDP$eval.at,      apply(absDP$bfhBin,  1,na.mean),type='l',   col='orange',lwd=2)
abline(h=0);abline(v=0)
dev.off()


# take a look at the parameter estimates. 
# could it just be lower power of SOA test?
picFile <- paste('trANOVA',extn,'_beta.eps',sep='')
eps.file(picFile, pdf=T, s=3, ratio=1)
xlim <- c(-10,250)
plot(  Intercept$eval.at,  apply( id(Intercept$parMat),1,na.mean),type='l', col='black', lwd=2,xlim=xlim,main=extn)
points(isiLog$eval.at,     apply( id(isiLog$parMat),   1,na.mean),type='l', col='blue',  lwd=2)
points(ampInd$eval.at,     apply( id(ampInd$parMat),   1,na.mean),type='l', col='red',  lwd=2)
abline(h=0)
dev.off()


# take a look at latency as a function of cortical depth
getLatency <- function(n,prd=Intercept){
  lat <- NA
  sigInd <- which(prd$bfhBin[,n]==T)
  if(length(sigInd)>0)
    lat <- prd$eval.at[ min(sigInd)]
  
  return( lat)
}
chdf$InterceptOnset <- sapply(1:dim(chdf)[1], getLatency, prd=Intercept)
chdf$SOAOnset       <- sapply(1:dim(chdf)[1], getLatency, prd=isiLog)
chdf$FrqOnset       <- sapply(1:dim(chdf)[1], getLatency, prd=ampInd)

# for each session extract most superficial MUA channel
getFirstMUADepth <- function(nd){ 
  frst   <- NA
  sigInd <- which(!is.na(chdf$InterceptOnset[nd]) )
  
  if(length(sigInd)>0)
    frst <- chdf$channel[nd][min(sigInd)]
  return(frst)
}
onsBySession <- tapply(1:dim(chdf)[1], chdf$session, getFirstMUADepth) 
chdf$firstMUADepth <- onsBySession[chdf$session]

#chdf$depth <- -chdf$channel+chdf$firstMUADepth                  # depth relative to first significant mua channel
#chdf$depth <- chdf$spacing * (-chdf$channel+chdf$Layer)         # depth relative to lfp inversion
sigChan    <- chdf$InterceptOnset > 0
mavOnsI     <- get.mav(chdf$InterceptOnset[ which(sigChan) ], chdf$depth[ which(sigChan) ],width=250, show=T, boot=T,nBoot=1000)
sigChan    <- chdf$FrqOnset > 0
mavOnsF     <- get.mav(chdf$FrqOnset[ which(sigChan) ], chdf$depth[ which(sigChan) ],width=250, show=T, boot=T,nBoot=1000)
sigChan    <- chdf$SOAOnset > 0
mavOnsS     <- get.mav(chdf$SOAOnset[ which(sigChan) ], chdf$depth[ which(sigChan) ],width=250, show=T, boot=T,nBoot=1000)


picFile <- paste('trANOVA',extn,'_OnsetByDepth.eps',sep='')
eps.file(picFile, pdf=T, s=4, ratio=1)

par(mar=c(5.1, 4.1, 4.1, 2.1))
ind <- 1:dim(chdf)[1] #c(1:24)+5*24
plot(   chdf$InterceptOnset[ind], jitter(chdf$depth[ind]), pch=20, col='black', xlim=c(-10,50),ylim=c(-2500,1500))
points( chdf$SOAOnset[ind],       chdf$depth[ind], pch=4,  col='blue')
points( chdf$FrqOnset[ind],       chdf$depth[ind], pch=3,  col='red')

points(mavOnsI$mav, mavOnsI$eval.at, type='l',col='black', lwd=2) 
points(mavOnsF$mav, mavOnsF$eval.at, type='l',col='red', lwd=2) 
points(mavOnsS$mav, mavOnsS$eval.at, type='l',col='blue', lwd=2) 
abline(h=0);abline(v=0)
dev.off()




# =====================================  onset density for intercept, soa and freq
sigChan    <- which( !is.na( chdf$InterceptOnset ))
dI         <- density( chdf$InterceptOnset[sigChan], bw=1   )

sigChan    <- which( !is.na( chdf$SOAOnset ))
dS         <- density( chdf$SOAOnset[sigChan], bw=1   )

sigChan    <- which( !is.na( chdf$FrqOnset ))
dF         <- density( chdf$FrqOnset[sigChan], bw=1   )

picFile <- paste('trANOVA',extn,'_OnsetLatencyDensity.eps',sep='')
eps.file(picFile, pdf=T, s=3, ratio=1)
plot(  dI$x, dI$y, col='black', lwd=2, type='l', xlim=c(-10,50) )
points(dS$x, dS$y, col='blue', lwd=2, type='l')
points(dF$x, dF$y, col='red', lwd=2, type='l')
dev.off()


## ============================================== plot average intercept by cortical depth
depthBrks <- c(10,5,0,-5,-10,-15,-20,-25)
depthBrks <- seq(1500, by=-500, to=-2500)
daxis     <- (depthBrks + my.lag(depthBrks,1))/2
scl       <- 100

picFile <- paste('trANOVA',extn,'_InterceptByDepth.eps',sep='')
eps.file(picFile, pdf=T, s=3, ratio=3/2)
plot(c(-1),c(-1), col='white', ylim=range( depthBrks ),xlim=c(-20,150))
for (i in 2:length(depthBrks)){
  tnd <- which(chdf$depth <= depthBrks[i-1] & chdf$depth>depthBrks[i])
  avg <- apply(Intercept$parMat[,tnd],1,mean)
  #points(Intercept$eval.at, mean(chdf$depth[tnd]) + 5*avg, type='l')
  mean.plot( daxis[i] + 7*scl*t(Intercept$parMat[,tnd]), Intercept$eval.at, add=T, widthFUN=na.se, col='black' )
  points(Intercept$eval.at, daxis[i] + 0*avg, type='l', lty=2)
}
abline(v=0); abline(v=c(50,100,150),lty=2)
#abline(v=c(7,12,17,29), lty=3, col='blue')
dev.off()


## ============================================== plot average isiLog parameter estimate by cortical depth
picFile <- paste('trANOVA',extn,'_isiLogByDepth.eps',sep='')
eps.file(picFile, pdf=T, s=3, ratio=3/2)

plot(c(-1),c(-1), col='white', ylim=range( depthBrks ),xlim=c(-20,150) )
for (i in 2:length(depthBrks)){
  tnd <- which(chdf$depth <= depthBrks[i-1] & chdf$depth>depthBrks[i])
  avg <- apply(Intercept$parMat[,tnd],1,mean)
  #points(Intercept$eval.at, mean(chdf$depth[tnd]) + 5*avg, type='l')
  mean.plot( daxis[i] + 20*scl*t(isiLog$parMat[,tnd]), isiLog$eval.at, add=T, widthFUN=na.se, col='blue' )
  points(Intercept$eval.at, daxis[i] + 0*avg, type='l', lty=2)
}
abline(v=0); abline(v=c(50,100,150),lty=2)
abline(v=c(7,12,17,29), lty=3, col='blue')
dev.off()


picFile <- paste('trANOVA_parameterEstimate',extn,'.eps',sep='')
eps.file(picFile, pdf=T, s=3, ratio=1)
par(mar=c(5.1, 4.1, 4.1, 2.1))
avg <- apply(ampInd$parMat,c(1),mean)
mean.plot( t(avg)/max(avg), Intercept$eval.at, add=F, widthFUN=na.se, col='red',lwd=2 , xlim=c(-20,200), ylim=c(-0.15,1.10))
avg <- apply(Intercept$parMat,1,mean)
mean.plot( t(Intercept$parMat)/max(avg), Intercept$eval.at, add=T, widthFUN=na.se, col='black',lwd=2 )
abline(h=0);abline(v=0)
avg <- apply(isiLog$parMat,1,mean)
mean.plot( t(isiLog$parMat)/max(avg), isiLog$eval.at, add=T, widthFUN=na.se, col='blue',lwd=2 )
dev.off()


picFile <- paste('trANOVA_parameterEstimate',extn,'_zoom.eps',sep='')
eps.file(picFile, pdf=T, s=3, ratio=1)
par(mar=c(5.1, 4.1, 4.1, 2.1))
avg <- apply(ampInd$parMat,c(1),mean)
mean.plot( t(avg)/max(avg), Intercept$eval.at, add=F, widthFUN=na.se, col='red',lwd=2 , xlim=c(-5,50), ylim=c(-0.15,1.10))
avg <- apply(Intercept$parMat,1,mean)
mean.plot( t(Intercept$parMat)/max(avg), Intercept$eval.at, add=T, widthFUN=na.se, col='black',lwd=2 )
abline(h=0);abline(v=0)
avg <- apply(isiLog$parMat,1,mean)
mean.plot( t(isiLog$parMat)/max(avg), isiLog$eval.at, add=T, widthFUN=na.se, col='blue',lwd=2 )
dev.off()




## ============================================== plot average frequency parameter estimate by cortical depth
picFile <- paste('trANOVA',extn,'_frequencyByDepth.eps',sep='')
eps.file(picFile, pdf=T, s=3, ratio=3/2)

plot(c(-1),c(-1), col='white', ylim=range( depthBrks ),xlim=c(-20,150) )
for (i in 2:length(depthBrks)){
  tnd <- which(chdf$depth <= depthBrks[i-1] & chdf$depth>depthBrks[i])
  
  avg <- apply(abs(pitchFreq$parMat[,,tnd]),c(1),mean)
  #avg <- avg - mean(avg[1:100])
  tmp <- apply( abs(pitchFreq$parMat[,,tnd]), c(1,3), mean )
  
  #points(Intercept$eval.at, mean(chdf$depth[tnd]) + 5*avg, type='l')
  mean.plot( daxis[i] + 10*scl*t( tmp - mean(avg[1:100]) ), Intercept$eval.at, add=T, widthFUN=na.se, col='red' )
  points(Intercept$eval.at, daxis[i] + 0*avg, type='l', lty=2)
}
abline(v=0); abline(v=c(50,100,150),lty=2)
abline(v=c(7,12,17,29), lty=3, col='blue')
dev.off()





picFile <- paste('trANOVA',extn,'byDepth_zoom.eps',sep='')
eps.file(picFile, pdf=T, s=3, ratio=3/2)
par(mar=c(5.1, 4.1, 4.1, 2.1))
scl       <- 300
plot(c(-1),c(-1), col='white', ylim=range( depthBrks ),xlim=c(-20,150), main=extn )
xlim <- c(-10,50)
for (i in 2:length(depthBrks)){
  tnd <- which(chdf$depth <= depthBrks[i-1] & chdf$depth>depthBrks[i])
  
  points(Intercept$eval.at,  daxis[i] + scl * rep(0,length(Intercept$eval.at)), lty=2, type='l')
  points(Intercept$eval.at,  daxis[i] + scl * rep(1,length(Intercept$eval.at)), lty=2, type='l')
  
  points(Intercept$eval.at,  daxis[i] + scl*apply(Intercept$bfhBin[,tnd],1,na.mean),type='l', col='black', lwd=2)
  points(pitchFreq$eval.at,  daxis[i] + scl*apply(pitchFreq$bfhBin[,tnd],1,na.mean),type='l', col='red',   lwd=2)
  points(isiLog$eval.at,     daxis[i] + scl*apply(isiLog$bfhBin[,tnd],   1,na.mean),type='l', col='blue',  lwd=2)
}
dev.off()



##========================================================================================  plot the average parameter estimate of deltaPitch and absDeltaPitch
picFile <- paste('trANOVA_parameterEstimate_absDeltaPitch',extn,'.eps',sep='')
eps.file(picFile, pdf=T, s=3, ratio=1)
par(mar=c(5.1, 4.1, 4.1, 2.1))
mean.plot( t(DP$parMat), Intercept$eval.at, add=F, widthFUN=na.se, col='black',lwd=2 , xlim=c(-50,250), ylim=c(-0.05,.05))
mean.plot( t(absDP$parMat), Intercept$eval.at, add=T, widthFUN=na.se, col='red',lwd=2 )
abline(h=0);abline(v=0)
dev.off()



## =================================================================================== SSA vs ND time-resolved
sigChan    <- which( !is.na( chdf$InterceptOnset ))

picFile <- paste('trANOVA_parameterEstimate_SSAvsND_sum',extn,'.eps',sep='')
eps.file(picFile, pdf=T, s=2, ratio=3/2)
mean.plot( t(Intercept$parMat), Intercept$eval.at, add=F, widthFUN=na.se, col='black',lwd=2 , xlim=c(-25,150), ylim=c(-0.05,.45))
#mean.plot( t(Intercept$parMat)+1*t(ND$parMat[,1,]), Intercept$eval.at, add=T, widthFUN=na.se, col='gray',lwd=2 )
mean.plot( t(Intercept$parMat)+1*t(ND$parMat[,2,]), Intercept$eval.at, add=T, widthFUN=na.se, col='blue',lwd=2 )
mean.plot( t(Intercept$parMat)+t(ND$parMat[,3,]), Intercept$eval.at, add=T, widthFUN=na.se, col='red',lwd=2 )
#mean.plot( t(ND$parMat[,1,]), Intercept$eval.at, add=T, widthFUN=na.se, col='gray',lwd=2 )
abline(h=0);abline(v=0)
dev.off()


mnDD  <- apply(ND$parMat[,3,sigChan]-ND$parMat[,2,sigChan],1,mean)
mnSSA <- apply(ND$parMat[,2,sigChan],1, mean)
picFile <- paste('trANOVA_parameterEstimate_SSAvsND',extn,'.eps',sep='')
eps.file(picFile, pdf=T, s=2, ratio=3/2)
mean.plot( (t(ND$parMat[,2,sigChan]))/max(mnSSA), Intercept$eval.at, add=F, widthFUN=na.se, col='blue',lwd=2 , xlim=c(-25,150), ylim=c(-0.1,1.1))
mean.plot( (t(ND$parMat[,3,sigChan])-t(ND$parMat[,2,sigChan]))/max(mnSSA), Intercept$eval.at, add=T, widthFUN=na.se, col='red',lwd=2 )
#mean.plot( t(ND$parMat[,1,]), Intercept$eval.at, add=T, widthFUN=na.se, col='gray',lwd=2 )
abline(h=0);abline(v=0)
dev.off()


mnISI <- apply(isiLog$parMat[,sigChan],1, mean)
picFile <- paste('trANOVA_parameterEstimate_SOAvsND',extn,'.eps',sep='')
eps.file(picFile, pdf=T, s=2, ratio=3/2)
mean.plot( t(isiLog$parMat[,sigChan])/max(mnISI), Intercept$eval.at, add=F, widthFUN=na.se, col='black',lwd=2 , xlim=c(-25,150), ylim=c(-0.05,1.1))
mean.plot( (t(ND$parMat[,3,sigChan])-t(ND$parMat[,2,sigChan]))/max(mnDD), Intercept$eval.at, add=T, widthFUN=na.se, col='red',lwd=2 )
#mean.plot( t(ND$parMat[,1,]), Intercept$eval.at, add=T, widthFUN=na.se, col='gray',lwd=2 )
abline(h=0);abline(v=0)
dev.off()


picFile <- paste('trANOVA_parameterEstimate_SOAvsNDvsSSA_3',extn,'.eps',sep='')
eps.file(picFile, pdf=T, s=2, ratio=3/2)
mean.plot( t(isiLog$parMat[,sigChan])/max(mnISI), Intercept$eval.at, add=F, widthFUN=na.se, col='black',lwd=2 , xlim=c(-25,250), ylim=c(-0.15,1.1))
mean.plot( (t(ND$parMat[,2,sigChan]))/max(mnSSA), Intercept$eval.at, add=T, widthFUN=na.se, col='blue',lwd=2)
mean.plot( (t(ND$parMat[,3,sigChan])-t(ND$parMat[,2,sigChan]))/max(mnDD), Intercept$eval.at, add=T, widthFUN=na.se, col='red',lwd=2 )
abline(h=0);abline(v=0)
dev.off()

picFile <- paste('trANOVA_parameterEstimate_SOAvsNDvsSSA_2',extn,'.eps',sep='')
eps.file(picFile, pdf=T, s=2, ratio=3/2)
mean.plot( t(isiLog$parMat[,sigChan])/max(mnISI), Intercept$eval.at, add=F, widthFUN=na.se, col='black',lwd=2 , xlim=c(-25,250), ylim=c(-0.15,1.1))
mean.plot( (t(ND$parMat[,2,sigChan]))/max(mnSSA), Intercept$eval.at, add=T, widthFUN=na.se, col='blue',lwd=2)
#mean.plot( (t(ND$parMat[,3,sigChan])-t(ND$parMat[,2,sigChan]))/max(mnDD), Intercept$eval.at, add=T, widthFUN=na.se, col='red',lwd=2 )
abline(h=0);abline(v=0)
dev.off()

picFile <- paste('trANOVA_parameterEstimate_SOAvsNDvsSSA_1',extn,'.eps',sep='')
eps.file(picFile, pdf=T, s=2, ratio=3/2)
mean.plot( t(isiLog$parMat[,sigChan])/max(mnISI), Intercept$eval.at, add=F, widthFUN=na.se, col='black',lwd=2 , xlim=c(-25,250), ylim=c(-0.15,1.1))
#mean.plot( (t(ND$parMat[,2,sigChan]))/max(mnSSA), Intercept$eval.at, add=T, widthFUN=na.se, col='blue',lwd=2)
#mean.plot( (t(ND$parMat[,3,sigChan])-t(ND$parMat[,2,sigChan]))/max(mnDD), Intercept$eval.at, add=T, widthFUN=na.se, col='red',lwd=2 )
abline(h=0);abline(v=0)
dev.off()



## ==========================================


## SSA and ND by cortical depth
depthBrks <- seq(1500, by=-250, to=-2500)
daxis     <- (depthBrks + my.lag(depthBrks,1))/2
scl       <- 1000

picFile <- paste('trANOVA',extn,'_NDByDepth.eps',sep='')
eps.file(picFile, pdf=T, s=3, ratio=3/2)
plot(c(-1),c(-1), col='white', ylim=range( depthBrks ),xlim=c(-20,250))
for (i in 2:length(depthBrks)){
  tnd <- which(chdf$depth <= depthBrks[i-1] & chdf$depth>depthBrks[i])
  avg <- apply(ND$parMat[,2,tnd],1,mean)
  mean.plot( daxis[i] + scl*t(ND$parMat[,2,tnd]), Intercept$eval.at, add=T, widthFUN=na.se, col='blue' )
  mean.plot( daxis[i] + scl*t(ND$parMat[,3,tnd]), Intercept$eval.at, add=T, widthFUN=na.se, col='red' )
  points(Intercept$eval.at, daxis[i] + 0*avg, type='l', lty=2)
}
abline(v=0); abline(v=c(50,100,150),lty=2)
abline(v=c(7,12,17,29), lty=3, col='blue')
dev.off()



# ## -------------- load sample test
# fileName <- paste(rda(),'glm_work/RSkernel.trStatsCar_',chdf$session[1],'_ch',chdf$channel[1],'_mua_random','.rda',sep='')
# load(fileName)
# 
# # define function to load results from the time-resolved ANOVA
# loadResults <- function(n,extn='_mua_random',what='Pmat',effect='Intercept'){
#   fileName <- paste(rda(),'glm_work/RSkernel.trStatsCar_',chdf$session[n],'_ch',chdf$channel[n],extn,'.rda',sep='')
#   load(fileName)
#   return(res[[what]][,effect])
# }
# 
# tst   <- sapply( 1:dim(chdf)[1], loadResults, extn='_mua_random',what='Pmat',effect='Intercept'    )
# taxis <- res$eval.at
# image( 1:dim(tst)[1], 1:dim(tst)[2], -log(tst)/log(10) )
# image( as.numeric(taxis), 1:dim(tst)[2], -log(tst)/log(10), zlim=c(0,10) )
# 
# 

