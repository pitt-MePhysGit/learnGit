baseDir <- '~/Dropbox/'
#source(paste(baseDir, 'ampOddClick/r-files/source_all_functions.R',sep=''))
source( paste(baseDir,'ampOddClick/r-files/functions/functions_AOCdual.R',sep='') )
source( paste(baseDir,'ampOddClick/r-files/functions/prepDFVprobe.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/expand2chdf.R',sep='') )
source( paste(baseDir,'ampOddClick/r-files/functions/function_VProbe_AOCscalar.R',sep='') )
source( paste(baseDir,'ampOddClick/r-files/functions/functions_ampOdd_presynapticPlasticity.R',sep='') )

library(lattice)

source( paste(baseDir,'ampOddClick/r-files/functions/getFeaturePopulation.R',sep='') )
source( paste(baseDir,'ampOddClick/r-files/functions/combineFeaturePopulation.R',sep='') )

source( paste(baseDir,'ampOddClick/r-files/figures/showTD.R',sep='') )
source( paste(baseDir,'ampOddClick/r-files/figures/showTD2.R',sep='') )
source( paste(baseDir,'ampOddClick/r-files/figures/showTD_delta.R',sep='') )

source( paste(baseDir,'r-compile/rUtils/R/rUtils.R',sep='') )


colClasses <- c('character',rep('numeric',2),'character',rep('numeric',5)  ,rep('character',3),rep('numeric',3) )
df         <- read.table(paste(baseDir,'ampOddClick/ampOddclickVProbe_link.txt',sep=''),header=T,colClasses=colClasses)
df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))
df$exp     <- df$session
df         <- prepDFVprobe(df)

existFeatureMat <- function(nx){ file.exists(paste('~/Dropbox/ampOddClick/Rdata/',df$session[nx],'/features.rda',sep='')) }
fME <- sapply(1:dim(df)[1],existFeatureMat )


animal <- 'Sam'
dataSetList <- c('g','r')
valnd <- which(df$drug=='N' & df$dataSet%in%dataSetList & df$animal==animal & fME)
dataSet <- dataSetList[1]

if (length(dataSetList)==2)
  dataSet <- paste(dataSetList[1],dataSetList[2],sep='')
if (length(dataSetList)==3)
  dataSet <- paste(dataSetList[1],dataSetList[2],dataSetList[3],sep='')

##  paths and such
#rda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
rda<- function(){'~/Dropbox/ampOddClick/Rdata/'}
trda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
spd  <- function(){'/Volumes/Drobo5D3/EEG/EEGLab/ampOddClick/'}
options(warn=0)


.GlobalEnv$pp    <- function(){paste(baseDir, 'ampOddClick/word/VProbe/plots/',dataSet,'_',animal,'/',sep='')}
pp <- .GlobalEnv$pp
if (!exists(pp()))
  system(paste('mkdir -p ', pp(), sep=''))

if (!exists(paste( pp(), '/STPSP/', sep='')))
  system(paste('mkdir -p ', pp(),'/STPSP/', sep=''))


slList     <- data.frame( what='H',    trng='t',  xtn='' , stringsAsFactors=FALSE)
slList[2,] <- data.frame( what='G',    trng='t',  xtn='' , stringsAsFactors=FALSE)
slList[3,] <- data.frame( what='B',    trng='t',  xtn='' , stringsAsFactors=FALSE)

slList[4,] <- data.frame( what='H',    trng='e',  xtn='' , stringsAsFactors=FALSE)
slList[5,] <- data.frame( what='G',    trng='e',  xtn='' , stringsAsFactors=FALSE)
slList[6,] <- data.frame( what='B',    trng='e',  xtn='' , stringsAsFactors=FALSE)
slList[7,] <- data.frame( what='A',    trng='e',  xtn='' , stringsAsFactors=FALSE)

slList[8,]  <- data.frame( what='H',   trng='i',   xtn=''  , stringsAsFactors=FALSE)
slList[9,]  <- data.frame( what='G',   trng='i',   xtn=''  , stringsAsFactors=FALSE)
slList[10,] <- data.frame( what='B',   trng='i',   xtn=''  , stringsAsFactors=FALSE)
slList[11,] <- data.frame( what='A',   trng='i',   xtn=''  , stringsAsFactors=FALSE)
slList[12,] <- data.frame( what='T',   trng='i',   xtn=''  , stringsAsFactors=FALSE)

slList[13,] <- data.frame( what='H',   trng='r',   xtn=''  , stringsAsFactors=FALSE)
slList[14,] <- data.frame( what='G',   trng='r',   xtn=''  , stringsAsFactors=FALSE)
slList[15,] <- data.frame( what='B',   trng='r',   xtn=''  , stringsAsFactors=FALSE)
slList[16,] <- data.frame( what='A',   trng='r',   xtn=''  , stringsAsFactors=FALSE)
slList[17,] <- data.frame( what='D',   trng='r',   xtn=''  , stringsAsFactors=FALSE)

slList[18,] <- data.frame( what='H',   trng='b',   xtn=''  , stringsAsFactors=FALSE)
slList[19,] <- data.frame( what='G',   trng='b',   xtn=''  , stringsAsFactors=FALSE)
slList[20,] <- data.frame( what='B',   trng='b',   xtn=''  , stringsAsFactors=FALSE)
slList[21,] <- data.frame( what='A',   trng='b',   xtn=''  , stringsAsFactors=FALSE)
slList[22,] <- data.frame( what='D',   trng='b',   xtn=''  , stringsAsFactors=FALSE)

slList[23,] <- data.frame( what='M',   trng='t',   xtn=''  , stringsAsFactors=FALSE)
slList[24,] <- data.frame( what='M',   trng='e',   xtn=''  , stringsAsFactors=FALSE)
slList[25,] <- data.frame( what='M',   trng='i',   xtn=''  , stringsAsFactors=FALSE)

slList <- slList[23:25,]
xppFeatureList <- lapply( 1:dim(slList)[1], function(fx){ combineFeaturePopulation( valnd,  df=df, removeBaseline=T, signal=slList$what[fx], epoch=slList$trng[fx] )} )


SOA_AMP_plot <- function(DV, SOA, AMP, ylim=NA, DVstr='DV',add=F){
  mnN1SOA <- tapply(DV, SOA, na.mean)
  seN1SOA <- tapply(DV, SOA, na.se)
  soax    <- tapply(SOA,SOA, mean)
  
  mnN1Amp <- tapply(DV, AMP, na.mean)
  seN1Amp <- tapply(DV, AMP, na.se)
  ampx    <- tapply(AMP,AMP, mean)
  
  if (is.na(ylim))
    ylim    <- range( c(mnN1SOA[1:6],mnN1Amp) ) + c(-.1,.1)*diff(range( c(mnN1SOA[1:6],mnN1Amp) ))
  
  errorplot(1:length(soax),mnN1SOA, ey=seN1SOA, pch=20, col='black', ylim=ylim, main=DVstr,add=add )
  errorplot(ampx,mnN1Amp, ey=seN1Amp, pch=20, add=T, col='red' )
}



addFeature <- function(fx,DVextn='', offsetCrit=c(.01,100),pCrit=c(10e-10), ylim=NA, add=F){
  DVstr     <- paste(slList$what[fx],'_',slList$trng[fx],DVextn,sep='')
  xpp       <- xppFeatureList[[fx]]
  subsBin   <- xpp$offset>offsetCrit[1] & xpp$offset<offsetCrit[2] & xpp$pValoff<pCrit[1] & xpp$valTrial
  
  if (sum(subsBin)>1){
    tone             <- list(xpp=xpp[xpp$offset>offsetCrit[1] & xpp$offset<offsetCrit[2] & xpp$pValoff<pCrit[1] & xpp$valTrial,], 
                             exp=paste(animal,dataSet,slList$what[fx],slList$trng[fx],sep='_'), 
                             DVstr=paste(slList$what[fx],'_',slList$trng[fx],DVextn,sep=''),
                             DVstrOff=paste(slList$what[fx],'_',slList$trng[fx],DVextn,'_off',sep=''))
    tone$xpp$exp     <- tone$DVstr
    SOA_AMP_plot( tone$xpp[[tone$DVstr]],tone$xpp$isiOct, tone$xpp$ampInd, DVstr=tone$DVstr, ylim=ylim, add=add )
  }
  
  if (sum(subsBin)<=1){
    plot(NA, xlim=c(0,1),ylim=c(0,1),main=DVstr)
  }
}


addFeaturePop <- function(fx,DVextn='', offsetCrit=c(.01,100),pCrit=c(10e-10,1,1), ylim=NA,mxUnits=64){
  DVstr     <- paste(slList$what[fx],'_',slList$trng[fx],DVextn,sep='')
  
  xpp       <- xppFeatureList[[fx]]
  xpp$unit  <- cumsum(xpp$firstTrial==1)
  ppp       <- xpp[which(xpp$firstTrial==1),]
  
  subsBin   <- ppp$offset>offsetCrit[1] & ppp$offset<offsetCrit[2] & ppp$pValoff<pCrit[1] & (ppp$pValISI<pCrit[2] | ppp$pValAmp<pCrit[3])
  subsInd   <- which(subsBin)
  
  addFeatureUnit <- function(ux){
    txpp <- xpp[ xpp$unit==ux & xpp$valTrial ,]
    SOA_AMP_plot( txpp[[DVstr]],txpp$isiOct, txpp$ampInd, DVstr=txpp$exp[1], ylim=ylim )
  }
  
  if (length(subsInd)>mxUnits)
    subsInd <- sample(subsInd,mxUnits, replace=F)
  
  Nunits <- length(subsInd)
  Nrow   <- 8
  Ncol   <- ceiling(Nunits/Nrow)
  
  par(mfrow=c(Ncol, Nrow), mar=c(3,2,2,1))
  disc <- lapply(subsInd, addFeatureUnit)
  
}



if (length(xppFeatureList)==25){
  eps.file('Figure5a_Scalars_SOA_AMP.eps',pdf=T,s=1,mf=c(5,5),ratio=1)
  par(mfcol=c(5,5))
  lapply(1:length(xppFeatureList), addFeature, DVextn='')
  dev.off()
  
  eps.file('Figure5a_Scalars_SOA_AMP_nrm.eps',pdf=T,s=1,mf=c(5,5),ratio=1)
  par(mfcol=c(5,5))
  lapply(1:length(xppFeatureList), addFeature, DVextn='_nrm',ylim=c(0.2,1.8))
  dev.off()
}

## ===================================================================
source( paste(baseDir,'ampOddClick/r-files/functions/functions_ampOdd_presynapticPlasticity.R',sep='') )

fitFeature <- function(fx, DVextn='_nrm', startFromScratch=F){
  
  xpp              <- xppFeatureList[[fx]]
  tone             <- list(xpp=xpp[xpp$offset> .01 & xpp$offset<100.04 & xpp$pValoff<10e-10 ,], 
                           exp=paste(animal,dataSet,slList$what[fx],slList$trng[fx],sep='_'), 
                           DVstr=paste(slList$what[fx],'_',slList$trng[fx],DVextn,sep=''),
                           DVstrOff=paste(slList$what[fx],'_',slList$trng[fx],DVextn,'_off',sep=''))
  
  tone$xpp$ampIndFactor <- as.factor(tone$xpp$ampInd)
  #tone        <- list(xpp=xpp[xpp$offset> .01 & xpp$offset<100.04 & xpp$pValoff<10e-10,], exp=paste(animal,this_signal,this_epoch,sep='_'), DVstr=paste(this_signal,this_epoch,'nrm',sep='_'),DVstrOff=paste(this_signal,this_epoch,'nrm_off',sep='_') )
  tone$xpp[[tone$DVstrOff]] <- tone$xpp[[tone$DVstr]]
  tone$xpp$exp <- as.factor(tone$exp)
  
  eps <- T
  # range of increasingly more complex and accurate models:
  
  #tst_D_O  <- callSTPSP(tone,   DVstr=tone$DVstrOff, lmOffset=T, gd=T, show=T, showCh=1, who='all', BYstr='exp', FUN=ampFUN_D, minISI=.190, maxISI=21, startFromScratch=startFromScratch, eps=F)
  tst_U1_O <- callSTPSP(tone,   DVstr=tone$DVstrOff, lmOffset=T, gd=T, show=T, showCh=1, who='all', BYstr='exp', FUN=ampFUN_U1, minISI=.190, maxISI=21, startFromScratch=startFromScratch, eps=eps)
  tst_U0_O <- callSTPSP(tone,   DVstr=tone$DVstrOff, lmOffset=T, gd=T, show=T, showCh=1, who='all', BYstr='exp', FUN=ampFUN_U0, minISI=.190, maxISI=21, startFromScratch=startFromScratch, eps=eps)
  
  # depletion independent of amplitude
  tstscl <- callSTPSP(tone,   DVstr=tone$DVstr, lmOffset=F, gd=T, show=T, showCh=1, who='all', BYstr='exp', FUN=ampFUNscl, minISI=.190, maxISI=21, startFromScratch=startFromScratch, eps=eps)
  tst    <- callSTPSP(tone,   DVstr=tone$DVstr, lmOffset=F, gd=T, show=T, showCh=1, who='all', BYstr='exp', FUN=ampFUN,    minISI=.190, maxISI=21, startFromScratch=startFromScratch, eps=eps)
  
  # adding an offset:
  tstsclO <- callSTPSP(tone,   DVstr=tone$DVstrOff, lmOffset=T, gd=T, show=T, showCh=1, who='all', BYstr='exp', FUN=ampFUNscl, minISI=.190, maxISI=21, startFromScratch=startFromScratch, eps=eps)
  tstO    <- callSTPSP(tone,   DVstr=tone$DVstrOff, lmOffset=T, gd=T, show=T, showCh=1, who='all', BYstr='exp', FUN=ampFUN,    minISI=.190, maxISI=21, startFromScratch=startFromScratch, eps=eps)
  
  
  eps.file(paste('/STPSP/STPSPfit_',tone$exp,'_','compareRSS_models.eps',sep=''),pdf=T,s=2,ratio=3/2)
  par(mfcol=c(1,1))
  plot(c(#tst_U0_O$modelFit$modelTest$RSS[2],
    tst_U1_O$modelFit$modelTest$RSS[2],
    tstscl$modelFit$modelTest$RSS[2],
    tst$modelFit$modelTest$RSS[2],
    tstsclO$modelFit$modelTest$RSS[2],
    tstO$modelFit$modelTest$RSS[2]),
    ylab='RSS', type='b', pch=20)
  dev.off()
  
  eps.file(paste('/STPSP/STPSPfit_',tone$exp,'_','compareRSS_ResDf_models.eps',sep=''),pdf=T,s=2,ratio=3/2)
  par(mfcol=c(1,1))
  plot(c(tst_U1_O$modelFit$modelTest$RSS[2]/tst_U1_O$modelFit$modelTest$Res.Df[2],
         tstscl$modelFit$modelTest$RSS[2]/tstscl$modelFit$modelTest$Res.Df[2],
         tst$modelFit$modelTest$RSS[2]/tst$modelFit$modelTest$Res.Df[2],
         tstsclO$modelFit$modelTest$RSS[2]/tstsclO$modelFit$modelTest$Res.Df[2],
         tstO$modelFit$modelTest$RSS[2]/tstO$modelFit$modelTest$Res.Df[2]),
       ylab='RSS/Res.Df', type='b', pch=20)
  dev.off()
  
} 

disc <- lapply(1:length(xppFeatureList), fitFeature, DVextn='_nrm', startFromScratch=T)



## ======================================== Figures for MUA separated by exicited vs inhibited


# 1) get a box-plot of significant vs non-significant mean across each DV


addOffsetHist <- function(fx,DVextn='', ylim=NA){
  xpp              <- xppFeatureList[[fx]]
  tone             <- list(xpp, exp=paste(animal,dataSet,slList$what[fx],slList$trng[fx],sep='_'), 
                           DVstr=paste(slList$what[fx],'_',slList$trng[fx],DVextn,sep=''),
                           DVstrOff=paste(slList$what[fx],'_',slList$trng[fx],DVextn,'_off',sep=''))
  ppp              <- xpp[which(xpp$firstTrial==1),]
  
  allDens <- hist( (ppp$offset%<=%2)%>=%-1, breaks=seq(-1,2,by=.1), plot=F)
  sigDens <- hist( (ppp$offset[ppp$pValoff<10e-10]%<=%2)%>=%-1, breaks=seq(-1,2,by=.1), plot=F)
  
  my.boxplot(allDens$counts, allDens$mids, dxfrac=2, xaxt='s', main=tone$DVstr)
  my.boxplot(sigDens$counts, sigDens$mids, dxfrac=2, add=T, col='red')
}

eps.file('FigureX_MUA_Offsets.eps',pdf=T,s=2,mf=c(3,1),ratio=1)
par(mfcol=c(1,3))
disc <- lapply(1:length(xppFeatureList), addOffsetHist, DVextn='')
dev.off()


eps.file('Figure5c_MUA_SOA_AMP.eps',pdf=T,s=2,mf=c(3,2),ratio=1)
par(mfrow=c(2,3))
lapply(1:length(xppFeatureList), addFeature, DVextn='', offsetCrit=c(.01,100),pCrit=c(10e-10))
lapply(1:length(xppFeatureList), addFeature, DVextn='', offsetCrit=c(-100,-.001),pCrit=c(10e-10))
dev.off()

eps.file('Figure5c_MUA_SOA_AMP_nrm.eps',pdf=T,s=2,mf=c(3,2),ratio=1)
par(mfrow=c(2,3))
lapply(1:length(xppFeatureList), addFeature, DVextn='_nrm', offsetCrit=c(.01,100),pCrit=c(10e-10), ylim=c(0,2))
lapply(1:length(xppFeatureList), addFeature, DVextn='_nrm', offsetCrit=c(-100,-.001),pCrit=c(10e-10), ylim=c(-2,0))
dev.off()


addFeaturePop(2,DVextn='_nrm', offsetCrit=c(-100,.001),pCrit=c(10e-10,1,1)        , ylim=NA, mxUnits=32 )
addFeaturePop(3,DVextn='_nrm', offsetCrit=c(-100,.001),pCrit=c(10e-10,10e-3,10e-3), ylim=NA, mxUnits=32 )

# this_signal <- 'M'
# this_epoch  <- 'e'
# xpp              <- xppFeatureList[[ which(slList$what==this_signal & slList$trng==this_epoch) ]]
# xpp$ampIndFactor <- as.factor(xpp$ampInd)
# #tone        <- list(xpp=xpp[xpp$offset> .01 & xpp$offset<100.04 & xpp$pValoff<10e-10,], exp=paste(animal,this_signal,this_epoch,sep='_'), DVstr=paste(this_signal,this_epoch,'nrm',sep='_'),DVstrOff=paste(this_signal,this_epoch,'nrm_off',sep='_') )
# tone$xpp[[paste(tone$DVstr,'off',sep='_')]] <- tone$xpp[[tone$DVstr]]
# tone$xpp$exp <- as.factor(tone$exp)



# 
# this_signal <- 'A'
# this_epoch <- 'e'
# xpp              <- combineFeaturePopulation(valnd, df=df, signal='A',epoch='e')
# xpp$ampIndFactor <- as.factor(xpp$ampInd)
# ppp              <- xpp[which(xpp$firstTrial==1),]
# 
# xpt              <- combineFeaturePopulation(valnd[1:39], df=df, signal=this_signal,epoch='t')
# xpt$ampIndFactor <- as.factor(xpt$ampInd)
# ppt              <- xpt[which(xpp$firstTrial==1),]
# 
# 
# plot( -log10(ppp$pValISI), -log10(ppp$pValAmp), xlim=c(0,10),ylim=c(0,10) )
# plot( sort(-log10(ppp$pValoff)), type='l')#,ylim=c(0,400) )
# plot( sort(ppp$offset[ppp$pValoff<10e-2]), type='l' )
# 
# 
# 
# tone             <- list(xpp=xpp[xpp$offset< -.01 & xpp$pValoff<10e-5,], exp=paste(animal,this_signal,this_epoch,sep='_'), DVstr=paste(this_signal,this_epoch,sep='_'),DVstrOff=paste(this_signal,this_epoch,'off',sep='_') )
# tone             <- list(xpp=xpp[xpp$offset> .01 & xpp$offset<100.04 & xpp$pValoff<10e-10,], exp=paste(animal,this_signal,this_epoch,sep='_'), DVstr=paste(this_signal,this_epoch,sep='_'),DVstrOff=paste(this_signal,this_epoch,'off',sep='_') )
# 
# tone             <- list(xpp=xpt[xpt$offset> .01 & xpp$offset<100.04 & xpt$pValoff<10e-10,], exp=paste(animal,this_signal,this_epoch,sep='_'), DVstr=paste(this_signal,this_epoch,sep='_'),DVstrOff=paste(this_signal,this_epoch,'off',sep='_') )
# tone$xpp$exp     <- tone$DVstr
# 
# 






tst_U1_O$modelFit$modelTest
tstscl$modelFit$modelTest
tst$modelFit$modelTest
tstsclO$modelFit$modelTest$RSS[2]
tstO$modelFit$modelTest$RSS[2]


## ========== compare it to log SOA linear models
y        <- tone$xpp[ tst$res$yind,]
y$logISI <- log2(y$ISI)


# ------------------ first model: ampInd + logISI
lmLL0 <- lm( y[[tone$DVstr]] ~ 1)
lmLL  <- lm( y[[tone$DVstr]] ~ y$ampInd + y$logISI)
lmLLA <- anova(lmLL0,lmLL)

showSTPSPFit( list(y=y,DVstr=tone$DVstr,lm1=lmLL, fileName='dummy', eps=eps ) )


par(mfcol=c(1,1))
plot(c(lmLLA$RSS[2],tst_U1_O$modelFit$modelTest$RSS[2],tstscl$modelFit$modelTest$RSS[2],tst$modelFit$modelTest$RSS[2],tstsclO$modelFit$modelTest$RSS[2],tstO$modelFit$modelTest$RSS[2]),
     ylab='RSS')








anova(tst$res$lm1,tst$res$lm2)
tst$res$lm2

tstO$res$modelTest$`Pr(>F)`[2]
tstsclO$res$modelTest$`Pr(>F)`[2]

tstO$res$modelTest
tstsclO$res$modelTest


tst$modelFit$par
tstO$modelFit$par

tstscl$res$modelTest




## no scaling with amplitude at all
cumFit <- callSTPSP(tone, DVstr='muanrm', gd=T, show=TRUE, showCh=1, who='random', BYstr='exp',
                    FUN=cumsum1, minISI=.190, maxISI=12.8, startFromScratch=T, eps=T)
cumFit$modelFit$par



## scale R/EEGamp with tone amplitude
sclFit <- callSTPSP(tone,                 DVstr='muanrm', gd=T, show=TRUE, showCh=1, who='random', BYstr='exp', 
                    FUN=cumsumscl, minISI=.190, maxISI=12.8, startFromScratch=T, eps=T)
sclFit$modelFit$par


## scale expended Fraction with tone amplitude
ampFit <- callSTPSP(tone,                 DVstr='muanrm', gd=T, show=TRUE, showCh=1, who='random', BYstr='exp', 
                    FUN=cumsumamp, minISI=.190, maxISI=12.8, startFromScratch=F, eps=T)
ampFit$modelFit$par

