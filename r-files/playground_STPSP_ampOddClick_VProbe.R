

baseDir <- '~/Dropbox/'
#source(paste(baseDir, 'ampOdd_click/r-files/source_all_functions.R',sep=''))
source( paste(baseDir,'ampOdd_click/r-files/functions/functions_AOCdual.R',sep='') )
source( paste(baseDir,'ampOdd_click/r-files/functions/prepDFVprobe.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/expand2chdf.R',sep='') )
source( paste(baseDir,'ampOdd_click/r-files/functions/function_VProbe_AOCscalar.R',sep='') )
source( paste(baseDir,'ampOdd_click/r-files/functions/functions_ampOdd_presynapticPlasticity.R',sep='') )

source( paste(baseDir,'r-compile/rUtils/R/rUtils.R',sep='') )


colClasses <- c('character',rep('numeric',2),'character',rep('numeric',5)  ,rep('character',3),rep('numeric',3) )
df         <- read.table(paste(baseDir,'ampOdd_click/ampOddclickVProbe_link.txt',sep=''),header=T,colClasses=colClasses)
df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))
df$exp     <- df$session
df         <- prepDFVprobe(df)

##  paths and such
rda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
trda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
spd  <- function(){'/Volumes/Drobo5D3/EEG/EEGLab/ampOddClick/'}
options(warn=0)


dataSet <- 'g'
animals <- c('Sam')

chdf     <- function_VProbe_AOCscalar(dataSet='g',animals=c('Sam'),who='all',cmpStr='MUA1',dpthStr='All', pthrAmp = 0.05, pthrISI = 0.05)
valnd    <- which( (chdf$pValAmp < .00001 | chdf$pValISI < .00001) & chdf$depth< 0 & chdf$depth> -2000  )
#valnd    <- valnd[1:5]

## read in all valid MUs and normalize their activity
#for (nd in valnd){
getData <- function(nd){
  
  ## ================================== load channel data
  chInd      <- chdf$chNr[nd]
  chStr      <- paste('ch',chInd,sep='')
  mua        <- loadAOC( chdf[nd,] , what='mua', chStr=chStr, trda=trda,blcorr=T  )
  
  mua$xpp$chdfnd <- c(nd)
  
  ## ================================== prep some additional information
  mua$xpp$ISI[ which( mua$xpp$ISI<0)]  <- 20
  mua$xpp$deltaISI                     <- c(NA, diff(log2(mua$xpp$ISI)))
  mua$xpp$absDeltaISI                  <- abs(mua$xpp$deltaISI)
  
  brks                 <- seq(-4.5,by=1,4.5)
  brks[1]              <- -10; brks[length(brks)]<- 10 
  mua$xpp$deltaISIbin  <- cut(mua$xpp$deltaISI,  breaks=brks)
  
  brks                    <- seq(-0.25,by=0.5,4.5)
  brks[length(brks)]      <- 10 
  mua$xpp$absDeltaISIbin  <- cut(mua$xpp$absDeltaISI, breaks=brks)
  
  ## ========================================================= p2p mua
  tmp          <- apply(mua$lfp,c(1,3),max) - apply(mua$lfp,c(1,3),min)
  mua$xpp$p2p  <- apply( tmp, 1, max )
  
  qtmp         <- quantile(mua$xpp$p2p,probs=c(0.5,0.75))
  mxmua        <- qtmp[1]+10*diff(qtmp)
  mua$xpp$val  <- mua$xpp$p2p < mxmua
  
  ## ================================================== extract and normalize mua
  tnd           <- which(mua$taxis > 25 & mua$taxis < 45)
  mua$xpp$mua   <- apply(mua$lfp[,tnd,1], c(1), mean)
  
  mnMUA          <- na.mean(mua$xpp$mua[which(mua$xpp$val)])
  sdMUA          <- na.sd(mua$xpp$mua[which(mua$xpp$val)])
  
  mua$xpp$muanrm <- (mua$xpp$mua - mnMUA)/sdMUA
  
  return(mua$xpp)
}


dataList <- lapply( valnd, getData)
xpp <- dataList[[1]]
for (i in 2:length(dataList))
  xpp <- rbind( xpp, dataList[[i]])

xpp$ampIndFactor <- as.factor(xpp$ampInd)
tone <- list(xpp=xpp, exp=dataSet)



mnN1SOA <- tapply(tone$xpp$muanrm, tone$xpp$isiOct, na.mean)
seN1SOA <- tapply(tone$xpp$muanrm, tone$xpp$isiOct, na.se)
soax    <- tapply(tone$xpp$ISI, tone$xpp$isiOct, mean)
errorplot(log2(soax),mnN1SOA, ey=seN1SOA, pch=20, xlim=c(-2,4) )



source( paste(baseDir,'ampOdd_click/r-files/functions/functions_ampOdd_presynapticPlasticity.R',sep='') )


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

