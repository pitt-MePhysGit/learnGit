rm(list=ls())

library(lattice)
baseDir <- '~/Dropbox/'
source( paste(baseDir,'r-compile/rUtils/R/rUtils.R',sep='') )
source( paste(baseDir,'r-compile/ePhysLab/R/ePhysLab.R',sep='') )
source('~/Dropbox/ampOdd_click/r-files/functions/functions_ampOdd_presynapticPlasticity.R')

pp   <- function(){'~/Dropbox/ampOdd_click/Simulations/'}
eps  <- FALSE


## =================================================
## parameters of the model 
tau    <- 0.5
U      <- 0.9
Off    <- 4
thrsh  <- c(0) 
thrsh  <- seq(0,1,length=31)
linear <- T

## ======================================================================
## settings of the experimental design
Ntr      <- 50000
stepBoth <- F # if set to true, both delta pitch and isi will step simultaneously

pp   <- function(){'~/Dropbox/ampOdd_click/Simulations/stepSeparate/'}
if (stepBoth)
  pp   <- function(){'~/Dropbox/ampOdd_click/Simulations/stepJoint/'}

scl     <- .4               # scaling factor that determines distance between tones
ampVal  <- 0.5 + scl * seq(-.4,.4, by=.2)

pChng      <- c(.1,.9)  # regular, random
pChngISI   <- c(.1,.9)  # equal chance of changing ISI in random and regular condition
Nblck      <- length(pChng)

minISI <- .25
tscl   <- 1.0
ISIlist  <- 2^seq(log2(minISI),by=tscl,length=6) 

## ============================================
## implementing the model parameter settings

FRthrsh <- 0.0
prm     <- c(tau, U, Off)


## ================================= contruct data frame
amplist <- 1:length(ampVal)

x            <- data.frame(trialType=rep(seq(1,to=Nblck,by=1), each=Ntr), rnd=runif(Nblck*Ntr), rnd2=runif(Nblck*Ntr), ISI=sample(ISIlist,Nblck*Ntr,replace=T), ampInd=sample(amplist,Nblck*Ntr,replace=T) )
x$pChng      <- pChng[x$trialType] 
x$pChngISI   <- pChngISI[x$trialType] 
stp          <- cumsum(c(1, x$rnd[2:dim(x)[1]] < x$pChng[2:dim(x)[1]] ))
x$ampInd    <- x$ampInd[stp]
stp2         <- stp
if (stepBoth==F){
  stp2 <- cumsum(c(1, x$rnd2[2:dim(x)[1]] < x$pChngISI[2:dim(x)[1]] ))
}
x$ISI       <- x$ISI[stp2]
x$logISI    <- log2(x$ISI)

x$deltaA     <- c(0,diff(x$ampInd))
x$deltaISI   <- c(0, diff( x$logISI  ))

x$ampType  <- x$deltaA == 0

x$ampDeviant <- abs(x$deltaA)>0 
x$isiDeviant <- abs(x$deltaISI)>0 

rpCnt             <- as.numeric(x$ampType)==1
tmp               <- get.epoch( as.numeric(x$toneType)==1  )
tmp$dur           <- tmp$offset - tmp$onset
rpCnt[tmp$offset] <- -tmp$dur  
x$repNr    <- cumsum(rpCnt)

brks   <- c(-0.5, 0.5, 1.5, 2.5, 6.5, 18.8, 100 )
x$repNrBlk <- cut(x$repNr, breaks=brks)


## ================================== simulate the model
#y               <- x
#y$ampInd        <- ampVal[x$ampInd]
#x$R             <- cumsumamp(y$ISI, y$ampInd, prm[1], prm[2])

x$R <- STPSPmodelFun_Amp(x, prm, thrsh=thrsh, linear=linear, returnStr='sum')
  
#disc <- runModelSTPSP(prm, tone, ybin, DVstr='x', BYstr='exp',shuffel='no', 
#                      returnModel=FALSE, show=T, eps=F,chNr=1,fileName='defaultFileName', FUN=cumsum1)


## ================================== preprocess data frame
IOfcn           <- function(x,thr=0){ (x-thr)%>=%0  }
x$FR            <- IOfcn(x$R, thr=FRthrsh) 
x$FR            <- (x$FR-mean(x$FR))/sd(x$FR)


x$ampIndFactor <- as.factor(x$ampInd)

x$Rn       <- lm(R~logISI+ampIndFactor, data=x)$residuals
x$Rp       <- lm(R~ampInd, data=x)$residuals

x$FRn       <- lm(FR~logISI+ampIndFactor, data=x)$residuals
x$FRp       <- lm(FR~ampInd, data=x)$residuals


## ================================================ analyse data
par(mfcol=c(1,4))

mnSOA <- tapply(x$FR, list(x$ampInd,x$ISI), na.mean)
seSOA <- tapply(x$FR, list(x$ampInd,x$ISI), na.se)

#image(ampVal, log2(ISIlist), mnSOA)

ylim <- c(-2.5,2.5) #c(0,2.5)
colList                <- bluered.colors(length(amplist))
xax <- ampVal
xnd <- 1:length(ampVal)
for (ix in 1:(dim(mnSOA)[2]))
  errorplot( xax[xnd], mnSOA[xnd,ix], ey=seSOA[xnd,ix],  type='b', lwd=2,pch=20, col=colList[ix],ylim=ylim, 
             add=ifelse(ix==1,F,T), xlab='Intensity [au]' )


ylim <- c(-2.5,2.5) #c(0,2.5)
colList                <- soa.colors(length(ISIlist))
xax <- log2(ISIlist)
xnd <- 1:length(ISIlist)
for (ix in 1:(dim(mnSOA)[1]))
  errorplot( xax[xnd], mnSOA[ix,xnd], ey=seSOA[ix,xnd],  type='b', lwd=2,pch=20, col=colList[ix],ylim=ylim, 
             add=ifelse(ix==1,F,T), xlab='SOA [au]' )



## ====================== by delta Amplitude
xax <-  tapply(x$deltaA, list(x$deltaA), na.mean)
mnDA <- tapply(x$FR,     list(x$deltaA,x$ampInd), na.mean)
seDA <- tapply(x$FR,     list(x$deltaA,x$ampInd), na.se)
colList <- bluered.colors(5)

for (ix in 1:(dim(mnSOA)[1]))
  errorplot( xax, mnDA[,ix], ey=seDA[,ix],  type='b', lwd=2,pch=20, col=colList[ix],ylim=ylim, 
             add=ifelse(ix==1,F,T), xlab='delta Amp [au]' )

#errorplot(xax, mnDA[,1], ey=seDA[,1], ylim=c(-1.75,1.75) )
#errorplot(xax, mnDA[,2], ey=seDA[,2], col='red',add=T )

## deltaAmp split by isiOct
N1        <- tapply(x$FR,              list(x$ISI,x$deltaA,x$ampInd), mean )[,,3]
N1se      <- tapply(x$FR,              list(x$ISI,x$deltaA,x$ampInd), na.se )[,,3]
colList   <- soa.colors(6)
vi <- 1:6
for (ix in vi ){
  errorplot( xax, N1[ix,], ey=N1se[ix,],  type='b', lwd=2,pch=20, col=colList[ix],ylim=ylim, 
             add=ifelse(ix==1,F,T), xlab='delta Amp [dB]' )
}



