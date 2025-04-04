## preliminaries
#args <- commandArgs( trailingOnly=TRUE )

baseDir <- '~/Dropbox/'
source(paste(baseDir, 'ampOdd_click/r-files/source_all_functions.R',sep=''))

## preliminaries
pp  <- function(){paste(baseDir, 'ampOdd_click/plots/',sep='')}
rda <- function(){'~/Dropbox/ampOdd_click/Rdata/'}


animalList <- c('Rockey','Jesse','Walter','Sam')
dataSet    <- 'g'
cmpStrList <- c('p1x','p1a','p1b','nx','px','n1','p2','n2')
#cmpStrList <- c('p1x.r','p1a.r','p1b.r','nx.r','px.r','n1.r','p2.r','n2.r')
#cmpStrList <- paste('z.',cmpStrList,sep='')

maxp2p     <- 450
minISI     <- 0.190
maxISI     <- 12.8
DVsdscale  <- 5


chStr       <- NA
getData <- function(subj){
  print(subj)
  tone        <- loadAmpOddClickAll( df[which(df$animal==subj&df$dataSet==dataSet),] , what='raw', chStr=chStr  )
  tone$exp    <- paste(subj,'_',dataSet,sep='')
  return(tone)
}
dataList <- lapply( animalList, getData)


### ============================================================================================
runOneFunctionOneSubject <- function(n, thisDVstr='n1.r', FUN=subjFigure_CmpByDeltaAmpAndAmp,...){
  tn <- dataList[[n]]
  animal <- strsplit(tn$exp,'_')[[1]][1]
  print(animal)
  if(animal=='Rockey') # & what %in% c('raw','mlp','mbp','bsp'))
    maxEOG <- 15000
  
  if(animal=='Walter')# & what %in% c('raw','mlp','mbp','bsp'))
    maxEOG <- 5000
  
  if(animal=='Jesse')# & what %in% c('raw','mlp','mbp','bsp'))
    maxEOG <- 2000
  
  if(animal=='Sam')# & what %in% c('raw','mlp','mbp','bsp'))
    maxEOG <- 15000
  
  ## do the subsetting before calling the fitting function
  if(is.character(tn$xpp$exp))
    tn$xpp$exp <- factor(tn$xpp$exp)
  
  tn$xpp$firstTrialInDay <- diff( c(0,as.numeric(tn$xpp$exp)) )
  
  # exclud trials based on several criteria, including one related to thisDVs
  tn$xpp$DV <- tn$xpp[[thisDVstr]]
  
  x         <- tn$xpp
  mDV       <- quantile( x$DV, pnorm(c(-.5,.5),m=0,sd=1) )# prob=c(0.25,0.5,0.75) )
  valRng    <- mDV[2] + c(-DVsdscale,DVsdscale) * diff(mDV) 
  propAccpt <- length( which( x$DV>valRng[1] &  x$DV<valRng[2] ) )/dim(x)[1]
  print(mDV)
  print(valRng)
  print(propAccpt)
  
  ybin <- ( x$p2p<maxp2p                  & 
              x$ISI>minISI                & 
              x$ISI<maxISI                & 
              x$powEOG<maxEOG             &
              x$firstTrialInDay==0        & 
              x$DV>valRng[1]        &  
              x$DV<valRng[2]        )
  
  tn        <- subs.data(tn,   ybin )
  
  res <- FUN(     tn, DVstr=thisDVstr,...)
  
  return(res)
}

### =========================================================================================== DeltaAmp and Amp
doIt <- function(DVstr){ 
  resList <- list()
  for (ind in 1:length(dataList))
    resList[[ind]] <- runOneFunctionOneSubject(  ind, thisDVstr=DVstr, FUN=subjFigure_CmpByDeltaAmpAndAmp,trialType=c(1),show=F,eps=F )

  disc <- gaFigure_CmpByDeltaAmpAndAmp(resList, DVstr=DVstr, show=T, eps=T)
  return( list(resList=resList) )
}
resAA <- lapply(cmpStrList, doIt)

getPval <- function(n,m=1){return(resAA[[m]]$resList[[n]]$lmresult$Pr[2]) }
getPvalCmp <- function(m){return(sapply(1:4,getPval,m=m)) }
pmat <- sapply(1:length(cmpStrList), getPvalCmp )
image(1:4,1:8,  (-1*log(pmat)/log(10))%<=a=0%log(0.05)/log(10), col=heat.colors(100),xaxt='n', yaxt='n',xlab='l',ylab='l' ) 
axis(1, at=1:4, labels=animalList)
axis(2, at=1:8, labels=cmpStrList)



### =========================================================
doIt <- function(DVstr){ 
  resList <- list()
  for (ind in 1:length(dataList))
    resList[[ind]] <- runOneFunctionOneSubject(  ind, thisDVstr=DVstr, FUN=subjFigure_CmpByDeltaAmp,trialType=c(1),show=F,eps=F )
  
  disc <- gaFigure_CmpByDeltaAmp(resList, DVstr=DVstr, show=T, eps=T)
  return( list(resList=resList) )
}
resA <- lapply(cmpStrList, doIt)

eps.file('GA_AllCmpByDeltaAmp.eps',s=3,ratio=3/2)
for (j in 1:length(cmpStrList) )
  gaFigure_CmpByDeltaAmp(resA[[j]][[1]],DVstr=cmpStrList[j], show=T, eps=F,col=componentColor2(cmpStrList[[j]]),add=ifelse(j==1,F,T))
dev.off()

extractCoeff <- function(n,m=1){ resA[[m]][1]$resList[[n]]$coeff['tn$xpp$deltaAmpInd'] }

getttestcmp <- function(m){ 
  tmp <- sapply(1:4, extractCoeff,m=m)
  tst <- t.test(tmp)
  print(tst)
  return(tst$p.val)
}
sapply(1:8, getttestcmp)


getPval <- function(n,m=1){return(resA[[m]]$resList[[n]]$lmresult$Pr[2]) }
getPvalCmp <- function(m){return(sapply(1:4,getPval,m=m)) }
pmat <- sapply(1:length(cmpStrList), getPvalCmp )
image(1:4,1:8,  (-1*log(pmat)/log(10))%<=a=0%2, col=heat.colors(100),xaxt='n', yaxt='n',xlab='l',ylab='l' ) 
axis(1, at=1:4, labels=animalList)
axis(2, at=1:8, labels=cmpStrList)

tst <- runOneFunctionOneSubject(  1, thisDVstr='n1', FUN=subjFigure_CmpByDeltaAmp,trialType=c(1,2),show=T,eps=F )



## ==================================== get normalized component amplitude for each animal by ampInd and isiOct
doIt <- function(DVstr){ 
  resList <- list()
  for (ind in 1:length(dataList))
    resList[[ind]] <- runOneFunctionOneSubject(  ind, thisDVstr=DVstr, FUN=subjFigure_CmpByAmpAndSOA,trialType=c(1),show=T,eps=T )
  
  disc <- gaFigure_CmpByAmpAndSOA(resList, DVstr=DVstr, show=T, eps=T)
  return( list(resList=resList) )
}
resA <- lapply(cmpStrList, doIt)


#tst <- runOneFunctionOneSubject(  1, thisDVstr='n1', FUN=subjFigure_CmpByAmpAndSOA,trialType=c(1),show=F,eps=F )




