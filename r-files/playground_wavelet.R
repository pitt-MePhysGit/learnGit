baseDir <- '~/Dropbox/'

library(R.matlab)
library(rgenoud)
library(lattice)
library(akima)
library(Rwave)
library(MASS)
library(fastICA)

source(paste(baseDir, 'ampOdd_click/r-files/functions/prepDFVprobe.R',sep=''))
source( paste(baseDir,'ampOdd_click/r-files/functions/FLA_VProbe_STPSP.R',sep='') )

source( paste(baseDir,'ampOdd_click/r-files/functions/functions_ampOddclick.R',sep='') )
source( paste(baseDir,'ampOdd_click/r-files/functions/functions_AOCdual.R',sep='') )
source( paste(baseDir,'ampOdd_click/r-files/functions/functions_VProbe_STPSP.R',sep='') )


source( paste(baseDir,'ampOdd_click/r-files/functions/getCWTevk.R',sep='') )
source( paste(baseDir,'ampOdd_click/r-files/functions/getCWT.R',sep='') )
source( paste(baseDir,'ampOdd_click/r-files/functions/showCWT.R',sep='') )
source( paste(baseDir,'ampOdd_click/r-files/functions/showICx.R',sep='') )
source( paste(baseDir,'ampOdd_click/r-files/functions/getSLix.R',sep='') )
source( paste(baseDir,'ampOdd_click/r-files/functions/analyzeCorrelationStructure.R',sep=''))
        
source( paste(baseDir,'ampOdd_click/r-files/figures/showTD.R',sep='') )
source( paste(baseDir,'ampOdd_click/r-files/figures/showTD2.R',sep='') )
source( paste(baseDir,'ampOdd_click/r-files/figures/showBandsByVar.R',sep='') )

source( paste(baseDir,'r-compile/rUtils/R/rUtils.R',sep='') )
source( paste(baseDir,'r-compile/ePhysLab/R/ePhysLab.R',sep='') )

source( paste(baseDir,'RSkernel/r-files/functions/tapered_ccf.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/lamiarLFPtheory/smooth2DCSD.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/call_fastICA.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/save_fastICA.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/lamiarLFPtheory/showICAcomponents.R',sep='') )


source( paste(baseDir,'RSkernel/r-files/functions/eCogFileStr.R',sep='') )
source( paste(baseDir,'aMMN_controlled/r-files/functions_eCog.R',sep='') )

# read the curated file
colClasses <- c('character',rep('numeric',2),'character',rep('numeric',5)  ,rep('character',3),rep('numeric',9) )
df         <- read.table(paste(baseDir,'ampOdd_click/ampOddclickVProbe_link.txt',sep=''),header=T,colClasses=colClasses)

df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))

df <- prepDFVprobe(df)

## preliminaries
pp   <- function(){paste(baseDir, 'ampOdd_click/word/VProbe/plots/',sep='')}
rda  <- function(){'~/Dropbox/ampOdd_click/Rdata/'}
trda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}


##====================================================== load example ICA component
ind   <- 48
ind   <- 50

startFromScratch <- F

Nica  <- 6
pp    <- function(){paste(baseDir, 'ampOdd_click/word/VProbe/plots/',df$exp[ind],'/ICA_Npca',Nica,'/',sep='')}
if (!exists(pp()))
    system(paste('mkdir -p ', pp(), sep=''))

chInd      <- 1:6
chStr      <- paste('ic',chInd,sep='')
icv        <- loadAOC( df[ind,] , what='ICA_evk_N6', chStr=chStr, trda=trda,blcorr=T,loadCC=T  )


calcSingleTrialWavelet <- T
if ( file.exists( paste(rda(),'/',df$exp[ind],'/features.rda',sep='') ) !startFromScratch )
  calcSingleTrialWavelet <- F
             
if (!calcSingleTrialWavelet)
  load( file=paste(rda(),'/',df$exp[ind],'/features.rda',sep='') )


## ========================= load MUA data
chInd      <- df$startVIndex[ind] - 1 + 1:df$numChan[ind]
chStr      <- paste('ch',chInd,sep='')
mua        <- loadAOC( df[ind,] , what='mua', chStr=chStr, trda=trda,blcorr=T  )


## =====================================================
if (1==0){
  chInd      <- 1:12
  chStr      <- paste('ic',chInd,sep='')
  eegic      <- loadAOC( df[ind,] , what='ICA_EEG_evk_N12', chStr=chStr, trda=trda,blcorr=T,loadCC=T  )
  
  eCogFile <- paste('~/Dropbox/eCog/',eCogFileStr(df$animal[ind]),'/positionList_ISI_eCog.txt',sep='')
  colClasses <- c('numeric','character',rep('numeric',4))
  ep             <- read.table(eCogFile, header=T, colClasses=colClasses)
  eegic$ep       <- ep
  
  showICAcomponents(eegic,srtICA=T,csdWeights=T, showFFT=T, showEvokedResponse=T, includeexp=F, eps=F)
  showScalpMap(eegic$weights,rep(1:12),eval.at=1:12,eCogFile=eCogFileStr(df$animal[ind]),mirror=F,flat=T, eps=F,pdf=T)
}

## ======================================================
llp          <- icv
tmp          <- apply(llp$lfp,c(1,3),max) - apply(llp$lfp,c(1,3),min)
llp$xpp$p2p  <- apply( tmp, 1, max )


## ============================== build logarithmic time-axis
toff                                   <- 0
llp$logtaxis                           <- log2( llp$taxis + toff)
llp$logtaxis[ (llp$taxis + toff <= 0)] <- .01*(llp$taxis[ which(llp$taxis + toff <= 0) ] - 1 + toff) #- log2(taxis[min(which(taxis>-toff))] )# + .01*taxis[ (llp$taxis + toff == 0)]
mua$logtaxis <- llp$logtaxis

brks <- 0.250 * 2^(0:7)
brks[1] <- 0.2#; brks[length(brks)]<- 22 
llp$xpp$isiOct  <- cut(llp$xpp$ISI,  breaks=brks)
llp$xpp$futureISI <- my.lag(llp$xpp$ISI,lag=-1)

mua$xpp$isiOct  <- cut(mua$xpp$ISI,  breaks=brks)
mua$xpp$futureISI <- my.lag(mua$xpp$ISI,lag=-1)

llp$srt   <- sort(apply(llp$weight,1,var),decreasing=T, index.return=T)$ix

llp$var <- apply(llp$weights,1,var)
plot( sort(llp$var,decreasing=T),type='b',pch=20, ylim=c(0,max(llp$var)) )



## =====================================================
freqList <- list(D=c(2,4), T=c(4,7),   A=c(7,15),   B=c(15,40),   G=c(40,85), H=c(85,250) )
trngList <- list(t=c(8,20),e=c(20,40), i=c(40,100), r=c(100,250), s=c(250,500))

if (calcSingleTrialWavelet)
  cwbnd    <- lapply(1:length(freqList), function(cx){ getCWT(llp, 1:dim(llp$lfp)[1], chx=cx, return.raw=F, freqList=freqList )})

## ========================================================= evoked power
cwt1 <- lapply(1:6, function(cx){getCWT(llp, which(llp$xpp$ISI>.250 &llp$xpp$ISI<.500 & llp$xpp$p2p<450),chx=cx)})
cwt2 <- lapply(1:6, function(cx){getCWT(llp, which(llp$xpp$ISI>.500 &llp$xpp$ISI<1.00 & llp$xpp$p2p<450),chx=cx)})
cwt3 <- lapply(1:6, function(cx){getCWT(llp, which(llp$xpp$ISI>1.00 &llp$xpp$ISI<2.00 & llp$xpp$p2p<450),chx=cx)})
cwt4 <- lapply(1:6, function(cx){getCWT(llp, which(llp$xpp$ISI>2.00 &llp$xpp$ISI<4.00 & llp$xpp$p2p<450),chx=cx)})
cwt5 <- lapply(1:6, function(cx){getCWT(llp, which(llp$xpp$ISI>4.00 &llp$xpp$ISI<8.00 & llp$xpp$p2p<450),chx=cx)})
cwt6 <- lapply(1:6, function(cx){getCWT(llp, which(llp$xpp$ISI>8.00 &llp$xpp$ISI<16.0 & llp$xpp$p2p<450),chx=cx)})

cwa1 <- lapply(1:6, function(cx){getCWT(llp, which(llp$xpp$ampInd==1 & llp$xpp$p2p<450),chx=cx)})
cwa2 <- lapply(1:6, function(cx){getCWT(llp, which(llp$xpp$ampInd==2 & llp$xpp$p2p<450),chx=cx)})
cwa3 <- lapply(1:6, function(cx){getCWT(llp, which(llp$xpp$ampInd==3 & llp$xpp$p2p<450),chx=cx)})
cwa4 <- lapply(1:6, function(cx){getCWT(llp, which(llp$xpp$ampInd==4 & llp$xpp$p2p<450),chx=cx)})
cwa5 <- lapply(1:6, function(cx){getCWT(llp, which(llp$xpp$ampInd==5 & llp$xpp$p2p<450),chx=cx)})

cwt <- list(cwt1=cwt1, cwt2=cwt2, cwt3=cwt3, cwt4=cwt4, cwt5=cwt5, cwt6=cwt6)
cwa <- list(cwa1=cwa1, cwa2=cwa2, cwa3=cwa3, cwa4=cwa4, cwa5=cwa5)


## ============================== components by SOA and SPL
eps.file('Figure1_S_SOA_SPL.eps',pdf=T,s=2,mf=c(6,2),ratio=1)
par(mfcol=c(2,6),mar=c(4,4,1,1))
tlim <-  llp$logtaxis[ which(llp$taxis%in%c(5,250)) ]
for (cx in 1:6){
  ylim <- c(-1,1) *na.max( c(abs(cwa[[5]][[llp$srt[cx]]]$mn),abs(cwt[[6]][[llp$srt[cx]]]$mn) ) )
  
  plot(NA, xlim=tlim, ylim=ylim, xaxt='n')
  taxat <- c(5,10,15,20,30,40,50,100,200,300)
  taxnd <- which( llp$taxis%in%c(taxat) )
  axis(1, at=llp$logtaxis[taxnd], taxat  )
  abline(v=llp$logtaxis[taxnd],lty=2,lwd=.5, col='gray')
  abline(h=0,lty=1,lwd=.5,col='black')
  for (tx in 1:6)
    points( llp$logtaxis, cwt[[tx]][[llp$srt[cx]]]$mn, col=soa.colors(7)[tx], type='l'  )
  #points( cwt[[tx]][[llp$srt[cx]]]$taxis, cwt[[tx]][[llp$srt[cx]]]$mn, col=soa.colors(7)[tx], type='l'  )
  
  plot(NA, xlim=tlim, ylim=ylim, xaxt='n')
  taxnd <- which( llp$taxis%in%c(taxat) )
  axis(1, at=llp$logtaxis[taxnd], taxat  )
  abline(v=llp$logtaxis[taxnd],lty=2,lwd=.5, col='gray')
  abline(h=0,lty=1,lwd=.5,col='black')
  for (ax in 1:5)
    points( llp$logtaxis, cwa[[ax]][[llp$srt[cx]]]$mn, col=soa.colors(6)[ax], type='l'  )
}
dev.off()



#------------------------ tapered cross-correlation
eps.file('Figure2_CC(S).eps',pdf=T,s=2,mf=c(6,6),ratio=1)
par(mfrow=c(6,6))
for (cx1 in 1:6){
  for (cx2 in 1:6){
    if(cx1<= cx2)
      disc <- tapered_ccf(llp$lfp[,,c(llp$srt[cx1],llp$srt[cx2])], lag.max=250)
    if(cx1 > cx2)
      plot(NA,xlim=c(-1,1),ylim=c(-1,1))
  }
}
dev.off()


## =============================== grand average time-frequency plots
showICAcomponents(llp,srtICA=T, csdWeights=T, showFFT=T, showEvokedResponse=T, includeexp=F, eps=F)

eps.file('Figure3a_Wavelet(S).eps',pdf=T,s=2,mf=c(6,1),ratio=1)
par(mfrow=c(1,length(llp$chStr)))
for (cx in 1:length(llp$chStr))
  disc <- showICx(llp,llp$srt[cx],logtime=T, xlim=c(3,700), subs=(llp$xpp$trialType==1&llp$xpp$futureISI>.750) )
dev.off()

eps.file('Figure3b_Wavelet(S).eps',pdf=T,s=2,mf=c(6,1),ratio=1)
par(mfrow=c(1,length(llp$chStr)))
for (cx in 1:length(llp$chStr))
  disc <- showICx(llp,llp$srt[cx],logtime=F, xlim=c(-150,700), subs=(llp$xpp$trialType==1&llp$xpp$futureISI>.750) )
dev.off()


if (calcSingleTrialWavelet){
  ## =============================== bands by VAR
  showBandsByVar(cwbnd, llp, logtime=T, varStr='isiOct',subs=llp$xpp$futureISI>.750&llp$xpp$trialType==1,
                 fileName='Figure4a_Bands(S)_SOA.eps', freqList=freqList, trngList=trngList)
  
  showBandsByVar(cwbnd, llp, logtime=T, varStr='ampInd',subs=llp$xpp$futureISI>.750&llp$xpp$trialType==1,
                 fileName='Figure4a_Bands(S)_AMP.eps', freqList=freqList, trngList=trngList)
  
  
  showBandsByVar(cwbnd, llp, logtime=F, varStr='isiOct',subs=llp$xpp$futureISI>.750&llp$xpp$trialType==1,
                 fileName='Figure4b_Bands(S)_SOA.eps', freqList=freqList, trngList=trngList)
  
  showBandsByVar(cwbnd, llp, logtime=F, varStr='ampInd',subs=llp$xpp$futureISI>.750&llp$xpp$trialType==1,
                 fileName='Figure4b_Bands(S)_AMP.eps', freqList=freqList, trngList=trngList)
}

showMUAByVar(mua, logtime=T, varStr='isiOct',subs=llp$xpp$futureISI>.750&llp$xpp$trialType==1,
               fileName='Figure4c_Bands(M)_SOA.eps', freqList=freqList, trngList=trngList)

showMUAByVar(mua, logtime=T, varStr='ampInd',subs=llp$xpp$futureISI>.750&llp$xpp$trialType==1,
             fileName='Figure4c_Bands(M)_AMP.eps', freqList=freqList, trngList=trngList)



if (calcSingleTrialWavelet){
  ## ============================= extract scalars
  freqList
  trngList
  slList     <- data.frame( what='H',    trng='t', xtn='a' , stringsAsFactors=FALSE)
  slList[2,] <- data.frame( what='G',    trng='t', xtn='a' , stringsAsFactors=FALSE)
  slList[3,] <- data.frame( what='B',    trng='t', xtn='a' , stringsAsFactors=FALSE)
  
  slList[4,] <- data.frame( what='H',    trng='e',  xtn='b' , stringsAsFactors=FALSE)
  slList[5,] <- data.frame( what='G',    trng='e',  xtn='b' , stringsAsFactors=FALSE)
  slList[6,] <- data.frame( what='B',    trng='e',  xtn='b' , stringsAsFactors=FALSE)
  
  slList[7,] <- data.frame( what='G',   trng='i',  xtn='c'  , stringsAsFactors=FALSE)
  slList[8,] <- data.frame( what='B',   trng='i',  xtn='c'  , stringsAsFactors=FALSE)
  slList[9,] <- data.frame( what='A',   trng='i',  xtn='a'  , stringsAsFactors=FALSE)
  
  slList[10,] <- data.frame( what='A',  trng='r',  xtn='b'  , stringsAsFactors=FALSE)
  slList[11,] <- data.frame( what='D',  trng='r',  xtn='a'  , stringsAsFactors=FALSE)
  
  
  sL <- getSLix( what=slList$what[1], when=slList$trng[1] ,xtn=slList$xtn[1], trngList=trngList, freqList=freqList )
  for (ix in 2:dim(slList)[1])
    sL <- rbind(sL, getSLix( what=slList$what[ix], when=slList$trng[ix],xtn=slList$xtn[ix], trngList=trngList, freqList=freqList ))
  
  extractScalars <- function(sx,bnd=cwbnd,scalarList=sL1){
    if ( scalarList$what[sx]!='ts' ){
      cx    <- scalarList$icx[sx]
      bx    <- which(names(freqList)%in%scalarList$what[sx])
      taxis <- bnd[[llp$srt[cx]]]$taxis
      tndx  <- which(taxis>=scalarList$tmin[sx] & taxis<=scalarList$tmax[sx])
      res   <- apply( bnd[[llp$srt[cx]]]$bands[,tndx,bx], 1, mean)
    }
    if ( scalarList$what[sx]=='ts' ){
      cx    <- scalarList$icx[sx]
      taxis <- llp$taxis
      tndx  <- which(taxis>=scalarList$tmin[sx] & taxis<=scalarList$tmax[sx])
      res   <- apply( llp$lfp[ , tndx, llp$srt[cx] ], 1, mean)
    }
    
    return(res)
  }
  
  predMat                <- sapply(1:dim(sL)[1], extractScalars, scalarList=sL, bnd=cwbnd)
  dimnames(predMat)[[2]] <- sL$cmpStrName
  
  xpp <- cbind(llp$xpp, predMat)
  xpp$valTrial <- xpp$p2p < 450 & 1:dim(xpp)[1]>5
  
  
  ## =================== extract average MUA per trial and time-window
  extractScalars <- function(td=mua,trng=c(8,12),ext='t'){
    taxis <- td$taxis
    tndx  <- which(taxis>=trng[1] & taxis<=trng[2])
    res   <- apply( td$lfp[,tndx,], c(1,3), mean)
    dimnames(res)[[2]] <- paste('M_',ext,'_',1:dim(td$lfp)[3],sep='')
    return(res)
  }
  
  predMatMua <-                   extractScalars(td=mua, trng=c( 8, 12), ext='t')
  predMatMua <- cbind(predMatMua, extractScalars(td=mua, trng=c(20, 40), ext='e'))
  predMatMua <- cbind(predMatMua, extractScalars(td=mua, trng=c(60,100), ext='i'))
  
  xpp <- cbind(xpp, predMatMua)
  
  ## ================================ save the feature map
  fileName  <- paste(rda(),'/', df$exp[ind], '/features.rda',sep='')
  save( file=paste(rda(),'/',df$exp[ind],'/features.rda',sep=''), xpp)
}


## Compare IC weights to MUA depth profile
muaProfile <- apply( extractScalars(td=mua, trng=c(20, 40), ext='e'), 2, mean)
muaProfile <- muaProfile/max(abs(muaProfile))
showICAcomponents(icv, muaProfile=muaProfile,srtICA=T,csdWeights=T, showFFT=T, showEvokedResponse=T, includeexp=F, eps=T, pp=pp, fxtn='_muaProfile')

tmpmat <- extractScalars(td=mua, trng=c(20, 40), ext='e')
muaProfile <- tapply(1:dim(mua$xpp)[1], mua$xpp$trialNr, function(tx){apply(tmpmat[tx,], 2, mean)} )

plot(   muaProfile[[1]],-1*(1:32), type='l', col='red')
points( muaProfile[[6]],-1*(1:32), type='l', col='gray')
points( muaProfile[[12]],-1*(1:32), type='l', col='blue')

plot( apply(tmpmat,1,mean) )



## ====================================== run model fits for each predictor
source( paste(baseDir,'ampOdd_click/r-files/functions/functions_ampOdd_presynapticPlasticity.R',sep='') )
tn     <- llp
tn$xpp <- xpp

#ampFit <- callSTPSP(tn,   DVstr='bspa1', gd=T, show=F, showCh=1, who='all', BYstr='exp', 
#                    FUN=cumsumamp, minISI=.190, maxISI=12.8, startFromScratch=F, eps=F)

FUN        <- cumsumamp
ampFitList <- lapply( sL$cmpStrName, function(dvs){
                      callSTPSP(tn,   DVstr=dvs, gd=T, show=F, showCh=1, who='all', BYstr='exp', FUN=FUN, minISI=.190, maxISI=21, startFromScratch=T, eps=F)})


extractResiduals   <- function(n){ampFitList[[n]]$res$lm1$residuals}
residMat           <- sapply(1:dim(sL)[1], extractResiduals)
dimnames(residMat) <- dimnames(predMat)


extractFitted    <- function(n){ampFitList[[n]]$res$lm1$fitted}
fitMat           <- sapply(1:dim(sL)[1], extractFitted)
dimnames(fitMat) <- dimnames(predMat)


extractpVal       <- function(n){ampFitList[[n]]$res$modelTest$`Pr(>F)`[2]}
pMat              <- sapply(1:dim(predMat)[2], extractpVal)

extractpar              <- function(n){ampFitList[[n]]$modelFit$par}
parMat                  <- sapply(1:dim(predMat)[2], extractpar)
dimnames(parMat)[[1]]   <- c('tau','U','mu','sd','Off') 

dfPar            <- data.frame(feature=dimnames(predMat)[[2]], t(parMat),pMat )



## --------------------- mua
#tst1 <- callSTPSP(tn,   DVstr='mua_evk_12', gd=T, show=T, showCh=1, who='all', BYstr='exp', 
#                                        FUN=cumsumamp, minISI=.190, maxISI=12.8, startFromScratch=T, eps=F)
                    
ampFitListMua <- lapply( dimnames(predMatMua)[[2]], function(dvs){
  callSTPSP(tn,   DVstr=dvs, gd=T, show=F, showCh=1, who='all', BYstr='exp', FUN=FUN, minISI=.190, maxISI=21, startFromScratch=F, eps=T)})


extractResiduals      <- function(n){ampFitListMua[[n]]$res$lm1$residuals}
residMatMua           <- sapply(1:dim(predMatMua)[2], extractResiduals)
dimnames(residMatMua) <- dimnames(predMatMua)

extractFitted       <- function(n){ampFitListMua[[n]]$res$lm1$fitted}
fitMatMua           <- sapply(1:dim(predMatMua)[2], extractFitted)
dimnames(fitMatMua) <- dimnames(predMatMua)


extractpVal               <- function(n){ampFitListMua[[n]]$res$modelTest$`Pr(>F)`[2]}
pMat                      <- sapply(1:dim(predMatMua)[2], extractpVal)

extractpar                 <- function(n){ampFitListMua[[n]]$modelFit$par}
parMatMua                  <- sapply(1:dim(predMatMua)[2], extractpar)
dimnames(parMatMua)[[1]]   <- c('tau','U','mu','sd','Off') 
dfParMua                   <- data.frame(feature=dimnames(predMatMua)[[2]], t(parMatMua),pMat )




## ===================================================== clusters in model parameters?
dfP         <- rbind(dfPar, dfParMua)
dfP$musd    <- (((dfP$mu-2.5)/dfP$sd)%<=%4)%>=%-4
dfP$taulog  <- log2( (log(2)/dfP$tau)%<=%5 )

extractWhat <- function(nx,sx=1){ substr(dfP$feature[nx],sx,sx) }
dfP$what    <- sapply(1:dim(dfP)[[1]], extractWhat, sx=1)
dfP$time    <- sapply(1:dim(dfP)[[1]], extractWhat, sx=3)
dfP$chan    <- sapply(1:dim(dfP)[[1]], extractWhat, sx=5)

bandColors <- c(H='red',G='orange',B='green',A='blue',D='black',M='purple')
timePch    <- c(t=18,e=20,i=17,r=15)

valFeat <- which( -log10(dfP$p) > 5 )

par(mfrow=c(4,4))
dns <- density(dfP$taulog)
plot(dns$x, dns$y*length(dfP$p), type='l', xlab='taulog')
dns <- density(dfP$taulog[valFeat])
points(dns$x, dns$y*length(valFeat), type='l', col='red')
plot( dfP$taulog[valFeat], dfP$U[valFeat],    pch=timePch[dfP$time[valFeat]], col=bandColors[dfP$what[valFeat]], xlab='taulog', ylab='U');abline(v=0)
plot( dfP$taulog[valFeat], dfP$musd[valFeat], pch=timePch[dfP$time[valFeat]], col=bandColors[dfP$what[valFeat]], xlab='taulog', ylab='musd');abline(v=0)
plot( dfP$taulog[valFeat], dfP$Off[valFeat],  pch=timePch[dfP$time[valFeat]], col=bandColors[dfP$what[valFeat]], xlab='taulog', ylab='Off');abline(v=0)

plot( dfP$U[valFeat], dfP$taulog[valFeat], pch=timePch[dfP$time[valFeat]], col=bandColors[dfP$what[valFeat]],xlab='U',ylab='taulog');abline(h=0)
dns <- density(dfP$U)
plot(dns$x, dns$y*length(dfP$p), type='l', xlab='U')
dns <- density(dfP$U[valFeat])
points(dns$x, dns$y*length(valFeat), type='l', col='red')
plot( dfP$U[valFeat], dfP$musd[valFeat],    pch=timePch[dfP$time[valFeat]], col=bandColors[dfP$what[valFeat]],xlab='U',ylab='musd');abline(h=0)
plot( dfP$U[valFeat], dfP$Off[valFeat],  pch=timePch[dfP$time[valFeat]], col=bandColors[dfP$what[valFeat]],xlab='U',ylab='Off')


plot( dfP$musd[valFeat], dfP$taulog[valFeat], pch=timePch[dfP$time[valFeat]], col=bandColors[dfP$what[valFeat]], xlab='musd',ylab='taulog');abline(h=0)
plot( dfP$musd[valFeat], dfP$U[valFeat],    pch=timePch[dfP$time[valFeat]], col=bandColors[dfP$what[valFeat]],xlab='musd',ylab='U');abline(v=0)
dns <- density(dfP$musd)
plot(dns$x, dns$y*length(dfP$p), type='l', xlab='muasd')
dns <- density(dfP$musd[valFeat])
points(dns$x, dns$y*length(valFeat), type='l', col='red')
plot( dfP$musd[valFeat], dfP$Off[valFeat],  pch=timePch[dfP$time[valFeat]], col=bandColors[dfP$what[valFeat]],xlab='musd',ylab='Off');abline(v=0)


plot( dfP$Off[valFeat], dfP$taulog[valFeat], pch=timePch[dfP$time[valFeat]], col=bandColors[dfP$what[valFeat]], xlab='Off',ylab='taulog');abline(h=0)
plot( dfP$Off[valFeat], dfP$U[valFeat],    pch=timePch[dfP$time[valFeat]], col=bandColors[dfP$what[valFeat]],xlab='Off',ylab='U')
plot( dfP$Off[valFeat], dfP$musd[valFeat],  pch=timePch[dfP$time[valFeat]], col=bandColors[dfP$what[valFeat]],xlab='Off',ylab='musd');abline(h=0)
dns <- density(dfP$Off)
plot(dns$x, dns$y*length(dfP$p), type='l', xlab='muasd')
dns <- density(dfP$Off[valFeat])
points(dns$x, dns$y*length(valFeat), type='l', col='red')


pcatst <- pca( dfP[valFeat,c('taulog','U','musd','Off','sd')], center=T, scale=T, ncomp=2  )
plot(pcatst$explained_variance )

plot(pcatst$x[,1],pcatst$x[,2], pch=timePch[dfP$time[valFeat]], col=bandColors[dfP$what[valFeat]] )


## ======================================= plot dependence on SOA and AMP for the scalars
eps.file('Figure5a_Scalars_SOA_AMP_Fit.eps',pdf=T,s=2,mf=c(6,dim(sL)[1]%/%6),ratio=1)
par( mfrow=c(dim(sL)[1]%/%6,6), mar=c(4,4,2,1) )
disc <- lapply(1:dim(sL)[1], showTD,  pM=predMat, fM=fitMat, xpp=xpp)
dev.off()

eps.file('Figure5a_Scalars_SOA_AMP.eps',pdf=T,s=2,mf=c(6,dim(sL)[1]%/%6),ratio=1)
par( mfrow=c(dim(sL)[1]%/%6,6), mar=c(4,4,2,1) )
disc <- lapply(1:dim(sL)[1], showTD,  pM=predMat, fM=NULL, xpp=xpp)
dev.off()

eps.file('Figure5a_Scalars2a_SOA_AMP.eps',pdf=T,s=2,mf=c(6,dim(sL)[1]%/%6),ratio=1)
par( mfrow=c(dim(sL)[1]%/%6,6), mar=c(4,4,2,1) )
disc <- lapply(1:dim(sL)[1], showTD2,  pM=predMat, xpp=xpp, fM=NULL, IVstr1='isiOct', IVstr2='ampInd', colFUN=soa.colors)
dev.off()

eps.file('Figure5a_Scalars2b_AMP_SOA.eps',pdf=T,s=2,mf=c(6,dim(sL)[1]%/%6),ratio=1)
par( mfrow=c(dim(sL)[1]%/%6,6), mar=c(4,4,2,1) )
disc <- lapply(1:dim(sL)[1], showTD2,  pM=predMat, xpp=xpp, fM=NULL, IVstr1='ampInd', IVstr2='isiOct', colFUN=bluered.colors)
dev.off()



eps.file('Figure5b_MUA_SOA_AMP_Fit.eps',pdf=T,s=2,mf=c(8,dim(predMatMua)[2]%/%8),ratio=1)
par( mfrow=c(dim(predMatMua)[2]%/%8,8), mar=c(4,4,2,1) )
disc <- lapply(1:dim(predMatMua)[2], showTD,  pM=predMatMua, fM=fitMatMua, xpp=xpp)
dev.off()

eps.file('Figure5b_MUA_SOA_AMP.eps',pdf=T,s=2,mf=c(8,dim(predMatMua)[2]%/%8),ratio=1)
par( mfrow=c(dim(predMatMua)[2]%/%8,8), mar=c(4,4,2,1) )
disc <- lapply(1:dim(predMatMua)[2], showTD,  pM=predMatMua, fM=NULL,xpp=xpp)
dev.off()

eps.file('Figure5b_MUA2a_SOA_AMP.eps',pdf=T,s=2,mf=c(6,dim(sL)[1]%/%6),ratio=1)
par( mfrow=c(dim(predMatMua)[2]%/%8,8), mar=c(4,4,2,1) )
disc <- lapply(1:dim(predMatMua)[2], showTD2,  pM=predMatMua, xpp=xpp, fM=NULL, IVstr1='isiOct', IVstr2='ampInd', colFUN=soa.colors)
dev.off()

eps.file('Figure5b_MUA2b_AMP_SOA.eps',pdf=T,s=2,mf=c(6,dim(sL)[1]%/%6),ratio=1)
par( mfrow=c(dim(predMatMua)[2]%/%8,8), mar=c(4,4,2,1) )
disc <- lapply(1:dim(predMatMua)[2], showTD2,  pM=predMatMua, xpp=xpp, fM=NULL, IVstr1='ampInd', IVstr2='isiOct', colFUN=bluered.colors)
dev.off()


## ================================================================ raw correlation matrix
## ============================================================================== noise correlation
disc <- analyzeCorrelationStructure(predMat   , predMat,    fileName='Figure6a_Cor(CSD)_raw.eps'  )
disc <- analyzeCorrelationStructure(predMat   , predMatMua, fileName='Figure6a_Cor(CSD,MUA)_raw.eps' )
disc <- analyzeCorrelationStructure(predMatMua, predMatMua, fileName='Figure6a_Cor(MUA)_raw.eps' )

disc <- analyzeCorrelationStructure(residMat   , residMat,   fileName='Figure6b_Cor(CSD)_noise.eps')
CC_SM_N <- analyzeCorrelationStructure(residMat   , residMatMua,fileName='Figure6b_Cor(CSD,MUA)_noise.eps', alpha=.001)
disc <- analyzeCorrelationStructure(residMatMua, residMatMua,fileName='Figure6b_Cor(MUA)_noise.eps')

disc <- analyzeCorrelationStructure(fitMat, fitMat,fileName='Figure6c_Cor(CSD)_fitted.eps',zlimCor=1)
disc <- analyzeCorrelationStructure(fitMat, fitMatMua,fileName='Figure6c_Cor(CSD,MUA)_fitted.eps',zlimCor=1)
disc <- analyzeCorrelationStructure(fitMatMua, fitMatMua,fileName='Figure6c_Cor(MUA)_fitted.eps',zlimCor=1)


## ========================================== SVM-based decoding
library(e1071)
tx          <- xpp[, names(xpp) %in% c(dimnames(predMat)[[2]], dimnames(predMatMua)[[2]]) ]

eps.file('Figure7_Decoding_all.eps',pdf=T,s=2,mf=c(2,1),ratio=1)
par( mfrow=c(1,2), mar=c(4,4,2,1) )

tx$category <- as.factor(xpp$ampInd)
decodeSVM(tx,subs=xpp$ISI<16)

tx$category <- as.factor(xpp$isiOct)
decodeSVM(tx,subs=xpp$ISI<16)
dev.off()




tset        <- sample( which(tx$valTrial), dim(tx)[1]/2 )
trainingSet <- tx[tset,c(43:206)]

tmpset  <- which(tx$valTrial)
eset    <- setdiff(tmpset,tset)
evalSet <- tx[eset,c(43:206)]


svmfit = svm(category ~ ., data = trainingSet, kernel = "sigmoid", cost = 1, scale = TRUE, cross=10)#, gamma=.0015)
print(svmfit)
trainingSet$prediction <- predict(svmfit, trainingSet)
trainingSet$correct    <- trainingSet$category == trainingSet$prediction
mean(trainingSet$correct)

evalSet$prediction <- predict(svmfit,evalSet)
evalSet$correct    <- evalSet$category == evalSet$prediction
mean(evalSet$correct)

table(category=evalSet$category, prediction=evalSet$prediction)
getCnfCat <- function(cx){ table(category=evalSet$category, prediction=evalSet$prediction)[cx,]/table(evalSet$category)[cx] }
confusionMat <- sapply(seq(length(levels(tx$category)),by=-1,to=1), getCnfCat)
confusionMat
image(confusionMat, col=gray.colors(100), zlim=c(0,1))



## ===========================================  try Chord plots
library(circlize)


thrsh <- .08



## ================= functional plot
mat <- abs(CC_SM_N$corSqMat) 
mat[ which(CC_SM_N$corPmatbf == 1 )] <- 0
mat <- mat[seq(dim(mat)[1],by=-1,to=1) ,]

extractWhat <- function(nx,d=1,sx=1){ substr(dimnames(mat)[[d]][nx],sx,sx) }
whatList1 <- sapply(1:dim(mat)[[1]], extractWhat, d=1, sx=1)
timeList1 <- sapply(1:dim(mat)[[1]], extractWhat, d=1, sx=3)
chanList1 <- sapply(1:dim(mat)[[1]], extractWhat, d=1, sx=5)

whatList2 <- sapply(1:dim(mat)[[2]], extractWhat, d=2, sx=1)
timeList2 <- sapply(1:dim(mat)[[2]], extractWhat, d=2, sx=3)
chanList2 <- sapply(1:dim(mat)[[2]], extractWhat, d=2, sx=5)

bandColors <- c(H='red',G='orange',B='green',A='blue',D='black')
timeColors <- c(t='red',e='orange',i='green',s='blue')


grid.colors = c( bandColors[whatList1], timeColors[timeList2] )
names(grid.colors) <- c( dimnames(mat)[[1]], dimnames(mat)[[2]])

chordDiagram(mat, annotationTrack = "grid", grid.col=grid.colors, big.gap=20,
             preAllocateTracks = list(list(track.height = max(strwidth(unlist(dimnames(mat))))),
                                      list(track.height=.1) ),
                                      list(track.height=.1)
             )

circos.track(track.index=1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important


circos.clear()




chordDiagram(mat, annotationTrack = "grid", preAllocateTracks=3)




chordDiagram( tst )
circos.clear()



tmp <- disc[1:10,1:10]
chordDiagram( tmp, annotationTrack='grid',grid.border='black',
              preAllocateTracks = list(
                list(track.height = 0.02)
              ) 
)
circos.clear




stop()


pmICA <- fastICA(predMat, 8, method='C', row.norm=T)
image(pmICA$K, col=bluered.colors(100),zlim=c(-1,1))

par( mfrow=c(1,dim(pmICA$S)[2]), mar=c(2,2,2,1) )
disc <- lapply(1:dim(pmICA$S)[2], showTD,  pM=pmICA$S)


what <- 'cwe'; # cwe or cw
cx   <- 3
flab <- c(1,2,3,10,20,30,50,100,200)
xlim <- c(-50,250)
ylim <- c(1,300)
zlim <- c(0,.4)
yoff <- 10

par(mfrow=c(2,5))
showCWT(cwt1[[llp$srt[cx]]], xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab,what=what)
showCWT(cwt2[[llp$srt[cx]]], xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab,what=what)
showCWT(cwt3[[llp$srt[cx]]], xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab,what=what)
showCWT(cwt4[[llp$srt[cx]]], xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab,what=what)
showCWT(cwt5[[llp$srt[cx]]], xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab,what=what)

showCWT(cwa1[[llp$srt[cx]]], xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab,what=what)
showCWT(cwa2[[llp$srt[cx]]], xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab,what=what)
showCWT(cwa3[[llp$srt[cx]]], xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab,what=what)
showCWT(cwa4[[llp$srt[cx]]], xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab,what=what)
showCWT(cwa5[[llp$srt[cx]]], xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab,what=what)








flab <- c(1,2,3,10,20,30,50,100,200)
xlim <- c(-150,750)
ylim <- c(1,300)
zlim <- c(0,.3)
yoff <- 10

flab <- c(10,20,30,50,100,200)
xlim <- c(-10,50)
ylim <- c(50,300)
zlim <- c(0,.2)
yoff <- 10





## ========================================== pull out time-series

addLine <- function(cw,frqrng=c(8,10),...){
  fndx <- which(cw$faxis[length(cw$faxis):-1:1]>=frqrng[1]&cw$faxis[length(cw$faxis):-1:1]<=frqrng[2] )
  mn   <- apply(Mod(cw$cw[,fndx]),1,mean)
  points(cw$taxis,mn,...)
}

xlim   <- c(-10,70)
ylim   <- c(0,.2)
frqrng <-c(150,250)

par(mfcol=c(1,1))
plot(NA, xlim=xlim, ylim=ylim)
addLine(cwt5, frqrng=frqrng,type='l')
addLine(cw5,  frqrng=frqrng,type='l', col='red')


xlim   <- c(-150,750)
ylim   <- c(0,.3)
frqrng <-c(7,15)

par(mfcol=c(1,1))
plot(NA, xlim=xlim, ylim=ylim)
addLine(cw1, frqrng=frqrng,type='l',col='black')
addLine(cw2, frqrng=frqrng,type='l',col='blue')
addLine(cw3, frqrng=frqrng,type='l',col='green')
addLine(cw4, frqrng=frqrng,type='l',col='orange')
addLine(cw5, frqrng=frqrng,type='l',col='red')



xlim   <- c(-10,70)
ylim   <- c(0,.25)
frqrng <-c(150,250)

xlim   <- c(-150,350)
frqrng <-c(7,18)

xlim   <- c(-150,750)
frqrng <-c(2,4)


par(mfrow=c(2,2))
plot(NA, xlim=xlim, ylim=ylim)
addLine(cw1, frqrng=frqrng,type='l',col='black')
addLine(cw2, frqrng=frqrng,type='l',col='blue')
addLine(cw3, frqrng=frqrng,type='l',col='green')
addLine(cw4, frqrng=frqrng,type='l',col='orange')
addLine(cw5, frqrng=frqrng,type='l',col='red')

plot(NA, xlim=xlim, ylim=ylim)
addLine(cwa1, frqrng=frqrng,type='l',col='black')
addLine(cwa2, frqrng=frqrng,type='l',col='blue')
addLine(cwa3, frqrng=frqrng,type='l',col='green')
addLine(cwa4, frqrng=frqrng,type='l',col='orange')
addLine(cwa5, frqrng=frqrng,type='l',col='red')

plot(NA, xlim=xlim, ylim=ylim)
addLine(cwt1, frqrng=frqrng,type='l',col='black')
addLine(cwt2, frqrng=frqrng,type='l',col='blue')
addLine(cwt3, frqrng=frqrng,type='l',col='green')
addLine(cwt4, frqrng=frqrng,type='l',col='orange')
addLine(cwt5, frqrng=frqrng,type='l',col='red')

plot(NA, xlim=xlim, ylim=ylim)
addLine(cwta1, frqrng=frqrng,type='l',col='black')
addLine(cwta2, frqrng=frqrng,type='l',col='blue')
addLine(cwta3, frqrng=frqrng,type='l',col='green')
addLine(cwta4, frqrng=frqrng,type='l',col='orange')
addLine(cwta5, frqrng=frqrng,type='l',col='red')








## two frequencies in one plot
xlim   <- c(-50,250)
frqrng1 <-c(7,18)
frqrng2 <-c(150,250)


plot(NA, xlim=xlim, ylim=ylim)
addLine(cw1, frqrng=frqrng1,type='l',col='black')
addLine(cw2, frqrng=frqrng1,type='l',col='blue')
addLine(cw3, frqrng=frqrng1,type='l',col='green')
addLine(cw4, frqrng=frqrng1,type='l',col='orange')
addLine(cw5, frqrng=frqrng1,type='l',col='red')

addLine(cw1, frqrng=frqrng2,type='l',col='black',lty=2)
addLine(cw2, frqrng=frqrng2,type='l',col='blue',lty=2)
addLine(cw3, frqrng=frqrng2,type='l',col='green',lty=2)
addLine(cw4, frqrng=frqrng2,type='l',col='orange',lty=2)
addLine(cw5, frqrng=frqrng2,type='l',col='red',lty=2)





## ========================================================== 
OneOverFcorrection <- array( rep(faxis[length(faxis):1],each=length(taxis)), dim(wmn1)   ) 
nwmn1 <- Mod(wmn1) * OneOverFcorrection
nwmn2 <- Mod(wmn2) * OneOverFcorrection
nwmn3 <- Mod(wmn3) * OneOverFcorrection

mxPow <- max(c(Mod(nwmn1[,30:108]), Mod(nwmn2[,30:108]), Mod(nwmn3[,30:108])) )
mxPow <- max(c(Mod(nwmn3)))

flab <- c(5,seq(10,by=10,to=70))
yConv <- function(x){ log(2*pi*x/(w0*nF))/log(2)   }
ylim <- yConv(c(3,100) )

par(mfrow=c(1,3))
image(taxis, faxlog , nwmn1[,dim(wmn1)[2]:1],xlim=c(-50,120),ylim=ylim,zlim=c(0,mxPow),col=jet.colors(256),useRaster=T,yaxt='n' )
axis(2,at=yConv(flab),  label=flab)
abline(h=yConv(c(10,20,30,40,50,60,70)), col='white',lty=2 )
abline(v=0,col='white')
#abline(v=c(8,23,35,40,55),col='white',lty=2)

image(taxis, faxlog , nwmn2[,dim(wmn2)[2]:1],xlim=c(-50,120),ylim=ylim,zlim=c(0,mxPow),col=jet.colors(256),useRaster=T,yaxt='n' )
axis(2,at=yConv(flab),  label=flab)
abline(h=yConv(c(10,20,30,40,50,60,70)), col='white',lty=2 )
abline(v=0,col='white')
#abline(v=c(8,23,35,40,55),col='white',lty=2)

image(taxis, faxlog , nwmn3[,dim(wmn3)[2]:1],xlim=c(-50,120),ylim=ylim,zlim=c(0,mxPow),col=jet.colors(256),useRaster=T,yaxt='n' )
axis(2,at=yConv(flab),  label=flab)
abline(h=yConv(c(10,20,30,40,50,60,70)), col='white',lty=2 )
abline(v=0,col='white')
#abline(v=c(8,23,35,40,55),col='white',lty=2)


## =======================================================================
##====================================================== Induced activity
noctave <- 7
nvoice  <- 6
sr      <- 500 # [1/sec]
nF      <- sr/2
w0      <- 2*pi

faxis  <- w0/(2*pi) * nF * 2^(seq( to=0, by=1/nvoice, length=nvoice*noctave ))
faxlog <- seq( to=0, by=1/nvoice, length=nvoice*noctave )

flab <- c(5,seq(10,by=10,to=70))
yConv <- function(x){ log(2*pi*x/(w0*nF))/log(2)   }
ylim <- yConv(c(3,100) )


getCWT2   <- function(n){ cwt(llp$lfp[n,,], noctave=noctave,  nvoice=nvoice, w0=w0, twoD=TRUE, plot=FALSE)}
tst      <- getCWT2(1)
getCWT <- function(ind){
  tmp      <- sapply(ind, getCWT2 )
  dim(tmp) <- c(dim(tst),length(ind))
  #res      <- Mod(apply( tmp/Mod(tmp), c(1,2), mean))
  res      <-     apply( Mod(tmp),     c(1,2), mean)
  return(res)
}
ind1 <- which(llp$xpp$ISI>.4 &llp$xpp$ISI<.8 & llp$xpp$trialType==1 & llp$xpp$p2p<450)
ind2 <- which(llp$xpp$ISI>.8 & llp$xpp$ISI<3.2 & llp$xpp$trialType==1& llp$xpp$p2p<450)
ind3 <- which(llp$xpp$ISI>3.2 & llp$xpp$ISI<6.4 & llp$xpp$trialType==1& llp$xpp$p2p<450)

wmn1 <- getCWT( sample( ind1, 500, replace=F ) )
wmn2 <- getCWT( sample( ind2, 500, replace=F ) )
wmn3 <- getCWT( ind3 )

mxPow    <- max(c(wmn1,wmn2,wmn3))
xlim <- c(-150,450)
par(mfrow=c(1,3))
image(taxis, faxlog , wmn1[,dim(wmn1)[2]:1],xlim=xlim,ylim=ylim,zlim=c(0,mxPow),col=jet.colors(256),useRaster=T,yaxt='n' )
axis(2,at=yConv(flab),  label=flab)
abline(h=yConv(c(10,20,30,40,50,60,70)), col='white',lty=2 )
abline(v=0,col='white')

image(taxis, faxlog , wmn2[,dim(wmn1)[2]:1],xlim=xlim,ylim=ylim,zlim=c(0,mxPow),col=jet.colors(256),useRaster=T,yaxt='n' )
axis(2,at=yConv(flab),  label=flab)
abline(h=yConv(c(10,20,30,40,50,60,70)), col='white',lty=2 )
abline(v=0,col='white')

image(taxis, faxlog , wmn3[,dim(wmn1)[2]:1],xlim=xlim,ylim=ylim,zlim=c(0,mxPow),col=jet.colors(256),useRaster=T,yaxt='n' )
axis(2,at=yConv(flab),  label=flab)
abline(h=yConv(c(10,20,30,40,50,60,70)), col='white',lty=2 )
abline(v=0,col='white')





## =====================================================
# -------------
noctave <- 10
nvoice  <- 5

sr      <- 1000 # [1/sec]
nF      <- sr/2
w0      <- 1*pi

faxis  <- w0/(2*pi) * nF * 2^(seq( to=0, by=1/nvoice, length=nvoice*noctave ))

chirp <- sin(2*pi * 10*taxis/1000) + sin(2*pi * 100*taxis/1000)
retChirp <- cwt(chirp, noctave=noctave, nvoice=nvoice, w0=w0, plot=FALSE)
image(taxis, faxis, Mod(retChirp[,dim(retChirp)[2]:1]),col=jet.colors(256),useRaster=F )
abline(h=10)


#
plot( Re(morlet(500, 150, 50, w0=2*pi)), type='l' )
plot( Mod(gabor(500, 250, 100,4)), type='l' )








## ===============================================================================================
## ==================================================== OLD
## ---------------------- SOA
chx <- 5
mxp2p <- 450
cw1 <- getCWTevk(llp, which(llp$xpp$ISI>.250 &llp$xpp$ISI<.500 & llp$xpp$p2p<mxp2p), chx=chx )
cw2 <- getCWTevk(llp, which(llp$xpp$ISI>.500 &llp$xpp$ISI<1.00 & llp$xpp$p2p<mxp2p), chx=chx )
cw3 <- getCWTevk(llp, which(llp$xpp$ISI>1.00 &llp$xpp$ISI<2.00 & llp$xpp$p2p<mxp2p), chx=chx )
cw4 <- getCWTevk(llp, which(llp$xpp$ISI>2.00 &llp$xpp$ISI<4.00 & llp$xpp$p2p<mxp2p), chx=chx )
cw5 <- getCWTevk(llp, which(llp$xpp$ISI>4.00 &llp$xpp$ISI<8.00 & llp$xpp$p2p<mxp2p), chx=chx )

## --------------------- AMP
cwa1 <- getCWTevk(llp, which(llp$xpp$ampInd==1 & llp$xpp$p2p<mxp2p), chx=chx )
cwa2 <- getCWTevk(llp, which(llp$xpp$ampInd==2 & llp$xpp$p2p<mxp2p), chx=chx )
cwa3 <- getCWTevk(llp, which(llp$xpp$ampInd==3 & llp$xpp$p2p<mxp2p), chx=chx )
cwa4 <- getCWTevk(llp, which(llp$xpp$ampInd==4 & llp$xpp$p2p<mxp2p), chx=chx )
cwa5 <- getCWTevk(llp, which(llp$xpp$ampInd==5 & llp$xpp$p2p<mxp2p), chx=chx )


flab <- c(10,20,30,50,100,200)
xlim <- c(-50,350)
ylim <- c(3,300)
zlim <- c(0,.9)
yoff <- 10

flab <- c(10,20,30,50,100,200)
xlim <- c(-10,50)
ylim <- c(50,300)
zlim <- c(0,.2)
yoff <- 10

par(mfrow=c(2,5))
showCWT(cw1, xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab)
showCWT(cw2, xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab)
showCWT(cw3, xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab)
showCWT(cw4, xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab)
showCWT(cw5, xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab)

showCWT(cwa1, xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab)
showCWT(cwa2, xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab)
showCWT(cwa3, xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab)
showCWT(cwa4, xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab)
showCWT(cwa5, xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab)


