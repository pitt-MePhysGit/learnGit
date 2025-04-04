function_wavelet <- function(ind){
  
  baseDir <- '~/Dropbox/'
  library(R.matlab)
  library(rgenoud)
  library(lattice)
  library(akima)
  library(Rwave)
  library(MASS)
  library(fastICA)
  library(e1071)
  library(circlize)
  
  
  freqList <- list(D=c(2,4), T=c(4,7),   A=c(7,15),   B=c(15,40),   G=c(40,85), H=c(85,250) )
  trngList <- list(t=c(8,20),e=c(20,40), i=c(40,100), r=c(100,250), s=c(250,500))
  
  
  source(paste(baseDir, 'ampOdd_click/r-files/functions/prepDFVprobe.R',sep=''))
  source( paste(baseDir,'ampOdd_click/r-files/functions/FLA_VProbe_STPSP.R',sep='') )
  
  source( paste(baseDir,'ampOdd_click/r-files/functions/functions_ampOddclick.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/functions/functions_AOCdual.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/functions/functions_VProbe_STPSP.R',sep='') )
  
  
  
  source( paste(baseDir,'r-compile/rUtils/R/rUtils.R',sep='') )
  source( paste(baseDir,'r-compile/ePhysLab/R/ePhysLab.R',sep='') )
  
  source( paste(baseDir,'RSkernel/r-files/functions/tapered_ccf.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/lamiarLFPtheory/smooth2DCSD.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/call_fastICA.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/save_fastICA.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/lamiarLFPtheory/showICAcomponents.R',sep='') )
  
  
  source( paste(baseDir,'RSkernel/r-files/functions/eCogFileStr.R',sep='') )
  source( paste(baseDir,'aMMN_controlled/r-files/functions_eCog.R',sep='') )
  
  
  source( paste(baseDir,'ampOdd_click/r-files/functions/getCWTevk.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/functions/getCWT.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/functions/showCWT.R',sep='') )
  
  source( paste(baseDir,'ampOdd_click/r-files/functions/showICx.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/functions/getSLix.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/functions/analyzeCorrelationStructure.R',sep=''))
  source( paste(baseDir,'ampOdd_click/r-files/functions/decodeSVM.R',sep=''))
  
  source( paste(baseDir,'ampOdd_click/r-files/figures/showTD.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/figures/showTD2.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/figures/showBandsByVar.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/figures/showMUAByVar.R',sep='') )
  
  
  
  # read the curated file
  colClasses <- c('character',rep('numeric',2),'character',rep('numeric',5)  ,rep('character',3),rep('numeric',9) )
  df         <- read.table(paste(baseDir,'ampOdd_click/ampOddclickVProbe_link.txt',sep=''),header=T,colClasses=colClasses)
  df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))
  df         <- prepDFVprobe(df)
  
  ## preliminaries
  #pp   <- function(){paste(baseDir, 'ampOdd_click/word/VProbe/plots/',sep='')}
  rda  <- function(){'~/Dropbox/ampOdd_click/Rdata/'}
  trda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
  .GlobalEnv$rda <- rda
  .GlobalEnv$trda <- trda
  
  
  ##====================================================== load example ICA component
  #ind   <- 48
  #ind   <- 52
  
  Nica  <- 6
  .GlobalEnv$pp    <- function(){paste(baseDir, 'ampOdd_click/word/VProbe/plots/',df$exp[ind],'/ICA_Npca',Nica,'/',sep='')}
  pp <- .GlobalEnv$pp
  if (!exists(pp()))
    system(paste('mkdir -p ', pp(), sep=''))
  
  if (!exists(paste( pp(), '/STPSP/', sep='')))
    system(paste('mkdir -p ', pp(),'/STPSP/', sep=''))
  
  chInd      <- 1:6
  chStr      <- paste('ic',chInd,sep='')
  icv        <- loadAOC( df[ind,] , what='ICA_evk_N6', chStr=chStr, trda=trda,blcorr=T,loadCC=T  )
  
  ## ========================= load MUA data
  chInd      <- df$startVIndex[ind] - 1 + 1:df$numChan[ind]
  chStr      <- paste('ch',chInd,sep='')
  mua        <- loadAOC( df[ind,] , what='mua', chStr=chStr, trda=trda,blcorr=T  )
  
  ## ======================================================
  llp          <- icv
  tmp          <- apply(llp$lfp,c(1,3),max) - apply(llp$lfp,c(1,3),min)
  llp$xpp$p2p  <- apply( tmp, 1, max )
  
  ## ============================== build logarithmic time-axis
  toff                                   <- 0
  llp$logtaxis                           <- log2( llp$taxis + toff)
  llp$logtaxis[ (llp$taxis + toff <= 1)] <- .01*(llp$taxis[ which(llp$taxis + toff <= 1) ] - 1 + toff) #- log2(taxis[min(which(taxis>-toff))] )# + .01*taxis[ (llp$taxis + toff == 0)]
  
  mua$logtaxis <- llp$logtaxis
  
  
  brks <- 0.250 * 2^(0:7)
  brks[1] <- 0.2#; brks[length(brks)]<- 22 
  llp$xpp$isiOct  <- cut(llp$xpp$ISI,  breaks=brks)
  llp$xpp$futureISI <- my.lag(llp$xpp$ISI,lag=-1)
  
  mua$xpp$isiOct  <- cut(mua$xpp$ISI,  breaks=brks)
  mua$xpp$futureISI <- my.lag(mua$xpp$ISI,lag=-1)
  
  
  llp$srt   <- sort(apply(llp$weight,1,var),decreasing=T, index.return=T)$ix
  
  llp$var <- apply(llp$weights,1,var)

  ## =====================================================
  cwbnd    <- lapply(1:length(freqList), function(cx){ getCWT(llp, 1:dim(llp$lfp)[1], chx=cx, return.raw=F, freqList=freqList )})
  
  
  ## ========================================================= evoked power
  cwt1 <- lapply(1:Nica, function(cx){getCWT(llp, which(llp$xpp$ISI>.250 &llp$xpp$ISI<.500 & llp$xpp$p2p<450),chx=cx)})
  cwt2 <- lapply(1:Nica, function(cx){getCWT(llp, which(llp$xpp$ISI>.500 &llp$xpp$ISI<1.00 & llp$xpp$p2p<450),chx=cx)})
  cwt3 <- lapply(1:Nica, function(cx){getCWT(llp, which(llp$xpp$ISI>1.00 &llp$xpp$ISI<2.00 & llp$xpp$p2p<450),chx=cx)})
  cwt4 <- lapply(1:Nica, function(cx){getCWT(llp, which(llp$xpp$ISI>2.00 &llp$xpp$ISI<4.00 & llp$xpp$p2p<450),chx=cx)})
  cwt5 <- lapply(1:Nica, function(cx){getCWT(llp, which(llp$xpp$ISI>4.00 &llp$xpp$ISI<8.00 & llp$xpp$p2p<450),chx=cx)})
  cwt6 <- lapply(1:Nica, function(cx){getCWT(llp, which(llp$xpp$ISI>8.00 &llp$xpp$ISI<16.0 & llp$xpp$p2p<450),chx=cx)})
  
  cwa1 <- lapply(1:Nica, function(cx){getCWT(llp, which(llp$xpp$ampInd==1 & llp$xpp$p2p<450),chx=cx)})
  cwa2 <- lapply(1:Nica, function(cx){getCWT(llp, which(llp$xpp$ampInd==2 & llp$xpp$p2p<450),chx=cx)})
  cwa3 <- lapply(1:Nica, function(cx){getCWT(llp, which(llp$xpp$ampInd==3 & llp$xpp$p2p<450),chx=cx)})
  cwa4 <- lapply(1:Nica, function(cx){getCWT(llp, which(llp$xpp$ampInd==4 & llp$xpp$p2p<450),chx=cx)})
  cwa5 <- lapply(1:Nica, function(cx){getCWT(llp, which(llp$xpp$ampInd==5 & llp$xpp$p2p<450),chx=cx)})
  
  cwt <- list(cwt1=cwt1, cwt2=cwt2, cwt3=cwt3, cwt4=cwt4, cwt5=cwt5, cwt6=cwt6)
  cwa <- list(cwa1=cwa1, cwa2=cwa2, cwa3=cwa3, cwa4=cwa4, cwa5=cwa5)
  
  llp$xpp$blockNr <- 1+(llp$xpp$trialNr-1)%/%2
  Nblock          <- length(unique(llp$xpp$blockNr))
  
  cwbf <- function(bx){ lapply(1:Nica, function(cx){getCWT(llp, which(llp$xpp$blockNr==bx & llp$xpp$p2p<450),chx=cx)}) }
  cwb  <- lapply(1:Nblock, cwbf)
  
  ## ============================== components by SOA and SPL and BLK
  eps.file('Figure1_S_SOA_SPL.eps',pdf=T,s=1,mf=c(6,3),ratio=1)
  par(mfcol=c(3,6),mar=c(4,4,1,1))
  lwd <- .7
  tlim <-  llp$logtaxis[ which(llp$taxis%in%c(5,450)) ]
  for (cx in 1:6){
    ylim <- c(-1,1) *na.max( c(abs(cwa[[5]][[llp$srt[cx]]]$mn),abs(cwt[[6]][[llp$srt[cx]]]$mn) ) )
    
    plot(NA, xlim=tlim, ylim=ylim, xaxt='n')
    taxat <- c(5,10,15,20,30,40,50,100,200,300)
    taxnd <- which( llp$taxis%in%c(taxat) )
    axis(1, at=llp$logtaxis[taxnd], taxat  )
    abline(v=llp$logtaxis[taxnd],lty=2,lwd=.5, col='gray')
    abline(h=0,lty=1,lwd=.5,col='black')
    for (tx in 1:6)
      points( llp$logtaxis, cwt[[tx]][[llp$srt[cx]]]$mn, col=soa.colors(7)[tx], type='l',lwd=lwd  )
    #points( cwt[[tx]][[llp$srt[cx]]]$taxis, cwt[[tx]][[llp$srt[cx]]]$mn, col=soa.colors(7)[tx], type='l'  )
    
    plot(NA, xlim=tlim, ylim=ylim, xaxt='n')
    taxnd <- which( llp$taxis%in%c(taxat) )
    axis(1, at=llp$logtaxis[taxnd], taxat  )
    abline(v=llp$logtaxis[taxnd],lty=2,lwd=.5, col='gray')
    abline(h=0,lty=1,lwd=.5,col='black')
    for (ax in 1:5)
      points( llp$logtaxis, cwa[[ax]][[llp$srt[cx]]]$mn, col=bluered.colors(5)[ax], type='l',lwd=lwd  )
    
    if(Nblock>1){
      plot(NA, xlim=tlim, ylim=ylim, xaxt='n')
      taxnd <- which( llp$taxis%in%c(taxat) )
      axis(1, at=llp$logtaxis[taxnd], taxat  )
      abline(v=llp$logtaxis[taxnd],lty=2,lwd=.5, col='gray')
      abline(h=0,lty=1,lwd=.5,col='black')
      for (bx in 1:Nblock)
        points( llp$logtaxis, cwb[[bx]][[llp$srt[cx]]]$mn, col=bluered.colors(Nblock)[bx], type='l',lwd=lwd  )
      
    }
  }
  dev.off()
  
  
  
  #------------------------ tapered cross-correlation
  eps.file('Figure2_CC(S).eps',pdf=T,s=1.25,mf=c(6,6),ratio=1)
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
  eps.file('Figure3a_Wavelet(S).eps',pdf=T,s=2,mf=c(6,1),ratio=1)
  par(mfrow=c(1,length(llp$chStr)))
  for (cx in 1:length(llp$chStr))
    disc <- showICx(llp,llp$srt[cx],logtime=T, xlim=c(3,700), subs=(llp$xpp$trialType==1&llp$xpp$futureISI>.750), freqList=freqList, trngList=trngList )
  dev.off()
  
  eps.file('Figure3b_Wavelet(S).eps',pdf=T,s=2,mf=c(6,1),ratio=1)
  par(mfrow=c(1,length(llp$chStr)))
  for (cx in 1:length(llp$chStr))
    disc <- showICx(llp,llp$srt[cx],logtime=F, xlim=c(-150,700), subs=(llp$xpp$trialType==1&llp$xpp$futureISI>.750), freqList=freqList, trngList=trngList )
  dev.off()
  
  
  
  ## =============================== bands by VAR
  showBandsByVar(cwbnd, llp, logtime=T, varStr='isiOct',subs=llp$xpp$futureISI>.750&llp$xpp$trialType==1,
                 fileName='Figure4a_Bands(S)_SOA.eps', freqList=freqList, trngList=trngList)
  
  showBandsByVar(cwbnd, llp, logtime=T, varStr='ampInd',subs=llp$xpp$futureISI>.750&llp$xpp$trialType==1,
                 fileName='Figure4a_Bands(S)_AMP.eps', freqList=freqList, trngList=trngList)
  
  
  showBandsByVar(cwbnd, llp, logtime=F, varStr='isiOct',subs=llp$xpp$futureISI>.750&llp$xpp$trialType==1,
                 fileName='Figure4b_Bands(S)_SOA.eps', freqList=freqList, trngList=trngList)
  
  showBandsByVar(cwbnd, llp, logtime=F, varStr='ampInd',subs=llp$xpp$futureISI>.750&llp$xpp$trialType==1,
                 fileName='Figure4b_Bands(S)_AMP.eps', freqList=freqList, trngList=trngList)

  
  showMUAByVar(mua, logtime=T, varStr='isiOct',subs=llp$xpp$futureISI>.750&llp$xpp$trialType==1,
               fileName='Figure4c_Bands(M)_SOA.eps', freqList=freqList, trngList=trngList)
  
  showMUAByVar(mua, logtime=T, varStr='ampInd',subs=llp$xpp$futureISI>.750&llp$xpp$trialType==1,
               fileName='Figure4c_Bands(M)_AMP.eps', freqList=freqList, trngList=trngList)
  
  ## ============================= extract scalars
  freqList
  trngList
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
  
  
  sL <- getSLix( what=slList$what[1], when=slList$trng[1] ,xtn=slList$xtn[1], freqList=freqList, trngList=trngList )
  for (ix in 2:dim(slList)[1])
    sL <- rbind(sL, getSLix( what=slList$what[ix], when=slList$trng[ix],xtn=slList$xtn[ix], freqList=freqList, trngList=trngList ))
  
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
  saveDir  <- paste(rda(),'/', df$exp[ind],sep='')
  
  if (!exists( saveDir ) )
    system(paste('mkdir -p ', saveDir, sep=''))
  
  save( file=paste(saveDir ,'/features.rda',sep=''), xpp)
  

  
  ## Compare IC weights to MUA depth profile
  muaProfile <- apply( extractScalars(td=mua, trng=c(20, 40), ext='e'), 2, mean)
  muaProfile <- muaProfile/max(abs(muaProfile))
  showICAcomponents(icv, muaProfile=muaProfile,srtICA=T,csdWeights=T, showFFT=T, showEvokedResponse=T, includeexp=F, eps=T, pp=pp, fxtn='_muaProfile',logtime=T)
  
  
  ## =========================== predictors without fit
  eps.file('Figure5a_Scalars_SOA_AMP.eps',pdf=T,s=2,mf=c(6,dim(sL)[1]%/%6),ratio=1)
  par( mfrow=c(dim(sL)[1]%/%6,6), mar=c(4,4,2,1) )
  disc <- lapply(1:dim(sL)[1], showTD,  pM=predMat, fM=NULL, xpp=xpp)
  dev.off()
  
  eps.file('Figure5a_Scalars2a_SOA_AMP.eps',pdf=T,s=2,mf=c(6,dim(sL)[1]%/%6),ratio=1)
  par( mfrow=c(dim(sL)[1]%/%6,6), mar=c(4,4,2,1) )
  disc <- lapply(1:dim(sL)[1], showTD2,  pM=predMat, fM=NULL, xpp=xpp, IVstr1='isiOct', IVstr2='ampInd', colFUN=soa.colors)
  dev.off()
  
  eps.file('Figure5a_Scalars2b_AMP_SOA.eps',pdf=T,s=2,mf=c(6,dim(sL)[1]%/%6),ratio=1)
  par( mfrow=c(dim(sL)[1]%/%6,6), mar=c(4,4,2,1) )
  disc <- lapply(1:dim(sL)[1], showTD2,  pM=predMat, fM=NULL, xpp=xpp, IVstr1='ampInd', IVstr2='isiOct', colFUN=bluered.colors)
  dev.off()
  
  
  eps.file('Figure5b_MUA_SOA_AMP.eps',pdf=T,s=2,mf=c(8,dim(predMatMua)[2]%/%8),ratio=1)
  par( mfrow=c(dim(predMatMua)[2]%/%8,8), mar=c(4,4,2,1) )
  disc <- lapply(1:dim(predMatMua)[2], showTD,  pM=predMatMua, fM=NULL, xpp=xpp)
  dev.off()
  
  
  eps.file('Figure5b_MUA2a_SOA_AMP.eps',pdf=T,s=2,mf=c(6,dim(sL)[1]%/%6),ratio=1)
  par( mfrow=c(dim(predMatMua)[2]%/%8,8), mar=c(4,4,2,1) )
  disc <- lapply(1:dim(predMatMua)[2], showTD2,  pM=predMatMua, fM=NULL, xpp=xpp, IVstr1='isiOct', IVstr2='ampInd', colFUN=soa.colors)
  dev.off()
  
  eps.file('Figure5b_MUA2b_AMP_SOA.eps',pdf=T,s=2,mf=c(6,dim(sL)[1]%/%6),ratio=1)
  par( mfrow=c(dim(predMatMua)[2]%/%8,8), mar=c(4,4,2,1) )
  disc <- lapply(1:dim(predMatMua)[2], showTD2,  pM=predMatMua, fM=NULL, xpp=xpp, IVstr1='ampInd', IVstr2='isiOct', colFUN=bluered.colors)
  dev.off()
  
  ## ====================================== run model fits for each predictor
  source( paste(baseDir,'ampOdd_click/r-files/functions/functions_ampOdd_presynapticPlasticity.R',sep='') )
  tn     <- llp
  tn$xpp <- xpp
  #ampFit <- callSTPSP(tn,   DVstr='bspa1', gd=T, show=F, showCh=1, who='all', BYstr='exp', 
  #                    FUN=cumsumamp, minISI=.190, maxISI=12.8, startFromScratch=F, eps=F)
  
  FUN        <- cumsumamp
  ampFitList <- lapply( dimnames(predMat)[[2]], function(dvs){
    callSTPSP(tn,   DVstr=dvs, gd=T, show=T, showCh=1, who='all', BYstr='exp', FUN=FUN, minISI=.190, maxISI=21, startFromScratch=F, eps=T)})
  
  
  extractResiduals   <- function(n){ampFitList[[n]]$res$lm1$residuals}
  residMat           <- sapply(1:dim(sL)[1], extractResiduals)
  dimnames(residMat) <- dimnames(predMat)
  
  
  extractFitted    <- function(n){ampFitList[[n]]$res$lm1$fitted}
  fitMat           <- sapply(1:dim(sL)[1], extractFitted)
  dimnames(fitMat) <- dimnames(predMat)
  
  ## --------------------- mua
  ampFitListMua <- lapply( dimnames(predMatMua)[[2]], function(dvs){
    callSTPSP(tn,   DVstr=dvs, gd=T, show=T, showCh=1, who='all', BYstr='exp', FUN=FUN, minISI=.190, maxISI=21, startFromScratch=T, eps=T)})
  
  
  extractResiduals      <- function(n){ampFitListMua[[n]]$res$lm1$residuals}
  residMatMua           <- sapply(1:dim(predMatMua)[2], extractResiduals)
  dimnames(residMatMua) <- dimnames(predMatMua)
  
  extractFitted       <- function(n){ampFitListMua[[n]]$res$lm1$fitted}
  fitMatMua           <- sapply(1:dim(predMatMua)[2], extractFitted)
  dimnames(fitMatMua) <- dimnames(predMatMua)
  
  
  ## ======================================= plot dependence on SOA and AMP for the scalars
  eps.file('Figure5a_Scalars_SOA_AMP_Fit.eps',pdf=T,s=2,mf=c(6,dim(sL)[1]%/%6),ratio=1)
  par( mfrow=c(dim(sL)[1]%/%6,6), mar=c(4,4,2,1) )
  disc <- lapply(1:dim(sL)[1], showTD,  pM=predMat, fM=fitMat, xpp=xpp)
  dev.off()
  
  eps.file('Figure5b_MUA_SOA_AMP_Fit.eps',pdf=T,s=2,mf=c(8,dim(predMatMua)[2]%/%8),ratio=1)
  par( mfrow=c(dim(predMatMua)[2]%/%8,8), mar=c(4,4,2,1) )
  disc <- lapply(1:dim(predMatMua)[2], showTD,  pM=predMatMua, fM=fitMatMua, xpp=xpp)
  dev.off()
  
  
  
  ## ================================================================ raw correlation matrix
  ## ============================================================================== noise correlation
  
  disc <- analyzeCorrelationStructure(predMat   , predMat,    fileName='Figure6a_Cor(CSD)_raw.eps'  )
  disc <- analyzeCorrelationStructure(predMat   , predMatMua, fileName='Figure6a_Cor(CSD,MUA)_raw.eps' )
  disc <- analyzeCorrelationStructure(predMatMua, predMatMua, fileName='Figure6a_Cor(MUA)_raw.eps' )
  
  disc <- analyzeCorrelationStructure(residMat   , residMat,   fileName='Figure6b_Cor(CSD)_noise.eps')
  disc <- analyzeCorrelationStructure(residMat   , residMatMua,fileName='Figure6b_Cor(CSD,MUA)_noise.eps')
  disc <- analyzeCorrelationStructure(residMatMua, residMatMua,fileName='Figure6b_Cor(MUA)_noise.eps')
  
  disc <- analyzeCorrelationStructure(fitMat, fitMat,fileName='Figure6c_Cor(CSD)_fitted.eps',zlimCor=1)
  disc <- analyzeCorrelationStructure(fitMat, fitMatMua,fileName='Figure6c_Cor(CSD,MUA)_fitted.eps',zlimCor=1)
  disc <- analyzeCorrelationStructure(fitMatMua, fitMatMua,fileName='Figure6c_Cor(MUA)_fitted.eps',zlimCor=1)
  
  ## ========================================== SVM-based decoding
  
  tx          <- xpp[, names(xpp) %in% c(dimnames(predMat)[[2]], dimnames(predMatMua)[[2]]) ]
  
  eps.file('Figure7_Decoding_all.eps',pdf=T,s=2,mf=c(2,1),ratio=1)
  par( mfrow=c(1,2), mar=c(4,4,2,1) )
  
  tx$category <- as.factor(xpp$ampInd)
  decodeSVM(tx,subs=xpp$ISI<16)
  
  tx$category <- as.factor(xpp$isiOct)
  decodeSVM(tx,subs=xpp$ISI<16)
  dev.off()
  
  
}