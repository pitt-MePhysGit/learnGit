FLA_SUA <- function(mua, lfp, cmpStr='n1',eps=TRUE, DVsdscale=0 ){

  doLFP  <- FALSE#!is.null(lfp)
  doPSTH <- FALSE#!is.null(mua$psth)
  
  doFigures   <- TRUE
  doTrStats   <- FALSE
  doSTPSP     <- FALSE
  doScalarSUA <- FALSE
  
  resList <- list()
  
  who <- 'all'
  if(!is.null(mua$who))
    who <- mua$who
  
  trialTypeIndex <- c(1)
  if (who=='all')
    trialTypeIndex <- c(1,2)
  if (who=='regular')
    trialTypeIndex <- c(2)
  
  # prepare directory for figures
  ppDir   <- paste(pp(),'/',mua$exp,'/ch_',mua$ch,'_',mua$cell,'/',mua$who,'/',sep='')
  #ppDir   <- paste(pp(),'/', mua$exp,'/',mua$who,'/',sep='')
  if(!exists(ppDir))
    system(paste('mkdir -p',ppDir) )
  
  ## 
  thisDir   <- paste(mua$exp,'/ch_',mua$ch,'_',mua$cell,'/',mua$who,sep='')
  baseName  <- paste( mua$exp,'/ch_',mua$ch,'_',mua$cell, '/', mua$who, '/', mua$exp,'_ch_',mua$ch,'_',mua$cell,sep='')
  
  #extn        <- paste(mua$animal,'/FLA/NoModel_',mua$exp,'_',cmpStr, sep='')
  resFileName <- paste('SUA_',mua$exp,'_ch_',mua$ch,'_',mua$cell,'_',cmpStr,'_', mua$who,'.rda',sep='')
  
  tmeanPSTH <- function(psth, FUN=mean){
    print(dim(psth))
    ## little helper function that returns valid mean psths even
    ## if the number of repetitions is one or zero
    res <- array(NA,c(dim(psth)[c(2,3)]))  
    
    if(is.array(psth)&length(dim(psth))==3 )
      res <- apply(psth,c(2,3),FUN)
    
    if( is.array(psth) & length(dim(psth))==2 & identical(FUN,mean))
      res <- psth
    
    return(res)
  }
  
  
  ## ================================== potential outliers
  mxmua             <- 10
  tmp               <- apply(mua$lfp,c(1,3),max) - apply(mua$lfp,c(1,3),min)
  mua$xpp$p2p       <- apply( tmp, 1, max )
  mua$xpp$valInd    <- mua$xpp$p2p < mxmua & mua$xpp$trialType%in%trialTypeIndex;
  mua$xpp$valInd[1] <- F
  valIndMUA         <- mua$xpp$valInd
   
  mxlfp             <- 650
  tmp               <- apply(lfp$lfp,c(1,3),max) - apply(lfp$lfp,c(1,3),min)
  lfp$xpp$p2p       <- apply( tmp, 1, max )
  lfp$xpp$valInd    <- mua$xpp$p2p < mxlfp & lfp$xpp$trialType%in%trialTypeIndex;
  lfp$xpp$valInd[1] <- F
  valIndLFP         <- lfp$xpp$valInd
  
  
  if(doFigures){
    xlim <- c(-50,200)
    
    mua$xpp$spk <- spike.count( mua, range=c(-1000,1000) , rate=F  )
    eps.file(paste(baseName,'_spikeCount.eps',sep=''), pdf=T,s=6, ratio=1/2)
    plot(mua$xpp$spk, col=1+mua$xpp$trialType%in%trialTypeIndex)
    abline(v=seq(1,by=100,length=20))
    dev.off()
    
    
    ##  ----------------------------------------------------------------------  master overview figure
    eps.file(paste(baseName,'_AOC_timeseries_MUA.eps',sep=''), pdf=T,s=4)
    disc <- figure_AOC_timeseries(mua, whoStr='lfp', rngT=c(-10,100))
    dev.off()

    
    eps.file(paste(baseName,'_AOC_timeseries_PSTH.eps',sep=''), pdf=T,s=4)
    disc <- figure_AOC_timeseries(mua, whoStr='psth', rngT=c(-10,100))
    dev.off()
    
    
    eps.file(paste(baseName,'_AOC_timeseries_LFP.eps',sep=''), pdf=T,s=4)
    disc <- figure_AOC_timeseries(lfp, whoStr='lfp', rngT=c(-25,250), chnd=2)
    dev.off()
    
    eps.file(paste(baseName,'_AOC_timeseries_CSD.eps',sep=''), pdf=T,s=4)
    disc <- figure_AOC_timeseries(lfp, whoStr='lfp', rngT=c(-25,250), chnd=4)
    dev.off()
    
    
    ## ==================================== GRAND AVERAGES SUA and MUA
    psthtaxis    <- 1000*as.numeric(dimnames(mua$psth)[[2]])
    
    resList$mnGAPSTH <- apply(mua$psth[which(valIndMUA),,1],2,mean)
    eps.file(paste(baseName,'_PSTH_GA.eps',sep=''), pdf=T,s=3)
    mean.plot(mua$psth[which(valIndMUA),,1],psthtaxis, widthFUN=na.se,add=F,lwd=2,xlim=xlim)
    abline(v=0);abline(h=0)
    dev.off()
    
    resList$mnGAEPSP <- apply(mua$epsp[which(valIndMUA),,1],2,mean)
    eps.file(paste(baseName,'_EPSP_GA.eps',sep=''), pdf=T,s=3)
    mean.plot(mua$epsp[which(valIndMUA),,1],psthtaxis, widthFUN=na.se,add=F,lwd=2,xlim=xlim)
    abline(v=0);abline(h=0)
    dev.off()
    
    resList$mnGAMUA <- apply(mua$lfp[which(valIndMUA),,1],2,mean)
    eps.file(paste(baseName,'_MUA_GA.eps',sep=''), pdf=T,s=3)
    mean.plot(mua$lfp[which(valIndMUA),,1],mua$taxis, widthFUN=na.se,add=F,lwd=2,xlim=xlim)
    abline(v=0);abline(h=0)
    dev.off()
    
    
    ## =================================== GRAND AVERAGES LFP and CSD
    
    resList$mnGALFP <- apply(lfp$lfp[which(valIndLFP),,2],2,mean)
    eps.file(paste(baseName,'_LFP_GA.eps',sep=''), pdf=T,s=3)
    mean.plot(lfp$lfp[which(valIndLFP),,2],lfp$taxis, widthFUN=na.se,add=F,lwd=2,xlim=xlim )
    abline(h=0);abline(v=0)
    dev.off()
    
    resList$mnGACSD <- apply(lfp$lfp[which(valIndLFP),,4],2,mean)
    eps.file(paste(baseName,'_CSD_GA.eps',sep=''), pdf=T,s=3)
    mean.plot(lfp$lfp[which(valIndLFP),,4],lfp$taxis, widthFUN=na.se,add=F,lwd=2,xlim=xlim)
    abline(h=0);abline(v=0)
    dev.off()
    
    
    ## ============================================= GRAND AVERAGE OVERLAY FIGURE
    eps.file(paste(baseName,'_ALL_GA.eps',sep=''), pdf=T,s=3)
    mean.plot(mua$epsp[which(valIndMUA),,1]/max(abs(resList$mnGAEPSP)),psthtaxis, 
              widthFUN=na.se,add=F,lwd=4, ylim=c(-1,1),xlim=xlim)
    
    mean.plot(mua$lfp[which(valIndMUA),,1]/max(abs(resList$mnGAMUA)),lfp$taxis, 
              widthFUN=na.se,add=T,lwd=2,col='black')
    
    mean.plot(lfp$lfp[which(valIndLFP),,2]/max(abs(resList$mnGALFP)),lfp$taxis, 
              widthFUN=na.se,add=T,lwd=2,col='gray')
    
    #mean.plot(lfp$lfp[which(valIndLFP),,4]/max(abs(resList$mnGACSD)),lfp$taxis, 
    #          widthFUN=na.se,add=T,lwd=2,col='red')
    
    mean.plot(mua$epsp[which(valIndMUA),,1]/max(abs(resList$mnGAEPSP)),psthtaxis, 
              widthFUN=na.se,add=T,lwd=2)
    
    abline(h=0);abline(v=0)
    dev.off()
    
    
    ## =======================================================
    print('Computing average ISIbin')
    extn <- 'ISIbin_'
    
    binKerISI       <- mua$xpp$isiOct 
    
    mn        <- tapply(1:dim(mua$psth)[1], list(binKerISI),function(ind){ tmeanPSTH(mua$psth[ ind[ind%in%which(valIndMUA)] ,,,drop=F]) })
    nullInd   <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
    for ( i in nullInd)
      mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
    
    Naverage        <- table(binKerISI)
    print(Naverage)
    
    ylim <- na.range(c(unlist(mn)))
    
    eps.file(paste(baseName,'_PSTH_ByISIOct.eps',sep=''), pdf=T,s=3)
    plot(c(-100),c(-100),ylim=ylim, xlim=c(0.5,7) )
    abline(h=0)
    for (i in seq(6,to=1,by=-1))
      mean.plot(mua$psth[which(as.numeric(mua$xpp$isiOct)==i&valIndMUA),,1],i+psthtaxis/1000, widthFUN=na.se,
                add=ifelse(i==6,T,T), col=soa.colors(7)[i], lwd=2 )
    dev.off()
    
    ## ================================= one-line TRF
    eps.file(paste(baseName,'_PSTH_ByAmp.eps',sep=''), pdf=T,s=3)
    plot(c(-100),c(-100), xlim=c(0.5, 6),ylim=ylim )
    abline(h=0)
    for (i in seq(1,to=5,by=1))
      mean.plot(mua$psth[which(as.numeric(mua$xpp$ampInd)==i&valIndMUA),,1],i+psthtaxis/1000, widthFUN=na.se,
                add=T, col=bluered.colors(5)[i], lwd=2 )
    dev.off()
    
    
    
    ## ==== MUA
    mn        <- tapply(1:dim(mua$lfp)[1], list(binKerISI),function(ind){ tmeanPSTH(mua$lfp[ ind[ind%in%which(valIndMUA)] ,,,drop=F]) })
    nullInd   <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
    for ( i in nullInd)
      mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
    
    Naverage        <- table(binKerISI)
    print(Naverage)
    
    ylim <- na.range(c(unlist(mn)))
    
    eps.file(paste(baseName,'_MUA_ByISIOct.eps',sep=''), pdf=T,s=3)
    plot(c(-100),c(-100), xlim=xlim,ylim=ylim);abline(h=0)
    for (i in seq(6,to=1,by=-1))
      mean.plot(mua$lfp[which(as.numeric(mua$xpp$isiOct)==i&valIndMUA),,1],mua$taxis, widthFUN=na.se,
                add=ifelse(i==6,T,T), col=soa.colors(6)[i], lwd=2 )
    dev.off()
    
    
    
    ## --------------- mua by amplitude
    mn        <- tapply(1:dim(mua$lfp)[1], list(binKerISI),function(ind){ tmeanPSTH(mua$lfp[ ind[ind%in%which(valIndMUA)] ,,,drop=F]) })
    nullInd   <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
    for ( i in nullInd)
      mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
    
    Naverage        <- table(binKerISI)
    print(Naverage)
    
    eps.file(paste(baseName,'_MUA_ByAmpInd.eps',sep=''), pdf=T,s=3)
    plot(c(-100),c(-100), xlim=xlim,ylim=ylim);abline(h=0)
    for (i in seq(5,to=1,by=-1))
      mean.plot(mua$lfp[which(as.numeric(mua$xpp$ampInd)==i&valIndMUA),,1],mua$taxis, widthFUN=na.se,
                add=ifelse(i==6,T,T), col=bluered.colors(5)[i], lwd=2 )
    dev.off()
    
    
    
    
    
    ## ================================= one-line TRF
    eps.file(paste(baseName,'_MUA_ByFreq.eps',sep=''), pdf=T,s=3)
    plot(c(-100),c(-100), xlim=c(0, 12),ylim=ylim );abline(h=0)
    for (i in seq(1,to=11,by=1))
      mean.plot(mua$lfp[which(as.numeric(mua$xpp$toneInd)==i&valIndMUA),,1],i+mua$taxis/1000, widthFUN=na.se,
                add=T, col='black', lwd=2 )
    dev.off()
    
    
  }# end doFigures
  
  
  muarng     <- c(0.015,0.070)  
  mua$xpp$DV         <- spike.count(mua,range=muarng,     rate=T)
  
  mua$xpp$FRbl       <- spike.count(mua,range=c(-0.100,0),rate=T)
  mua$xpp$FR         <- spike.count(mua,range=muarng,     rate=T)
  mua$xpp$FRnrm      <- mua$xpp$FR - mua$xpp$FRbl
  
  correctForSlowFluctuations <- TRUE
  if(correctForSlowFluctuations){
    mavDV <- get.mav(mua$xpp$FR, 1:dim(mua$xpp)[1], eval.at=seq(1,dim(mua$xpp)[1],by=1 ), show=T, distribution='Gaussian',boot=F, width=100)
    mua$xpp$dmFR <- mua$xpp$FR - mavDV$mav + mean(mavDV$mav)
    
    mavDV <- get.mav(mua$xpp$FRnrm, 1:dim(mua$xpp)[1], eval.at=seq(1,dim(mua$xpp)[1],by=1 ), show=T, distribution='Gaussian',boot=F, width=100)
    mua$xpp$dmFRnrm <- mua$xpp$FRnrm - mavDV$mav + mean(mavDV$mav)
    
    mavDV <- get.mav(mua$xpp$FRbl, 1:dim(mua$xpp)[1], eval.at=seq(1,dim(mua$xpp)[1],by=1 ), show=T, distribution='Gaussian',boot=F, width=100)
    mua$xpp$dmFRbl <- mua$xpp$FRbl - mavDV$mav + mean(mavDV$mav)
  }
  
  
  
  if (doTrStats){
    ## ================================================================== time-resolved test
    getRSTrBaseline <- function(input,subs=NULL){ 
      width     <- .020
      eval.at   <- seq(-0.010, 0.950,by=.002)
      variables <- c('sp','isiOct')
      #subs      <- mua$xpp$p2p<400 & !is.na(mua$xpp$isiOct) & mua$xpp$trialType==1 & mua$xpp$ISI>1
      
      this.test <- function(input){
        df  <- input$xpp[,variables]
        lm0 <- lm(sp ~  -1, data=df)
        lm1 <- lm(sp ~   1, data=df)
        
        res <- anova(lm0,lm1)
        return(res)
      }
      tmp <- tres.test(input, main='ISI', eval.at=eval.at, test=this.test, width=width, subs=subs, blcor=T)
      
      sigEpochs        <- get.epoch( tmp$P<1e-4 )
      sigEpochs$length <- eval.at[sigEpochs$offset] - eval.at[sigEpochs$onset]
      sigEpochs        <- subset(sigEpochs, sigEpochs$length>0.010 & eval.at[sigEpochs$onset]>0 )
      tmp$sigEpochs    <- sigEpochs
      
      return(tmp)
    }
    ## ===================================================================== extract relevant component
    # find response onset and duration. do this independent of condition to enable comparison
    trResList         <- getRSTrBaseline(mua,subs=mua$xpp$p2p<mxmua)
    resList$trResList <- trResList
    
    if(eps)
      eps.file( paste(baseName,'_tresMeanFR.eps',sep=''),ratio=1/.66, s=2 )
    show.tres.test( trResList)
    if(eps)
      dev.off()
    
    #muarng    <- c(0.005,0.014)
    
    if (!is.null(dim(trResList$sigEpochs))){ 
      if(!dim(trResList$sigEpochs)[1]==0){
        tmpRng    <- trResList$eval.at[ c(unlist(trResList$sigEpochs[1,c('onset','offset')])) ]
        if (tmpRng[1]<0.09)
          muarng    <- (tmpRng%>=%0.015)%<=%0.100
      }
    }
    resList$muarng <- muarng
    
    #mua$xpp$FRbl       <- spike.count(mua,range=c(-0.100,0),rate=T)
    #mua$xpp$FR         <- spike.count(mua,range=muarng,     rate=T)
    #mua$xpp$FRnrm      <- mua$xpp$FR - mua$xpp$FRbl
    
    #correctForSlowFluctuations <- TRUE
    #if(correctForSlowFluctuations){
    #  mavDV <- get.mav(mua$xpp$FR, 1:dim(mua$xpp)[1], eval.at=seq(1,dim(mua$xpp)[1],by=1 ), show=T, distribution='Gaussian',boot=F, width=100)
    #  mua$xpp$dmFR <- mua$xpp$FR - mavDV$mav + mean(mavDV$mav)
    #  
    #  mavDV <- get.mav(mua$xpp$FRnrm, 1:dim(mua$xpp)[1], eval.at=seq(1,dim(mua$xpp)[1],by=1 ), show=T, distribution='Gaussian',boot=F, width=100)
    #  mua$xpp$dmFRnrm <- mua$xpp$FRnrm - mavDV$mav + mean(mavDV$mav)
    #  
    #  mavDV <- get.mav(mua$xpp$FRbl, 1:dim(mua$xpp)[1], eval.at=seq(1,dim(mua$xpp)[1],by=1 ), show=T, distribution='Gaussian',boot=F, width=100)
    #  mua$xpp$dmFRbl <- mua$xpp$FRbl - mavDV$mav + mean(mavDV$mav)
    #}
    
    # use FR as the standin for DV. 
    mua$xpp$DV         <- mua$xpp$FR
    
    ## time-resolved general linear model 
    mua$xpp$isiLog <- log(mua$xpp$ISI)
    
    
    #variables=c('sp','isiLog','pitchFreq','deltaPitchOct','absDeltaPitchOct','trialNr','seqNr')
    #resList$trStats <- trStats.RSkernel(mua,width=.030,eval.at=mua$taxis,variables=variables,
    #                                    subs=valIndMUA,extn='',show=T, fix=FALSE, what='spk')
    
    if(1==1){# 
      variables=c('sp','isiLog','pitchFreq','deltaPitchOct','absDeltaPitchOct','trialNr','seqNr','toneType')
      resList$trStats <- trStats.RSkernel(mua,width=.050,eval.at=mua$taxis,variables=variables,
                                          subs=valIndMUA,extn=paste('_cell_',mua$cell,'_epsp_',who,sep=''),show=T, fix=FALSE, what='epsp')
      
      if(eps)
        eps.file( paste(baseName,'_tres_Pmap.eps',sep=''),ratio=1/.66, s=2)
      show.tres.glm(resList$trStats,showPar=F,bfh=T,yaxt=T)
      if(eps)
        dev.off()
    }
  }
  
  ## subset the data
  maxISI    <- 25
  minISI    <-  0.195
  maxp2p    <-  10
  DVsdscale <- NA
  maxPowEOG <- Inf
  sfc       <- TRUE
  
  if(mua$animal=='Rockey') # & what %in% c('raw','mlp','mbp','bsp'))
    maxPowEOG <- 15000
  
  if(mua$animal=='Walter')# & what %in% c('raw','mlp','mbp','bsp'))
    maxPowEOG <- 5000
  
  if(mua$animal=='Jesse')# & what %in% c('raw','mlp','mbp','bsp'))
    maxPowEOG <- 2000
  
  if(mua$animal=='Sam')# & what %in% c('raw','mlp','mbp','bsp'))
    maxPowEOG <- 15000
  
  x <- mua$xpp
  if (!is.na(DVsdscale)){ 
    mDV       <- quantile( x$DV, prob=c(0.25,0.5,0.75) )
    valRng    <- mDV[2] + c(-DVsdscale,DVsdscale) * diff(mDV) 
    print(mDV)
  }
  
  if (is.na(DVsdscale)){ 
    valRng    <- quantile( x$DV, prob=c(0,1) )
  }
  
  mua$xpp$firstTrialInDay <- diff( c(0,rep(1,length(mua$xpp$ISI)) ))
  
  ybin <- ( mua$xpp$p2p<maxp2p                     &   # present in STPSP
              mua$xpp$ISI>minISI                   &   # present in STPSP
              mua$xpp$ISI<maxISI                   &   # present in STPSP
              mua$xpp$firstTrialInDay==0           &   # present in STPSP exclude the very first tone of each day
              mua$xpp$trialType%in%trialTypeIndex  &   # present in STPSP
              !is.na(mua$xpp$binISI)               &   
              mua$xpp$DV>=valRng[1]                 &   # present in STPSP
              mua$xpp$DV<=valRng[2]                 &   # present in STPSP
              mua$xpp$powEOG<maxPowEOG)                # present in STPSP
  
  print(sum(ybin)/length(ybin))
  
  # subset from STPSP  
  y           <- subset(mua$xpp, ybin)
  
  ## ===================================================================== do rasters
  spkcnt <- get.mav( spike.count(mua,range=c(-1,1),rate=F), 1:length(mua$spk), show=F)
  muacnt <- get.mav( spike.count(mua,range=c(-1,1),rate=F, whoStr='mua'), 1:length(mua$spk), show=F)
  
  eps.file(paste(baseName,'_SPKCountByTrialNr.eps',sep=''), pdf=T,s=3,ratio=1/0.67)
  mav.show(muacnt, points=F,col='black', ylim=c(0,na.max(muacnt$mav)) );abline(h=0)
  mav.show(spkcnt, points=F,col='red',add=T)
  dev.off()
  
  ## =========== example rasters ISI and Frequency
  ## ================================ 
  
  eps.file(paste(baseName,'_sampleRaster.eps',sep=''), pdf=T,s=3,ratio=1/0.67)
  raster.plot(mua$spk[sample( which(ybin),1000)], xlim=c(-.05,.25) ,xlab='time after tone onset',col='black',yaxt='n',range=c(0,1))
  abline(v=0)
  dev.off()
  
  sampleValidIndeces <- function(nd,nsmp=20){ 
    valnd <- nd[ which(ybin[nd]) ]
    return(sample(valnd,min(nsmp,length(valnd))))
  }
  
  tst    <- tapply( 1:dim(mua$xpp)[1] , mua$xpp$isiOct,  sampleValidIndeces,nsmp=100)
  allInd <- unlist(tst[1:6]) 
  tlim <- c(-0.05, 0.25)
  eps.file(paste(baseName,'_sampleRaster_ISI.eps',sep=''), pdf=T,s=4,ratio=1)
  raster.plot(mua$spk[allInd], xlim=tlim, gind=mua$xpp$isiOct[allInd],gsort=T,xlab='time after tone onset',
              col=soa.colors(7),ylab='SOA [sec]',yaxt='n',range=c(0,1),cex=1.0)
  #axis(2,at=squeeze(fat,from.range=frng,to.range=c(0,1)),labels=fpretty)
  abline(v=0)
  dev.off()
  
  tst    <- tapply( 1:dim(mua$xpp)[1], mua$xpp$ampInd,  sampleValidIndeces, nsmp=100)
  allInd <- unlist(tst) 
  eps.file(paste(baseName,'_sampleRaster_AMP.eps',sep=''), pdf=T,s=4,ratio=1)
  raster.plot(mua$spk[allInd], xlim=c(-.05,.250), gind=mua$xpp$ampInd[allInd],gsort=T,xlab='time after tone onset',
              col=bluered.colors(5),ylab='Amplitude',yaxt='n',range=c(0,1),cex=1.0)
  #axis(2,at=squeeze(fat,from.range=frng,to.range=c(0,1)),labels=fpretty)
  abline(v=0)
  dev.off()
  
  ### ======================================================================================= set of analyses using a scalar variable such as mean mua activity in specific time-window
  if(doScalarSUA){
    
    cmpStrDF     <- data.frame( what='sua', cmpStrName='SUA0',  tmin=.010,  tmax=.018,  stringsAsFactors=FALSE) 
    cmpStrDF[2,] <- data.frame( what='sua', cmpStrName='SUA1',  tmin=.025,  tmax=.045,  stringsAsFactors=FALSE)
    cmpStrDF[3,] <- data.frame( what='sua', cmpStrName='SUA1a', tmin=.020,  tmax=.030,  stringsAsFactors=FALSE)
    cmpStrDF[4,] <- data.frame( what='sua', cmpStrName='SUA1b', tmin=.030,  tmax=.040,  stringsAsFactors=FALSE)
    cmpStrDF[5,] <- data.frame( what='sua', cmpStrName='SUA1c', tmin=.040,  tmax=.050,  stringsAsFactors=FALSE)
    cmpStrDF[6,] <- data.frame( what='sua', cmpStrName='SUA2' , tmin=.065,  tmax=.150, stringsAsFactors=FALSE)
    cmpStrDF[7,] <- data.frame( what='sua', cmpStrName='SUA3' , tmin=.150, tmax=.250, stringsAsFactors=FALSE)
    cmpStrDF[8,] <- data.frame( what='sua', cmpStrName='SUAbl', tmin= -.100, tmax=0, stringsAsFactors=FALSE)
    cmpStrDF[9,] <- data.frame( what='sua', cmpStrName='SUAsig',tmin=resList$muarng[1], tmax=resList$muarng[2], stringsAsFactors=FALSE)
    
    cmpStrDF[10,]   <- data.frame(what='nsua', cmpStrName='nSUA0',  tmin=.010,  tmax=.018,  stringsAsFactors=FALSE) 
    cmpStrDF[11,] <- data.frame(  what='nsua',   cmpStrName='nSUA1',  tmin=.025,  tmax=.045,  stringsAsFactors=FALSE)
    cmpStrDF[12,] <- data.frame(  what='nsua',   cmpStrName='nSUA1a', tmin=.020,  tmax=.030,  stringsAsFactors=FALSE)
    cmpStrDF[13,] <- data.frame(  what='nsua',   cmpStrName='nSUA1b', tmin=.030,  tmax=.040,  stringsAsFactors=FALSE)
    cmpStrDF[14,] <- data.frame(  what='nsua',   cmpStrName='nSUA1c', tmin=.040,  tmax=.050,  stringsAsFactors=FALSE)
    cmpStrDF[15,] <- data.frame(  what='nsua',   cmpStrName='nSUA2' , tmin=.065,  tmax=.150, stringsAsFactors=FALSE)
    cmpStrDF[16,] <- data.frame(  what='nsua',   cmpStrName='nSUA3' , tmin=.0150, tmax=.250, stringsAsFactors=FALSE)
    cmpStrDF[17,] <- data.frame(  what='nsua',   cmpStrName='nSUAbl', tmin= -.100, tmax=0, stringsAsFactors=FALSE)
    cmpStrDF[18,] <- data.frame(  what='nsua',   cmpStrName='nSUAsig',tmin=resList$muarng[1], tmax=resList$muarng[2], stringsAsFactors=FALSE)
    
    for (cx in 1:dim(cmpStrDF)[1]){# loop that runs across all potential measures
      dat     <- mua
      DVstr   <- paste(mua$chStr, '_', mua$cell, '_', cmpStrDF$cmpStrName[cx], sep='')
      print(DVstr)
      
      trng <- c( cmpStrDF[cx,3],cmpStrDF[cx,4]) 
      dat$xpp$TMP   <- spike.count(mua,range=trng,     rate=T)
      dat$xpp$DV    <- dat$xpp$TMP 
      
      nd                  <- which(names(dat$xpp)=='TMP')
      names(dat$xpp)[nd]  <- DVstr
      
      dat$xpp$pitchInd       <- dat$xpp$toneInd   #(dat$xpp$pitch-1)%%11
      if ( cmpStrDF$cmpStrName[cx] %in% c('SUAbl','nSUAbl') ){     # shift pitchInd to reflect pitch of last tone
        dat$xpp$pitchInd         <- my.lag( dat$xpp$pitchInd) #                 (dat$xpp$pitch-1)%%11, 1)
        dat$xpp$deltaPitchInd    <- my.lag( dat$xpp$deltaPitchInd,1)
        dat$xpp$absDeltaPitchInd <- my.lag( dat$xpp$absDeltaPitchInd,1)
      }
      
      if( cmpStrDF$what[cx] %in% c('nsua')){
        mavDV            <- get.mav(dat$xpp[[DVstr]], 1:dim(dat$xpp)[1], eval.at=seq(1,dim(dat$xpp)[1],by=1 ), show=T, distribution='Gaussian',boot=F, width=100)
        dat$xpp[[DVstr]] <- dat$xpp[[DVstr]] - mavDV$mav + mean(mavDV$mav)
      }
      
      dat$xpp$pitchIndFactor <- as.factor(dat$xpp$pitchInd)
      
      # --------------- run analyzeRSscalar
      valBin <- getValidTrialsRS(dat, who=who,DVstr=DVstr, DVsdscale=Inf)
      
      if (max(dat$xpp[[DVstr]])==0)# quick and dirty way of preventing rank deficient linear models if there is not a single spike
        dat$xpp[[DVstr]][ which(valBin)[4] ] <- 1
      
      disc   <- analyzeRSscalar(dat, ybin=valBin, who=who, cmpStr=DVstr)
    }
  }# end doScalarSUA
  
  
  ## ========================================================================================
  if(doSTPSP){
    
    # fit model to non-baseline corrected component
    tmpDir <- paste(pp(),thisDir,'/STPSP',sep='' )
    if(!exists(tmpDir))
      system(paste('mkdir',tmpDir) )
    
    cellSpecifier <- paste('ch_',mua$ch,'_',mua$cell,'_',cmpStr,sep='')
    fe <- paste(cellSpecifier,'',sep='')
    
    sfc <- T
    
    print('running STPSP...')
    
    STPSPFR      <- callSTPSP( mua, gd=T, startFromScratch=sfc, DVstr='FR', who=who, minFrq=500, toneBy=1,NOct=3,fileExtn=fe,    
                               minISI=minISI, maxISI=maxISI, DVsdscale=DVsdscale, show=T, eps=T,baseName=thisDir)
    STPSPFRnrm   <- callSTPSP( mua, gd=T, startFromScratch=sfc, DVstr='FRnrm', who=who, minFrq=500, toneBy=1,NOct=3,fileExtn=fe,    
                               minISI=minISI, maxISI=maxISI, DVsdscale=DVsdscale, show=T, eps=T,baseName=thisDir)
    STPSPFRbl   <- callSTPSP( mua, gd=T, startFromScratch=sfc, DVstr='FRbl', who=who, minFrq=500, toneBy=1,NOct=3,fileExtn=fe,    
                              minISI=minISI, maxISI=maxISI, DVsdscale=DVsdscale, show=T, eps=T,baseName=thisDir)
    
    STPSPdmFR      <- callSTPSP( mua, gd=T, startFromScratch=sfc, DVstr='dmFR', who=who, minFrq=500, toneBy=1,NOct=3,fileExtn=fe,    
                                 minISI=minISI, maxISI=maxISI, DVsdscale=DVsdscale, show=T, eps=T,baseName=thisDir)
    STPSPdmFRnrm   <- callSTPSP( mua, gd=T, startFromScratch=sfc, DVstr='dmFRnrm', who=who, minFrq=500, toneBy=1,NOct=3,fileExtn=fe,    
                                 minISI=minISI, maxISI=maxISI, DVsdscale=DVsdscale, show=T, eps=T,baseName=thisDir)
    STPSPdmFRbl   <- callSTPSP( mua, gd=T, startFromScratch=sfc, DVstr='dmFRbl', who=who, minFrq=500, toneBy=1,NOct=3,fileExtn=fe,    
                                minISI=minISI, maxISI=maxISI, DVsdscale=DVsdscale, show=T, eps=T,baseName=thisDir)
    
    
    resList$STPSPparFR    <- STPSPFR$modelFit$par
    resList$STPSPparFRnrm <- STPSPFRnrm$modelFit$par
    resList$STPSPparFRbl  <- STPSPFRbl$modelFit$par
    
    resList$STPSPpardmFR    <- STPSPdmFR$modelFit$par
    resList$STPSPpardmFRnrm <- STPSPdmFRnrm$modelFit$par
    resList$STPSPpardmFRbl  <- STPSPdmFRbl$modelFit$par
    
    mua$xpp$toneIndFactor <- as.factor(mua$xpp$toneInd)
    
    print('running STPSPfreq')
    
    fe <- paste(cellSpecifier,'_frq3par',sep='')
    STPSPFrqFR       <- callFrqSTPSP( mua, gd=T, startFromScratch=sfc,predFUN=STPSPmodelFun_freq, DVstr='FR', who=who,fileExtn=fe,    
                                      minISI=minISI, maxISI=maxISI, DVsdscale=DVsdscale, show=T, eps=T,baseName=thisDir)
    STPSPFrqFRnrm    <- callFrqSTPSP( mua, gd=T, startFromScratch=sfc,predFUN=STPSPmodelFun_freq, DVstr='FRnrm', who=who,fileExtn=fe,    
                                      minISI=minISI, maxISI=maxISI, DVsdscale=DVsdscale, show=T, eps=T,baseName=thisDir)
    
    STPSPFrqdmFR       <- callFrqSTPSP( mua, gd=T, startFromScratch=sfc,predFUN=STPSPmodelFun_freq, DVstr='dmFR', who=who,fileExtn=fe,    
                                        minISI=minISI, maxISI=maxISI, DVsdscale=DVsdscale, show=T, eps=T,baseName=thisDir)
    STPSPFrqdmFRnrm    <- callFrqSTPSP( mua, gd=T, startFromScratch=sfc,predFUN=STPSPmodelFun_freq, DVstr='dmFRnrm', who=who,fileExtn=fe,    
                                        minISI=minISI, maxISI=maxISI, DVsdscale=DVsdscale, show=T, eps=T,baseName=thisDir)
    
    # for the baselineFR to be selective to a particular frequency, it needs to be selective to the 
    # pitch of the last sound, not the current sound.
    muatmp                      <- mua
    muatmp$xpp$toneInd          <- my.lag(mua$xpp$toneInd,1)
    muatmp$xpp$toneIndFactor    <- my.lag(mua$xpp$toneIndFactor,1)
    muatmp$xpp$deltaPitchOct    <- my.lag(mua$xpp$deltaPitchOct,1)
    muatmp$xpp$absDeltaPitchOct <- my.lag(mua$xpp$absDeltaPitchOct,1)
    
    STPSPFrqFRbl     <- callFrqSTPSP( muatmp, gd=T, startFromScratch=sfc,predFUN=STPSPmodelFun_freq, DVstr='FRbl', who=who,fileExtn=fe,    
                                      minISI=minISI, maxISI=maxISI, DVsdscale=DVsdscale, show=T, eps=T,baseName=thisDir)
    
    STPSPFrqdmFRbl     <- callFrqSTPSP( muatmp, gd=T, startFromScratch=sfc,predFUN=STPSPmodelFun_freq, DVstr='dmFRbl', who=who,fileExtn=fe,    
                                        minISI=minISI, maxISI=maxISI, DVsdscale=DVsdscale, show=T, eps=T,baseName=thisDir)
    
    resList$STPSPFrqparFR    <- STPSPFrqFR$modelFit$par
    resList$STPSPFrqparFRnrm <- STPSPFrqFRnrm$modelFit$par
    resList$STPSPFrqparFRbl  <- STPSPFrqFRbl$modelFit$par
    
    resList$STPSPFrqpardmFR    <- STPSPFrqdmFR$modelFit$par
    resList$STPSPFrqpardmFRnrm <- STPSPFrqdmFRnrm$modelFit$par
    resList$STPSPFrqpardmFRbl  <- STPSPFrqdmFRbl$modelFit$par
  }
  
  ## ================================================= save results
  save(file=paste(rda(),resFileName,sep=''), resList)
  
}
