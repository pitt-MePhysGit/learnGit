
FLA_VProbe_STPSP <- function(j,df=df,baseDir=baseDir,trda=trda,ind=ind){
  
  i <- ind[j]
  print(paste('starting process', i))
  print(df$session[i])
  
  library(R.matlab)
  library(rgenoud)
  library(lattice)
  library(akima)
  library(Rwave)
  library(MASS)
  
  
  #source(paste(baseDir, 'RSkernel/r-files/functions/source_all_functions_parallel.R',sep=''))
  source( paste(baseDir,'ampOdd_click/r-files/functions/functions_ampOddclick.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/functions/functions_AOCdual.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/functions/functions_VProbe_STPSP.R',sep='') )
  #source( paste(baseDir,'RSkernel/r-files/functions/FLA_VProbe_STPSP.R',sep='') )
  
  source( paste(baseDir,'r-compile/rUtils/R/rUtils.R',sep='') )
  source( paste(baseDir,'r-compile/mav.test/R/mav.test.R',sep='') )
  source( paste(baseDir,'r-compile/ePhysLab/R/ePhysLab.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/calcCSD.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/showCSDimage.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/showCSD.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/multiVplot.R',sep='') )
  
  source( paste(baseDir,'ampOdd_click/r-files/functions/prepSTPSPVProbeMUA.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/functions/functions_ampOdd_presynapticPlasticity.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/functions/trStats_ampOddClick.R',sep='') )
  
  source( paste(baseDir,'ampOdd_click/r-files/functions/getSessionAverage.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/functions/getValidTrialsAOC.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/functions/analyzeAOCscalar.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/functions_correctForRapitTransitions.R',sep='') )
  
  
  options(warn=0)
  
  mxmua <- 40
  mxlfp <- 400
  
  colList = c('black','blue','cyan','green','orange','red')
  
  doSessionAverage <- F  # simple averages across ISI, pitch, etc  
  CSDplots         <- F  # a couple of simple csd plots 
  doTrLM           <- F  # time-resolved analysis of variance
  doSTPSPbyChan    <- F  # STPSP model for one single outcome variable per channel
  doSTPSP          <- F  # STPSP model for one single outcome varialbe per probe
  doScalarMUA      <- T  # scalar values for MUA
  
  #print(df[i,])
  #print(i)
  
  baseName <- paste(df$session[i],'/',sep='')
  
  #pp   <- function(){paste(baseDir, 'RSkernel/word/VProbe/plots/',sep='')}
  #.GlobalEnv$rda  <- function(){'~/Dropbox/RSkernel/Rdata/'}
  .GlobalEnv$rda  <- trda
  
  dxt <- ifelse(df$drug[i]=='N','', paste('_',df$drug[i],sep=''))
  pxt <- ifelse(df$probe[i]==1 ,'', paste('_',df$probe[i],sep=''))
  
  .GlobalEnv$pp   <- function(){paste(baseDir, 'ampOdd_click/word/VProbe/plots/',df$session[i],dxt,pxt,'/',sep='')}
  print(pp())
  
  interpChanInd <- c();#c(16)
  
  if (!file.exists(pp()))
    system( paste('mkdir', pp()) )
  
  #source(paste(baseDir, 'RSkernel/r-files/functions/source_all_functions.R',sep=''))
  # read the data
  #browser()
  chInd      <- df$startVIndex[i] + 1:df$numChan[i]
  chStr      <- paste('ch',chInd,sep='')
  muaR       <- loadAOC( df[i,] , what='mua', chStr=chStr, trda=trda,blcorr=F  )
  mua        <- loadAOC( df[i,] , what='mua', chStr=chStr, trda=trda,blcorr=T  )
  lfp        <- loadAOC( df[i,] , what='raw', chStr=chStr, trda=trda  )
  eeg        <- loadAOC( df[i,] , what='raw', chStr=paste('ch',1:33,sep=''), trda=trda  )
  eeg$what   <- 'eeg'
  
  muaR$exp <- paste( df$session[i], dxt, pxt, sep='')
  mua$exp  <- paste( df$session[i], dxt, pxt, sep='')
  lfp$exp  <- paste( df$session[i], dxt, pxt, sep='')
  eeg$exp  <- paste( df$session[i], dxt, pxt, sep='')

  muaR$dataSet <- 'g'
  mua$dataSet  <- 'g'
  lfp$dataSet  <- 'g'
  eeg$dataSet  <- 'g'
  
  mua$xpp$ISI[ which( mua$xpp$ISI<0)]  <- 20
  muaR$xpp$ISI[which(muaR$xpp$ISI<0)]  <- 20
  lfp$xpp$ISI[ which( lfp$xpp$ISI<0)]  <- 20
  eeg$xpp$ISI[ which( eeg$xpp$ISI<0)]  <- 20
  
  ## ========================================================= p2p mua
  tmp          <- apply(mua$lfp,c(1,3),max) - apply(mua$lfp,c(1,3),min)
  mua$xpp$p2p  <- apply( tmp, 1, max )
  
  qtmp         <- quantile(mua$xpp$p2p,probs=c(0.5,0.75))
  mxmua        <- qtmp[1]+10*diff(qtmp)
  valIndMUA    <- mua$xpp$p2p < mxmua;

  
  ## ========================================================= p2p lfp
  tmp          <- apply(lfp$lfp,c(1,3),max) - apply(lfp$lfp,c(1,3),min)
  lfp$xpp$p2p  <- apply( tmp, 1, max )
  
  qtmp         <- quantile(lfp$xpp$p2p,probs=c(0.5,0.75))
  mxlfp        <- qtmp[1]+10*diff(qtmp)
  valIndLFP    <- lfp$xpp$p2p < mxlfp;
  
  ## calculate single-trial csd
  lag     <- 1
  smooth  <- T
  sd      <- 100
  spacing <- df[i,]$spacing
  varCor  <- F
  
  csd      <- calcSingleTrialCSD(lfp, spacing=spacing, lag=lag,smooth=T, sd=sd)
  csd$what <- 'csd'
  
  # ----------------
  tmp          <- apply(csd$lfp,c(1,3),max) - apply(csd$lfp,c(1,3),min)
  csd$xpp$p2p  <- apply( tmp, 1, max )
  
  qtmp         <- quantile(csd$xpp$p2p,probs=c(0.5,0.75))
  mxcsd        <- qtmp[1]+10*diff(qtmp)
  valIndCSD    <- csd$xpp$p2p < mxcsd;

  ### ======================================================================================= set of analyses using a scalar variable such as mean mua activity in specific time-window
  if(doScalarMUA){
    ## =================================================================================================== analyses per session (mostly CSD)
    
    cmpStrDF     <- data.frame( what='csd', cmpStrName='CSD0',  tmin=8,  tmax=13,  dmin= -1500, dmax= -900, stringsAsFactors=FALSE) 
    cmpStrDF[2,] <- data.frame( what='csd', cmpStrName='CSD1',  tmin=15,  tmax=40,  dmin= -500,  dmax=  500, stringsAsFactors=FALSE)
    cmpStrDF[3,] <- data.frame( what='csd', cmpStrName='CSD2',  tmin=45,  tmax=55,  dmin= -500,  dmax=  500, stringsAsFactors=FALSE)
    cmpStrDF[4,] <- data.frame( what='csd', cmpStrName='CSD3',  tmin=65,  tmax=110, dmin= -500,  dmax=  500, stringsAsFactors=FALSE)
    cmpStrDF[5,] <- data.frame( what='csd', cmpStrName='CSD4',  tmin=100, tmax=150,  dmin= -500, dmax=  500, stringsAsFactors=FALSE)
    cmpStrDF[6,] <- data.frame( what='csd', cmpStrName='CSD5',  tmin=140, tmax=240,  dmin= -500, dmax=  500, stringsAsFactors=FALSE)
    
    #cmpStrDF <- cmpStrDF[c(2,3),]
    
    for (cx in 1:dim(cmpStrDF)){
      #print(cx)
      if (identical(cmpStrDF$what[cx],'mua') )
        dat <- mua
      if (identical(cmpStrDF$what[cx],'muaR') )
        dat <- muaR
      if (identical(cmpStrDF$what[cx],'csd') )
        dat <- csd
      if (identical(cmpStrDF$what[cx],'lfp') )
        dat <- lfp
      
      #dat$xpp$pitchInd       <- dat$xpp$toneInd #(dat$xpp$pitch-1)%%11
      #dat$xpp$pitchIndFactor <- as.factor(dat$xpp$pitchInd)
      
      
      tnd        <- which(dat$taxis >= cmpStrDF$tmin[cx] & dat$taxis < cmpStrDF$tmax[cx])
      
      dpthrng       <- c( cmpStrDF[cx,c('dmin','dmax')] )
      chDepth       <- (df$Layer[i] - 1:df$numChan[i]) * df$spacing[i]
      chnd          <- which(chDepth>=dpthrng[1] & chDepth<=dpthrng[2])
      
      if(length(chnd)>0){# channels present in depth range?
        DVstr         <- cmpStrDF$cmpStrName[cx]
        
        dat$xpp$TMP   <- apply(dat$lfp[,tnd,chnd], c(1), mean)
        if (identical(cmpStrDF$what[cx],'csd'))
          dat$xpp$TMP   <- apply(abs(dat$lfp[,tnd,chnd]), c(1), mean)
        
        dat$xpp$DV    <- dat$xpp$TMP 
        
        nd                  <- which(names(dat$xpp)=='TMP')
        names(dat$xpp)[nd]  <- DVstr
        
        dat$xpp$ampIndFactor <- as.factor(dat$xpp$ampInd)
        
        # --------------- run analyzeRSscalar
        valBin <- getValidTrialsAOC(dat, who='all')
        disc   <- analyzeAOCscalar(dat, ybin=valBin, who='all', cmpStr=DVstr)
        
        valBin <- getValidTrialsAOC(dat, who='regular')
        disc   <- analyzeAOCscalar(dat, ybin=valBin, who='regular', cmpStr=DVstr)
        
        valBin <- getValidTrialsAOC(dat, who='random')
        disc   <- analyzeAOCscalar(dat, ybin=valBin, who='random', cmpStr=DVstr)
      }# if channels present in depth range
    }# end cmpStrList loop
    
    
    ## ==================================================================================================== analyses by channel (mostly MUA)
    cmpStrDF     <- data.frame( what='mua', cmpStrName='MUA0',  tmin=8,  tmax=13,  stringsAsFactors=FALSE) 
    cmpStrDF[2,] <- data.frame( what='mua', cmpStrName='MUA1',  tmin=18,  tmax=35,  stringsAsFactors=FALSE)
    cmpStrDF[3,] <- data.frame( what='mua', cmpStrName='MUA1a', tmin=15,  tmax=25,  stringsAsFactors=FALSE)
    cmpStrDF[4,] <- data.frame( what='mua', cmpStrName='MUA1b', tmin=25,  tmax=35,  stringsAsFactors=FALSE)
    cmpStrDF[5,] <- data.frame( what='mua', cmpStrName='MUA1c', tmin=35,  tmax=45,  stringsAsFactors=FALSE)
    cmpStrDF[6,] <- data.frame( what='mua', cmpStrName='MUA2' , tmin=55,  tmax=125, stringsAsFactors=FALSE)
    cmpStrDF[7,] <- data.frame( what='mua', cmpStrName='MUA3' , tmin=160, tmax=250, stringsAsFactors=FALSE)
    
    cmpStrDF[8,]  <- data.frame( what='csd', cmpStrName='CSD0',  tmin=15,   tmax=19,  stringsAsFactors=FALSE)
    cmpStrDF[9,]  <- data.frame( what='csd', cmpStrName='CSD1',  tmin=25,   tmax=45,  stringsAsFactors=FALSE)
    cmpStrDF[10,] <- data.frame( what='csd', cmpStrName='CSD2',  tmin=50,   tmax=60,  stringsAsFactors=FALSE)
    cmpStrDF[11,] <- data.frame( what='csd', cmpStrName='CSD3',  tmin=65,   tmax=110, stringsAsFactors=FALSE)
    cmpStrDF[12,] <- data.frame( what='csd', cmpStrName='CSD4',  tmin=100,  tmax=150, stringsAsFactors=FALSE)
    cmpStrDF[13,] <- data.frame( what='csd', cmpStrName='CSD5',  tmin=160,  tmax=250, stringsAsFactors=FALSE)
    
    cmpStrDF[14,] <- data.frame( what='lfp', cmpStrName='LFP0',  tmin=10,  tmax=17, stringsAsFactors=FALSE)
    cmpStrDF[15,] <- data.frame( what='lfp', cmpStrName='LFP1',  tmin=20,  tmax=30, stringsAsFactors=FALSE)
    
    cmpStrDF[16,] <- data.frame( what='muaR', cmpStrName='MUAbl',  tmin= -50,  tmax=0,   stringsAsFactors=FALSE)
    cmpStrDF[17,] <- data.frame( what='muaR', cmpStrName='MUAR1',  tmin=18,    tmax=35,  stringsAsFactors=FALSE)
    cmpStrDF[18,] <- data.frame( what='muaR', cmpStrName='MUAR1a', tmin=15,    tmax=25,  stringsAsFactors=FALSE)
    cmpStrDF[19,] <- data.frame( what='muaR', cmpStrName='MUAR1b', tmin=25,    tmax=35,  stringsAsFactors=FALSE)
    cmpStrDF[20,] <- data.frame( what='muaR', cmpStrName='MUAR1c', tmin=35,    tmax=45,  stringsAsFactors=FALSE)
    cmpStrDF[21,] <- data.frame( what='muaR', cmpStrName='MUAR2' , tmin=55,    tmax=125, stringsAsFactors=FALSE)
    cmpStrDF[22,] <- data.frame( what='muaR', cmpStrName='MUAR3' , tmin=160,   tmax=250, stringsAsFactors=FALSE)
    
    #cmpStrDF <- cmpStrDF[c(2,6),]#,16,17,21),]
    
    for (chnd in 1:df$numChan[i]){
      for (cx in 1:dim(cmpStrDF)[1]){
        
        if (identical(cmpStrDF$what[cx],'mua') )
          dat <- mua
        if (identical(cmpStrDF$what[cx],'muaR') )
          dat <- muaR
        if (identical(cmpStrDF$what[cx],'csd') )
          dat <- csd
        if (identical(cmpStrDF$what[cx],'lfp') )
          dat <- lfp
        
        tnd        <- which(dat$taxis >= cmpStrDF$tmin[cx] & dat$taxis < cmpStrDF$tmax[cx])
        
        DVstr         <- paste('ch', df$startVIndex[i]+chnd, '_', cmpStrDF$cmpStrName[cx], sep='')
        dat$xpp$TMP   <- apply(dat$lfp[,tnd,chnd], c(1), mean)
        dat$xpp$DV    <- dat$xpp$TMP 
        
        nd                  <- which(names(dat$xpp)=='TMP')
        names(dat$xpp)[nd]  <- DVstr
        
        dat$xpp$pitchInd       <- dat$xpp$toneInd   #(dat$xpp$pitch-1)%%11
        if (identical( cmpStrDF$cmpStrName,'MUAbl') ){   # shift pitchInd to reflect pitch of last tone
          dat$xpp$ampInd      <- my.lag( dat$xpp$ampInd) #                 (dat$xpp$pitch-1)%%11, 1)
          dat$xpp$deltaAmpInd <- my.lag( dat$xpp$deltaAmpInd,1)
        }
        
        
        if (identical(cmpStrDF$what[cx],'muaR')){
          fileName <- paste( rda(), df$session[i], '/MUAbl_brkpoints.rda',sep='' )
          
          if (file.exists(fileName)){# use previously saved break points to remove sharp artifacts
            print(fileName)
            load(file=fileName)# load variable 'brks'
            
            tb   <- unique( c(0,brks,length(dat$xpp[[DVstr]])))
            x    <- 1:length(dat$xpp[[DVstr]])
            xbin <- cut(x, tb )
            
            lm0 <- lm(dat$xpp[[DVstr]] ~ 1)
            lm1 <- lm(dat$xpp[[DVstr]] ~ 1 + xbin)
            pVal <- anova(lm0,lm1)$P[2]           
            
            if( pVal<0.001 )
              dat$xpp[[DVstr]] <- lm1$residuals
          }
          
          # correct for slow fluctuations
          if (1==1){ #(df$corMUAbl[i]==2){
            print('correcting for slow fluctuations')
            mavDV            <- get.mav(dat$xpp[[DVstr]], 1:dim(dat$xpp)[1], eval.at=seq(1,dim(dat$xpp)[1],by=1 ), show=T, distribution='Gaussian',boot=F, width=100)
            dat$xpp[[DVstr]] <- dat$xpp[[DVstr]] - mavDV$mav + mean(mavDV$mav)
          }
        }
        
        
        dat$xpp$ampIndFactor <- as.factor(dat$xpp$ampInd)
        
        # --------------- run analyzeRSscalar
        valBin <- getValidTrialsAOC(dat, who='all',DVstr=DVstr)
        disc   <- analyzeAOCscalar(dat, ybin=valBin, who='all', cmpStr=DVstr)
        
        valBin <- getValidTrialsAOC(dat, who='regular')
        disc   <- analyzeAOCscalar(dat, ybin=valBin, who='regular', cmpStr=DVstr)
        
        valBin <- getValidTrialsAOC(dat, who='random')
        disc   <- analyzeAOCscalar(dat, ybin=valBin, who='random', cmpStr=DVstr)
      }# end cmpStrList loop
    }# end channel loop
    
  }# end RSscalar
  
  
  
  
  
  if(doSessionAverage){
    #browser()
    
    fltxtn <- ''
   
    getSessionAverage(lfp, maxp2p=mxlfp, who='all', fextn=fltxtn)
    gc()
    
    getSessionAverage(eeg, maxp2p=450, who='all', fextn=fltxtn)
    gc()
    
     getSessionAverage(csd, maxp2p=mxcsd, who='all', fextn=fltxtn)
    gc()
    
    getSessionAverage(mua, maxp2p=mxmua, who='all', fextn=fltxtn)
    gc()
    
    ## random trialType==1
    vnd <- which(lfp$xpp$trialType==1)
    if (length(vnd)>100){

      getSessionAverage(subs.data(lfp,lfp$xpp$trialType==1), maxp2p=mxlfp, who='random', fextn=fltxtn)
      gc()
      
      getSessionAverage(subs.data(eeg,eeg$xpp$trialType==1), maxp2p=450, who='random', fextn=fltxtn)
      gc()
      
      getSessionAverage(subs.data(csd,csd$xpp$trialType==1), maxp2p=mxcsd, who='random', fextn=fltxtn)
      gc()
      
      getSessionAverage(subs.data(mua,mua$xpp$trialType==1), maxp2p=mxmua, who='random', fextn=fltxtn)
      gc()
    }

    ## regular trialType==2
    vnd <- which(lfp$xpp$trialType==2)
    if (length(vnd)>100){
      getSessionAverage(subs.data(lfp,lfp$xpp$trialType==2), maxp2p=mxlfp, who='regular', fextn=fltxtn)
      gc()
      
      getSessionAverage(subs.data(eeg,eeg$xpp$trialType==2), maxp2p=450, who='regular', fextn=fltxtn)
      gc()
      
      getSessionAverage(subs.data(csd,csd$xpp$trialType==2), maxp2p=mxcsd, who='regular', fextn=fltxtn)
      gc()
      
      getSessionAverage(subs.data(mua,mua$xpp$trialType==2), maxp2p=mxmua, who='regular', fextn=fltxtn)
      gc()
      }
  }
  
  ## =======================================  Grand Average CSD
  lfpscl  <- 4
  csdscl  <- 300#500#1300
  
  csdGA   <- calcCSD(lfp, mua, lag=lag, mxmua=mxmua,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sd,varCor=varCor,spacing=spacing)
  
  eps.file('GA_CSDImage.eps',pdf=T,s=5.80)
  disc <- showCSDImage(csdGA,trng=c(-10,250),lfpscl=lfpscl,csdscl=csdscl,tat=NA,addcsd=T)
  dev.off()
  
  if(CSDplots){
    print('CSDplots')
    ## =========================================  CSD by SOA
    isicsdList <- list()
    for (j in 1:6){
      subs <- as.numeric(lfp$xpp$isiOct)==j 
      subs[which(is.na(subs))] <- F
      isicsdList[[j]] <- calcCSD(subs.data(lfp,subs),subs.data(mua,subs),lag=lag, mxmua=mxmua,mxlfp=mxlfp, 
                                 interpChanInd==interpChanInd,smooth=smooth,sd=sd,correction=csdGA$correction,spacing=spacing )
    }
    
    ## all in one image
    eps.file('CSDbySOA.eps',pdf=T,s=5.80)
    disc <- showCSD(isicsdList[[1]],lfpscl=20, csdscl=.1, muascl=.75,trng=c(-10,250),addpoly=F, col=colList[1],add=F)
    for (j in 2:6)
      disc <- showCSD(isicsdList[[j]],lfpscl=20, csdscl=.1, muascl=.75,trng=c(-10,250),addpoly=F, col=colList[j],add=T)
    dev.off()
    
    # csdImage for each SOA separately
    for (j in 1:6){
      eps.file(paste('CSDImageBySOA_',j,'.eps',sep=''),pdf=T,s=5.80)
      disc <- showCSDImage(isicsdList[[j]],lfpscl=lfpscl, csdscl=csdscl, muascl=100,trng=c(-10,250), addcsd=T)
      dev.off()
    }
    
    ## AVREC
    colList <- c('black','blue','cyan','green','orange','red')
    eps.file('AVRECBySOA.eps',pdf=T,s=2.5,ratio=3/2)
    disc <- plot(isicsdList[[6]]$txlfp, isicsdList[[6]]$avrec,xlim=c(-25,250), col=colList[1],type='l',lwd=2)
    for (j in 1:6)
      disc <- points(isicsdList[[1]]$txlfp,isicsdList[[j]]$avrec, col=colList[j],type='l',lwd=2)
    dev.off()
    
    
    ## =========================================  CSD by Pitch
    pitchcsdList <- list()
    for (j in 1:5){
      subs            <- as.numeric(lfp$xpp$ampInd)==j 
      subs[which(is.na(subs))] <- F
      pitchcsdList[[j]] <- calcCSD(subs.data(lfp,subs),subs.data(mua,subs),lag=lag, mxmua=mxmua,mxlfp=mxlfp, 
                                   interpChanInd==interpChanInd,smooth=smooth,sd=sd,correction=csdGA$correction,spacing=spacing )
    }
    
    ## all in one image
    eps.file('CSDbyAmp.eps',pdf=T,s=5.80)
    disc <- showCSD(pitchcsdList[[1]],lfpscl=20, csdscl=.1, muascl=.75,trng=c(-10,250),addpoly=F, col=colList[1],add=F)
    for (j in 2:5)
      disc <- showCSD(pitchcsdList[[j]],lfpscl=20, csdscl=.1, muascl=.75,trng=c(-10,250),addpoly=F, col=colList[j],add=T)
    dev.off()
    
    # csdImage for each SOA separately
    for (j in 1:5){
      eps.file(paste('CSDImageByAmp_',j,'.eps',sep=''),pdf=T,s=5.80)
      disc <- showCSDImage(pitchcsdList[[j]],lfpscl=lfpscl, csdscl=csdscl, muascl=100,trng=c(-10,250),addcsd=T)
      dev.off()
    }
  
    
    ## AVREC
    colList <- c('black','blue','cyan','green','orange','red')
    eps.file('AVRECByAmp.eps',pdf=T,s=2.5,ratio=3/2)
    disc <- plot(pitchcsdList[[5]]$txlfp, pitchcsdList[[5]]$avrec,xlim=c(-25,250), col=colList[1],type='l',lwd=2)
    for (j in 1:5)
      disc <- points(pitchcsdList[[1]]$txlfp,pitchcsdList[[j]]$avrec, col=colList[j],type='l',lwd=2)
    dev.off()
  }
  
  if (doTrLM){
    print('time-resolved linear model')
    #browser()
    mua$xpp$isiLog  <- log(mua$xpp$ISI)
    muaR$xpp$isiLog <- log(muaR$xpp$ISI)
    csd$xpp$isiLog  <- log(csd$xpp$ISI)
    lfp$xpp$isiLog  <- log(lfp$xpp$ISI)
    eeg$xpp$isiLog  <- log(eeg$xpp$ISI)
    
    mua$xpp$absDeltaAmpInd  <- abs( mua$xpp$deltaAmpInd )
    muaR$xpp$absDeltaAmpInd <- abs( muaR$xpp$deltaAmpInd )
    csd$xpp$absDeltaAmpInd  <- abs( csd$xpp$deltaAmpInd )
    lfp$xpp$absDeltaAmpInd  <- abs( lfp$xpp$deltaAmpInd )
    eeg$xpp$absDeltaAmpInd  <- abs( eeg$xpp$deltaAmpInd )
    
    variables  <- c('sp','isiLog','ampInd','deltaAmpInd', 'absDeltaAmpInd','repNr','trialType','trialNr')
    variables2 <- c('sp','isiLog','ampInd','deltaAmpInd', 'absDeltaAmpInd','repNr',            'trialNr')

    ## run time-resolved general linear model 
    for (chNr in 2:dim(eeg$lfp)[3]){
      print(chNr)
      
      eeg$ch <- chNr
      disc   <- trStats.ampOddClick(eeg,width=.050,eval.at=eeg$taxis,variables=variables,what='lfp',chNr=chNr,
                                 subs=valIndCSD,extn='_eeg_all',show=F, fix=FALSE)
      disc   <- trStats.ampOddClick(eeg,width=.050,eval.at=eeg$taxis,variables=variables2,what='lfp',chNr=chNr,
                                 subs=valIndCSD&eeg$xpp$trialType==1,extn='_eeg_random',show=F, fix=FALSE)
      disc   <- trStats.ampOddClick(eeg,width=.050,eval.at=eeg$taxis,variables=variables2,what='lfp',chNr=chNr,
                                 subs=valIndCSD&eeg$xpp$trialType==2,extn='_eeg_regular',show=F, fix=FALSE)
      
    }
    
    for (chNr in 1:df$numChan[i]){
      print(chNr)
      mua$ch <- chNr
      disc   <- trStats.ampOddClick(mua,width=.050,eval.at=mua$taxis,variables=variables,what='lfp',chNr=chNr,
                                 subs=valIndMUA,extn='_mua_all',show=F, fix=FALSE)
      disc   <- trStats.ampOddClick(mua,width=.050,eval.at=mua$taxis,variables=variables2,what='lfp',chNr=chNr,
                                 subs=valIndMUA&mua$xpp$trialType==1,extn='_mua_random',show=F, fix=FALSE)
      disc   <- trStats.ampOddClick(mua,width=.050,eval.at=mua$taxis,variables=variables2,what='lfp',chNr=chNr,
                                 subs=valIndMUA&mua$xpp$trialType==2,extn='_mua_regular',show=F, fix=FALSE)
      
      muaR$ch <- chNr
      disc    <- trStats.ampOddClick(muaR,width=.050,eval.at=mua$taxis,variables=variables,what='lfp',chNr=chNr,
                                  subs=valIndMUA,extn='_muaR_all',show=F, fix=FALSE)
      disc   <- trStats.ampOddClick(muaR,width=.050,eval.at=muaR$taxis,variables=variables2,what='lfp',chNr=chNr,
                                 subs=valIndMUA&muaR$xpp$trialType==1,extn='_muaR_random',show=F, fix=FALSE)
      disc   <- trStats.ampOddClick(muaR,width=.050,eval.at=muaR$taxis,variables=variables2,what='lfp',chNr=chNr,
                                 subs=valIndMUA&muaR$xpp$trialType==2,extn='_muaR_regular',show=F, fix=FALSE)
      
      Nchan   <- dim(mua$lfp)[3]
      csd$ch <- chNr
      if(chNr>1 & chNr<Nchan){
        disc   <- trStats.ampOddClick(csd,width=.050,eval.at=csd$taxis,variables=variables,what='lfp',chNr=chNr,
                                   subs=valIndCSD,extn='_csd_all',show=F, fix=FALSE)
        disc   <- trStats.ampOddClick(csd,width=.050,eval.at=csd$taxis,variables=variables2,what='lfp',chNr=chNr,
                                   subs=valIndCSD&csd$xpp$trialType==1,extn='_csd_random',show=F, fix=FALSE)
        disc   <- trStats.ampOddClick(csd,width=.050,eval.at=csd$taxis,variables=variables2,what='lfp',chNr=chNr,
                                   subs=valIndCSD&csd$xpp$trialType==2,extn='_csd_regular',show=F, fix=FALSE)
      }
    }# chNr
  }
  
  
  # this is being used for STPSPbyChan and STPSP
  trng    <- c(20,70)
  blrng   <- c(-50,0)
  minISI  <- 0.190
  maxISI  <- 12.8
  DVscl   <- 5
  sfc     <- 1
  Nchan   <- dim(mua$lfp)[3]
  muach   <- 1:Nchan
  
  
  if (doSTPSPbyChan){
    print('do STPSP by Channel')
    
    ## ================================================================================= analyses by channel (mostly MUA)
    cmpStrDF     <- data.frame( what='mua', cmpStrName='MUA0',  tmin=8 ,  tmax=13,  stringsAsFactors=FALSE) 
    cmpStrDF[2,] <- data.frame( what='mua', cmpStrName='MUA1',  tmin=18,  tmax=35,  stringsAsFactors=FALSE)
    cmpStrDF[3,] <- data.frame( what='mua', cmpStrName='MUA1a', tmin=15,  tmax=25,  stringsAsFactors=FALSE)
    cmpStrDF[4,] <- data.frame( what='mua', cmpStrName='MUA1b', tmin=25,  tmax=35,  stringsAsFactors=FALSE)
    cmpStrDF[5,] <- data.frame( what='mua', cmpStrName='MUA1c', tmin=35,  tmax=45,  stringsAsFactors=FALSE)
    cmpStrDF[6,] <- data.frame( what='mua', cmpStrName='MUA2' , tmin=55,  tmax=125, stringsAsFactors=FALSE)
    cmpStrDF[7,] <- data.frame( what='mua', cmpStrName='MUA3' , tmin=160, tmax=250, stringsAsFactors=FALSE)
    
    cmpStrDF[8,]  <- data.frame( what='csd', cmpStrName='CSD0',  tmin=15,   tmax=19,  stringsAsFactors=FALSE)
    cmpStrDF[9,]  <- data.frame( what='csd', cmpStrName='CSD1',  tmin=25,   tmax=45,  stringsAsFactors=FALSE)
    cmpStrDF[10,] <- data.frame( what='csd', cmpStrName='CSD2',  tmin=50,   tmax=60,  stringsAsFactors=FALSE)
    cmpStrDF[11,] <- data.frame( what='csd', cmpStrName='CSD3',  tmin=65,   tmax=110, stringsAsFactors=FALSE)
    cmpStrDF[12,] <- data.frame( what='csd', cmpStrName='CSD4',  tmin=100,  tmax=150, stringsAsFactors=FALSE)
    cmpStrDF[13,] <- data.frame( what='csd', cmpStrName='CSD5',  tmin=160,  tmax=250, stringsAsFactors=FALSE)
    
    cmpStrDF[14,] <- data.frame( what='lfp', cmpStrName='LFP0',  tmin=8,  tmax=13, stringsAsFactors=FALSE)
    cmpStrDF[15,] <- data.frame( what='lfp', cmpStrName='LFP1',  tmin=20,  tmax=30, stringsAsFactors=FALSE)
    
    cmpStrDF[16,] <- data.frame( what='muaR', cmpStrName='MUAbl',  tmin= -50,  tmax=0,  stringsAsFactors=FALSE)
    cmpStrDF[17,] <- data.frame( what='muaR', cmpStrName='MUAR1',  tmin=18,    tmax=35, stringsAsFactors=FALSE)
    cmpStrDF[18,] <- data.frame( what='muaR', cmpStrName='MUAR1a', tmin=15,  tmax=25,  stringsAsFactors=FALSE)
    cmpStrDF[19,] <- data.frame( what='muaR', cmpStrName='MUAR1b', tmin=25,  tmax=35,  stringsAsFactors=FALSE)
    cmpStrDF[20,] <- data.frame( what='muaR', cmpStrName='MUAR1c', tmin=35,  tmax=45,  stringsAsFactors=FALSE)
    cmpStrDF[21,] <- data.frame( what='muaR', cmpStrName='MUAR2' , tmin=55,  tmax=125, stringsAsFactors=FALSE)
    cmpStrDF[22,] <- data.frame( what='muaR', cmpStrName='MUAR3' , tmin=160, tmax=250, stringsAsFactors=FALSE)
    
    cmpStrDF[23,] <- data.frame( what='mua',  cmpStrName='MUA'  , tmin=20, tmax=70, stringsAsFactors=FALSE)
    cmpStrDF[24,] <- data.frame( what='muaR', cmpStrName='MUAR' , tmin=20, tmax=70, stringsAsFactors=FALSE)
    
    #cmpStrDF <- cmpStrDF[c(16,17,21),]   #cmpStrDF[c(1,2,6,9,15,17,21),]
    
    for (cx in 1:dim(cmpStrDF)[1]){
      
      if (identical(cmpStrDF$what[cx],'mua') )
        dat <- mua
      if (identical(cmpStrDF$what[cx],'muaR') )
        dat <- muaR
      if (identical(cmpStrDF$what[cx],'lfp') )
        dat <- lfp
      if (identical(cmpStrDF$what[cx],'csd') ){
        dat <- csd
        Nchan <- dim(dat$lfp)[3]
        dat$lfp[,,1]      <- dat$lfp[,,2]   
        dat$lfp[,,Nchan]  <- dat$lfp[,,Nchan-1]# avoid model crashes
      }
      
      if(identical(cmpStrDF$cmpStrName[cx],'MUAbl')  ){
        dat$xpp$toneInd          <- my.lag(muaR$xpp$toneInd,1)
        #  dat$xpp$ISI              <- my.lag(mua$xpp$ISI,1)
        #  dat$xpp$ISIOct           <- my.lag(mua$xpp$ISIOct,1)
      }
      
      ## ================= call STPSP
      parline        <- cmpStrDF[cx,]
      STPSPmua       <- prepSTPSPVProbeMUA(dat, minISI=minISI,maxISI=maxISI,DVscl=DVscl,sfc=sfc,who='all'    , parline=parline, dfline=df[i,])
      STPSPmuaRnd    <- prepSTPSPVProbeMUA(dat, minISI=minISI,maxISI=maxISI,DVscl=DVscl,sfc=sfc,who='random' , parline=parline, dfline=df[i,])# crashes for LFP ch32 regular
      STPSPmuaReg    <- prepSTPSPVProbeMUA(dat, minISI=minISI,maxISI=maxISI,DVscl=DVscl,sfc=sfc,who='regular', parline=parline, dfline=df[i,])
    } 
  }
  
  
  
  if (doSTPSP){# not yet adjusted to ampOddClick
    
    print('computing STPSP model fits')
    ## =========================================== STPSP model for 
    ## =========================================== MUA and STPSP
    #muach           <- 1:24
    muaind          <- which( mua$taxis>=trng[1] & mua$taxis<=trng[2] )
    mua$xpp$mua     <- apply( mua$lfp[,muaind,muach],1,mean)
    
    # prepare baseline data
    muaR$xpp$toneInd          <- my.lag(mua$xpp$toneInd,1)
    muaR$xpp$deltaPitchOct    <- my.lag(mua$xpp$deltaPitchOct,1)
    muaR$xpp$absDeltaPitchOct <- my.lag(mua$xpp$absDeltaPitchOct,1)
    
    
    muaind             <- which( mua$taxis>=blrng[1] & mua$taxis<=blrng[2] )
    muaR$xpp$muabl     <- apply( muaR$lfp[,muaind,muach],1,mean)
    
    thisDir <- ''    
    #thisDir   <- paste(mua$exp,sep='')
    tmpDir <- paste(pp(),thisDir,'/STPSP',sep='' )
    if(!exists(tmpDir))
      system(paste('mkdir',tmpDir) )
    
    DVstr      <- 'mua'
    STPSPfit   <- callSTPSP(mua, gd=T, startFromScratch=sfc, DVstr=DVstr,who='all', minFrq=500, toneBy=1,NOct=3,showCh=1,
                            minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=T, eps=T, fileExtn='allch',baseName=thisDir)
    print( c(log(2)/STPSPfit$modelFit$par[1],STPSPfit$model$par[2]) )
    
    STPSPfitRnd   <- callSTPSP(mua, gd=T, startFromScratch=sfc, DVstr=DVstr,who='random', minFrq=500, toneBy=1,NOct=3,showCh=1,
                               minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=T, eps=T, fileExtn='allch',baseName=thisDir)
    print( c(log(2)/STPSPfitRnd$modelFit$par[1],STPSPfitRnd$model$par[2]) )
    
    STPSPfitReg   <- callSTPSP(mua, gd=T, startFromScratch=sfc, DVstr=DVstr,who='regular', minFrq=500, toneBy=1,NOct=3,showCh=1,
                               minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=T, eps=T, fileExtn='allch',baseName=thisDir)
    print( c(log(2)/STPSPfitReg$modelFit$par[1],STPSPfitReg$model$par[2]) )
    
    
    
    STPSPfitFrq   <- callFrqSTPSP(mua, gd=T, startFromScratch=sfc,predFUN=STPSPmodelFun_freq, DVstr=DVstr,who='all',
                                  minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=T, eps=T, fileExtn='allch_frq3par',baseName=thisDir)
    
    STPSPfitFrqRnd   <- callFrqSTPSP(mua, gd=T, startFromScratch=sfc,predFUN=STPSPmodelFun_freq, DVstr=DVstr,who='random', 
                                     minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=T, eps=T, fileExtn='allch_frq3par',baseName=thisDir)
    
    STPSPfitFrqRnd   <- callFrqSTPSP(mua, gd=T, startFromScratch=sfc,predFUN=STPSPmodelFun_freq, DVstr=DVstr,who='regular',
                                     minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=T, eps=T, fileExtn='allch_frq3par',baseName=thisDir)
    
    
    
    ## baseline 
    DVstr      <- 'muabl'
    STPSPfit   <- callSTPSP(muaR, gd=T, startFromScratch=sfc, DVstr=DVstr,who='all', minFrq=500, toneBy=1,NOct=3,showCh=1,
                            minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=T, eps=T, fileExtn='allch',baseName=thisDir)
    print( c(log(2)/STPSPfit$modelFit$par[1],STPSPfit$model$par[2]) )
    
    STPSPfitRnd   <- callSTPSP(muaR, gd=T, startFromScratch=sfc, DVstr=DVstr,who='random', minFrq=500, toneBy=1,NOct=3,showCh=1,
                               minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=T, eps=T, fileExtn='allch',baseName=thisDir)
    print( c(log(2)/STPSPfitRnd$modelFit$par[1],STPSPfitRnd$model$par[2]) )
    
    STPSPfitReg   <- callSTPSP(muaR, gd=T, startFromScratch=sfc, DVstr=DVstr,who='regular', minFrq=500, toneBy=1,NOct=3,showCh=1,
                               minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=T, eps=T, fileExtn='allch',baseName=thisDir)
    print( c(log(2)/STPSPfitReg$modelFit$par[1],STPSPfitReg$model$par[2]) )
    
    
    
    STPSPfitFrq   <- callFrqSTPSP(muaR, gd=T, startFromScratch=sfc,predFUN=STPSPmodelFun_freq, DVstr=DVstr,who='all',
                                  minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=T, eps=T, fileExtn='allch_frq3par',baseName=thisDir)
    
    STPSPfitFrqRnd   <- callFrqSTPSP(muaR, gd=T, startFromScratch=sfc,predFUN=STPSPmodelFun_freq, DVstr=DVstr,who='random', 
                                     minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=T, eps=T, fileExtn='allch_frq3par',baseName=thisDir)
    
    STPSPfitFrqRnd   <- callFrqSTPSP(muaR, gd=T, startFromScratch=sfc,predFUN=STPSPmodelFun_freq, DVstr=DVstr,who='regular',
                                     minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=T, eps=T, fileExtn='allch_frq3par',baseName=thisDir)
    
    
    
    eps.file(paste('GA_MUAbySOA.eps',sep=''),pdf=T,s=5.80)
    soax   <- log2(tapply(mua$xpp$ISI, mua$xpp$isiOctFine, mean))
    
    N      <- table(mua$xpp$isiOctFine)
    valInd <- which(as.numeric(N)>10)
    
    mnMUA <- tapply( mua$xpp$mua[valIndMUA], list(mua$xpp$isiOctFine[valIndMUA],mua$xpp$trialType[valIndMUA]), mean)
    NMUA  <- tapply( mua$xpp$mua[valIndMUA], list(mua$xpp$isiOctFine[valIndMUA],mua$xpp$trialType[valIndMUA]), length)
    ylim  <- c(0.98*na.min(mnMUA), 1.02*na.max(c(mnMUA)) ) 
    
    plot(   soax, tapply( mua$xpp$mua[valIndMUA], mua$xpp$isiOctFine[valIndMUA], mean), type='l', col='black',lwd=2,ylim=ylim)
    addsebar <- function(n){sebar(mua$xpp$mua[valIndMUA & as.numeric(mua$xpp$isiOctFine)==n],at=soax[n],lwd=2,col='black') }
    disc     <- sapply(1:length(N),addsebar)
    
    points(   soax, mnMUA[,1], type='l', col='blue',lwd=2)
    addsebar <- function(n,m=1){sebar(mua$xpp$mua[valIndMUA & as.numeric(mua$xpp$isiOctFine)==n&mua$xpp$trialType==m],at=soax[n],lwd=2,col='blue') }
    disc     <- sapply(1:length(N),addsebar)
    
    points(   soax, mnMUA[,2], type='l', col='red',lwd=2)
    addsebar <- function(n,m=2){sebar(mua$xpp$mua[valIndMUA & as.numeric(mua$xpp$isiOctFine)==n&mua$xpp$trialType==m],at=soax[n],lwd=2,col='red') }
    disc     <- sapply(1:length(N),addsebar)
    dev.off()
    
    
    ## ===============================================================
    ## do diff LFP above and below putative layer 2/3 source/sink pair (not present in all animals, but use default anyway)
    if(df$Layer[i]==0)
      df$Layer[i] <- 10
    
    tipChan <- (ceiling( df$Layer[i] )+2) %<=% 24
    topChan <- (floor(   df$Layer[i] )-2) %>=% 1
    diff    <- lfp$lfp[,,topChan] - lfp$lfp[,,tipChan]  
    #plot(lfp$taxis,apply(diff,2,mean),type='l',xlim=c(-10,100) )
    p35rng  <- c(20,50)
    p35ind  <- which(lfp$taxis>=p35rng[1] & lfp$taxis<=p35rng[2])
    mua$xpp$p35 <- apply(diff[,p35ind],1,mean)
    
    
    minISI     <- .275   # seems to be enhanced for ISIs below 275 ms
    DVstr      <- 'p35'
    STPSPfit   <- callSTPSP(mua, gd=T, startFromScratch=sfc, DVstr=DVstr,who='all', minFrq=500, toneBy=1,NOct=3,showCh=1,
                            minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=F, eps=T, fileExtn='')
    
    STPSPfitRnd   <- callSTPSP(mua, gd=T, startFromScratch=sfc, DVstr=DVstr,who='random', minFrq=500, toneBy=1,NOct=3,showCh=1,
                               minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=F, eps=F, fileExtn='')
    
    STPSPfitReg   <- callSTPSP(mua, gd=T, startFromScratch=sfc, DVstr=DVstr,who='regular', minFrq=500, toneBy=1,NOct=3,showCh=1,
                               minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=F, eps=F, fileExtn='')
    
    print( c(log(2)/STPSPfit$modelFit$par[1],STPSPfit$model$par[2]) )
    print( c(log(2)/STPSPfitRnd$modelFit$par[1],STPSPfitRnd$model$par[2]) )
    print( c(log(2)/STPSPfitReg$modelFit$par[1],STPSPfitReg$model$par[2]) )
  }
}