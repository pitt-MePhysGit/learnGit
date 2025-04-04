##20211021 - SK: made most of the csd calculations conditional on a flag; saves a lot of time

FLA_VProbe_STPSP <- function(j,df=df,baseDir=baseDir,trda=trda,ind=ind){
  
  ## PRELIMINARIES                ====
  i <- ind[j]
  print(df$session[i])
  
  project <- df$project[i]#'ampOddClick'
  baseName <- paste(df$session[i],'/',sep='')

  ## SET PATHS                    ====
    
  #pp   <- function(){paste(baseDir, 'RSkernel/word/VProbe/plots/',sep='')}
  #.GlobalEnv$rda  <- function(){'~/Dropbox/RSkernel/Rdata/'}
  #.GlobalEnv$rda  <- function(){ return(paste(trda(),df$session[i],'rda/',sep='/')) }
  .GlobalEnv$rda  <- trda
  
  if (!file.exists(rda()))
    system( paste('mkdir', rda() ))
  
  dxt <- ifelse(df$drug[i]=='N','', paste('_',df$drug[i],sep=''))
  pxt <- ifelse(df$probe[i]==1 ,'', paste('_',df$probe[i],sep=''))
  
  .GlobalEnv$pp   <- function(){paste(baseDir, 'ampOddClick/word/VProbe/plots/',df$session[i],dxt,pxt,'/',sep='')}
  print(pp())
  
  interpChanInd <- c();#c(16)
  
  if (!file.exists(pp()))
    system( paste('mkdir', pp()) )
  
  
  ## LOAD LIBRARIES               =================================
  library(R.matlab)
  library(rgenoud)
  library(lattice)
  library(akima)
  library(Rwave)
  library(MASS)
  library(fastICA)
  
  #source(paste(baseDir, 'RSkernel/r-files/functions/source_all_functions_parallel.R',sep=''))
  source( paste(baseDir,'ampOddClick/r-files/functions/functions_ampOddclick.R',sep='') )
  source( paste(baseDir,'ampOddClick/r-files/functions/functions_AOCdual.R',sep='') )
  source( paste(baseDir,'ampOddClick/r-files/functions/functions_VProbe_STPSP.R',sep='') )
  #source( paste(baseDir,'RSkernel/r-files/functions/FLA_VProbe_STPSP.R',sep='') )
  
  source( paste(baseDir,'RSkernel/r-files/functions/removeCommonReferenceNoise_ica.R',sep='') )
  
  source( paste(baseDir,'r-compile/rUtils/R/rUtils.R',sep='') )
  source( paste(baseDir,'r-compile/mav.test/R/mav.test.R',sep='') )
  source( paste(baseDir,'r-compile/ePhysLab/R/ePhysLab.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/calcCSD.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/showCSDimage.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/showCSD.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/multiVplot.R',sep='') )
  
  source( paste(baseDir,'ampOddClick/r-files/functions/prepSTPSPVProbeMUA.R',sep='') )
  source( paste(baseDir,'ampOddClick/r-files/functions/functions_ampOdd_presynapticPlasticity.R',sep='') )
  
  ##Need OKEEG to stay FALSE in line ~846 because this sourcing requires a package that I can't install on tinyone
  #  source( paste(baseDir,'ampOddClick/r-files/functions/trStats_ampOddClick.R',sep='') )
  
  source( paste(baseDir,'ampOddClick/r-files/functions/getSessionAverage.R',sep='') )
  source( paste(baseDir,'ampOddClick/r-files/functions/getValidTrialsAOC.R',sep='') )
  source( paste(baseDir,'ampOddClick/r-files/functions/analyzeAOCscalar.R',sep='') )
  source( paste(baseDir,'ampOddClick/r-files/functions/trStats_ampOddClick.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/functions_correctForRapitTransitions.R',sep='') )
  
  
  source( paste(baseDir,'RSkernel/r-files/functions/tapered_ccf.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/lamiarLFPtheory/smooth2DCSD.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/call_fastICA.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/save_fastICA.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/lamiarLFPtheory/showICAcomponents.R',sep='') )
  
  source('~/Dropbox/r-compile/signalProcessoR/R/signalProcessoR.R')
  
  options(warn=0)
  
  mxmua <- 40
  mxlfp <- 400
  
  colList = c('black','blue','cyan','green','orange','red')
  
  ## SELECT TASKS TO RUN          ====
  startFromScratch <- F  # rerun analyses that have alrady been performed?
  doSessionAverage <- F  # simple averages across ISI, pitch, etc  
  doCSDplots       <- F  # a couple of simple csd plots 
  doTrLM           <- F  # time-resolved analysis of variance
  doSTPSPbyChan    <- F  # STPSP model for one single outcome variable per channel
  doSTPSP          <- F  # STPSP model for one single outcome variable per probe
  doScalarMUA      <- F  # scalar values for MUA
  doBSP            <- F  # create CSD plots for BSP filtered data 
  
  doCalcICA        <- F  # calculate ICA components from evoked activity
  #doICAWavelet     <- T  # calculate single-trial wavelet transform
  
  ## CHECK if already processed   ============
  CSDfile       <- paste(pp(), '/CSDImageByPitch_1.pdf', sep='')
  doneCSDplots  <-  file.exists( CSDfile ) & !startFromScratch
  doCSDplots    <- doCSDplots & !doneCSDplots
  
  ICAfile      <- paste(pp(), 'ICA_summary_Npca8_evk.pdf',sep='')
  doneCalcICA  <-  file.exists( ICAfile ) & !startFromScratch
  doCalcICA    <- doCalcICA & !doneCalcICA
  
  STPSPDir         <- paste(rda(), '/STPSP/', df$session[i], '/', sep='')
  doneSTPSPbyChan  <-  file.exists( STPSPDir ) & !startFromScratch
  doSTPSPbyChan    <- doSTPSPbyChan & !doneSTPSPbyChan
  
  RSscalarDir        <- paste('/Volumes/EEGCyno/RES/ampOddClick/RSscalar/', df$session[i], '/RSscalar*.rda', sep='')
  if (df$probe[i]>1)
    RSscalarDir        <- paste('/Volumes/EEGCyno/RES/ampOddClick/RSscalar/', df$session[i],'_',df$probe[i], '/RSscalar*.rda', sep='')
  
  NRSfiles           <- length( system(paste('ls ', RSscalarDir, sep=''), intern=T))
  doneScalarMUA      <-  NRSfiles > 100 & !startFromScratch
  doScalarMUA        <- doScalarMUA & !doneScalarMUA
  
  doAny     <- doCalcICA | doScalarMUA | doSessionAverage | doCSDplots | doTrLM | doSTPSPbyChan | doSTPSP
  loadMUA   <-             doScalarMUA | doSessionAverage | doCSDplots | doTrLM | doSTPSPbyChan | doSTPSP
  
  if (!doAny){
    print('Nothing to be done')
    return(NULL)
  }
  
  ### LOAD MUA                       ====
  chInd      <- df$startVIndex[i] - 1 + 1:df$numChan[i]
  chStr      <- paste('ch',chInd,sep='')
  
  if (loadMUA){
    muaR       <- loadSignalProcessor(  df[i,],chStr=chStr, project=project, epoch='mua',blrng=NA )
    muaR$daxis <- -1*(1:dim(muaR$lfp)[3] - df$Layer[i]) * df$spacing[ i ]
    
    mua          <- loadSignalProcessor(  df[i,],chStr=chStr, project=project, epoch='mua',blrng=c(-50,0) )
    mua$daxis    <- -1*(1:dim(mua$lfp)[3] - df$Layer[i]) * df$spacing[ i ]
    
    muaC         <- removeCommonReferenceNoise_ica( mua,   Npca=3, ica_cutoff=1/4, show=F, check=F)
    tmp          <- apply(muaC$lfp,c(1,3),max) - apply(muaC$lfp,c(1,3),min)
    muaC$xpp$p2p  <- apply( tmp, 1, max )
    
    icaCleanMUA  <- 0
    
    if (icaCleanMUA){
      mua <- muaC
    }
  }
  
  
  ### LOAD LFP (RAW)                 ====
  lfpR          <- loadSignalProcessor(  df[i,],chStr=chStr, project=project, epoch='raw', blrng=NA )
  lfpR$daxis    <- -1*(1:dim(lfpR$lfp)[3] - df$Layer[i]) * df$spacing[ i ]
  
  lfp           <- lfpR
  tmp           <- baseLineCorrection(lfpR$lfp,c(-50,0), lfpR$taxis)  
  lfp$lfp       <- tmp$res
  
  
  ### Identify TYPES                 =============================
  tblTrTy <- c(NA,NA)
  tblTrTy[1] <- length(which(lfp$xpp$trialType==1))
  tblTrTy[2] <- length(which(lfp$xpp$trialType==2))
  tmp     <- tblTrTy > (length(lfp$xpp$trialType)/4)
  
  doTrTy     <- c( tmp[1] & tmp[2], tmp)
  names(doTrTy) <- c('all','random','regular')
  
  
  
  ### LOAD EEG                       ====
  eeg          <- loadSignalProcessor(  df[i,] ,chStr=paste('ch',1:33,sep=''), project=project, epoch='raw'  )
  eeg$daxis    <- 1:dim(eeg$lfp)[3]
  
  
  ### set units for CSD plots        ====
  lfpUnit    <- 40    # [uV]      per channel-distance
  muaUnit    <- 5     # [uV]      per channel-distance
  csdUnit    <-.5     # [mV/mm2] per channel-distance 
  ffrUnit    <- 2#4     # [uV]      per channel-distance
  csdffrUnit <- 0.01  # [mV/mm2] per channel-distance 
  
  lfpscl <- 1/lfpUnit * df$spacing[ind] # [uV]
  csdscl <- 1/csdUnit * df$spacing[ind] # [mV/mm2]
  muascl <- 1/muaUnit * df$spacing[ind] # [uV]
  
  ffrscl <- 1/ffrUnit * df$spacing[ind] # [uV]
  csdffrscl <- 1/csdffrUnit * df$spacing[ind] # [mV/mm2]
  
  # for making plots that match FFR
  lfpUnit250    <- 30    # [uV]      per channel-distance
  muaUnit250    <- 1     # [uV]      per channel-distance
  csdUnit250    <-.1     # [mV/mm2] per channel-distance 
  
  lfpscl250 <- 1/lfpUnit250 * df$spacing[ind] # [uV]
  csdscl250 <- 1/csdUnit250 * df$spacing[ind] # [mV/mm2]
  muascl250 <- 1/muaUnit250 * df$spacing[ind] # [uV]
  
  
  ## ========================================================= p2p mua
  #qtmp         <- quantile(mua$xpp$p2p,probs=c(0.5,0.75))
  #mxmua        <- qtmp[1]+10*diff(qtmp)
  #valIndMUA    <- mua$xpp$p2p < mxmua;
  
  ## ========================================================= p2p muaC
  #qtmp         <- quantile(muaC$xpp$p2p,probs=c(0.5,0.75))
  #mxmuaC        <- qtmp[1]+10*diff(qtmp)
  #valIndMUAC   <- muaC$xpp$p2p < mxmuaC;
  
  ## ========================================================= p2p lfp
  #qtmp         <- quantile(lfp$xpp$p2p,probs=c(0.5,0.75))
  #mxlfp        <- qtmp[1]+10*diff(qtmp)
  #valIndLFP    <- lfp$xpp$p2p < mxlfp;
  
  
  
  ### SINGLE TRIAL CSD               =====================================
  lag     <- 1
  smooth  <- T
  sd      <- 100
  spacing <- df[i,]$spacing
  varCor  <- F
  
  # new version from RSkernel
  csd      <- calcSingleTrialCSD(lfpR, spacing=spacing, lag=lag,smooth=T, sd=sd)
  
  tmp          <- apply(csd$lfp,c(1,3),max) - apply(csd$lfp,c(1,3),min)
  csd$xpp$p2p  <- apply( tmp, 1, max )
  
  #qtmp         <- quantile(csd$xpp$p2p,probs=c(0.5,0.75))
  #mxcsd        <- qtmp[1]+10*diff(qtmp)
  #valIndCSD    <- csd$xpp$p2p < mxcsd
  
  
  ### CALCULATE AVREC                ====================================================
  avRec <- csd
  borderChannels <- c(1,dim(csd$lfp)[3] )
  avRec$lfp <- apply(csd$lfp[,,-borderChannels], c(1,2), function(x){ mean(abs(x)) } )
  dim(avRec$lfp) <- c(dim(csd$lfp)[c(1,2)],1)
  avRec$channelPosition <- csd$channelPosition[1,]
  avRec$chStr <- 'avRec'
  
  tmp       <- baseLineCorrection( avRec$lfp, c(-50,0), avRec$taxis )
  avRec$lfp <- tmp$res
  
  tmp            <- apply(avRec$lfp,c(1,3),max) - apply(avRec$lfp,c(1,3),min)
  avRec$xpp$p2p  <- apply( tmp, 1, max )
  
  qtmp              <- quantile(avRec$xpp$p2p,probs=c(0.5,0.75))
  mxavrec           <- qtmp[1]+10*diff(qtmp)
  #valIndAVRec       <- avRec$xpp$p2p < mxavrec
  avRec$xpp$valTone <- avRec$xpp$p2p < mxavrec
  
  
  
  
  # old CSD version from ampOddClick
  # llp          <- loadSignalProcessor(  df[i,],chStr=chStr, project=project, epoch='llp' )
  # 
  # ## ========================================================= p2p llp
  # qtmp         <- quantile(llp$xpp$p2p,probs=c(0.5,0.75))
  # mxlfp        <- qtmp[1]+10*diff(qtmp)
  # valIndLLP    <- llp$xpp$p2p < mxlfp;
  # 
  # ## ========================================================= p2p eeg
  # qtmp         <- quantile(eeg$xpp$p2p,probs=c(0.5,0.75))
  # mxeeg        <- qtmp[1]+10*diff(qtmp)
  # valIndEEG    <- eeg$xpp$p2p < mxeeg;
  # 
  
  

  ####   RUN FUNCTIONS ====
  if (doScalarMUA){
    ## =================================================================================================== analyses per session (mostly CSD)
    
    cmpStrDF     <- data.frame( what='csd', cmpStrName='ARC0',  tmin=8,  tmax=13,  dmin= -1500, dmax= -900, stringsAsFactors=FALSE) 
    cmpStrDF[2,] <- data.frame( what='csd', cmpStrName='ARC1',  tmin=15,  tmax=40,  dmin= -500,  dmax=  500, stringsAsFactors=FALSE)
    cmpStrDF[3,] <- data.frame( what='csd', cmpStrName='ARC2',  tmin=45,  tmax=55,  dmin= -500,  dmax=  500, stringsAsFactors=FALSE)
    cmpStrDF[4,] <- data.frame( what='csd', cmpStrName='ARC3',  tmin=65,  tmax=110, dmin= -500,  dmax=  500, stringsAsFactors=FALSE)
    cmpStrDF[5,] <- data.frame( what='csd', cmpStrName='ARC4',  tmin=100, tmax=150,  dmin= -500, dmax=  500, stringsAsFactors=FALSE)
    cmpStrDF[6,] <- data.frame( what='csd', cmpStrName='ARC5',  tmin=140, tmax=240,  dmin= -500, dmax=  500, stringsAsFactors=FALSE)
    
    #cmpStrDF <- cmpStrDF[c(2,3),]
    
    for (cx in 1:dim(cmpStrDF)[1]){
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
      #chDepth       <- (df$Layer[i] - 1:df$numChan[i]) * df$spacing[i]
      chDepth       <- (df$Layer[i] - dat$chDepthInd) * df$spacing[i]
      
      chnd          <- which(chDepth>=dpthrng[1] & chDepth<=dpthrng[2])
      if (is.na(df$Layer[i]))
        chnd <- 1:dim(dat$lfp)[3]
      
      if(length(chnd)>0){# channels present in depth range?
        DVstr         <- cmpStrDF$cmpStrName[cx]
        
        dat$xpp$TMP   <- apply(dat$lfp[,tnd,chnd], c(1), mean)
        if (identical(cmpStrDF$what[cx],'csd'))# CSD should default to avRec across channels, not mean
          dat$xpp$TMP   <- apply(abs(dat$lfp[,tnd,chnd]), c(1), function(x){mean(abs(x))})
        
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
    
    ## ==================================================================================================== analyze by channel (mostly MUA)
    cmpStrDF     <- data.frame( what='mua', cmpStrName='MUA0',  tmin=8,  tmax=13,  stringsAsFactors=FALSE) 
    cmpStrDF[2,] <- data.frame( what='mua', cmpStrName='MUA1',  tmin=18,  tmax=35,  stringsAsFactors=FALSE)
    cmpStrDF[3,] <- data.frame( what='mua', cmpStrName='MUA1a', tmin=15,  tmax=25,  stringsAsFactors=FALSE)
    cmpStrDF[4,] <- data.frame( what='mua', cmpStrName='MUA1b', tmin=25,  tmax=35,  stringsAsFactors=FALSE)
    cmpStrDF[5,] <- data.frame( what='mua', cmpStrName='MUA1c', tmin=35,  tmax=45,  stringsAsFactors=FALSE)
    cmpStrDF[6,] <- data.frame( what='mua', cmpStrName='MUA2' , tmin=55,  tmax=125, stringsAsFactors=FALSE)
    cmpStrDF[7,] <- data.frame( what='mua', cmpStrName='MUA3' , tmin=160, tmax=250, stringsAsFactors=FALSE)
    cmpStrDF[8,] <- data.frame( what='mua', cmpStrName='MUAS' , tmin=1, tmax=70, stringsAsFactors=FALSE)
    
    cmpStrDF[9,]  <- data.frame( what='csd', cmpStrName='CSD0',  tmin=15,   tmax=19,  stringsAsFactors=FALSE)
    cmpStrDF[10,] <- data.frame( what='csd', cmpStrName='CSD1',  tmin=25,   tmax=45,  stringsAsFactors=FALSE)
    cmpStrDF[11,] <- data.frame( what='csd', cmpStrName='CSD2',  tmin=50,   tmax=60,  stringsAsFactors=FALSE)
    cmpStrDF[12,] <- data.frame( what='csd', cmpStrName='CSD3',  tmin=65,   tmax=110, stringsAsFactors=FALSE)
    cmpStrDF[13,] <- data.frame( what='csd', cmpStrName='CSD4',  tmin=100,  tmax=150, stringsAsFactors=FALSE)
    cmpStrDF[14,] <- data.frame( what='csd', cmpStrName='CSD5',  tmin=160,  tmax=250, stringsAsFactors=FALSE)
    
    cmpStrDF[15,] <- data.frame( what='lfp', cmpStrName='LFP0',  tmin=10,  tmax=17, stringsAsFactors=FALSE)
    cmpStrDF[16,] <- data.frame( what='lfp', cmpStrName='LFP1',  tmin=20,  tmax=30, stringsAsFactors=FALSE)
    
    cmpStrDF[17,] <- data.frame( what='muaR', cmpStrName='MUAbl',  tmin= -50,  tmax=0,   stringsAsFactors=FALSE)
    cmpStrDF[18,] <- data.frame( what='muaR', cmpStrName='MUAR1',  tmin=18,    tmax=35,  stringsAsFactors=FALSE)
    cmpStrDF[19,] <- data.frame( what='muaR', cmpStrName='MUAR1a', tmin=15,    tmax=25,  stringsAsFactors=FALSE)
    cmpStrDF[20,] <- data.frame( what='muaR', cmpStrName='MUAR1b', tmin=25,    tmax=35,  stringsAsFactors=FALSE)
    cmpStrDF[21,] <- data.frame( what='muaR', cmpStrName='MUAR1c', tmin=35,    tmax=45,  stringsAsFactors=FALSE)
    cmpStrDF[22,] <- data.frame( what='muaR', cmpStrName='MUAR2' , tmin=55,    tmax=125, stringsAsFactors=FALSE)
    cmpStrDF[23,] <- data.frame( what='muaR', cmpStrName='MUAR3' , tmin=160,   tmax=250, stringsAsFactors=FALSE)
    
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
        
        #DVstr         <- paste('ch', df$startVIndex[i]+chnd, '_', cmpStrDF$cmpStrName[cx], sep='')
        DVstr         <- paste(chStr[chnd], '_', cmpStrDF$cmpStrName[cx], sep='')
        
        if(chnd <= dim(dat$lfp)[3]){ # csd has fewer channels on the d-probe...
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
            if (1==0){ #(df$corMUAbl[i]==2){
              print(paste0('correcting for slow fluctuations; Ch. ', chnd, '; Epoch:', cmpStrDF$cmpStrName[cx]))
              mavDV            <- get.mav(dat$xpp[[DVstr]], 1:dim(dat$xpp)[1], eval.at=seq(1,dim(dat$xpp)[1],by=1 ), show=T, distribution='Gaussian',boot=F, width=100)
              dat$xpp[[DVstr]] <- dat$xpp[[DVstr]] - mavDV$mav + mean(mavDV$mav)
            }
          }
          
          dat$xpp$ampIndFactor <- as.factor(dat$xpp$ampInd)
          
          # --------------- run analyzeRSscalar
          #getValidTrialsAOC is returning logical(0)
          valBin <- getValidTrialsAOC(dat, who='all',DVstr=DVstr)
          disc   <- analyzeAOCscalar(dat, ybin=valBin, who='all', cmpStr=DVstr)
          
          valBin <- getValidTrialsAOC(dat, who='regular')
          disc   <- analyzeAOCscalar(dat, ybin=valBin, who='regular', cmpStr=DVstr)
          
          valBin <- getValidTrialsAOC(dat, who='random')
          disc   <- analyzeAOCscalar(dat, ybin=valBin, who='random', cmpStr=DVstr)
        }# end if chnd < 
      }# end cmpStrList loop
    }# end channel loop
    
  }              # linear models of scalar MUA values
  if (doCalcICA) {
    ## ================================================== calc ICA if requested
    Npca <- 6
    ica6 <- call_fastICA(lfp, Npca=Npca, subset=valIndLFP, valChan=1:df$numChan[i], show=F, restoreAll=T)
    disc <- save_fastICA(ica6, extn='evk_',trda=trda)
    showICAcomponents(ica6,srtICA=T,csdWeights=T, showFFT=T, showEvokedResponse=T, logtime=F, subset=valIndLFP, includeexp=F,eps=T, pp=pp, fxtn='_evk')
    showICAcomponents(ica6,srtICA=T,csdWeights=T, showFFT=T, showEvokedResponse=T, logtime=T, subset=valIndLFP, includeexp=F,eps=T, pp=pp, fxtn='_evkLog', trng=c(3,500))
    
    Npca <- 8
    ica8 <- call_fastICA(lfp, Npca=Npca, subset=valIndLFP, valChan=1:df$numChan[i], show=F, restoreAll=T)
    disc <- save_fastICA(ica8, extn='evk_',trda=trda)
    showICAcomponents(ica8,srtICA=T,csdWeights=T, showFFT=T, showEvokedResponse=T, subset=valIndLFP, includeexp=F,eps=T, pp=pp, fxtn='_evk')
    showICAcomponents(ica8,srtICA=T,csdWeights=T, showFFT=T, showEvokedResponse=T, subset=valIndLFP, includeexp=F,eps=T, pp=pp, fxtn='_evklog',logtime=T, trng=c(3,500) )
    
    Npca <- 6
    ica6 <- call_fastICA(llp, Npca=Npca, subset=valIndLLP, valChan=1:df$numChan[i], show=F, restoreAll=T)
    disc <- save_fastICA(ica6, extn='llp_',trda=trda)
    showICAcomponents(ica6,srtICA=T,csdWeights=T, showFFT=T, showEvokedResponse=T, subset=valIndLFP, includeexp=F,eps=T, pp=pp, fxtn='_llp')
    showICAcomponents(ica6,srtICA=T,csdWeights=T, showFFT=T, showEvokedResponse=T, subset=valIndLFP, includeexp=F,eps=T, pp=pp, fxtn='_llplog', logtime=T, trng=c(3,500))
    
    #Npca <- 6
    #icaffr6 <- call_fastICA(ffr, Npca=Npca, subset=valIndLFP, valChan=1:df$numChan[i], show=F, restoreAll=T)
    #disc <- save_fastICA(icaffr6, extn='evk_ffr_',trda=trda)
    #showICAcomponents(icaffr6,srtICA=T,csdWeights=T, showFFT=T, showEvokedResponse=T, subset=valIndLFP, includeexp=F,eps=T, pp=pp, fxtn='_FFR')
    
    sumEEG <- sum(eeg$xpp$p2p<450)
    OKEEG <- sumEEG>dim(eeg$lfp)[1]/3
    
    if (OKEEG){
      Npca <- 8
      eegica8 <- call_fastICA(eeg, Npca=Npca, subset=valIndEEG, valChan=1:33, show=F, restoreAll=T)
      disc    <- save_fastICA(eegica8, extn='EEG_evk_',trda=trda)
      showICAcomponents(eegica8,srtICA=T,csdWeights=F, showFFT=T, showEvokedResponse=T, subset=valIndEEG, includeexp=F,eps=T, pp=pp, fxtn='_EEG')
      
      Npca <- 12
      eegica12 <- call_fastICA(eeg, Npca=Npca, subset=valIndEEG, valChan=1:33, show=F, restoreAll=T)
      disc <- save_fastICA(eegica12, extn='EEG_evk_',trda=trda)
      showICAcomponents(eegica12,srtICA=T,csdWeights=F, showFFT=T, showEvokedResponse=T, subset=valIndEEG, includeexp=F,eps=T, pp=pp, fxtn='_EEG')
      
      Npca <- 24
      eegica24 <- call_fastICA(eeg, Npca=Npca, subset=valIndEEG, valChan=1:33, show=F, restoreAll=T)
      disc <- save_fastICA(eegica24, extn='EEG_evk_',trda=trda)
      showICAcomponents(eegica32,srtICA=T,csdWeights=F, showFFT=T, showEvokedResponse=T, subset=valIndEEG, includeexp=F,eps=T, pp=pp, fxtn='_EEG')
    }
  }               # calculate ICA 
  if (doCSDplots){
    print('CSDplots')
    
    ## ==============================================  GA CSD 
    csdGA      <- calcCSD(lfp,  mua, lag=lag, mxmua=mxmua,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sd,varCor=varCor,spacing=spacing)
    csdGAC     <- calcCSD(lfp, muaC, lag=lag, mxmua=mxmuaC,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sd,varCor=varCor,spacing=spacing)
    
    eps.file('GAC_CSDImage.eps',pdf=T,s=5.50)
    disc <- showCSDImage(csdGAC,trng=c(-10,250),muascl=muascl,lfpscl=lfpscl,csdscl=csdscl,zlim=c(-2,2)*lfpUnit,demeanLFP=F,tat=NA,addcsd=T)
    dev.off()
    
    eps.file('GAR_CSDImage.eps',pdf=T,s=5.50)
    disc <- showCSDImage(csdGA,trng=c(-10,250),muascl=muascl,lfpscl=lfpscl,csdscl=csdscl,zlim=c(-2,2)*lfpUnit,demeanLFP=F,tat=NA,addcsd=T)
    dev.off()
    
    eps.file('GA_CSDImage.eps',pdf=T,s=5.50)
    disc <- showCSDImage(csdGA,trng=c(-10,250),muascl=muascl,lfpscl=lfpscl,csdscl=csdscl,zlim=c(-2,2)*lfpUnit,demeanLFP=F,tat=NA,addcsd=T)
    dev.off()
    
    eps.file('GA_CSDImage_log.eps',pdf=T,s=5.50)
    disc <- showCSDImage(csdGA,trng=c(3,500),muascl=muascl,lfpscl=lfpscl,csdscl=csdscl,zlim=c(-2,2)*lfpUnit,demeanLFP=F,tat=NA,addcsd=T, logtime=T)
    dev.off()
    
    eps.file('GA_CSDImage_demean.eps',pdf=T,s=5.50)
    disc <- showCSDImage(csdGA,trng=c(-10,250),muascl=muascl,lfpscl=lfpscl,csdscl=csdscl,zlim=c(-2,2)*lfpUnit,demeanLFP=T,tat=NA,addcsd=T)
    dev.off()
    
    eps.file('GA_CSDImage_log_demean.eps',pdf=T,s=5.50)
    disc <- showCSDImage(csdGA,trng=c(3,500),muascl=muascl,lfpscl=lfpscl,csdscl=csdscl,zlim=c(-2,2)*lfpUnit,demeanLFP=T,tat=NA,addcsd=T, logtime=T)
    dev.off()
    
    csdGA      <- calcCSD(lfp,  mua, lag=lag, mxmua=mxmua,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sd,varCor=varCor,spacing=spacing)
    csdGA250   <- calcCSD(lfp, mua, lag=lag, mxmua=mxmua,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=250,varCor=varCor,spacing=spacing)
    eps.file('GA_CSDImage250.eps',pdf=T,s=5.50)
    disc <- showCSDImage(csdGA250,trng=c(-10,350),muascl=muascl250,lfpscl=lfpscl250,csdscl=csdscl250,zlim=c(-1,1)*lfpUnit250,demeanLFP=F,lwd=.5,addcsd=F)
    dev.off()
    
    eps.file('GA_CSDImage250_demean.eps',pdf=T,s=5.50)
    disc <- showCSDImage(csdGA250,trng=c(-10,350),muascl=muascl250,lfpscl=lfpscl250,csdscl=csdscl250,zlim=c(-1,1)*lfpUnit250,demeanLFP=T,lwd=.5,addcsd=F)
    dev.off()
    
    
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
      eps.file(paste('CSDImageBySOA_',j,'.eps',sep=''),pdf=T,s=5.50)
      disc <- showCSDImage(isicsdList[[j]],lfpscl=lfpscl, csdscl=csdscl, muascl=muascl,zlim=c(-2,2)*lfpUnit, logtime=T, trng=c(-10,250), addcsd=T)
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
      eps.file(paste('CSDImageByAmp_',j,'.eps',sep=''),pdf=T,s=5.50)
      disc <- showCSDImage(pitchcsdList[[j]],lfpscl=lfpscl, csdscl=csdscl, muascl=muascl,zlim=c(-2,2)*lfpUnit,logtime=T,trng=c(-10,250),addcsd=T)
      dev.off()
    }
    
    
    ## AVREC
    colList <- c('black','blue','cyan','green','orange','red')
    eps.file('AVRECByAmp.eps',pdf=T,s=2.5,ratio=3/2)
    disc <- plot(pitchcsdList[[5]]$txlfp, pitchcsdList[[5]]$avrec,xlim=c(-25,250), col=colList[1],type='l',lwd=2)
    for (j in 1:5)
      disc <- points(pitchcsdList[[1]]$txlfp,pitchcsdList[[j]]$avrec, col=colList[j],type='l',lwd=2)
    dev.off()
  }               # CSD grand average plots
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
    sumEEG <- sum(eeg$xpp$p2p<450)
    OKEEG <- sumEEG>dim(eeg$lfp)[1]/3
    OKEEG <- F
    if (OKEEG){
      for (chNr in 2:dim(eeg$lfp)[3]){
        print(chNr)
        
        eeg$ch <- chNr
        disc   <- trStats.ampOddClick(eeg,width=.050,eval.at=eeg$taxis,variables=variables,what='lfp',chNr=chNr,
                                      subs=valIndEEG,extn='_eeg_all',show=F, fix=FALSE)
        disc   <- trStats.ampOddClick(eeg,width=.050,eval.at=eeg$taxis,variables=variables2,what='lfp',chNr=chNr,
                                      subs=valIndEEG&eeg$xpp$trialType==1,extn='_eeg_random',show=F, fix=FALSE)
        disc   <- trStats.ampOddClick(eeg,width=.050,eval.at=eeg$taxis,variables=variables2,what='lfp',chNr=chNr,
                                      subs=valIndEEG&eeg$xpp$trialType==2,extn='_eeg_regular',show=F, fix=FALSE)
      }
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
      
      #Nchan   <- dim(mua$lfp)[3]
      csd$ch <- chNr
      if(chNr>1 & chNr<=dim(csd$lfp)[3]){
        disc   <- trStats.ampOddClick(csd,width=.050,eval.at=csd$taxis,variables=variables,what='lfp',chNr=chNr,
                                      subs=valIndCSD,extn='_csd_all',show=F, fix=FALSE)
        disc   <- trStats.ampOddClick(csd,width=.050,eval.at=csd$taxis,variables=variables2,what='lfp',chNr=chNr,
                                      subs=valIndCSD&csd$xpp$trialType==1,extn='_csd_random',show=F, fix=FALSE)
        disc   <- trStats.ampOddClick(csd,width=.050,eval.at=csd$taxis,variables=variables2,what='lfp',chNr=chNr,
                                      subs=valIndCSD&csd$xpp$trialType==2,extn='_csd_regular',show=F, fix=FALSE)
      }
    }# chNr
  }                   # time-resolved linar models
  if (doSTPSPbyChan | doSTPSP){
    trng    <- c(20,70)
    blrng   <- c(-50,0)
    minISI  <- 0.190
    maxISI  <- 12.8
    DVscl   <- 5
    sfc     <- 1
    Nchan   <- dim(mua$lfp)[3]
    muach   <- 1:Nchan
  }  # set STPSP parameters
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
        dat               <- csd
        Nchan             <- dim(dat$lfp)[3]
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
      #function(vmua,minISI=0.190, maxISI=25,     DVscl=5,    sfc=1, who='all',      parline=NA,      dfline=NA)
      STPSPmuaRnd    <- prepSTPSPVProbeMUA(dat, minISI=minISI,maxISI=maxISI,DVscl=DVscl,sfc=sfc,who='random' , parline=parline, dfline=df[i,])# crashes for LFP ch32 regular
      STPSPmuaReg    <- prepSTPSPVProbeMUA(dat, minISI=minISI,maxISI=maxISI,DVscl=DVscl,sfc=sfc,who='regular', parline=parline, dfline=df[i,])
    } 
  }            # STPSP for individual channels
  if (doSTPSP){# ===================================================== not yet adjusted to ampOddClick
    
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
  }                  # STPSP for whole probe, e.g., CSD
  if (doBSP){
    
    bsp          <- loadSignalProcessor(  df[i,],chStr=chStr, project=project, epoch='bsp' )
    hdmua        <- loadSignalProcessor(  df[i,],chStr=chStr, project=project, epoch='hdmua' )
    hdmuaC       <- removeCommonReferenceNoise_ica( hdmua,   Npca=3, ica_cutoff=1/4, show=F, check=F)
    
    
    bsp$daxis       <- -1*( bsp$chDepthInd    - df$Layer[i]) * df$spacing[ i ]
    hdmua$daxis     <- -1*( hdmua$chDepthInd  - df$Layer[i]) * df$spacing[ i ]
    hdmuaC$daxis    <- -1*( hdmuaC$chDepthInd - df$Layer[i]) * df$spacing[ i ]
    
    bsp$exp      <- paste( df$session[i], dxt, pxt, sep='')
    bsp$dataSet  <- 'g'
    
    
    sdFFR   <- 250
    
    #muascl  <- 30
    #lfpscl  <- 4
    #csdscl  <- 300 #500#1300
    
    csdGA     <- calcCSD(subs.data(bsp,bsp$xpp$ampInd >0), subs.data(hdmua,hdmua$xpp$ampInd >0), lag=lag, mxmua=mxlfp,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sdFFR,varCor=varCor,spacing=spacing)
    csdGAC    <- calcCSD(subs.data(bsp,bsp$xpp$ampInd >0), subs.data(hdmuaC,hdmuaC$xpp$ampInd >0), lag=lag, mxmua=mxlfp,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sdFFR,varCor=varCor,spacing=spacing)
    
    
    #disc <- showCSDImage(csdGA,trng=c(-5,25), muascl=5,lfpscl=100,csdscl=800,tat=NA,addcsd=F)
    #muascl <- 5
    #lfpscl <- 40#20
    #csdscl <- 1200#400
    
    trng   <- c(-1,25) 
    
    eps.file('ABRR_GA_CSDImage.eps',pdf=T,s=5.80)
    disc <- showCSDImage(csdGA,trng=trng, muascl=muascl,lfpscl=ffrscl,csdscl=csdffrscl,lwd=0.5,zlim=c(-2,2)*ffrUnit,demeanLFP=T )
    dev.off()
    
    eps.file('ABRC_GA_CSDImage.eps',pdf=T,s=5.80)
    disc <- showCSDImage(csdGAC,trng=trng, muascl=muascl,lfpscl=ffrscl,csdscl=csdffrscl,lwd=0.5,zlim=c(-2,2)*ffrUnit,demeanLFP=T )
    dev.off()
    
    icaCleanMUAhd <- TRUE
    if (icaCleanMUAhd){
      hdmua <- hdmuaC
    }
    
    csdGA     <- calcCSD(subs.data(bsp,bsp$xpp$ampInd >0), subs.data(hdmua,hdmua$xpp$ampInd >0), lag=lag, mxmua=mxlfp,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sdFFR,varCor=varCor,spacing=spacing)
    
    eps.file('ABR_GA_CSDImage.eps',pdf=T,s=5.80)
    disc <- showCSDImage(csdGA,trng=trng, muascl=muascl,lfpscl=ffrscl,csdscl=csdffrscl,lwd=0.5,zlim=c(-2,2)*ffrUnit,demeanLFP=T )
    dev.off()
    
    ## CSDs by click intensity:
    csdGA1   <- calcCSD(subs.data(bsp,bsp$xpp$ampInd==1), subs.data(hdmua,hdmua$xpp$ampInd==1), lag=lag, mxmua=mxlfp,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sdFFR,varCor=varCor,spacing=spacing)
    csdGA2   <- calcCSD(subs.data(bsp,bsp$xpp$ampInd==2), subs.data(hdmua,hdmua$xpp$ampInd==2), lag=lag, mxmua=mxlfp,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sdFFR,varCor=varCor,spacing=spacing)
    csdGA3   <- calcCSD(subs.data(bsp,bsp$xpp$ampInd==3), subs.data(hdmua,hdmua$xpp$ampInd==3), lag=lag, mxmua=mxlfp,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sdFFR,varCor=varCor,spacing=spacing)
    csdGA4   <- calcCSD(subs.data(bsp,bsp$xpp$ampInd==4), subs.data(hdmua,hdmua$xpp$ampInd==4), lag=lag, mxmua=mxlfp,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sdFFR,varCor=varCor,spacing=spacing)
    csdGA5   <- calcCSD(subs.data(bsp,bsp$xpp$ampInd==5), subs.data(hdmua,hdmua$xpp$ampInd==5), lag=lag, mxmua=mxlfp,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sdFFR,varCor=varCor,spacing=spacing)
    
    
    eps.file('ABR_CSDImageByAmp_1.eps',pdf=T,s=5.80)
    disc <- showCSDImage(csdGA1,trng=trng, muascl=muascl,lfpscl=ffrscl,csdscl=csdffrscl,zlim=c(-2,2)*ffrUnit,addcsd=F,demeanLFP=T)
    dev.off()
    eps.file('ABR_CSDImageByAmp_2.eps',pdf=T,s=5.80)
    disc <- showCSDImage(csdGA2,trng=trng, muascl=muascl,lfpscl=ffrscl,csdscl=csdffrscl,zlim=c(-2,2)*ffrUnit,addcsd=F,demeanLFP=T)
    dev.off()
    eps.file('ABR_CSDImageByAmp_3.eps',pdf=T,s=5.80)
    disc <- showCSDImage(csdGA3,trng=trng, muascl=muascl,lfpscl=ffrscl,csdscl=csdffrscl,zlim=c(-2,2)*ffrUnit,addcsd=F,demeanLFP=T)
    dev.off()
    eps.file('ABR_CSDImageByAmp_4.eps',pdf=T,s=5.80)
    disc <- showCSDImage(csdGA4,trng=trng, muascl=muascl,lfpscl=ffrscl,csdscl=csdffrscl,zlim=c(-2,2)*ffrUnit,addcsd=F,demeanLFP=T)
    dev.off()
    eps.file('ABR_CSDImageByAmp_5.eps',pdf=T,s=5.80)
    disc <- showCSDImage(csdGA5,trng=trng, muascl=muascl,lfpscl=ffrscl,csdscl=csdffrscl,zlim=c(-2,2)*ffrUnit,addcsd=F,demeanLFP=T)
    dev.off()
    
    
    ## CSDs by ICI:
    bsp$xpp$ICIoct[is.na(bsp$xpp$ICIoct)] <- levels(bsp$xpp$ICIoct)[1]
    hdmua$xpp$ICIoct[is.na(hdmua$xpp$ICIoct)] <- levels(bsp$xpp$ICIoct)[1]
    csdICI2   <- calcCSD(subs.data(bsp,c(bsp$xpp$ICIoct)==2), subs.data(hdmua,c(hdmua$xpp$ICIoct)==2), lag=lag, mxmua=mxlfp,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sdFFR,varCor=varCor,spacing=spacing)
    csdICI3   <- calcCSD(subs.data(bsp,c(bsp$xpp$ICIoct)==3), subs.data(hdmua,c(hdmua$xpp$ICIoct)==3), lag=lag, mxmua=mxlfp,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sdFFR,varCor=varCor,spacing=spacing)
    csdICI4   <- calcCSD(subs.data(bsp,c(bsp$xpp$ICIoct)==4), subs.data(hdmua,c(hdmua$xpp$ICIoct)==4), lag=lag, mxmua=mxlfp,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sdFFR,varCor=varCor,spacing=spacing)
    csdICI5   <- calcCSD(subs.data(bsp,c(bsp$xpp$ICIoct)==5), subs.data(hdmua,c(hdmua$xpp$ICIoct)==5), lag=lag, mxmua=mxlfp,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sdFFR,varCor=varCor,spacing=spacing)
    csdICI6   <- calcCSD(subs.data(bsp,c(bsp$xpp$ICIoct)==6), subs.data(hdmua,c(hdmua$xpp$ICIoct)==6), lag=lag, mxmua=mxlfp,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sdFFR,varCor=varCor,spacing=spacing)
    csdICI7   <- calcCSD(subs.data(bsp,c(bsp$xpp$ICIoct)==7), subs.data(hdmua,c(hdmua$xpp$ICIoct)==7), lag=lag, mxmua=mxlfp,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sdFFR,varCor=varCor,spacing=spacing)
    
    
    
    eps.file('ABR_CSDImageBySOA_2.eps',pdf=T,s=5.80)
    disc <- showCSDImage(csdICI2,trng=trng, muascl=muascl,lfpscl=ffrscl,csdscl=csdffrscl,zlim=c(-2,2)*ffrUnit,addcsd=F,demeanLFP=T)
    dev.off()
    
    eps.file('ABR_CSDImageBySOA_3.eps',pdf=T,s=5.80)
    disc <- showCSDImage(csdICI3,trng=trng,muascl=muascl,lfpscl=ffrscl,csdscl=csdffrscl,zlim=c(-2,2)*ffrUnit,addcsd=F,demeanLFP=T)
    dev.off()
    
    eps.file('ABR_CSDImageBySOA_4.eps',pdf=T,s=5.80)
    disc <- showCSDImage(csdICI4,trng=trng,muascl=muascl,lfpscl=ffrscl,csdscl=csdffrscl,zlim=c(-2,2)*ffrUnit,addcsd=F,demeanLFP=T)
    dev.off()
    
    
    eps.file('ABR_CSDImageBySOA_5.eps',pdf=T,s=5.80)
    disc <- showCSDImage(csdICI5,trng=trng, muascl=muascl,lfpscl=ffrscl,csdscl=csdffrscl,zlim=c(-2,2)*ffrUnit,addcsd=F,demeanLFP=T)
    dev.off()
    
    eps.file('ABR_CSDImageBySOA_6.eps',pdf=T,s=5.80)
    disc <- showCSDImage(csdICI6,trng=trng, muascl=muascl,lfpscl=ffrscl,csdscl=csdffrscl,zlim=c(-2,2)*ffrUnit,addcsd=F,demeanLFP=T)
    dev.off()
    
    eps.file('ABR_CSDImageBySOA_7.eps',pdf=T,s=5.80)
    disc <- showCSDImage(csdICI7,trng=trng, muascl=muascl,lfpscl=lfpscl,csdscl=csdscl,zlim=c(-2,2)*ffrUnit,tat=NA,addcsd=F)
    dev.off()
    
    
    
    #lfp$daxis    <- -1*( lfp$chDepthInd - df$Layer[i]) * df$spacing[ i ]
    #mua$daxis    <- lfp$daxis
    #muaR$daxis   <- lfp$daxis
    #llp$daxis    <- lfp$daxis
    #eeg$daxis    <- 1:dim(eeg$lfp)[3]
    
    
    #muaR$exp <- paste( df$session[i], dxt, pxt, sep='')
    #mua$exp  <- paste( df$session[i], dxt, pxt, sep='')
    #lfp$exp  <- paste( df$session[i], dxt, pxt, sep='')
    #llp$exp  <- paste( df$session[i], dxt, pxt, sep='')
    ##ffr$exp  <- paste( df$session[i], dxt, pxt, sep='')
    #eeg$exp  <- paste( df$session[i], dxt, pxt, sep='')
    
    muaR$dataSet <- 'g'
    mua$dataSet  <- 'g'
    lfp$dataSet  <- 'g'
    llp$dataSet  <- 'g'
    #ffr$dataSet  <- 'g'
    eeg$dataSet  <- 'g'
    
    #mua$xpp$ISI[ which( mua$xpp$ISI<0)]  <- 20
    #muaR$xpp$ISI[which(muaR$xpp$ISI<0)]  <- 20
    #lfp$xpp$ISI[ which( lfp$xpp$ISI<0)]  <- 20
    #llp$xpp$ISI[ which( lfp$xpp$ISI<0)]  <- 20
    ##ffr$xpp$ISI[ which( ffr$xpp$ISI<0)]  <- 20
    #eeg$xpp$ISI[ which( eeg$xpp$ISI<0)]  <- 20
    
    ## ============================== build logarithmic time-axis
    #if (1==0){
    #  toff                                   <- 0
    #  lfp$logtaxis                           <- log2( lfp$taxis + toff)
    #  #lfp$logtaxis[ (lfp$taxis + toff <= 0)] <- .01*(lfp$taxis[ which(lfp$taxis + toff <= 0) ] - 1 + toff) #- log2(taxis[min(which(taxis>-toff))] )# + .01*taxis[ (llp$taxis + toff == 0)]
    #  lfp$logtaxis[ (lfp$taxis + toff <= 1)] <- .01*(lfp$taxis[ which(lfp$taxis + toff <= 1) ] - 1 + toff) #- log2(taxis[min(which(taxis>-toff))] )# + .01*taxis[ (llp$taxis + toff == 0)]
    #  
    #  mua$logtaxis  <- lfp$logtaxis
    #  muaR$logtaxis <- lfp$logtaxis
    #  llp$logtaxis  <- lfp$logtaxis
    #  eeg$logtaxis  <- eeg$logtaxis
    #}
    
    
    ## calculate single-trial csd ====
    
    lfpnbl <- lfp
    tmp            <- baseLineCorrection(lfp$lfp,range(lfp$taxis), lfp$taxis)
    lfpnbl$lfp <- tmp$res
    csd      <- calcSingleTrialCSD(lfpnbl, spacing=spacing, lag=lag,smooth=T, sd=sd)
    csd$what <- 'csd'
    
    avrec          <- apply( abs(csd$lfp), c(1,2), sum)
    tmp            <- baseLineCorrection(avrec,range(csd$taxis), csd$taxis)
    avrec          <- tmp$res
    csd$lfp[,,dim(csd$lfp)[3]] <- avrec
    
    tmp          <- apply(csd$lfp,c(1,3),max) - apply(csd$lfp,c(1,3),min)
    csd$xpp$p2p  <- apply( tmp, 1, max )
    
    qtmp         <- quantile(csd$xpp$p2p,probs=c(0.5,0.75))
    mxcsd        <- qtmp[1]+10*diff(qtmp)
    valIndCSD    <- csd$xpp$p2p < mxcsd;
    
    #plot(csd$taxis,apply(avrec[valIndLFP,],2,mean), type='l' )
    ## CSDs by block number
    ## =========================================  CSD by SOA
    lfp$xpp$blockNr <- 1+(lfp$xpp$trialNr-1)%/%2
    Nblock <- length(unique(lfp$xpp$blockNr))
    if(Nblock>1){
      blockcsdList <- list()
      for (j in 1:Nblock){
        subs <- as.numeric(lfp$xpp$blockNr)==j 
        subs[which(is.na(subs))] <- F
        blockcsdList[[j]] <- calcCSD(subs.data(lfp,subs),subs.data(mua,subs),lag=lag, mxmua=mxmua,mxlfp=mxlfp, 
                                     interpChanInd==interpChanInd,smooth=smooth,sd=sd,spacing=spacing )
      }
      
      ## all in one image
      colListBlk <- bluered.colors(Nblock)
      eps.file('CSDbyBLK.eps',pdf=T,s=5.80)
      disc <- showCSD(blockcsdList[[1]],lfpscl=50, csdscl=.4, muascl=2,trng=c(-10,250),addpoly=F, col=colListBlk[1],add=F)
      for (j in 2:Nblock)
        disc <- showCSD(blockcsdList[[j]],lfpscl=50, csdscl=.4, muascl=2,trng=c(-10,250),addpoly=F, col=colListBlk[j],add=T)
      dev.off()
      
      # csdImage for each Block separately
      #csdscl <- 250
      #muascl <- 30
      #lfpscl <- 4
      for (j in 1:Nblock){
        eps.file(paste('CSDImageByBLK_',j,'.eps',sep=''),pdf=T,s=5.80)
        disc <- showCSDImage(blockcsdList[[j]],lfpscl=lfpscl, csdscl=csdscl, muascl=muascl,trng=c(-10,400),zlim=c(-2,2)*lfpUnit, addcsd=T,logtime=T)
        dev.off()
      }
      
      ## AVREC
      colList <- c('black','blue','cyan','green','orange','red')
      ylim <- na.range( blockcsdList[[1]]$avrec) * c(0,1.5)
      eps.file('AVRECByBLK.eps',pdf=T,s=2.5,ratio=3/2)
      disc <- plot(blockcsdList[[1]]$txlfp, blockcsdList[[1]]$avrec,xlim=c(-25,250),ylim=ylim, col=colListBlk[1],type='l',lwd=1)
      for (j in 1:Nblock)
        disc <- points(blockcsdList[[1]]$txlfp,blockcsdList[[j]]$avrec, col=colListBlk[j],type='l',lwd=1)
      dev.off()
      
    }
  }                    # CSD plots for BSP-filtered data  
  if(doSessionAverage){
    #browser()
    
    fltxtn <- ''
    
    getSessionAverage(lfp, maxp2p=mxlfp, who='all', fextn=fltxtn)
    gc()
    
    sumEEG <- sum(eeg$xpp$p2p<mxeeg)
    OKEEG <- sumEEG>dim(eeg$lfp)[1]/3
    #OKEEG <- F
    if (OKEEG){
      getSessionAverage(eeg, maxp2p=mxeeg, who='all', fextn=fltxtn)
      gc()
    }
    getSessionAverage(csd, maxp2p=mxcsd, who='all', fextn=fltxtn)
    gc()
    
    getSessionAverage(mua, maxp2p=mxmua, who='all', fextn=fltxtn)
    gc()
    
    ## random trialType==1
    vnd <- which(lfp$xpp$trialType==1)
    if (length(vnd)>100){
      
      getSessionAverage(subs.data(lfp,lfp$xpp$trialType==1), maxp2p=mxlfp, who='random', fextn=fltxtn)
      gc()
      
      if (OKEEG){
        getSessionAverage(subs.data(eeg,eeg$xpp$trialType==1), maxp2p=mxeeg, who='random', fextn=fltxtn)
        gc()
      }
      
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
      
      if(OKEEG){
        getSessionAverage(subs.data(eeg,eeg$xpp$trialType==2), maxp2p=mxeeg, who='regular', fextn=fltxtn)
        gc()
      }
      
      getSessionAverage(subs.data(csd,csd$xpp$trialType==2), maxp2p=mxcsd, who='regular', fextn=fltxtn)
      gc()
      
      getSessionAverage(subs.data(mua,mua$xpp$trialType==2), maxp2p=mxmua, who='regular', fextn=fltxtn)
      gc()
    }
  }          # calculate session averages; obsolete: see averager.R  
}