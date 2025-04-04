getSessionAverage <- function(ttn,maxp2p=450,who='all',maxPowEOG=Inf,fextn='',logVar=F){
  ## a function that computes full-channel count averages
  ## for a number of main effect second order analyses
  
  # A) simple average across all trials                          mean
  # B) split by ISI bin                                          mean_ISIbin
  # C) split by tone number                                      mean_Tone
  # D) split by tone bin                                         mean_Tonebin
  # E) split by time in block                                    mean_SeqNr
  # F) split by tone type ('deviant' vs 'standard')              mean_ToneType
  
  # the variable 'who' is a string that specifies which data of the session was used for the analysis
  # 'all'     : all data from both trial types (regular and random)
  # 'random'  : only tones from random environment
  # 'regular' : only tones from regular environment
  # other extensions are possible depending on subsetting of the trials
  
  doAll <- TRUE
  doVar <- TRUE
  
  doISI           <- TRUE
  doAmp           <- TRUE
  doISIAmp        <- TRUE
  doSeqNr         <- TRUE
  doType          <- TRUE
  doDeltaAmp      <- TRUE
  doBlock         <- TRUE
  
  if (doAll){
    doISI           <- TRUE
    doAmp           <- TRUE
    doSeqNr         <- TRUE
    doType          <- TRUE
    doDeltaAmp      <- TRUE
    doBlock         <- TRUE
  }
  
  print(ttn$exp)
  
  ## check input
  if( length(dim(ttn$lfp))==2 )
    stop('only use data-sets with all channels')
  
  if( length(dim(ttn$lfp))==3 & dim(ttn$lfp)[3]<10  )
    stop('only use data-sets with all channels')
  
  ## define paths
  animal <- strsplit(ttn$exp,'_')[[1]][1]
  
  if (animal=='Rockey')
    eCogFile <- 'rockey'
  
  if (animal=='Jesse')
    eCogFile <- 'jessesym'
  
  if (animal=='Walter')
    eCogFile <- 'waltersym'
  
  if (animal=='Sam')
    eCogFile <- 'samsym'
  
  resFolder <- paste(rda(), ttn$exp,'/', sep='')
  if (!file.exists(resFolder))
    dir.create(resFolder, recursive=T)
  #system(paste('mkdir', resFolder))
  
  ## prepare data
  valInd <- which( ttn$xpp$p2p<maxp2p & ttn$xpp$powEOG<maxPowEOG& ttn$xpp$ISI>0.1 & ttn$xpp$ISI<13)
  if (logVar)
    valInd <- which( ttn$xpp$p2p<maxp2p & ttn$xpp$powEOG<maxPowEOG& ttn$xpp$ISI>0.1 & ttn$xpp$ISI<32)
  
  lfp    <- ttn$lfp[valInd,,]
  ## ======================== define helper function
  tmeanPSTH <- function(psth, FUN=mean){
    #print(dim(psth))
    ## little helper function that returns valid mean psths even
    ## if the number of repetitions is one or zero
    res <- array(NA,c(dim(lfp)[c(2,3)]))  
    
    if(is.array(psth)&length(dim(psth))==3 )
      res <- apply(psth,c(2,3),FUN)
    
    if( is.array(psth) & length(dim(psth))==2 & identical(FUN,mean))
      res <- psth
    
    return(res)
  }
  
  
  # define the sums of square function
  SSQFUN  <- function(x,n=NULL){  return(sum( (x-mean(x))^2 )) }
  #     if (is.null(n))
  #       ret <- sum( (x-mean(x))^2 )
  #     if (!is.null(n))
  #       ret <- sum( n* (x-mean(x))^2 )
  #     return( ret ) 
  #   }
  
  
  ## ==================================================== Compute and save the averages
  # A) simple average across all trials 
  print('Computing grand average')
  mn       <- apply(lfp, c(2,3), mean )
  Naverage <- dim(lfp)[1]
  xax      <- c(1)
  brks     <- NA
  
  save(mn, Naverage, xax, brks, file=paste(resFolder,who,'_mean_',ttn$what,fextn,'.rda',sep='') )
  
  if (doVar){
    ## =======================================================
    ## prepare the computation of time-resolved variance
    print('Computing time-resolved context-induced variance')
    
    # and remove baseline-correction
    rawlfp    <- subNormLFP(lfp, apply(lfp, c(1,3), mean) )
    
    mndotdot  <- apply(rawlfp, c(2,3), mean )  # calculate grand average across all conditions
    ssqTMat   <- apply(rawlfp, c(2,3), SSQFUN )  # calculate total sums of squares
    
    ssq <- ssqTMat
    vr  <- ssq/( dim(lfp)[1]-1 )
    
    save(ssq, Naverage, xax, brks, file=paste(resFolder,who,'_ssq_',ttn$what,fextn,'.rda',sep='') )
    save(vr, Naverage, xax, brks, file=paste(resFolder,who,'_var_',ttn$what,fextn,'.rda',sep='') )  
  }
  gc()
  
  if (doISI){
    # B) split by ISI bin
    rm(mn, xax, brks, Naverage, ssq, vr)
    print('Computing average ISIbin')
    extn <- 'ISIbin_'
    brks <- exp( log(2) * seq(log(0.2)/log(2),by=1/3,to=3.7) )
    
    #brks               <- exp( log(2) * seq(log(0.2)/log(2),by=.5,to=4.5) )
    brks[ 1 ]          <- 0.195
    brks[length(brks)] <- 13
    
    #if(logVar)
    #  brks <- exp( log(2) * seq(log(0.25)/log(2),by=1/2,to=5) )
    
    binKerISI       <- cut( ttn$xpp$ISI[valInd], breaks=brks )
    
    mn        <- tapply(1:dim(lfp)[1], list(binKerISI),function(ind){tmeanPSTH(lfp[ind,,]) })
    nullInd   <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
    for ( i in nullInd)
      mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
    
    
    Naverage        <- table(binKerISI)
    print(Naverage)
    xax             <- tapply(ttn$xpp$ISI[valInd], binKerISI, mean)
    names(mn)       <- xax
    save(mn, Naverage, xax, brks, file=paste(resFolder,who,'_mean_',extn,ttn$what,fextn,'.rda',sep=''))
    
    
    ## this function computes a matrix with the sums of squares for the Effect described in mn
    getSSQ    <- function(tmn,tNaverage){
      if(! (length(tmn)==length(tNaverage)) )
        browser()
      
      getSSQCnd <- function(n,mndd=mndotdot,Nav=tNaverage){ 
        res <- array(0, dim(mndotdot) )
        
        if (Nav[n]>0)
          res <- Nav[n]*(tmn[[n]]-mndd)^2
        
        return( res  )  
      }
      
      disc    <- lapply(1:length(tmn), getSSQCnd, mndd=mndotdot )
      dmpMat  <- array( unlist(disc), c(dim(disc[[1]]),length(disc) )  ) 
      ssqMat  <- apply(dmpMat,c(1,2), na.sum)
      return(ssqMat)
    }
    ## =========================================================================================
    if (doVar){
      mn       <- tapply(1:dim(lfp)[1], list(binKerISI),function(ind){tmeanPSTH(rawlfp[ind,,]) })
      nullInd  <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
      for ( i in nullInd)
        mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
      
      ssq   <- getSSQ(mn,Naverage)
      vr    <- ssq/(sum(Naverage)-length(Naverage[ which(Naverage>0) ]))
      eta   <- ssq/ssqTMat
      eta[which(is.na(eta)) ] <- 0;
      
      save(ssq , Naverage, xax, brks, file=paste(resFolder,who,'_ssq_',extn, ttn$what,fextn,'.rda',sep=''))
      save(vr  , Naverage, xax, brks, file=paste(resFolder,who,'_var_',extn,ttn$what,fextn,'.rda',sep=''))
      save(eta , Naverage, xax, brks, file=paste(resFolder,who,'_eta_',extn,ttn$what,fextn,'.rda',sep=''))  
    }
    
    # Bii) split by larger ISI bins, for visualization purposes ----------------------------------
    rm(mn, xax, brks, Naverage, ssq, vr,eta)
    gc()
    
    print('Computing average ISI large bins')
    extn <- 'ISIbinOct_'
    #brks <- exp( log(2) * seq(log(0.2)/log(2),by=1,to=3.7) )
    #brks[ 1 ]          <- 0.195
    #brks[length(brks)] <- 13
    
    brks <- exp( log(2) * seq(log2(0.250),by=1,to=4) )
    brks[ 1 ]          <- 0.245
    brks[length(brks)] <- 21
    
    #if(logVar)
    #  brks <- exp( log(2) * seq(log(0.25)/log(2),by=1,to=5) )
    
    
    binKerISI       <- cut( ttn$xpp$ISI[valInd], breaks=brks )
    
    mn       <- tapply(1:dim(lfp)[1], list(binKerISI),function(ind){tmeanPSTH(lfp[ind,,]) })
    nullInd         <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
    for ( i in nullInd)
      mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
    
    
    Naverage        <- table(binKerISI)
    print(Naverage)
    xax             <- tapply(ttn$xpp$ISI[valInd], binKerISI, mean)
    names(mn)       <- xax
    save(mn, Naverage, xax, brks, file=paste(resFolder,who,'_mean_',extn,ttn$what,fextn,'.rda',sep=''))
    
    if (doVar){
      # do variance
      mn       <- tapply(1:dim(lfp)[1], list(binKerISI),function(ind){tmeanPSTH(rawlfp[ind,,]) })
      nullInd         <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
      for ( i in nullInd)
        mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
      
      ssq   <- getSSQ(mn,Naverage)
      vr    <- ssq/(sum(Naverage)-length(Naverage[ which(Naverage>0) ]))
      eta   <- ssq/ssqTMat
      eta[which(is.na(eta)) ] <- 0;
      
      save(ssq , Naverage, xax, brks, file=paste(resFolder,who,'_ssq_',extn, ttn$what,fextn,'.rda',sep=''))
      save(vr  , Naverage, xax, brks, file=paste(resFolder,who,'_var_',extn,ttn$what,fextn,'.rda',sep=''))
      save(eta , Naverage, xax, brks, file=paste(resFolder,who,'_eta_',extn,ttn$what,fextn,'.rda',sep=''))  
    }
  }# end doISI
  
  if (doAmp){
    # C) split by tone number
    rm(mn, xax, brks, Naverage, ssq, vr, eta)
    gc()
    
    print('Computing average Tone')
    extn <- 'Amp_'
    ampInd       <- ttn$xpp$ampInd  #1 + (ttn$xpp$pitch-1)%%11
    
    brks           <- seq(0.5,by=1,to=5.5)
    xax            <- (brks[-1]+brks[-length(brks)])/2  
    binKerTone     <- cut( ampInd[valInd], breaks=brks )
    
    mn              <- tapply(1:dim(lfp)[1], list(binKerTone),function(ind){tmeanPSTH(lfp[ind,,]) })
    nullInd         <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
    for ( i in nullInd)
      mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
    
    
    Naverage        <- table(binKerTone)
    print(Naverage)
    xax             <- sort( unique(ttn$xpp$dB) )
    names(mn)       <- xax  
    brks            <- NA
    
    save(mn, Naverage, xax, brks, file=paste(resFolder,who,'_mean_',extn,ttn$what,fextn,'.rda',sep='') )
    
    if (doVar){
      # do the ssq and variance
      mn    <- tapply(1:dim(lfp)[1], list(binKerTone),function(ind){tmeanPSTH(rawlfp[ind,,]) })
      nullInd         <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
      for ( i in nullInd)
        mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
      
      ssq   <- getSSQ(mn,Naverage)
      vr    <- ssq/(sum(Naverage)-length(Naverage))
      eta   <- ssq/ssqTMat
      eta[which(is.na(eta)) ] <- 0;
      
      save(ssq , Naverage, xax, brks, file=paste(resFolder,who,'_ssq_',extn, ttn$what,fextn,'.rda',sep=''))
      save(vr  , Naverage, xax, brks, file=paste(resFolder,who,'_var_',extn,ttn$what,fextn,'.rda',sep=''))
      save(eta , Naverage, xax, brks, file=paste(resFolder,who,'_eta_',extn,ttn$what,fextn,'.rda',sep=''))  
      
      
      ## calculate sum of squares of the residuals within tones as an estimate of history-depndent variance
      # subtract the ssq explained by tones from total ssq
      ssq <- ssqTMat - ssq
      vr  <- ssq/(dim(lfp)[1])
      
      save(ssq, Naverage, xax, brks, file=paste(resFolder,who,'_ssqR_',ttn$what,fextn,'.rda',sep='') )
      save(vr, Naverage, xax, brks, file=paste(resFolder,who,'_varR_',ttn$what,fextn,'.rda',sep='') )  
    }
  }
  
  
  if (doISIAmp){
    # C) split by tone number
    rm(mn, xax, brks, Naverage, ssq, vr, eta)
    gc()
    
    print('Computing average Tone')
    extn <- 'ISIAmp_'
    
    # isi first
    brks               <- exp( log(2) * seq(log2(0.250),by=1,to=4) )
    brks[ 1 ]          <- 0.245
    brks[length(brks)] <- 21
    binKerISI          <- cut( ttn$xpp$ISI[valInd], breaks=brks )
    
    # amp second
    ampInd         <- ttn$xpp$ampInd  #1 + (ttn$xpp$pitch-1)%%11
    brks           <- seq(0.5,by=1,to=5.5)
    xax            <- (brks[-1]+brks[-length(brks)])/2  
    binKerTone     <- cut( ampInd[valInd], breaks=brks )
    
    
    mn              <- tapply(1:dim(lfp)[1], list(binKerISI,binKerTone),function(ind){tmeanPSTH(lfp[ind,,]) })
    nullInd         <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
    for ( i in nullInd)
      mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
    
    
    Naverage        <- table(binKerISI,binKerTone)
    print(Naverage)
    
    xax1                 <- tapply(ttn$xpp$ISI[valInd], binKerISI, mean)
    xax2                 <- sort( unique(ttn$xpp$dB) )
    dimnames(mn)[[1]]    <- xax1  
    dimnames(mn)[[2]]    <- xax2  
    brks                 <- NA
    
    save(mn, Naverage, xax, brks, file=paste(resFolder,who,'_mean_',extn,ttn$what,fextn,'.rda',sep='') )
    
    if (doVar){
      # do the ssq and variance
      mn    <- tapply(1:dim(lfp)[1], list(binKerISI,binKerTone),function(ind){tmeanPSTH(rawlfp[ind,,]) })
      nullInd         <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
      for ( i in nullInd)
        mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
      
      ssq   <- getSSQ(mn,Naverage)
      vr    <- ssq/(sum(Naverage)-length(Naverage))
      eta   <- ssq/ssqTMat
      eta[which(is.na(eta)) ] <- 0;
      
      save(ssq , Naverage, xax, brks, file=paste(resFolder,who,'_ssq_',extn, ttn$what,fextn,'.rda',sep=''))
      save(vr  , Naverage, xax, brks, file=paste(resFolder,who,'_var_',extn,ttn$what,fextn,'.rda',sep=''))
      save(eta , Naverage, xax, brks, file=paste(resFolder,who,'_eta_',extn,ttn$what,fextn,'.rda',sep=''))  
      
      
      ## calculate sum of squares of the residuals within tones as an estimate of history-depndent variance
      # subtract the ssq explained by tones from total ssq
      ssq <- ssqTMat - ssq
      vr  <- ssq/(dim(lfp)[1])
      
      save(ssq, Naverage, xax, brks, file=paste(resFolder,who,'_ssqR_',ttn$what,fextn,'.rda',sep='') )
      save(vr, Naverage, xax, brks, file=paste(resFolder,who,'_varR_',ttn$what,fextn,'.rda',sep='') )  
    }
  }
  
  
  
  
  
  if (doSeqNr){
    # E) split by time in block
    rm(mn, xax, brks, Naverage, ssq, vr, eta)
    gc()
    
    print('Computing average SeqNrbin')
    extn <- 'SeqNrbin_'
    ##brks       <- seq(0,by=5,to=600)
    #brks       <- c(0, 5, 10, 20, 30, 40, 50, 100, 200, 300, 600)
    brks       <- c(0, 25, 50, 100, 150, 200, 250, 300, 650)
    
    if (logVar)
      brks       <- seq(0, by=20, to=200)
    
    xax        <- (brks[-1]+brks[-length(brks)])/2
    
    binSeq     <- cut( ttn$xpp$seqNr[valInd], breaks=brks )
    
    ##mn                <- by(  1:length(binSeq), binSeq, function(nd){  apply( lfp[nd,,],c(2,3),mean) } )
    mn                <- tapply(1:dim(lfp)[1], list(binSeq),function(ind){tmeanPSTH(lfp[ind,,]) })
    nullInd         <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
    for ( i in nullInd)
      mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
    
    
    Naverage          <- table(binSeq)
    print(Naverage)
    
    #xax               <- brks[-length(brks)] + 2.5
    names(mn)         <- xax
    
    save(mn, Naverage, xax, brks, file=paste(resFolder,who,'_mean_',extn,ttn$what,fextn,'.rda',sep='') )
    
    if (doVar){
      # do the ssq and variance
      mn    <- tapply(1:dim(lfp)[1], list(binSeq),function(ind){tmeanPSTH(rawlfp[ind,,]) })
      nullInd         <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
      for ( i in nullInd)
        mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
      
      ssq   <- getSSQ(mn,Naverage)
      vr    <- ssq/(sum(Naverage)-length(Naverage))
      eta   <- ssq/ssqTMat
      eta[which(is.na(eta)) ] <- 0;
      
      save(ssq , Naverage, xax, brks, file=paste(resFolder,who,'_ssq_',extn, ttn$what,fextn,'.rda',sep=''))
      save(vr  , Naverage, xax, brks, file=paste(resFolder,who,'_var_',extn,ttn$what,fextn,'.rda',sep=''))
      save(eta , Naverage, xax, brks, file=paste(resFolder,who,'_eta_',extn,ttn$what,fextn,'.rda',sep=''))  
    }
  }
  if (doType){
    # F) split by tone type ('deviant' vs 'standard')
    rm(mn, xax, brks, Naverage, ssq, vr, eta)
    gc()
    
    print('Computing average by tone type')
    extn <- 'ToneType_'
    toneType        <- ttn$xpp$toneType[valInd]
    
    brks           <- seq(-0.5,by=1,to=1.5)
    xax            <- (brks[-1]+brks[-length(brks)])/2  
    toneType       <- cut( as.numeric(ttn$xpp$toneType[valInd]), breaks=brks )
    
    
    ##mn              <- by(  1:length(toneType), toneType, function(nd){  apply( lfp[nd,,],c(2,3),mean) } )
    mn              <- tapply(1:dim(lfp)[1], list(toneType),function(ind){tmeanPSTH(lfp[ind,,]) })
    nullInd         <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
    for ( i in nullInd)
      mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
    
    Naverage        <- table(toneType)
    print(Naverage)
    
    xax             <- 1:2
    brks            <- NA
    
    save(mn, Naverage, xax, brks, file=paste(resFolder,who,'_mean_',extn,ttn$what,fextn,'.rda',sep=''))
    
    if (doVar){
      # do the ssq and variance
      mn    <- tapply(1:dim(lfp)[1], list(toneType),function(ind){tmeanPSTH(rawlfp[ind,,]) })
      nullInd         <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
      for ( i in nullInd)
        mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
      
      
      ssq   <- getSSQ(mn,Naverage)
      vr    <- ssq/(sum(Naverage)-length(Naverage))
      eta   <- ssq/ssqTMat
      eta[which(is.na(eta)) ] <- 0;
      
      save(ssq , Naverage, xax, brks, file=paste(resFolder,who,'_ssq_',extn, ttn$what,fextn,'.rda',sep=''))
      save(vr  , Naverage, xax, brks, file=paste(resFolder,who,'_var_',extn,ttn$what,fextn,'.rda',sep=''))
      save(eta , Naverage, xax, brks, file=paste(resFolder,who,'_eta_',extn,ttn$what,fextn,'.rda',sep=''))  
    }
  }
  
  if (doDeltaAmp){
    # G) split by pitch difference
    rm(mn, xax, brks, Naverage, ssq, vr, eta)
    gc()
    
    print('Computing average by tone difference')
    extn <- 'binDeltaAmp_'
    #binKerTone      <- 1 + (ttn$xpp$pitch-1)%%11
    #ttn$xpp$deltaPitch <- c(NA,diff(binKerTone))
    
    brks           <- seq(-4.5,by=1,to=4.5)
    xax            <- (brks[-1]+brks[-length(brks)])/2  
    binDeltaPitch  <- cut( ttn$xpp$deltaAmpInd[valInd], breaks=brks )
    
    ##mn              <- by(  1:length(toneType), toneType, function(nd){  apply( lfp[nd,,],c(2,3),mean) } )
    mn              <- tapply(1:dim(lfp)[1], list(binDeltaPitch),function(ind){tmeanPSTH(lfp[ind,,]) })
    nullInd         <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
    for ( i in nullInd)
      mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
    
    
    Naverage        <- table(binDeltaPitch)
    
    save(mn, Naverage, xax, brks, file=paste(resFolder,who,'_mean_',extn,ttn$what,fextn,'.rda',sep=''))
    
    if (doVar){
      # do the ssq and variance
      mn    <- tapply(1:dim(lfp)[1], list(binDeltaPitch),function(ind){tmeanPSTH(rawlfp[ind,,]) })
      nullInd         <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
      for ( i in nullInd)
        mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
      
      
      ssq   <- getSSQ(mn,Naverage)
      vr    <- ssq/(sum(Naverage)-length(Naverage))
      eta   <- ssq/ssqTMat
      eta[which(is.na(eta)) ] <- 0;
      
      save(ssq , Naverage, xax, brks, file=paste(resFolder,who,'_ssq_',extn, ttn$what,fextn,'.rda',sep=''))
      save(vr  , Naverage, xax, brks, file=paste(resFolder,who,'_var_',extn,ttn$what,fextn,'.rda',sep=''))
      save(eta , Naverage, xax, brks, file=paste(resFolder,who,'_eta_',extn,ttn$what,fextn,'.rda',sep=''))  
    }
  }
  
  if (doBlock){
    # I) split by block
    rm(mn, xax, brks, Naverage, ssq, vr, eta)
    gc()
    
    print('Computing average block')
    extn <- 'Block_'
    
    brks       <- c(0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5)
    if(logVar)
      brks       <- c(0.5, 4.5, 8.5, 12.5)
    
    xax        <- (brks[-1]+brks[-length(brks)])/2
    
    binBlock     <- cut(ttn$xpp$trialNr[valInd], breaks=brks )
    
    mn              <- tapply(1:dim(lfp)[1], list(binBlock),function(ind){tmeanPSTH(lfp[ind,,]) })
    nullInd         <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
    for ( i in nullInd)
      mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
    
    Naverage          <- table(binBlock)
    print(Naverage)
    
    #xax               <- brks[-length(brks)] + 2.5
    names(mn)         <- xax
    
    save(mn, Naverage, xax, brks, file=paste(resFolder,who,'_mean_',extn,ttn$what,fextn,'.rda',sep='') )
    
    if (doVar){
      # do the ssq and variance
      mn    <- tapply(1:dim(lfp)[1], list(binBlock),function(ind){tmeanPSTH(rawlfp[ind,,]) })
      nullInd         <- which(sapply( 1:length(mn), function(n){is.null(mn[[n]])}))
      for ( i in nullInd)
        mn[[ i ]] <- array(NA,c(dim(lfp)[c(2,3)]))
      
      ssq   <- getSSQ(mn,Naverage)
      vr    <- ssq/(sum(Naverage)-length(Naverage) )
      eta   <- ssq/ssqTMat
      eta[which(is.na(eta)) ] <- 0;
      
      save(ssq , Naverage, xax, brks, file=paste(resFolder,who,'_ssq_',extn, ttn$what,fextn,'.rda',sep=''))
      save(vr  , Naverage, xax, brks, file=paste(resFolder,who,'_var_',extn,ttn$what,fextn,'.rda',sep=''))
      save(eta , Naverage, xax, brks, file=paste(resFolder,who,'_eta_',extn,ttn$what,fextn,'.rda',sep=''))  
    }
  }
}
