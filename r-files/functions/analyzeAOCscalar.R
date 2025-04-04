analyzeAOCscalar <- function(tone,ybin=NA,who='all',cmpStr='DV'){
  
  # built on the RSkernel version analyzeRSscalar
  # replaces pitchInd with ampInd
  # replaces deltaPitchInd with deltaAmpInd
  
  ## DV  : raw dependent variable 
  ## DVz : z-scored dependent variable
  ## DVn : variance-normalized residulals of DV ~ 1 + I(log(y$ISI)) 
  ## DVo : variance-normalized residulals of DV ~ 1 + I(log(y$ISI)) + pitchInd
  ## DVp : variance-normalized residulals of DV ~ 1                 + pitchInd
  ## DVq : variance-normalized residulals of DV ~ 1 + I(log(y$ISI)) + pitchInd + deltaPitchBin
  
  ## 1) difference in RF between random and regular                       DVn ~ pitchInd  | trialType   ( only for 'all' )
  ##                                                                      DVn ~ pitchInd
  ## 2) RF specificity as a function of SOA                               DVz ~ pitchInd  | SOA       
  ## 3) Difference between regular and random as a function of SOA        DVz ~ trialType | SOA         ( only for 'all' )
  ## 4) Difference in RF between standard and deviant                     DVn ~ pitchInd  | toneType    ( compare this measure between reg and rnd )
  ## 5) Difference between deviant and standard as a function of SOA      DVz ~ toneType  | SOA         ( compare this measure between reg and rnd )
  ## 6) Joint effect of RepNr and SOA                                     DVp ~ repNr     | SOA         ( only for 'all' and 'regular')
  ##                                                                      DVn ~ pitchInd  | repBlk
  ## 7) Difference between toneType * trialType by SOA                    DVz ~ toneType*trialType | SOA (only for 'all' )
  ## 8) Difference in RF between toneType * trialType                     DVn ~ pitchInd  | toneType*trialType    ( only for 'all' )
  ## 9) absDeltaPitch                                                     DVq ~ absDeltaPitchBin
  ##10) absDeltaPitch by trialType                                        DVq ~ absDeltaPitchBin  | trialType
  ##11) by ToneType and TrialType                                         DVznpq ~ toneType * trialType
  ##12) by SeqNrBlock (number of trials in block)                         DVn ~ seqNrBlck
  ##13) by deltaISI                                                       DVq ~ deltaISIbin
  ##                                                                      DVq ~ deltaISIbin | trialType    
  ##                                                                      DVq ~ absDeltaISIbin    
  ##                                                                      DVq ~ absDeltaISIbin | trialType    
  
  
  ######################
  ###20211021 - SK: lots of variables simply don't exist, I'm adding provisional versions of them here 
  # ##if seqNr does not exist in xpp; function fails
  ### 2023 03 08 - TT: all of these variables play key roles. you can't just replace them with dummy variables to prevent the function from crashing. 
  ###                   you need to make sure that they are created appropriately.
  # if(is.null(tone$xpp$seqNr)){
  #   tone$xpp$seqNr<-rep(1, dim(tone$xpp)[1])
  # }
  # if(is.null(tone$xpp$deltaAmpInd)){
  #   tone$xpp$deltaAmpInd<-c(0, diff(tone$xpp$ampInd))
  # }
  # if(is.null(tone$xpp$toneType)){
  #   tone$xpp$toneType<-rep(1, dim(tone$xpp)[1])
  # }
  # if(is.null(tone$xpp$deltaISIbin)){
  #   tone$xpp$deltaISIbin<-c(0, diff(tone$xpp$ISI))
  # }
  # if(is.null(tone$xpp$absDeltaISIbin)){
  #   tone$xpp$absDeltaISIbin<-abs(tone$xpp$ISI)
  # }
  # if(is.null(tone$xpp$repNr)){
  #   tone$xpp$repNr<-rep(1, dim(tone$xpp)[1])
  # }
  # 
  
  # backwards compatibility with switch in naming convention from ISI to ICI
  if(is.null(tone$xpp$ISI))
    tone$xpp$ISI <- tone$xpp$ICI
  
  if(is.null(tone$xpp$deltaISIbin))
    tone$xpp$deltaISIbin <- tone$xpp$deltaICIbin
  
  if(is.null(tone$xpp$absDeltaISIbin))
    tone$xpp$absDeltaISIbin <- tone$xpp$absDeltaICIbin
  
  
  ## preliminaries
  tone$xpp$DV <- tone$xpp[[cmpStr]]
  resList     <- list()
  
  tone$xpp$seqNrBin <- cut(tone$xpp$seqNr, breaks=c(0,10,50,100,200,300,400,500,600,2000) )
  tone$xpp$index    <- 1:length(tone$xpp$seqNr)
  tone$xpp$index    <- tone$xpp$index - mean(tone$xpp$index)
  
  #resDir <- paste(pp(),tone$animal,'/FLA/',sep='')
  #if (!file.exists(resDir))
  #  system(  paste('mkdir -p ', resDir, sep=''), intern=T )
  
  #resDir1 <- paste(rda(),tone$exp,'/rda/RSscalar/',sep='')
  #if (!file.exists(resDir1))
  #  system(  paste('mkdir -p ', resDir1, sep=''), intern=T )
  
  resPath      <- paste( '/Volumes/EEGCyno/RES/ampOddClick/RSscalar/', tone$exp,'/',sep='')
  if (tone$dfline$probe>1)
    resPath      <- paste( '/Volumes/EEGCyno/RES/ampOddClick/RSscalar/', tone$exp, '_',tone$dfline$probe,'/',sep='')
  
  if (!file.exists(resPath))
    system(  paste('mkdir -p ', resPath, sep=''), intern=T )
    
  
  #extn        <- paste(tone$animal,'/FLA/NoModel_',tone$exp,'_',cmpStr, sep='')
  #resFileName <- paste(tone$exp,'/RSscalar_',cmpStr,'_',who,'.rda',sep='')
  #resFileName <- paste('/RSscalar/RSscalar_',cmpStr,'_',who,'.rda',sep='')
  #resFileName <- paste(tone$exp,'/rda/RSscalar/RSscalar_',cmpStr,'_',who,'.rda',sep='')
  resFileName <- paste('/RSscalar_',cmpStr,'_',who,'.rda',sep='')
  
  
  # no subsetting within the routine. Assume that subsetting has already been taken care of before calling this function
  
  ##SK: this was erroring; had to add any() and length(ybin)==0
  if (any(is.na(ybin)) | length(ybin)==0)
    ybin <- rep(T,length(tone$xpp$DV))
    
  ybin[1:3]     <- FALSE
  y                <- subset(tone$xpp, ybin)
 
  # mean and standard deviation of average mua
  resList$mnDV  <- mean(y$DV)
  resList$sdDV  <- sd(y$DV)
  resList$DV <- y$DV

  ## ======================================= linear models
  lm0   <- lm(y$DV ~  1  + y$index                  )
  lm1a  <- lm(y$DV ~  1  + y$index + I(log(y$ISI))  )
  lm1b  <- lm(y$DV ~  1  + y$index                + y$ampIndFactor  )
  lm2   <- lm(y$DV ~ -1  + y$index + I(log(y$ISI)) + y$ampIndFactor )
  lm3   <- lm(y$DV ~ -1  + y$index + I(log(y$ISI)) + y$ampIndFactor + y$deltaAmpInd)
  
  
  # ======================================== the normalized dependent variables  
  y$DVz  <- (y$DV-mean(y$DV))/sd(y$DV)
  y$DVi  <- lm1a$residuals/sd(lm1a$residuals)   # regress out logISI
  y$DVa  <- lm1b$residuals/sd(lm1b$residuals)   # regress out ampInd
  y$DVai <- lm2$residuals/sd(lm2$residuals)     # regress out logISI, ampInd
  y$DVq  <- lm3$residuals/sd(lm3$residuals)     # regress out logISI, ampInd, deltaAmpInd
  
  # ======================================== 
  tmp             <- anova(lm1b,lm2)
  resList$pValISI <- tmp$`Pr(>F)`[2]
  resList$betaISI <- lm2$coefficients[2]
  
  tmp               <- anova(lm1a,lm2)
  resList$pValAmp <- tmp$`Pr(>F)`[2]
  resList$betaAmp <- lm2$coefficients[3:length(lm2$coefficients)]
  
  tmp               <- anova(lm2,lm3)
  resList$pValDA <- tmp$`Pr(>F)`[2]
  resList$betaDA <- lm3$coefficients[length(lm3$coefficients)]
  
  if (identical(who,'all') ){
    lm2ToTr      <- lm(y$DV ~ -1 + y$index + I(log(y$ISI)) + y$ampIndFactor + y$toneType + y$trialType)
    lm2To        <- lm(y$DV ~ -1 + y$index + I(log(y$ISI)) + y$ampIndFactor + y$toneType              )
    lm2Tr        <- lm(y$DV ~ -1 + y$index + I(log(y$ISI)) + y$ampIndFactor              + y$trialType)
    
    tmp                    <- anova(lm2To, lm2ToTr)
    resList$pValTrialType  <- tmp$`Pr(>F)`[2]
    resList$betaTrialType  <- lm2ToTr$coefficients[length(lm2ToTr$coefficients)]
    
    tmp                   <- anova(lm2Tr, lm2ToTr)
    resList$pValToneType  <- tmp$`Pr(>F)`[2]
    resList$betaToneType  <- lm2ToTr$coefficients[ length(lm2ToTr$coefficients)-1 ]
  }
  
  if (!identical(who,'all') ){
    lm2To      <- lm(y$DV ~ -1 + y$index + I(log(y$ISI)) + y$ampIndFactor + y$toneType )
    
    tmp                   <- anova(lm2, lm2To)
    resList$pValToneType  <- tmp$`Pr(>F)`[2]
    resList$betaToneType  <- lm2To$coefficients[ length(lm2To$coefficients) ]
  }
  
  
  
  ## =====================================================================================================================================
  ## 1) effect of ampInd between random and regular                       DVi ~ ampInd  (| trialType)   ( only for 'all' )
  resList$DVbyPitch               <- tapply(y$DVi, list(y$ampInd), mean)
  resList$sdbyPitch               <- tapply(y$DVi, list(y$ampInd), na.se)
  
  if (identical(who,'all' )){
    resList$DVbyPitchandTrialType <- tapply( y$DVi,          list(y$ampInd,y$trialType),  mean)
    resList$DVbyPitchandTrialType <- tapply( y$DVi,          list(y$ampInd,y$trialType),  na.se)
  }
  
  ## 2) Effect of SOA
  resList$DVbySOA <- tapply( y$DVa,          list(y$isiOct),  mean)
  resList$sdbySOA <- tapply( y$DVa,          list(y$isiOct),  na.se)
  
  
  ## 3) effect of ampInd as a function of SOA
  resList$DVbyAmpAndSOA <- tapply( y$DVz,          list(y$ampInd,y$isiOct),  mean)
  resList$sdbyAmpAndSOA <- tapply( y$DVz,          list(y$ampInd,y$isiOct),  na.se)
  

  ## 4) Effect of SOA by trialType        DVa ~ trialType (| SOA)         ( only for 'all' )
  if (identical(who,'all' )){
    resList$DVbyTrialTypeAndSOA <- tapply( y$DVa,          list(y$trialType,y$isiOct),  mean)
    resList$sdbyTrialTypeAndSOA <- tapply( y$DVa,          list(y$trialType,y$isiOct),  na.se)
  }
  
  
  ## 5) Effect of amp by trialType                     DVn ~ pitchInd  | toneType    ( compare this measure between reg and rnd )
  resList$DVbyAmp <- tapply( y$DVi,          list(y$ampInd),  mean)
  resList$sdbyAmp <- tapply( y$DVi,          list(y$ampInd),  na.se)
  
  if (identical(who,'all' )){
    resList$DVbyTrialTypeAndAmp <- tapply( y$DVi,          list(y$trialType,y$ampInd),  mean)
    resList$sdbyTrialTypeAndAmp <- tapply( y$DVi,          list(y$trialType,y$ampInd),  na.se)
  }
  
  ## 6) Effect of ISI and Amp by toneType
  resList$DVbyToneTypeAndSOA <- tapply( y$DVa,          list(y$toneType,y$isiOct),  mean)
  resList$sdbyToneTypeAndSOA <- tapply( y$DVa,          list(y$toneType,y$isiOct),  na.se)
  
  resList$DVbyToneTypeAndAmp <- tapply( y$DVi,          list(y$toneType,y$ampInd),  mean)
  resList$sdbyToneTypeAndAmp <- tapply( y$DVi,          list(y$toneType,y$ampInd),  na.se)
  
  

  ## 7) Joint effect of RepNr and SOA                                     DVz ~ repNr     | SOA         ( only for 'all' )
  if (identical(who,'all' ) | identical(who,'regular')){
    brks   <- c(-0.5, 0.5, 1.5, 2.5, 6.5, 18.8, 100 )
    repBlk <- cut(y$repNr, breaks=brks)
    
    xax                 <- tapply(c(y$repNr), repBlk, mean)
    yax                 <- tapply(y$ISI, y$isiOct, mean)
    mny                 <- tapply(y$DVa, list(repBlk, y$isiOct), mean )
    sey                 <- tapply(y$DVa, list(repBlk, y$isiOct), na.se )
    Ny                  <- tapply(y$DVa, list(repBlk, y$isiOct), length )
    
    dimnames(mny)[[1]] <- xax
    dimnames(mny)[[2]] <- yax
    
    indSOA        <- tapply(y$DVa, list(y$repNr,y$isiOct))
    ZeroIndex     <- seq(1,by=length(table(y$repNr)), length=length(table(y$isiOct)))
    mnyminus1     <- rep(NA, length(table(y$isiOct)) ) 
    seyminus1     <- mnyminus1
    for( i in 1:length(mnyminus1)){
      mnyminus1[i]     <- mean(y$DVa[ which(indSOA==ZeroIndex[i])-1 ] )
      seyminus1[i]     <- na.se(y$DVa[ which(indSOA==ZeroIndex[i])-1 ] )
    }
    mnyadd         <- rbind(mnyminus1, mny)
    seyadd         <- rbind(seyminus1, sey)
    Nyadd          <- rbind(Ny[1,],Ny ) 
    xaxadd         <- c(-1,xax)
    
    resList$DVbyRepNrAndSOA <- mnyadd
    resList$NbyRepNrAndSOA  <- Nyadd
    
    #colList                <- soa.colors(length(yax))
    #for (ix in 1:(dim(mny)[2]-1))
    #  errorplot( 1:length(xaxadd), mnyadd[,ix], ey=seyadd[,ix],  type='b', lwd=2,pch=20, col=colList[ix],ylim=c(-1,1), add=ifelse(ix==1,F,T) )
   
    ## ===============================================================  effect of repetition Nr on ampInd
    resList$DVbyAmpAndRepNr <- tapply(y$DVi, list(y$ampInd, repBlk), mean )
    
  }
  
  ## 7) Difference between toneType * trialType by SOA                    DVz ~ toneType*trialType | SOA (only for 'all' )
  if (identical(who,'all')){
    xax               <- tapply(c(y$toneType),  y$toneType, mean)
    yax               <- tapply(c(y$trialType), y$trialType, mean)
    zax               <- tapply(y$ISI, y$isiOct, mean)
    
    mny               <- tapply(y$DVz, list(y$toneType,y$trialType,y$isiOct), mean)
    sey               <- tapply(y$DVz, list(y$toneType,y$trialType,y$isiOct), na.se)
    
    dimnames(mny)[[1]] <- xax
    dimnames(mny)[[2]] <- yax
    dimnames(mny)[[3]] <- zax
    
    resList$DVbyToneTypeTrialTypeAndSOA <- mny
  }
  
  
  ## 8) Difference in RF between toneType * trialType                     DVn ~ pitchInd  | toneType*trialType    ( only for 'all' )
  if (identical(who,'all')){
    xax               <- tapply(y$ampInd, y$ampIndFactor, mean)
    yax               <- tapply(c(y$toneType), y$toneType, mean)
    zax               <- tapply(c(y$trialType), y$trialType, mean)
    
    mny               <- tapply(y$DVi, list(y$ampInd,y$toneType,y$trialType), mean)
    sey               <- tapply(y$DVi, list(y$ampInd,y$toneType,y$trialType), na.se)
    
    dimnames(mny)[[1]] <- xax
    dimnames(mny)[[2]] <- yax
    dimnames(mny)[[3]] <- zax
    
    resList$DVbyAmpToneTypeAndTrialType <- mny
  }
  
  ## 9) deltaAmp                                                     DVq ~ absDeltaPitchBin
  xax               <- tapply(y$deltaAmp, y$deltaAmpInd, mean)
  
  mny               <- tapply(y$DVai, list(y$deltaAmpInd), mean)
  sey               <- tapply(y$DVai, list(y$deltaAmpInd), na.se)
  
  dimnames(mny)[[1]] <- xax
  
  resList$DVbyDeltaAmp <- mny
  
  
  ## 9b) effect of deltaAmp by ampInd
  resList$DVzbyDeltaAmpAndAmp <- tapply( y$DVz,          list(y$deltaAmpInd,y$ampInd),  mean)
  resList$sdzbyDeltaAmpAndAmp <- tapply( y$DVz,          list(y$deltaAmpInd,y$ampInd),  na.se)
  
  
  
  ##10) absDeltaPitch by trialType                                        DVq ~ absDeltaPitchBin  | trialType
  if (identical(who,'all')){
    xax               <- tapply(y$deltaAmp, y$deltaAmpInd, mean)
    yax               <- tapply(c(y$trialType), y$trialType, mean)
    
    mny               <- tapply(y$DVai, list(y$deltaAmpInd,y$trialType), mean)
    sey               <- tapply(y$DVai, list(y$deltaAmpInd,y$trialType), na.se)
    
    dimnames(mny)[[1]] <- xax
    dimnames(mny)[[2]] <- yax
    
    resList$DVbyAbsDeltaPitchAndTrialType <- mny
  }
  
  
  ##11) by ToneType and TrialType                                         DVznpq ~ toneType * trialType
  if (identical(who,'all')){
    
    resList$DVbyToneTypeAndTrialType  <- tapply(y$DV,  list(y$toneType,y$trialType), mean)
    resList$DVzbyToneTypeAndTrialType <- tapply(y$DVz, list(y$toneType,y$trialType), mean)
    resList$DVabyToneTypeAndTrialType <- tapply(y$DVa, list(y$toneType,y$trialType), mean)
    resList$DVibyToneTypeAndTrialType <- tapply(y$DVi, list(y$toneType,y$trialType), mean)
    resList$DVaibyToneTypeAndTrialType <- tapply(y$DVai, list(y$toneType,y$trialType), mean)
  }
  
  #12) by SeqNrBin
  resList$DVbySeqNrBin  <- tapply(y$DVai,  list(y$seqNrBin), mean)
  if (identical(who,'all'))
    resList$DVbySeqNrBinAndTrialType  <- tapply(y$DVai,  list(y$seqNrBin,y$trialType), mean)
  
  
  
  # 13) by diffISI
  resList$DVbyDeltaISIbin              <- tapply(y$DVai,  list(y$deltaISIbin), mean)
  resList$DVbyAbsDeltaISIbin           <- tapply(y$DVai,  list(y$absDeltaISIbin), mean)
  if (identical(who,'all')){
    resList$DVbyDeltaISIbinAndTrialType    <- tapply(y$DVai,  list(y$deltaISIbin, y$trialType), mean)
    resList$DVbyAbsDeltaISIbinAndTrialType <- tapply(y$DVai,  list(y$absDeltaISIbin,y$trialType), mean)
  }
  
  save(file=paste(resPath,resFileName,sep=''), resList)
  
  return(resList)
}
