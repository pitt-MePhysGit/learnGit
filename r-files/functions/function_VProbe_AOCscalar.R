####
####20211025 - SK: put most of the file in a makeFig control; rnd and reg are lists not matrices, and seem to 
####break errorbar plots so I added the correct formulation (around line 215) - will need to be added to all other plots


function_VProbe_AOCscalar <- function(dataSet='g',animals=c('Sam'),who='All',cmpStr='MUA1c',dpthStr='MUA',useInd=NULL,df=NULL,
                                      cleanOnly=FALSE,valBin=NA,shft=NA,pthrISI=0.05, pthrAmp=0.05, makeFigs = TRUE){
  
  
  efpar <- par()
  
  baseDir <- '~/Dropbox/'
  source(paste(baseDir, 'RSkernel/r-files/functions/source_all_functions.R',sep=''))
  source( paste(baseDir,'ampOddClick/r-files/functions/functions_AOCdual.R',sep='') )
  source( paste(baseDir,'ampOddClick/r-files/functions/prepDFVprobe.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/expand2chdf.R',sep='') )
  
  source('~/Dropbox/r-compile/signalProcessoR/R/signalProcessoR.R')
  
  #colClasses <- c('character',rep('numeric',2),'character',rep('numeric',5)  ,rep('character',3),rep('numeric',3) )
  #df         <- read.table(paste(baseDir,'ampOddClick/ampOddclickVProbe_link.txt',sep=''),header=T,colClasses=colClasses)
  #df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))
  #df$exp     <- df$session
  #df$project <- 'ampOddClick'
  
  if(is.null(df)){
    df       <- prepDFsignalProcessor('ampOddClick','ampOddclickVProbe_link.txt')
    df$exp   <- df$session
  }
  
  ##  paths and such
  #rda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
  rda <- function(){'/Volumes/EEGCyno/RES/ampOddClick/RSscalar/'}
  trda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
  spd  <- function(){'/Volumes/Drobo5D3/EEG/EEGLab/ampOddClick/'}
  options(warn=0)
  
  ## ===============================================
  badDays <- c(1) 
  df$use  <- c(1)
  
  animalStr <- '_all'
  if (length(animals)==1)
    animalStr <- paste('_',animals,sep='')
  
  
  cleanStr  <- ''
  if (cleanOnly){
    df$use[badDays] <- 0
    cleanStr <- '_clean'
  }
  
  ## ====================================== subset recording sessions
  useInd <- which(df$drug=='N' &  df$dataSet==dataSet & df$animal%in%animals )#& df$use==1)# & df$probeString=='d')
  
  df         <- prepDFVprobe(df[useInd,],getTrTy=T)
  df <- df[which(!is.na(df$Nrand+df$Nrand)),]
  
  if (who == 'All')
    useInd <- which(df$Nrand > 200 & df$Nreg > 200)
  
  chdf   <- expand2chdf(df)
  
  ## ====================================== subset channels within session
  if(identical( substr(cmpStr,1,3),'CSD'))
    chdf <- chdf[chdf$useCSD,]
  
  # ================================
  uIdir <- paste('g', animalStr, cleanStr,'/', who,'/',cmpStr, '/', dpthStr,sep='')
  
  .GlobalEnv$pp   <- function(){paste(baseDir, '/ampOddClick/word/VProbe/plots/',uIdir,'/',sep='')}
  if (!file.exists(pp()))
    system( paste('mkdir -p ', pp()) )
  
  print(pp())
  
  df[useInd,]
  
  ## ================================= find depth
  layers     <- data.frame( cmpStrName='L2',      dmin=  200, dmax= 1000, stringsAsFactors=FALSE) 
  layers[2,] <- data.frame( cmpStrName='L3',      dmin= -700,  dmax=200, stringsAsFactors=FALSE)
  layers[3,] <- data.frame( cmpStrName='L4',      dmin= -1100,  dmax= -700, stringsAsFactors=FALSE)
  layers[4,] <- data.frame( cmpStrName='L5',      dmin= -1500,  dmax= -1100, stringsAsFactors=FALSE)
  layers[5,] <- data.frame( cmpStrName='L6',      dmin= -2000,  dmax= -1500, stringsAsFactors=FALSE)
  layers[6,] <- data.frame( cmpStrName='SGSo',    dmin=   40,  dmax= 500, stringsAsFactors=FALSE)
  layers[7,] <- data.frame( cmpStrName='SGSi',    dmin= -500,  dmax=  -40, stringsAsFactors=FALSE)
  layers[8,] <- data.frame( cmpStrName='GSi' ,    dmin= -1100,  dmax= -700, stringsAsFactors=FALSE)
  layers[9,] <- data.frame( cmpStrName='MUA' ,    dmin= -1500,  dmax= -300, stringsAsFactors=FALSE)
  layers[10,] <- data.frame( cmpStrName='SG' ,    dmin= -700,  dmax= 1000, stringsAsFactors=FALSE)
  layers[11,] <- data.frame( cmpStrName='IG' ,    dmin= -2000,  dmax= -1100, stringsAsFactors=FALSE)
  layers[12,] <- data.frame( cmpStrName='ThAf' ,  dmin= -2500,  dmax= -300, stringsAsFactors=FALSE)
  layers[13,] <- data.frame( cmpStrName='All' ,   dmin= -Inf,  dmax= +Inf, stringsAsFactors=FALSE)
  
  lnd     <- which( layers$cmpStrName %in% dpthStr )
  if( length(lnd)==1 )
    dpthRng <- c( unlist(layers[lnd, c('dmin','dmax')] ))
  
  if(length(lnd)==0)
    dpthRng <- c(-Inf, Inf)
  
  #chdf <- chdf[ which(chdf$depth > dpthRng[1] & chdf$depth < dpthRng[2]),]
  
  
  ## ==================================================
  readRSscalar <- function(ix,who='all',cmpStr='MUA1'){
    resFileName <- paste(rda(),chdf$exp[ix],'/RSscalar_ch',chdf$chNrCons[ix],'_',cmpStr,'_',who,'.rda',sep='')
    print(resFileName)
    load(resFileName)  
    return(resList)
  }
  
  ##if you want to know which sessions are failing; I'm leaving this on to exclude non-processed
  if(1==1){
    feList<-list()
    for(i in 1:dim(chdf)[1]){
      print(i)
      feList[[i]]<- file.exists(paste(rda(),chdf$exp[i],'/RSscalar_ch',chdf$chNrCons[i],'_',cmpStr,'_',who,'.rda',sep=''))
    }
    print('missing Sessions:')
    chdf$preprocessed <- unlist(feList)
    table(chdf$session, chdf$preprocessed)
    chdf <- chdf[which(unlist(feList)),]
  }
  
  #unique(chdf[errorc, 1])
  
  resList    <- lapply(1:dim(chdf)[1], readRSscalar, who=who, cmpStr=cmpStr)
  #tst <- readRSscalar(16, who=who, cmpStr=cmpStr)
  
  chdf$pValISI   <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$pValISI ) } )
  chdf$pValAmp   <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$pValAmp ) } )
  chdf$pValTrial <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$pValTrial ) } )
  chdf$pValTone  <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$pValTone ) } )
  
  chdf$pValDA    <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$pValDA ) } )
  chdf$betaISI   <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$betaISI ) } )
  chdf$betaTrial <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$betaTrial ) } )
  chdf$betaTone  <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$betaTone ) } )
  
  
  if (is.na(valBin[1])){
    print('calculating valBin')
    valBin           <- chdf$pValISI<pthrISI  | chdf$pValAmp<pthrAmp
  }
  
  valBinDepth <- chdf$depth > dpthRng[1] & chdf$depth < dpthRng[2]
  
  vlnd       <- which(valBin & valBinDepth)
  avgBetaISI <- mean(chdf$betaISI[vlnd])
  ##sgn         <- +1 #sign(avgBetaISI)       # remove once analyzeRSscalar has been rerun
  
  
  if(makeFigs){
    
    if (identical(who,'all')){
      if (length(vlnd)>2){
        ## ============================================================== MMN and ND  DVz
        tt        <- lapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVzbyToneTypeAndTrialType ) }) 
        ttMat     <- array(unlist(tt), c(dim(tt[[1]]),length(tt)) )
        
        xax                <- c(2,1)
        mnTT               <- apply(ttMat[,,vlnd], c(1,2), na.mean)
        seTT               <- apply(ttMat[,,vlnd], c(1,2), na.se)
        
        picFile <- paste('MMN_NDz.eps',sep='')
        eps.file(picFile, pdf=T, s=2, ratio=3/2)
        errorplot( xax, mnTT[,1], ey=seTT[,1],  type='b', lwd=2,pch=20, col='blue',ylim=c(-.2,.2),add=F,xlim=c(0.5,2.5))  # Rnd Std vs Dev
        errorplot( xax, mnTT[,2], ey=seTT[,2],  type='b', lwd=2,pch=20, col='black',add=T )                # Reg Std vs dEv
        dev.off()
        
        picFile <- paste('MMN_NDz2.eps',sep='')
        eps.file(picFile, pdf=T, s=2, ratio=3/2)
        errorplot( 1:3, mnTT[c(4,1,3)], ey=seTT[c(4,1,3)],  type='p', lwd=2,pch=20, col=c('black','blue','red'),
                   ylim=c(-.2,.2),xlim=c(0,4),add=F,wd=0.45)  # Rnd Std vs Dev
        dev.off()
        
        ## ============================================================== MMN and ND  DVai
        tt        <- lapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVaibyToneTypeAndTrialType ) }) 
        ttMat     <- array(unlist(tt), c(dim(tt[[1]]),length(tt)) )
        
        xax                <- c(2,1)
        mnTT               <- apply(ttMat[,,vlnd], c(1,2), na.mean)
        seTT               <- apply(ttMat[,,vlnd], c(1,2), na.se)
        
        picFile <- paste('MMN_NDai.eps',sep='')
        eps.file(picFile, pdf=T, s=2, ratio=3/2)
        errorplot( xax, mnTT[,1], ey=seTT[,1],  type='b', lwd=2,pch=20, col='blue',ylim=c(-.2,.2),add=F,xlim=c(0.5,2.5))  # Rnd Std vs Dev
        errorplot( xax, mnTT[,2], ey=seTT[,2],  type='b', lwd=2,pch=20, col='black',add=T )                # Reg Std vs dEv
        dev.off()
        
        picFile <- paste('MMN_NDai2.eps',sep='')
        eps.file(picFile, pdf=T, s=2, ratio=3/2)
        errorplot( 1:3, mnTT[c(4,1,3)], ey=seTT[c(4,1,3)],  type='p', lwd=2,pch=20, col=c('black','blue','red'),
                   ylim=c(-.2,.2),xlim=c(0,4),add=F,wd=0.45)  # Rnd Std vs Dev
        dev.off()
        
        
        
        ## ============================================================== MMN and ND  as a function of SOA
        tt        <- lapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVbyToneTypeTrialTypeAndSOA ) }) 
        ttMat     <- array(unlist(tt), c(dim(tt[[1]]),length(tt)) )
        
        xax                 <- as.numeric(unlist(dimnames(resList[[1]]$DVbyToneTypeTrialTypeAndSOA)[3] ))
        mnTT                <- apply(ttMat[,,,vlnd], c(1,2,3), na.mean)
        seTT                <- apply(ttMat[,,,vlnd], c(1,2,3), na.se)
        
        soand <- 1:6
        
        picFile <- paste('MMN_NDbySOA.eps',sep='')
        eps.file(picFile, pdf=T, s=2, ratio=3/2)
        errorplot( log2(xax[soand]), mnTT[2,1,soand], ey=seTT[1,1,soand],  type='b', lwd=2,pch=20, col='gray',ylim=c(-.6,.6),add=F)  # Std Rnd
        errorplot( log2(xax[soand]), mnTT[2,2,soand], ey=seTT[1,2,soand],  type='b', lwd=2,pch=20, col='black',add=T )                # Std Reg
        errorplot( log2(xax[soand]), mnTT[1,1,soand], ey=seTT[2,1,soand],  type='b', lwd=2,pch=20, col='blue',add=T)                  # Dev Rnd
        errorplot( log2(xax[soand]), mnTT[1,2,soand], ey=seTT[2,2,soand],  type='b', lwd=2,pch=20, col='red',add=T)                   # Dev Reg 
        dev.off()
        
        ## =================================================================== compare SOA dependence between regular and random condition
        rnd        <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVbyTrialTypeAndSOA[1,] ) }) 
        reg        <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVbyTrialTypeAndSOA[2,] ) }) 
        
        dsc <- sapply( 1:dim(chdf)[1], function(ix){ return( sum(dim( resList[[ix]]$DVbyTrialTypeAndSOA))) })
        
        xax                 <- 1:dim(resList[[1]]$DVbyTrialTypeAndSOA)[2] #as.numeric(unlist(dimnames(resList[[1]]$DVbyTrialTypeAndSOA)[2] ))
        mnReg               <- apply(reg[,vlnd], 1, na.mean)
        seReg               <- apply(reg[,vlnd], 1, na.se)
        
        mnRnd               <- apply(rnd[,vlnd], 1, na.mean)
        seRnd               <- apply(rnd[,vlnd], 1, na.se)
        
        picFile <- paste('TrialTypeAndSOA.eps',sep='')
        eps.file(picFile, pdf=T, s=2, ratio=3/2)
        errorplot( xax, mnReg[soand], ey=seReg[soand],  type='b', lwd=2,pch=20, col='black',ylim=c(-.6,.6), add=F )
        errorplot( xax, mnRnd[soand], ey=seRnd[soand],  type='b', lwd=2,pch=20, col='blue',add=T)
        dev.off()
        
        
        ## ============================================================================================================== check out deltaISIbin as a function of trial type
        rnd        <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVbyDeltaISIbinAndTrialType[,1] ) }) 
        reg        <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVbyDeltaISIbinAndTrialType[,2] ) }) 
        
        xaxlab              <- unlist(dimnames(resList[[1]]$DVbyDeltaISIbinAndTrialType)[[1]] )
        xax                 <- 1:length(xaxlab); xax <- xax-mean(xax)
        
        if(is.null(dim(reg))){
          mnReg               <- apply(matrix(unlist(reg[vlnd]), ncol=length(xax), byrow = TRUE), 2, na.mean)
          seReg               <- apply(matrix(unlist(reg[vlnd]), ncol=length(xax), byrow = TRUE), 2, na.se)
          
          mnRnd               <- apply(matrix(unlist(rnd[vlnd]), ncol=length(xax), byrow = TRUE), 2, na.mean)
          seRnd               <- apply(matrix(unlist(rnd[vlnd]), ncol=length(xax), byrow = TRUE), 2, na.se)
        }else{
          mnReg               <- apply(reg[,vlnd], 1, na.mean)
          seReg               <- apply(reg[,vlnd], 1, na.se)
          
          mnRnd               <- apply(rnd[,vlnd], 1, na.mean)
          seRnd               <- apply(rnd[,vlnd], 1, na.se)
        }
        
        picFile <- paste('DeltaISIAndTrialType.eps',sep='')
        eps.file(picFile, pdf=T, s=2, ratio=3/2)
        errorplot( xax, mnReg, ey=seReg,  type='b', lwd=2,pch=20, col='black',ylim=c(-.4,.4), add=F )
        errorplot( xax, mnRnd, ey=seRnd,  type='b', lwd=2,pch=20, col='blue',add=T)
        dev.off()
        
        
        ## =============================================================================================================== check out absDeltaISIbin as a function of trial type
        rnd        <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVbyAbsDeltaISIbinAndTrialType[,1] ) }) 
        reg        <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVbyAbsDeltaISIbinAndTrialType[,2] ) }) 
        
        xaxlab              <- unlist(dimnames(resList[[1]]$DVbyAbsDeltaISIbinAndTrialType)[[1]] )
        xax                 <- (0:(length(xaxlab)-1)) / 2
        mnReg               <- apply(reg[,vlnd], 1, na.mean)
        seReg               <- apply(reg[,vlnd], 1, na.se)
        
        mnRnd               <- apply(rnd[,vlnd], 1, na.mean)
        seRnd               <- apply(rnd[,vlnd], 1, na.se)
        
        picFile <- paste('absDeltaISIAndTrialType.eps',sep='')
        eps.file(picFile, pdf=T, s=2, ratio=3/2)
        errorplot( xax, mnReg, ey=seReg,  type='b', lwd=2,pch=20, col='black',ylim=c(-.4,.4), add=F )
        errorplot( xax, mnRnd, ey=seRnd,  type='b', lwd=2,pch=20, col='blue',add=T)
        dev.off()
        
      }# if vlnd>2
    }# if all
    
    ## ============================================ compare pitch selectivity between random and regular condition
    #chdf$bestPitchInd   <- sapply( 1:dim(chdf)[1], function(ix){ return(which( resList[[ix]]$betaPitch==max(resList[[ix]]$betaPitch))[1]) } )
    
    if (length(vlnd)>2) {
      
      
      ## ====================================================== Amp by SOA
      tmpsoa     <- lapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVbyAmpAndSOA ) }) 
      soa        <- array(unlist(tmpsoa), c(dim(tmpsoa[[1]]),length(tmpsoa)) )
      
      
      soaax                 <- (seq(log2(.2),by=1,length=6) + seq(log2(.4),by=1,length=6))/2
      dbax                  <- seq(62,by=6,length=5)  
      mnSOA               <- apply(soa[,,vlnd], c(1,2), na.mean)
      seSOA               <- apply(soa[,,vlnd], c(1,2), na.se)
      
      picFile <- paste('SOAandAmp.eps',sep='')
      eps.file(picFile, pdf=T, s=2, ratio=3/2)
      colList                <- soa.colors(6)
      xnd <- 1:5
      for (ix in 1:(dim(mnSOA)[2]))
        errorplot( dbax[xnd], mnSOA[xnd,ix], ey=seSOA[xnd,ix],  type='b', lwd=2,pch=20, col=colList[ix],ylim=c(-.5,1), add=ifelse(ix==1,F,T) )
      dev.off()
      
      picFile <- paste('SOAandAmp2.eps',sep='')
      eps.file(picFile, pdf=T, s=2, ratio=3/2)
      colList                <- bluered.colors(5)
      xnd <- 1:6
      for (ix in 1:(dim(mnSOA)[1]))
        errorplot( soaax[xnd], mnSOA[ix,xnd], ey=seSOA[ix,xnd],  type='b', lwd=2,pch=20, col=colList[ix],ylim=c(-.5,1), add=ifelse(ix==1,F,T) )
      dev.off()
      
      
      
      picFile <- paste('SOAandAmpImg.eps',sep='')
      eps.file(picFile, pdf=T, s=2, ratio=3/2)
      image(mnSOA%>=% -.5, zlim=c(-.4,.6), col=jet.colors(100), useRaster=T)
      dev.off()
      
      soand <- 1:6
      xnd <- 1:5
      fitISIslope <- function(pnd){
        res <- NA
        tdf <- data.frame( y=c(soa[pnd,soand,vlnd]), soand=rep(soand, length(vlnd)), unit=rep(1:length(vlnd),each=length(soand)) );
        
        try( lm0 <- lm(tdf$y ~ tdf$soand ) )
        try( res <- lm0$coeff[2])
        return(res)
      }
      
      picFile <- paste('SOAandPitchSlope.eps',sep='')
      eps.file(picFile, pdf=T, s=2, ratio=3/2)
      soaslope <- sapply( xnd, fitISIslope)
      plot(dbax[xnd],soaslope, type='b', pch=20)
      dev.off()
      
      
      tmpdAA     <- lapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVzbyDeltaAmpAndAmp ) }) 
      dAA        <- array(unlist(tmpdAA), c(dim(tmpdAA[[1]]),length(tmpdAA)) )
      
      #xax                 <- 1:5
      xax                 <- seq(-24,by=6,to=24)
      
      mndAA               <- apply(dAA[,,vlnd], c(1,2), na.mean)
      sedAA               <- apply(dAA[,,vlnd], c(1,2), na.se)
      
      picFile <- paste('deltaAmpandAmp.eps',sep='')
      eps.file(picFile, pdf=T, s=2, ratio=3/2)
      colList                <- bluered.colors(5)
      xnd <- 1:9
      for (ix in 1:(dim(mndAA)[2]))
        errorplot( xax[xnd], mndAA[xnd,ix], ey=sedAA[xnd,ix],  type='b', lwd=2,pch=20, col=colList[ix],ylim=c(-.5,1), add=ifelse(ix==1,F,T) )
      dev.off()
      
      
      
      if (identical(who,'all')){
        # ==================================================== RF by Reg and Rnd
        rnd       <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVbyTrialTypeAndAmp[1,] ) })
        reg       <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVbyTrialTypeAndAmp[2,] ) })
        
        
        xax                 <- 1:5
        mnReg               <- apply(reg[,vlnd], 1, na.mean)
        seReg               <- apply(reg[,vlnd], 1, na.se)
        
        mnRnd               <- apply(rnd[,vlnd], 1, na.mean)
        seRnd               <- apply(rnd[,vlnd], 1, na.se)
        
        picFile <- paste('TrialTypeAndAmp.eps',sep='')
        eps.file(picFile, pdf=T, s=2, ratio=3/2)
        errorplot( xax, mnReg, ey=seReg,  type='b', lwd=1.5,pch=20, col='black',ylim=c(-.6,.8), add=F,xlim=c(0,6) )
        errorplot( xax, mnRnd, ey=seRnd,  type='b', lwd=1.5,pch=20, col='blue',add=T)
        dev.off()
        
        dlta                 <- rnd - reg
        mnDlta               <- apply(dlta[,vlnd], 1, na.mean)
        seDlta               <- apply(dlta[,vlnd], 1, na.se)
        picFile <- paste('TrialTypeAndPitchDelta.eps',sep='')
        eps.file(picFile, pdf=T, s=2, ratio=3/2)
        errorplot( xax, mnDlta, ey=seDlta,  type='b', lwd=2,pch=20, col='black',ylim=c(-0.1,.5), add=F )
        abline(h=0);abline(v=0)
        dev.off()
        
        
        # ## ===================================================================== RF by ND and MMN
        # tmp     <- lapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVbyAmpAndSOA ) }) 
        # ndmmn   <- array(unlist(tmpsoa), c(dim(tmpsoa[[1]]),length(tmpsoa)) )
        # 
        # 
        # xax                 <- (-10:10)/4
        # mnNDMMN               <- apply(ndmmn[,,,vlnd], c(1,2,3), na.mean)
        # seNDMMN               <- apply(ndmmn[,,,vlnd], c(1,2,3), na.se)
        # 
        # picFile <- paste('MMN_NDbyRF.eps',sep='')
        # eps.file(picFile, pdf=T, s=2, ratio=3/2)
        # errorplot( xax, mnNDMMN[,1,1], ey=seNDMMN[,1,1],  type='b', lwd=2,pch=20, col='gray',ylim=c(-0.6,.8), add=F )
        # errorplot( xax, mnNDMMN[,1,2], ey=seNDMMN[,1,2],  type='b', lwd=2,pch=20, col='black', add=T )
        # errorplot( xax, mnNDMMN[,2,1], ey=seNDMMN[,2,1],  type='b', lwd=2,pch=20, col='blue', add=T )
        # errorplot( xax, mnNDMMN[,2,2], ey=seNDMMN[,2,2],  type='b', lwd=2,pch=20, col='red', add=T )
        # dev.off()
      }
      
      if (identical(who,'all') ){# | identical(who,'regular')){  
        # ## ===================================================================== RF by RepNr
        # xax                 <- (-10:10)/4
        # 
        # mnRepNr               <- apply(repnr[,,vlndp], c(1,2), na.mean)
        # seRepNr               <- apply(repnr[,,vlndp], c(1,2), na.se)
        # 
        # picFile <- paste('RFbyRepNr.eps',sep='')
        # eps.file(picFile, pdf=T, s=2, ratio=3/2)
        # colList                <- bluered.colors(6)
        # xnd <- 4:18
        # for (ix in 1:(dim(mnRepNr)[2]-0))
        #   errorplot( xax[xnd], mnRepNr[xnd,ix], ey=seRepNr[xnd,ix],  type='b', lwd=2,pch=20, col=colList[ix],ylim=c(-1,1), add=ifelse(ix==1,F,T) )
        # dev.off()
        # 
        # picFile <- paste('RFbyRepNrImg.eps',sep='')
        # eps.file(picFile, pdf=T, s=2, ratio=3/2)
        # image(mnRepNr[,6:1]%>=% -.5, zlim=c(-.5,.75), col=jet.colors(100), useRaster=T)
        # dev.off()
        # 
        # 
        # soand <- 1:6
        # fitISIslope <- function(pnd){
        #   res <- NA
        #   tdf <- data.frame( y=c(repnr[pnd,soand,vlndp]), soand=rep(soand, length(vlndp)), unit=rep(1:length(vlndp),each=length(soand)) );
        #   
        #   try( lm0 <- lm(tdf$y ~ tdf$soand ) )
        #   try( res <- lm0$coeff[2])
        #   return(res)
        # }
        # 
        # picFile <- paste('RFbyRepNrSlope.eps',sep='')
        # eps.file(picFile, pdf=T, s=2, ratio=3/2)
        # repnrslope <- sapply( xnd, fitISIslope)
        # plot(xax[xnd],repnrslope, type='b', pch=20)
        # dev.off()
        # 
        
      }
    }
    if (identical(who,'all') ){
      if(length(vlnd)>2){
        ## ============================================================================  DV by RepNr and SOA
        rep        <- lapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVbyRepNrAndSOA ) }) 
        repMat     <- array(unlist(rep), c(dim(rep[[1]]),length(rep)) )
        Nt         <- lapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$NbyRepNrAndSOA ) }) 
        NMat       <- array(unlist(Nt), c(dim(Nt[[1]]),length(Nt)) )
        repMat[which(NMat<20)] <- NA
        
        mnRep <- apply(repMat[,,vlnd], c(1,2), na.mean)
        seRep <- apply(repMat[,,vlnd], c(1,2), na.se)
        xax   <- 0:(dim(repMat)[1]-1)
        
        picFile <- paste('RepNrAndSOA.eps',sep='')
        eps.file(picFile, pdf=T, s=2, ratio=3/2)
        colList                <- soa.colors(length(xax))
        for (ix in 1:(dim(mnRep)[2]))
          errorplot( 1:length(xax), mnRep[,ix], ey=seRep[,ix],  type='b', lwd=2,pch=20, col=colList[ix],ylim=c(-.6,.6), add=ifelse(ix==1,F,T) )
        dev.off()
      }
    }
    ## =========================================================================== DV by absDeltaPitchOct
    adp        <- sapply( 1:dim(chdf)[1], function(ix){ return( resList[[ix]]$DVbyDeltaAmp ) }) 
    print(length(adp))
    print(length(adp[[1]]))
    xax                 <- -4:4
    mnADP               <- apply(adp[,vlnd], 1, na.mean)
    seADP               <- apply(adp[,vlnd], 1, na.se)
    
    picFile <- paste('absDeltaAmp.eps',sep='')
    eps.file(picFile, pdf=T, s=2, ratio=3/2)
    errorplot( xax, mnADP, ey=seADP,  type='b', lwd=2,pch=20, col='black',ylim=c(-.1,.15), add=F )
    text( xax[1], 0.4, labels=paste(length(vlnd),'/',sum(valBinDepth), '/', length(valBinDepth),sep=' '), adj=c(0,1) )
    dev.off()
    
    
    if (identical(who,'all')){
      xax                 <- -4:4
      
      tmp        <- lapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVbyAbsDeltaPitchAndTrialType ) }) 
      tmpMat     <- array(unlist(tmp), c(dim(tmp[[1]]),length(tmp)) )
      if( length(unlist(tmp)) != length(unlist(tmpMat)) )
        error('something is wrong with absDeltaAmpAndTrialType')
      
      mnADP               <- apply(tmpMat[,,vlnd], c(1,2), na.mean)
      seADP               <- apply(tmpMat[,,vlnd], c(1,2), na.se)
      
      
      picFile <- paste('absDeltaAmpByTrialType.eps',sep='')
      eps.file(picFile, pdf=T, s=2, ratio=3/2)
      errorplot( xax, mnADP[,1], ey=seADP[,1],  type='b', lwd=2,pch=20, col='blue',ylim=c(-.15,.15), add=F )
      errorplot( xax, mnADP[,2], ey=seADP[,2],  type='b', lwd=2,pch=20, col='black', add=T )
      abline(v=0)
      dev.off()
      
      
      # ========================== seqNrBin by TrialType
      tmp        <- lapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVbySeqNrBinAndTrialType ) }) 
      tmpMat     <- array(unlist(tmp), c(dim(tmp[[1]]),length(tmp)) )
      
      if( length(unlist(tmp)) != length(unlist(tmpMat)) )
        error('something is wrong with SeqNrBinAndTriallType')
      
      
      brks <- c(0,10,50,100,200,300,400,500,600,700)
      xax <- 1/2*(brks[1:(length(brks)-1)] + brks[2:length(brks)])
      soand <- 1:5
      
      mnSNR               <- apply(tmpMat[,,vlnd], c(1,2), na.mean)
      seSNR               <- apply(tmpMat[,,vlnd], c(1,2), na.se)
      
      picFile <- paste('SeqNrBinByTrialType.eps',sep='')
      eps.file(picFile, pdf=T, s=2, ratio=3/2)
      errorplot( xax[soand], mnSNR[soand,1], ey=seSNR[soand,1],  type='b', lwd=2,pch=20, col='blue',ylim=c(-.15,.15), add=F )
      errorplot( xax[soand], mnSNR[soand,2], ey=seSNR[soand,2],  type='b', lwd=2,pch=20, col='black', add=T )
      dev.off()
    }
    
  }
  
  chdf$valBin <- valBin
  return(chdf)
  
}

