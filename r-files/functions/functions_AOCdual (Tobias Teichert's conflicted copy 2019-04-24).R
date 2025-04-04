loadAOC <- function( dfline, chStr='frontoCentral',showCh=1,show=F,exp=NA, from=NA, to=NA, what='llp2', suaMethod='nev',
                    minFrq=800,dFrq=0.21783, trda=NA, extn='',xppextn='',blcorr=T, suarng=c(-1,1), byBlock=T,loadCC=FALSE){
  
  readSUA <- exists('cell',dfline) & exists('channel',dfline)
  
  if(!is.na(trda))
    rda <- trda
  
  print('loadRS')
  print(rda())
  
  if(!is.na(dfline)){
    exp    <- dfline$session
    from   <- ifelse(is.na(dfline$minTR),NA, dfline$minTR)
    to     <- ifelse(is.na(dfline$maxTR),NA, dfline$maxTR)
    minFrq <- dfline$minOct
    dFrq  <- dfline$stepOct
  }
  
  tone       <- list()
  tone$xpp   <- read.table( paste(rda(),'/',exp,'/', exp, xppextn,'.txt',sep=''), header=T ) 
  if(nchar(extn)>0)
    print(extn)
  
  raw   <- identical(what,'raw')
  llp   <- identical(what,'llp')
  mua   <- identical(what,'mua')
  mbp   <- identical(what,'mbp')
  bsp   <- identical(what,'bsp')
  
  chextn <- ''
  if (llp)
    chextn <- extn
  if (raw)
    chextn <- paste('_raw',extn,sep='')
  if (mua)
    chextn <- paste('_mua',extn,sep='')
  
  tone$taxis <- NULL
  tone$lfp   <- NULL
  
  if (!is.na(chStr) ){# get lfps only if chStr is specified
    
    if (!loadCC){
      print(paste(rda(),exp,'/',exp,'_',chStr[1],chextn,'.mat',sep=''))
      tmpLFP   <- t( readMat( paste(rda(), exp, '/', exp,'_',chStr[1],chextn,'.mat',sep=''))[[1]] )
    }
    
    if (loadCC){ # new naming convention introduced for set2R and cross-correlation. Ultimately transition to this convention
      print(paste(rda(), exp, '/', exp,'_',what,'_',chStr[1],'.mat',sep=''))
      tmpLFP   <- t( readMat(paste(rda(), exp, '/', exp,'_',what,'_',chStr[1],'.mat',sep=''))[[1]] )
    }
    
    tone$lfp <- array(NA, c(dim(tmpLFP)[1],dim(tmpLFP)[2], length(chStr)) )
    tone$lfp[,,1] <- tmpLFP 
    if(length(chStr)>1){
      for (chnd in 2:length(chStr)){
        if (!loadCC){
          print(paste(rda(),exp,'_',chStr[chnd],chextn,'.mat',sep=''))
          tone$lfp[,,chnd]   <- t( readMat( paste(rda(), exp,'/', exp,'_',chStr[chnd],chextn,'.mat',sep=''))[[1]] )
        }
        
        if (loadCC){
          print(paste(rda(), exp, '/', exp,'_',what,'_',chStr[chnd],'.mat',sep=''))
          #browser()
          tone$lfp[,,chnd]   <- t( readMat( paste(rda(), exp, '/', exp,'_',what,'_',chStr[chnd],'.mat',sep='') )[[1]] )
        }
        
      }
    }
    
    textn <- ''
    if (raw)
      textn <- '_raw'
    if (mua)
      textn <- '_mua' 
    
    if (!loadCC){
      print(paste(rda(), exp,'/', exp,'_taxis',textn,'.mat',sep=''))
      tone$taxis <- readMat( paste(rda(),exp,'/', exp,'_taxis',textn,'.mat',sep=''))[[1]]
    }
    
    if (loadCC){
      print(paste(rda(), exp, '/', exp,'_',what,'_taxis.mat',sep=''))
      #browser()
      tone$taxis   <- readMat( paste(rda(), exp, '/', exp,'_',what,'_taxis.mat',sep='') )[[1]]
    }
    
    
    # non-baseline-corrected lfp
    lfpnbl <- tone$lfp
    
    tmp            <- baseLineCorrection(tone$lfp,c(-50,0), tone$taxis)
    if (blcorr==1){
      tone$lfp       <- tmp$res
      tone$xpp$blcor <- tmp$correction
    }
  }
  
  shp <- NULL
  if (readSUA){
    
    if (suaMethod=='nev'){
      spktimeN  <- read.table(  paste( spd(), dfline$session,'/ch_', dfline$channel,'/',dfline$cell,'_nev_',dfline$fextn,'_time.txt',sep=''))$V1
      spktimeM  <- spktimeN#read.table(  paste( spd(), dfline$session,'/ch_', dfline$channel,'/time',ifelse(is.na(dfline$thrsh),'',paste('_',dfline$thrsh,sep='')),'_nev.txt',sep=''))[[1]]
      
      
      shapeFile <- paste( spd(), dfline$session,'/ch_', dfline$channel,'/',dfline$cell,'_nev_',dfline$fextn,'_shape.mat',sep='')
      if (file.exists(shapeFile))
        shp       <- readMat(shapeFile)$shp
      #browser()
    }
    
    # read in the event Time file
    eventFile     <-  paste( spd(), dfline$session, '/events.txt',sep='');
    events        <-  read.table(  eventFile )
    
    spk <- list()
    mua <- list()
    for (tnd in 1:dim(events)[1]){
      relSpkTimeN <- spktimeN - events$V2[tnd]
      spk[[tnd]] <- relSpkTimeN[ which(relSpkTimeN>suarng[1] & relSpkTimeN<suarng[2]) ]
      
      relSpkTimeM <- spktimeM - events$V2[tnd]
      mua[[tnd]] <- relSpkTimeM[ which(relSpkTimeM>suarng[1] & relSpkTimeM<suarng[2]) ]
    }
    tone$spk <- spk
    tone$mua <- mua
  }
  
  
  tone$xpp$exp <- rep(exp,length(tone$xpp$seqNr))
  
  ## write out some info about the data contained =========================
  tone$exp     <- exp
  tone$chStr   <- chStr
  tone$what    <- what
  
  ## prepare xpp =====================================
  amplitude              <- seq(from=62,by=6,length=5) 
  tone$xpp$dB            <- amplitude[tone$xpp$ampInd]
  
  pastAmpInd             <- my.lag(tone$xpp$ampInd,1)
  pastAmpInd[1]          <- NA
  tone$xpp$deltaAmpInd   <- tone$xpp$ampInd - pastAmpInd 
  
  #pastISI               <- my.lag(tone$xpp$ISI,1)
  #pastISI[c(1,2)]       <- NA
  #tone$xpp$deltaISI     <- tone$xpp$ISI - pastISI 
  #tone$xpp$deltaISIfrac <- (pastISI+tone$xpp$deltaISI)/(pastISI)
  
  tone$xpp$ISI[ which( tone$xpp$ISI<0)]  <- 20
  
  tone$xpp$deltaISI     <- c(NA, diff(log2(tone$xpp$ISI)))
  tone$xpp$absDeltaISI  <- abs(tone$xpp$deltaISI)
  
  brks    <- seq(-4.5,by=1,4.5)
  brks[1] <- -10; brks[length(brks)]<- 10 
  tone$xpp$deltaISIbin  <- cut(tone$xpp$deltaISI,  breaks=brks)
  
  brks    <- seq(-0.25,by=0.5,4.5)
  brks[length(brks)]<- 10 
  tone$xpp$absDeltaISIbin  <- cut(tone$xpp$absDeltaISI, breaks=brks)
  
  tone$xpp$ampDev  <- abs(tone$xpp$deltaAmpInd )     >  0       
  tone$xpp$timeDev <- abs(tone$xpp$deltaISI)         >= 0.1  # time Deviant if ISI differs by more than 5% from previous ISI
  
  tttmp    <- tone$xpp$ampDev | tone$xpp$timeDev
  tttmp[1] <- TRUE
  
  tone$xpp$toneType <- as.numeric(!tttmp)
  
  rpCnt             <- as.numeric(tone$xpp$toneType)==1
  tmp               <- get.epoch( as.numeric(tone$xpp$toneType)==1  )
  tmp$dur           <- tmp$offset - tmp$onset
  rpCnt[tmp$offset] <- -tmp$dur  
  tone$xpp$repNr    <- cumsum(rpCnt)
  
  
  ## subset range of valid tones if specified ======================
  if (byBlock)
    trialRng <- na.range(c(tone$xpp$trialNr))
  
  if (!byBlock)
    trialRng <- c(1,length(tone$xpp$trialNr))
  
  if(!is.na(from))
    trialRng[1] <- from
  
  if(!is.na(to))
    trialRng[2] <- to
  
  if (!byBlock)
    tone <- subs.data(tone, 1:length(tone$xpp$trialNr) %in% trialRng[1]:trialRng[2])
  
  if (byBlock)
    tone <- subs.data(tone, tone$xpp$trialNr %in% trialRng[1]:trialRng[2])
  
  if(!is.null(tone$lfp)){
    #  tone$lfp <- tone$lfp[trialRng[1]:trialRng[2],,]
    if(length(dim(tone$lfp))==2)
      dim(tone$lfp) <- c(dim(tone$lfp), 1)
  }
  
  
  brks            <- exp( log(2) * seq(log(0.2)/log(2),by=1,to=log(12.800)/log(2)) )
  brks[1]         <- 0.190
  brks[length(brks)] <- 21
  tone$xpp$isiOct <- cut( tone$xpp$ISI, breaks=brks )
  table(tone$xpp$isiOct)
  
  brks                <- exp( log(2) * seq(log(0.2)/log(2),by=1/2,to=log(12.800)/log(2)) )
  brks[1]             <- 0.190
  brks[length(brks)] <- 21
  tone$xpp$isiOctFine <- cut( tone$xpp$ISI, breaks=brks )
  table(tone$xpp$isiOctFine)
 
  tone$shp <- shp
  
  if (show){
    showAOC(tone)
  }
  return(tone) 
}


showAOC <- function(tone,showCh=1){  
  ## -----------------------------------------------
  
  if( !is.null(tone$lfp) ){
    
    par(mfrow=c(2,2))
    lwd <- 2    
    
    if(is.null(tone$xpp$powEOG))
      tone$xpp$powEOG <- rep(0,length(dim(tone$xpp)[1]))
    
    valInd <- which( tone$xpp$p2p<450 & tone$xpp$toneType==1  )
    if(length(valInd)==0)
      valInd <- which( tone$xpp$p2p<450 & tone$xpp$toneType==2  )
    
    ylim <- c(-15,20)
    xlim <- range(tone$taxis)
    
    if(tone$what=='mua'){
      ylim <- c(-.1, 1)
    }
    
    disc   <- mean.plot( tone$lfp[valInd,,showCh], tone$taxis, add=F, sd=T, traces=F, posFUN=na.mean, 
                         widthFUN=na.se, range=NULL, col=c('black'), lwd=lwd, ylim=ylim, xlim=xlim)
    
    abline(h=0);abline(v=0)
    if(length(unique(tone$xpp$toneType))>1){
      valInd <- which( tone$xpp$p2p<450 & tone$xpp$toneType==2  )
      disc   <- mean.plot( tone$lfp[valInd,,showCh], tone$taxis, add=T, sd=T, traces=F, posFUN=na.mean, 
                           widthFUN=na.se, range=NULL, col=c('red'), lwd=lwd)
    }
    
   
    
    plot( sort(tone$xpp$p2p),1:dim(tone$xpp)[1]/dim(tone$xpp)[1], type='l', ylab='',xlab='peak to peak', xlim=c(0,1000) )
    plot( sort(tone$xpp$powEOG),1:dim(tone$xpp)[1]/dim(tone$xpp)[1], type='l', ylab='',xlab='powEOG', xlim=c(0,10000) )
    
    par(mfrow=c(1,2))
    N1  <- tapply(tone$xpp$mmn, tone$xpp$pitchOct, mean)
    xax  <- tapply(tone$xpp$pitchFreq, tone$xpp$pitchOct, mean)
    plot(xax,N1, type='b',log='x',pch=20,xlab='Pitch [Oct]')
    
    N1   <- tapply(tone$xpp$mmn, tone$xpp$isiOct, mean)
    xax  <- tapply(tone$xpp$ISI, tone$xpp$isiOct, mean)
    plot(xax,N1, type='b',log='x',xlab='ISI',pch=20)
    
    ## new plots
    if ( length(unique(tone$xpp$trialType))>1 & length(unique(tone$xpp$toneType))>1 ){
      randomVsPredictable(tone)
    }
  }
  return(tone)
}


loadAOCall <- function(tdf, what='llp2', chStr=chStr, show=F,showCh=1,rda=rda,extn='',xppextn='',trda=NA, byBlock=T){
  
  print(tdf$session[1])
  tone    <- loadAOC( tdf[1,],chStr=chStr, show=show, showCh=showCh, what=what,extn=extn,xppextn=xppextn,trda=trda,byBlock=byBlock )
  toneexp <- tone$exp
  
  if(dim(tdf)[1]>1){
    for (i in 2:dim(tdf)[1]){
      print(tdf$session[i])
      ttn   <- loadAOC( tdf[i,], chStr=chStr, show=show,showCh=showCh, what=what,extn=extn,xppextn=xppextn,trda=trda, byBlock=byBlock )
      tone <- append.data(tone,ttn)
    }
  }
  tone$animal  <- tdf[1,]$animal
  tone$xpp$exp <- as.factor(tone$xpp$exp)
  tone$exp     <- toneexp# paste('RSflex_',what,sep='')
  tone$chStr   <- chStr
  tone$what    <- what
  
  return(tone)
}

