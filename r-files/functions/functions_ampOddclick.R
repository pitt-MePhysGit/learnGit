loadAmpOddClick <- function( dfline, chStr='frontoCentral',showCh=1,show=F,exp=NA, from=NA, to=NA, what='llp2', 
                       minFrq=800,dFrq=0.21783,trda=NA,extn='',xppextn=''){
  
  if(!is.na(trda))
    rda <- trda
  
  print('loadAmpOddClick')
  print(rda())
  
  if(!is.na(dfline)){
    exp    <- dfline$session
    from   <- ifelse(is.null(dfline$from),NA, dfline$from)
    to     <- ifelse(is.null(dfline$to),NA, dfline$to)
    #minFrq <- dfline$minOct
    #dFrq  <- dfline$stepOct
  }
  
  tone       <- list()
  
  tone$xpp   <- read.table( paste(rda(),exp,xppextn,'.txt',sep=''), header=T ) 
  if(nchar(extn)>0)
    print(extn)
  
  raw   <- identical(what,'raw')
  raw5k <- identical(what,'raw5k')
  llp   <- identical(what,'llp2')
  mlp   <- identical(what,'mlp')
  mbp   <- identical(what,'mbp')
  bsp   <- identical(what,'bsp')
  
  
  chextn <- ''
  if (llp)
    chextn <- extn
  if (raw)
    chextn <- paste('_raw',extn,sep='')
  if (raw5k)
    chextn <- paste('_raw5k',extn,sep='')
  if (mlp)
    chextn <- paste('_mlp',extn,sep='')
  if (mbp)
    chextn <- paste('_mbp',extn,sep='')
  if (bsp)
    chextn <- paste('_bsp',extn,sep='')
  
  tone$taxis <- NULL
  tone$lfp   <- NULL
  if (!is.na(chStr) ){# get lfps only if chStr is specified
    
    print(paste(rda(),exp,'_',chStr[1],chextn,'.mat',sep=''))
    tmpLFP   <- t( readMat( paste(rda(),exp,'_',chStr[1],chextn,'.mat',sep=''))[[1]] )
    tone$lfp <- array(NA, c(dim(tmpLFP)[1],dim(tmpLFP)[2], length(chStr)) )
    tone$lfp[,,1] <- tmpLFP 
    if(length(chStr)>1){
      for (chnd in 2:length(chStr)){
        tone$lfp[,,chnd]   <- t( readMat( paste(rda(),exp,'_',chStr[chnd],chextn,'.mat',sep=''))[[1]] )
        print(paste(rda(),exp,'_',chStr[chnd],chextn,'.mat',sep=''))
      }
    }
    
    textn <- ''
    if (mlp)
      textn <- '_mlp'
    if (mbp)
      textn <- '_mbp' 
    if (bsp)
      textn <- '_bsp'
    if (raw5k)
      textn <- '_raw5k'
    
    print(paste(rda(),exp,'_taxis',textn,'.mat',sep=''))
    #if(!raw5k)
    tone$taxis <- readMat( paste(rda(),exp,'_taxis',textn,'.mat',sep=''))[[1]]
    
    if(raw5k)# new correction for raw5k time-axis
      tone$taxis   <- seq(-100,by=0.2, length=dim(tone$lfp)[2] )
    
    # non-baseline-corrected lfp
    lfpnbl <- tone$lfp
    
    tmp            <- baseLineCorrection(tone$lfp,c(-50,0), tone$taxis)
    tone$lfp       <- tmp$res
    tone$xpp$blcor <- tmp$correction
  }
  
  tone$xpp$exp <- rep(exp,length(tone$xpp$seqNr))
  
  ## write out some info about the data contained =========================
  tone$exp     <- exp
  tone$chStr   <- chStr
  tone$what    <- what
  
  ## prepare xpp =====================================
  # already implemented 
  #newTrInd                    <- which(diff(c(0,tone$xpp$trialNr))==1)
  #tone$xpp$toneType[newTrInd] <- 1
  
  #amplitude              <- seq(from=56,by=8,length=5) 
  amplitude              <- seq(from=62,by=6,length=5) 
  tone$xpp$dB            <- amplitude[tone$xpp$ampInd]
  
  pastAmpInd             <- my.lag(tone$xpp$ampInd,1)
  pastAmpInd[1]          <- NA
  tone$xpp$deltaAmpInd   <- tone$xpp$ampInd - pastAmpInd 
  
  pastISI               <- my.lag(tone$xpp$ISI,1)
  pastISI[c(1,2)]       <- NA
  tone$xpp$deltaISI     <- tone$xpp$ISI - pastISI 
  tone$xpp$deltaISIfrac <- (pastISI+tone$xpp$deltaISI)/(pastISI)
  
  
  rpCnt             <- as.numeric(tone$xpp$toneType)==1
  tmp               <- get.epoch( as.numeric(tone$xpp$toneType)==1  )
  tmp$dur           <- tmp$offset - tmp$onset
  rpCnt[tmp$offset] <- -tmp$dur  
  tone$xpp$repNr    <- cumsum(rpCnt)
  
  tone$xpp$ampDev  <- abs(tone$xpp$deltaAmpInd )     >  0       
  tone$xpp$timeDev <- abs(tone$xpp$deltaISIfrac-1)   >= 0.05  # time Deviant if ISI differs by more than 5% from previous ISI
  
  ## subset range of valid tones if specified ======================
  trialRng <- c(1,dim(tone$xpp)[1])
  if(!is.na(from))
    trialRng[1] <- from
  
  if(!is.na(to))
    trialRng[2] <- to
  
  tone$xpp <- tone$xpp[trialRng[1]:trialRng[2],]
  
  if(!is.null(tone$lfp)){
    tone$lfp <- tone$lfp[trialRng[1]:trialRng[2],,]
    if(length(dim(tone$lfp))==2)
      dim(tone$lfp) <- c(dim(tone$lfp), 1)
  }
  
  
  brks            <- exp( log(2) * seq(log(0.2)/log(2),by=1,to=log(25.600)/log(2)) )
  brks[1]         <- 0.190
  tone$xpp$isiOct <- cut( tone$xpp$ISI, breaks=brks )
  table(tone$xpp$isiOct)
  
  brks                <- exp( log(2) * seq(log(0.2)/log(2),by=1/2,to=log(25.600)/log(2)) )
  brks[1]             <- 0.190
  tone$xpp$isiOctFine <- cut( tone$xpp$ISI, breaks=brks )
  table(tone$xpp$isiOctFine)
  
  ## create difference components
  tone$xpp$dn1p1 <- tone$xpp$n1 - tone$xpp$p1
  tone$xpp$dn1px <- tone$xpp$n1 - tone$xpp$px
  
  ## create normalized component amplitude
  cmpStrList <- c('p1x.r','p1a.r','p1b.r','nx.r','px.r','n1.r','p2.r','n2.r')
  cmpStrListZ <- paste('z.',cmpStrList,sep='')
  tmpdf       <- tone$xpp[,cmpStrList] 
  names(tmpdf)<- cmpStrListZ
  for (j in 1:dim(tmpdf)[2] )
    tmpdf[,cmpStrListZ[j]] <- robustZscore( tone$xpp[,cmpStrList[j]],subs=tone$xpp$p2p<450&(!is.na(tone$xpp$isiOct))  )
  
  tone$xpp <- data.frame(tone$xpp,tmpdf)
  tone$xpp$mmn       <- tone$xpp$n1.r
  
  print( apply(tone$xpp[,cmpStrListZ],2,median) )
  ## -----------------------------------------------
  if(show & !is.null(tone$lfp) ){
    par(mfrow=c(2,2))
    lwd <- 2    
    if(is.null(tone$xpp$powEOG))
      tone$xpp$powEOG <- rep(0,length(dim(tone$xpp)[1]))
      
    valInd <- which( tone$xpp$p2p<450 & tone$xpp$toneType==1  )
    if(length(valInd)==0)
      valInd <- which( tone$xpp$p2p<450 & tone$xpp$toneType==2  )
    
    ylim <- c(-15,20)
    xlim <- range(tone$taxis)
    
    if(mlp){
      ylim <- c(-5, 5)
      xlim <- c(min(tone$taxis), 55 )
    }
    
    disc   <- mean.plot( tone$lfp[valInd,,showCh], tone$taxis, add=F, sd=T, traces=F, posFUN=na.mean, 
                         widthFUN=na.se, range=NULL, col=c('black'), lwd=lwd, ylim=ylim, xlim=xlim)
    
    abline(h=0);abline(v=0)
    if(length(unique(tone$xpp$toneType))>1){
      valInd <- which( tone$xpp$p2p<450 & tone$xpp$toneType==2  )
      disc   <- mean.plot( tone$lfp[valInd,,showCh], tone$taxis, add=T, sd=T, traces=F, posFUN=na.mean, 
                           widthFUN=na.se, range=NULL, col=c('red'), lwd=lwd)
    }
    
    normData <- lfpnbl[,,showCh] #$baseLineCorrection(tone$lfp[,,showCh],c(-151,301), tone$taxis)
    ##normData <- baseLineCorrection(tone$lfp[,,showCh],c(-151,301), tone$taxis)$res
    
    valInd <- which( tone$xpp$p2p<450 & tone$xpp$toneType==1  )
    disc   <- mean.plot( normData[valInd,],    tone$taxis, add=F, sd=F, traces=F, posFUN=na.sd,   widthFUN=na.se, range=NULL, col='black',lwd=lwd, ylim=c(0,50))
    if(length(which(tone$xpp$toneType==2))>1){
      valInd <- which( tone$xpp$p2p<450 & tone$xpp$toneType==2  )
      disc   <- mean.plot( normData[valInd,],    tone$taxis, add=T, sd=F, traces=F, posFUN=na.sd,   widthFUN=na.se, range=NULL, col='red',lwd=lwd)  
    }
    abline(v=0)
    
    plot( sort(tone$xpp$p2p),1:dim(tone$xpp)[1]/dim(tone$xpp)[1], type='l', ylab='',xlab='peak to peak', xlim=c(0,1000) )
    plot( sort(tone$xpp$powEOG),1:dim(tone$xpp)[1]/dim(tone$xpp)[1], type='l', ylab='',xlab='powEOG', xlim=c(0,10000) )
    
    par(mfrow=c(1,2))
    N1  <- tapply(tone$xpp$mmn, tone$xpp$ampInd, mean)
    xax  <- tapply(tone$xpp$dB, tone$xpp$ampInd, mean)
    plot(xax,N1, type='b',log='x',pch=20,xlab='Amplitude [dB]')
    
    N1   <- tapply(tone$xpp$mmn, tone$xpp$isiOct, mean)
    xax  <- tapply(tone$xpp$ISI, tone$xpp$isiOct, mean)
    plot(xax,N1, type='b',log='x',xlab='ISI',pch=20)
    
    ## new plots
    #if ( length(unique(tone$xpp$trialType))>1 & length(unique(tone$xpp$toneType))>1 ){
    #  randomVsPredictable(tone)
    #}
  }
  return(tone)
}





loadAmpOddClickAll <- function(tdf, what='llp2', chStr=chStr, show=F,showCh=1,rda=rda,extn='',xppextn='',trda=NA){
  
  print(tdf$session[1])
  tone    <- loadAmpOddClick( tdf[1,],chStr=chStr, show=show, showCh=showCh, what=what,extn=extn,xppextn=xppextn,trda=trda )
  toneexp <- tone$exp
  
  if(dim(tdf)[1]>1){
    for (i in 2:dim(tdf)[1]){
      print(tdf$session[i])
      ttn   <- loadAmpOddClick( tdf[i,], chStr=chStr, show=show,showCh=showCh, what=what,extn=extn,xppextn=xppextn,trda=trda )
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



showGlobalPower <- function(tn, extn='vanilla', eps=T, tlim=c(-25,250), show=T,cmpStrList=NULL,electrodeIndex=NULL,
                            showCompPoly=F,ylim=NULL,highlight=c(4,6),highCol=c('green','blue'),... ){
  
  dotargs <- list(...)
  
  animal <- strsplit(tn$exp,'_')[[1]][1]
  if (animal=='Rockey')
    eCogFile <- 'rockey'
  
  if (animal=='Jesse')
    eCogFile <- 'jessesym'
  
  if (animal=='Walter')
    eCogFile <- 'waltersym'
  
  if (animal=='Sam')
    eCogFile <- 'samsym'
  
  ep             <- read.table(paste(baseDir,'/eCog/',eCogFile,'/positionList_ISI_eCog.txt',sep=''),header=T)  
  
  if(is.null(electrodeIndex))
    electrodeIndex <- ep$Number[which(ep$use==1)] ## tmp[sind]
  
  valInd <- 1:dim(tn$lfp)[1]
  if (exists('xpp',where=tn))
    valInd <- which(tn$xpp$p2p<450 & tn$xpp$ISI>0 & tn$xpp$ISI<15.5 )
  
  sdat   <- tn$lfp[valInd,,]
  
  mnMat <- apply(sdat,c(2,3),na.mean)  
  powerMat <- apply(mnMat[,electrodeIndex], 1, na.sd)

  if(show){
    ## -------------------------------------------
    if(!eps)
      par(mfrow=c(1,1))
    
    if(eps)
      eps.file( paste('showModel_erp_',extn,'.eps',sep=''),ratio=0.67 )
    
    if(is.null(ylim))
      ylim <- c(-1,1)*na.max(c(mnMat[,electrodeIndex]))
   
    tfn <- function(x){squeeze(x, to.range=c(1,2),from.range=ylim  )}
    #if(is.null(cmpStrList))
    plot(tn$taxis, tfn(mnMat[,electrodeIndex[1]]) ,type='l', ylim=c(0,2), xlim=tlim,..., xaxt='n',yaxt='n', bty='n')
    points( c(0,0), tfn(ylim), type='l',lwd=2 )
    axis(1) 
    axis(2, at=tfn(ylim), label=ylim)
    abline( h=tfn(0),lwd=2 )
    
    if(!is.null(cmpStrList)){
      compSuperList <- ploygonERPComponents( cmpStrList, show=showCompPoly, add=T,ylim=tfn(ylim),tlim=tlim, animal=animal )
      ##points( tn$taxis,mnMat[,1],type='l',...)
    }
    
    for(i in electrodeIndex )
      points(tn$taxis,tfn(mnMat[,i]),type='l',...)
    
    if(!is.null(highlight) ){
      tlwd <- 3
      if(!is.null(dotargs$lwd) )
        tlwd <- 2* dotargs$lwd
        
      for(i in 1:length(highlight)){
        points(tn$taxis, tfn(mnMat[,highlight[i]]), col=rep(highCol,length(highlight) )[i],type='l',lwd=tlwd ) 
      }
    }
      
    
    if (!is.null(cmpStrList) & !showCompPoly){
  
      addCmpLines <- function(n){
        tind  <- which(resFrame$taxis >= compSuperList[[n]]$trng[1] & resFrame$taxis <= compSuperList[[n]]$trng[2])
        chInd <- compSuperList[[n]]$chInd
        
        for(i in chInd  )
          points(tn$taxis[tind],tfn(mnMat[tind,i]),type='l',col=(n+1), lwd=2)
      }
      disc <- lapply(1:length(compSuperList), addCmpLines )
    }
#     if(eps)
#       dev.off()
#     
#     ## -------------------------------------------
#     if(eps)
#       eps.file( paste('showModel_globalPower_',extn,'.eps',sep=''),ratio=0.67 )  
#     
    ## second part of the plot, the actual global power
    ylim <- c(0,ceiling(max(c(powerMat))))
    
    tfn <- function(x){squeeze(x, to.range=c(0,.8),from.range=ylim  )}
    ##if(is.null(cmpStrList) )
    points(tn$taxis,tfn(powerMat),type='l',xlim=tlim,ylim=ylim,...)  
    points( c(0,0), tfn(c(0,ylim[2])), type='l',lwd=2 )
    axis(2, at=tfn(ylim), label=ylim)


    if(!is.null(cmpStrList)){
      compSuperList <- ploygonERPComponents( cmpStrList, show=T, add=T,tlim=tlim, ylim=tfn(ylim), animal=animal )
      ##points( tn$taxis,powerMat,type='l',...)  
    }
    abline( h=tfn(0),lwd=2 )
    
    if(eps)
      dev.off()
  }
  
  return(list(powerMat=powerMat, mnMat=mnMat))
  
}



showModelGlobalPower <- function(parm, tn, extn='vanilla', eps=T, tlim=c(-25,250) ){
  
  animal <- strsplit(tn$exp,'_')[[1]][1]
  if (animal=='Rockey')
    eCogFile <- 'rockey'
  
  if (animal=='Jesse')
    eCogFile <- 'jessesym'
  
  ep             <- read.table(paste(baseDir,'/eCog/',eCogFile,'/positionList_ISI_eCog.txt',sep=''),header=T)  
  electrodeIndex <- ep$Number[which(ep$use==1)] ## tmp[sind]
  
  x <- tn$xpp
  
  if ( !is.na(parm[1]) ){
    lmd <- parm[1]
    sgm <- parm[2]
    
    memnl  <- FALSE     # memory non-linearity. If true, pass tpKer through non-linearity
    if(length(parm)>2){
      memnl <- TRUE
      mem <- parm[3]
    }
    
    ## prepare the history kernels
    Nhist      <- 50
    
    ## define the pitch-tuning function
    tmpfn <- function(zz){ dnorm(zz,m=0,sd=1/abs(sgm)) / dnorm(0,m=0,sd=1/abs(sgm)) }
    tdn <- function(z){ 
      
      if(is.array(z) )
        res <- array(NA, dim(z) )
      if(!is.array(z) )
        res <- rep(NA,length(z))
      
      res[ is.finite(z) ] <- 1
      
      if(sgm>0)
        res <- tmpfn(z)
      
      if(sgm<0)
        res <- 1+tmpfn(4)-tmpfn( z ) 
      
      return(res)
    }
    histRes    <- getHist(tn,Nhist=Nhist,useOct=TRUE)
    timeHist   <- histRes$timeHist
    toneHist   <- histRes$toneHist
    
    x$tKer  <- apply( exp(  lmd*timeHist),1,na.sum )
    x$pKer  <- apply( tdn(toneHist), 1,na.sum )
    x$tpKer <- apply( exp( lmd*timeHist) * tdn(toneHist), 1,na.sum )
  
    
    ybin <- ( x$p2p<450 & x$ISI>0 & x$ISI<15.5 & is.finite(x$tpKer) )
    yind <- which(ybin)
  
    ## subset the data set
    y           <- subset(x, ybin)
    y$tpKerNorm <- y$tpKer - mean(y$tpKer)
  
    ## create the binned Kernel 
    Nbin <- 6
    binX<- quantile(y$tpKer, prob=seq(1/Nbin,by=1/Nbin,length=Nbin-1) )
  
    binKer <- rep(1, dim(y)[1] )
    for (bnd in 1:length(binX)){
      binKer <- binKer + (y$tpKer>binX[bnd])
    }
    binKer <- binX[binKer]
  }
  
  ## create groups based on ISI and pitch 
  #if ( is.na(parm[1]) ){
  #  #browser()
    
  #  ybin <- ( x$p2p<450 & x$ISI>0 & x$ISI<15.5 )
  #  yind <- which( ybin )
    
    ## subset the data set
  #  y           <- subset(x, ybin)
    
    brks          <- exp( seq(-1.7,by=.25,to=2.75) )
    binKerISI        <- cut( y$ISI, breaks=brks )
    binKerdP        <- y$deltaPitch
  #}
  
  sdat <- tn$lfp[yind,,electrodeIndex]
  
  mnLstC  <- lapply( sort(unique(binKer)) , function(x){apply( sdat[which(binKer==x),,],c(2,3),mean) } )
  mnMatC <- array(unlist(mnLstC), c(dim(mnLstC[[1]]),length(mnLstC)) );
  
  mnLstCpitch  <- lapply( sort(unique(binKerdP)) , function(x){apply( sdat[which(binKerdP==x),,],c(2,3),mean) } )
  mnMatCpitch <- array(unlist(mnLstCpitch), c(dim(mnLstCpitch[[1]]),length(mnLstCpitch)) );
  
  mnLstCISI  <- lapply( sort(unique(binKerISI)) , function(x){apply( sdat[which(binKerISI==x),,],c(2,3),mean) } )
  mnMatCISI <- array(unlist(mnLstCISI), c(dim(mnLstCISI[[1]]),length(mnLstCISI)) );
  
  
  mnMat <- apply(sdat,c(2,3),mean) 
  
  adaptMatC      <- apply( apply(mnMatC,     c(1,2),sd), 1, mean)
  adaptMatCISI   <- apply( apply(mnMatCISI,  c(1,2),sd), 1, mean)
  adaptMatCpitch <- apply( apply(mnMatCpitch,c(1,2),sd), 1, mean)
  
  adaptMat  <- apply( apply(sdat,c(2,3),sd), 1, mean)
  
  
  
  powerMat <- apply(mnMat, 1, sd)
  
  if(!eps)
    par(mfrow=c(3,1))
  
  
  ## -------------------------------------------
  if(eps)
    eps.file( paste('showModel_erp_',extn,'.eps',sep=''),ratio=0.67 )
  plot(tn$taxis,mnMat[,1],type='l', ylim=range(mnMat),xlim=tlim)
  for(i in 1:dim(mnMat)[2])
    points(tn$taxis,mnMat[,i],type='l')
  if(eps)
    dev.off()
  
  ## -------------------------------------------
  if(eps)
    eps.file( paste('showModel_globalPower_',extn,'.eps',sep=''),ratio=0.67 )  
  plot(tn$taxis,powerMat,type='l',xlim=tlim,lwd=2)  
  if(eps)
    dev.off()
  
  ## -------------------------------------------
  if(eps)
    eps.file( paste('showModel_adaptation_',extn,'.eps',sep=''),ratio=0.67 )
  
  plot(tn$taxis,adaptMatC,type='l',xlim=tlim,lwd=2, ylim=c(0,max(adaptMatC,adaptMatCpitch,adaptMatCISI))) 
  lines(tn$taxis,adaptMatCpitch,lwd=2,col='blue')
  lines(tn$taxis,adaptMatCISI,lwd=2,col='red')
  
  if(eps)
    dev.off()
  
  
  return(list(powerMat=powerMat, adaptMat=adaptMatC))
}





showISIbyKet <- function(animal,dataSet='l',speaker='TB',what='llp2',eps=F,pdf=T,baseName='ketamineByISI',
                         doseListCut=c(0.01,0.2),chStr=chStr,df=tdf){
  
  doseListCut <- c(-.01, doseListCut, Inf)
  table(df$animal,df$ket,df$dataSet)
  indList <- list()
  
  for (j in 1:(length(doseListCut)-1) ){
    indList[[j]] <- which(df$animal==animal & df$dataSet==dataSet & df$speaker==speaker & df$ket>doseListCut[j] &  df$ket<doseListCut[j+1] )
    print( df[indList[[j]],] )
  }
  doseList <- unlist( lapply(indList, function(nd){ mean(df$ket[nd]) }))
  cR <- colorRampPalette(c('black','blue','green','orange','red'), bias = 1, space ="Lab",interpolate="spline")
  
  xlim <- c(-10,150)
  ylim <- c(-75,50)
  
  if(what=='mlp'){
    xlim <- c(-10,70)
    ylim <- c(-15,15)    
  }
  
  if(what=='raw5k'){
    xlim <- c(-10,150)
    ylim <- c(-75,50)    
  }
  
  lwd <- 1
  
  if(eps)
    eps.file( paste(animal,'_',dataSet,'/',what,'/', baseName, what,'.eps',sep=''),mf=c(1,3),ratio=1, s=2 )
  
  par(mfcol=c(2               ,1),mar=c(3,2,2,1)+.1)
  
  datList <-  list()
  for (dind in 1:length(doseList) ){
    print( df[indList[[dind]],] )
    dat             <- loadAmpOddClickAll( df[indList[[dind]],] , chStr=chStr, what=what)
    
    dat$animal      <- df[ indList[[dind]][1], 'animal']
    datList[[dind]] <- dat
    datList[[dind]]$exp <- paste(animal,dataSet,sep='_')
    
    dat <- subs.data(dat, dat$xpp$trialNr>4 & dat$xpp$p2p<450 )
    
    tbl   <- table(dat$xpp$isiOct)
    print(tbl)
    
    mnISI <- tapply(dat$xpp$ISI, dat$xpp$isiOct, mean )
    TBind <- which(mnISI>0.2 & mnISI<16)
    print(tbl[TBind])
    
    if (dind == 1)
      palette(cR(length(TBind)))
    
    for(n in 1:(length(TBind)-1) ){
      
      valInd <- which( dat$xpp$isiOct == levels(dat$xpp$isiOct)[TBind[n]]  )
      
      if(length(valInd)>20)
        disc   <- mean.plot( dat$lfp[valInd,,], dat$taxis, add=ifelse(n==1,F,T), sd=TRUE, traces=F, 
                             posFUN=na.mean, widthFUN=na.se, range=NULL, col=n,lwd=lwd, ylim=ylim,
                             xlim=xlim)
    }
    abline(h=0)
  }
  if(eps)
    dev.off()
  
  return(datList)
}




runAllTests_ampOddClick <- function(cmpStr,thisall=all){
  
  ## ====================================================================
  ## Only saline days, first order analysis (combining data from all days)
  
  thisall$xpp$DV     <- thisall$xpp[[cmpStr]]
  thisall$xpp$isiLog <- thisall$xpp$isiLog - log2(mxISI) 
  
  y         <- subset(thisall$xpp, thisall$xpp$drug==FALSE )
  lmSalFull <- aov( DV ~ 1 + isiLog + ampInd + trialNr, data=y)
  resSal    <- Anova(lmSalFull, type=2)
  
  
  ## ===================================================================  
  ## Saline and MK-801 days, first order analysis (combining data from all days)
  y         <- thisall$xpp
  lmAllFull <- aov( DV ~ 1 + drug + isiLog + ampInd + trialNr + isiLog:drug + ampInd:drug , data=y)
  resAll    <- Anova(lmAllFull, type=2)
  Anova(lmAllFull,type=2)
  
  
  ## ======================================================================================================
  mnTDexp     <- t(tapply(y$DV, list(y$isiOct,y$exp), mean  ))
  mnLDexp     <- t(tapply(y$DV, list(y$ampInd,y$exp), mean  ))
  
  drg         <- tapply(y$drug, y$exp, mean)
  idataT      <- data.frame(IV=1:5)
  
  lm1DT  <- lm(  mnTDexp ~ 1 + drg )
  lm1T   <- lm(  mnTDexp ~ 1)
  lm0T   <- lm(  mnTDexp ~ 0)
  
  anova(lm1DT, M=~IV, X= ~1,       idata=idataT, test="Spherical")
  anova(lm1DT,lm1T, lm0T, M=~IV, X= ~1,       idata=idataT, test="Spherical")
  anova(lm1DT, M=~IV, X= ~1,       idata=idataT, test="Spherical")
  
  
  # non-repeated analysis
  dfTD <- data.frame( exp=dimnames(mnTDexp)[[1]], IV=rep(1:5, each=dim(mnTDexp)[1]), drug=c(drg), DV=c(mnTDexp)  )
  dfLD <- data.frame( exp=dimnames(mnLDexp)[[1]], IV=rep(1:5, each=dim(mnLDexp)[1]), drug=c(drg), DV=c(mnLDexp)  )
  dfTD$drug <- as.factor(dfTD$drug)
  dfLD$drug <- as.factor(dfLD$drug)
  
  lmAllExpTD   <- aov( DV ~ 1 + drug + IV + drug:IV + exp, data=dfTD )
  resAllExpTD  <- Anova(lmAllExpTD,type=2)
  Anova(lmAllExpTD, type=2)
  
  lmAllExpLD     <- aov( DV ~ 1 + drug + IV + drug:IV, data=dfLD )
  resAllExpLD    <- Anova(lmAllExpLD)
  Anova(lmAllExpLD)

  # extrace one slope value per day for LDAEP and TDAEP
  Nexp   <- length(unique(y$exp))
  
  lmtmp  <- aov( DV ~ 1 + exp + isiLog:exp + ampInd:exp, data=y)
  dfLT   <- data.frame(drug=drg, expNr=1:Nexp, TD=lmtmp$coefficients[Nexp+(1:Nexp)], LD=lmtmp$coefficients[2*Nexp+(1:Nexp)], LTI=lmtmp$coefficients[3*Nexp+(1:Nexp)] )
  
  # test effect of loudness and time in saline condition
  slopeTDSal <- Anova( aov( TD~1+expNr, data=subset(dfLT,dfLT$drug==0)), type=3 )
  slopeLDSal <- Anova( aov( LD~1+expNr, data=subset(dfLT,dfLT$drug==0)), type=3 )
  
  # test effects of drug
  slopeTD <- Anova( aov( TD~drug+expNr, data=dfLT),type=2 )
  slopeLD <- Anova( aov( LD~drug+expNr, data=dfLT),type=2 )
  
  ## FIGURES
  ## =================================================== TDAEP
  xax      <- tapply(log(y$ISI), list(y$isiOct,y$drug), mean  )
  mnN1     <- tapply(y$DV, list(y$isiOct,y$drug), mean  )
  seN1     <- tapply(y$DV, list(y$isiOct,y$drug), na.se  )
  N1fit    <- tapply(lmAllFull$fitted, list(y$isiOct,y$drug), mean  )
  
  valOct <- 1:5 
  nrmFctr <- N1fit[max(valOct),1]# lmExp$coefficient[1] - lmExp$coefficients[3]
  mnN1    <- mnN1[valOct,]/nrmFctr
  seN1    <- seN1[valOct,]/nrmFctr
  N1fit <- N1fit[valOct,]/nrmFctr
  xax   <- xax[valOct,]
  
  errorplot(xax, mnN1, ey=seN1, type='p', wd=.05, main=cmpStr)
  points(xax[,1], mnN1[,1], type='p', pch=20, col='black')
  points(xax[,2], mnN1[,2], type='p', pch=20, col='red')
  points(xax[,1], N1fit[,1], type='l')
  points(xax[,2], N1fit[,2], type='l',col='red')
  abline(h=c(0,1))  
  mnTDAEP <- mnN1
  
  ## =====================================================  LDAEP
  lax      <- tapply(y$dB, list(y$ampInd,y$drug), mean  )
  mnN1     <- tapply(y$DV, list(y$ampInd,y$drug), mean  )
  seN1     <- tapply(y$DV, list(y$ampInd,y$drug), na.se  )
  N1fit    <- tapply(lmAllFull$fitted, list(y$ampInd,y$drug), mean  )
  
  valOct <- 1:5 
  nrmFctr <- N1fit[max(valOct),1]# lmExp$coefficient[1] - lmExp$coefficients[3]
  mnN1    <- mnN1[valOct,]/nrmFctr
  seN1    <- seN1[valOct,]/nrmFctr
  N1fit   <- N1fit[valOct,]/nrmFctr
  lax     <- lax[valOct,]
  
  errorplot(lax, mnN1, ey=seN1, type='p', wd=1, main=cmpStr)
  points(lax[,1], mnN1[,1], type='p', pch=20, col='black')
  points(lax[,2], mnN1[,2], type='p', pch=20, col='red')
  points(lax[,1], N1fit[,1], type='l')
  points(lax[,2], N1fit[,2], type='l',col='red')
  abline(h=c(0,1))  
  mnLDAEP     <- mnN1
  
  return(list(lmSalFull=lmSalFull$coefficients, resSal=resSal, lmAllFull=lmAllFull$coefficients, resAll=resAll, mnTDAEP=mnTDAEP, mnLDAEP=mnLDAEP, 
              slopeTDSal=slopeTDSal, slopeLDSal=slopeLDSal, slopeTD=slopeTD, slopeLD=slopeLD,dfLT=dfLT))
}

