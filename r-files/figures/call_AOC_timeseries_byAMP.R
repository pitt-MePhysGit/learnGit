call_AOCtimeseries_byAMP <- function(what='mua',who='all', dpthStr='MUA',ganrm=T,dxdn=NA,...){
  ## ============================================ by SOA
  
  SOAndx <- 1
  AMPndx <- 1:5
  dBVal  <- seq(62, by=6, length=5)
  Ncnd   <- length(SOAndx) * length(AMPndx)
  
  
  if(!is.na(dxdn))
    tmpSub <- loadAllSessions( useInd, extn='mean_amp',who='all',what=what,xnd=dxdn, ynd=NA, ganrm=ganrm)
  
  
  rm(lfp)
  for (ax in AMPndx){
    tmp <- loadAllSessions( useInd, extn='mean_amp',who='all',what=what,xnd=ax, ynd=NA, ganrm=ganrm)
    if(!is.na(dxdn))
      tmp <- tmp - tmpSub
    
    if (ax == 1)   
      lfp <- tmp
    if(ax > 1)
      lfp <- cbind(lfp,tmp)
  }
  
  # ================================================= select valid channels and create xpp file
  chdf          <- expand2chdf(tdf)
  if(what=='csd')
    chdf          <- chdf[chdf$useCSD,]
  
  Nchan  <- dim(chdf)[1]   #dim(tst)[2]
  
  
  
  layers     <- data.frame( cmpStrName='L2',      dmin=  200, dmax= 1000, stringsAsFactors=FALSE) 
  layers[2,] <- data.frame( cmpStrName='L3',      dmin= -700,  dmax=200, stringsAsFactors=FALSE)
  layers[3,] <- data.frame( cmpStrName='L4',      dmin= -1100,  dmax= -700, stringsAsFactors=FALSE)
  layers[4,] <- data.frame( cmpStrName='L5',      dmin= -1500,  dmax= -1100, stringsAsFactors=FALSE)
  layers[5,] <- data.frame( cmpStrName='L6',      dmin= -2000,  dmax= -1500, stringsAsFactors=FALSE)
  layers[6,] <- data.frame( cmpStrName='SGSo',    dmin=   40,  dmax= 500, stringsAsFactors=FALSE)
  layers[7,] <- data.frame( cmpStrName='SGSi',    dmin= -500,  dmax=  -40, stringsAsFactors=FALSE)
  layers[8,] <- data.frame( cmpStrName='GSi' ,    dmin= -1100,  dmax= -700, stringsAsFactors=FALSE)
  layers[9,] <- data.frame( cmpStrName='MUA' ,    dmin= -1500,  dmax= -300, stringsAsFactors=FALSE)
  layers[10,] <- data.frame( cmpStrName='SG' ,     dmin= -700,  dmax= 1000, stringsAsFactors=FALSE)
  layers[11,] <- data.frame( cmpStrName='IG' ,    dmin= -2000,  dmax= -1100, stringsAsFactors=FALSE)
  layers[12,] <- data.frame( cmpStrName='ThAf' ,  dmin= -2500,  dmax= -300, stringsAsFactors=FALSE)
  layers[13,] <- data.frame( cmpStrName='All' ,  dmin= -Inf,  dmax= +Inf, stringsAsFactors=FALSE)
  
  lnd     <- which( layers$cmpStrName %in% dpthStr )
  if( length(lnd)==1 )
    dpthRng <- c( unlist(layers[lnd, c('dmin','dmax')] ))
  
  if(length(lnd)==0)
    dpthRng <- c(-Inf, Inf)
  
  chdf$use <- chdf$depth >= dpthRng[1] & chdf$depth <= dpthRng[2]
  
  tmpxpp     <- expand.grid(isiOct=SOAndx, ampInd=AMPndx)
  tmpxpp$ISI <- tmpxpp$isiOct
  tmpxpp$dB  <- dBVal[tmpxpp$ampInd]

  xpp     <- tmpxpp[rep(1:dim(tmpxpp)[1], each=Nchan), ]
  xpp$use <- rep(chdf$use,Ncnd)
  
  # ===================================== create ephys data frame
  tone       <- list()
  tone$xpp   <- xpp[xpp$use,]
  tone$lfp   <- t(lfp[,xpp$use]); dim(tone$lfp) <- c(dim(tone$lfp),1)
  tone$taxis <- seq(-150, by=1, length=900)
  
  #mean.plot(tone$lfp, tone$taxis,add=F)
  fileName <- paste('AOC_timeseries_byAMP_',what,'_',who,'_',dpthStr,'.eps',sep='') 
  if (!is.na(dxdn))
    fileName <- paste('AOC_timeseries_byAMP_Delta_',dxdn,'_',what,'_',who,'_',dpthStr,'.eps',sep='') 
  
  eps.file(fileName, pdf=T, s=2, ratio=3/2)
  disc <- figure_AOC_timeseries_byAMP(tone, whoStr='lfp', ...)
  dev.off()
  
}
