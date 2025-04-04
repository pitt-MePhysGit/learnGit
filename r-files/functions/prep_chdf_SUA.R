prep_chdf_SUA <- function(chdf){
  
  extractProbeNr    <- function(n){
    dfline    <- chdf[n,]
    numProbes <- nchar(dfline$ElectrodeConfig)-1 
    probeStrings <- rep('x',numProbes)
    for (px in 1:numProbes)
      probeStrings[px] <- substr(dfline$ElectrodeConfig,1+px,1+px)
    
    numChanProbe <- c(24,32,32); names(numChanProbe) <- c('v','t','d')
    probeCutoff  <- cumsum( numChanProbe[probeStrings]  )
    probe        <- min( which(dfline$channel<=probeCutoff)  ) 
    return(probe)
  }
  chdf$probe   <- sapply(1:dim(chdf)[1],extractProbeNr)
  
  extractProbeStr    <- function(n){
    dfline <- chdf[n,]
    substr(dfline$ElectrodeConfig,dfline$probe+1,dfline$probe+1)
  }
  chdf$probeString   <- sapply(1:dim(chdf)[1],extractProbeStr)
  
  channelNrInProbe <- function(n){
    dfline    <- chdf[n,]
    numProbes <- nchar(dfline$ElectrodeConfig)-1 
    probeStrings <- rep('x',numProbes)
    for (px in 1:numProbes)
      probeStrings[px] <- substr(dfline$ElectrodeConfig,1+px,1+px)
    
    numChanProbe <- c(24,32,32); names(numChanProbe) <- c('v','t','d')
    probeCutoff  <- cumsum( numChanProbe[probeStrings]  )
    probe        <- min( which(dfline$channel<=probeCutoff)  )
    relativeChannel <- dfline$channel - c(0,probeCutoff)[probe]
    return(relativeChannel)
  }
  chdf$channelNrInProbe   <- sapply(1:dim(chdf)[1],channelNrInProbe)
  
  lyr          <- df$Layer; names(lyr) <- df$session
  chdf$Layer   <- lyr[chdf$session]
  spc          <- df$spacing; names(spc) <- df$session
  chdf$spacing <- spc[chdf$session]
  
  chdfDepth <- function(n){# a function that determines relative recording depth based on electrod config and channel contact
    dfline     <- chdf[n,]
    chDepthInd <- 1:32
    if(dfline$probeString == 'd')
      chDepthInd <- c(1:12,rep(13:20,each=2), 21:24)
    
    chDepth <-  dfline$spacing * (-chDepthInd[dfline$channelNrInProbe] + dfline$Layer) 
    return(chDepth)
  }
  chdf$depth   <- sapply(1:dim(chdf)[1],chdfDepth)
  chdf$animal  <- unlist( lapply(1:dim(chdf)[1], function(n){    animal <- strsplit(chdf$session[n],'_')[[1]][1] }  ))
  
  numChan      <- c(24,32,32); names(numChan) <- c('v','t','d')
  chdf$numChan <- numChan[chdf$probeString] 
  
  
  # hard wired for now... fix as soon as using more than one probe or changing number of eeg channels
  chdf$startVindex <- 34
  
  
  
  return(chdf)
}