prepDFVprobe <- function(df, getTrTy=F){
  
  # if probe information is not provided include all probes here
  if (!exists('probe',df)){
    numProbe      <- function(n){return(nchar(as.character(df$ElectrodeConfig[n]))-1)}
    df$numProbe   <- sapply(1:dim(df)[1],numProbe)
    sessionNr     <- rep(1:dim(df)[1],df$numProbe)
    tmp           <- df[rep(1:dim(df)[1],df$numProbe),]
    tmp$probe     <- c(1)
    sameSession   <- 1-diff( c(0,sessionNr)) 
    tmp$probe     <- 1 + sameSession
    df            <- tmp
  }
  
  nChanVprobe      <- c(24,32,32,32,32); names(nChanVprobe)<- c('v','t','s','d','z')
  extractProbeStr  <- function(n){dfline <- df[n,];substr(dfline$ElectrodeConfig,dfline$probe+1,dfline$probe+1)}
  df$probeString   <- sapply(1:dim(df)[1],extractProbeStr)
  df$numChan       <- nChanVprobe[df$probeString]
  df$startVIndex   <- 34
  extractProbeStrM <- function(n,m=2){dfline <- df[n,];substr(dfline$ElectrodeConfig,m,m)}
  
  if (max(df$probe)>1){
    for(rx in 2:max(df$probe) ){
      ndx                 <- which(df$probe == rx)
      df$startVIndex[ndx] <- df$startVIndex[ndx] + nChanVprobe[ extractProbeStrM(ndx, m=rx) ]  
    }
  }
  
  for (i in 1:length(df$session) ){
    dxt <- ifelse(df$drug[i]=='N','', paste('_',df$drug[i],sep=''))
    pxt <- ifelse(df$probe[i]==1 ,'', paste('_',df$probe[i],sep=''))
    df$exp[i] <- paste( df$session[i], dxt, pxt, sep=''  )
  }

  
  if (getTrTy){
    # try to read xpp file and determine the number of random/regular blocks
    getData <- function(i){ loadSignalProcessor(  df[i,], project=df$project[i], chStr=NA )} 
    getTable <- function(i){   
      tblTrTy <- c(NA,NA)
      xpp     <- NULL
      try( xpp <- getData(i))
      if (!is.null(xpp)){
        tblTrTy[1] <- length(which(xpp$xpp$trialType==1))
        tblTrTy[2] <- length(which(xpp$xpp$trialType==2))
      }
      return(tblTrTy) 
    }
    tblTrTy <- sapply(1:dim(df)[1], getTable)
    df$Nrand <- tblTrTy[1,]
    df$Nreg  <- tblTrTy[2,]
    
    plot(df$Nreg, col='red')
    points(df$Nrand, col='blue')
    
  }
  
  
  
  df$spacing[df$probeString=='d'] <- 150
  
  return(df)
}