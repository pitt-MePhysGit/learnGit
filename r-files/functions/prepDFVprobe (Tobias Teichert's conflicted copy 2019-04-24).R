prepDFVprobe <- function(df){
  
  nChanVprobe      <- c(24,32,32); names(nChanVprobe)<- c('v','t','s')
  extractProbeStr  <- function(n){dfline <- df[n,];substr(dfline$ElectrodeConfig,dfline$probe+1,dfline$probe+1)}
  df$probeString   <- sapply(1:dim(df)[1],extractProbeStr)
  df$numChan       <- nChanVprobe[df$probeString]
  df$startVIndex   <- 33
  extractProbeStrM <- function(n,m=2){dfline <- df[n,];substr(dfline$ElectrodeConfig,m,m)}
  
  for(rx in 2:max(df$probe) ){
    ndx                 <- which(df$probe == rx)
    df$startVIndex[ndx] <- df$startVIndex[ndx] + nChanVprobe[ extractProbeStrM(ndx, m=rx) ]  
  }
  
  
  for (i in 1:length(df$session) ){
    dxt <- ifelse(df$drug[i]=='N','', paste('_',df$drug[i],sep=''))
    pxt <- ifelse(df$probe[i]==1 ,'', paste('_',df$probe[i],sep=''))
    df$exp[i] <- paste( df$session[i], dxt, pxt, sep=''  )
  }
  
  
  return(df)
}