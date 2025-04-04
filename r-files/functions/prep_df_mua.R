
## load df and chdf
colClasses <- c('character',rep('numeric',2),'character',rep('numeric',5)  ,rep('character',3),rep('numeric',6) )
df         <- read.table(paste(baseDir,'ampOddClick/ampOddClickVProbe_link.txt',sep=''),header=T,colClasses=colClasses)
df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))
df$exp     <- paste(df$session,df$probe,sep='_')


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
