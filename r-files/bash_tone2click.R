## master FFR tone2click
rm(list=ls())
defpar <- par()
library(parallel)

args <- list('animal_Ben','date_20201007')

Nargs <- length(args)
objectList <- list()
valueList  <- list()
for (ax in 1:Nargs){
  objectList[ax] <- strsplit(args[[ax]],'_')[[1]][1]
  valueList[ax]  <- strsplit(args[[ax]],'_')[[1]][2]
}
#argDF <- data.frame(object=objectList, value=valueList)
print(objectList)
print(valueList)


baseDir <- '~/Dropbox/'
#source( paste(baseDir,'ampOdd_click/r-files/functions/function_wavelet.R',sep='') )


# read the curated file
colClasses <- c('character',rep('numeric',11),rep('character',3),rep('numeric',4),rep('character',8))
df        <- read.table(paste(baseDir,'FFR/FFR_tone2click_link.txt',sep=''),header=T,colClasses=colClasses)
df$animal  <- unlist( lapply(1:dim(dfC)[1], function(n){    animal <- strsplit(dfC$session[n],'_')[[1]][1] }  ))

# identify recordings that are separated into two or more files
getDateTime    <- function(n,m=1,tdf=df){ strsplit(tdf$session[n],'_')[[1]][m] }
df$animal      <- sapply(1:dim(df)[1],getDateTime, m=1, tdf=df)
df$date        <- sapply(1:dim(df)[1],getDateTime, m=2, tdf=df)
df$time        <- sapply(1:dim(df)[1],getDateTime, m=3, tdf=df)
df$depthStr    <- sapply(1:dim(df)[1],getDateTime, m=4, tdf=df)

getDepth       <- function(n,tdf=df){ neg <- (substr(tdf$depthStr[n],2,2)=='m'); res <- ifelse(neg, -1 * as.numeric(substr(tdf$depthStr[n],3,10)), as.numeric(substr(tdf$depthStr[n],2,10)));return(res) }
df$depth       <- sapply(1:dim(df)[1],getDepth, tdf=df)
df$day         <- paste(df$animal, df$date, sep='_')



## preliminaries
pp   <- function(){paste(baseDir, 'ampOdd_click/word/VProbe/plots/',sep='')}
rda  <- function(){'~/Dropbox/ampOdd_click/Rdata/'}
trda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}


## ------------ extract list of valid sessions
valMat <- array(0, c(Nargs,dim(df)[1]))
for (ax in 1:Nargs)
  valMat[ax,] <- df[[ objectList[[ax]] ]] == valueList[[ax]]

for (ax in 1:Nargs)
  valMat[ax,grep( valueList[[ax]], df[[ objectList[[ax]] ]])] <- 1 


valBin <- apply(valMat, 2, sum)
valInd <- which(valBin==Nargs)
print(df[valInd,])


Nclus             <- 12

run_deconvolution <- T


## ==================================================================
cluster               <- (Nclus > 1) & (length(valInd)>1)

# ====================================== only one session
if (length(valInd)==1){
  if (run_deconvolution)
    function_ffr_click_deconvolution( valInd[1] )
}

# ==================================== many sessions, no multi-core
if (length(ind)>1 & !cluster & run_deconvolution)
  disc <- lapply(valInd, function_ffr_click_deconvolution )

# =================================== many sessions multi-core
if (length(ind)>1 & cluster){
  nClus <- min( c(length(ind),  Nclus, parallel::detectCores()-2))
  parallelCluster <- parallel::makeCluster(nClus,outfile='')
  print(parallelCluster)
  
  if (run_deconvolution){
    tryCatch(
      disc <- parallel::parLapply(parallelCluster, ind, function_ffr_click_deconvolution), error = function(e) print(e)
    )
  }
  
  
  # Shutdown cluster neatly
  if(!is.null(parallelCluster)) {
    parallel::stopCluster(parallelCluster)
    parallelCluster <- c()
  }
}