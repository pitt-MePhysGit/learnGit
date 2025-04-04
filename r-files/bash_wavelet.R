## master VProbe STPSP
rm(list=ls())
defpar <- par()
library(parallel)

baseDir <- '~/Dropbox/'
source(paste(baseDir, 'ampOdd_click/r-files/functions/prepDFVprobe.R',sep=''))
source( paste(baseDir,'ampOdd_click/r-files/functions/function_wavelet.R',sep='') )
source( paste(baseDir,'ampOdd_click/r-files/functions/function_aoc_featureExtraction.R',sep='') )
source( paste(baseDir,'ampOdd_click/r-files/functions/function_aoc_featureFitting.R',sep='') )

source( paste(baseDir,'ampOdd_click/r-files/functions/function_aoc_featureDecoding.R',sep='') )

# read the curated file
colClasses <- c('character',rep('numeric',2),'character',rep('numeric',5)  ,rep('character',3),rep('numeric',9) )
df         <- read.table(paste(baseDir,'ampOddClick/ampOddclickVProbe_link.txt',sep=''),header=T,colClasses=colClasses)
df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))
df         <- prepDFVprobe(df)

existFeatureMat <- function(nx){ file.exists(paste('~/Dropbox/ampOdd_click/Rdata/',df$session[nx],'/features.rda',sep='')) }
fME <- sapply(1:dim(df)[1],existFeatureMat )

## preliminaries
pp   <- function(){paste(baseDir, 'ampOdd_click/word/VProbe/plots/',sep='')}
rda  <- function(){'~/Dropbox/ampOdd_click/Rdata/'}
trda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
 
ind <- which(df$drug=='N' & df$dataSet=='r' )#& !fME)# & df$probeString=='d' & df$animal%in%c('Walter'))
#ind <- which(df$drug=='N' & df$dataSet=='g' & fME)
#ind <- ind[31:length(ind)]
#ind <- 52
print(df[ind,])

Nclus                 <- 12

run_featureExtraction <- T
run_featureFitting    <- T
run_featureDecoding   <- F
run_icaAnalysis       <- F


## ==================================================================
#function_wavelet( ind[1] )
#disc <- lapply(ind, function_wavelet )


cluster               <- (Nclus > 1) & (length(ind)>1)


if (length(ind)==1){
  if (run_featureExtraction)
    function_aoc_featureExtraction( ind[1] )
  
  if (run_featureFitting)
    function_aoc_featureFitting( ind[1] )
  
  if (run_featureDecoding)
    function_aoc_featureDecoding( ind[1] )
  
  if (run_icaAnalysis)
    function_aoc_icaAnalysis( ind[1] )
}


if (length(ind)>1 & !cluster & run_featureExtraction)
  disc <- lapply(ind, function_aoc_featureExtraction )

if (length(ind)>1 & !cluster & run_featureFitting)
  disc <- lapply(ind, function_aoc_featureFitting )

if (length(ind)>1 & !cluster & run_featureDecoding)
  disc <- lapply(ind, function_aoc_featureDecoding )

if (length(ind)>1 & !cluster & run_icaAnalysis)
  disc <- lapply(ind, function_aoc_icaAnalysis )


if (length(ind)>1 & cluster){
  nClus <- min( c(length(ind),  Nclus, parallel::detectCores()-2))
  parallelCluster <- parallel::makeCluster(nClus,outfile='')
  print(parallelCluster)
  
  if (run_featureExtraction){
    tryCatch(
      disc <- parallel::parLapply(parallelCluster, ind, function_aoc_featureExtraction), error = function(e) print(e)
    )
  }
  
  if (run_featureFitting){
    tryCatch(
      disc <- parallel::parLapply(parallelCluster, ind, function_aoc_featureFitting), error = function(e) print(e)
    )
  }
  
  if (run_featureDecoding){
    tryCatch(
      disc <- parallel::parLapply(parallelCluster, ind, function_aoc_featureDecoding), error = function(e) print(e)
    )
  }
  
  if (run_icaAnalysis){
    tryCatch(
      disc <- parallel::parLapply(parallelCluster, ind, function_aoc_icaAnalysis), error = function(e) print(e)
    )
  }
  
  # Shutdown cluster neatly
  if(!is.null(parallelCluster)) {
    parallel::stopCluster(parallelCluster)
    parallelCluster <- c()
  }
}