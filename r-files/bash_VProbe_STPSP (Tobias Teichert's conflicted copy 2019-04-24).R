## master VProbe STPSP
rm(list=ls())
defpar <- par()
library(parallel)

baseDir <- '~/Dropbox/'
source(paste(baseDir, 'RSkernel/r-files/functions/prepDFVprobe.R',sep=''))
#source(paste(baseDir, 'RSkernel/r-files/functions/source_all_functions.R',sep=''))
#source( paste(baseDir,'RSkernel/r-files/functions/functions_RSdual.R',sep='') )
#source( paste(baseDir,'RSkernel/r-files/functions/functions_VProbe_STPSP.R',sep='') )
source( paste(baseDir,'ampOdd_click/r-files/functions/FLA_VProbe_STPSP.R',sep='') )


colClasses <- c('character',rep('numeric',2),'character',rep('numeric',5)  ,rep('character',3),rep('numeric',3) )
df         <- read.table(paste(baseDir,'ampOdd_click/ampOddclickVProbe_link.txt',sep=''),header=T,colClasses=colClasses)
df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))

df <- prepDFVprobe(df)

## preliminaries
##pp <- function(){paste(baseDir, 'RSkernel/Figures/',sep='')}
pp   <- function(){paste(baseDir, 'ampOdd_click/word/VProbe/plots/',sep='')}
rda  <- function(){'~/Dropbox/ampOdd_click/Rdata/'}
trda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}

ComputerName <- 'tinyone' #system(paste('scutil --get ComputerName'),intern=T)


if(identical(ComputerName,'Fry')){
  ind     <- c(1)
}


if( identical(ComputerName,'zoidberg') ){
  ind <- which(df$drug=='N')# & df$animal%in%c('Sam'))
}

if(identical(ComputerName,'Farnsworth')){
  ind <- which(df$drug=='N')# & df$animal%in%c('Walter'))
}

ind <- which(df$drug=='N' & df$dataSet=='g' & df$animal%in%c('Sam','Walter'))
#ind <- ind[15]
#ind <- c(10)#,7)
print(df[ind,])


## ==================================================================
if (length(ind)==1)
  FLA_VProbe_STPSP(1, df=df, baseDir=baseDir,trda=trda, ind=ind )

cluster <- T


if (length(ind)>1 & !cluster)
  disc <- lapply(2:length(ind), FLA_VProbe_STPSP, df=df, baseDir=baseDir,trda=trda, ind=ind )
  

if (length(ind)>1 & cluster){
  nClus <- min( c(length(ind),parallel::detectCores()-1))
  parallelCluster <- parallel::makeCluster(nClus,outfile='')
  print(parallelCluster)
  
  tryCatch(
    disc <- parallel::parLapply(parallelCluster, 1:length(ind), FLA_VProbe_STPSP,df=df, baseDir=baseDir,trda=trda, ind=ind ), error = function(e) print(e)
  )
  
  # Shutdown cluster neatly
  if(!is.null(parallelCluster)) {
    parallel::stopCluster(parallelCluster)
    parallelCluster <- c()
  }
}