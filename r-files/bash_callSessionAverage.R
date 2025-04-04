## preliminaries
args <- commandArgs( trailingOnly=TRUE )

baseDir <- '~/Dropbox/'
source(paste(baseDir, 'ampOdd_click/r-files/source_all_functions.R',sep=''))
#rda <- function(){return('/Volumes/rawData/EEG/ampOddClick/rda/')}

dataSet    <- 'g'
speaker    <- 'TB'
cmpStrList <- c('p1x.r','p1a.r','p1b.r','nx.r','p1x.r','p1.r','n1.r','p2.r')

if (length(args)>=2)
  dataSet <- as.character(args[2])

if (length(args)>=3)
  speaker <- as.character(args[3])

useInd <- which(df$speaker == speaker & df$dataSet==dataSet & df$animal==as.character(args[1]) )

if(length(args)>=4)
  useInd <- args[4:length(args)]

useInd
disc <- lapply(useInd,function(n){ callMaster_session(df[n,], whatList=c('raw'),fltxtn='' ) } )

#disc <- lapply(useInd,function(n){ callMaster_session(df[n,], whatList=c('raw','raw5k'),fltxtn='' ) } )
#disc <- lapply(useInd,function(n){ callMaster_session(df[n,], whatList=c('llp2','mlp','mbp','bsp'),fltxtn='_s' ) } )



