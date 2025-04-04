## master VProbe STPSP
rm(list=ls())
defpar <- par()
library(parallel)

baseDir <- '~/Dropbox/'
#source(paste(baseDir, 'RSkernel/r-files/functions/prepDFVprobe.R',sep=''))
source(paste(baseDir, 'ampOddClick/r-files/functions/prepDFVprobe.R',sep=''))
#source(paste(baseDir, 'RSkernel/r-files/functions/source_all_functions.R',sep=''))feat
#source( paste(baseDir,'RSkernel/r-files/functions/functions_RSdual.R',sep='') )
#source( paste(baseDir,'RSkernel/r-files/functions/functions_VProbe_STPSP.R',sep='') )
source( paste(baseDir,'ampOddClick/r-files/functions/FLA_VProbe_STPSP.R',sep='') )

curated <- T
if (!curated){
  ## read un-curated file
  colClasses <- c('character',rep('numeric',9),rep('character',3),rep('numeric',1) )# ,rep('character',3),rep('numeric',3) )
  df         <- read.table(paste(baseDir,'ampOddClick/ampOddclickdual_link.txt',sep=''),header=T,colClasses=colClasses)
  #df$probe   <- 1
  df$Layer   <- 10
  df$spacing <- 100
  df$drug    <- 'N'
}

if (curated){
  # read the curated file
  colClasses <- c('character',rep('numeric',2),'character',rep('numeric',5)  ,rep('character',3),rep('numeric',9) )
  df         <- read.table(paste(baseDir,'ampOddClick/ampOddclickVProbe_link.txt',sep=''),header=T,colClasses=colClasses)
}

df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))
df         <- prepDFVprobe(df)

## preliminaries
##pp <- function(){paste(baseDir, 'RSkernel/Figures/',sep='')}
pp   <- function(){paste(baseDir, 'ampOddClick/word/VProbe/plots/',sep='')}
rda  <- function(){'~/Dropbox/ampOddClick/Rdata/'}
#trda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
trda <- function(){'/Volumes/EEGCyno/EEG/ampOddClick/'}



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

##first check for raw, llp, mua preprcoessing
existraw<-function(x){file.exists(paste0("/Volumes/Drobo5D3/EEG/ampOddClick/rda/", df$session[x], "/", df$session[x], "_ch34_raw.mat"))}
#existllp<-function(x){file.exists(paste0("/Volumes/Drobo5D3/EEG/ampOddClick/rda/", df$session[x], "/", df$session[x], "_ch34_llp.mat"))}
existmua<-function(x){file.exists(paste0("/Volumes/Drobo5D3/EEG/ampOddClick/rda/", df$session[x], "/", df$session[x], "_ch34_mua.mat"))}

#df[which(unlist(lapply(1:dim(df)[1], existraw))+unlist(lapply(1:dim(df)[1], existmua))!=2),]


existRscalar<-function(x){file.exists(paste0("/Volumes/Drobo5D3/EEG/ampOddClick/rda/", df$session[x], "/RSscalar_ch34_MUA1_all.rda"))}
ind<-intersect(which(!unlist(lapply(1:dim(df)[1], existRscalar))), which(unlist(lapply(1:dim(df)[1], existraw))+unlist(lapply(1:dim(df)[1], existmua))==2))

#ind<-which(unlist(lapply(1:dim(df)[1], existraw))+unlist(lapply(1:dim(df)[1], existmua))==2)


ind <- which(df$drug=='N' & df$dataSet=='r' )#& !df$probeString=='d')# & df$animal%in%c('Walter'))
#ind <- which(df$drug=='N' & df$dataSet=='g' & df$animal%in%c('Sam') & df$preProc==1)
ind <- which(df$drug=='N' & df$dataSet%in%c('g','r') & df$animal%in%c('Walter'))
#ind <- which(df$numProbe==2 & df$dataSet%in%c('g','r'))# & df$animal%in%c('Walter'))
#ind <- ind[15:length(ind)]
ind <- c(55) # nice recording Sam_20190318
ind <- c(82) # 82 85 86last Walter recording
#ind <- 39:46
#ind <- ind[40:length(ind)]
#ind <- c(10)#,7)

ind <- which(df$drug=='N' & df$dataSet%in%c('g','r') & df$animal%in%c('Fuzzy'))

print(df[ind,])

##new only
#ind<-ind[32:length(ind)]
#ind<-c(86, 87, 88, 89, 90, 93)

## ==================================================================
if (length(ind)==1)
  FLA_VProbe_STPSP(1, df=df, baseDir=baseDir,trda=trda, ind=ind )

cluster <- T

if (length(ind)>1 & !cluster)
  disc <- lapply(1:length(ind), FLA_VProbe_STPSP, df = df, baseDir = baseDir, trda = trda, ind = ind )


#errorlist<-c()
#for(i in 1:length(ind)){
#  tryCatch(FLA_VProbe_STPSP(i, df = df, baseDir = baseDir, trda = trda, ind = ind), error = function(x){errorlist<<-c(errorlist, i)})
#}


if (length(ind)>1 & cluster){
  nClus <- min( c(length(ind),  6, parallel::detectCores()-2))
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