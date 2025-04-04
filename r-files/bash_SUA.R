## master_noModel
## this script calls all the important first-level model free analyses
args <- commandArgs( trailingOnly=TRUE )

efpar <- par()

baseDir <- '~/Dropbox/'
source(paste(baseDir, 'ampOdd_click/r-files/source_all_functions.R',sep=''))
source(paste(baseDir, 'ampOdd_click/r-files/functions/call_FLA_SUA.R',sep=''))

.GlobalEnv$spd       <- function(){paste('/Volumes/Drobo5D3/EEG/sortedZJ_ampOddClick/',sep='')}
.GlobalEnv$trda      <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}


## load df and chdf
colClasses <- c('character',rep('numeric',2),'character',rep('numeric',5)  ,rep('character',3),rep('numeric',3) )
df         <- read.table(paste(baseDir,'ampOdd_click/ampOddclickVProbe_link.txt',sep=''),header=T,colClasses=colClasses)
df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))

suaMethod        <- 'nev'
startFromScratch <- TRUE

colClasses <- c('character',rep('numeric',11),rep('character',3),'numeric','character','numeric','character','character','character')
chdf       <- read.table(paste(spd(),'ampOddClick_SUA_ZJ_link.txt',sep=''),header=T,colClasses=colClasses)
chdf       <- prep_chdf_SUA(chdf)

fileExists      <- function(n){return(file.exists(paste( spd(), chdf$session[n],'/ch_', chdf$channel[n],'/',chdf$cell[n],'_nev_',chdf$fextn[n],'_shape.mat',sep='') )) }
chdf$fileExists <- unlist( lapply(1:dim(chdf)[1], fileExists ))


## ========================================================
alreadyRun <- function(n){ return( file.exists( paste( rda(),'/STPSP/STPSPfit_',chdf$session[n],'_all_FR_ch_',chdf$channel[n],'_',chdf$cell[n],'_',chdf$fextn[n],'_sig1_gd.rda',sep='') )) }
chdf$alreadyRun <- unlist( lapply(1:dim(chdf)[1], alreadyRun ))
if (startFromScratch)
  chdf$alreadyRun <- c(FALSE)

chdf <- chdf[ which(chdf$dataSet=='g' & !chdf$alreadyRun & chdf$fileExists),]
options(warn=0)

basepp     <- function(){paste(baseDir, '/ampOdd_click/word/SUA/plots/',sep='')}
cmpStrList <- c('sig1')

dataSet    <- 'g' 
whoList    <- c('all')#,'regular','random')

ind     <- 1:length(chdf$session)
cluster <- T

# ========================= test-call without cluster
#call_FLA_SUA(3 ,df=chdf, baseDir=baseDir, ind=ind, whoList=whoList, dataSet=dataSet, cmpStrList=cmpStrList, suaMethod=suaMethod )


if (!cluster)
  disc <- lapply(1:length(ind), call_FLA_SUA,df=chdf, baseDir=baseDir, ind=ind, whoList=whoList, dataSet=dataSet, cmpStrList=cmpStrList, suaMethod=suaMethod )

if(cluster){
  parallelCluster <- parallel::makeCluster( min( c(46, length(ind), parallel::detectCores()-2) ) ,outfile='')
  print(parallelCluster)
  
  tryCatch(
    disc <- parallel::clusterApplyLB(parallelCluster, 1:length(ind), call_FLA_SUA,df=chdf, baseDir=baseDir, ind=ind, whoList=whoList, dataSet=dataSet, cmpStrList=cmpStrList, suaMethod=suaMethod ), error = function(e) print(e)
  )
  
  # Shutdown cluster neatly
  if(!is.null(parallelCluster)) {
    parallel::stopCluster(parallelCluster)
    
    parallelCluster <- c()
  }
}
