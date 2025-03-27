library(R.matlab)
library(lattice)#for eps.file
try(library(car))   # was commented out... not sure why. necessary for averagR; can't be installed on tinyone
library(signal)
library(akima)# for interp
library(Rwave)# for cwt
library(scales)

repDir <- '~/pitt-MePhysGit'

try(library(fastICA))
try( source( paste(repDir,'/RSkernel/r-files/functions/call_fastICA.R',sep='') ) )
try( source( paste(repDir,'/RSkernel/r-files/functions/save_fastICA.R',sep='') ) )
try( source( paste(repDir,'/RSkernel/r-files/functions/removeICAcomponent.R',sep='') ) )

try( source( paste(repDir,'/RSkernel/r-files/lamiarLFPtheory/showICAcomponents.R',sep='') ) )

try(library(xlsx))
#try(library(xlsx2)) # why is xslx2 not being loaded?

#source(repDir,'/r-compile/psyFit/R/psyfit.R')
source(repDir,'/r-compile/psyfit/R/psyfit.R')

source(repDir,'/r-compile/rUtils/R/rUtils.R')
source(repDir,'/r-compile/ePhysLab/R/ePhysLab.R')
source(repDir,'/r-compile/mav.test/R/mav.test.R')
source(repDir,'/r-compile/stereotax/R/stereotax.R')


source(repDir,'/RSkernel/r-files/functions/tapered_fft.R')
source(repDir,'/RSkernel/r-files/functions/tapered_ccf.R')

# files to be included in the signalProcessor pacakge:
source( paste(repDir,'/changeDetection/r-files/functions/plotAnatomy.R',sep='') )

#source(repDir,'/r-compile/signalProcessor/R/wrapper_figures/wrapper_TRF_HeatPlot.R')
source(repDir,'/r-compile/SignalProcessoR/R/wrapper_figures/wrapper_TRF_HeatPlot.R')

#source(repDir,'/r-compile/signalProcessor/R/wrapper_figures/wrapper_STRF_HeatPlot.R')
source(repDir,'/r-compile/SignalProcessoR/R/wrapper_figures/wrapper_STRF_HeatPlot.R')

#source(repDir,'/r-compile/signalProcessor/R/wrapper_figures/wrapper_CWT_HeatPlot.R')
source(repDir,'/r-compile/SignalProcessoR/R/wrapper_figures/wrapper_CWT_HeatPlot.R')



## Change Log
## TT 20211019 prepXPP changeDetection getCSOA: removed a factor 0.001 to cause CSOA to be reported in seconds
## TT 20211019 loadSignalProcessor_session eventStrobe reports time in ms, events in seconds. added a factor of 001 to eventStrobe[,1]


## NchanBankList
#NchanBankList        <- c( 33,  24,  32,  32,  64, 96,   8)
#names(NchanBankList) <- c('e',  'v', 't', 's', 'c', 'u', 'y')

NchanBankList        <- c( 33,  24,  32,  32,  32,  64, 64,  96,   8,  128,  0,  32, 16)
names(NchanBankList) <- c('e',  'v', 't', 'd', 'c', 's', 'S', 'u', 'y', 'b','q', 'z', 'Z')

prepDFsignalProcessor <- function(project,linkFileName=NULL){
  
  if (is.null(linkFileName))
    linkFileName    <- paste(project,'_link.txt',sep='')
  
  projectDir  <- paste(repDir,'/',project,'/',sep='')
  
  print( paste(projectDir,linkFileName,sep='') )   
  df         <- read.table(paste(projectDir,linkFileName,sep=''),header=T)
  df$project <- project
  df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(as.character(df$session)[n],'_')[[1]][1] }  ))
  df$link    <- linkFileName#strsplit(linkFileName,'_link')[[1]][1]
  df$session <- as.character(df$session)
  
  # identify recordings that are separated into two or more files
  getDateTime    <- function(n,m=1,tdf=df){ strsplit(as.character(tdf$session[n]),'_')[[1]][m] }
  df$animal      <- sapply(1:dim(df)[1],getDateTime, m=1, tdf=df)
  df$date        <- sapply(1:dim(df)[1],getDateTime, m=2, tdf=df)
  df$time        <- sapply(1:dim(df)[1],getDateTime, m=3, tdf=df)
  df$depthStr    <- sapply(1:dim(df)[1],getDateTime, m=4, tdf=df)
  
  getYear <- function(n,tdf=df){ substr(tdf$date[n],1,4) }
  df$year <- sapply(1:dim(df)[1], getYear, tdf=df)
  
  getMonth <- function(n,tdf=df){ substr(tdf$date[n],5,6) }
  df$month <- sapply(1:dim(df)[1], getMonth, tdf=df)
  
  
  getDepth       <- function(n,tdf=df){ neg <- (substr(tdf$depthStr[n],2,2)=='m'); res <- ifelse(neg, -1 * as.numeric(substr(tdf$depthStr[n],3,10)), as.numeric(substr(tdf$depthStr[n],2,10)));return(res) }
  df$depth       <- sapply(1:dim(df)[1],getDepth, tdf=df)
  df$day         <- paste(df$animal, df$date, sep='_')
  df$index       <- 1:dim(df)[1]
  
  if (!'probe'%in%names(df))# is.null(df$probe))
    df$probe <- 1
  
  df$probe[is.na(df$probe) ] <- 1
  
  df$ElectrodeConfig <- as.character(df$ElectrodeConfig)
  
  extractProbeStr  <- function(n){dfline <- df[n,];substr(dfline$ElectrodeConfig,dfline$probe+1,dfline$probe+1)}
  df$probeString   <- sapply(1:dim(df)[1],extractProbeStr)
  df$numChan       <- NchanBankList[df$probeString]
  
  df$startVIndex <- ifelse(df$animal=="Rockey", 22, 34)
  extractProbeStrM <- function(n,m=2){dfline <- df[n,];substr(dfline$ElectrodeConfig,m,m)}
  
  if (max(df$probe)>1){
    for(rx in 2:max(df$probe) ){
      ndx                 <- which(df$probe == rx)
      df$startVIndex[ndx] <- df$startVIndex[ndx] + NchanBankList[ extractProbeStrM(ndx, m=rx) ]  
    }
  }
  
  
  
  df$prelimStartMIndex <- ifelse(df$animal=="Rockey", 22, 34)
  
  df$ElectrodeConfig[which(is.na(df$ElectrodeConfig)) ] <- 'e'

  if(!is.null(df$ElectrodeConfig)){
    tv1<-c("v", "t", "d", "s","S", "c", "u", "z","Z", "a", "b", "e", "y", "x", "?", "Q")
    tv2<-c(24,  32,  32,  64,  64, 32,  96,  32,  16, 32,  128, 0,   0,   0,   0,    0)
    
    
    getNumProbes<-function(n){
      tmp<-0
      for(i in 1:nchar(df$ElectrodeConfig[n])){
        tmp <- tmp + ifelse(substr(df$ElectrodeConfig[n], i, i)%in%c("v", "t", "d", "s", "c", "u", "z","Z"), 1, 0) 
      }
      return(tmp)
    }
    df$numProbes <- sapply(1:dim(df)[1], getNumProbes)

    tf<-function(x){
      tmp<-0
      for(i in 1:nchar(x)){
        tmp<-tmp+tv2[which(substr(x, i, i)==tv1)]
      }
      return(tmp)
    }
    df$prelimEndMIndex <- df$prelimStartMIndex+unlist(lapply(as.character(df$ElectrodeConfig), tf))-1
    
    #df$startVIndex <- df$prelimStartMIndex
    #df$numChan <- df$prelimEndMIndex - df$prelimStartMIndex + 1
    
    
    mesoInd <- which(df$animal=='Soleil' & df$ElectrodeConfig%in%c('ebbby','ebbbbbbby'))
    df$prelimStartMIndex[mesoInd] <- 129
    df$startVIndex[mesoInd] <- 1
    df$numChan[mesoInd] <- 513
    df$prelimEndMIndex[mesoInd] <- 513
  }
  
  return(df)
}

# list pos possible locations to search if remote_rda is not specified:
remote_rda_list    <- list()
remote_rda_list[1] <- '/Volumes/EEGCyno/'
remote_rda_list[2] <- '/Volumes/Drobo5D3/'
remote_rda_list[3] <- '~/localDATA/'

checkSignalProcessor <- function(df){
  
  # Stage A ==== [check raw data]
  # Define locations to look for raw data =====
  pathsA <- c('/Volumes/rawDataCyno/EEG/', '/Volumes/rawData/EEG/', '/Volumes/EEGLab/rawDataOverflow')
  names(pathsA) <- c('rawCyno','rawDrobo','rawOverflow')
  resultsA <- matrix(data = NA, nrow = dim(df)[1], ncol = length(pathsA))
  
  for (m in 1:length(pathsA)){
    for (p in 1:dim(df)[1]){
      resultsA[p,m]   <- dir.exists( paste(pathsA[m], '/', df$animal[p], '/', df$animal[p],'_', df$date[p],'/', df$session[p], sep = ''))
    }
  }
  Adf            <- data.frame(resultsA)
  names (Adf)    <- names(pathsA)
  getRawNdx      <- function(n){res <- which(resultsA[n,]); if(length(res)>1){res <- res[1]}; if(length(res)==0){res=0}; return(res) }
  Adf$rawDataDir <- c('noRaw',names(pathsA))[ 1+unlist(sapply(1:dim(Adf)[1], getRawNdx )) ]
  df             <- data.frame(df, Adf  )
  
  # Define locations to look for EEG data =====
  pathsA <- c('/Volumes/EEGCyno/EEG/', '/Volumes/Drobo5D3/EEG/')
  names(pathsA) <- c('EEGCyno','EEGDrobo')
  resultsA <- matrix(data = NA, nrow = dim(df)[1], ncol = length(pathsA))
  
  for (m in 1:length(pathsA)){
    for (p in 1:dim(df)[1]){
      resultsA[p,m]   <- dir.exists( paste(pathsA[m], '/',df$project[p],'/', df$session[p], sep = ''))
    }
  }
  Adf            <- data.frame(resultsA)
  names (Adf)    <- names(pathsA)
  Adf$EEGDir     <- c('noEEG',names(pathsA))[ 1+unlist(sapply(1:dim(Adf)[1], getRawNdx )) ]
  df             <- data.frame(df, Adf  )
  
  # Look for Error Messages ======
  errorLogDirList  <- c( '/Volumes/EEGCyno/EEG/', '/Volumes/EEGLab/EEGLab/' ); names(errorLogDirList) <- c('Cyno','Drobo')
  errorLogNameList <- c('errorLogR','errorLogC','errorLogM','errorLogA','errorLogV')
  
  existErrorLog  <- function(n, dir='/Volumes/EEGCyno/EEG/', name='errorLogR'){ file.exists(paste(dir,'/',df$project[n],'/', df$session[n],'/',name, '.txt', sep='' ) )  }
  existErrorLog2 <- function(mstr, dr='/Volumes/EEGCyno/EEG/'){ res <- sapply( 1:dim(df)[1], existErrorLog, name=mstr, dir=dr )}
  
  existErrorLogByFolder <- function(fx, nameList=c('errorLogR','errorLogC','errorLogM','errorLogA','errorLogV'), whatStr='Error'){
    tst <- sapply( nameList, existErrorLog2, dr=errorLogDirList[fx] )
    rdf <- data.frame(tst,apply(tst,1, sum))
    
    dimnames(rdf)[[2]] <- list(paste( names(errorLogDirList)[fx], c(nameList,paste('any',whatStr,sep='')), sep='_'  ) )[[1]]
    return(rdf)
  }
  
  errLogC <- existErrorLogByFolder( 1, nameList = errorLogNameList )
  errLogD <- existErrorLogByFolder( 2, nameList = errorLogNameList )
  
  
  ## look for successX.txt files
  successLogNameList <- c('success','successC','successM','successA','successV')
  successLogC        <- existErrorLogByFolder( 1, nameList = successLogNameList, whatStr='success' )
  
  
  df <- data.frame(df, errLogC, errLogD, successLogC)
  
  # look for bhv file
  existbhvFile  <- function(n, extn='.bhv'){ file.exists(paste('/Volumes/rawDataCyno/EEG/bhv/', df$session[n],extn, sep='' ) )  }
  df$existbhv1 <- sapply( 1:dim(df)[1], existbhvFile, extn='.bhv')
  df$existbhv2 <- sapply( 1:dim(df)[1], existbhvFile, extn='.bhv2')
  df$existbhvm <- sapply( 1:dim(df)[1], existbhvFile, extn='.mat')
  
  df$existbhv <- df$existbhv1 | df$existbhv2 | df$existbhvm
  
  
  if (1==0){
    # Stage B =====
    #Define locations to looks for signal processor data 
    StageB <- c('EEGLab','EEG','VIS','AVG')
    StageBloc <- c('Cyno','Drobo')
    pathsB <- data.frame(matrix(data = NA, nrow = length(StageB), ncol =length(StageBloc)))
    names(pathsB) <- StageBloc
    row.names(pathsB) <- StageB
    pathsB$Cyno  <- c(NA,
                      paste('/Volumes/EEGCyno/EEG/', projectName, sep = ''),
                      NA, 
                      paste('/Volumes/EEGCyno/AVG/', projectName, '/', sep, '/' =''))
    pathsB$Drobo <- c(paste('/Volumes/EEGLab/EEGLab/', projectName, '/', sep = ''), 
                      paste('/Volumes/Drobo5D3/EEG/', projectName, '/', sep = ''),
                      paste('/Volumes/Drobo5D3/signalProcessorPlots/', projectName, '/', sep = ''),
                      NA)
    print(pathsB)
    #Sort df$session into matrices to store results (list, 1 matrix per location)
    resultsB <- vector(mode ="list", length = length(StageBloc))
    names(resultsB) <- StageBloc
    for(a in 1:length(resultsB)){
      resultsB[[a]] <- data.frame(matrix(data = NA, ncol = length(StageB), nrow =dim(df)[1]))
      row.names(resultsB[[a]]) <- df$session
      names(resultsB[[a]]) <- StageB
    }
    #Loop through all paths and report if files exist
    for(m in 1:dim(pathsB)[2]){
      for(n in 1:dim(pathsB)[1]){
        for(p in 1:dim(df)[1]){
          .GlobalEnv$directory   <- function(){paste(pathsB[n,m], df[p,1], sep = '')}
          resultsB[[m]][p,n] <- dir.exists(path = directory())
          if(is.na(pathsB[n,m]))
            resultsB[[m]][p,n] <- NA
        }
      }
    }
    print(resultsB)
    
  }
  
  return(df)
}

loadSignalProcessor <- function(tdf, project=NA, epoch='llp', analysis=NULL, chStr='frontoCentral',showCh=1, show=F, nike=F, 
                                local_rda='/Volumes/Drobo5D3/EEGLab/', remote_rda=NULL, eeglab_rda='/Volumes/EEGLab/',checkFiles=F,
                                extn='',xppextn='',blrng=c(-50,0), suarng=c(-1,1), suaMethod='nev', what=NULL, chBlkN=16, appendByTrial=T, chan2var=F){
  
  if (is.na(project))
    project <- tdf$project[1]
  
  passChStr <- chStr
  if (is.na(chStr[1]))
    chStr <- paste('ch',tdf$chNr,sep='')
  
  secondLevelFunction <- loadSignalProcessor_session
  if (!is.null(analysis))
    secondLevelFunction <- loadSignalProcessor_avg
  
  if (!is.null(analysis) & nike)
    secondLevelFunction <- loadSignalProcessor_nike
  
  
  # check which of the possible rda diretories are mounted:
  if (is.null(remote_rda)){# try to identify the correct remote_rda from remote_rda_list
    rda_access_list <- sapply( remote_rda_list, function(xstr){file.access(xstr)} )
    
    if (sum(rda_access_list==0)==0)# no remote directory mounted:
      stop('none of the potential data directories is mounted')
    
    if (sum(rda_access_list==0)==1)# only one remote directory mounted:
      remote_rda <- remote_rda_list[[ which(rda_access_list==0) ]]
    
    if (sum(rda_access_list==0)>1){# more than one potential data directory mounted
      folderExist <- function(n){# check in which folder, if any, the EEG data resides 
        return(sum(lapply( paste(remote_rda_list[n],'/EEG/', tdf$project, '/', tdf$session, sep = ''), function(xstr){file.access(xstr)})==0))  
      }
      if (!is.null(analysis))# look for data in the AVG folder instead
        folderExist <- function(n){# check in which folder, if any, the EEG data resides 
          return(sum(lapply( paste(remote_rda_list[n],'/AVG/', tdf$project, '/', tdf$session, sep = ''), function(xstr){file.access(xstr)})==0))  
        }
      
      NExist <- sapply( which(rda_access_list==0), folderExist )
      if(max(NExist)==0)
        stop('data was not found in any of the directories')
      remote_rda <- remote_rda_list[[ which(NExist==max(NExist))[1] ]]
    }
  }
  print(remote_rda)
  
  # check for the existence of relevant files and exclude sessions that have not yet been processed
  if (!file.access(remote_rda)==0)
    stop('data directory not mounted')
  
  # first-pass check to see if the requested sessions have been processed:
  # fileExist <- function(n){ file.exists( paste(local_rda,'/',project,'/',tdf$session[n],'/eventStrobe.txt',sep='')) }
  # fileExist2 <- function(n){ file.exists( paste(remote_rda, '/EEGLab/', project,'/', tdf$session[n],'/eventStrobe.txt',sep='')) }
  #These functions seem never to point to the eventStrobe file, and are identically named;
  #they are crashing every time; are they necessary?
  
  fileExist <- function(n){# looking for the EEGLab folder: 
    
    if (remote_rda == remote_rda_list[[1]])
      dirName   <- paste(remote_rda,'/EEG/', project, '/', tdf$session[n], sep = '')
    
    if (remote_rda == remote_rda_list[[2]])
      dirName   <- paste(remote_rda,'/EEGLab/', project, '/', tdf$session[n], sep = '')
    
    if (remote_rda == remote_rda_list[[3]])
      dirName   <- paste(remote_rda,'/EEG/', project, '/', tdf$session[n], sep = '')
    
    #if(!dir.exists(dirName)){
    #  paste0('\n This directory does not exist: \n', dirName)
    #  dirName   <- paste(eeglab_rda,'/EEGLab/', project, '/', tdf$session[n], sep = '')
    #}
    
    #if(!dir.exists(dirName)){
    #  paste0('\n This directory does not exist: \n', dirName)
    #  dirName <- paste0('/Volumes/EEGLab/EEGLab/', project, '/', tdf$session[n])
    #}
    
    #if(!dir.exists(dirName)){
    #  paste0('\n This directory does not exist: \n', dirName)
    #  dirName <- paste0(remote_rda,'/EEG/', project, '/', tdf$session[n])
    #}
    
    if(!dir.exists(dirName)){
      paste0('\n This directory does not exist: \n ', dirName)
      warning('\n No EEGLab folder to load the data from. \n')
    }
    
    paste0('Working data directory: ', dirName)
    
    fileName1  <- paste(dirName, '/eventStrobe.txt', sep = '')
    fileName1b <- paste(remote_rda, '/EEG/',project, '/',tdf$session[n], '/eventStrobe.txt', sep = '')
    fileName2  <- paste(dirName, '/success.txt', sep = '')
    
    res1 <- file.exists( fileName1 ) | file.access( fileName1 )==0 | file.exists(fileName1b)
    res2 <- file.exists( fileName2 ) | file.access( fileName2 )==0
    
    if (tdf$project[n]%in%c('resting','tDCS'))
      res1 <- res2                          # resting has no monkey-logic bhv file. use res2 instead
    
    if (!res1)
      warning(paste(fileName1,' does not exist in the data directory.\n'))
    if (!res2)
      warning(paste(fileName2,' does not exist in the data directory.\n'))
    return(list(res1, res2))
  }
  
  if (!is.null(analysis))
    fileExist <- function(n){
      dirName   <- paste(remote_rda,'/AVG/', project, '/', tdf$session[n],'/',epoch,'/',analysis, sep = '')
      if(!dir.exists(dirName)){
        paste0('\n This directory does not exist: \n ', dirName)
        warning('\n No AVG folder to load the data from. \n')
      }
      
      paste0('AVG data directory: ', dirName)
      fileName1  <- paste(dirName, '/xpp.rda', sep = '')
      fileName2  <- paste(dirName, '/',                     chStr[1],'.rda', sep = '')
      fileName2b  <- paste(dirName, '/',                    chStr[1],'.mat', sep = '')
      fileName2c  <- paste(dirName, '/', tdf$animal[n],'_', chStr[1],'.mat', sep = '')
      
      res1  <- file.exists( fileName1 )  | file.access( fileName1 )==0 # for xpp
      res2  <- file.exists( fileName2 )  | file.access( fileName2 )==0 # for regular channels
      res2b <- file.exists( fileName2b ) | file.access( fileName2b )==0 # for nike import
      res2c <- file.exists( fileName2c ) | file.access( fileName2c )==0 # for nike import
      
      res <- res1 | res2 | res2b | res2c
      
      return(list(res,res))
    }
  
  fE        <- sapply(1:dim(tdf)[1], fileExist)

  missingSessions <- NULL
  if (sum(unlist(fE[1,]))<dim(fE)[2]){
    warning('some sessions have not yet been processed:')
    missingSessions <- tdf[!unlist(fE[1,]),] 
  }
  tdf <- tdf[unlist(fE[1,]),]
  if(dim(tdf)[1]==0){
    stop("No valid sessions.  Exiting.")
  }
  if(sum(unlist(fE[2,]))<dim(tdf)[1]){
    #    stop("Code has been added which fails if no eventStrobe.txt file exists for any of the specified sessions.  Exiting.")
    warning("A session was input with the EEGLab folder missing success.txt.  IF USING NEW PREPROCESSOR: Errors may be due to data not having been converted to R-readable format:")
    tdf[!unlist(fE[2,]),]
  }
  
  if (checkFiles)
    return(fE)
  
  # second pass test to see if the requested sessions have been averaged: ====
  if (!is.null(analysis)){
    avgExist <- function(n){ 
      thisChStr <- chStr
      if(is.null(chStr) & 'chStr'%in%names(tdf))
        thisChStr <- tdf$chStr[n]
      avgFileName1 <- paste(remote_rda,'/AVG/',project,'/',tdf$session[n],'/',epoch,'/',analysis,'/',                 thisChStr[length(thisChStr)],'.rda',sep='')
      avgFileName2 <- paste(remote_rda,'/AVG/',project,'/',tdf$session[n],'/',epoch,'/',analysis,'/',                  thisChStr[length(thisChStr)],'.mat',sep='')
      avgFileName3 <- paste(remote_rda,'/AVG/',project,'/',tdf$session[n],'/',epoch,'/',analysis,'/',tdf$animal[n],'_',thisChStr[length(thisChStr)],'.mat',sep='')
      
      
      #print(avgFileName)
      res1 <- file.exists( avgFileName1 ) | file.access( avgFileName1 )==0
      res2 <- file.exists( avgFileName2 ) | file.access( avgFileName2 )==0
      res3 <- file.exists( avgFileName3 ) | file.access( avgFileName3 )==0
      
      res <- res1 | res2 | res3
      
      avgFileName <- avgFileName1
      if (res2)
        avgFileName <- avgFileName2
      if (res3)
        avgFileName <- avgFileName3
      
      
      if (!res)
        warning(paste(avgFileName,' does not exist\n'))
      if (res)
        warning(paste(avgFileName,' does exist\n'))
      return(res)
    }
    
    fE <- sapply(1:dim(tdf)[1], avgExist)
    if (sum(fE)<length(fE)){
      warning('some averages have not yet been calculated...')
      print(tdf[!fE,])
    }
    tdf <- tdf[fE,]
  }
  
  print(tdf$session[1])
  print(remote_rda)
  
  ## READ EACH SESSION ====
  tone    <- secondLevelFunction( tdf[1,],project=project, epoch=epoch, analysis=analysis, chStr=passChStr,showCh=showCh,show=show, 
                                  local_rda=local_rda, remote_rda=remote_rda,eeglab_rda=eeglab_rda, 
                                  extn=extn,xppextn=xppextn,blrng=blrng, suarng=suarng, suaMethod=suaMethod, what=what,chBlkN=chBlkN)
  tone$xpp$sessionID       <- 1
  tone$trialInfo$sessionID <- 1
  
  if (chan2var)
    tone <- chan2varSignalProcessor( tone )
  
  #toneexp <- tone$exp
  
  if(dim(tdf)[1]>1){
    for (i in 2:dim(tdf)[1]){
      print(tdf$session[i])
      ttn   <- secondLevelFunction( tdf[i,], project=project, epoch=epoch, analysis=analysis, chStr=passChStr,showCh=showCh,show=show, 
                                    local_rda=local_rda, remote_rda=remote_rda,eeglab_rda=eeglab_rda, 
                                    extn=extn,xppextn=xppextn,blrng=blrng, suarng=suarng, suaMethod=suaMethod, what=what, chBlkN=chBlkN )
      ttn$xpp$sessionID       <- i
      ttn$trialInfo$sessionID <- i
      
      if (!appendByTrial)
        tone <- append.channels(tone,ttn)
      
      if (appendByTrial){
        if (chan2var)
          ttn <- chan2varSignalProcessor( ttn )
        
        
        if (!is.null(tone$lfp) & !is.null(ttn$lfp) )
          if ( dim(tone$lfp)[2]!=dim(ttn$lfp)[[2]] ){# very bizare. this seems to have been using == instead of !=, but it never crashed until now...?
            
            warning('read data data sets have different time-axis. Selecting time-points common to both data sets')
            
            validTonetaxis   <- which(tone$taxis %in% ttn$taxis )
            validTTNtaxis    <- which(ttn$taxis %in% tone$taxis )
            
            if ( min( c(length(validTTNtaxis),length(validTonetaxis)) )==0)
              browser()
            
            ttn$taxis <- ttn$taxis[validTTNtaxis]
            ttn$lfp   <- ttn$lfp[,validTTNtaxis,]
            
            tone$taxis <- tone$taxis[validTonetaxis]
            tone$lfp   <- tone$lfp[,validTonetaxis,]
          }
          
        tone <- append.data(tone,ttn)
      }
      
    }
  }
  
  if (chan2var ){
    tone$xpp$chIDbyDay <- tone$xpp$chID
    tone$xpp$chID      <- tapply(1:dim(tone$xpp)[1], list(tone$xpp$chIDbyDay, tone$xpp$sessionID))
  }
  
  tone$epoch <- epoch
  tone$project <- project
  
  if (!is.null(missingSessions)){
    print('missing sessions: ')
    print( missingSessions )
  }
  
  return(tone)
}


loadSignalProcessor_avg <- function( dfline, project='TRF', epoch='llp',analysis='sessionID', chStr='frontoCentral',showCh=1,show=F, 
                                     local_rda='~/Dropbox/', remote_rda='/Volumes/EEGCyno/',eeglab_rda=eeglab_rda, 
                                     extn='',xppextn='',blrng=c(-50,0), suarng=c(-1,1), suaMethod='nev', what=NULL,chBlkN=16){
  print('loading averages')
  # for backward compatibility
  if (!is.null(what))
    epoch <- what
  
  
  if (is.na(chStr[1]))
    chStr <- paste('ch',dfline$chNr,sep='')
  
  print(remote_rda)
  print(dfline)
  if (!file.access(remote_rda)==0){
    stop('data directory not mounted')
  }
  
  # enable loading of individual channels from chdf-type data ====
  if (is.null(chStr) & 'chStr'%in%names(dfline))
    chStr <- dfline$chStr
  
  # write out some info about the data contained ======
  tone             <- list()
  tone$project     <- project
  tone$chStr       <- chStr
  tone$epoch       <- epoch
  tone$analysis    <- analysis
  tone$dfline      <- dfline
  tone$exp         <- dfline$sessions
  # data directory ====
  rda            <- paste(remote_rda,'/AVG/',project,'/',dfline$session,'/',epoch,'/',analysis,'/',sep='')
  
  print(rda)
  if (! (file.access(rda) | file.exists(rda) ))
    stop('directory has not yet been processed')
  
  ### load xpp      ====  
  load( paste(rda,'xpp.rda',sep=''))
  tone$xpp <- xpp
  
  ### load time-axis ====
  load( paste(rda,'taxis.rda',sep=''))
  tone$taxis <- taxis
  
  tone$dt       <- mean(diff(tone$taxis)) # temporal resolution in ms
  tone$SR       <- 1000/tone$dt
  
  ### load channels ====
  loadChannel <- function(cx){
    load( paste(rda,chStr[cx],'.rda',sep=''))
    if (!is.array(chdata))
      dim(chdata) <- c(1,length(chdata),1)
    if (is.array(chdata) & length(dim(chdata))==2)
      dim(chdata) <- c(dim(chdata),1)
    return(chdata)
  }
  chDataList <- lapply(1:length(chStr), loadChannel )
  
  #  old ====
  #myenv<-new.env()
  #load(paste(rda,'xpp.rda',sep=''), env = myenv)
  #load("/Volumes/Drobo5D3/AVG/changeDetection/Jesse_20211007_1150_zm9999/llp/seqNr/xpp.rda", env = myenv)
  #ls(envir = myenv)
  #dim(get("xpp", envir = myenv))
  #pxpp<-get("xpp", envir = myenv)
  #dim(myenv$xpp)
  #myenv<-new.env()
  #load(paste(rda, chStr[cx], '.rda', sep = ''), env = myenv)
  #load("/Volumes/Drobo5D3/AVG/changeDetection/Jesse_20211007_1150_zm9999/llp/seqNr/ch2.rda", env = myenv)  #ls(envir = myenv)
  #dim(get("chdata", envir = myenv))
  
  #  for(i in 1:33){
  #    load(paste("/Volumes/Drobo5D3/AVG/changeDetection/Jesse_20211007_1150_zm9999/llp/seqNr/", chStr[i], '.rda', sep = ''), env = myenv)
  #    print(dim(get("chdata", envir = myenv)))
  #  }
  
  # re-shape matrix ====
  lfpdim    <- dim(chDataList[[1]])
  if (length(lfpdim)==3)
    lfpdim[3] <- length(chDataList)
  if (!length(lfpdim)==3)
    stop('something is wrong in loadSignalprocessor_avg')
  
  lfp      <- array( unlist(chDataList), lfpdim  )
  tone$lfp <- lfp
  
  #browser()
  
  lfpnbl <- tone$lfp
  if (length(blrng)==2){
    print('baseline correction')
    tmp            <- baseLineCorrection(tone$lfp,blrng, tone$taxis)  
    tone$lfp       <- tmp$res
    tone$xpp$blcor <- tmp$correction
  }
  
  dimnames(tone$lfp) <- list(trial=1:dim(tone$lfp)[1],time=tone$taxis,channel=chStr)
  
  
  # create log-timeaxis
  toff                                     <- 1
  suppressWarnings(tone$logtaxis           <- log2( tone$taxis + toff))
  tone$logtaxis[ (tone$taxis + toff <= 1)] <- .01*(tone$taxis[ which(tone$taxis + toff <= 1) ] - 1 + toff) 
  
  tone$xpp$valTone <- T
  tone$xpp$animal <- tone$dfline$animal
  
  
  tone$channelPosition <- data.frame(xpos=1+(0:(length(chStr)-1))%/%chBlkN, ypos=1+(0:(length(chStr)-1))%%chBlkN,  zpos=rep(0,length(chStr)) )
  tone$chDepthInd      <- 1:length(chStr)
  
  return(tone)
}


loadSignalProcessor_nike <- function( dfline, project='TRF', epoch='ffr',analysis='toneNr', chStr='frontoCentral',showCh=1,show=F, 
                                      local_rda='~/Dropbox/', remote_rda='/Volumes/EEGCyno/',eeglab_rda=eeglab_rda, 
                                      extn='',xppextn='',blrng=c(-50,0), suarng=c(-1,1), suaMethod='nev', what=NULL,chBlkN=16){
  
  # for backward compatibility
  if (!is.null(what))
    epoch <- what
  
  print(remote_rda)
  print(dfline)
  if (!file.access(remote_rda)==0){
    stop('data directory not mounted')
  }
  
  # enable loading of individual channels from chdf-type data ====
  if (is.null(chStr) & 'chStr'%in%names(dfline))
    chStr <- dfline$chStr
  
  # data directory ====
  #rda            <- paste(remote_rda,'/AVG/',project,'/',dfline$session,'/',sep='')
  rda            <- paste(remote_rda,'/AVG/',project,'/',dfline$session,'/',epoch,'/',analysis,'/',sep='')
  
  addSubjPrefix <- TRUE
  subjPrefix    <- paste(dfline$animal,'_',sep='')
  
  #if (chStr=='all'){
  if (length(chStr)==1){
    #chLst <- system(paste('ls ', rda, dfline$animal, '*.mat',sep=''), intern=T)
    chLst1 <- system(paste('ls ', rda, dfline$animal,'_', chStr ,'.mat',sep=''), intern=T)
    chLst2 <- system(paste('ls ', rda,                    chStr ,'.mat',sep=''), intern=T)
    
    # assuming subject prefix <- T
    extracChStr <- function(nx){ 
      tmp <- strsplit(chLst1[[nx]],paste(dfline$animal,'_',sep=''))[[1]][3]
      tmp2 <- strsplit(tmp, '.mat')
      return(tmp2)
    }
    chLst <- chLst1
    
    # what if no subject prefix:
    if ( length(chLst1)<length(chLst2) ){
      addSubjPrefix <- F
      subjPrefix <- ''
      
      
      extracChStr <- function(nx){ 
        tmp <- strsplit(chLst2[[nx]],rda)[[1]][2]
        tmp2 <- strsplit(tmp, '.mat')
        return(tmp2)
      }
      chLst <- chLst2
    }
    chStr <- unlist(sapply(1:length(chLst),extracChStr))
  }
  # write out some info about the data contained ======
  tone             <- list()
  tone$project     <- project
  tone$chStr       <- chStr
  tone$epoch       <- epoch
  tone$analysis    <- analysis
  tone$dfline      <- dfline
  tone$exp         <- dfline$sessions
  
  print(rda)
  if (! (file.access(rda) | file.exists(rda) ))
    stop('directory has not yet been processed')
  
  ### create xpp      ====  
  #load( paste(rda,'xpp.rda',sep=''))
  xpp <- data.frame(subject=dfline$session, toneNr=c(1,2,3,4))
  tone$xpp <- xpp
  
  ### load time-axis ====
  
  tst <- readMat( paste(rda, subjPrefix,chStr[1],'.mat',sep='') )
  convfctr <- ifelse(max(tst$t)<2,1000,1)
  
  tone$taxis <- convfctr*tst$t
  
  tone$dt       <- mean(diff(tone$taxis)) # temporal resolution in ms
  
  tone$SR       <- 1000/tone$dt
  
  ### load channels ====
  loadChannel <- function(cx){
    chdata <- t(readMat( paste(rda, subjPrefix,chStr[cx],'.mat',sep='') )$avgffr)
    
    if (!is.array(chdata))
      dim(chdata) <- c(1,length(chdata),1)
    if (is.array(chdata) & length(dim(chdata))==2)
      dim(chdata) <- c(dim(chdata),1)
    return(chdata)
  }
  chDataList <- lapply(1:length(chStr), loadChannel )
  
  # re-shape matrix ====
  lfpdim    <- dim(chDataList[[1]])
  if (length(lfpdim)==3)
    lfpdim[3] <- length(chDataList)
  if (!length(lfpdim)==3)
    stop('something is wrong in loadSignalprocessor_avg')
  
  lfp      <- array( unlist(chDataList), lfpdim  )
  tone$lfp <- lfp
  
  #browser()
  
  lfpnbl <- tone$lfp
  if (length(blrng)==2){
    print('baseline correction')
    tmp            <- baseLineCorrection(tone$lfp,blrng, tone$taxis)  
    tone$lfp       <- tmp$res
    tone$xpp$blcor <- tmp$correction
  }
  
  dimnames(tone$lfp) <- list(trial=1:dim(tone$lfp)[1],time=tone$taxis,channel=chStr)
  
  
  # create log-timeaxis
  toff                                     <- 1
  suppressWarnings(tone$logtaxis           <- log2( tone$taxis + toff))
  tone$logtaxis[ (tone$taxis + toff <= 1)] <- .01*(tone$taxis[ which(tone$taxis + toff <= 1) ] - 1 + toff) 
  
  tone$xpp$valTone <- T
  tone$xpp$animal <- tone$dfline$session
  
  
  tone$channelPosition <- data.frame(xpos=1+(0:(length(chStr)-1))%/%chBlkN, ypos=1+(0:(length(chStr)-1))%%chBlkN,  zpos=rep(0,length(chStr)), chStr=chStr )
  tone$chDepthInd      <- 1:length(chStr)
  
  tone <- resampleSignalProcessor( tone, seq(-50,by=.1,to=499.9) )
  tone$lfp[is.na(tone$lfp)] <- 0
  
  tone <- prepLFP(tone)
  
  return(tone)
}


loadSignalProcessor_session <- function( dfline, project='TRF', epoch='llp', analysis=NULL, chStr='frontoCentral',showCh=1,show=F, 
                                         local_rda='~/Dropbox/', remote_rda='/Volumes/Drobo5D3/',eeglab_rda=eeglab_rda, 
                                         extn='',xppextn='',blrng=c(-50,0), suarng=c(-1,1), suaMethod='nev', what=NULL,chBlkN=16){
  # for backward compatibility
  if (!is.null(what))
    epoch <- what
  
  ##  preliminaries                                                    ====
  print('loadSignalProcessor')
  
  print(remote_rda)
  if (!file.access(remote_rda)==0){
    stop('data directory not mounted')
  }
  
  # write out some info about the data contained                       ====
  tone          <- list()
  tone$project  <- project
  tone$chStr    <- chStr
  tone$epoch    <- epoch
  tone$dfline   <- dfline
  tone$exp      <- dfline$session
  
  exp    <- dfline$session
  from   <- ifelse(is.na(dfline$minTR),NA, dfline$minTR)
  to     <- ifelse(is.na(dfline$maxTR),NA, dfline$maxTR)
  
  #rda <- paste(remote_rda,'/EEG/',project,'/rda/',exp,'/',sep='')
  #if (!file.exists(rda))
  #  rda <- paste(remote_rda,'/EEG/',project,'/',exp,'/',sep='')
  
  
  rda <- paste(remote_rda,'/EEG/',project,'/',exp,'/',sep='')
  if (!file.exists(rda))
    rda <- paste(remote_rda,'/EEG/',project,'/rda/',exp,'/',sep='')
  
  
  print(paste( 'rda:', rda))
  
  # check files                                                        ====
  spd                <- paste(remote_rda,'/EEGLab/',project,'/',sep='')
  readSUA            <- exists('cell',dfline) & exists('channel',dfline)
  
  eventStrobeFile    <- paste(remote_rda, '/EEGLab/', project,'/', exp,'/eventStrobe.txt',sep='')
  eventTimingFile    <- paste(remote_rda, '/EEGLab/', project,'/', exp,'/events.txt',sep='')
  trialStrobeFile    <- paste(remote_rda, '/EEGLab/', project,'/', exp, '/eventStrobeTrial.txt',sep='')
  
  if(!file.exists(trialStrobeFile)){
    trialStrobeFile    <- paste(local_rda,              project,'/Rdata/',exp,'/eventStrobeTrial.txt',sep='')
  }
  
  
  # new EEGLab directory in volumes                                    ====
  if(!file.exists(eventStrobeFile)){
    eventStrobeFile    <- paste('/Volumes/EEGLab/EEGLab/', project, '/', exp, '/eventStrobe.txt', sep = '')
  }
  if(!file.exists(eventTimingFile)){
    eventTimingFile    <- paste('/Volumes/EEGLab/EEGLab/', project,'/', exp,'/events.txt',sep='')
  }
  if(!file.exists(trialStrobeFile)){
    trialStrobeFile    <- paste('/Volumes/EEGLab/EEGLab/', project, '/', exp, '/eventStrobeTrial.txt', sep = '')
  }
  
  # new Cyno EEG layout has ML files in same directory as intan files  ====
  if(!file.exists(eventStrobeFile)){
    eventStrobeFile    <- paste(remote_rda,'/EEG/', project, '/', exp, '/eventStrobe.txt', sep = '')
  }
  if(!file.exists(eventTimingFile)){
    eventTimingFile    <- paste(remote_rda,'/EEG/', project,'/', exp,'/events.txt',sep='')
  }
  if(!file.exists(trialStrobeFile)){
    trialStrobeFile    <- paste(remote_rda,'/EEG/', project, '/', exp, '/eventStrobeTrial.txt', sep = '')
  }
  
  
  ###mounting bug
  #if(!file.exists(eventStrobeFile)){
  #  eventStrobeFile    <- paste('/Volumes/EEGLab-1/EEGLab/', project, '/', exp, '/eventStrobe.txt', sep = '')
  #}
  #if(!file.exists(eventTimingFile)){
  #  eventTimingFile    <- paste('/Volumes/EEGLab-1/EEGLab/', project,'/', exp,'/events.txt',sep='')
  #}
  #if(!file.exists(trialStrobeFile)){
  #  trialStrobeFile    <- paste('/Volumes/EEGLab-1/EEGLab/', project, '/', exp, '/eventStrobeTrial.txt', sep = '')
  #}
  
  
  ## the Strobe files previously didn't live on Drobo, but this is easier than changing all the calls to signalprocessor
  if(!file.exists(eventStrobeFile)){
    eventStrobeFile    <- paste('~/Dropbox/', project,'/RData/', exp,'/eventStrobe.txt',sep='')
  }
  if(!file.exists(trialStrobeFile)){
    trialStrobeFile    <- paste('~/Dropbox/', project,'/RData/', exp, '/eventStrobeTrial.txt',sep='')
  }
  
  print('File Locations:')
  print(eventStrobeFile)
  print(eventTimingFile)
  print(trialStrobeFile)
  
  
  if (!file.exists(eventTimingFile))
    print(paste('missing:', eventTimingFile) )
  if (!file.exists(eventStrobeFile))
    print(paste('missing:', eventStrobeFile) )
  if (!file.exists(trialStrobeFile))
    print(paste('missing:', trialStrobeFile) )
  
  ## read MLtrials data if present                                     ====
  events <- array(1,c(1,2))
  if (file.exists(eventTimingFile)){
    print(eventTimingFile)
    events            <- read.table( eventTimingFile, header=F )
  }
  if (file.exists(eventStrobeFile)){
    print(eventStrobeFile)
    eventStrobe       <- read.table( eventStrobeFile, header=T )
  }
  if (file.exists(trialStrobeFile)){
    print(trialStrobeFile)
    trialStrobe       <- read.table( trialStrobeFile, header=T )
  }
  
  if ( dfline$project%in%c('resting','tDCS')){# resting does not have eventStrobe data. Create dummy version
    eventStrobe            <- events
    eventStrobe$eventTimes <- events$V2  
  }
  
  # wrong sampling rate in trodes software
  if (dim(events)[1] == dim(eventStrobe)[1]){
    lm0 <- lm( I( 1000*diff(events$V2) ) ~ -1 + I(diff(eventStrobe$eventTimes)) )
    if (lm0$coeff<.8){
      print('correcting wrong sampling rate in spike gadgets')
      events$V2 <- events$V2 * 1/lm0$coeff # for test purposes only
    }
  }
  
  #plot(diff(events$V2[1:1316]), diff(eventStrobe$eventTimes[1:1316]) )
  
  # fix the FFR bug                                                    ====
  # by accident, the start of each trial has an event marker of 14. In FFR_16Tone, 14 is a valid tone marker, hence trial onsets need to be manually removed here:
  # this is causing some confusion: it seems that eventStrobe2==99 should be removed for tone2click, but it is not...
  # why can't it be removed for everybody? Should it be NOT in FFR16Tone or FFR_OA. would make sense but why, FFR_OA?
  
  if (identical(project,'FFR') & dfline$link[1]%in%c('FFR_16Tone_link.txt','FFR_OA_link.txt') )
    eventStrobe <- eventStrobe[-which(eventStrobe$eventStrobe2==99) ,]
  
  if (identical(project,'FFR') & dfline$link[1]%in%c('FFR_tone2click_link.txt','FFR_link_tone.txt') )
    eventStrobe <- eventStrobe[-which(eventStrobe$eventStrobe2==99) ,]
  
  # only for old data sets...
  #if (identical(project,'instrumentalConditioning')&dfline$dataSet%in%c('c','d') )
  #  eventStrobe <- eventStrobe[-which(eventStrobe$eventStrobe4==35) ,]
  
  
  successLogTFileName    <- paste(remote_rda,'/EEG/', project, '/', exp, '/successT.txt', sep = '')
  errorLogTFileName   <- paste(remote_rda,'/EEG/', project, '/', exp, '/errorLogT.txt', sep = '')
  
  pp <- function(){ return(paste(remote_rda,'/EEG/', project, '/', exp, '/',sep='')) }
  successLogTFigureName <- paste( 'successLogTfigure.eps', sep = '')
  errorLogTFigureName   <- paste( 'errorLogTfigure.eps', sep = '')
  
  
  if (dim(eventStrobe)[1] == dim(events)[[1]]){
    
    system( paste('touch ', successLogTFileName, sep = '') )
    
    if (file.exists(errorLogTFileName))
      system( paste('rm ', errorLogTFileName, sep = '') )
  }
  ### check for consistency between timing and eventStrobes file       ====
  if (file.exists(eventStrobeFile) & file.exists(eventTimingFile) & all(!is.na(chStr)) ){
    
    minTr <- min(dim(events)[1], dim(eventStrobe)[1])  
    ##I think we can/should take the tolerance on the initial diff'ed deviation down; maybe 10ms
    ##Never mind; we should be able to, but for some reason the systems can start substantially de-synced
    ##e.g., Jesse TRF 0903 has a 53ms deviation for the first event, falling to very low numbers thereafter 
    
    tmpdiff <- abs(diff(eventStrobe$eventTimes[1:minTr]) - 1000*diff(events$V2[1:minTr]))
    brk <- NA
    
    mxTmpDiff <- 85
    
    # in vocalizations the trigger is not at the beginning of the sound file but somewhere in the middle, this can add to jitter
    if (project %in% c('vocalizations','SMART'))
      mxTmpDiff <- 150
    
    if (identical(project,'SSA')){ # not sure what is going on with SSA, but at the begninning of each block, the isi of ML is way shorter than it actually is. larger discrepancy for larger ISIs...?
      mxTmpDiff <- rep(mxTmpDiff, minTr-1)
      
      blockBreaks <- c( seq(1,250,to=minTr), seq(250,250,to=minTr) )
      #tmpdiffBlockBreaks <- tmpdiff[blockBreaks]
      
      mxTmpDiff[blockBreaks] <- 260
    }
    
    #which(diff(eventStrobe$eventTimes)>3000)
    #plot(diff(events$V2)[1:40], ylim=c(0,2))
    #plot(diff(eventStrobe$eventTimes)[1:40]/1000, ylim=c(0,2))
    
    #plot( abs(diff(eventStrobe$eventTimes[1:minTr]) - 1000*diff(events$V2[1:minTr])), xlim=c(0,500) )
    
    if (any(tmpdiff>mxTmpDiff))
      brk <- min(which( tmpdiff > mxTmpDiff))
    
    if ((dim(eventStrobe)[1] == dim(events)[1])){
      system( paste('rm ', errorLogTFileName, sep = '') )
      system( paste('rm ', pp(), errorLogTFigureName, sep = '') )
      
      system( paste('touch ', successLogTFileName, sep = '') )
      
      fileConn<-file(successLogTFileName)
      writeLines( paste( dim(eventStrobe)[1],'/',dim(events)[1], sep=''), fileConn)
      close(fileConn)
      
      
      eps.file2(successLogTFigureName, s=6, ratio=1/2, pp=pp())
      plot(   eventStrobe$eventTimes[-length(eventStrobe$eventTimes)]/60000,     diff(eventStrobe$eventTimes)/1000 , col='black', type='l', xlab='time[minutes]', ylab='SOA [minutes]')#, xlim=c(0,6000))
      points(  events$V2[-length(events$V2)]/60, diff(events$V2), col='red', type='l' )
      dev.off()
    }
    
       
    if ((dim(eventStrobe)[1] != dim(events)[1])){ 
      #browser()
      #break
      system( paste('touch ', errorLogTFileName, sep = '') )
    
      fileConn<-file(errorLogTFileName)
      writeLines( paste( dim(eventStrobe)[1],'/',dim(events)[1], '/',brk, sep=''), fileConn)
      close(fileConn)
      
      if ((dim(eventStrobe)[1] < dim(events)[1]))
        warning('too many traces')
      
      if ((dim(eventStrobe)[1] > dim(events)[1]))
        warning('missing traces.')
      
      #cor(eventStrobe$eventTimes[1:minTr], events$V2[1:minTr])
      #minTr <- min(dim(events)[1], dim(eventStrobe)[1])  
      
      #plot( diff(eventStrobe$eventTimes[1:minTr]) - 1000*diff(events$V2[1:minTr]), col='black', type='l', ylim=c(-100,100))
      #brk <- min(which( abs(diff(eventStrobe[1:minTr, 1]) - 1000*diff(events$V2[1:minTr])) > mxTmpDiff))
      
      plot(   eventStrobe$eventTimes[-length(eventStrobe$eventTimes)]/60000,     diff(eventStrobe$eventTimes)/1000 , col='black', type='l', xlim=c(0,3))
      points(  events$V2[-length(events$V2)]/60, diff(events$V2), col='red', type='l' )
      abline(v=eventStrobe$eventTimes[brk]/60000, col='green', lty=2)
      
      plot(   diff(eventStrobe$eventTimes)/1000 , col='black', type='l')#, xlim=c(0,300))
      points(  diff(events$V2), col='red', type='l' )
      
      eps.file2(errorLogTFigureName, s=6, ratio=1/2, pp=pp())
      plot(   eventStrobe$eventTimes[-length(eventStrobe$eventTimes)]/60000,     diff(eventStrobe$eventTimes)/1000 , col='black', type='l', xlab='time[minutes]', ylab='SOA [minutes]')#, xlim=c(0,6000))
      points(  events$V2[-length(events$V2)]/60, diff(events$V2), col='red', type='l' )
      abline(v=eventStrobe$eventTimes[brk]/60000, col='green', lty=2)
      dev.off()
      
      
      #plot(        diff(eventStrobe$eventTimes[150:200]) , col='black', type='l', ylim=c(0,4500))
      #points( 1000*diff(events$V2[150:200]), col='red', type='l' )
      
      
      if(brk>1 & dim(events)[1] > dim(eventStrobe)[1]){ 
        #This can also occur (I'm assuming one machine got interrupted/turned off)
        #Check best way to handle
        warning('eventStrobe is missing data at the end.')
        #events <- events[1:minTr,]
        events <- events[1:brk,]
      }
      if(brk>1 & dim(events)[1] < dim(eventStrobe)[1]){ 
        # this is the harmless case when one or two trials are missing at the end. Probably related to either the recording system being turned off
        # prematurely, or the SignalProcessor Toolbox not adequately handling the last minute of the recording  
        warning('RHD seems to have missed the last trials...')
        #eventStrobe <- eventStrobe[1:minTr,]
        eventStrobe <- eventStrobe[1:brk,]
      }
    }
    if (1==0){
      if(is.finite(brk)){
        #browser()
        
        #eventStrobe$V6 =      diff( c(0, eventStrobe$V1));
        #events$V5      = 1000*diff( c(0, events$V2));
        
        #plot(eventStrobe$V6)
        #points(events$V5,col='green')
        
        #        eventStrobe[brk + seq(-5,5),]
        #        events[brk + seq(-5,5),]
        #        diff(eventStrobe$V1[ brk + seq(-5,5) ])
        #        1000*diff(events$V2[ brk + seq(-5,5) ])
        
        #        plot(eventStrobe$V6, xlim=c(brk+c(-20,20)))
        #        points(events$V5,col='green')
        
        #        cleanedEvents <- events[-brk,]
        #        cleanedEvents$V5      = 1000*diff( c(0, cleanedEvents$V2));
        
        #        plot(eventStrobe$V6)#, xlim=c(brk+c(-20,20)))
        #        points(cleanedEvents$V5,col='green')
        
        #pad needs to be big enough to account for the ffr first/last trial in a set offset+all other spurious
        pad<-30
        ##the maximum allowable difference in time between the two signals (in s), before we try to 
        ##correct for a presumably missing or inserted click signal
        max_deviation<-0.075
        diff_deviation<-0.018
        te1<-eventStrobe[,1]/1000-eventStrobe[1,1]/1000
        te2<-events$V2-events$V2[1]
        
        d1<-length(te1)
        d2<-length(te2)
        if(d1>d2){
          te2<-c(te2, rep(0, d1-d2))
        }
        if(d1<d2){
          te1<-c(te1, rep(0, d2-d1))
        }
        rmv<-c()
        
        ##we need an aligned starting point for the vectors, because we want the first event at time = 0 for both
        ##make sure the inital pad length is a multiple of 4
        startCheck<-10
        #      endpad<-round(startCheck/4)
        
        #      tt<-toeplitz(c(rev(diff(te1)[1:startCheck]), rep(0, endpad)))
        #      te1off<-tt[1:(2*endpad),(startCheck+endpad):(startCheck-endpad)]
        #      te2off<-matrix(diff(te2)[1:(2*endpad+1)], nrow = 2*endpad, ncol = 2*endpad+1, byrow = TRUE)
        #      te1off<-matrix(diff(te1)[rep(1:startCheck, startCheck)+rep((1:startCheck)-1, each = startCheck)], nrow = startCheck)
        #      te2off<-matrix(diff(te2)[1:startCheck], nrow = startCheck, ncol = startCheck, byrow = TRUE)
        Dmat<-array(0, dim = c(startCheck, startCheck, startCheck))
        for(i in 1:startCheck){
          for(j in 1:startCheck){
            for(k in 1:startCheck){
              Dmat[i, j, k]<-diff(te1)[i+k-1]-diff(te2)[j+k-1]
            }
          }
        }
        diff_mat<-apply(abs(Dmat)<diff_deviation, c(1, 2), sum)
        startind<-which(diff_mat==max(diff_mat), arr.ind = TRUE)[1,]
        #      diff_mat<-abs(te1off-te2off)<diff_deviation
        #      te1si<-which.max(rowSums(diff_mat))
        te1si<-startind[1]
        if(te1si>1){
          rmv<-rbind(rmv, cbind(1, 1:(te1si-1)))
        }
        #      if(prod(diff_mat[te1si,])){
        #        te2si<-1
        #      }else{
        #        te2si<-max(which(diff_mat[te1si,]=="FALSE"))+1
        #        rmv<-rbind(rmv, cbind(2, 1:(te2si-1)))
        #      }      
        te2si<-startind[2]
        if(te2si>1){
          rmv<-rbind(rmv, cbind(2, 1:(te2si-1)))
        }
        
        if((te1si*te2si)>1){
          te1<-te1-te1[te1si]
          te2<-te2-te2[te2si]
          
          rm1i<-rmv[rmv[,1]==1, 2]
          rm2i<-rmv[rmv[,1]==2, 2]
          remain1i<-setdiff(1:(length(te1)-length(rm2i)), rm1i)
          remain2i<-setdiff(1:(length(te2)-length(rm1i)), rm2i)
          brkall<-which(abs(te2[remain2i]-te1[remain1i])>max_deviation)
          if(length(brkall)>0){
            brki1<-min(brkall)+length(rm1i)
            brki2<-min(brkall)+length(rm2i)
          }
        }else{      
          brkall<-which(abs(te2-te1)>max_deviation)
          brki1<-brki2<-min(brkall)
        }
        te1<-c(te1, rep(0, pad))
        te2<-c(te2, rep(0, pad))
        
        seqerr<-1
        kick<-1
        ##this loop checks if removing one event time brings the sequences back within tolerance
        ##if not, we add one to our putative error sequence length (seqerr), and look again
        ##
        ##Currently we may remove from either source, to correspond to a potential missing signal from either, 
        ##imagining that monkeylogic sent a signal, but for some reason no tone was played, we would remove 
        ##the monkeylogic event
        ##
        ##However, we never impute a signal; this assumes that removing something is sufficient to fix a mismatch
        
        rm1i<-rmv[rmv[,1]==1, 2]
        rm2i<-rmv[rmv[,1]==2, 2]
        remain1i<-setdiff(1:(length(te1)-length(rm2i)), rm1i)
        remain2i<-setdiff(1:(length(te2)-length(rm1i)), rm2i)
        
        if(length(brkall)>0){
          kick<-0
          diff1<-abs(te2[brki2]-te1[brki1+seqerr])
          diff2<-abs(te2[brki2+seqerr]-te1[brki1])
        }
        #      browser()
        
        while(kick==0){
          #removing diff1 means we can no longer assume a spurious ML event, 
          #that is, we can now only remove from the Intan events 
          #  diff1<-abs(te2[brki]-te1[brki+seqerr])
          #  diff2<-abs(te2[brki+seqerr]-te1[brki])
          
          ##This first #if is for FFR; single mismatches that return to within tolerance after the point 
          ##are "left alone", these are the final tone in a sequence translated over to the next sequence
          ##(visible in matlab ML data), though they are definitely messing up downstream analysis
          ##Will have to figure this out later
          
          if(((brki1+seqerr)<=d1) & ((brki2+seqerr)<=d2) & (abs(te2[brki2+seqerr]-te1[brki1+seqerr])<max_deviation)){
            #  if(((brki+seqerr)<=min(d1, d2)) & (abs(te2[remain2i][brki+seqerr]-te1[remain1i][brki+seqerr])<max_deviation)){
            rmv<-rbind(rmv, cbind(3, seq(from = brki1, to = brki1+seqerr-1)), cbind(4, seq(from = brki2, to = brki2+seqerr-1)))
            kick<-1
          }else if(diff1<max_deviation){
            rmv<-rbind(rmv, cbind(1, seq(from = brki1, to = brki1+seqerr-1)))
            kick<-1
          }else if(diff2<max_deviation){
            rmv<-rbind(rmv, cbind(2, seq(from = brki2, to = brki2+seqerr-1)))
            kick<-1
          }else{
            seqerr<-seqerr+1
            #        diff1<-abs(te2[remain2i][brki]-te1[remain1i][brki+seqerr])
            #        diff2<-abs(te2[remain2i][brki+seqerr]-te1[remain1i][brki])
            diff1<-abs(te2[brki2]-te1[brki1+seqerr])
            diff2<-abs(te2[brki2+seqerr]-te1[brki1])
          }
          if(kick==1){
            rm1i<-rmv[(rmv[,1]==1 | rmv[,1]==3), 2]
            rm2i<-rmv[(rmv[,1]==2 | rmv[,1]==4), 2]
            remain1i<-setdiff(1:(length(te1)-length(rm2i)), rm1i)
            remain2i<-setdiff(1:(length(te2)-length(rm1i)), rm2i)
            #?
            #          remain1i<-setdiff(1:length(te1), rm1i)
            #          remain2i<-setdiff(1:length(te2), rm2i)
            brkall<-which(abs(te2[remain2i]-te1[remain1i])>max_deviation)
            if(length(brkall)>0){
              brki1<-min(brkall)+length(rm1i)
              brki2<-min(brkall)+length(rm2i)
              if(brki1<(d1+2) & brki2<(d2+2)){
                seqerr<-1
                kick<-0
                #      diff1<-abs(te2[remain2i][brki]-te1[remain1i][brki+seqerr])
                #      diff2<-abs(te2[remain2i][brki+seqerr]-te1[remain1i][brki])
                diff1<-abs(te2[brki2]-te1[brki1+seqerr])
                diff2<-abs(te2[brki2+seqerr]-te1[brki1])
              }
            }
          }
          if(seqerr>10){
            print('Event sequences grossly misaligned; data needs manual review.')
            break
          }
        }
        
        #  browser()
        
        overflow<-which(rmv[,2]>c(d1, d2)[1+(rmv[,1]%%2)])
        if(length(overflow)>0){
          rmv<-rmv[-overflow,]
        }
        
        if(sum(rmv[,1]==1)>0){
          eSrmvI<-rmv[which(rmv[,1]==1), 2]
          numBegineSrmv<-max(which(cumsum((1:length(eSrmvI))-eSrmvI)==0))
          if(min(eSrmvI)==1){
            if(numBegineSrmv==sum(rmv[,1]==1)){
              warning(paste("We are removing", numBegineSrmv, "event(s) from the beginning of eventStrobe, 
              appearing as spurious (but really corresponding to missing elements of events), likely 
              because of a slight delay in the start of Intan's recording"))
            }
          }   
          if((min(eSrmvI)!=1) | (numBegineSrmv!=sum(rmv[,1]==1))){
            warning("There are extra events in eventStrobe, not solely at the beginning of the session, 
                    which should be removed to ideally align.  Consider reviewing manually.")
          }
          eventStrobe<-eventStrobe[-rmv[which(rmv[,1]==1), 2],]
        }
        if(sum(rmv[,1]==2)>0){
          events<-events[-rmv[which(rmv[,1]==2), 2],]
          events$V5<-1000*diff(c(0, events$V2));
        }
      }
    }
  }
  
  #browser()
  ##the events.txt and eventStrobe.txt files didn't always have column names
  ##!!!!!let's hope they always kept the same ordering/identity
  if (file.exists(eventStrobeFile)){
    #tone$xpp          <- data.frame( e.MLonset=eventStrobe[,1], strobe1=eventStrobe[,2], strobe2=eventStrobe[,3], strobe3=eventStrobe[,4], trialNr=eventStrobe[,5])
    
    if (!exists('eventStrobe4', eventStrobe))
      eventStrobe$eventStrobe4 <- c(1)
    if (!exists('eventStrobe5', eventStrobe))
      eventStrobe$eventStrobe5 <- c(1)
    
    
    tone$xpp          <- data.frame( e.MLonset=eventStrobe$eventTimes, strobe1=eventStrobe$eventStrobe1, strobe2=eventStrobe$eventStrobe2, strobe3=eventStrobe$eventStrobe3, strobe4=eventStrobe$eventStrobe4, strobe5=eventStrobe$eventStrobe5, trialNr=eventStrobe$trialNr)
    
    
    if (file.exists(eventTimingFile) & dim(tone$xpp)[1]==dim(events)[1] ){
      tone$xpp$e.onset=events[,2]
    } else{
      tone$xpp$e.onset=.001*eventStrobe[,1]
    }
  }
  
  #if (!file.exists(eventStrobeFile)){
  #  tone$xpp          <- data.frame( e.onset=events[,2] )
  #}
  
  #  tone$xpp          <- data.frame( e.MLonset=eventStrobe$eventTimes, e.onset=events$V2, strobe1=eventStrobe$eventStrobe1, strobe2=eventStrobe$eventStrobe2, strobe3=eventStrobe$eventStrobe3, trialNr=eventStrobe$trialNr  )
  #tone$xpp          <- data.frame( e.MLonset=eventStrobe$V1, e.onset=events$V2, strobe1=eventStrobe$V2, strobe2=eventStrobe$V3, strobe3=eventStrobe$V4, trialNr=eventStrobe$V5  )
  
  #browser()
  
  tone$trialInfo <- NULL
  if (file.exists(trialStrobeFile))
    tone$trialInfo          <- trialStrobe
  
  if (!file.exists(eventStrobeFile)){# if ml file is not readable, create random strobe values to keep the program from crashing
    warning('simulating factors for strobe; do not use this data for real analyses')
    tone$xpp          <- data.frame( e.onset=events[,2], strobe1=events[,1], strobe2=rbinom(dim(events)[1],1,.5), strobe3=rbinom(dim(events)[1],1,.5), strobe4=rbinom(dim(events)[1],1,.5) )
    #tone$xpp          <- data.frame( e.onset=events$V2, strobe1=events$V1)
  }
  tone$xpp$toneNr   <- 1:dim(tone$xpp)[1]
  tone$project      <- project
  tone$dfline       <- dfline
  
  # from here on, assume that MLtrials and the recording system have identified the same events
  
  ##  read time-continuous data if requested ====
  if(nchar(extn)>0)
    print(extn)
  
  chextn <- paste('_',epoch,extn,sep='')
  print(chextn)
  
  #if (identical(epoch,'llp'))
  #  chextn <- extn
  
  tone$taxis <- NULL
  tone$lfp   <- NULL
  
  #  browser()
  if (!epoch=='dmu'){ # make sure we are not requesting digital mu
    if ( all(!is.na(chStr)) ){ 
      
      if (!exists('numChan',dfline))
        dfline$numChan <- length(chStr)
      if (!exists('startVIndex',dfline))
        dfline$startVIndex <- 1
      if (!exists('probeString',dfline))
        dfline$probeString <- substr(dfline$ElectrodeConfig,1,1)
      
      if (is.na(dfline$numChan))
        dfline$numChan <- length(chStr)
      if (is.na(dfline$startVIndex))
        dfline$startVIndex <- 1
      
      chDepthInd <- 1:dfline$numChan
      csdUseInd  <- rep(T, dfline$numChan)
      
      if (dfline$probeString=='d'){
        chDepthInd <- c(1:12,rep(13:20,each=2), 21:24)
        csdUseInd[seq(14,by=2,length=8)] <- F
      }
      chNum           <- unlist(lapply(1:length(chStr),function(n){ as.numeric(strsplit(chStr[n],'h')[[1]][2]) })) - dfline$startVIndex + 1
      #    tone$chDepthInd <- chDepthInd[chNum]
      tone$chDepthInd <- c(chNum[chNum<1], chDepthInd[chNum>0])
      if(sum(chNum<1)>0 & sum(chNum>0)>0){
        warning("Element chDepthInd mixes valid and invalid depths; don't use this structure for CSD.")
      }
      
      tone$csdUseInd       <- csdUseInd
      tone$channelPosition <- data.frame(xpos=1+(0:(length(chStr)-1))%/%chBlkN, ypos=1+(0:(length(chStr)-1))%%chBlkN,  zpos=rep(0,length(chStr)) )
      
      
      # read time-axis ====
      signalDirExtn <- ''
      taxFile <-             paste(rda,'/', signalDirExtn,exp,'_taxis',chextn,'.mat',sep='')
      if (!file.exists(taxFile)){# File may not exist because it is in epoch-specific subfolder
        nynew <- T  
        signalDirExtn <- paste(epoch,'/',sep='')
      }
      taxFile <-             paste(rda,'/', signalDirExtn,exp,'_taxis',chextn,'.mat',sep='')
      print(                 taxFile )
      
      
      #tone$taxis <- readMat( paste(rda,'/', exp,'_taxis',chextn,'.mat',sep=''))[[1]]
      tone$taxis <- readMat( taxFile )[[1]]
      tone$taxis <- c(tone$taxis)
      
      tone$dt       <- mean(diff(tone$taxis)) # temporal resolution in ms
      tone$SR       <- 1000/tone$dt
      
      # create log-timeaxis ====
      toff                                     <- 1
      suppressWarnings(tone$logtaxis           <- log2( tone$taxis + toff))
      tone$logtaxis[ (tone$taxis + toff <= 1)] <- .01*(tone$taxis[ which(tone$taxis + toff <= 1) ] - 1 + toff) 
      
      # read channels ====
      
      print(                  paste(rda,'/',signalDirExtn,exp,'_',chStr[1],chextn,'.mat',sep=''))
      tmpLFP   <- t( readMat( paste(rda,'/',signalDirExtn,exp,'_',chStr[1],chextn,'.mat',sep=''))[[1]] )
      tone$lfp <- array(NA, c(dim(tmpLFP)[1],dim(tmpLFP)[2], length(chStr)) )
      tone$lfp[,,1] <- tmpLFP 
      
      if(length(chStr)>1){
        for (chnd in 2:length(chStr)){
          tone$lfp[,,chnd]   <- t( readMat( paste(rda,'/',signalDirExtn, exp,'_',chStr[chnd],chextn,'.mat',sep=''))[[1]] )
          print(                            paste(rda,'/',signalDirExtn, exp,'_',chStr[chnd],chextn,'.mat',sep=''))
        }
      }
      gc()
      if(exists('rmv')){
        if(sum(rmv[,1]==2)>0){
          if (length(dim(tone$lfp))==3){
            tone$lfp <- tone$lfp[-rmv[rmv[,1]==2,2],,]
          }
          if (length(dim(tone$lfp))==2){
            tone$lfp <- tone$lfp[-rmv[rmv[,1]==2,2],]
          }
        }
      }
      
      # Not sure what this is for...
      if ('trialInfo' %in% names(tone)){
        if ( dim(tone$lfp)[1] < dim(tone$xpp)[1] & dim(tone$lfp)[1]==dim(tone$trialInfo)[1] ){
          print(' swapping xpp and trialInfo')
          xppBackup    <- tone$xpp
          tone$xpp     <- tone$trialInfo
          tone$snp1    <- list()
          tone$snp2    <- list()
          tone$snp3    <- list()
          tone$snp4    <- list()
          tone$snpTime <- list()
          
          for (tx in 1:dim(tone$xpp)[1]){
            tone$snp1[[tx]]    <- xppBackup$strobe1[   which(xppBackup$trialNr==tx) ]
            tone$snp2[[tx]]    <- xppBackup$strobe2[   which(xppBackup$trialNr==tx) ]
            tone$snp3[[tx]]    <- xppBackup$strobe3[   which(xppBackup$trialNr==tx) ]
            tone$snp4[[tx]]    <- xppBackup$strobe4[   which(xppBackup$trialNr==tx) ]
            tone$snpTime[[tx]] <- xppBackup$e.MLonset[ which(xppBackup$trialNr==tx) ]
          }
        }
      }
      
      #tst <- compareEventTimes( .001*tone$xpp$e.MLonset, tone$xpp$e.onset)
      
      # removing trials if xpp and lfp don't match ====
      if (dim(tone$lfp)[1]> dim(tone$xpp)[1] ){
        warning('too many trials in tone$lfp. Removing trials.  This should be manually reviewed!')
        
        ##Error in events$V3 : $ operator is invalid for atomic vectors
        ##likely means that eventTimingFile doesn't exist
        #table(events$V3)
        
        if (length(dim(tone$lfp))==3)
          tone$lfp <- tone$lfp[1:dim(tone$xpp)[1],,]
        if (length(dim(tone$lfp))==2)
          tone$lfp <- tone$lfp[1:dim(tone$xpp)[1],]
      }
      
      if ( dim(tone$lfp)[1] < dim(tone$xpp)[1]  ){
        warning('missing traces. Removing lines from xpp. You should probably re-run signal processor to get all trials')
        tone$xpp <- tone$xpp[1:dim(tone$lfp)[1],]
      }
      
      # baseline-correction                        ====
      lfpnbl <- tone$lfp # store non-baseline correct lfp for later use
      if (length(blrng)==2){
        print('baseline correction')
        tmp            <- baseLineCorrection(tone$lfp,blrng, tone$taxis)  
        tone$lfp       <- tmp$res
        tone$xpp$blcor <- tmp$correction
        rm(tmp); gc()
      }
      
      if(length(dim(tone$lfp))==2){
        tone$lfp<-array(tone$lfp, dim = c(dim(tone$lfp), 1))
      }
      
      dimnames(tone$lfp) <- list(trial=1:dim(tone$lfp)[1],time=tone$taxis,channel=chStr)
    }
  }
  
  ## read digital multi unit data. Unsorted threshold crossings ====
  if (epoch=='dmu'){
    if ( all(!is.na(chStr)) ){ # get lfps only if chStr is specified
      
      dmu <- list()
      for (chx in 1:length(chStr)){
        dmuFileName <-             paste(rda, 'spk/',chStr[chx],'/time_m30.txt',sep='')
        spktime <- read.table(dmuFileName, header=T) 
        SR      <- as.numeric(strsplit(names(spktime),'X' )[[1]][2])
        spktimeM <- spktime$X30000 / SR
        
        tdmu <- list()
        for (tnd in 1:dim(events)[1]){
          relSpkTimeM <- spktimeM - events$V2[tnd]
          tdmu[[tnd]] <- relSpkTimeM[ which(relSpkTimeM>suarng[1] & relSpkTimeM<suarng[2]) ]
        }
        dmu[[chx]] <- tdmu
      }# loop chStr
    }
    tone$dmu <- dmu
  }
  
  
  ## read spike data if requested            ====
  if (readSUA){
    ##will need a new section for the new format!
    
    print('reading threshold crossings')
    if (suaMethod=='nev'){
      ##This is all broken, spd() doesn't exist as a fxn, _nev.txt files don't exist in the directory, 
      #the txt files have other rows with time since previous and next stim. ONLY, no raw spike times 
      #updating to actually read the old data, e.g., RSk Sam_20160425, now returning channel MUA+hash too
      #SK 20220408
      
      # read in the nev spike timing file for the unit in question 
      #      print( paste( spd, dfline$session,'/ch_', dfline$channel,'/time_m30_nev.txt',sep=''))
      #      spktimeN  <- read.table(  paste( spd, dfline$session,'/ch_', dfline$channel,'/time_m30_nev.txt',sep=''))$V1
      #      spktimeM  <- spktimeN
      
      ##This will be old nev reading (in terms of soa information) ====
      if(project == "RSkernel"){
        print(paste(spd, dfline$session, '/ch_', dfline$channel, '/soa_', dfline$cell, '.txt', sep=''))
        
        spktimeN  <- readMat(paste( spd, dfline$session,'/ch_', dfline$channel, '/soa_', dfline$cell, '.mat', sep = ''))$soa
        
        chls<-list.files(paste( spd, dfline$session,'/ch_', dfline$channel, sep = ''))
        mil<-which(unlist(lapply(chls, function(x){gregexpr(paste(dfline$thrsh, '.mat', sep = ""), x)[[1]][1]!=-1})))
        
        if(length(mil)>1){
          spktimeMt  <- rbind(readMat(paste(spd, dfline$session, '/ch_', dfline$channel, '/', chls[mil[1]], sep = ''))$soa, 1)
          for(i in 2:length(mil)){
            spktimeMt  <- cbind(spktimeMt, rbind(readMat(paste(spd, dfline$session, '/ch_', dfline$channel, '/', chls[mil[i]], sep = ''))$soa, i))
          }
          tmp<-spktimeMt[,order(spktimeMt[1,])]
          spktimeM<-tmp[,order(tmp[3,])]
        }else{
          spktimeM  <- spktimeN
        }
        spktimeH  <- read.table(paste( spd, dfline$session,'/ch_', dfline$channel, '/time_', dfline$thrsh, '.txt', sep = ''))
        #        spktimeN  <- readMat(  paste( spd, dfline$session,'/ch_', dfline$channel,'/soa_n2_m30.mat',sep=''))$soa
        
        absspktimeN<-c(rep(events$V2[1], sum(spktimeN[3,]==0)), events$V2[spktimeN[3,]]) + c(spktimeN[2, which(spktimeN[3,]==0)], spktimeN[1, which(spktimeN[3,]>0)])
        absspktimeM<-c(rep(events$V2[1], sum(spktimeM[3,]==0)), events$V2[spktimeM[3,]]) + c(spktimeM[2, which(spktimeM[3,]==0)], spktimeM[1, which(spktimeM[3,]>0)])
        absspktimeH<-spktimeH[,1]
        
      }
      
      ##new reading process, keeping directories and whatever else possible, the same ====
      if(project == "ampOddClick"){
        
        #        browser()
        
        ##fp0 is for raw conversion to v6 format by the bottom of the standalone NEV2mat function
        fp0<-paste( spd, dfline$session,'/ch_', dfline$channel, '/sortS1C.mat', sep = '')
        fp1<-paste( spd, dfline$session,'/ch_', dfline$channel, '/sortS1t.mat', sep = '')
        fp2<-paste( spd, dfline$session,'/ch_', dfline$channel, '/sortS1.mat', sep = '')
        if(file.exists(fp0)){
          spktimeMt  <- readMat(fp0)$NEV[,,1]$NEV[,,1]$Data[,,1]$Spikes[,,1]
          print(paste0('Reading: ', fp0))
        }else if(file.exists(fp1)){
          spktimeMt  <- readMat(fp1)$NEV[,,1]$Data[,,1]$Spikes[,,1]
          print(paste0('Reading: ', fp1))
        }else{
          spktimeMt  <- readMat(fp2)$NEV[,,1]$Data[,,1]$Spikes[,,1]
          print(paste0('Reading: ', fp2))
        }
        
        findis<-table(spktimeMt$Unit)
        
        hashi<-as.numeric(names(which.max(findis)))
        muai<-setdiff(1:max(spktimeMt$Unit), hashi)
        suai<-muai[dfline$cell]
        
        #tone$SR is wrong is piggybacking on other epochs, find where that would be kept
        #in general, could it be != 30000?
        absspktimeN<-spktimeMt$TimeStamp[1, which(spktimeMt$Unit==suai)]/30000
        absspktimeM<-spktimeMt$TimeStamp[1, which(spktimeMt$Unit%in%muai)]/30000
        absspktimeH<-spktimeMt$TimeStamp[1,]/30000
        
      }
      
      #shapeFile <- paste( spd(), dfline$session,'/ch_', dfline$channel,'/',dfline$cell,'_nev_',dfline$fextn,'_shape.mat',sep='')
      #if (file.exists(shapeFile))
      #  shp       <- readMat(shapeFile)$shp
    }# end nev method
    
    ##no idea how the ws data would be structured, but safe bet this would now crash when 
    #calculating relative times
    if (suaMethod=='ws'){
      # read in the timing file and the wavesorter selection file 
      spktime   <- read.table(  paste( spd, dfline$session,'/ch_', dfline$channel,'/time.txt',sep='') )
      tmp       <- read.table(  paste( spd, dfline$session,'/ch_', dfline$channel,'/', dfline$cell,'.txt',sep='') )
      cellInd   <- tmp[-1,1] 
      if (length(cellInd) != length(spktime$V1)) 
        browser()
      
      spktimeM   <- spktime[                   ,1]
      spktimeN   <- spktime[ which(cellInd==1) ,1]
      
      spk <- list()
      mua <- list()
      for (tnd in 1:dim(events)[1]){
        relSpkTimeN <- spktimeN - events$V2[tnd]
        spk[[tnd]] <- relSpkTimeN[ which(relSpkTimeN>suarng[1] & relSpkTimeN<suarng[2]) ]
        
        relSpkTimeM <- spktimeM - events$V2[tnd]
        mua[[tnd]] <- relSpkTimeM[ which(relSpkTimeM>suarng[1] & relSpkTimeM<suarng[2]) ]
      }
    }# end ws method
    
    if(dim(events)[1]<dim(eventStrobe)[1]){
      eT<-eventStrobe$eventTimes/1000+events$V2[1]-eventStrobe$eventTimes[1]/1000
    }else{
      eT<-events$V2
    }
    
    spk  <- list()
    mca  <- list()
    hash <- list()
    for (tnd in 1:length(eT)){
      #find the list of eventtimes
      #all spikes within +-t of each 
      
      #      relSpkTimeN <- spktimeN - events$V2[tnd]
      #      spk[[tnd]] <- relSpkTimeN[ which(relSpkTimeN>suarng[1] & relSpkTimeN<suarng[2]) ]
      
      relSpkTimeN <- absspktimeN - eT[tnd]
      spk[[tnd]] <- relSpkTimeN[which(relSpkTimeN>suarng[1] & relSpkTimeN<suarng[2])]
      
      #      relSpkTimeM <- spktimeM - events$V2[tnd]
      #      mua[[tnd]] <- relSpkTimeM[ which(relSpkTimeM>suarng[1] & relSpkTimeM<suarng[2]) ]
      relSpkTimeM <- absspktimeM - eT[tnd]
      mca[[tnd]] <- relSpkTimeM[which(relSpkTimeM>suarng[1] & relSpkTimeM<suarng[2])]
      
      relSpkTimeH <- absspktimeH - eT[tnd]
      hash[[tnd]] <- relSpkTimeH[which(relSpkTimeH>suarng[1] & relSpkTimeH<suarng[2])]
    }
    tone$spk <- spk
    tone$mca <- mca
    tone$hash <- hash
    
  }# end readSUA
  
  # subset range of valid tones if specified ======================
  #trialRng <- c(1,length(tone$xpp$trialNr))
  #if(!is.na(from))
  #  trialRng[1] <- from
  
  #if(!is.na(to))
  #  trialRng[2] <- to
  
  #tone <- subs.data(tone, 1:length(tone$xpp$trialNr) %in% trialRng[1]:trialRng[2])
  
  # make sure lfp has three dimensions       ====
  if(!is.null(tone$lfp)){
    if(length(dim(tone$lfp))==2)
      dim(tone$lfp) <- c(dim(tone$lfp), 1)
  }
  gc()
  
  # prepXPP                                  ====  
  tone              <- prepXPP(tone)
  # return                                   ====
  return(tone) 
}


chan2varSignalProcessor <- function(ttn){
  
  tst     <- ttn
  tst$xpp <- ttn$xpp[ rep(1:dim(ttn$xpp)[1], each=dim(ttn$lfp)[3]) ,]
  lfp <- aperm(tst$lfp,c(3,1,2)) 
  #lfp <- aperm(tst$lfp,c(1,3,2)) 
  
  dim(lfp) <- c(dim(tst$lfp)[1]*dim(tst$lfp)[3], dim(tst$lfp)[2],1 )
  tst$lfp <- lfp
  
  tst$channelPosition <- data.frame(xpos=1,ypos=1,zpos=1)
  tst$xpp$chStr <- rep(tst$chStr)
  tst$xpp$chID <- rep(1:length(chStr))
  tst$chStr <- 'mix'
  
  return(tst)
}

nrmByChanSignalProcessor <- function(tone, FUN=function(x){max(abs(x))}, frml=~chID){
  
  tndx <- 1:length(tone$taxis)
  if (!is.null(trng))
    tndx <- which( tone$taxis >= trng[1] & tone$taxis <= trng[2] )
  
  mn         <- avgSignalProcessor( tone, frml=frml )
  
  if (dim(mn$lfp)[3]>1)
    mn$xpp$nrm <- apply( mn$lfp[,tndx,], c(1,3), FUN )
  if (dim(mn$lfp)[3]==1)
    mn$xpp$nrm <- apply( mn$lfp[,tndx,], c(1), FUN )
  
  frmStr <- deparse(frml[[length(frml)]])
  tst         <- tone
  tst$xpp$nrm <- NA
  for (cx in 1:dim(tone$lfp)[1]){
    tst$xpp$nrm[cx] <- mn$xpp$nrm[ tst$xpp[[frmStr]][cx] ]
    tst$lfp[cx,,]   <- tone$lfp[cx,,] / mn$xpp$nrm[ tst$xpp[[frmStr]][cx] ]
  }
  return(tst)
}

# this is the default version that removes mean and normalizes sd to 1, or both
nrmSignalProcessor <- function(tone, demean=T, devar=T, trng=NULL){
  
  tndx <- 1:length(tone$taxis)
  if ( !is.null(trng) )
    tndx <- which(tone$taxis >= trng[1] & tone$taxis <= trng[2])
  
  # remove mean for all channels (jointly across all trials) 
  tLFP    <- three2twoDimSignalProcessor( tone$lfp)
  if (demean){
    mnLFP   <- array( rep( apply(tLFP,2,mean), each=dim(tLFP)[1]), dim(tLFP)) 
    tLFP    <- tLFP - mnLFP
  }
  
  if (devar){
    vrLFP   <- array( rep( apply(tLFP,2,var), each=dim(tLFP)[1]), dim(tLFP)) 
    tLFP  <- tLFP / sqrt(vrLFP)
  }
  tone$lfp <- two2threeDimSignalProcessor(tLFP, dm=dim(tone$lfp))
  return(tone)
}
nrmByTrialSignalProcessor <- function(tone, demean=T, devar=T){
  
  # remove mean for all channels (jointly across all trials) 
  tLFP    <- three2twoDimSignalProcessor( tone$lfp)
  if (demean){
    mnLFP   <- array( rep( apply(tLFP,2,mean), each=dim(tLFP)[1]), dim(tLFP)) 
    tLFP    <- tLFP - mnLFP
  }
  
  if (devar){
    vrLFP   <- array( rep( apply(tLFP,2,var), each=dim(tLFP)[1]), dim(tLFP)) 
    tLFP  <- tLFP / sqrt(vrLFP)
  }
  tone$lfp <- two2threeDimSignalProcessor(tLFP, dm=dim(tone$lfp))
  return(tone)
}


checkISISignalProcessor <- function(tn, show=F, mxJitter=28){
  
  # e.MLonset provides timing from the MonkeyLogic Computer
  # e.onset provides timing for the same events from the Intan recording
  # they should be identical other than an offset, a scaling factor of 1000 and a jitter of approximately 1/frame rate (18 ms)
  # 
  # approach: calculate time difference between adjacent events for both measures of time
  #           fit a linar model, and check model fit. If max error is small, the data set is good
  
  lm1     <- lm( diff(tn$xpp$e.MLonset) ~ diff(tn$xpp$e.onset) )
  
  if (show)
    plot(diff(tn$xpp$e.MLonset), diff(tn$xpp$e.onset) )
  
  
  # compare intervals   
  #somethingWrong <-   max(abs(lm1$residuals)) > mxJitter
  somethingWrong <-   max(abs(lm1$residuals[-1])) > mxJitter | abs(lm1$residuals[1])>80
  
  
  if (abs(lm1$residuals[1])>mxJitter)
    warning('first Trigger seems to be off by more than mxJitter, but below 80 ms. Ignoring the warning because first trial timing can be more variable than other triggers')
  
  if (somethingWrong){
    
    badNdx <- which(abs(lm1$residuals) > mxJitter)
    
    N1Error <- F
    # CASE 1: simple case of one wrong e.onset
    if ( length( badNdx )  == 2){
      if (diff(badNdx) == 1  )   {
        
        N1Error <- T
        errorNdx <- min( badNdx )
        
        # step 1: estimate correct e.onset for the wrong e.onset
        lm2 <- lm(tn$xpp$e.onset ~ tn$xpp$e.MLonset )
        tn$xpp$e.onset[ badNdx ] <- lm2$fitted.values[ badNdx ]
        
        
        # Step 2: flag the error trial with a p2p of Inf
        tn$xpp$p2p[errorNdx]     <- Inf
        tn$xpp$valTone[errorNdx] <- F
        
        warning('One wrong e.onset was found and corrected. ')
        
      }
    }
    
    # CASE 2: everything else. Not handled yet
    if ( !N1Error  )   {
      # unrecoverable error
      stop(' timing between intan and ML is off. Most likely caused by missing or false triggers.')
    }
  }
  
  return(tn)
}


prepLFP <- function(tone, trng=NULL){
  
  if (exists('taxis',tone)){
    if (abs(mean(diff(tone$taxis))-1)<10e-10){
      tone$taxis <- round(tone$taxis)
      warning('setting sampling rate to 1ms')
    }
  }
  
  if (exists('lfp',tone)){
    
    tone$Ntrial  <- dim(tone$lfp)[1]
    tone$Ntst    <- dim(tone$lfp)[2]
    tone$Nchan   <- dim(tone$lfp)[3]
    
    tndx <- 1:tone$Ntst
    if (!is.null(trng) )
      tndx <- which(tone$taxis >= trng[1] & tone$taxis<=trng[2])
    
    tone$p2pByChan    <- apply(tone$lfp[,tndx,],c(1,3),max) - apply(tone$lfp[,tndx,],c(1,3),min)
    tone$xpp$p2p      <- apply( tone$p2pByChan, 1, median )
    
    qtmp               <- quantile(tone$xpp$p2p,probs=c(0.5,0.75))
    mxp2p              <- qtmp[1]+5*diff(qtmp)
    tone$xpp$valTone   <- tone$xpp$p2p < mxp2p;
    
    if ( tone$Nchan > 10 ){
      
      # analyze over signal quality of individual channels
      tone$channelQmed <- apply(tone$p2pByChan,2,median)
      tone$channelQmin <- apply(tone$p2pByChan,2,min)
      
      if (1==0){
        plot( tone$channelQmed, type='l' )
        points( tone$channelQmin, type='l', col='red' )
      }
      
      ## reject an entire channel if the minum or medial p2p is well above the range of other channels
      #rngQmed <- quantile(tone$channelQmed, probs=c(.25,.75) )
      #mxQmed  <- rngQmed[1] + 3*diff(rngQmed)
      #which(tone$channelQmed>mxQmed)
      
      #rngQmin <- quantile(tone$channelQmin, probs=c(.25,.75) )
      #mxQmin  <- rngQmin[1] + 3*diff(rngQmin)
      #which(tone$channelQmin>mxQmin)
      
      # run rejection criteria for each channel separately
      getValToneChannel <- function(cx){
        qtmp               <- quantile(tone$p2pByChan[,cx],probs=c(0.25,0.75))
        mxp2p              <- qtmp[1]+5*diff(qtmp)
        valTone            <- tone$p2pByChan[,cx] <= mxp2p
        return(valTone)
      }
      
      tone$valp2pByChan <- sapply(1:tone$Nchan, getValToneChannel)
      chQ               <- apply(tone$valp2pByChan,2,mean)
      
      if (1==0){
        plot(chQ, ylab='Proportion of accepted trials') 
        plot(sort(chQ), ylab='Proportion of accepted trials')
      }
      
      
      criterion        <- 0.90                           # at least 90 percent of trials pass criterion for a channel
      flexcriterion    <- quantile( chQ, probs=c(.1) )   # reject 10 percent of channels at most
      tone$cleanChan   <- chQ >  min(criterion,flexcriterion)
      
      avgVal  <- apply(tone$valp2pByChan[,which(tone$cleanChan)],1,mean)
      
      if(1==0)
        plot(sort(avgVal))
      
      tone$xpp$valTone  <- avgVal == 1 
      if (dim(tone$lfp)[3]>400)
        tone$xpp$valTone  <- avgVal >= (dim(tone$lfp)[3]-2)/dim(tone$lfp)[3] # allow two bad channels for each trial 
      
    }
  }
  gc()
  return(tone)
}


prepXPP <- function(tone){
  
  tone$xpp$p2p         <- 0
  tone$xpp$randomSplit <- rbinom(dim(tone$xpp)[1],1,.5)
  
  if (exists('lfp',tone)){
    
    trng <- NULL
    
    if (tone$project %in% c('rewardAnticipation'))
      trng <- c(-250,1800)
    
    tone <- prepLFP(tone, trng=trng)
    #tmp           <- apply(tone$lfp,c(1,3),max) - apply(tone$lfp,c(1,3),min)
    #tone$channelQ <- apply(tmp,2,mean)
    #tone$xpp$p2p  <- apply( tmp, 1, median )
    
    #qtmp               <- quantile(tone$xpp$p2p,probs=c(0.5,0.75))
    #mxp2p              <- qtmp[1]+5*diff(qtmp)
    #tone$xpp$valTone   <- tone$xpp$p2p < mxp2p;
  }
  
  #tone$xpp$p2p                 <- apply(tone$lfp,1, max) - apply(tone$lfp,1, min)
  tone$xpp$animal <- tone$animal
  
  ##this may be removed in the future depending on the format of link files
  if(!exists('link', tone$dfline)){
    tone$dfline$link<-tone$project
  }

  ## STRF                 ====
  if ( tone$project %in% c('STRF')){
    
    tone$xpp$freqInd <- tone$xpp$strobe2
    # ------------ frequency
    freqList          <- tone$dfline$minOct*exp(log(2)*seq(0,by=tone$dfline$stepOct, length=tone$dfline$NFrq))
    tone$xpp$freqHz   <- freqList[ tone$xpp$freqInd ]
    tone$xpp$freqOct  <- log2( freqList[ tone$xpp$freqInd ])
  }
  
  ## TRF                  ====
  if ( tone$project %in% c('TRF','activeTRF','amRF','noiseTRF','narrowBandNoise','clickRF','reactivationTRF') & !tone$dfline$link=='DTRF'){
    
    tone <- checkISISignalProcessor( tone )
    
    tone$xpp$freqInd <- tone$xpp$strobe2
    
    
    # old version
    #tone$xpp$ampInd  <- tone$xpp$strobe3
    
    if (tone$project=='reactivationTRF'){
      tone$xpp$typeInd <- tone$xpp$strobe1
      tndx             <- which(tone$xpp$typeInd==4)
      wndx             <- which(tone$xpp$typeInd %in% c(3,5,6))
      
      if ( !sum(tndx+1==wndx)==length(tndx) )
        stop('tones and masks not alternating')
      
      tone$xpp$freqInd        <- tone$xpp$strobe2
      tone$xpp$freqInd[wndx]  <- tone$xpp$freqInd[tndx]
      
      tone$xpp$strobe3[wndx]  <- tone$xpp$strobe3[tndx]
      tone$xpp$stimType       <- tone$xpp$strobe3
      
      brksTmp <- seq(1.5, by=1, to=max(tone$xpp$freqInd)+0.5)
      brks                   <- c(0, brksTmp, max(tone$xpp$freqInd)+2.5 )# add lower bound, separation to silence and upper bound
      
      NfrqInd             <- length(unique(tone$xpp$freqInd))
      tone$xpp$freqIndAll <- tone$xpp$freqInd + NfrqInd * (tone$xpp$stimType-1)
      
      tone$xpp$category     <- cut( tone$xpp$freqIndAll, brks )
      tone$xpp$catIndx      <- as.numeric( tone$xpp$category  )
      
      #levels(mua$xpp$category) <- seq(1:length(unique(mua$xpp$category)))
    }
    
    
    Namp <- length(unique(tone$xpp$strobe3))-1
    if(tone$project=='narrowBandNoise')
      Namp <- 5
    
    tone$xpp$ampInd   <- 1+(tone$xpp$strobe3-1)%%Namp
    tone$xpp$stimType <- 1+(tone$xpp$strobe3-1)%/%Namp
    tone$xpp$ampInd[tone$xpp$strobe3==max(tone$xpp$strobe3)] <- Namp+1
    
    #
    tone$xpp$toneNr   <- 1:dim(tone$xpp)[1]
    
    # ------------ frequency
    freqList          <- tone$dfline$minOct*exp(log(2)*seq(0,by=tone$dfline$stepOct, length=tone$dfline$NFrq))
    tone$xpp$freqHz   <- freqList[ tone$xpp$freqInd ]
    tone$xpp$freqOct  <- log2( freqList[ tone$xpp$freqInd ])
    
    # -------------- amplitude
    ampList           <- tone$dfline[grep(  'dB' , names(tone$dfline))] 
    tone$xpp$ampdB    <- c(unlist(ampList[ tone$xpp$ampInd ]))
    
    # sequence number
    tonePerTrial   <- table(tone$xpp$trialNr)
    tone$xpp$seqNr <- c(unlist( sapply(tonePerTrial, function(n){1:n}) ))
  }
  
  ## RFmap_noise          ====
  if ( tone$project == 'RFmap'){
    tone$xpp$ampInd <- tone$xpp$strobe2
  }
  
  
  
  
  ## DRF                  ====
  if ( tone$project == 'DRF'){
    tone$xpp$freqInd <- tone$xpp$strobe2
    tone$xpp$durInd  <- tone$xpp$strobe3
    #
    tone$xpp$toneNr   <- 1:dim(tone$xpp)[1]
    
    # ------------ frequency
    freqList          <- tone$dfline$minOct*exp(log(2)*seq(0,by=tone$dfline$stepOct, length=tone$dfline$NFrq))
    tone$xpp$freqHz   <- freqList[ tone$xpp$freqInd ]
    tone$xpp$freqOct  <- log2( freqList[ tone$xpp$freqInd ])
    
    # -------------- duration
    durList           <- tone$dfline[grep(  'dur' , names(tone$dfline))] 
    tone$xpp$durMs    <- c(unlist(durList[ tone$xpp$durInd ]))
    
    # sequence number
    tonePerTrial   <- table(tone$xpp$trialNr)
    tone$xpp$seqNr <- c(unlist( sapply(tonePerTrial, function(n){1:n}) ))
  }
  
  ## DTRF                 ====
  if (tone$dfline$link=='DTRF'){
    
    tone$xpp$freqInd <- tone$xpp$strobe2
    tone$xpp$freq2Ind  <- tone$xpp$strobe3
    
    # ------------ frequency 1
    freqList          <- tone$dfline$minOct*exp(log(2)*seq(0,by=tone$dfline$stepOct, length=tone$dfline$NFrq))
    tone$xpp$freqHz   <- freqList[ tone$xpp$freqInd ]
    tone$xpp$freqOct  <- log2( freqList[ tone$xpp$freqInd ])
    
    # ------------ frequency 2
    freqList2          <- tone$dfline$minOct2*exp(log(2)*seq(0,by=tone$dfline$stepOct2, length=tone$dfline$NFrq2))
    tone$xpp$freq2Hz   <- freqList2[ tone$xpp$freq2Ind ]
    tone$xpp$freq2Oct  <- log2( freqList2[ tone$xpp$freq2Ind ])
    
    # sequence number
    tonePerTrial   <- table(tone$xpp$trialNr)
    tone$xpp$seqNr <- c(unlist( sapply(tonePerTrial, function(n){1:n}) ))
    
  }
  
  ## sweepRF              ====
  if ( tone$project %in% c('sweepRF') ){
    
    tone$xpp$freqInd <- tone$xpp$strobe2
    
    # old version
    #tone$xpp$ampInd  <- tone$xpp$strobe3
    tone$xpp$sweepInd   <- 1+(tone$xpp$strobe3-1)%%5
    tone$xpp$stimType <- 1+(tone$xpp$strobe3-1)%/%5
    tone$xpp$sweepInd[tone$xpp$strobe3==max(tone$xpp$strobe3)] <- 6
    
    #
    tone$xpp$toneNr   <- 1:dim(tone$xpp)[1]
    
    # ------------ frequency
    freqList          <- tone$dfline$minOct*exp(log(2)*seq(0,by=tone$dfline$stepOct, length=tone$dfline$NFrq))
    tone$xpp$freqHz   <- freqList[ tone$xpp$freqInd ]
    tone$xpp$freqOct  <- log2( freqList[ tone$xpp$freqInd ])
    
    # -------------- amplitude
    ampList               <- tone$dfline[grep(  'dB' , names(tone$dfline))] 
    tone$xpp$sweepAmount  <- c(unlist(ampList[ tone$xpp$sweepInd ]))
    
    # sequence number
    tonePerTrial   <- table(tone$xpp$trialNr)
    tone$xpp$seqNr <- c(unlist( sapply(tonePerTrial, function(n){1:n}) ))
  }
  
  ## MFF                  ====
  if ( tone$project == 'MFF'){
    
    tone$xpp$freqInd <- tone$xpp$strobe2
    tone$xpp$ampInd  <- tone$xpp$strobe3
    #
    tone$xpp$toneNr   <- 1:dim(tone$xpp)[1]
    
    # ------------ frequency
    freqList          <- tone$dfline$minOct*exp(log(2)*seq(0,by=tone$dfline$stepOct, length=tone$dfline$NFrq))
    tone$xpp$freqHz   <- freqList[ tone$xpp$freqInd ]
    tone$xpp$freqOct  <- log2( freqList[ tone$xpp$freqInd ])
    
    # -------------- amplitude
    ampList           <- tone$dfline[grep(  'dB' , names(tone$dfline))] 
    tone$xpp$ampdB    <- c(ampList[ tone$xpp$ampInd ])
    
    # sequence number
    tonePerTrial   <- table(tone$xpp$trialNr)
    tone$xpp$seqNr <- c(unlist( sapply(tonePerTrial, function(n){1:n}) ))
  }
  
  ## clickGroups          ====
  if ( tone$project %in% c('clickGroups')){
    
    tone$xpp$deviant        <- tone$xpp$strobe2
    tone$xpp$toneType       <- tone$xpp$strobe3
  }
  
  ## harmonicity          ====
  if ( tone$project %in% c('harmonicity')){
    
    tone$xpp$freqInd    <- tone$xpp$strobe2
    tone$xpp$jitterInd  <- tone$xpp$strobe3
    
    if (tone$dfline$link[1] == 'harmonicity_coo_link.txt'){
      tone$xpp$toneID    <- tone$xpp$strobe2
      tone$xpp$freqInd    <- tone$xpp$strobe3
      tone$xpp$jitterInd  <- tone$xpp$strobe4
      tone$xpp$jitterInd[tone$xpp$toneID>=16] <- length(unique(tone$xpp$jitterInd))+1
      tone$xpp$freqInd[tone$xpp$toneID>=16]   <- length(unique(tone$xpp$freqInd))+1
    }
    
    tone$xpp$toneNr   <- 1:dim(tone$xpp)[1]
    
    # ------------ frequency
    freqList          <- tone$dfline$minOct*exp(log(2)*seq(0,by=tone$dfline$stepOct, length=tone$dfline$NFrq))
    tone$xpp$freqHz   <- freqList[ tone$xpp$freqInd ]
    tone$xpp$freqOct  <- log2( freqList[ tone$xpp$freqInd ])
    
    # -------------- amplitude
    ampList            <- tone$dfline[grep(  'dB' , names(tone$dfline))] 
    tone$xpp$jitter    <- unlist(c(ampList[ tone$xpp$jitterInd ]))
    
    # sequence number
    tonePerTrial   <- table(tone$xpp$trialNr)
    tone$xpp$seqNr <- c(unlist( sapply(tonePerTrial, function(n){1:n}) ))
  }
  
  
  ## textures             ====
  if ( tone$project %in% c('textures')){
    tone$xpp$texID    <- tone$xpp$strobe2
    tone$xpp$cndID    <- tone$xpp$strobe3
  }
  
  
  ## Change Detection     ====
  if ( tone$project %in% c('changeDetection','detectTone')){
    
    tone$trialInfo$standardNumber  <- tone$trialInfo$trialStrobe2
    tone$trialInfo$absSig          <- tone$trialInfo$trialStrobe3
    tone$trialInfo$toneDuration    <- tone$trialInfo$trialStrobe4
    if ( any(unique(tone$trialInfo$toneDuration)== -1) )
      tone$trialInfo$toneDuration <- tone$dfline$toneDur
    
    tone$trialInfo$signal           <- (tone$trialInfo$absSig-1)%%5 + 1
    tone$trialInfo$signalSigned     <- c(0,1,2,3,4,0,-1,-2,-3,-4)[tone$trialInfo$absSig]
    
    sigValues            <- tone$dfline[1,c('S1','S2','S3','S4','S5')]
    sigValues            <- cbind(sigValues, -1*sigValues)
    tone$trialInfo$signalStrength   <- t(sigValues[tone$trialInfo$signal])
    tone$trialInfo$signalStrength   <- t(sigValues[tone$trialInfo$absSig])
    
    ## detectTone ====
    if (tone$project == 'detectTone'){
      dBsum <- function(xdB){   20*log10( sum( 10^(xdB/20)) )   }
      dBinv <- function(xdB){  10^(xdB/20)}
      
      
      dBWN   <- 60  #40dB with A-level, 60dB with C-level
      sigdB  <- 74
      
      getSignalToNoise <- function(sig, noisedB=dBWN,sinedB=sigdB){ 20*log10( dBinv( sinedB + 20*log10( sig ) )   /  dBinv( noisedB  ) ) }
      sapply(sigValues,getSignalToNoise,noisedB=dBWN)
      
      sigValuesAdj     <- sigValues
      #sigValuesAdj[1]  <- min(sigValues[sigValues>0])/
      sigValuesdB      <- sapply(sigValues,getSignalToNoise,noisedB=dBWN)
      sigValuesdBAdj   <- sigValuesdB
      sigValuesdBAdj[1] <- sigValuesdB[2] - 1.5*diff(sigValuesdB[c(2,3)])
      
      #-------------------------------------------------------------
    }
    
    # identify catch trials 
    tone$trialInfo$catchTrial <- tone$trialInfo$signal == 1
    
    
    tone$trialInfo$trialDuration    <- tone$trialInfo$trialEndTime - tone$trialInfo$firstToneTime 
    
    ## determine if trial ended with a lever release or not
    tone$trialInfo$leverRelease                                                <- rep(NA, dim(tone$trialInfo)[1])
    
    # signal present trials ended with lever release if trialerror is 0, 2, 5, 6, 7
    tone$trialInfo$leverRelease[tone$trialInfo$trialError%in%c(0,5,6,7) & tone$trialInfo$signal>1.5] <- TRUE
    tone$trialInfo$leverRelease[tone$trialInfo$trialError%in%c(2,4)     & tone$trialInfo$signal>1.5] <- FALSE
    
    # correct signal absent trials ended without lever release
    tone$trialInfo$leverRelease[tone$trialInfo$trialError==0          & tone$trialInfo$signal<1.5]   <- FALSE
    tone$trialInfo$leverRelease[tone$trialInfo$trialError%in%c(5,6,7) & tone$trialInfo$signal<1.5]   <- TRUE
    
    
    ## 
    tone$trialInfo$leverReleaseTime                                       <- rep(NA, dim(tone$trialInfo)[1])
    tone$trialInfo$leverReleaseTime[ which(tone$trialInfo$leverRelease) ] <- tone$trialInfo$trialEndTime[ which(tone$trialInfo$leverRelease) ]
    #tone$trialInfo$lastToneTime                                           <- apply(times,1,max) 
    
    #
    tone$trialInfo$RTalt                                       <- tone$trialInfo$trialEndTime - tone$trialInfo$lastToneTime 
    #tone$trialInfo$RTtar                                       <- tone$trialInfo$trialEndTime - tone$trialInfo$targetTime 
    
    tone$trialInfo$RTalt[which(!tone$trialInfo$leverRelease)]             <- NA
    #tone$trialInfo$RTtar[which(!xpp$leverRelease)]             <- NA
    tone$trialInfo$RT[   which(!tone$trialInfo$leverRelease)]             <- NA
    
    #vnd <- which(xpp$signal>1.5 & xpp$trialError==0)
    
    ## ======================= check for unrealistically short SOAs if they exist
    tmpSOA                <- c(NA,diff(tone$xpp$e.onset))
    tmpSOA[1]             <- 20
    rmIndx <- which(tmpSOA < 0.001 )
    
    if (length(rmIndx)>0){
      warning(paste('removing trials with zero SOA: ', tone$dfline$session, 'trial# ', rmIndx, sep=' '))
      #print(tone$xpp[rmIndx,])
      
      # first adjust seqLen entry in trialInfo
      for (rx in rmIndx){
        trialNr <- tone$xpp$trialNr[rx]
        tiL     <- which(tone$trialInfo$trialNr == trialNr)
        tone$trialInfo$seqLen[tiL] <- tone$trialInfo$seqLen[tiL] - 1  
      }
      tone <- subs.data( tone, !tmpSOA<.001, saveRejected=F)
    }
    
    ## 
    trialNr <- rep(tone$trialInfo$trialNr[1],tone$trialInfo$seqLen[1])
    seqNr <- c()
    #seqNr   <- 1:tone$trialInfo$seqLen[1]
    for (i in 1:dim(tone$trialInfo)[1]){
      trialNr <- c(trialNr, rep(tone$trialInfo$trialNr[i],tone$trialInfo$seqLen[i]) )
      if (tone$trialInfo$seqLen[i]>0)
        seqNr   <- c(seqNr, 1:tone$trialInfo$seqLen[i])
    }
    tone$xpp$seqNr <- seqNr
    seqNrBinBrks <- c(0.5, 1.5, 2.5, 4.5, 7.5, 14.5)
    #seqNrBinBrks <- c(0.5, 1.5, 2.5, 4.5, 6.8, 8.5, 13.5)
    tone$xpp$seqNrBin <- cut(tone$xpp$seqNr, breaks=seqNrBinBrks)
    
    seqNrBinBrksCoarse      <- c(0, 1.5, 3.5, 10.5, 14)
    tone$xpp$seqNrBincoarse <- cut(tone$xpp$seqNr, breaks=seqNrBinBrksCoarse)
    
    seqNrBinBrksCoarse      <- c(0, 1.5, 3.5, 10.5, 14)
    tone$xpp$seqNrBincoarse <- cut(tone$xpp$seqNr, breaks=seqNrBinBrksCoarse)
    
    seqNrBinBrks2 <- c(0, 1.5, 2.5, 14)
    tone$xpp$seqBin2 <- cut(tone$xpp$seqNr, breaks=seqNrBinBrks2)
    
    seqNrBinBrks3 <- c(0, 1.5, 2.5, 3.5, 14)
    tone$xpp$seqBin3 <- cut(tone$xpp$seqNr, breaks=seqNrBinBrks3)
    
    
    disc               <- as.numeric(tone$xpp$strobe2==1)
    if (tone$project=='detectTone')
      disc               <- as.numeric(tone$xpp$strobe2>11)
    
    
    tone$xpp$target     <- diff(c(0,disc))==1 
    tone$xpp$postTarget <- tone$xpp$strobe2==1 & !tone$xpp$target
    
    tone$xpp$signal                        <- rep(1,dim(tone$xpp)[1])
    tone$xpp$signal[which(tone$xpp$target)] <- tone$trialInfo$signal[tone$xpp$trialNr[which(tone$xpp$target)]]
    
    tone$xpp$signalStrength                         <- rep(0,dim(tone$xpp)[1])
    tone$xpp$signalStrength[which(tone$xpp$target)] <- tone$trialInfo$signalStrength[tone$xpp$trialNr[which(tone$xpp$target)]]
    
    if (tone$project=='detectTone')
      tone$xpp$signalStrengthdB <- sigValuesdBAdj[tone$xpp$signal]
    
    
    tone$xpp$deltaPitchStndrd                          <- rep(0,dim(tone$xpp)[1])
    tone$xpp$deltaPitchStndrd[which(!tone$xpp$target)] <- tone$trialInfo$signalStrength[tone$xpp$trialNr[which(!tone$xpp$target)]]
    
    tone$xpp$catchTrial                                <- tone$trialInfo$catchTrial[ tone$xpp$trialNr ]
    tone$xpp$standardNumber                            <- tone$trialInfo$standardNumber[tone$xpp$trialNr]
    tone$xpp$targetNumber                              <- tone$trialInfo$standardNumber[tone$xpp$trialNr]
    # not sure why I originally used 'standardNumber' should be 'targetNumber' for all I can tell
    # leaving 'standardNumber' identical to 'targetNumber' for backwards compatibility
    
    tone$xpp$toneDuration                     <- tone$trialInfo$toneDuration[tone$xpp$trialNr]
    tone$xpp$ampInd                           <- tone$trialInfo$toneDuration[tone$xpp$trialNr]
    
    tone$xpp$toneDurationBin <- as.factor(tone$xpp$toneDuration)
    tone$xpp$ampIndBin       <- as.factor(tone$xpp$ampInd)
    
    tone$xpp$responded           <- rep(FALSE,dim(tone$xpp)[1])
    tone$xpp$leverReleaseTime    <- tone$trialInfo$leverReleaseTime[tone$xpp$trialNr]
    tone$xpp$trialError          <- tone$trialInfo$trialError[tone$xpp$trialNr]
    
    tone$xpp$deltaRT            <- tone$xpp$leverReleaseTime - tone$xpp$e.MLonset
    tone$xpp$responded          <- tone$xpp$deltaRT > 100 & tone$xpp$deltaRT < 700
    tone$xpp$responded[which(is.na(tone$xpp$leverReleaseTime))] <- 0
    
    #    browser()
    
    #tone$signalStrengthdB <- sigValuesdB[tone$signal]
    tone$xpp$SOA                <- c(NA,diff(tone$xpp$e.onset))
    tone$xpp$ISI                <- tone$xpp$SOA - 0.001*tone$xpp$toneDuration
    tone$xpp$SOA[tone$xpp$seqNr==1] <- NA
    tone$xpp$ISI[tone$xpp$seqNr==1] <- NA
    tone$xpp$logSOA             <- log2(tone$xpp$SOA)     
    tone$xpp$logISI             <- log2(tone$xpp$ISI)     
    tone$xpp$SOAlag1            <- my.lag(tone$xpp$SOA,1)
    tone$xpp$ISIlag1            <- my.lag(tone$xpp$ISI,1)
    tone$xpp$futureSOA          <- my.lag(tone$xpp$SOA,-1)
    
    tone$xpp$logSOAlag          <- log2(tone$xpp$SOAlag) 
    tone$xpp$logISIlag          <- log2(tone$xpp$ISIlag) 
    
    
    # frequency of the 11 target tones
    #targetPitchOct      <- seq( log2(900), log2(2750), length=11) 
    targetPitchOct      <- seq( log2(as.numeric(tone$dfline$minF)), by=as.numeric(tone$dfline$deltaF), length=11)
    
    # pitch of each tone in log 2
    tone$xpp$tonePitchOct   <- targetPitchOct[ tone$xpp$targetNumber ] + tone$xpp$deltaPitchStndrd
    #tone$tonePitchOct   <- targetPitchOct[ tone$targetNumber ] - tone$deltaPitchStndrd # wrong way
    
    # pitch standard tone in log 2 for each trial
    tone$xpp$standardPitchOct <- tone$xpp$tonePitch
    tone$xpp$standardPitchOct[which(tone$xpp$target)] <- tone$xpp$standardPitchOct[which(tone$xpp$target)-1]
    
    tone$xpp$standardInRng <- tone$xpp$standardPitch >= min(targetPitchOct) & tone$xpp$standardPitch <= max(targetPitchOct) 
    
    # repeat the whole thing with mean-normalized values
    targetPitchOctNrm   <- targetPitchOct - mean(targetPitchOct)
    tone$xpp$tonePitchNrm   <- targetPitchOctNrm[ tone$xpp$targetNumber ] + tone$xpp$deltaPitchStndrd
    
    tone$xpp$standardPitchNrm <- tone$xpp$tonePitchNrm
    tone$xpp$standardPitchNrm[which(tone$xpp$target)] <- tone$xpp$standardPitchNrm[which(tone$xpp$target)-1]
    
    
    brks                  <- .250*2^seq(0,6,by=1)
    brks[length(brks)]    <- 20.000
    tone$xpp$soaBin       <- cut(tone$xpp$SOA,brks)
    tone$xpp$soaBinLag1   <- my.lag(tone$xpp$soaBin,1)
    tone$xpp$futureSOAbin <- my.lag(tone$xpp$soaBin,-1)
    
    
    brks                <- .250*2^seq(0,6,by=1/2)
    brks[length(brks)]  <- 20.000
    tone$xpp$soaBinfine <- cut(tone$xpp$SOA,brks)
    
    brks                  <- .250*2^seq(0,6,by=2)
    brks[length(brks)]    <- 20.000
    tone$xpp$soaBincoarse <- cut(tone$xpp$SOA,brks)
    
    
    brks                <- 0.0625*2^seq(0,9,by=1)
    brks[length(brks)]  <- 20.000; brks[1] <- .050
    tone$xpp$isiBin     <- cut(tone$xpp$ISI,brks)
    tone$xpp$isiBinLag1 <- my.lag(tone$xpp$isiBin,1)
    
    
    
    
    ##What is this doing here?
    vnd <- which(tone$xpp$responded==T)
    
    minTarget <- min(tone$xpp$seqNr[which(tone$xpp$target)])
    tone$xpp$preTarget <- tone$xpp$seqNr < minTarget-.5
    
    
    ## ================ 
    # identify the distribution of reaction times when there is only one valid option to choose from
    tmp           <- tapply(tone$xpp$responded, tone$xpp$trialNr, sum)
    singleTrialNr <- as.numeric(names(which(tmp==1))) 
    snd           <- which(tone$xpp$trialNr%in%singleTrialNr & tone$xpp$responded)
    valDist       <- sort(tone$xpp$deltaRT[snd])
    
    # identify all trials that have more than one valid RT option
    doubleTrialNr <- as.numeric(names(which(tmp>1))) 
    print(doubleTrialNr)
    
    
    # function to resolve which RT is more likely based on the valDist
    getLikelyResponseTime <- function(n){
      RToption   <- tone$xpp$deltaRT[which(tone$xpp$trialNr==n & tone$xpp$responded) ]
      #sgnlOption <- tone$signal[which(tone$trialNr==n & tone$responded) ]
      dn         <- dnorm(RToption, m=mean(valDist), s=sd(valDist))/dnorm(mean(valDist), m=mean(valDist), s=sd(valDist) )
      #dn         <- dn*(sgnlOption)
      return(  c(2*as.numeric(dn==max(dn))-1)  )
    }
    
    for (j in doubleTrialNr){
      print(j)
      vj                     <- which(tone$xpp$trialNr==j & tone$xpp$responded)
      tone$xpp$responded[vj] <- getLikelyResponseTime(j)
    }
    tone$xpp$rawResponded <- tone$xpp$responded
    tone$xpp$responded    <- tone$xpp$responded %>=% 0
    
    
    ## ----------------------------------------
    getSeqNrTarget <- function(n){
      snt <- NA
      valnd <- which(tone$xpp$trialNr == n)
      if (length(valnd)>0){
        if(sum(tone$xpp$target[valnd]==1)>0){
          snt <- which( tone$xpp$target[valnd] == 1 )
        }
      }
      return(snt)
    }
    
    tone$xpp$seqNrTarget <- sapply(1:dim(tone$xpp)[1], function(n){getSeqNrTarget(tone$xpp$trialNr[n])} )
    tone$xpp$seqNrToTarget <- tone$xpp$seqNr - tone$xpp$seqNrTarget
    
    # identify three types of trial: catchTrialNr: 0) signal present trials, 1) catch trials, 3) aborted pre-target
    #tone$catchTrialNdx <- 
    
    getSeqNrResponse <- function(n){
      snt <- NA
      valnd <- which(tone$xpp$trialNr == n)
      if (length(valnd)>0){
        if(sum(tone$xpp$responded[valnd]==1)>0){
          snt <- which( tone$xpp$responded[valnd] == 1 )
        }
      }
      return(snt)
    }
    
    tone$xpp$seqNrResponse <- sapply(1:dim(tone$xpp)[1], function(n){getSeqNrResponse(tone$xpp$trialNr[n])} )
    tone$xpp$seqNrToResponse <- tone$xpp$seqNr - tone$xpp$seqNrResponse
    
    tone$xpp$postResponse <- FALSE
    tone$xpp$postResponse[which(tone$xpp$seqNrToResponse > 0.5)] <- TRUE
    
    
    getCSOA <- function(n){
      csoa <- c()
      valnd <- which(tone$xpp$trialNr == n)
      if (length(valnd)>0){
        tmp    <- tone$xpp$SOA[valnd]
        tmp[1] <- 0
        csoa <- cumsum( tmp )
      }
      return(csoa)
    }
    tone$xpp$CSOA <- unlist( sapply(unique(tone$xpp$trialNr), getCSOA  ))
    
    
    valTrial  <- which( tone$trialInfo$trialError %in% c(0,2,5,6,7) & tone$trialInfo$seqLen>0 )#& tone$trialInfo$trialNr >= trNrRng[1] & tone$trialInfo$trialNr <= trNrRng[2])
    valTones  <- which(tone$trialNr %in% valTrial)
    
    tone <- subs.data(tone, tone$xpp$trialNr %in% valTrial)
    
  }
  
  ## temporal Discounting ====
  if ( tone$project == 'temporalDiscounting'){
    tone$xpp$leverReleaseTime    <- tone$trialInfo$trialEndTime[tone$xpp$trialNr]
    tone$xpp$trialError          <- tone$trialInfo$trialError[tone$xpp$trialNr]
    tone$xpp$blockNumber         <- tone$trialInfo$BlockNumber[tone$xpp$trialNr]
    tone$xpp$conditionNumber     <- tone$trialInfo$ConditionNumber[tone$xpp$trialNr]
    
    delayValues            <- tone$dfline[1,c('delay1','delay2','delay3','delay4','delay5')]
    
    
    # red target is high reward
    tone$xpp$redTargetPos    <- 1 +  (tone$xpp$conditionNumber-1)%%2        # left  or right of fixation
    tone$xpp$redCuePos       <- 1 + ((tone$xpp$conditionNumber-1)%/%2)%%2   # above or below of fixation
    
    # infer response
    tone$xpp$numReward     <-   tone$trialInfo$numReward[tone$xpp$trialNr]
    #tone$xpp$choiceType    <- 1 + as.numeric(tone$trialInfo$numReward==max(tone$xpp$numReward))
    tone$xpp$choiceType    <- 1 + as.numeric(tone$trialInfo$trialError[tone$xpp$trialNr]==0)
    
    tone$xpp$choiceDir                                   <- tone$xpp$redTargetPos
    tone$xpp$choiceDir[ tone$xpp$choiceType==1 ]         <- 3 - tone$xpp$choiceDir[ tone$xpp$choiceType==1 ]
    
    #table(tone$xpp$redTargetPos, tone$xpp$choiceDir, tone$xpp$blockNumber)
    #table(tone$xpp$choiceType, tone$xpp$blockNumber)
    
    tone$xpp$RT                 <- tone$trialInfo$RT[tone$xpp$trialNr]
    tone$xpp$deltaRT            <- tone$xpp$leverReleaseTime - tone$xpp$e.MLonset
    
    
    
    tone$xpp$delay <- t(delayValues[tone$xpp$blockNumber])
    
    tone$xpp$choiceType[ !tone$xpp$trialError%in%c(0,6) ] <- NA
    tone$xpp$choiceDir[  !tone$xpp$trialError%in%c(0,6) ] <- NA
    tone$xpp$RT[         !tone$xpp$trialError%in%c(0,6) ] <- NA
    
    tone$xpp$smallReward <- tone$xpp$choiceType == 2
    tone$xpp$delayCue    <- tone$xpp$strobe1 == 4
    tone$xpp$responseCue <- tone$xpp$strobe1 == 5
    
    #tapply(tone$xpp$RT, tone$xpp$choiceDir, na.mean)
    #tapply(tone$xpp$RT, tone$xpp$choiceType, na.mean)
    #tapply(tone$xpp$RT, list(tone$xpp$choiceType,tone$xpp$blockNumber), na.mean)
    
  }
  
  ## ampOddClick          ====
  if ( tone$project == 'ampOddClick'){
    
    #    browser()
    #I believe we want the intan record here, for accuracy, but that is already 
    #on the second scale, not ms like ML
    #    tone$xpp$ICI                 <- c(20,.001*diff(tone$xpp$e.onset))
    tone <- checkISISignalProcessor( tone )
    
    tone$xpp$ampInd              <- ((tone$xpp$strobe2-1) %% 11)+1
    tone$xpp$trialType <- tone$trialInfo$trialStrobe2[ tone$xpp$trialNr ]
    
    amplitude              <- seq(from=62,by=6,length=5) 
    tone$xpp$dB            <- amplitude[tone$xpp$ampInd]
    
    pastAmpInd             <- my.lag(tone$xpp$ampInd,1)
    pastAmpInd[1]          <- NA
    tone$xpp$deltaAmpInd   <- tone$xpp$ampInd - pastAmpInd 
    
    tone$xpp$ICI                 <- c(20, diff(tone$xpp$e.onset))
    tone$xpp$ICI[tone$xpp$ICI<0] <- c(20)
    tone$xpp$ICI                 <- tone$xpp$ICI %<=% 20
    
    brks <- 2^c(-3,-1,0,1,2,3,5)
    tone$xpp$ICIoct              <- cut(tone$xpp$ICI, breaks=brks)
    
    tone$xpp$ISI    <- tone$xpp$ICI
    tone$xpp$isiOct <- tone$xpp$ICIoct
    tone$xpp$iciInd <- as.numeric(tone$xpp$ICIoct)
    
    #tone$xpp$logISI <- log2(tone$xpp$ISI)
    
    tone$xpp$deltaICI     <- c(NA, diff(log2(tone$xpp$ICI)))
    tone$xpp$absDeltaICI  <- abs(tone$xpp$deltaICI)
    
    brks    <- seq(-4.5,by=1,4.5)
    brks[1] <- -10; brks[length(brks)]<- 10 
    tone$xpp$deltaICIbin  <- cut(tone$xpp$deltaICI,  breaks=brks)
    
    brks    <- seq(-0.25,by=0.5,4.5)
    brks[length(brks)]<- 10 
    tone$xpp$absDeltaICIbin  <- cut(tone$xpp$absDeltaICI, breaks=brks)
    
    tone$xpp$ampDev  <- abs(tone$xpp$deltaAmpInd )     >  0       
    tone$xpp$timeDev <- abs(tone$xpp$deltaICI)         >= 0.1  # time Deviant if ISI differs by more than 5% from previous ISI
    
    
    tttmp             <- tone$xpp$ampDev | tone$xpp$timeDev
    tttmp[1]          <- TRUE
    tone$xpp$toneType <- as.numeric(!tttmp)
    
    rpCnt             <- as.numeric(tone$xpp$toneType)==1
    tmp               <- get.epoch( as.numeric(tone$xpp$toneType)==1  )
    tmp$dur           <- tmp$offset - tmp$onset
    rpCnt[tmp$offset] <- -tmp$dur  
    tone$xpp$repNr    <- cumsum(rpCnt)
    
    # calculate seqNr
    tmpSeqNr       <- 1:length(tone$xpp$trialNr)
    numToneBlock   <- table(tone$xpp$trialNr)
    subtractList   <- c(0, cumsum(numToneBlock))
    tone$xpp$seqNr <- tmpSeqNr - subtractList[tone$xpp$trialNr]
    
    #tone$xpp$p2p                 <- apply(tone$lfp,1, max) - apply(tone$lfp,1, min)
    
  }
  
  
  ## RSkernel             ====
  if ( tone$project == 'RSkernel'){
    
    tone$xpp$ICI                 <- c(20,diff(tone$xpp$e.onset))
    tone$xpp$ICI[tone$xpp$ICI<0] <- c(20)
    tone$xpp$ICI                 <- tone$xpp$ICI %<=% 20
    brks <- 2^c(-3,-1,0,1,2,3,5)
    
    tone$xpp$ICIoct              <- cut(tone$xpp$ICI, breaks=brks)
    
    tone$xpp$ISI    <- tone$xpp$ICI
    tone$xpp$isiOct <- tone$xpp$ICIoct
    
    tone$xpp$isiLog  <- log(tone$xpp$ISI)
    
    #tone$xpp$deltaISI         <- c(NA, diff(log2(tone$xpp$ISI)))
    #tone$xpp$absDeltaISI      <- abs(tone$xpp$deltaISI)
    tone$xpp$ISI[ which( tone$xpp$ISI>20)]  <- 20
    tone$xpp$ISI[ which( tone$xpp$ISI<0)]   <- 20
    tone$xpp$deltaISI                       <- c(NA, diff(log2(tone$xpp$ISI)))
    tone$xpp$absDeltaISI                    <- abs(tone$xpp$deltaISI)
    
    brks                  <- seq(-4.5,by=1,4.5)
    brks[1]               <- -10; brks[length(brks)]<- 10 
    tone$xpp$deltaISIbin  <- cut(tone$xpp$deltaISI,  breaks=brks)
    
    brks    <- seq(-0.25,by=0.5,4.5)
    brks[length(brks)]<- 10 
    tone$xpp$absDeltaISIbin  <- cut(tone$xpp$absDeltaISI, breaks=brks)
    
    
    tone$xpp$p2p                 <- apply(tone$lfp,1, max) - apply(tone$lfp,1, min)
    
    tone$xpp$freqInd              <- ((tone$xpp$strobe2-1) %% 11)+1
    tone$xpp$toneInd              <- tone$xpp$freqInd # for backwards compatibility with loadRS
    tone$xpp$pitchInd             <- tone$xpp$freqInd # for backwards compatibility with loadRS
    
    tone$xpp$trialType <- tone$trialInfo$trialStrobe2[ tone$xpp$trialNr ]
    
    # calculate seqNr
    tmpSeqNr <- 1:length(tone$xpp$trialNr)
    numToneBlock   <- table(tone$xpp$trialNr)
    subtractList <- c(0, cumsum(numToneBlock))
    
    tone$xpp$seqNr <- tmpSeqNr - subtractList[tone$xpp$trialNr]
    
    
    # prepare xpp (from functions_RSdual) =====================================
    pitch         <- tone$dfline$minOct*exp(log(2)*seq(0,by=tone$dfline$stepOct, length=11))
    tone$xpp$toneCount        <- 1:dim(tone$xpp)[1]
    
    pastPitchInd              <- my.lag(tone$xpp$toneInd,1)
    pastPitchInd[1]           <- NA
    tone$xpp$deltaPitchInd    <- tone$xpp$toneInd - pastPitchInd 
    tone$xpp$absDeltaPitchInd <- abs(tone$xpp$deltaPitchInd)
    tone$xpp$deltaPitch       <- tone$xpp$deltaPitchInd # for compatibility with FLA_VProbe_STPSP 
    
    tone$xpp$pitchIndFactor <- as.factor(tone$xpp$pitchInd)
    
    tone$xpp$pitchFreq        <- pitch[ tone$xpp$toneInd ]
    tone$xpp$pitchOct         <- log(pitch[ tone$xpp$toneInd  ])/log(2)
    pastPitchOct              <- my.lag(tone$xpp$pitchOct,1)
    
    pastPitchOct[1]           <- NA
    tone$xpp$deltaPitchOct    <- signif( tone$xpp$pitchOct - pastPitchOct, 4) 
    tone$xpp$absDeltaPitchOct <- abs(tone$xpp$deltaPitchOct)
    
    ## tone Type
    tone$xpp$toneType    <- as.numeric( (tone$xpp$strobe2 > 11 & tone$xpp$deltaPitchInd!=0) + 1 )
    tone$xpp$toneType[1] <- 3
    tone$xpp$toneType <- as.factor(tone$xpp$toneType)
    
    #tone$xpp$mmn         <- tone$xpp$n1
    
    ## create a repetition counter that gets reset every time tone Type is equal to 1, ie. each time there is a pitch deviant 
    rpCnt             <- as.numeric(tone$xpp$toneType)==1
    tmp               <- get.epoch( as.numeric(tone$xpp$toneType)==1  )
    tmp$dur           <- tmp$offset - tmp$onset
    rpCnt[tmp$offset] <- -tmp$dur  
    if ( tone$xpp$toneType[length(rpCnt)]==1 ) # handle false end-of-series marker at the end of the experimental session
      rpCnt[length(rpCnt)] <- 1
    
    tone$xpp$repNr    <- cumsum(rpCnt)
    
    
    
    tone$xpp$toneDeviant <- abs(tone$xpp$deltaPitchOct) > .01
    tone$xpp$timeDeviant <- abs(tone$xpp$deltaISI)      > .02
    tone$xpp$anyDeviant  <- tone$xpp$toneDeviant | tone$xpp$timeDeviant
    
    ## create a repetition counter that gets reset every time there is a time deviant 
    rm(rpCnt)
    tone$xpp$timeDeviant[1] <- TRUE
    rpCnt             <- as.numeric(!tone$xpp$timeDeviant)
    tmp               <- get.epoch( !tone$xpp$timeDeviant  )
    tmp$dur           <- tmp$offset - tmp$onset
    rpCnt[tmp$offset] <- -tmp$dur  
    tone$xpp$repNrTimeDeviant    <- cumsum(rpCnt)
    
    
    # create a new variable that codes trialType and toneType in one factor: ND
    # ND1 : regular, standard
    # ND2 : random , standard
    # ND3 : random , deviant
    # ND4 : regular, deviant
    # trialType:             1: random  ; 2: regular
    # toneType:              1: standard; 2: deviant
    
    #tmp             <- rep('ND1', dim(tone$xpp)[1])
    #tst             <- tapply( tone$xpp$absDeltaPitchInd, list(tone$xpp$trialType,tone$xpp$toneType) )# 1:rnd,std; 2:reg,std; 3:rnd,dev; 4: reg:dev
    #tone$xpp$NDind  <- c(2,1,3,4)[tst]
    
    tone$xpp$NDind <- c(NA)
    tone$xpp$NDind[ tone$xpp$trialType == 1 & tone$xpp$toneType==1] <- 1
    tone$xpp$NDind[ tone$xpp$trialType == 2 & tone$xpp$toneType==1] <- 2
    tone$xpp$NDind[ tone$xpp$trialType == 1 & tone$xpp$toneType==2] <- 3
    tone$xpp$NDind[ tone$xpp$trialType == 2 & tone$xpp$toneType==2] <- 4
    
    tone$xpp$NDfct  <- factor(tone$xpp$NDind, labels=c('ND1','ND2','ND3','ND4')[ sort(unique(tone$xpp$NDind)) ])    #,ordered=T)
    
  }
  
  
  ## FFR                  ====
  if ( tone$project == 'FFR'){
    tone$xpp$stimNr   <- tone$xpp$strobe1
    tone$xpp$toneNr   <- 1+(tone$xpp$strobe1-1)%/%2
    tone$xpp$polarity <- 1+(tone$xpp$strobe1-1)%%2
    
    tone$xpp$block     <- tone$trialInfo$BlockNumber[tone$xpp$trialNr]
    tone$xpp$condition <- tone$trialInfo$ConditionNumber[tone$xpp$trialNr]
    tone$xpp$context   <- 1 + ((tone$xpp$condition-1)%%3)
    
    tone$xpp$half  <- 1 + (0:(length(tone$xpp$trialNr)-1))%/%2000
    tone$xpp$quart <- 1 + (0:(length(tone$xpp$trialNr)-1))%/%1000
    
    if (tone$dfline$link[1] == 'FFR_16Tone_link.txt'){
      tone$xpp$stimNr   <- tone$trialInfo$ConditionNumber[tone$xpp$trialNr]
      tone$xpp$toneNr   <- 1+(tone$xpp$stimNr-1)%/%4
      tone$xpp$polarity <- 1+(tone$xpp$strobe1-tone$xpp$stimNr)%%2
      
      tone$xpp$speakerID<- 1+(tone$xpp$stimNr-1)%%4
      tone <- subs.data(tone, tone$xpp$strobe2<18)
    }
    if (tone$dfline$link[1] == 'FFR_OA_link.txt'){
      #  tone <- subs.data(tone, tone$xpp$strobe2<8)
    }
  }
  
  ## rewardAnticipation   ====
  if ( tone$project == 'rewardAnticipation'){
    
    tone$xpp$rewardID <- tone$xpp$strobe2
    tone$xpp$toneID   <- tone$xpp$strobe2
    
    tone$xpp$block                     <- tone$trialInfo$trialStrobe3[1:dim(tone$xpp)[1]]   # index 6 - 64
    tone$xpp$toneID[tone$xpp$block==2] <- c(2,1,4,3)[ tone$xpp$rewardID[tone$xpp$block==2] ]
    
  }
  
  ## vocalization         ====
  if ( tone$project == 'vocalizations'){
    
    NSpc   <- tone$dfline$NSpc
    NType  <- tone$dfline$NCall
    NID    <- tone$dfline$NID
    
    tone$xpp$stimulusNr  <- tone$xpp$strobe2
    tone$xpp$callID      <- tone$xpp$strobe3
    tone$xpp$callType    <- tone$xpp$strobe4
    tone$xpp$callSpecies <- tone$xpp$strobe5
    
    #tmpNID <- length(unique(tone$xpp$callID))
    
    tone$xpp$callID[   tone$xpp$callSpecies==(NSpc+1) ] <- NID+1
    tone$xpp$callType[ tone$xpp$callSpecies==(NSpc+1) ] <- NType+1
    
  }
  
  ## chait                ====
  if ( tone$project == 'chait'){
    
    tone$xpp$trialType   <- tone$xpp$trialStrobe1 - 64
    tone$xpp$snpNr       <- tone$xpp$trialStrobe2 
    tone$xpp$BlockNumber <- tone$xpp$trialStrobe3 
  }
  
  ## cyclicNoise2         ====
  if ( tone$project == 'cyclicNoise2'){
    tone$xpp$trialType   <- tone$xpp$trialStrobe1 -64 
    tone$xpp$snpNr       <- tone$xpp$trialStrobe2
    tone$xpp$BlockNumber <- tone$xpp$trialStrobe3
  }
  
  ## continuity           ====
  if ( tone$project == 'continuity'){
    tone$xpp$trialType   <- tone$xpp$strobe2
  }
  ## instrumentalConditioning ====
  if (tone$project == 'instrumentalConditioningOld'){
    
    tone$xpp$active  <- tone$xpp$strobe1 == 4
    tone$xpp$reward  <- tone$xpp$strobe1 == 3
    tone$xpp$condID  <- tone$trialInfo$ConditionNumber[ tone$xpp$trialNr ]
    tone$xpp$blockID <- tone$trialInfo$BlockNumber[ tone$xpp$trialNr ]
    
    tone$xpp$soundID <- tone$xpp$strobe1
    tone$xpp$soundID[tone$xpp$strobe1 == 4] <- 5 + (tone$xpp$condID[tone$xpp$strobe1 == 4]-1)%%4%/%2
    
    # fix the passive sounds:
    switchInd                  <- which(tone$xpp$condID %in% c(3,4,7,8,11,12) & tone$xpp$strobe1%in%c(5,6))
    tone$xpp$soundID[switchInd] <- c(1,2,3,4,6,5)[tone$xpp$strobe1[switchInd]] 
    
  }
  if (tone$project == 'instrumentalConditioning'){
    
    
    # lowRewPos <-  layout of the choice targets
    # choiceDir <-  direction of the choice
    # choiceID  <- identity of the choice
    # soundID   <- identity of the sound played; identical to choiceDir, except in oddball trials 
    
    tone$xpp$blockID    <- tone$trialInfo$BlockNumber[ tone$xpp$trialNr ]
    tone$xpp$trialError <- tone$trialInfo$trialError[ tone$xpp$trialNr ]
    
    
    # there are 6 conditions grouped into 3 types: normal, invert outcome, single choice
    # on the invert outcome trials, the sound is switched and the reward magnitude is switched. not sure about trialError...
    # on single choice trials, the animal only has the small option to choose from
    tone$xpp$condID    <- tone$trialInfo$ConditionNumber[ tone$xpp$trialNr ]
    tone$xpp$trialType <- 1 + (tone$xpp$condID-1)%/%2 
    tone$xpp$lowRewPos <- 1 + (tone$xpp$condID-1)%%2 
    
    tone$trialInfo$lowRewPos <- 1 + (tone$trialInfo$ConditionNumber-1)%%2 
    

    # eventType0 : lever movement and sound presentation
    # eventType1 : reward
    # eventType2 : passive sound or silence
    
    tone$xpp$eventType  <- as.numeric(tone$xpp$strobe1 %in% c(5,9,10,11)) + as.numeric(tone$xpp$strobe1%in%c(9,10,11))
    
    # eventStrobe1 + eventType 0 : identity of the sound played during movement
    # in most cases this is identical to the direction of the movement
    # eventStrobe1 + eventType 2 : identity of the passive sound
    tone$xpp$soundID                       <- tone$xpp$strobe1
    tone$xpp$soundID[tone$xpp$strobe1 > 6] <- tone$xpp$strobe1[tone$xpp$strobe1 > 6] - 8
    

    # direction of the choice 
    tone$xpp$choiceDir <- c(NA)
    tone$xpp$choiceDir[tone$xpp$strobe1 %in% c(1,2)] <- tone$xpp$strobe1[tone$xpp$strobe1 %in% c(1,2)]
    # get choice info on a per trial basis
    ndx <- which( !is.na(tone$xpp$choiceDir ))
    tone$trialInfo$choiceDir <- NA
    tone$trialInfo$choiceDir[ tone$xpp$trialNr[ndx]   ] <- tone$xpp$choiceDir[ndx]
    # bring the choice info back into the xpp frame for all events
    tone$xpp$choiceDir <- tone$trialInfo$choiceDir[ tone$xpp$trialNr]
    
    
    # strobe1 always indicates the movement direction, even on inverted trials
    # invert movement direction for the inverted trials
    #tone$xpp$choiceDir[ tone$xpp$eventType==0 & tone$xpp$trialType == 3] <- 3-tone$xpp$strobe1[tone$xpp$eventType==0 & tone$xpp$trialType == 3]
    
    
    # identity of the choice target
    tone$trialInfo$choiceID <- as.numeric( tone$trialInfo$choiceDir == 3-tone$trialInfo$lowRewPos)
    
    tone$xpp$choiceID <- c(NA)
    tone$xpp$choiceID <- tone$trialInfo$choiceID[ tone$xpp$trialNr]
    
    
    
    # fix the passive sounds:
    #switchInd                  <- which(tone$xpp$condID %in% c(3,4,7,8,11,12) & tone$xpp$strobe1%in%c(5,6))
    #tone$xpp$soundID[switchInd] <- c(1,2,3,4,6,5)[tone$xpp$strobe1[switchInd]] 
    
    tone$xpp$reward  <- tone$xpp$strobe1 == 3
    
    
  }
  ## SSA                                                ===
  if (tone$project == 'SSA'){
    tone$xpp$deviant     <- tone$xpp$strobe2 == 2
    tone$xpp$deviant[which(tone$xpp$deviant)-1] = 2
    tone$xpp$deltaF      <- tone$xpp$strobe4
    
    tone$xpp$condID  <- tone$trialInfo$ConditionNumber[ tone$xpp$trialNr ]
    tone$xpp$blockID <- tone$trialInfo$BlockNumber[ tone$xpp$trialNr ]
    
  }
  
  ## SMART task           ====
  if ( tone$project == 'SMART' & (tone$dfline$link =='SMART_link.txt' | tone$dfline$link =='SMART_amwn_link.txt'| tone$dfline$link =='SMART_sweep_link.txt') ){
    
    
    catID <- c(tone$dfline$minF, tone$dfline$deltaF)
    
    # remove and ignore trials that were not initiated (trialError 3 and 4)
    # this is important when looking at information about previous trials. 
    # For this purpose, non-initiated trials can be ignored.
    valNdx                        <- which(tone$trialInfo$trialError%in%c(0,1,2,5,6))

    ndxmap         <- rep(NA, dim(tone$trialInfo)[1])
    ndxmap[valNdx] <- 1:length(valNdx)
  
    tone$trialInfo                <- subset(tone$trialInfo, tone$trialInfo$trialError%in%c(0,1,2,5,6))
    
    #ndxmap <- 1:dim(tone$trialInfo)
    #names(ndxmap) <- tone$trialInfo$trialNr
    
    # COND: direction of the visual target
    tone$trialInfo$COND            <- 1 + (tone$trialInfo$ConditionNumber-1)%%2
    tone$trialInfo$dirID           <- 1 + (tone$trialInfo$ConditionNumber-1)%%2
    tone$trialInfo$congruent       <- tone$trialInfo$BlockNumber %in% c(1,3,4) # block 2 is the incongruent block
    
    tone$trialInfo$trialErrorLag1  <- my.lag(tone$trialInfo$trialError,1, circular=F)
    tone$trialInfo$trialErrorLag2  <- my.lag(tone$trialInfo$trialError,2, circular=F)
    tone$trialInfo$trialErrorLead1 <- my.lag(tone$trialInfo$trialError,-1, circular=F)
    
    tone$trialInfo$congruentLag1   <- my.lag(tone$trialInfo$congruent,1,  circular=F)
    tone$trialInfo$congruentLag2   <- my.lag(tone$trialInfo$congruent,2,  circular=F)
    
    tone$trialInfo$congruentLead1  <- my.lag(tone$trialInfo$congruent,-1, circular=F)
    
    
    tone$trialInfo$BlockNumberLag1   <- my.lag(tone$trialInfo$BlockNumber,1,  circular=F)
    tone$trialInfo$BlockNumberLag2   <- my.lag(tone$trialInfo$BlockNumber,2,  circular=F)
    tone$trialInfo$BlockNumberLead1  <- my.lag(tone$trialInfo$BlockNumber,-1, circular=F)
    
    
    tone$trialInfo$correct                                                  <- tone$trialInfo$trialError == 0
    tone$trialInfo$correct[which( !tone$trialInfo$trialError %in% c(0,6) )] <- NA
    
    
    tone$trialInfo$targetLocation <- (tone$trialInfo$ConditionNumber-1)%%2
    
    tone$xpp$trialNrNew          <- ndxmap[tone$xpp$trialNr]
    tone$xpp$RT                  <- tone$trialInfo$RT[tone$xpp$trialNrNew]
    tone$xpp$trialError          <- tone$trialInfo$trialError[tone$xpp$trialNrNew]
    tone$xpp$trialErrorLag1      <- tone$trialInfo$trialErrorLag1[tone$xpp$trialNrNew]
    tone$xpp$trialErrorLag2      <- tone$trialInfo$trialErrorLag2[tone$xpp$trialNrNew]
    tone$xpp$trialErrorLead1     <- tone$trialInfo$trialErrorLead1[tone$xpp$trialNrNew]
    tone$xpp$targetLocation      <- tone$trialInfo$targetLocation[tone$xpp$trialNrNew]
    
    
    tone$xpp$leverReleaseTime    <- tone$trialInfo$trialEndTime[tone$xpp$trialNrNew]
    tone$xpp$blockID             <- tone$trialInfo$BlockNumber[tone$xpp$trialNrNew]
    
    tone$xpp$blockIDLag1             <- tone$trialInfo$BlockNumberLag1[tone$xpp$trialNrNew]
    tone$xpp$blockIDLag2             <- tone$trialInfo$BlockNumberLag2[tone$xpp$trialNrNew]
    tone$xpp$blockIDLead1            <- tone$trialInfo$BlockNumberLead1[tone$xpp$trialNrNew]
    
    tone$xpp$condID              <- tone$trialInfo$ConditionNumber[tone$xpp$trialNrNew]
    tone$xpp$COND                <- tone$trialInfo$COND[tone$xpp$trialNrNew]
    tone$xpp$dirID               <- tone$trialInfo$dirID[tone$xpp$trialNrNew]
    
    tone$xpp$congruent           <- tone$trialInfo$congruent[tone$xpp$trialNrNew]
    tone$xpp$congruentLag1       <- tone$trialInfo$congruentLag1[tone$xpp$trialNrNew]
    tone$xpp$congruentLag2       <- tone$trialInfo$congruentLag2[tone$xpp$trialNrNew]
    
    tone$xpp$congruentLead1      <- tone$trialInfo$congruentLead1[tone$xpp$trialNrNew]
    
    tone$xpp$correct             <- tone$trialInfo$correct[tone$xpp$trialNrNew]
    
    
    tone$xpp$deltaRT            <- tone$xpp$leverReleaseTime - tone$xpp$e.MLonset
    
    
    
    #tone$xpp$correct                                            <- tone$xpp$trialError == 0
    #tone$xpp$correct[which( !tone$xpp$trialError %in% c(0,6) )] <- NA
    
    #tone$xpp$congruent <- tone$xpp$COND %in% c(1,3)
    
    tone$xpp$seqID       <- tone$xpp$strobe2   
    tone$xpp$toneID      <- tone$xpp$strobe3 - 4
    tone$xpp$catID       <- tone$xpp$strobe4 - 4
    #browser()
    
    # based on the tone identity. gives inverted results for inconguent trials
    #tone$xpp$choiceDir                                       <- tone$xpp$catID == catID[2] # 0 is left; 1 is right
    #tone$xpp$choiceDir[ tone$xpp$trialError%in%c(1,2,3,4,5)] <- NA # set to NA if no response was given
    #tone$xpp$choiceDir[tone$xpp$trialError%in%c(6)]          <- tone$xpp$catID[tone$xpp$trialError%in%c(6)] == catID[1]# invert for errors
    
    tone$xpp$choiceDir                                       <- tone$xpp$targetLocation # 0 is left; 1 is right
    tone$xpp$choiceDir[ tone$xpp$trialError%in%c(1,2,3,4,5)] <- NA # set to NA if no response was given
    tone$xpp$choiceDir[tone$xpp$trialError%in%c(6)]          <- !tone$xpp$choiceDir[tone$xpp$trialError%in%c(6)]# invert for errors
    
    
  }
  
  if ( tone$project == 'SMART' & tone$dfline$link =='irfbat_SSA_all_link.txt'){
    tone$xpp$devID  <- tone$xpp$strobe2   
    tone$xpp$toneID <- tone$xpp$strobe3 - 10
    tone$xpp$catID  <- tone$xpp$strobe4 - 4
  }
  ## return               ====
  gc()
  return(tone)
}


getSC96pos <- function(){
  
  sc96pos                              <- read.table('~/Dropbox/eCog/jessesym/positionList_ISI_SC96.txt', header=T)
  getElectrodeStr                      <- function(ex){ tmp <- paste('00',sc96pos$Electrode[ex],sep=''); return(substr(tmp,nchar(tmp)-1,nchar(tmp))) }
  sc96pos$ElectrodeStr                 <- sapply(1:dim(sc96pos)[1], getElectrodeStr)
  sc96pos$bank                         <- c('B','C')[1+as.numeric(sc96pos$Number>64)] 
  sc96pos$Intan                        <- paste( sc96pos$bank,sc96pos$ElectrodeStr,sep='')
  sc96pos$Electrode                    <- sc96pos$Electrode + 1
  sc96pos$Electrode[sc96pos$Number>64] <- sc96pos$Electrode[sc96pos$Number>64] + 64
  sc96pos$order                        <- sort(sc96pos$Intan, index.return=T)$ix
  
  sc96pos$xpos     <- sc96pos$ml1
  sc96pos$ypos     <- -1*sc96pos$ap1
  sc96pos$zpos     <- c(0)
  
  sc96pos <- sc96pos[sc96pos$order,]
  
  return(sc96pos)
}


getSC224pos <- function(animal='Cirque',implant='SC224',sc96ExtnNr='_Intan'){
  
  masterXLS <- paste('~/Dropbox/admin/',animal,implant,'/',animal,implant,'_masterSheet',sc96ExtnNr,'.xlsx',sep='')
  ##Something about my R version on Tinyone doesn't allow reading xlsx, had to change a copy to xls format
  if(Sys.info()[["user"]]=="spencer"){
    #This directory doesn't work any more with the Dropbox permissions bug -SK 20220224
    #masterXLS <- paste('~/Dropbox/Spencer/',animal,implant,'_masterSheet',sc96ExtnNr,'.xls',sep='')
    masterXLS <- paste('/Volumes/Drobo5D3/Transfers/admin/',animal,implant,'/',animal,implant,'_masterSheet',sc96ExtnNr,'.xls',sep='')
  }
  dfImp        <- read.xls(masterXLS, sheet='Impedance', header=T)
  sc96pos      <- dfImp[,c('Electrode','Intan')]
  sc96pos$Number <- sc96pos$Electrode
  sc96pos$Electrode <- 1:dim(sc96pos)[1]
  sc96pos$xpos <- dfImp$Column
  sc96pos$ypos <- -1*dfImp$Row
  sc96pos$zpos <- c(0)
  #sc96pos                              <- read.table('~/Dropbox/eCog/cirquesym/positionList_ISI_SC224.txt', header=T)
  getElectrodeStr                      <- function(ex){ tmp <- paste('00',sc96pos$Electrode[ex],sep=''); return(substr(tmp,nchar(tmp)-1,nchar(tmp))) }
  sc96pos$ElectrodeStr                 <- sapply(1:dim(sc96pos)[1], getElectrodeStr)
  
  #sc96pos$order                        <- sort(sc96pos$Intan, index.return=T)$ix
  
  #sc96pos$xpos     <- sc96pos$ml1
  #sc96pos$ypos     <- -1*sc96pos$ap1
  #sc96pos$zpos     <- c(0)
  
  #sc96pos <- sc96pos[sc96pos$order,]
  
  return(sc96pos)
}


getEEGpos <- function(eCogFile,xscale=.5,yscale=1,flat=T){
  
  if (eCogFile == 'Ben')
    eCogFile <- 'benjaminsym'
  
  ep             <- read.table(paste('~/Dropbox/eCog/',eCogFile,'/positionList_ISI_eCog.txt',sep=''),header=T)
  srt <- sort(ep$Number, index.return=T)$ix
  ep <- ep[srt,]
  
  xscale <- .5
  yscale <- 1
  ep$xpos <- ep$ml1 * xscale
  ep$ypos <- ep$ap1 * yscale
  ep$zpos <- ep$dv1
  if (flat){
    ep$zpos <- ep$zpos - max(ep$zpos)
    ep$xpos <- ep$xpos - sign(ep$xpos) * 0.9 * ep$zpos * ifelse(abs(ep$xpos)>.1,1,0) # * abs(ep$xpos)/max(abs(ep$xpos))
    ep$ypos <- ep$ypos - sign(ep$ypos) * 0.0 * ep$zpos
  }
  return(ep)
}


getMePhyspos <- function(chStr, type='exact', measured=F, ref='skull', layout='tight', returnAll=F){
  
  # type:    'exact' uses the exact depth, '' uses the index as depth
  # layout:  '3rows': 3 rows with equal space for each slice no matter how many electrodes are in that slice
  #          'tight': one row with one space between each slice  
  # measured: default is FALSE; 
  #           if TRUE, returns the locations of the contacts based on the post implantation CT and the spline interpolation. If electrodes have not been implanted, will return NA for ml1, ap1, and dv1
  # ref:      reference frame. default is skull
  #           ref=='mesoscope' will transform ap1, ml1 and dv1 to mesoscope space before calculating electrode layout
  
  # this function provides the x y and z-coordinates of each contact in the mesoscope based on the chStr
  # Step 1:
  # read xls file and calculate position of all conta
  # needs ml1 ap1 dv1 [3D coordinates of electrodes] as well as  xpos ypos zpos [flattened coordinates], and chStr
  
  # Note that this function returns position in skull space. 
  # Note that this function uses the requested length of the blank space, not the actual length
  
  ## EEG electrode position from Ben unless specified                           ====
  eegChannelPosition               <- getEEGpos('ben2022')
  if (measured)
    eegChannelPosition               <- getEEGpos('soleil2022')
  
  eegChannelPosition$chStr         <- paste('ch',eegChannelPosition$Number,sep='') 
  eegChannelPosition$xpos          <- 0.9 * eegChannelPosition$xpos 
  eegChannelPosition$ypos          <- 1.5 * eegChannelPosition$ypos + 1
  eegChannelPosition$Yx            <- 0
  eegChannelPosition$Xx            <- 0
  eegChannelPosition$Zx            <- 0
  eegChannelPosition$Xgrid         <- 0
  eegChannelPosition$Ygrid         <- 0
  
  
  # get channel Number from channel String
  chNr <- sapply(chStr, function(cs){ as.numeric(strsplit(cs,'ch')[[1]][2] ) })
  
  chNrtmp <- chNr-34
  chNrtmp <- chNrtmp[chNrtmp>=0]
  chNrtmp <- 3 + floor(chNrtmp/16)
  
  
  # ## load the excel file with the info about all contacts                       ====
  # folderName <- '~/Dropbox/VolumeEEG/Solei/Plexon/'
  # fileName <- 'ShaftLength_ForPlexon_newLayout_FINAL.xlsx'
  # 
  # old <- FALSE
  # if (old){
  #   # read info about each shaft:
  #   plxtmp       <- read.xls( paste(folderName,fileName,sep=''), sheet='forPlexon', header=T, blank.lines.skip=1 )
  #   plxtmp$ml1   <- plxtmp$ml1/10
  #   plxtmp$ap1   <- plxtmp$ap1/10
  #   
  #   # this is probably what they were like:
  #   plxtmp$Xx    <- sapply( 1:dim(plxtmp)[1], function(ix){ as.numeric(substr(plxtmp$Grid.Label[ix],2,2) ) } ) 
  #   plxtmp$Yx    <- sapply( 1:dim(plxtmp)[1], function(ix){ as.numeric(substr(plxtmp$Grid.Label[ix],5,nchar(as.character(plxtmp$Grid.Label[ix]))) ) } ) 
  #   
  #   
  #   fieldNames   <-  c('ShaftNr','ml1','ap1','top','tip','b','a','c','s','Xx','Yx')  
  #   
  #   #sort from left to righ and back to front
  #   plx <- plxtmp[ sort(plxtmp$ap+plxtmp$ml/20, index.return=T)$ix, fieldNames]
  #   old <- plx
  # }# the old routine using read.xls
  # 
  # # new routine using read.xlsx2
  # colClasses <- c( rep('numeric',13), rep('character',3), rep('numeric',7))
  # plxtmp     <- read.xlsx2( paste(folderName,fileName,sep=''), sheetName='forPlexon', header=T, startRow = 2, colClasses = colClasses, stringsAsFactors=FALSE)
  # plxtmp     <- subset(plxtmp, !is.nan(plxtmp$ShaftNr))
  # 
  # plxtmp$ml1   <- plxtmp$ml1/10
  # plxtmp$ap1   <- plxtmp$ap1/10
  # 
  # # not sure what happened here. These were the lines on 06/18 but they can't be correct. I commented them out and tried fixing them below
  # ##plxtmp$Xx    <- sapply( 1:dim(plxtmp)[1], function(ix){ as.numeric(substr(plxtmp$Grid.Label[ix],2,2) ) } ) 
  # #plxtmp$Y    <- sapply( 1:dim(plxtmp)[1], function(ix){ as.numeric(substr(plxtmp$Grid.Label[ix],5,nchar(as.character(plxtmp$Grid.Label[ix]))) ) } ) 
  # 
  # # this is probably what they were like:
  # plxtmp$Xx    <- sapply( 1:dim(plxtmp)[1], function(ix){ as.numeric(substr(plxtmp$Grid.Label[ix],2,2) ) } ) 
  # plxtmp$Yx    <- sapply( 1:dim(plxtmp)[1], function(ix){ as.numeric(substr(plxtmp$Grid.Label[ix],5,nchar(as.character(plxtmp$Grid.Label[ix]))) ) } ) 
  # 
  # 
  # fieldNames   <-  c('ShaftNr','ml1','ap1','top','tip','b','a','c','s','Xx','Yx')  
  # 
  # #sort from left to righ and back to front
  # plx <- plxtmp[ sort(plxtmp$ap+plxtmp$ml/20, index.return=T)$ix, fieldNames]
  
  
  ## load the excel file with the info about all contacts 
  folderName <- '~/Dropbox/VolumeEEG/Solei/Plexon/'
  fileName <- 'ShaftLength_ForPlexon_newLayout_FINAL.xlsx'
  
  old <- FALSE
  
  if (old) {
    
    # Original routine using read.xls
    plxtmp <- read.xls(paste(folderName, fileName, sep=''), sheet='forPlexon', header=TRUE, blank.lines.skip=1)
    plxtmp$ml1 <- plxtmp$ml1 / 10
    plxtmp$ap1 <- plxtmp$ap1 / 10
    
    plxtmp$Xx <- sapply(1:dim(plxtmp)[1], function(ix) { as.numeric(substr(plxtmp$Grid.Label[ix], 2, 2)) }) 
    plxtmp$Yx <- sapply(1:dim(plxtmp)[1], function(ix) { as.numeric(substr(plxtmp$Grid.Label[ix], 5, nchar(as.character(plxtmp$Grid.Label[ix])))) }) 
    
    fieldNames <- c('ShaftNr', 'ml1', 'ap1', 'top', 'tip', 'b', 'a', 'c', 's', 'Xx', 'Yx')  
    
    plx <- plxtmp[sort(plxtmp$ap + plxtmp$ml / 20, index.return=TRUE)$ix, fieldNames]
    old <- plx
  } 
  
  # New routine using read.xlsx2 (if available)
  if (exists("read.xlsx2")) {
    colClasses <- c(rep('numeric', 13), rep('character', 3), rep('numeric', 7))
    plxtmp <- read.xlsx2(paste(folderName, fileName, sep=''), sheetName='forPlexon', header=TRUE, startRow=2, colClasses=colClasses, stringsAsFactors=FALSE)
  
    } else {
    
    # Alternative using openxlsx when read.xlsx2 is not avaliable
    library(openxlsx)
    
    plxtmp <- openxlsx::read.xlsx(paste(folderName, fileName, sep=''), sheet='forPlexon', colNames=TRUE, startRow=2)
    
    # Convert columns manually to match colClasses because openxlsx doesn't do that automatically
    colClasses <- c(rep('numeric', 13), rep('character', 3), rep('numeric', 7))
    for (i in seq_along(colClasses)) {
      if (i <= ncol(plxtmp)) {
        if (colClasses[i] == "numeric") {
          plxtmp[[i]] <- suppressWarnings(as.numeric(plxtmp[[i]]))
        } else if (colClasses[i] == "character") {
          plxtmp[[i]] <- as.character(plxtmp[[i]])
        }
      }
    }
  }
  
  plxtmp <- subset(plxtmp, !is.nan(plxtmp$ShaftNr))
  
  plxtmp$ml1 <- plxtmp$ml1 / 10
  plxtmp$ap1 <- plxtmp$ap1 / 10
  
  plxtmp$Xx <- sapply(1:dim(plxtmp)[1], function(ix) { as.numeric(substr(plxtmp$Grid.Label[ix], 2, 2)) }) 
  plxtmp$Yx <- sapply(1:dim(plxtmp)[1], function(ix) { as.numeric(substr(plxtmp$Grid.Label[ix], 5, nchar(as.character(plxtmp$Grid.Label[ix])))) }) 
  
  fieldNames <- c('ShaftNr', 'ml1', 'ap1', 'top', 'tip', 'b', 'a', 'c', 's', 'Xx', 'Yx')  
  
  plx <- plxtmp[sort(plxtmp$ap + plxtmp$ml / 20, index.return=TRUE)$ix, fieldNames]
  
  plx$planned <- 
    c(0,0,           # posterior
      0,1,1,0,  
      0,1,1,1,1,  
      0,1,1,1,1,1,  
      0,1,1,1,1,1,  
      0,1,1,1,1,1,  
      0,1,1,1,1,1,1,  
      0,1,1,1,1,1,1,  
      1,1,1,1,1,1,
      1,1,1,1,1,1,
      1,1,1,1,1,1,
      1,1,1,1,1,
      1,1,1,1,
      1,1)            #anterior
  
  plx$use <- c(0,0,           # posterior
               0,2,2,0,  
               0,1,1,1,1,  
               0,1,1,1,1,1,  
               0,1,1,1,1,1,  
               0,1,1,1,1,2,  
               0,1,1,1,1,1,1,  
               0,9,1,1,1,2,1,  
               9,1,1,1,1,1,
               9,1,1,1,1,1,
               9,1,1,1,1,1,
               9,1,1,1,2,         # is X5_Y12 broken or 
               1,1,1,2,           # is X4_Y13 broken, or both???
               2,1)               #anterior
  
  #plx$oNr <- c(NA,NA,           # posterior
  #             NA, 3, 4,NA,  
  #             NA, 5, 6, 7, 8,  
  #             NA, 9,10,11,12,16,  
  #             NA,13,14,15,19,20,  
  #             NA,17,18,22,23,24,  
  #             NA,21,28,27,26,25,32,  
  #             NA,29,30,31,34,35,36,  
  #             33,37,38,39,40,44,
  #             41,42,43,46,47,48,
  #             45,49,50,51,52,56,
  #             53,54,55,59,60,
  #             57,58,63,64,
  #             61,62)            #anterior
  
  
  # with channel order inverted for all left-hand boards. Not sure why, but this is the right way.
  plx$oNr <- c(NA,NA,           # posterior
               NA, 4, 3,NA,  
               NA, 5, 6, 7, 8,  
               NA,12,11,10, 9,16,  
               NA,13,14,15,18,17,  
               NA,20,19,22,23,24,  
               NA,21,28,27,26,25,32,
               NA,29,30,31,35,34,33,  
               36,37,38,39,40,41,
               44,43,42,46,47,48,
               45,52,51,50,49,56,
               53,54,55,58,57,
               60,59,63,64,
               61,62)            #anterior
  
  # Same as above, but accounting for a misplugging: X2 Y11 was plugged into board 6 left first channel instead of board 7 left first 
  # channel/ Swapped 52 with 44
  plx$oNr <- c(NA,NA,           # posterior
               NA, 4, 3,NA,  
               NA, 5, 6, 7, 8,  
               NA,12,11,10, 9,16,  
               NA,13,14,15,18,17,  
               NA,20,19,22,23,24,  
               NA,21,28,27,26,25,32,
               NA,29,30,31,35,34,33,  
               36,37,38,39,40,41,
               52,43,42,46,47,48,
               45,44,51,50,49,56,
               53,54,55,58,57,
               60,59,63,64,
               61,62)            #anterior
  
  
  
  #plx$side <- c(NA,NA,           # posterior
  #             NA,1,1,NA,  
  #             NA,2,2,2,2,  
  #             NA,1,1,1,1,2,  
  #             NA,2,2,2,1,1,  
  #             NA,1,1,2,2,2,  
  #             NA,2,1,1,1,1,2,  
  #             NA,2,2,2,1,1,1,  
  #             1, 2,2,2,2,1,
  #             1, 1,1,2,2,2,
  #             2, 1,1,1,1,2,
  #             2, 2,2,1,1,
  #             1, 1,2,2,
  #             2, 2)            #anterior
  
  
  tplx <- subset(plx, plx$planned==1 & plx$oNr %in% chNrtmp)
  
  #plot(eegChannelPosition$ml1, eegChannelPosition$ap1, col='red')
  #points(tplx$ml1, tplx$ap1-1)
  
  
  ## create the xposition of the display                                        ====
  if (layout == 'tight'){
    tplx$xpos <- -7.00 - 1:dim(tplx)[1] - cumsum(as.numeric(diff(c(tplx$ap1[1],tplx$ap1))>0))
    tplx$addToYpos <- 0
  }
  
  if (layout == '3rows'){# Xx is the ml index; Yx is the ap index
    slicePerRow    <- 3
    tplx$xpos      <- 4 + -1*(tplx$Xx+1) - 8*(tplx$Yx-2)%%slicePerRow
    tplx$addToYpos <- 12 * (tplx$Yx-2)%/%slicePerRow
  }
  
  if (layout == '4rows'){# Xx is the ml index; Yx is the ap index
    slicePerRow    <- 4
    tplx$xpos      <- 4 + -1*(tplx$Xx+1) - 8*(tplx$Yx-2)%%slicePerRow
    tplx$addToYpos <- 12 * (tplx$Yx-2)%/%slicePerRow
  }
  
  if (layout == '5rows'){# Xx is the ml index; Yx is the ap index
    
    Yx_shft       <-   0
    Xgrid_scl     <-  -8 # how much space for each slice in ML
    Ygrid_scl     <-  12 # how much space for each slice in DV
    
    slicePerRow    <- 5
    #tplx$xpos      <- 4 + -1*(tplx$Xx+1) - 8*(tplx$Yxtmp-2)%%slicePerRow
    #tplx$addToYpos <- 12 * (tplx$Yxtmp-2)%/%slicePerRow
    
    tplx$Xgrid <- Xgrid_scl * (tplx$Yx - Yx_shft)%%slicePerRow
    tplx$Ygrid <- Ygrid_scl * (tplx$Yx - Yx_shft)%/%slicePerRow
    
    tplx$xpos      <- 4 + 1*(tplx$Xx+1) + tplx$Xgrid
    tplx$addToYpos <-                      tplx$Ygrid
    
  }
  
  
  
  ## create the data frame for each contact                                     ====
  cntct    <- tplx[ rep(1:dim(tplx)[1],each=16),   ]
  cntct$Zx <- rep(1:16,dim(tplx)[1])
  
  #cntct$Ex <- 1:dim(cntct)[1]
  cntct$Ex    <- 33 + 16 * (cntct$oNr-3) + cntct$Zx
  cntct$chStr <- paste('ch',cntct$Ex,sep='')
  
  
  ## calculate actual depth in cm of each   
  cntct$dv1 <-(45 + 2.5 - cntct$b - (cntct$Zx-1)*cntct$s/1000 )/10
  
  ## Load measured electrode locations if requested                             ====
  if (measured == TRUE){
    ## load the estimated electrode locations from CT and spline fit
    xlsDir     <- '~/Dropbox/VolumeEEG/Solei/plexon/'
    rdaPosFile <- 'actualShaftLocations.rda'
    load( paste(xlsDir,rdaPosFile,sep='') )
    
    # link acutalElectrodePositions in ECs to estimated electrode positions in cntct
    findMatchingIndx <- function(ix){
      mx <- NA
      tmp <- which( ECs$Xx == cntct$Xx[ix] & ECs$Yx == cntct$Yx[ix] & ECs$Zx == cntct$Zx[ix])
      if (length(tmp)==1)
        mx <- tmp
      
      return(mx)
    }
    mx <- sapply(1:dim(cntct)[1], findMatchingIndx )
    
    ## debug cntc vs ECs position [ECs is in mesoscope space and in mm, cntct is in skull space and cm]
    mes2skl <- read.table( '~/Dropbox/data/monkey.highres/stereotax_coordinate_frame/solei/soleil2022/mesoscope2soleil_skull.mat')
    
    A_mes    <- as.matrix(ECs[,1:3])
    A_skl    <- applyFSLtransform( A_mes, mes2skl)
    
    if (identical(ref, 'skull') ){
      cntct[!is.na(mx),c('ml1','ap1','dv1')] <- A_skl/10
      cntct[ is.na(mx),c('ml1','ap1','dv1')] <- NA
    }
    
    if (identical(ref, 'mesoscope') ){
      cntct[!is.na(mx),c('ml1','ap1','dv1')] <- A_mes/10
      cntct[ is.na(mx),c('ml1','ap1','dv1')] <- NA
    }
      
    
    debug <- FALSE
    if (debug){
      browser()
      #A  <- array( c(0,0,44.5), c(1,3) ) # example point in cm relative to zero
      #A  <- array( rep(c(0,0,44.5),each=2), c(2,3) ) # example point in cm relative to zero
      tst   <- applyFSLtransform( A, mes2skl)
      
      p1 <- c( 0, 0,   44.5) # origin
      p2 <- c(20, 9.5, 44.5)
      p3 <- c(0, 34.5, 44.5)
      p4 <- c(0, -30,  45.5)
      
      t1 <- c( 2, 12.5,  53.5)
      t2 <- c(22, 22,    53)
      t3 <- c(2,  47,    53)
      t4 <- c(2, -17.5,  54.5)
      
      #mes2skl[1,4] <- 0
      #skl2mes[1,4] <- 0  
      applyFSLtransform( rbind(p1,p2,p3, p4), mes2skl) - rbind(t1,t2,t3,t4)
      applyFSLtransform( rbind(t1,t2,t3, t4), skl2mes) - rbind(p1,p2,p3, p4)
      
      # there is a large discrepancy in depth for the most anterior contacts.
      # check out slice 13 for an extreme example.
      # 
      for (tx in 2:14){# coronal slices
        
        tx1 <- which(cntct$Yx==tx)
        plot(cntct$ml1[tx1], cntct$dv1[tx1], xlim=c(0,3.5), ylim=c(-2,4), main = paste('Slice #:', tx))
        
        tx2 <- which(ECs$Yx==tx)
        #points(ECs$X[tx2]/10, ECs$Z[tx2]/10, col='red')
        points(tst[tx2,1]/10, tst[tx2,3]/10, col='green')
      }  
      for (tx in 1:7){# sagital slices
        
        tx1 <- which(cntct$Xx==tx)
        plot(cntct$ap1[tx1], cntct$dv1[tx1], xlim=c(-2,5), ylim=c(-2,4), main = paste('Slice #:', tx))
        
        tx2 <- which(ECs$Xx==tx)
        #points(ECs$X[tx2]/10, ECs$Z[tx2]/10, col='red')
        points(tst[tx2,2]/10, tst[tx2,3]/10, col='green')
      }  
    }
    
  }
  
  ## convert to mesoscope coordinates if reqested                               ====
  if (measured == FALSE & identical(ref, 'mesoscope')){
    skl2mes <- read.table( '~/Dropbox/data/monkey.highres/stereotax_coordinate_frame/solei/soleil2022/soleil_skull2mesoscope.mat')
    
    A_skl    <- as.matrix( cntct[,c('ml1','ap1','dv1')] )
    A_mes    <- applyFSLtransform( 10 * A_skl, skl2mes)
    
    cntct[,c('ml1','ap1','dv1')] <- A_mes/10
  }
  
  
  ## add depth information                                                      ====
  # use actual depth for display 
  cntct$ypos <- cntct$addToYpos + 2.5*cntct$dv1
  
  if (type == 'index')
    cntct$ypos <- cntct$addToYpos + 6.5 - (cntct$Zx-1)*.75
  
  cntct             <- cntct[ sort(cntct$Ex,index.return=T)$ix ,]
  #  plot(cntct$Xx[cntct$Yx==3], cntct$dv[cntct$Yx==3] )
  
  
  ## combine eeg and TProbe data                                                ====
  valFieldStr <- c('ml1','ap1','dv1','xpos','ypos','chStr','use','Yx','Xx','Zx','Xgrid','Ygrid')
  eegChannelPosition$use <- c(1)
  res                    <- rbind( eegChannelPosition[,valFieldStr],cntct[,valFieldStr] )
  
  #print(dim(res))
  
  # create an index so lines from cntct can be called by chStr
  tmp <- 1:dim(res)[1]; names(tmp) <- res$chStr
  
  # subset chStr to include only channels that have use==1
  #chStr <- chStr[ res$use==1 ]
  
  tres <- res[ tmp[chStr], valFieldStr] 
  
  if (!returnAll)
    tres <- tres[tres$use==1,]
  
  tres$slice <- as.factor(tres$Yx)
  
  ## return the results                                                         ===== 
  tres$Xgrid <- as.numeric(tres$Xgrid)
  tres$Ygrid <- as.numeric(tres$Ygrid)
  
  return( tres )
}


getMePhysAnatomy <- function( chPos, ref='skull' ){
  
  chPos <- chPos[ which(chPos$Yx > 0),]
  # load the corresponding nifti anatomy file (only D99 works so far), determine the relevant slices, and extract gray-matter outlines using contourLines
  GMClist <- NULL
  
  ovl     <- 'D99'
  imgDir  <- '~/Dropbox/data/monkey.highres/stereotax_coordinate_frame/solei/soleil2022/tmp/'
  imgName <- paste(ovl, '_t1mask_', ref, '.nii.gz',sep='')
  img     <- read.nii( paste(imgDir, imgName, sep='') )
  
  # flip left-right to get to a right-right convention
  img$data <-  img$data[ seq(dim(img$data)[1],to=1,by=-1) ,,,,drop=F] 
  
  
  #ctName  <- paste('soleil', '_ct1_', ref, '.nii.gz',sep='')
  #ct     <- read.nii( paste(imgDir, ctName, sep='') )
  if (identical(ovl, 'D99') ){
    zlim <- c( 1,180)
    img$data <- img$data%%256
  }
  
  xbx <- c(55:185)
  ybx <- c(45:220)
  zbx <- c(20:145)
  
  ## prepare Overlay by computing slices with electrodes ==
  zero                <- c(120,100,50)
  # dimensions of the images in cm
  xax <- ((0:(img$header$dim[2]-1))-zero[1])*img$header$pixdim[2]/10
  yax <- ((0:(img$header$dim[3]-1))-zero[2])*img$header$pixdim[2]/10
  zax <- ((0:(img$header$dim[4]-1))-zero[3])*img$header$pixdim[2]/10
  
  # find the closest slice index for each coronal slice ===
  mnYslc   <- round(20*tapply(chPos$ap1,chPos$Yx, na.mean))/20
  Yx_list  <- as.numeric(names(mnYslc))
  SliceNdx <- which(yax %in% mnYslc)
  
  print(Yx_list)
  print( yax[SliceNdx] )
  
  # 
  getGrayMatterContours <- function(Yx){
    
    Yxi <- which(Yx_list==Yx)
    slx  <- SliceNdx[Yxi]
    cpx <- which(chPos$Yx==Yx)
    
    slc <- img$data[,slx,,1]
    slc[which(xax<0),] <- 0
    CLs <- contourLines( xax[xbx], zax[zbx], slc[xbx,zbx], levels=c(70,130) ) 
    
    getLength <- function(cx){ length(CLs[[cx]]$x) }
    CLl <- sapply(1:length(CLs), getLength)
    return(CLs[ CLl > 20])
  }
  GMClist               <- sapply( sort(unique(chPos$Yx)), getGrayMatterContours )
  names(GMClist)        <- sort(unique(chPos$Yx))
  
  
  ## this is the option that loads leo's old gray-matter outlines
  #load(paste(baseDir,'eCog/leo_20130122/brainOutline.rda',sep=''))
  #bOnrm <- brainOutline
  #for (cx in 1:length(brainOutline)){
  #  bOnrm[[cx]]$x <- brainOutline[[cx]]$x/10
  #  bOnrm[[cx]]$y <- brainOutline[[cx]]$y/10 - 2 # to match the downward shift in showSignalProcessorMap
  #}
  #GMClist[[1]] <- bOnrm
  
  ## this is the option that uses gray-matter contours through 
  slc <- apply(img$data[,,76:120,1]>70,c(1,2), sum)
  slc[which(abs(xax)<0.07),] <- 0
  CLs <- contourLines( xax[xbx], yax[ybx]-2, slc[xbx,ybx], levels=c(2.5) ) 
  getLength <- function(cx){ length(CLs[[cx]]$x) }
  CLl <- sapply(1:length(CLs), getLength)
  
  
  ## read cusomized sulci from brainOutline directory
  sulciFileNameList <- system( paste('ls ', '~/Dropbox/eCog/soleil2022/brainOutline/*.txt' ), intern = T ) 
  readSulci <- function(fx, mlt=1){
    xyz <- read.table(sulciFileNameList[fx], header=T)
    return(list(x=mlt*xyz$x/10, y=xyz$y/10 - 2) )
  }
  sulciList  <- lapply(1:length(sulciFileNameList),readSulci)
  sulciListR <- lapply(1:length(sulciFileNameList),readSulci, mlt=-1)
  
  tmp <- CLs[CLl>200]
  tmp[ (length(tmp)+1) : (length(tmp)+length(sulciList))] <- sulciList
  tmp[ (length(tmp)+1) : (length(tmp)+length(sulciList))] <- sulciListR
  
  
  # Since slice 2 is not being presented, swap out 
  resList                        <- list( tmp )
  #resList[1]                     <- tmp
  resList[2:(length(GMClist)+1)] <- GMClist
  names(resList) <- c("0", names(GMClist))
  
  return(resList)
}


baselineCorrectionSignalProcessor <- function(tone, blrng=c(-50,0), deVar=F){
  tmp            <- baseLineCorrection(tone$lfp,blrng, tone$taxis, sdnrm=deVar)  
  tone$lfp       <- tmp$res
  tone$xpp$blcor <- tmp$correction
  return(tone)
}


reRefSignalProcessor <- function(tone, reRefChan=1, refActImport=NULL,dim=3){
  
  # tone$lfp 3dimensional vector with lfps [trial, time, chNr]
  # reRefChan: list of numbers of channels to be used as reference instead 
  
  if (dim==1)
    tone$lfp <- aperm(tone$lfp, c(3,2,1))
  
  if (is.null(refActImport)){
    if (length(reRefChan)==1)
      refAct   <- tone$lfp[,,reRefChan]
    
    if (length(reRefChan)>1)
      refAct <- apply(tone$lfp[,,reRefChan, drop=F],c(1,2),na.mean)
    
    if (dim(tone$lfp)[1]==1)# dropped dimension fix
      dim(refAct) <- c(dim(tone$lfp)[1:2])
  }
  
  if (!is.null(refActImport)){
    refAct <- refActImport
  }
  
  if ( !identical( dim(refAct), dim(tone$lfp)[c(1,2)]) )
    browser()
  
  
  rerefLFP  <- tone$lfp
  
  for(chind in 1:dim(tone$lfp)[3])
    rerefLFP[,,chind] <- tone$lfp[,,chind] - refAct
  
  tone$lfp <- rerefLFP
  
  if (dim==1)
    tone$lfp <- aperm(tone$lfp, c(3,2,1))
  
  
  gc()
  
  return(tone)
}


filterSignalProcessor <- function(tn, hpHz=NA, lpHz=NA,filterFUN=butter, order=10, 
                                  dataStr='lfp', taxisStr='taxis', Npad=200, Ntaper=150){

  taxis <- tn[[taxisStr]]
  sr <- 1/mean( diff(.001*taxis) ) # sampling rate in Hz
  
  # ===     set up parameters
  nF       <- sr/2
  lpFrc    <- lpHz/nF 
  hpFrc    <- hpHz/nF 
  
  if ( !is.na(hpFrc) &  !is.na(lpFrc) ){
    type <- 'pass'
    cutF <- c(hpFrc,lpFrc)
  }
  
  if ( !is.na(hpFrc) &  is.na(lpFrc) ){
    type <- 'high'
    cutF <- hpFrc
  }
  
  if ( is.na(hpFrc) &  !is.na(lpFrc) ){
    type <- 'low'
    cutF <- lpFrc
  }
  
  if ( is.na(hpFrc) &  is.na(lpFrc) )
    stop('Need to specify at least one filter')
  
  bf       <- filterFUN(order, cutF, type=type)
  lfp      <- tn[[dataStr]]
  Nchan    <- dim(lfp)[3]
  Ntst     <- dim(lfp)[2]
  
  padLFP    <- lfp
  centerTax <- 1:Ntst
  taper     <- rep(1,Ntst)

  if (Npad > 0){# pad the data
    Npad <- min(Npad, round(dim(lfp)[2]/2) )
    padLFP <- zbind( lfp[,seq(Npad,by=-1,to=1),,drop=F], lfp, lfp[,seq(Ntst,by=-1,length=Npad),,drop=F], dimension = 2)
    centerTax <- seq(Npad+1, by=1, length=Ntst)
    
    taper     <- rep(1,dim(padLFP)[2])
    if (Ntaper>0){
      Ntaper   <- min(Ntaper, Npad)
      taper    <- cosine.window( dim(padLFP)[2], Ntaper )
    }
  }
  
  for (cx in 1:Nchan){
    tmp           <- sapply( 1:dim(padLFP)[1], function(nx){ filtfilt(bf, taper*padLFP[nx,,cx]) })
    lfp[,,cx]    <- t(tmp[centerTax,,drop=F])
  }
  #dim(tn$lfp) <- c(dim(tn$lfp),1)
  
  tn[[dataStr]] <- lfp
  return(tn)
}


fftSignalProcessor <- function(tone,trng=NULL,taperms=NULL){
  ##careful pre-specifying taperms, if leading to greater than half the series, expect crash
  if (is.null(trng))
    trng <- range(tone$taxis)
  
  tone$dt       <- mean(diff(tone$taxis)) # temporal resolution in ms
  tone$SR       <- 1000/tone$dt
  
  if(is.null(taperms) & (25/tone$dt < dim(tone$lfp)[2]/4)){
    taperms<-25
  }
  if(is.null(taperms)){
    taperms<-(dim(tone$lfp)[2]*tone$dt)/20
  }
  
  print(paste("Taper length (one side), in ms: ", taperms, sep = ""))
  
  Ntaper <- round(taperms/tone$dt)
  
  tnd <- which(tone$taxis>=trng[1]&tone$taxis<=trng[2])
  if ( (length(tnd)%%2)==1 ){
    warning('use even number of time-points. dropping last time point')
    tnd <- tnd[-length(tnd)]
  }
  getSpecByChan <- function(cx){
    tfft      <- t(tapered_fft( t(tone$lfp[,tnd,cx]), Ntaper=Ntaper))
    #spc      <- 2*Mod( spc )/dim(spc)[2] 
  }
  
  tst       <- lapply( 1:dim(tone$lfp)[3], getSpecByChan)
  tone$tfft <- array( unlist(tst), c(dim(tst[[1]]), length(tst)) ) 
  
  Nefrq             <- 1+dim(tone$tfft)[2]/2
  tone$faxis        <- seq(0, by = 1/(.001*tone$dt*dim(tone$tfft)[2]), length = Nefrq)
  tone$faxislog     <- log(tone$faxis)
  tone$faxislog[1]  <- tone$faxislog[2] - 2*diff(tone$faxislog[c(2,3)]) 
  tone$spc          <- 2*Mod(tone$tfft[,1:Nefrq,,drop=F])#^2#/dim(tone$tfft)[2]
  # divided by number of time-points is already done in tapered_fft
  # ^2 is the wrong thing to do...
  if(length(dim(tone$spc))==2)
    dim(tone$spc)<- c(dim(tone$spc),1)
  
  tone$spc[,1,]     <- 0.5*tone$spc[,1,]
  tone$spc[,Nefrq,] <- 0.5*tone$spc[,Nefrq,]
  tone$logspc       <- log(tone$spc)
  
  return(tone)
}


cwtSignalProcessor <- function(tone,noctave=6,nvoice=3,w0=pi, save=F, saveBaseDir='/Volumes/EEGCyno/EEG/'){
  
  tone$dt       <- mean(diff(tone$taxis)) # temporal resolution in ms
  tone$SR       <- 1000/tone$dt
  nF            <- tone$SR/2  
  
  # old version with inverted frequency axis
  #faxis  <- w0/(2*pi) * nF * 2^(seq( to=0, by=1/nvoice, length=nvoice*noctave ))
  #faxlog <- seq( to=0, by=1/nvoice, length=nvoice*noctave )
  
  faxlog <- seq( from=0, by=-1/nvoice, length=nvoice*noctave )
  faxis  <- w0/(2*pi) * nF * 2^(faxlog)
  yConv  <- function(x){ log(2*pi*x/(w0*nF))/log(2)   }
  
  ## wavelet transform for one trial and one channel
  getCWT2   <- function(n,mn=NA,cx=1){ 
    rtrn <- cwt(tone$lfp[n,,cx]-mn, noctave=noctave,  nvoice=nvoice, w0=w0, twoD=TRUE, plot=F)
    return(rtrn)
  }
  
  mnRsp    <- apply(tone$lfp,c(2),mean)
  tst      <- getCWT2(1,mn=0*mnRsp)
  
  # cwt for a range of trails, fixed channel number
  getCWT     <- function(ind,mn=NA,cx=1){
    tmp      <- sapply(ind, getCWT2, mn=mn,cx=cx )
    dim(tmp) <- c(dim(tst),length(ind))
    return(tmp)
  }
  
  tone$cwt <- array(c(0), c(dim(tone$lfp),dim(tst)[2] ) )
  chind <- 1:dim(tone$lfp)[3]
  for (cx in 1:length(chind) ){
    print( paste(tone$exp, cx, '...', sep=' ') )
    
    #mnRsp             <- apply(tone$lfp[,,cx],c(2),mean)
    singleTrialCWT    <- getCWT( 1:dim(tone$xpp)[1], mn=0*mnRsp, cx=cx )
    tone$cwt[,,cx,] <- aperm(singleTrialCWT, c(3,1,2) )
    rm(singleTrialCWT)
    gc()
  }
  #tone$cwt        <- Mod(tone$cwt)
  
  tone$yConv       <- yConv
  tone$faxiscwt    <- faxis
  tone$faxislogcwt <- faxlog
  
  
  if (save & dim(tone$dfline)[1]==1){
    
    print('saving CWT results to file...')
    saveDir <- paste(saveBaseDir,'/', tone$project, '/', tone$dfline$session, '/', tone$epoch, '_cwt_',noctave,'_',nvoice,'/',sep='')
    if (!file.exists( saveDir ))
      system(paste('mkdir ', saveDir, sep='') )
    
    xpp   <- tone$xpp
    taxis <- tone$taxis
    faxis <- tone$faxiscwt
    
    save(xpp,   file=paste(saveDir, 'xpp.rda',  sep='') )
    save(taxis, file=paste(saveDir, 'taxis.rda',sep='') )
    save(faxis, file=paste(saveDir, 'faxis.rda',sep='') )
    
    saveChannel <- function(cx,tlfp=tone$cwt){
      print( tone$chStr[cx] )
      chdata    <- tlfp[,,cx,, drop=F]
      save(chdata, file=paste(saveDir,tone$chStr[cx],'.rda',sep=''))
    }
    disc <- lapply(1:dim(tone$cwt)[3], saveChannel)  
    
  }
  
  return(tone)
}

# tst <- compareRepositories(checkPath='/Volumes/rawData/EEG/', refPath='/Volumes/rawDataCyno/EEG/',animal='Ben',checkRecordings=T, checkFiles=F)
compareRepositories <- function(checkPath='/Volumes/rawData/EEG/', refPath='/Volumes/rawDataCyno/EEG/',animal='*',checkRecordings=T, checkFiles=F, fileType='rhd'){
  
  #
  # functions checks if checkPath is redundant relative to refPath
  # if all sessions, recordings and data files in checkPath are found in refPath, checkPath is redundant and can be removed
  #
  path1 <- checkPath
  path2 <- refPath
  
  sessionList1 <- system(paste( 'ls ', path1, '/',animal,'/',sep=''), intern=T)
  sessionList2 <- system(paste( 'ls ', path2, '/',animal,'/',sep=''), intern=T)
  
  missingFrom1_ndx <- which(!sessionList2 %in% sessionList1)   
  missingFrom2_ndx <- which(!sessionList1 %in% sessionList2)  
  
  if ( length(missingFrom2_ndx)>0){
    print('checkPath is NOT redundant. Do not remove checkPath before copying the following sessions:')
    print('Sessions from checkPath that are missing from refPath:')
    print(sessionList1[missingFrom2_ndx])
  }
  
  #if ( length(missingFrom1_ndx)>0){
  #  print('Sessions from refPath that are missing Sessions from checkPath')
  #  print(sessionList2[missingFrom1_ndx])
  #}
  
  if ( length(missingFrom2_ndx)==0){
    
    print('All sessions from checkPath are present in refPath. checkPath is redundant at the session level!')
    
    if (checkRecordings){
      
      print('checking redundancy at the recording level:')
      
      recordingList1 <- system(paste( 'ls -1 ', path1, '/',animal,'/*/',sep=''), intern=T)
      recordingList2 <- system(paste( 'ls -1 ', path2, '/',animal,'/*/',sep=''), intern=T)
      
      valSubDir <- function(ix, rL=recordingList1){ !substr(rL[ix],1,1)%in%c('/','') }
      
      recordingList1 <- recordingList1[ sapply(1:length(recordingList1),valSubDir, rL=recordingList1 ) ]
      recordingList2 <- recordingList2[ sapply(1:length(recordingList2),valSubDir, rL=recordingList2 ) ]
      
      missingRecordingFrom2_ndx <- which(!recordingList1 %in% recordingList2)
      missingRecordingFrom1_ndx <- which(!recordingList2 %in% recordingList1)
      
      if (length(missingRecordingFrom2_ndx)>0){
        print('checkPath is NOT redundant. Do not remove checkPath before copying the following directories:')
        print('Recordings from checkPath that are missing in refPath:')
        print(recordingList1[missingRecordingFrom2_ndx])
      }
      
      #if (length(missingRecordingFrom1_ndx)>0){
      #  print('Recordings from path 1 that are missing in path 2:')
      #  print(recordingList2[missingRecordingFrom1_ndx])
      #}
      if ( length(missingRecordingFrom2_ndx)==0){
        print('All recordings from checkPath are present in refPath. checkPath is redundant at the recording level!')
        
        if (checkFiles){
          print('checking Files')
          
          fileListCheck <- rep(F, length(sessionList1) )
          for (sx in 1:length(sessionList1) ){
            fileList1 <- system(paste( 'ls ', path1, '/',animal,'/',sessionList1[sx],'/*/*.',fileType,sep=''), intern=T)
            fileList2 <- system(paste( 'ls ', path2, '/',animal,'/',sessionList1[sx],'/*/*.',fileType,sep=''), intern=T)
            
            extractFileName <- function(pathFile){ tmp=strsplit(pathFile,'/')[[1]];return(tmp[length(tmp)])  }
            fileNameList1 <- sapply(fileList1, extractFileName ) 
            fileNameList2 <- sapply(fileList2, extractFileName )
            names(fileNameList1) <- 1:length(fileNameList1)
            names(fileNameList2) <- 1:length(fileNameList2)
            
            #fileListCheck[sx] <- identical(fileNameList1, fileNameList2)
            fileListCheck[sx] <- length(which(!fileNameList1 %in% fileNameList2))==0
            
            
            if (fileListCheck[sx])
              print(paste(sessionList1[sx],'check') )
            
            if (!fileListCheck[sx])
              print(paste(sessionList1[sx],'fail') )
            
            
          }# check each session (all recordings of a session)
          
          missingFileIdx <- which(!fileListCheck)
          
          if (length(missingFileIdx)==0)
            print('checkPath is redundant at the file level.')
          
          if(length(missingFileIdx)>0){
            print('checkPath is NOT redundant. Do not remove checkPath before checking data files in the following sessions:')
            print('Sessions with data files from checkPath that are missing in refPath:')
            
            print(sessionList1[missingFileIdx])
          }
          
        } # check File
      } # recordings are Fine
      
      
      
    } # check Recordings
  } # Sessions are fine
  
  return(list(missingFrom1_ndx=missingFrom1_ndx, missingFrom2_ndx=missingFrom2_ndx))  
}



showSignalProcessorCWT <- function(tone,cx=1,subs=NULL,zlim=NULL,trng=NULL,absval=T,...){
  
  if (is.null(subs))
    subs <- rep(T, dim(tone$lfp)[1])
  
  if( is.null(tone$xpp$valTrial))
    tone$xpp$valTrial <- T
  
  tndx <- 1:dim(tone$cwt)[2]
  if (!is.null(trng))
    tndx <- tone$taxis>=trng[1] & tone$taxis<=trng[2]
    
  subs <- subs & tone$xpp$valTrial
  
  dimsx <- c(2,3)
  if (length(which(subs))==1)
    dimsx <- c(1,2)
  
  if (absval)
    mnPow    <- apply(Mod(tone$cwt[which(subs),tndx,cx,]),dimsx, mean)
  
  if (!absval)
    mnPow    <- apply(tone$cwt[which(subs),tndx,cx,],dimsx, mean)
  
  mxPow    <- max(c(mnPow))
  xlim     <- range(tone$taxis[tndx])
  ylim     <- range(tone$faxislogcwt)
  nf       <- max(tone$faxiscwt)
  
  flab <- c(1,2,3,4,5,10,20,30,40,50,100,200,300,400,500)
  
  if(is.null(zlim))
    zlim <- c(0,mxPow)
  
  image(tone$taxis[tndx], tone$faxislogcwt[dim(mnPow)[2]:1], mnPow[,dim(mnPow)[2]:1],xlim=xlim,ylim=ylim,zlim=zlim,
        col=potential.colors(256),useRaster=T, yaxt='n',xlab='Time [ms]', ylab='Frequency [Hz]',... )
  #points(tone$taxis,faxlog[5] + mnRsp/20, type='l')
  
  axis(2,at=tone$yConv(flab),  label=flab)
  #abline(h=yConv(c(10,20,30,40,50,60,70)), col='white',lty=2 )
  abline(v=0,col='white')
  
  return(zlim)
}


showSignalProcessorGraph <- function(tone,dataStr='graph', chIndx=NULL){
  
  if (is.null(chIndx))
    chIndx <- 1:dim(tone$lfp)[3]
  
  plot(NA, xlim=range(tone$channelPosition$x[1:length(chIndx)])+c(-.1,1), ylim=range(tone$channelPosition$y[1:length(chIndx)])+c(-1,1),
       xlab='',ylab='',xaxt='n',yaxt='n',bty='n')#, xlab='time', ylab='Channel Nr',bty='n',yaxt='n',xaxt='n')
  points( tone$channelPosition$x[chIndx], tone$channelPosition$y[chIndx], pch=20)
  
  
  graph <- tone[[dataStr]]
  print(dim(graph))
  
  
  addLine <- function(cx, edgedf=edgeIndx,...){
    cx1 <- edgedf$Px[cx]
    cx2 <- edgedf$Py[cx]
    
    points( tone$channelPosition$x[c(cx1,cx2)],tone$channelPosition$y[c(cx1,cx2)], type='l',... )
  }
  
  edgeIndx <- data.frame(Px=1 + (which(graph==-1)-1) %% dim(graph)[1],Py=1 + (which(graph==-1)-1) %/% dim(graph)[1])
  disc <- lapply(1:dim(edgeIndx)[1], addLine, edgedf=edgeIndx,col='blue' )
  
  
  edgeIndx <- data.frame(Px=1 + (which(graph==+1)-1) %% dim(graph)[1],Py=1 + (which(graph==+1)-1) %/% dim(graph)[1])
  disc <- lapply(1:dim(edgeIndx)[1], addLine, edgedf=edgeIndx,col='red' )
  
  disc <- showBrainOutline(scaleX=.14,scaleY=.15, lwd=2)
  
}


frml2str <- function(frml){
  Nfrml       <- length(frml)
  tmp         <- as.character(frml)[[Nfrml]]
  tmpLst      <- strsplit(tmp,' ')
  IVlist      <- tmpLst[[1]][ !tmpLst[[1]]%in%c('+','*',':','/') ]
  frmlStr <- paste(IVlist,collapse='_')
  return(frmlStr)
}

three2twoDimSignalProcessor <- function(lfp){
  tmp1      <- aperm(lfp, c(3,2,1))
  dim(tmp1) <- c(dim(tmp1)[1], dim(tmp1)[2]*dim(tmp1)[3] )
  gc()
  return( t(tmp1) )
}
two2threeDimSignalProcessor <- function(lfp, dm=NULL, Ntrial=NULL, Ntst=NULL ){
  
  if (is.null(dm) & is.null(Ntrial) & is.null(Ntst) )
    stop('specify at least one of the following: dm, Ntrial or Ntst')

  if (is.null(dm)){
    dm <- c(NA, NA, dim(lfp)[2])
    if ( !is.null(Ntst) )
      Ntrial <- dim(lfp)[1]/Ntst
    
    if ( !is.null(Ntrial) )
      Ntst <- dim(lfp)[1]/Ntrial
    
    dm[1] <- Ntrial
    dm[2] <- Ntst
  }
  
  tmp1      <- aperm(lfp, c(2,1))
  dim(tmp1) <- c(dm[3], dm[2], dm[1] )
  
  rm(lfp)
  gc()
  
  return( aperm(tmp1,c(3,2,1)) )
}


trial2contSignalProcessor <- function(tone){
  
  
  copyStr <- c('chStr','GMC','channelPosition','xpp','dfline','lfp','taxis','filter')
  copyNdx <- which(copyStr %in% names(tone))
  

  
  cnt          <- tone[copyStr[copyNdx]] # cnt will contain a single trial with all trials concatenated
  cnt$xpp      <- tone$xpp[1,]
  cnt$lfp      <- three2twoDimSignalProcessor(tone$lfp)
  dim(cnt$lfp) <- c(1,dim(tone$lfp)[1]*dim(tone$lfp)[2],dim(tone$lfp)[3] )
  
  
  sr           <- mean(diff(tone$taxis)) 
  cnt$taxis    <- sr *(1:dim(cnt$lfp)[2]-1)
  
  return(cnt)  
}

avgSignalProcessor <- function(tone, frml=~dummy, IVndx=NULL, trng=NULL,returnDF=T, FUN=mean,
                               saveData=F,saveBaseDir='/Volumes/EEGCyno/AVG/',analysisDir=NULL,
                               chIndx=NULL,subsExpr=NULL,dataStr='lfp',taxisStr='taxis',
                               Nrepmax=NA,...){
  
  # Nrepmax is the maximal number of trials to include in the averaging process
  # this value can be helpful when comparing signal-to-noise of the averaged data 
  # between diferent days.
  # If Nrepmax is set to 50, the program will find the first 50 instances for each of 
  # the conditions and only use these for averaging

  
  tone$xpp$dummy <- c(1)
  if ('rejected'%in%names(tone))
    tone$rejected <- c(1)
  
  if (!is.null(subsExpr))
    subs <- with(tone$xpp, eval(subsExpr) )
  
  if (is.null(subsExpr))
    subs <- rep(T, dim(tone[[dataStr]])[1] )
  
  subs[is.na(subs)] <- FALSE

  # subset the data
  ttn   <- subs.data(tone,subs)
  
  # handle Nrepmax if not set to NA
  if(!is.na(Nrepmax) ){
    
    xpp   <- ttn$xpp
    IV    <- model.frame(frml, data = xpp )
    
    mxIVtable <- max(table(IV))
    
    # this breaks it down in case of more than two IVs
    cnd    <- tapply(IV, IV)
    cndlst <- unique(sort(cnd))
    
    # initiate the condition counter as 0
    ttn$xpp$NrepCnd <- 0
    for (cix in cndlst){
      tcix <- which(cnd == cix)
      ttn$xpp$NrepCnd[tcix] <- 1:length(tcix)
    }
    # subset to the first Nrepmax trials
    ttn   <- subs.data(ttn, ttn$xpp$NrepCnd <= Nrepmax)
  }
  
  # at this point $lfp is not usefull anymore if dataStr is unqual to 'lfp'
  # hence replace $lfp with [[dataStr]]
  if (!identical(dataStr,'lfp'))
    ttn$lfp <- ttn[[dataStr]]
  
  # extract taxis and lfp based on taxisStr and dataStr
  taxis <- ttn[[taxisStr]]
  lfp   <- ttn[[dataStr]]
  if (identical(dataStr,'cwt'))
    lfp <- Mod(lfp)
  
  xpp   <- ttn$xpp
  IV    <- model.frame(frml, data = xpp )
  
  mxIVtable <- max(table(IV))
  doAVG     <- mxIVtable > 1
  
  Ntrial <- dim(lfp)[1]
  Ntst   <- dim(lfp)[2]
  Nchan  <- dim(lfp)[3]
  
  ## ======================== define helper function
  tmeanPSTH <- function(psth, FUN=mean){
    #print(dim(psth))
    ## little helper function that returns valid mean psths even
    ## if the number of repetitions is one or zero
    
    Ntrial <- dim(lfp)[1]
    Ntst   <- dim(lfp)[2]
    Nchan  <- dim(lfp)[3]
    
    res <- array(NA,c(Ntst,Nchan))  
    
    if(is.array(psth)&length(dim(psth))==3 )
      res <- apply(psth,c(2,3),FUN)
    
    if( is.array(psth) & length(dim(psth))==2 & identical(FUN,mean))
      res <- psth
    
    return(res)
  }
  
  tmeanPSTH_4D <- function(psth, FUN=mean){
    #print(dim(psth))
    ## little helper function that returns valid mean psths even
    ## if the number of repetitions is one or zero
    
    #Ntrial <- dim(lfp)[1]
    Ntst   <- dim(psth)[2]
    Nchan  <- dim(psth)[3]
    Nfreq  <- dim(psth)[4]
    
    
    res <- array(NA,c(Ntst,Nchan,Nfreq))  
    
    if(is.array(psth)&length(dim(psth))==4 )
      res <- apply(psth,c(2,3,4),FUN)
    
    if( is.array(psth) & length(dim(psth))==3 & identical(FUN,mean))
      res <- psth
    
    return(res)
  }
  # define the sums of square function
  SSQFUN  <- function(x,n=NULL){  return(sum( (x-mean(x))^2 )) }
  
  if (!doAVG){# if not averaging needs to be done, just return the subsetted input data
    rxpp <- xpp
    rlfp <- lfp
    rtn  <- ttn
  }
  
  if (doAVG){
    tapplyList <- lapply(1:dim(IV)[2], function(ix)IV[,ix]) # drops levels with zero entries
    
    if (length(dim(lfp))==3)
      mn         <- tapply(1:dim(lfp)[1], tapplyList, function(ind){tmeanPSTH(lfp[ind,,,drop=F],FUN=FUN) })
    
    if (length(dim(lfp))==4)
      mn         <- tapply(1:dim(lfp)[1], tapplyList, function(ind){tmeanPSTH_4D(lfp[ind,,,,drop=F],FUN=FUN) })
    
    Ntrial     <- tapply(1:dim(lfp)[1], tapplyList, length)
    mnF        <- mn[which(!is.na(Ntrial))]
    
    print(Ntrial)
    
    #rlfp <- array(unlist(mn), c(length(mn),Ntst, Nchan) )
    if (length(dim(lfp))==3){
      tmplfp <- array(unlist(mn), c(Ntst, Nchan,length(mnF)) )
      rlfp   <- aperm(tmplfp, c(3,1,2) )
      
      dimnames(rlfp)[1] <- list(1:dim(rlfp)[1])
      dimnames(rlfp)[2] <- dimnames(tone[[dataStr]])[2]
      dimnames(rlfp)[3] <- dimnames(tone[[dataStr]])[3]
      names(dimnames(rlfp))    <- names(dimnames(tone[[dataStr]]))
      names(dimnames(rlfp))[1] <- 'condition' 
    }
    
    if (length(dim(lfp))==4){
      Nfreq <- dim(lfp)[4]
      tmplfp <- array(unlist(mn), c(Ntst, Nchan, Nfreq, length(mnF)) )
      rlfp   <- aperm(tmplfp, c(4,1,2,3) )
      
      dimnames(rlfp)[1] <- list(1:dim(rlfp)[1])
      dimnames(rlfp)[2] <- dimnames(tone[[dataStr]])[2]
      dimnames(rlfp)[3] <- dimnames(tone[[dataStr]])[3]
      dimnames(rlfp)[4] <- dimnames(tone[[dataStr]])[4]
      
      names(dimnames(rlfp))    <- names(dimnames(tone[[dataStr]]))
      names(dimnames(rlfp))[1] <- 'condition' 
    }
    
    gc()
    ## =============================== get the new xpp 
    useBin       <- sapply( 1:dim(xpp)[2], function(cx){is.numeric(xpp[,cx]) } )
    useIndx      <- which(useBin)
    rxpp         <- aggregate( xpp[,useIndx], by=tapplyList, mean)  
    rxpp$Ntrials <- aggregate( xpp[,useIndx[1]], by=tapplyList, length)$x
    rxpp$valTone <- c(T)
    
    rtn             <- ttn
    rtn$xpp         <- rxpp
    rtn$xpp$valTone <- T
    rtn[[dataStr]]  <- rlfp
    # return averaged results also in $lfp, just in case...
    if (!identical(dataStr,'lfp'))
      rtn$lfp <- rlfp
    
  }
  
  if(frml==~seqNr){
    #browser()
  }
  
  if (saveData){
    print('saving...')
    
    if (is.null(tone$project))
      tone$project <- tone$dfline$project[1]
    
    if (is.null(analysisDir))
      analysisDir <- frml2str(frml)
    
    epochExtn <- ''
    if (!dataStr=='lfp')
      epochExtn <- paste('_',dataStr, sep='')
    
    saveDir <- paste(saveBaseDir,tone$project,'/',as.character(tone$exp),'/',tone$epoch,epochExtn,'/',analysisDir,'/',sep='')
    
    if (!is.na(Nrepmax))    
      saveDir <- paste(saveBaseDir,tone$project,'/',as.character(tone$exp),'/',tone$epoch,epochExtn,'/',analysisDir,'_Nrepmax',Nrepmax,'/',sep='')
    
    ##this was #exists() and getting remade every time
    if (!file.exists(saveDir))
      system( paste('mkdir -p ', saveDir,sep=''), intern=T )
    
    print(c(dim(xpp)[1], dim(rxpp)[1]))
    
    xpp <- rxpp
    save(xpp, file=paste(saveDir,'xpp.rda',sep='') )
    save(taxis, file=paste(saveDir,'taxis.rda',sep='') )
    saveChannel <- function(cx,tlfp=rlfp){
      chdata <- tlfp[,,cx]
      dim(chdata) <- c( dim(tlfp)[1:2],1)
      save(chdata, file=paste(saveDir,tone$chStr[cx],'.rda',sep=''))
    }
    disc <- lapply(1:dim(lfp)[3], saveChannel)  
  }
  
  
  res <- rlfp
  if (returnDF)
    res <- rtn
  
  gc()
  return( res )
  
}


extractMaxSignalProcessor <- function(tone,trng=NULL,FUN=max,dataStr='lfp',taxisStr='taxis',returnTime=T){
  
  lfp <- tone[[dataStr]]
  taxis <- tone[[taxisStr]]
  
  if (is.null(trng))
    trng <- range(taxis)
  
  tndx    <- which(taxis>=trng[1] & taxis<=trng[2])
  mxIndx  <- tndx[1] + apply( lfp[,tndx,,drop=F], c(1,3), function(x){which(x==FUN(x))} )
  timeMat <- array( taxis[mxIndx], dim=dim(mxIndx) )
  maxMat  <- apply( lfp[,tndx,,drop=F], c(1,3), FUN)
  
  res <- timeMat
  if (!returnTime)
    res <- maxMat
  
  return(timeMat)
}


csdSignalProcessor <- function(lfp, spacing=100,lag=1,smooth=F,sd=150){
  
  if (exists('csdUseInd',where=lfp)){
    lfp$lfp        <- lfp$lfp[,,lfp$csdUseInd,drop=F]
    lfp$chDepthInd <- lfp$chDepthInd[lfp$csdUseInd,drop=F]
  }
  
  Nchan         <- dim(lfp$lfp)[3]
  chanInd       <- lfp$chDepthInd;
  
  #csd            <- list()
  #csd$daxis      <- -1*spacing * (lfp$chDepthInd-1) # seq(0,by=-1*spacing,length=length(chanInd)) #electrode spacing in microns
  #mnlfp     <- lfp$lfp
  
  #if(exists('logtaxis',where=lfp))
  #  csd$logtaxis <- lfp$logtaxis
  
  if (smooth){
    
    sker <- dnorm(spacing*seq(-3,by=1,to=3),m=0,sd=sd)
    sker <- sker/sum(sker)
    print(sker)
    
    tmp       <- array(0, dim(lfp$lfp)+c(0,0,6) )
    smLFP     <- array(0, dim(lfp$lfp)+c(0,0,6) ) # set up smoothed LFP
    
    # extrapolate channels for smoothing 
    tmp[,,3+(1:Nchan)] <- lfp$lfp
    
    for (k in c(3,2,1) )
      tmp[,,k] <- tmp[,,k+1] - (tmp[,,k+2]-tmp[,,k+1])
    
    for (k in Nchan+3+(1:3) )
      tmp[,,k] <- tmp[,,k-1] - (tmp[,,k-2]-tmp[,,k-1])
    
    
    #gpuKer <- gpuMatrix( array(sker, c(length(sker),1))) ;
    for (k in (1:Nchan)+3){
      for (tx in 1:dim(tmp)[2] ){
        smLFP[,tx,k] = tmp[,tx,seq(k-3,by=1,to=k+3)]%*%sker
        #gpuA         <- gpuMatrix( tmp[,tx,seq(k-3,by=1,to=k+3)], dim(tmp[,tx,seq(k-3,by=1,to=k+3)]) )
        #smLFP[,tx,k] <- gpuA%*%gpuKer
      }
    }
    lfp$lfp <- smLFP[,,(1:Nchan)+3,drop=F] # replace standard lfp with spatially smoothed version for consistency
  }
  
  lfp$csd <- array(0, dim(lfp$lfp))
  tmp     <- lfp$lfp
  for (i in (1+lag):(Nchan-lag))
    lfp$csd[,,i]  =  .001 * (2*tmp[,,i] - tmp[,,i-lag] - tmp[,,i+lag]) / ((lag*spacing/1000)^2);  
  
  lfp$avrec <- apply(abs(lfp$csd),1,mean)
  
  return(lfp)
}


showSignalProcessor <- function(tone, frml=~dummy, colFUN=bluered.colors, colLst=NULL, IVndx=NULL, mnFUN=mean,
                                scl='global',sclfctr=1.5, chIndx=NULL, trng=NULL,bty='n',tcex=1,butterfly=F,muaplot=F,muaColLst=NULL,
                                add=F,addyaxis=T,logT=F,vlineat=NULL,mxNIV=23,invertIV=F,yshft=0,showSE=F,axlab=T,atms=NULL,
                                subsExpr=NULL,taxisStr='taxis',dataStr='lfp',posStr='channelPosition',chStrLabel='Intan',...){
  
  if (! posStr %in% names(tone) ){
    Nchan <- dim(tone$lfp)[3]
    tone$channelPosition <- data.frame(xpos=1+(0:(Nchan-1))%/%8, ypos=1+(0:(Nchan-1))%%8,  zpos=rep(0,Nchan), chStr=paste('ch',1:Nchan,sep=''))
  }
  
  if (!'chStr'%in% names(tone[[posStr]] ))
    tone[[posStr]]$chStr <- tone$chStr
  
  tone$xpp$dummy <- c(1)
  par(mfcol=c(1,1),mar=c(.1,.1,.1,.1))
  
  if (!is.null(frml)){
    mn             <- avgSignalProcessor(tone, frml=frml, IVndx=IVndx, trng=trng,returnDF=T, chIndx=chIndx, subsExpr=subsExpr, dataStr=dataStr, taxisStr=taxisStr, FUN=mnFUN,...)
    # if avg does actual averaging, it returns the results as lfp and taxis, regardless of the dataStr and taxisStr values.
    # if no averaging is necessary, it returns the original data frame...
    # this needs to be fixed in avg, otherwise the results differ 
    # 
    mn[[dataStr]]  <- mn$lfp
    mn[[taxisStr]] <- mn$taxis
    }
  else{
    mn <- tone
    mn$xpp$Ntrials==1
  }
  gc()
  
  se <- NULL
  if (showSE & min(mn$xpp$Ntrials)>1 & !is.null(frml)){# only try to calculate se if at least two trials in each condition.
    se <- avgSignalProcessor(tone, frml=frml, IVndx=IVndx, trng=trng,returnDF=T, chIndx=chIndx, subsExpr=subsExpr, dataStr=dataStr, taxisStr=taxisStr, FUN=na.se,...)
  }
  
  Ncnd  <- dim(mn[[dataStr]])[1]
  Ntnt  <- dim(mn[[dataStr]])[2]
  
  if (butterfly & Ncnd>1){
    ga     <- avgSignalProcessor(tone, frml=~dummy, IVndx=IVndx, trng=trng,returnDF=T, chIndx=chIndx, subsExpr=subsExpr, dataStr=dataStr, taxisStr=taxisStr,...)
    ga$lfp <- array( rep(ga$lfp, each=dim(mn$lfp)[1] ), dim(mn$lfp) )
    mn$lfp <- mn$lfp - ga$lfp
  }
  
  
  if (is.null(chIndx))
    chIndx          <- 1:dim(mn[[dataStr]])[3]
  
  taxis <- 1:Ntnt
  if (!is.null(taxisStr))
    taxis <- tone[[taxisStr]]
  
  if (logT & taxisStr=='taxis'& exists('logtaxis', tone))
    taxis <- tone$logtaxis
  
  xlim <- range(taxis)
  if (!is.null(trng))
    xlim <- c(max(trng[1],min(taxis)),min(trng[2],max(taxis)))
  
  channelPosition <- tone[[posStr]][chIndx,]
  
  tnd <- which(taxis >= xlim[1] & taxis <= xlim[2])
  if (length(tnd)==0){
    warning('showSignalProcessor: no valid times found; using the entire time axis')
    tnd <- 1:length(taxis)
    xlim <- range(taxis)
  }
  tsqueezrng <- c(0,0.85)
  taxis <- squeeze(taxis[tnd],from.range=xlim,to.range=tsqueezrng)
  
  lfp   <- mn[[dataStr]][,tnd,chIndx, drop=F]
  Ncnd  <- dim(mn[[dataStr]])[1]
  Ntnt  <- length(tnd)
  Nchan <- length(chIndx)
  
  #lfp   <- lfp[,tnd,]
  dim(lfp) <- c(Ncnd,length(tnd),Nchan)
  
  NIV <- min( c( dim(lfp)[1], mxNIV) )
  if (is.null(colLst)){
    colLst <- colFUN(NIV)
    if(NIV==1)
      colLst <- c('black')
  }
  
  if (is.null(muaColLst))
    muaColLst <- colLst
  
  if (!add){
    plot(NA, xlim=range(channelPosition$xpos)+c(-.1,1), ylim=range(channelPosition$ypos)+c(-.75,.75), xlab='', ylab='',bty='n',yaxt='n',xaxt='n',...)
    #axis(1, at=1+squeeze(atms, from.range=xlim,to.range=tsqueezrng), labels=atms)
  }
  #atms <- pretty(tone$taxis[tnd])
  
  if (!is.null(atms))
    atms <- squeeze( atms, from.range=xlim,to.range=tsqueezrng)
  
  if (is.null(atms)){
    atms <- pretty(taxis)
    atms <- atms[-c(1,length(atms))]
  }
  #browser()
  channelPosition$ty <- channelPosition$y - channelPosition$x/(max(channelPosition$x)+1)
  cxbl               <- which( channelPosition$ty == min(channelPosition$ty) )
  
  
  for (cx in 1:dim(lfp)[3]){
    
    sclx <- switch(scl,
                   'channel' =  sclfctr*na.max(abs(c(lfp[,,cx]))),
                   'numeric' =  sclfctr                          ,
                   'global'  =  sclfctr*na.max(abs(c(lfp[,,]  ))) 
    )
    
    atmv <- pretty( (c(-.5,.5)-yshft) *sclx,n=5)
    
    
    if ( (bty!='n' | cx==cxbl) & !add ){
      
      # left box
      if (bty %in% c('L','l','o') | cx==cxbl)
        points( rep(min(taxis),2) + channelPosition$x[cx], c(-0.45,0.45)+ channelPosition$y[cx], type='l',lwd=.5 )
      
      # right box
      if (bty %in% c('o'))
        points( rep(max(taxis),2) + channelPosition$x[cx], c(-0.45,0.45)+ channelPosition$y[cx], type='l',lwd=.5 )
      
      # top box
      if (bty %in% c('o'))
        points( taxis[c(1,length(taxis))] + c(.05,-.05) + channelPosition$x[cx],  rep(channelPosition$y[cx]+.5,2), col='black', type='l',lwd=.5 )
      
      # bottom box
      if (bty %in% c('o','l','L') )#| cx==cxbl)
        points( taxis[c(1,length(taxis))] + c(.05,-.05) + channelPosition$x[cx],  rep(channelPosition$y[cx]-.5,2), col='black', type='l',lwd=.5 )
      
      addXTick <- function(atx){
        points(rep( squeeze(atx, from.range=xlim,to.range=tsqueezrng),2)+ channelPosition$x[cx], c(-0.5,-0.53)+ channelPosition$y[cx], type='l',lwd=.5  )
        if (cx==cxbl & axlab)
          text(       squeeze(atx, from.range=xlim,to.range=tsqueezrng)   + channelPosition$x[cx],   -0.58      + channelPosition$y[cx], atx, cex=tcex )
      }
      if (bty %in% c('o','l','L') | cx==cxbl)
        disc <- lapply( atms[-c(1,length(atms))], addXTick )
      
      addYTick <- function(aty){
        points( c(-.03,0)+ channelPosition$x[cx], rep(aty/sclx,2)+ channelPosition$y[cx]+ yshft, type='l',lwd=.5  )
        if (cx == cxbl & axlab)
          text(    -.05+ channelPosition$x[cx], aty/sclx + channelPosition$y[cx]+ yshft, aty, adj=1, cex=tcex )
      }
      if (bty %in% c('o','l','L') | cx==cxbl)
        disc <- lapply( atmv[ -c(1,length(atmv)) ], addYTick )
      
      #if(cx==1)
      #not sure what this is doing; does not seem to be correct
      #        text(squeeze(c(0),from.range=xlim,to.range=tsqueezrng) + channelPosition$x[cx], aty[2]/sclx + channelPosition$y[cx], aty[2],adj=1.2, cex=.75 )
    }
    
    # add x-axis ticks for all panels
    addXTick <- function(atx){
      points(rep(atx,2)+ channelPosition$x[cx], c(-0,-0.03)+ channelPosition$y[cx], type='l',lwd=.5  )
    }
    disc <- lapply( atms, addXTick )
    
    NIVorder <- 1:NIV
    if (invertIV)
      NIVorder <- seq(NIV, by=-1,to=1)
    
    if (showSE & !is.null(se)){
      for (ix in NIVorder ){
        
        tmpcol <- col2rgb( colLst[ix])
        
        cll=rgb(red=tmpcol[1]/255, green=tmpcol[2]/255, blue=tmpcol[3]/255, alpha=0.2)
        errorband( taxis+channelPosition$x[cx], lfp[ix,,cx]/sclx + channelPosition$y[cx] + yshft, 2*se[[dataStr]][ix,tnd,cx]/sclx, col=cll, border=F )
      }
    }# error bands
    
    for (ix in NIVorder ){
      ##added an all to avoid warnings (do we ever get NA values for only some channels?)
      if (!all(!is.na(lfp[ix,1,]))){
        warning('Some or all channels are generating NA values when conditioned on the IV of interest.')
      }
      
      if (!muaplot)
        points(      taxis+channelPosition$x[cx], lfp[ix,,cx]/sclx + channelPosition$y[cx] + yshft, col=colLst[ix], type='l',... )
      
      if (muaplot)
        csdplot(     taxis+channelPosition$x[cx], lfp[ix,,cx]/sclx + channelPosition$y[cx]+ yshft, col=c(colLst[ix], muaColLst[ix]), thr=channelPosition$y[cx], type='l',... )
      
    }# averages
    
    
    aty <- pretty(c(-.5,.5)*sclx, n=5,min.n=5 )
    stp <- diff(aty[c(1,2)])
    aty <- c(-1,1)*stp
    points( taxis[c(1,length(taxis))]+ channelPosition$x[cx],  rep(channelPosition$y[cx],2)+ yshft, col='black', type='l',lwd=.3 )
    if (addyaxis)
      points( squeeze(c(0,0),from.range=xlim,to.range=tsqueezrng) + channelPosition$x[cx], ( c(-.25,.25) + yshft )%>=%-.5%<=%.5 + channelPosition$y[cx] + yshft, type='l',lwd=.3 )
    
    addVLine <- function(tm){points( squeeze(c(tm,tm),from.range=xlim,to.range=tsqueezrng) + channelPosition$x[cx], c(-.5,.5) + channelPosition$y[cx] + yshft, type='l',lwd=.3 )}
    if (!is.null(vlineat)){
      disc <- lapply(vlineat, addVLine)
    }# add vertical marker lines 
    
  }
  if (!is.null(chStrLabel)){
    if (chStrLabel%in%names(channelPosition))
      text(channelPosition$x-.05, channelPosition$y+.2, channelPosition[[chStrLabel]], cex = 0.6)
  }
}


showSignalProcessorWrapper <- function(tn,plotFUN=NULL,to.rangeX=c(0,0.8),to.rangeY=c(0,0.8),annotate=T,zlim=NULL,...){
  
  par(mfcol=c(1,1))
  Ncnd  <- dim(tn$lfp)[1]
  Ntst  <- dim(tn$lfp)[2]
  Nchan <- dim(tn$lfp)[3]
  
  par(mfcol=c(1,1),mar=c(.1,.1,.1,.1))
  
  channelPosition <- tn$channelPosition
  plot(NA, xlim=range(channelPosition$x)+c(-.1,1), ylim=range(channelPosition$y)+c(-.1,.5), xlab='', ylab='',bty='n',yaxt='n',xaxt='n')#,...)
  
  
  channelFunction <- function(cx){plotFUN(cx, tn, add=T, to.rangeX=to.rangeX, to.rangeY=to.rangeY,zlim=zlim,...)}
  res <- lapply(1:Nchan, channelFunction )
  
  if (annotate){
    #text( channelPosition$x-0, channelPosition$y+0.6, tn$chStr, cex=.5)
    text( channelPosition$x-0, channelPosition$y-0, channelPosition$chStr, cex=.5)
  }
  return(res)
}


showSignalProcessorXYC <- function(tone, xaxis='xpos',yaxis='ypos', caxis='dummy', colFUN=bluered.colors, colLst=NULL, IVndx=NULL, 
                                   scl='global',sclfctr=1.5, chIndx=NULL, trng=NULL,  tsqueezrng=c(0,0.80),
                                   add=F,addyaxis=T,logT=F,vlineat=NULL,mxNIV=10,muaplot=F,
                                   subsExpr=NULL,taxisStr='taxis',dataStr='lfp',dontSetmfcol=F,...){
  
  tone$xpp$dummy <- c(1)
  if (!dontSetmfcol)
    par(mfcol=c(1,1))
  
  if (!is.null(subsExpr))
    subs <- with(tone$xpp, eval(subsExpr) )
  
  if (is.null(subsExpr))
    subs <- rep(T, dim(tone$lfp)[1] )
  
  subs[is.na(subs)] <- FALSE
  
  ttn   <- subs.data(tone,subs)
  
  mn <- ttn
  
  lfp   <- mn[[dataStr]]
  Ncnd  <- dim(lfp)[1]
  Ntnt  <- dim(lfp)[2]
  Nchan <- dim(lfp)[3]
  
  taxis <- 1:Ntnt
  if (!is.null(taxisStr))
    taxis <- mn[[taxisStr]]
  
  if (logT & taxisStr=='taxis'& exists('logtaxis', tone))
    taxis <- tone$logtaxis
  
  xlim <- range(taxis)
  if (!is.null(trng))
    xlim <- c(max(trng[1],min(taxis)),min(trng[2],max(taxis)))
  
  
  # X Y and C can vary either as a function of a variable in XPP or a position in channelPosition
  # here we set up a lookup table that specifies X Y and C for each condition and each channe
  
  getLayout <- function(axstr){
    chk <- 0
    XLayout <- array(1, c(Ncnd, Nchan))
    
    if (axstr %in% names(mn$channelPosition)){
      chk     <- 1
      Xval    <- mn$channelPosition[[axstr]]
      XLayout <- array( rep(Xval, each=Ncnd),  c(Ncnd, Nchan) )
    }
    
    if (axstr %in% names(mn$xpp)){
      chk     <- 1
      Xval    <- as.numeric(as.factor(mn$xpp[[axstr]]))
      XLayout <- array( rep(Xval, Nchan),  c(Ncnd, Nchan) )
    }
    
    if (chk==0)
      warning(paste(axstr, ' not found. Using dummy variable.'))
    
    return(XLayout)
  }
  if (is.null(chIndx))
    chIndx          <- 1:Nchan
  
  XLayout <- getLayout(xaxis)
  YLayout <- getLayout(yaxis)
  Ctmp    <- as.factor(getLayout(caxis))
  Ctmp <- as.integer(as.character.factor(Ctmp))
  CLayout <- array(Ctmp, dim(XLayout))
  CLayout <- CLayout - min(CLayout[,chIndx]) + 1
  
  tnd <- which(taxis >= xlim[1] & taxis <= xlim[2])
  if (length(tnd)==0){
    warning('showSignalProcessor: no valid times found; using the entire time axis')
    tnd <- 1:length(taxis)
    xlim <- range(taxis)
  }
  
  if (is.null(colLst)){
    colLst <- colFUN( max(CLayout[,chIndx]) )
    if(max(CLayout)==1)
      colLst <- c('black')
  }
  
  taxis <- squeeze(taxis[tnd],from.range=xlim,to.range=tsqueezrng)
  lfp   <- lfp[,tnd,]
  dim(lfp) <- c(Ncnd,length(tnd),Nchan)
  
  if (!add){
    plot(NA, xlim=range(c(XLayout[,chIndx]))+c(-.1,1), ylim=range(YLayout[,chIndx])+c(-1,1), xlab='', ylab='',bty='n',yaxt='n',xaxt='n',...)
    atms <- pretty(taxis[tnd])
    #axis(1, at=1+squeeze(atms, from.range=xlim,to.range=tsqueezrng), labels=atms)
  }
  
  for (cx in chIndx){
    print(paste( 'channel', cx, sep='') )
    for (ix in 1:Ncnd ){
      
      ##added an all to avoid warnings (do we ever get NA values for only some channels?)
      if (!all(!is.na(lfp[ix,1,]))){
        warning('Some or all channels are generating NA values when conditioned on the IV of interest.')
      }else{
        sclx <- switch(scl,
                       'channel'   =  sclfctr*na.max(abs(c(lfp[,,cx]))),
                       'numeric'   =  sclfctr                          ,
                       'condition' =  sclfctr*na.max(abs(c(lfp[ix,,]))),
                       'global'    =  sclfctr*na.max(abs(c(lfp[,,]  )))
        )
        
        
        if (!muaplot){
          points( taxis+ XLayout[ix,cx], lfp[ix,,cx]/sclx + YLayout[ix,cx], col=colLst[CLayout[ix,cx]], type='l',... )
          points( taxis+ XLayout[ix,cx], lfp[ix,,cx]*0 + YLayout[ix,cx], col='gray', type='l', lwd=.5 )
        }
        if (muaplot)
          csdplot(taxis+XLayout[ix,cx], lfp[ix,,cx]/sclx + YLayout[ix,cx], col=c(colLst[CLayout[ix,cx]], 'gray'), thr=YLayout[ix,cx], type='l',... )
        
        
        
        
        if (!is.null(vlineat)){
          addvline <- function(vl){
            points( rep(squeeze(vl,from.range=xlim,to.range=tsqueezrng) ,2)+XLayout[ix,cx], c(0,.8) + YLayout[ix,cx], col='gray', type='l',lwd=.5)
          }
          disc <- sapply(vlineat, addvline)
        }
      }
      
      if (addyaxis){
        aty <- pretty(c(-.5,.5)*sclx, n=5,min.n=5 )
        stp <- diff(aty[c(1,2)])
        aty <- c(-1,1)*stp
        
        points( squeeze(c(0,0),from.range=xlim,to.range=tsqueezrng) + XLayout[ix,cx], aty[c(2,length(aty)-1)]/sclx+ YLayout[ix,cx], type='l',lwd=.3 )
        points( taxis[c(1,length(taxis))]+ XLayout[ix,cx],  rep(YLayout[ix,cx],2), col=CLayout[ix,cx][ix], type='l',lwd=.3 )
        
        #if(cx==1)
        #not sure what this is doing; does not seem to be correct
        #        text(squeeze(c(0),from.range=xlim,to.range=tsqueezrng) + channelPosition$x[cx], aty[2]/sclx + channelPosition$y[cx], aty[2],adj=1.2, cex=.75 )
        
      }
    }
  }
  #text(channelPosition$x-.05, channelPosition$y+.2, chStr, cex = 0.6)
}



showSignalProcessorTrial <- function(tone,trialNr,...){
  
  trial <- subs.data( tone, 1:dim(tone$xpp)[1]==trialNr)
  showSignalProcessor( trial,... )
  
}

old_showSignalProcessorTrial <- function(tone, trNr, scl='global', chIndx=1,colFUN=bluered.colors, trng=NULL, logT=F, spread=T,...){
  
  taxis <- tone$taxis
  if (logT & exists('logtaxis', tone))
    taxis <- tone$logtaxis
  
  
  xlim <- range(taxis)
  if (!is.null(trng))
    xlim <- trng
  
  if (logT & !is.null(trng))
    xlim <- range( taxis[ tone$taxis>=trng[1] & tone$taxis<=trng[2]] )
  
  tnd <- which(taxis >= xlim[1] & taxis <= xlim[2])
  
  colList <- colFUN(length(trNr))
  
  sprd <- 1
  if (!spread)
    sprd <- 0
  
  
  if (spread)
    plot(NA, xlim=xlim, ylim=c(-length(trNr)-1,.5), xlab='time', ylab='IV',...)
  
  if (!spread)
    plot(NA, xlim=xlim, ylim=c(-.9,.9), xlab='time', ylab='IV',...)
  
  for (trx in trNr ){
    sclx <- switch(scl,
                   'channel' =  max(abs(tone$lfp[trNr[trx],,chIndx])),
                   'global'  =  max(abs(tone$lfp[trNr,,chIndx]))
    )
    points( taxis, tone$lfp[trNr[trx],,chIndx]/sclx - sprd*trx, col=colList[trx], type='l' )
  }
  abline(v=0)
  
}


showSignalProcessorMap <- function(tone,eval.at=NULL,mirror=F,flip=F,zlim=NULL,earloberef=F,flat=F,chIndx=NULL,trialIndx=1,showTraces=F,trng=NULL,colmap=jet.colors,colList=bluered.colors,colorPallete=NULL,addContour=F,
                                   fileName='scalpMap.eps',eps=F,pdf=T,mfrow=NULL,setmfrow=T, isolines=F, isostp=5, taxisStr='taxis',dataStr='lfp',posStr='channelPosition',nrmVal=NULL,cexMult=NULL,
                                   dotscl=.01,sclDotSize=T,dot.colors=jet.colors,style='map',directColorMap=F,normByTimePoint=F,contourAlpha=.5,...){
  lfp   <- tone[[dataStr]]
  if (length(dim(lfp))==3)
    lfp <- lfp[trialIndx,,]
  
  # check if these are ica or pca weights
  # if so, check if fft has been run on the time-courses of the ica/pca
  # if so, prepare spectra to be displayed instead of AEPs if showTraces is set to TRUE
  lfpTraces    <- lfp
  taxTraces    <- tone$taxis
  showFFT      <- FALSE

  if (is.null(chIndx))
    chIndx <- 1:dim(lfp)[2]
  
  chIndxTraces <- chIndx
  
  if (dataStr%in%c('weights','nrmWeights') & showTraces){
    showTraces <- F
    if ('spc'%in% names(tone)){
      showTraces <- T
      showFFT    <- T
      lfpTraces <- apply( tone$spc, c(2,3), mean )
      for (cx in 1:dim(lfpTraces)[2])
        lfpTraces[,cx] <- lfpTraces[,cx]/ max(abs(lfpTraces[,cx]))

      taxTraces <- tone$faxislog
      chIndxTraces <- eval.at
    }
  }
  
  
  electrodeIndex <- chIndx[chIndx<=dim(lfp)[2]]
  
  if (!is.null(taxisStr))
    taxis <- tone[[taxisStr]]
  if (is.null(taxisStr) | identical(taxisStr,'index'))
    taxis <- 1:dim(lfp)[1]
  
  if (style=='dots'){
    tpx <- which(tone$chStr %in% paste('ch',34:1025,sep='') )
    if (length(tpx)){
      tone$channelPosition$xpos[tpx] <- tone$channelPosition$ml1[tpx] + tone$channelPosition$Xgrid[tpx]/2
      tone$channelPosition$ypos[tpx] <- tone$channelPosition$dv1[tpx] + tone$channelPosition$Ygrid[tpx]/2
    }
    
    eex <- which(tone$chStr %in% paste('ch',1:33,sep='') )
    if (length(eex)>0){
      tone$channelPosition$xpos[eex] <- tone$channelPosition$ml1[eex] + tone$channelPosition$Xgrid[eex]/2
      tone$channelPosition$ypos[eex] <- tone$channelPosition$ap1[eex] + tone$channelPosition$Ygrid[eex]/2 - 2
    } 
  }
  
  xposE <- ifelse(flip,-1,1)*tone[[posStr]]$xpos[electrodeIndex]
  yposE <-                   tone[[posStr]]$ypos[electrodeIndex]
  
  xlim <- range(xposE)
  ylim <- range(yposE)
  
  if (showTraces){
    traceZero <- ylim[2] + 0.2 * diff(ylim)
    traceSpace <- 0.2*diff(ylim)
    #ylim <- ylim + c(0,ylim[2]+.5*diff(ylim) )
    ylim <- ylim + c(0,2*traceSpace )
    
  }
  
  if(is.null(eval.at))
    eval.at <- pretty(taxis,n=5)
  
  eval.at <- eval.at[which(eval.at>=min(taxis)&eval.at<=max(taxis))]
  #print(eval.at)
  
  ## convert time to index
  findClosestIndex <- function(t){return( which(abs(taxis-t)== min( abs(taxis-t)))) }
  eval.at.index <- sapply(eval.at,findClosestIndex)
  
  ## ================================== make the plot
  addTimePointMap <- function(tind,addMain=T, dotscl=.01, normByTimePoint=normByTimePoint){
    map <- NULL
    
    if(length(tind)>1)
      map  <- apply( lfp[tind, electrodeIndex ], 2, na.mean )
    
    if(length(tind)==1)
      map  <- lfp[tind, electrodeIndex ]
    
    if(!is.null(map)){
      
      xo   <- seq(min(xposE), max(xposE), length = 100)
      yo   <- seq(min(yposE), max(yposE), length = 100)
      imap <- akima::interp(xposE,yposE,map,xo=xo,yo=yo,linear=T, extrap=F)
      
      if(is.null(zlim))
        zlim <- c(-1,1)*na.max(c(abs(imap[[3]])))
      
      imapmat <- imap[[3]]%>=%zlim[1]
      imapmat <- imapmat%<=%zlim[2]
      
      imapmatunclip <- imap[[3]]
      
      meanTime <- signif( taxis[tind], 3 ) 
      main     <- paste(meanTime,'ms')
      if (!addMain)
        main <- ''
      
      image(xo,yo,imapmat,col=colmap(256), zlim=zlim, main=main,cex.main=2,xaxt='n',yaxt='n',bty='n',
            xlim=xlim, ylim=ylim, useRaster=TRUE)
      
      if (showTraces){
        if (is.null(trng))
          trng <- range(tone$taxis)
        
        tndx <- which(tone$taxis>=trng[1] & tone$taxis<=trng[2])
        if (is.null(nrmVal))
          nrmVal <- na.max(c(abs(lfp[tndx,chIndx])))
        
        sqzfun <- function(x){ squeeze(x, from.range=trng, to.range=xlim) }
        tax <- sqzfun(taxis)
        if (is.null(colorPallete))
          colorPallete <- colList(length(chIndx))
        
        for (cx in 1:length(chIndx))
          points( tax, traceZero + traceSpace * lfp[,chIndx[cx]]/nrmVal,type='l', col=colorPallete[cx] )
        points( rep(tax[tind],2), traceZero + c(-1,1)*traceSpace, col='red',type='l', lwd=2 )
      }
      
      if (isolines){
        levels <- pretty(2*zlim,18)
        levels <- seq(-30,to=30, by=2)
        levels <- seq( min(zlim),by=isostp,to=max(zlim) )
        if (identical(zlim, c(-1,1)) )
          levels <- seq(-1,1,by=0.1)
        
        contour(xo,yo, imapmatunclip,add=T, levels=levels, col='black')
      }
    }
    
    if(is.null(map))
      plot(c(-100),c(-1000), xlim=xlim, ylim=ylim, xaxt='n',yaxt='n',bty='n' )
    
    scaleX <- 0.1
    scaleY <- 0.1
    if(flat)
      scaleX <- 0.12
    
    #brainOutline <- showBrainOutline('leo_20130122',scaleX=scaleX, scaleY=scaleY,add=T)
    #xpos <- ep$ml1[ electrodeIndex ]
    #ypos <- ep$ap1[ electrodeIndex ]
    points(xposE,yposE,pch=20,...)
  }#end addTimePoint
  
  addTimePointDots <- function(tind,addMain=T,dotscl=.01, normByTimePoint=normByTimePoint){
    map <- NULL
    
    if(length(tind)>1)
      map  <- apply( lfp[tind, electrodeIndex ], 2, na.mean )
    
    if(length(tind)==1)
      map  <- lfp[tind, electrodeIndex ]
    
    meanTime <- signif( taxis[tind], 3 ) 
    main     <- paste(meanTime,'ms')
    
    if (identical(taxisStr,'icaNr'))
      main     <- paste('ICA #',meanTime)
    
    if ('what' %in% names(tone))
      main     <- paste(tone$what, ' #', meanTime)
    
    
    if (!addMain)
      main <- ''
    
    # old approach
    #map     <- map/max(abs(zlim))
    #cexMult <- (sqrt(abs(map)))%>=%(.05*dotscl)
    #colNdx  <- round( (100*(map+1)/2)%>=%1%<=%100  )
    
    # handel zlim NULL
    if (is.null(zlim))
      zlim <- c(-1,1) * quantile(abs(map), probs=.95)
    
    # normalize by time point if requested
    if (normByTimePoint)
      map    <- (map-mean(map))/sd(map)
    
    # new approach
    map    <- (map%>=%zlim[1])%<=%zlim[2]
    colNdx <- map
    if (!directColorMap)
      colNdx <- round(squeeze(map, from.range = zlim, to.range = c(1,100) ))
    
    map     <- map/max(abs(zlim))
    
    if (is.null(cexMult))
      cexMult <- (sqrt(abs(map)))%>=%(.05*dotscl)
        
    if (!sclDotSize)
      cexMult <- 1
    
    if(!is.null(map))
      plot(xposE,yposE, cex=dotscl*cexMult, col=dot.colors(100)[ colNdx ], 
           pch=20, xlim=xlim, ylim=ylim, xaxt='n',yaxt='n',bty='n', main=main, xlab='',ylab='' )
    
    #plot(xposE,yposE, cex=dotscl*sqrt(abs(map)), col=c('blue','red')[1+(map>0)], pch=20, xlim=xlim, ylim=ylim, xaxt='n',yaxt='n',bty='n' )
    
    if(is.null(map))
      plot(c(-100),c(-1000), xlim=xlim, ylim=ylim, xaxt='n',yaxt='n',bty='n' )
    
    if (addContour & 'GMC'%in%names(tone)){
      GMC_M <- tone$GMC
      addContourSlice <- function(Yx){
        xg <- tapply(tone$channelPosition$Xgrid, tone$channelPosition$Yx, mean)[paste(Yx)]
        yg <- tapply(tone$channelPosition$Ygrid, tone$channelPosition$Yx, mean)[paste(Yx)]
        
        as.numeric(names(GMC_M))
        CN <- GMC_M[[ paste(Yx) ]]
        Nc <- length(CN)
        for (nx in 1:Nc){
          points( CN[[nx]]$x + xg/2, CN[[nx]]$y + yg/2, type='l', col=alpha(rgb(0,0,0), contourAlpha), lwd=.5 )
        }
        
      }
      sapply( unique(tone$channelPosition$slice), addContourSlice)
      
      if ( sum(tone$channelPosition$Yx==0)>0 ){
        # add AP slices to the EEG contour 
        # x and y origin of the EEG overlay
        xg <- tapply(tone$channelPosition$Xgrid, tone$channelPosition$Yx, mean)[["0"]]
        yg <- tapply(tone$channelPosition$Ygrid, tone$channelPosition$Yx, mean)[["0"]]
        
        # meso=electro shaft positions for one slice 
        #MEP_ml <- tapply(tone$channelPosition$ml1, list(tone$channelPosition$Xx,tone$channelPosition$Yx),mean)
        MEP_ap <- tapply(tone$channelPosition$ap1, tone$channelPosition$slice, mean)
        for (sx in 2:length(MEP_ap))
          points( c(0,3.5), MEP_ap[c(sx,sx)]-2, type='l', lty=1, col=alpha(rgb(0,0,0), 0.5),lwd=.5 )
      }
    }
    
    if (showTraces){
    
      if (is.null(trng))
        trng <- range(taxTraces)
      
      tndx   <- which(taxTraces>=trng[1] & taxTraces<=trng[2])
      
      if (is.null(nrmVal))
        nrmVal <- na.max(c(abs(lfpTraces[tndx,chIndxTraces])))
      
      sqzfun <- function(x){ squeeze(x, from.range=trng, to.range=xlim) }
      tax <- sqzfun(taxTraces)
      if (is.null(colorPallete))
        colorPallete <- colList(length(chIndxTraces))
      
      if (showFFT){
        #for (cx in 1:length(chIndxTraces))
        #  points( tax, traceZero + traceSpace * lfpTraces[,chIndxTraces[cx]]/nrmVal,type='l', col='gray' )
        
        points( tax, traceZero + traceSpace * lfpTraces[,tind]/nrmVal,type='l', col='black' )
        # x-axis
        points( tax, traceZero + traceSpace * 0*lfpTraces[,chIndxTraces[1]],type='l', col='black' )
        # x-axis ticks
        xatick <- c(1/30, 1/20, 1/10, 1, 2.5, 5,10,20,40, 80)
        xatick <- xatick[ log(xatick)>=trng[1] & log(xatick)<=trng[2]]
        addTick <- function(tx){points( rep(sqzfun( log(xatick[tx]) ),2), traceZero + c(-.1,0)*traceSpace, col='black',type='l', lwd=1 )}
        disc <- sapply(1:length(xatick), addTick )
      }
      
      if (!showFFT){
        for (cx in 1:length(chIndxTraces))
          points( tax[tndx], traceZero + traceSpace * lfpTraces[tndx,chIndxTraces[cx]]/nrmVal,type='l', col=colorPallete[cx] )
        
        # eval at indicator bar
        points( rep(tax[tind],2), traceZero + c(-1,1)*traceSpace, col='red',type='l', lwd=2 )
        
        # y-axis
        points( rep(sqzfun(0),2), traceZero + c(-1,1)*traceSpace, col='black',type='l', lwd=1 )
        # x-axis
        points( tax, traceZero + traceSpace * 0*lfp[,chIndxTraces[1]],type='l', col='black' )
        
        # color-palette
        points( rep( sqzfun(0) ,50), traceZero +  traceSpace * seq(zlim[1]/nrmVal,to=zlim[2]/nrmVal,length=50), col=dot.colors(50)[1:50],type='p', pch=20, cex=1 )
        
        # y-axis labels
        text(  sqzfun(0), traceZero + traceSpace*zlim[1]/nrmVal,labels=paste(zlim[1]), adj=1.2 )
        text(  sqzfun(0), traceZero + traceSpace*zlim[2]/nrmVal,labels=paste(zlim[2]), adj=1.2 )
      }
    }
    
    
    if (!showTraces | (showTraces & showFFT) ){
      
      # color-palette
      minX <- min(xlim)
      minY <- min(ylim)
      rngX <- diff(range(xlim))
      rngY <- diff(range(ylim))
      
      tst <- density( colNdx[!is.na(colNdx)], from=1, to=100, n=50 )
      
      xc <- seq( minX+.01*rngX, to=minX+.2*rngX, length=50)
      xc <- seq( minX+.35*rngX, to=minX+.5*rngX, length=50)
      
      yo <- rep(minY - .02 * rngX, 50)
      yc <- yo + .15*rngY*tst$y/max(tst$y) 
      
      points( xc, yo, col='black', type='l' )
      points( xc, yc, col=dot.colors(50)[1:50],type='p', pch=20, cex=1 )
      
    }
  }
  
  
  if (style=='map')
    addTimePoint <- addTimePointMap
  
  if (style=='dots')
    addTimePoint  <- addTimePointDots
  
  if (is.null(mfrow))
    mfrow <- c(1+((length(eval.at)-1)%/%5) , min(5,length(eval.at)))
  
  ##mf = c(min(5,length(eval.at)), 1+((length(eval.at)-1)%/%5) )
  ratio <- 1.3
  if (showTraces)
    ratio <- 1.8
  
  if(eps)
    eps.file(fileName, s=1.75 ,mf=mfrow[c(2,1)], ratio=ratio, pdf=pdf)
  
  if (setmfrow)
    par(mfrow=mfrow, mar=c(1,1,2,1))
  
  addMain <- TRUE
  
  if( mfrow[1]>1 & mfrow[2]==1){
    par(mfrow=mfrow, mar=c(0,1,0,1))
    addMain <- F
  }
  
  disc <- lapply(eval.at.index,addTimePoint, addMain=addMain, dotscl=dotscl, normByTimePoint=normByTimePoint)
  
  if(eps)
    dev.off()
}

showSignalProcessorGIF <- function( tone, eval.at=seq(-1,by=1,to=80), gifName='test',style='dots', dotscl=2, zlim=c(-400,400), 
                                    showTraces = T, trng=c(-25,150), dot.colors=potential.colors,addContour=T,trialIndx=1, 
                                    mfrow=c(1,1),setmfrow=T,colList=function(n)(rep('gray',n)), 
                                    ani.width=600, ani.height=600, interval=.1 ){
  
  doOneTimePoint <- function(tndx,...){   showSignalProcessorMap(tone,eval.at=tndx, ...)}
  saveGIF(
    {
      disc <- sapply( eval.at, doOneTimePoint, style=style, dotscl=dotscl, zlim=zlim, 
                      showTraces = showTraces, trng=trng, dot.colors=dot.colors,addContour=addContour, 
                      mfrow=mfrow,setmfrow=setmfrow,colList=colList, trialIndx=trialIndx  )
    },
    movie.name = paste(pp(),gifName,".gif",sep=''), interval = interval, nmax = 30, 
    ani.width = ani.width, ani.height = ani.height, loop=1
  )
}


showSignalProcessorDPMap <- function(tone,IVstr='freqInd', IVbyStr='dummy', nrm=T, lognrm=T, trng=c(0,500),zlim=NULL,FUN=function(x){sqrt(sum(x^2))},...){
  
  tone$xpp$IV <- as.numeric(as.factor( tone$xpp[[IVstr]]))
  IVval       <- levels(tone$xpp$IV)
  NIV         <- length(unique(tone$xpp$IV))
  
  
  tone$xpp$dummy <- c(1)
  tone$xpp$IVby <- tone$xpp$dummy
  
  if (IVbyStr %in% names(tone$xpp))
    tone$xpp$IVby <- as.numeric( as.factor( tone$xpp[[IVbyStr]]) )
  IVbyval       <- levels(tone$xpp$IVby)
  NIVby         <- length(unique(tone$xpp$IVby))
  
  
  for (tx in 1:NIVby){
    
    tone$xpp$valTone <- tone$xpp$IVby == tx
    mntn <- avgSignalProcessor( tone, frml=~IV, subsExpr = expression( valTone ))
    tndx <- which(tone$taxis> trng[1] & tone$taxis<trng[2])
    
    gaCS          <- mntn
    gaCS$taxis    <- c( 1:(dim(mntn$lfp)[1]+1))
    gaCS$lfp      <- gaCS$lfp[1,1:length(gaCS$taxis),]
    
    for (sx in 1:NIV)
      gaCS$lfp[sx,]      <-  apply(mntn$lfp[sx,tndx,],c(2), FUN)# / apply(gaAll$lfp[,tndx,],c(2), function(x){sqrt(mean(x^2))}))
    
    gaCS$lfp[NIV+1,]     <- apply( gaCS$lfp[1:NIV,], 2, mean)
    
    if (nrm){
      if (lognrm){
        for (sx in 1:NIV)
          gaCS$lfp[sx,]      <-  log2(gaCS$lfp[sx,] / gaCS$lfp[ NIV+1 ,])
      }
      if (!lognrm){
        for (sx in 1:NIV)
          gaCS$lfp[sx,]      <-  gaCS$lfp[sx,] - gaCS$lfp[ NIV+1 ,]
      }
    }
    
    if (is.null(zlim))
      zlim <- c(-.9,.9) * max(abs(gaCS$lfp[1:NIV,]))
    
    showSignalProcessorMap( gaCS,eval.at=1:NIV,zlim=zlim, ... )
  }
}


reconstructICompSignalProcessor <- function( ica, cmpNr=NULL ){
  
  ## reconstruct the data X using the time-courses S and weights W specified in ic
  # by defult use all components, if compNr is not null, use only the onse specified here
  # INPUT
  # ic: output of call_fastICA
  # cmpNr: vector of integers specifying the components to be used
  
  
  # PRELIMINARIES ====
  Nica   <- dim(ica$weights)[1]
  Nchan  <- dim(ica$weights)[2]
  Ntpt   <- length(ica$taxis)
  Ntrial <- dim(ica$lfp)[1]
  
  if (is.null(cmpNr))
    cmpNr <- 1:Nica
  
  print(paste('reconstructing ICA components', cmpNr, sep=''))
  
  S <- three2twoDimSignalProcessor(ica$lfp)
  X <- S[,cmpNr,drop=F] %*% ica$weights[cmpNr,,drop=F]
  
  LFP <- two2threeDimSignalProcessor(X, Ntrial=Ntrial )
  return(LFP)
}

removeICompSignalProcessor <- function( tone, noiseComponent, show=F, check=F){
  
  # the old version of removeICAcomponent did remove components, but instead reconstructed
  # the signal using all other components. This can have a very different effect depending
  # on how much variance the components actually explain.
  #
  # The new version of the fuction fixes this problem
  
  # check if ica data frame exists
  if ('ica' %in% names(tone)){
    
    
    if ('noiseComponent'%in%names(tone)){
      noiseComponent <- setdiff(noiseComponent, tone$noiseComponent)
      if (length(tst)==0)
        break('noise component has already been removed')
    }
    
    # PRELIMINARIES ====
    ica    <- tone$ica
    Nica   <- dim(ica$weights)[1]
    Nchan  <- dim(ica$weights)[2]
    Ntpt   <- length(ica$taxis)
    Ntrial <- dim(ica$lfp)[1]
    
    if (!identical( dim(tone$lfp), c(Ntrial, Ntpt, Nchan)))
      break('dimension mismatch in removeICompSignalProcessor')
    
    print(paste('removing ICA ', noiseComponent, sep=''))
    keepComponent <- which( ! (1:Nica %in% noiseComponent))  
    
    # creates an SP-like data frame with the to-be-removed components reconstructed
    Nsp <- reconstructICompSignalProcessor( ica, cmpNr=noiseComponent )
    tone$lfp <- tone$lfp - Nsp$lfp
    
    old <- F
    if (old){
      # first convert 3D representation of component time-courses into two-dimensional representation (trial*time)*component
      tmp1      <- aperm(ica$lfp, c(3,2,1))
      dim(tmp1) <- c(Nica, dim(ica$lfp)[1]*dim(tmp1)[2] )
      S        <- t(tmp1)
      
      # multiply component time-courses with weights-matrix to recover channel time-series
      #X    <- S[,keepComponents] %*% ica$weights[keepComponents,]
      Xnot <- S[,noiseComponent,drop=F] %*% ica$weights[noiseComponent,,drop=F]
      
      # convert back to 3D representation
      tmp1      <- aperm(Xnot, c(2,1))
      dim(tmp1) <- c(Nchan, Ntpt, Ntrial )
      #ica       <- list(lfp=aperm(tmp1,c(3,2,1)),  xpp=lfp$xpp[trlst,], weights=a$A, K=a$K, W=a$W,srt=srt, daxis=daxis, taxis=lfp$taxis)
      tone$lfp     <- tone$lfp - aperm(tmp1,c(3,2,1))
    }
    
    
    if ('noiseComponent'%in%names(tone))
      tone$noiseComponent <- c(tone$noiseComponent, noiseComponent)
    
    if (!'noiseComponent'%in%names(tone))
      tone$noiseComponent <- noiseComponent
    
    
    return(tone)
  }
}

applyPCAsignalprocessor <- function(tone, pca, whatStr='pca'){
  # tone is a signalProcessor data frame
  # pca is a pca result from call_fastPCA
  
  ica <- tone[ setdiff( names(tone), c('lfp','cwt','spc','logspc','pca','ica','tfft','faxis','faxislog')) ]
  
  A       <- pca$weights
  LFP     <- three2twoDimSignalProcessor( tn$lfp )
  RFP     <- t(A%*%t(LFP)) #%*%A

  # normalize to standard deviation of 1
  sdRFP <- array( rep( apply(RFP,2,sd), each=dim(RFP)[1]), dim(RFP)) 
  RFP   <- RFP / sdRFP
  
  ica$lfp <- two2threeDimSignalProcessor(RFP,  Ntrial=dim(tn$lfp)[1])
  
  ica$icaPosition <- pca$icaPosition
  ica$icaNr       <- pca$icaNr
  ica$weights     <- pca$weights
  
  tone[[whatStr]] <- ica
  
  return(tone)
}

ch2regSignalProcessor <- function(tone, subs=NULL, seedList=NULL, what='pca', Nica=2, mxNica=Inf, mnNica=1, chBlkN=8, demean=F){
  # creates a new SP-like data frame tone$reg that shares properties with $ica or $pca
  # however, it calculate separate pca/ica for each subset of channels in each elemet of seedList
  # It then selects a certain number of the pca's or ica's to include $reg as representative elemets of the seed region
  # it also creates a weight matrix that resembles the $pca$weight and can be used as a projection to extract regional time-courses 
  # for new data sets.
  
  # seedList contains either numeric or character vectors. If they are numeric, it is assmed that they are indexes into $lfp
  # if they are characters, they are assumed to be of the type 'ch45', in which case the approppiate index will be extracted
  # from tone$chStr
  
  Nreg    <- length(seedList)
  if (is.null(names(seedList)))
    names(seedList) <- paste('reg',1:Nreg,sep='')
  
  getChNdx    <- function(cs){ which(tone$chStr %in% cs)}
  
  if (is.null(subs))
    subs <- rep(T, dim(tone$lfp)[1] )
  
  extractFcn <- id
  if (is.character(seedList[[1]]) )
    extractFcn <- getChNdx
  
  NchList  <- sapply(1:length(seedList), function(sx){length(seedList[[sx]]) })
  
  if (all(Nica > 0) ){ # use Nica unless there are fewer than Nica channels in a region
    if (length(Nica)==1)
      NicaList <- apply(  cbind( rep( Nica, length(seedList) ),  NchList), 1, min)
    
    if (length(Nica)==length(seedList))
      NicaList <- Nica %>=%1
    
  }
  
  if (all(Nica < 0) )# use a fractino of the channels determined by Nica: Nchan * Nica 
    NicaList <-  round( -1*Nica*NchList ) %>=% mnNica %<=% mxNica
  
  if ( length(unique(sign(Nica)))>1 )
    stop( 'Nica elements must either all be larger than zero or less than zero. Can t yet mix signs=.' )
  
  
  # compute the pca/ica for each seed region with the requested number of components
  getRegions     <- function(sx, tn=tone){
    print(names(seedList)[sx])
    reg        <- subs.data( tn, rep(T,dim(tn$xpp)[1]), chIndx = getChNdx( seedList[[sx]] ) )
    reg$chIndx <- getChNdx( seedList[[sx]] )
    reg$regStr <- paste( names(seedList)[sx], 1:NicaList[sx], sep='#' )
    reg$ica    <- call_fastICA(reg, Npca=NicaList[sx], show=F, restoreAll=T, dosrt=F,subset=subs, demean = demean)
    reg$pca    <- recode_ICA(lfp=reg$lfp, ica=reg$ica, what='pca')
    return(reg)  
  }
  regList        <- lapply(1:Nreg, getRegions, tn=tone )   
  names(regList) <- names(seedList)
  
  # combine all of the regional pca/ica into a large data frame that resembles the $ica or $pca data frames
  regStr <- c( unlist( sapply(1:length(regList), function(sx){regList[[sx]]$regStr} ) ))

  pca     <- regList[[1]]$pca
  pca$lfp <- regList[[1]]$pca$lfp
  
  pca$weights <- array(0, c(length(regStr),dim(tone$lfp)[3]  ) )
  pca$lfp     <- array(0, c(dim(tone$lfp)[1],dim(tone$lfp)[2], length(regStr) ) )
  
  
  
  pca$chStr                <- tone$chStr
  pca$channelPosition      <- tone$channelPosition
  pca$icaPosition          <- data.frame(xpos=1+(0:(length(regStr)-1))%/%chBlkN, ypos=1+(0:(length(regStr)-1))%%chBlkN,  zpos=rep(0,length(regStr)) )
  pca$icaPosition$chStr    <- regStr  
  pca$icaNr                <- 1:sum(NicaList)
  
  # same for ica:
  ica         <- pca
  
  for (rx in 1:Nreg){

    tochx            <- 1:NicaList[rx] + cumsum(c(0,NicaList))[rx]
    
    for(cx in 1:NicaList[rx]){
      # set the sign of the component timecourses:
      # each region has a reference channel, always the first seed listed
      # find the sign of the weight of that channel for this component
      # multiply the time-course and weigths with that sign.
      
      rfx <- which( regList[[rx]]$chStr %in% seedList[[rx]][1] )# this is index of the first channel in the list
      
      # for the pca:
      thisSignP                      <- sign( regList[[rx]]$pca$weights[cx,rfx] )
      regList[[rx]]$pca$lfp[,,cx]    <- thisSignP * regList[[rx]]$pca$lfp[,,cx]
      regList[[rx]]$pca$weights[cx,] <- thisSignP * regList[[rx]]$pca$weights[cx,]
      
      # for the ica:
      thisSignP                      <- sign( regList[[rx]]$ica$weights[cx,rfx] )
      regList[[rx]]$ica$lfp[,,cx]    <- thisSignP * regList[[rx]]$ica$lfp[,,cx]
      regList[[rx]]$ica$weights[cx,] <- thisSignP * regList[[rx]]$ica$weights[cx,]

    }
      

    pca$lfp[,,tochx] <- regList[[rx]]$pca$lfp 
    ica$lfp[,,tochx] <- regList[[rx]]$ica$lfp 
    
    pca$weights[ tochx, regList[[rx]]$chIndx ] <- regList[[rx]]$pca$weights
    ica$weights[ tochx, regList[[rx]]$chIndx ] <- regList[[rx]]$ica$weights
  }
  
  #showSignalProcessor( pca, trng=c(-50,250), frml=~ampInd)
  #image( cor( t(pca$weights) ) )
  
  tone$pca <- pca
  tone$ica <- ica
  
  return(tone)
}



# tmeanPSTH <- function(psth, FUN=mean){
#   #print(dim(psth))
#   ## little helper function that returns valid mean psths even
#   ## if the number of repetitions is one or zero
#   res <- array(NA,c(dim(lfp)[c(2,3)]))  
#   
#   if(is.array(psth)&length(dim(psth))==3 )
#     res <- apply(psth,c(2,3),FUN)
#   
#   if( is.array(psth) & length(dim(psth))==2 & identical(FUN,mean))
#     res <- psth
#   
#   return(res)
# }


tmeanPSTH <- function(psth, FUN=mean){
  #print(dim(psth))
  ## little helper function that returns valid mean psths even
  ## if the number of repetitions is one or zero
  
  Ntrial <- dim(lfp)[1]
  Ntst   <- dim(lfp)[2]
  Nchan  <- dim(lfp)[3]
  
  res <- array(NA,c(Ntst,Nchan))  
  
  if(is.array(psth)&length(dim(psth))==3 )
    res <- apply(psth,c(2,3),FUN)
  
  if( is.array(psth) & length(dim(psth))==2 & identical(FUN,mean))
    res <- psth
  
  return(res)
}


downsampleSignalProcessor <- function( tone, newTaxis=NULL, widthMS=5  ){
  
  sr <- 1000 * mean( diff(tone$taxis) )
  if (is.null(newTaxis) )
    newTaxis <- seq( min(tone$taxis)+widthMS*sr/1000, by=widthMS*sr/1000,to=max(tone$taxis - widthMS*sr/1000) ) 
  
  lfp <- array(data = 0, dim=c(dim(tone$lfp)[3],dim(tone$lfp)[1],length(newTaxis) ))
  eval.at <- which( tone$taxis %in% newTaxis)
  for(ch in 1:dim(tone$lfp)[3]){
    lfp[ch,1:dim(tone$lfp)[1],1:length(newTaxis)] <- t(sapply( 1:(dim(tone$lfp)[1]),getEnvStmMN, widthMS=widthMS, eval.at=eval.at,ttn=tone, chndx=ch) )
  }
  
  dntone       <- tone
  dntone$taxis <- newTaxis
  dntone$lfp   <- aperm(lfp, c(2,3,1) )
  return(dntone)
}


padSignalProcessor <- function(tone, Npad=100, what=0){
  tnpad         <- tone
  tnpad$lfp     <- zbind(tone$lfp, array(what, c(dim(tone$lfp)[1], Npad, dim(tone$lfp)[3]) ) , dimension = 2)
  SR <- mean(diff(tone$taxis))
  tnpad$taxis   <- c(tone$taxis, max(tone$taxis)+ SR*(1:Npad) )
  
  return(tnpad)
}

corSignalProcessor  <- function(tone, Ntaper=NULL, type='correlation', chIndx=NULL){
  # calculcates cross correlation of the refChan to all other channels in the array
  
  Ntrial <- dim(tone$lfp)[1]
  Ntst   <- dim(tone$lfp)[2]  
  Nchan  <- dim(tone$lfp)[3]
  
  if (is.null(chIndx))
    chIndx <- 1:Nchan
  
  # calculate the taper window and multiply data with it
  if ( is.numeric(Ntaper) ){
    taper    <- cosine.window( Ntst, Ntaper )
    tmp      <- array( taper, c(Ntst, Ntrial, Nchan)  )
    tapermat <- aperm(tmp,c(2,1,3))
    tone$lfp <- tone$lfp * tapermat
  }
  
  LFP         <- three2twoDimSignalProcessor( tone$lfp )
  ccMat       <- cor(LFP)
  tone$corMat <- ccMat
  return(tone)
}
acSignalProcessor   <- function(tone, lag.max=200, Ntaper=50, jitter=F, type='correlation', chIndx=NULL){
  # calculcates cross correlation of the refChan to all other channels in the array
  
  Ntrial <- dim(tone$lfp)[1]
  Ntst   <- dim(tone$lfp)[2]  
  Nchan  <- dim(tone$lfp)[3]
  
  if (is.null(chIndx))
    chIndx <- 1:Nchan
  
  # calculate the taper window and multiply data with it
  if ( is.numeric(Ntaper) ){
    taper    <- cosine.window( Ntst, Ntaper )
    tmp      <- array( taper, c(Ntst, Ntrial, Nchan)  )
    tapermat <- aperm(tmp,c(2,1,3))
  }
  tone$lfp <- tone$lfp * tapermat
  
  LFP <- three2twoDimSignalProcessor( tone$lfp )
  
  # the cross-correlation function for a specific channel pair
  getccf<- function(cx, refChn=refChn, lag.max=lag.max){ 
    tmpac <- ccf( LFP[,cx], LFP[,cx], lag.max=lag.max, type=type, main = "", plot=F)
    return(c(tmpac$acf))
  }
  #tst <- getccf(1, refChn=23,lag.max=lag.max)
  
  # get cross correlation for each channel pair
  ccMat <- sapply(chIndx, getccf, refChn=refChan, lag.max=lag.max)
  
  ccDF                 <- tone
  ccDF$xpp             <- ccDF$xpp[1,]
  ccDF$lfp             <- ccMat
  ccDF$taxis           <- -lag.max:lag.max
  dim(ccDF$lfp)        <- c(1,dim(ccMat))
  
  return(ccDF)
}
ccSignalProcessor   <- function(tone, seedChIndx=1, lag.max=200, Ntaper=50, jitter=F, type='correlation', chIndx=NULL){
  # calculcates cross correlation of the refChan to all other channels in the array
  
  Ntrial <- dim(tone$lfp)[1]
  Ntst   <- dim(tone$lfp)[2]  
  Nchan  <- dim(tone$lfp)[3]
  
  if (is.null(chIndx))
    chIndx <- 1:Nchan
  
  # calculate the taper window and multiply data with it
  if ( is.numeric(Ntaper) ){
    taper    <- cosine.window( Ntst, Ntaper )
    tmp      <- array( taper, c(Ntst, Ntrial, Nchan)  )
    tapermat <- aperm(tmp,c(2,1,3))
  }
  tone$lfp <- tone$lfp * tapermat
  
  LFP <- three2twoDimSignalProcessor( tone$lfp )
  
  
  if (length(seedChIndx)==1)
    refSig <-LFP[,seedChIndx]
  
  if (length(seedChIndx)>1)
    refSig <- apply( LFP[,seedChIndx],1,mean )
  
  # the cross-correlation function for a specific channel pair
  getccf<- function(cx, refChn=refChn, lag.max=lag.max){ 
    tmpac <- ccf( LFP[,cx], refSig, lag.max=lag.max, type=type, main = "", plot=F)
    return(c(tmpac$acf))
  }
  #tst <- getccf(1, refChn=23,lag.max=lag.max)
  
  # get cross correlation for each channel pair
  ccMat <- sapply(chIndx, getccf, refChn=refChan, lag.max=lag.max)
  
  ccDF                 <- tone
  ccDF$xpp             <- ccDF$xpp[1,]
  ccDF$lfp             <- ccMat
  ccDF$taxis           <- -lag.max:lag.max
  dim(ccDF$lfp)        <- c(1,dim(ccMat))
  
  return(ccDF)
}
pcorSignalProcessor <- function(tone, Ntaper=NULL, type='correlation', returnP=FALSE,
                                seedChIndx=NULL, chIndx=NULL, pcorChIndx=NULL, trIndx=NULL){
 
  if (is.null(seedChIndx))
    stop('Specify at least one seed channel.')
  
  if (is.null(pcorChIndx)){
    warning('without a seed channel, corSignalProcessor is probably a better choice...')
  }
  
  Ntrial <- dim(tone$lfp)[1]
  Ntst   <- dim(tone$lfp)[2]  
  Nchan  <- dim(tone$lfp)[3]
  
  if (is.null(chIndx))
    chIndx <- 1:Nchan
  
  if (is.null(trIndx))
    trIndx <- 1:Ntrial

  # this uses the first principal component of the seed region
  # convert the list of channels back into strings of the type: ch45
  if (length(seedChIndx)> 1){
    seedList <- list( sd=tone$channelPosition$chStr[seedChIndx] )
    rg      <- ch2regSignalProcessor( tone, seedList = seedList, Nica=1)$pca
  
    # make sure that strongest weight is positive
    #if (max(abs(rg$weights))== max(-1*rg$weights)){
    if (mean(abs(rg$weights)) < 0){
      rg$weights <- -1 * rg$weights
      rg$lfp <- -1 * rg$lfp
    }
  }
  
  # calculate the taper window and multiply data with it
  if ( is.numeric(Ntaper) ){
    taper    <- cosine.window( Ntst, Ntaper )
    tmp      <- array( taper, c(Ntst, Ntrial, Nchan)  )
    tapermat <- aperm(tmp,c(2,1,3))
    tone$lfp <- tone$lfp * tapermat
     
    rg$lfp <- rg$lfp * tapermat[1,,]
  }
  
  
  LFP         <- three2twoDimSignalProcessor( tone$lfp[trIndx,,,drop=F] )
  
  if (length(seedChIndx)==1)
    refSig <-LFP[,seedChIndx]
  
  # this uses the average activity in the seed region
  #if (length(seedChIndx)>1)
  #  refSig <- apply( LFP[,seedChIndx],1,mean )
  
  # this uses the first principal component
  if (length(seedChIndx)>1){
    RG         <- three2twoDimSignalProcessor( rg$lfp[trIndx,,,drop=F] )
    refSig <- RG[,1]
  }
  
  getPcorChannel <- function(chx){
    
    sameseed <- F    
    if (length(seedChIndx)==1 & chx == seedChIndx[1] ){
      rs <- 1#cor( cbind(refSig, LFP[,chx]) )
      if (returnP)
        rs <- 0
      
      sameseed <- T
    }
    
    if (!sameseed){
      #tPcorx <- setdiff(pcorChIndx,chx)
      sameref <- chx %in% pcorChIndx
      rs <- 0
      if (returnP)
        rs <- 1
      
      if (!sameref){
         tmp <- pcor( cbind(refSig, LFP[,c(chx,pcorChIndx)] ))
         rs  <- tmp$estimate[1,2]
         if (returnP)
           rs  <- tmp$p.value[1,2]
      }
    }
  
    return(rs)#$estimate[1,2])
    }
  pcorVec      <- sapply(chIndx, getPcorChannel)
  tone$pcorVec <- t(pcorVec)
  return(tone)
}


resampleSignalProcessor <- function(tone, newTaxis){
  lfp   <- tone$lfp
  taxis <- tone$taxis
  
  newLFP <- array(c(0), c(dim(lfp)[1], length(newTaxis), dim(lfp)[3] ))
  
  for (chx in 1:dim(tone$lfp)[3]){
    for (tx in 1:dim(tone$lfp)[1]){
      interpFUN   <- approxfun(taxis, lfp[tx,,chx])
      newLFP[tx, ,chx] <- interpFUN( newTaxis )
    }
  }
  tone$taxis <- newTaxis
  tone$lfp <- newLFP
  return(tone)
}



trAOVSignalProcessor <- function(tone,frml=~dummy, lmName='',alpha=0.05,
                                 chIndx=NULL,subsExpr=expression(valTone),dataStr='lfp',taxisStr='taxis',
                                 fix=FALSE, trng=c(-25,250)){
  ## preliminaries                         ====
  tone$xpp$dummy <- sample(2,dim(tone$xpp)[1], replace=T)
  if ('rejected'%in%names(tone))
    tone$rejected <- c(1)
  
  if (!is.null(subsExpr))
    subs <- with(tone$xpp, eval(subsExpr) )
  
  if (is.null(subsExpr))
    subs <- rep(T, dim(tone$lfp)[1] )
  
  subs[is.na(subs)] <- FALSE
  
  teps  <- 10e-10 # fudge factor to account for inaccurate time-axis
  ttn   <- subs.data(tone,subs)
  tndx  <- which(ttn[[taxisStr]]>=(trng[1]-teps) &  ttn[[taxisStr]]<=(trng[2]+teps)) 
  taxis <- ttn[[taxisStr]][tndx]
  lfp   <- ttn[[dataStr]][,tndx,,drop=F]
  IV    <- model.frame(frml, data = ttn$xpp )
  
  Ntrial <- dim(lfp)[1]
  Ntst   <- dim(lfp)[2]
  Nchan  <- dim(lfp)[3]
  
  mnIVtable <- min( unlist( apply(IV,2,table) ))
  doAVG     <- mnIVtable < 2  # need at least 
  if (mnIVtable < 2)
    warning(paste('there may not be enough trials for a meaningful aov analysis'))
  
  ## create data from fro run_trAOV_SP     ====
  # create a reduced data frame for run_trAOV_SP, oavIn excludes unwanted trials, and time-points
  aovIn       <- ttn
  aovIn$xpp   <- IV
  aovIn$lfp   <- lfp
  aovIn$taxis <- taxis
  
  
  ## normalize numeric predictor variables ====
  for (vx in 1:dim(aovIn$xpp)[2]){
    if(is.numeric(aovIn$xpp[,vx])){
      if (na.sd(aovIn$xpp[,vx])>teps )
        aovIn$xpp[,vx]           <- (aovIn$xpp[,vx]           - na.mean(aovIn$xpp[,vx] ))/na.sd(aovIn$xpp[,vx])
    }
  }
  
  
  ## define the test that is handed to run_trAOV_SP
  this.test <- function(subinput){
    
    if(max(subinput$xpp$DV)==0 && min(subinput$xpp$DV)==0) ## not so nice trick to avoid a singular system
      subinput$xpp$DV[1] <- 1
    
    res         <- list()
    res$aovFull <- aov(DV ~ 1+., data=subinput$xpp)
    res$car     <- Anova(res$aovFull,type=2)
    
    lm0           <- aov(DV ~ -1, data=subinput$xpp)
    res$lm1       <- aov(DV ~  1, data=subinput$xpp)
    res$blcar     <- anova(lm0,res$lm1)
    return( res )
  }
  
  ## run the test using run_trAOV_SP       ====
  res <- run_trAOV_SP(aovIn,this.test,frml=frml,show=show,alpha=alpha,lmName=lmName)
  
  ## save the results ====
  #save(res, file=filename)
  
  return(res)
  
}




run_trAOV_SP <- function(input,thisTest,frml=frml,show=T,alpha=.05,
                         rw=NULL,subs=NULL, lmName=''){
  
  if (lmName == 'default')
    lmName <- frml2str(frml)
  
  ## paths and file names ====
  saveDir <- paste('/Volumes/EEGCyno/AVG/',input$project,'/',input$exp,'/',input$epoch,'/trAOV/' ,sep='')
  if (!lmName=='')
    saveDir <- paste('/Volumes/EEGCyno/AVG/',input$project,'/',input$exp,'/',input$epoch,'/trAOV_',lmName,'/' ,sep='')
  
  
  saveDirPar <- paste('/Volumes/EEGCyno/AVG/',input$project,'/',input$exp,'/',input$epoch,'/trAOV_parm/' ,sep='')
  if (!lmName=='')
    saveDirPar <- paste('/Volumes/EEGCyno/AVG/',input$project,'/',input$exp,'/',input$epoch,'/trAOV_parm_',lmName,'/' ,sep='')
  
  
  if (!file.exists(saveDir))
    system( paste('mkdir -p ', saveDir,sep=''), intern=T )
  
  if (!file.exists(saveDirPar))
    system( paste('mkdir -p ', saveDirPar,sep=''), intern=T )
  
  
  eval.at <- 1:length(input$taxis)
  
  ## loop by channel ====
  runChannel <- function(chx){
    
    filename <- input$chStr[chx]
    print(filename)
    
    do.it <- function(tndx){
      input$xpp$DV <- input$lfp[,tndx,chx]
      res          <- thisTest(input)
      return(res)
    }
    
    result <- tapply(eval.at,eval.at,do.it)
    
    Ntst <- length(result)
    Nvar <- dim(result[[1]]$car)[1]
    
    if(is.null(rw))
      rw <- dim(result[[1]]$car)[2]
    
    
    ## separate p-value analysis for the intercept =====
    if( exists('blcar',where=result[[1]]) ){
      
      extract.P <- function(n,rw=6,vStr='car'){return(result[[n]][[vStr]][,rw][2])}  
      blP <- t(matrix( unlist( tapply(eval.at,eval.at,extract.P,rw=6,vStr='blcar')), c(1,Ntst))) ## Pvalue
      blF <- t(matrix( unlist( tapply(eval.at,eval.at,extract.P,rw=5,vStr='blcar')), c(1,Ntst))) ## Fvalue
      blS <- t(matrix( unlist( tapply(eval.at,eval.at,extract.P,rw=3,vStr='blcar')), c(1,Ntst))) ## degrees of freedom
      
      extract.par <- function(n){return(result[[n]]$lm1$coeff)}
      blPar       <- t(matrix( unlist( lapply(eval.at,extract.par)), c(1,Ntst)))
      
    }
    ## p-values for all the other variables ====
    extract.P <- function(n,rw=6,vStr='car'){
      return(result[[n]][[vStr]][,rw])
    }
    
    Pmat   <- t(matrix( unlist( tapply(eval.at,eval.at,extract.P,rw=rw))  , c(Nvar,Ntst)))
    Fmat   <- t(matrix( unlist( tapply(eval.at,eval.at,extract.P,rw=rw-1)), c(Nvar,Ntst)))
    SSQmat <- t(matrix( unlist( tapply(eval.at,eval.at,extract.P,rw=1))   , c(Nvar,Ntst)))
    
    degfred    <- extract.P(1,rw=2)
    effectName <- dimnames(result[[1]]$car)[[1]]
    
    ## extract parameter estimates ====
    extract.par <- function(n){
      return(result[[n]]$aovFull$coeff)
    }
    Npar     <- 1 + sum(degfred[-length(degfred)]%>=%1)  ## intercept plus each degree of freedom
    parMat   <- t(matrix( unlist( lapply(1:length(eval.at),extract.par)), c(Npar,Ntst)))
    coeffName <- names(extract.par(1))
    
    
    if (exists('blcar',where=result[[1]])){
      Pmat       <- cbind(blP, Pmat)
      Fmat       <- cbind(blF, Fmat)
      SSQmat     <- cbind(blS, SSQmat)
      parMat[,1] <- blPar
      degfred    <- c(1,degfred)
      effectName <- c('Intercept', effectName)
      #coeffName  <- c('Intercept',coeffName)
    }
    
    
    dimnames(Pmat)[[2]] <- effectName
    dimnames(Pmat)[[1]] <- eval.at
    
    dimnames(Fmat)   <- dimnames(Pmat)
    dimnames(SSQmat) <- dimnames(Pmat)
    
    dimnames(parMat)[[1]] <- eval.at
    dimnames(parMat)[[2]] <- coeffName
    
    Pmat[,dim(Pmat)[2]] <- 1
    ##Pmat   <-   Pmat[,-dim(Pmat)[2]]
    ##Fmat   <-   Fmat[,-dim(Fmat)[2]]
    ##SSQmat <- SSQmat[,-dim(SSQmat)[2]]
    
    final <- list(allign=input$allign,eval.at=eval.at,Pmat=t(Pmat),Fmat=t(Fmat), SSQmat=t(SSQmat),
                  degfred=degfred,parMat=t(parMat),filename='show.tres.glm.eps')
    
    #save(final, file=paste( rda(),fileName,'.rda', sep=''))
    
    return(final)
  }
  Nchan   <- dim(input$lfp)[3]
  Ntst    <- dim(input$lfp)[2]
  
  resList <- lapply( 1:Nchan, runChannel)
  Nvar    <- dim(resList[[1]]$Pmat)[1]
  
  ## create signalProcessor style result for p-values ====
  effectName  <- dimnames(resList[[1]]$Pmat)[[1]]
  aovxpp.P    <- data.frame( effect=effectName, df=resList[[1]]$degfred  )
  extractPMat <- function(chx){ return(resList[[chx]]$Pmat) }
  Pmat        <- matrix( unlist( tapply(1:Nchan,1:Nchan,extractPMat)), c(Nvar,Ntst,Nchan))
  dim(Pmat)   <- c(Nvar,Ntst,Nchan)
  
  ptn <- list(xpp=aovxpp.P, lfp=Pmat, taxis=input$taxis,chStr=input$chStr, epoch=input$epoch, exp=input$exp, channelPosition=input$channelPosition)
  
  ## save signalProcessor style results for p-values ====
  xpp   <- ptn$xpp
  taxis <- ptn$taxis
  save(xpp, file=paste(saveDir,'xpp.rda',sep='') )
  save(taxis, file=paste(saveDir,'taxis.rda',sep='') )
  saveChannel <- function(cx){
    chdata <- ptn$lfp[,,cx]
    dim(chdata) <- c( dim(ptn$lfp)[1:2],1)
    save(chdata, file=paste(saveDir,ptn$chStr[cx],'.rda',sep=''))
  }
  disc <- lapply(1:dim(ptn$lfp)[3], saveChannel)  
  
  #if(show)
  #  try(show.tres.glm(final,alpha=alpha))
  
  ## Create signalProcessor stype results for parameter estimates ====
  
  aovxpp.par <- data.frame( effect=dimnames(resList[[1]]$parMat)[[1]] )
  extractParMat <- function(chx){ return(resList[[chx]]$parMat) }
  Nbeta <- dim(aovxpp.par)[1]
  parMat        <- matrix( unlist( tapply(1:Nchan,1:Nchan,extractParMat)), c(Nbeta,Ntst,Nchan))
  dim(parMat)   <- c(Nbeta,Ntst,Nchan)
  
  partn <- list(xpp=aovxpp.par, lfp=parMat, taxis=input$taxis,chStr=input$chStr, epoch=input$epoch, exp=input$exp, channelPosition=input$channelPosition)
  
  ## save signalProcessor style results for beta-values ====
  xpp   <- partn$xpp
  taxis <- partn$taxis
  save(xpp, file=paste(saveDirPar,'xpp.rda',sep='') )
  save(taxis, file=paste(saveDirPar,'taxis.rda',sep='') )
  saveChannel <- function(cx){
    chdata <- partn$lfp[,,cx]
    dim(chdata) <- c( dim(partn$lfp)[1:2],1)
    save(chdata, file=paste(saveDirPar,partn$chStr[cx],'.rda',sep=''))
  }
  disc <- lapply(1:dim(partn$lfp)[3], saveChannel)  
  
  
  return(ptn)
}

# 
# # ------------ frequency
# freqList          <- dfline$minOct*exp(log(2)*seq(0,by=dfline$stepOct, length=dfline$NFrq))
# tone$xpp$freqHz   <- freqList[ tone$xpp$freqInd ]
# tone$xpp$freqOct  <- log2( freqList[ tone$xpp$freqInd ])
# 
# # -------------- amplitude
# ampList           <- dfline[grep(  'dB' , names(dfline))] 
# tone$xpp$ampdB    <- c(ampList[ tone$xpp$ampInd ])
# 
# # sequence number
# tonePerTrial   <- table(tone$xpp$trialNr)
# tone$xpp$seqNr <- c(unlist( sapply(tonePerTrial, function(n){1:n}) ))
# 
# 

