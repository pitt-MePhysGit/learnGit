batchFLA_ampOddClick <- function(projectName, linkFileName, ndx){
  
  #projectName<-'ampOddClick'
  #linkFileName<-'ampOddclickVProbe_link.txt'
  #ndx<-183
  #ndx<-170
  
  source('~/Dropbox/r-compile/signalProcessoR/R/signalProcessoR.R')
  sourceDir( paste( '~/Dropbox/',projectName,'/r-files/FLA/',sep='')  )
  source('~/Dropbox/ampOddClick/r-files/functions/FLA_VProbe_STPSP.R')
  
  df          <- prepDFsignalProcessor(projectName, linkFileName)
  dfline      <- df[ndx,]
  print(dfline)
  
  ## determine which functions to run ====
  run_FLAVProbe      <- T
  run_spectroLaminar <- F
  
  vProbeSession <- grepl('v', dfline$ElectrodeConfig) | grepl('d', dfline$ElectrodeConfig) | grepl('t', dfline$ElectrodeConfig) | grepl('z',  dfline$ElectrodeConfig)
  
  
  ## start the requested functions ====
  
  if (run_FLAVProbe & vProbeSession){
    trda     <- function(){'/Volumes/EEGCyno/EEG/ampOddClick/'}
    baseDir  <- '~/Dropbox/'
    FLA_VProbe_STPSP(1, dfline, baseDir=baseDir,trda=trda,ind=1)
  }
  
  if (run_spectroLaminar & vProbeSession){
    trda     <- function(){'/Volumes/EEGCyno/EEG/ampOddClick/'}
    baseDir  <- '~/Dropbox/'
    function_VProbe_SpectroLaminar(1, dfline, baseDir=baseDir,trda=trda,ind=1)
  }
  

}