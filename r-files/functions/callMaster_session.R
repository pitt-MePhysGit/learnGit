callMaster_session <- function(dfline, whatList=c('raw'), fltxtn=''){
  disc <- lapply( whatList, callMaster_session_what, dfline=dfline, fltxtn=fltxtn  )
}


callMaster_session_what <- function(what, dfline=NA, fltxtn=''){
  
  print(paste( 'callMaster_session_what', what))
  
  ## define paths
  animal <- strsplit(dfline$session,'_')[[1]][1]
  
  print('loading data')
  chStr <- paste('ch',1:21,sep='')
  if (animal%in%c('Jesse','Walter','Sam'))
    chStr <- paste('ch',1:33,sep='')
    
  print(chStr)
  print(fltxtn)
  maxp2p    <- 450
  maxPowEOG <- Inf
  
  if(animal=='Rockey' & what %in% c('raw','raw5k','mlp','mbp','bsp'))
    maxPowEOG <- 5000
    
  if(animal=='Walter' & what %in% c('raw','raw5k','mlp','mbp','bsp'))
    maxPowEOG <- 5000
  
  if(animal=='Jesse' & what %in% c('raw','raw5k','mlp','mbp','bsp'))
    maxPowEOG <- 2000
  
  if(animal=='Sam' & what %in% c('raw','raw5k','mlp','mbp','bsp'))
    maxPowEOG <- 15000
  
  
  ## load data      
  tone   <- loadAmpOddClick( dfline, chStr=chStr,show=FALSE, what=what,trda=function(){'/Volumes/rawData/EEG//AmpOddClick/rda/'},extn=fltxtn )  
  ## call the calculations
  gc()
  
  getSessionAverage(tone, maxp2p=maxp2p, maxPowEOG=maxPowEOG, who='all', fextn=fltxtn)
  gc()
  
  print('random environment...')
  random <- subs.data(tone,tone$xpp$trialType==1,comment='nn')
  getSessionAverage(random, maxp2p=maxp2p, maxPowEOG=maxPowEOG, who='random', fextn=fltxtn)
  gc()
  
  #print('regular environment...')
  #regular <- subs.data(tone,tone$xpp$trialType==2,comment='nn')
  #getSessionAverage(regular, maxp2p=maxp2p, maxPowEOG=maxPowEOG, who='regular', fextn=fltxtn)
  #gc()
  
  ## ========== long ISIs to avoid big impact of past and upcoming sounds
  #print('long ISIs...')
  
  isi       <- tone$xpp$ISI
  futureisi <- c(tone$xpp$ISI[-1],12)
  
  longISI <- subs.data(tone,isi>.300 & futureisi>.300,comment='nn')
  getSessionAverage(longISI, maxp2p=maxp2p, maxPowEOG=maxPowEOG, who='longISI', fextn=fltxtn)
  
  rm(tone)
  gc() 
}
