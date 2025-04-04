batchFLA_ampOddClick_Spencer <- function(projectName, linkFileName, ndx, modality = 'LMECF', ovr = TRUE){
  
  ##in future, will want ability to pass the modalities FLA is wanted for
  ##E: EEG, N: electrodes without CSD, C: electrodes with CSD, M: MUA
  ##ENM would be allowable, but the index should be appropriate
  #
  #projectName<-'ampOddClick'
  #linkFileName<-'ampOddclickVProbe_link.txt'
  #ndx<-61
  #modality<-'C'
  
  #linkFileName<-'ampOddclick_link.txt'
  #ndx<-111
  #ndx<-116
  #ndx<-57
  #ndx<-109
  
  #batchFLA_ampOddClick('ampOddClick', 'ampOddclick_link.txt', 111)
  
  
  source('~/Dropbox/r-compile/signalProcessoR/R/signalProcessoR.R')
  sourceDir( paste( '~/Dropbox/',projectName,'/r-files/FLA/',sep='')  )
  
  df          <- prepDFsignalProcessor(projectName, linkFileName)
  dfline      <- df[ndx,]
  
  if(!grepl('v', dfline$ElectrodeConfig) & !grepl('d', dfline$ElectrodeConfig) & !grepl('t', dfline$ElectrodeConfig) & grepl('C', modality)){
    modality<-gsub('C', '', modality)
  }
  
  ##need to change function_aoc_featureExtraction to accept linkFileName
  run_aocfeatures  <-TRUE
  run_LRMfit       <-TRUE
  
  
  
  
  
  lep<-''
  eep<-''
  cep<-''
  
  ##need to make sure featureExtraction can handle the sc96 (maybe different function?)
  if(grepl('N', modality)){
    if(run_aocfeatures){
#      function_aoc_featureExtraction(ndx, Nica, linkFileName = linkFileName)
    }
  }

  if(grepl('L', modality)){
    if(run_aocfeatures){
      if(any(unlist(lapply(c('v', 't', 'd', 's', 'c', 'u'), function(x){grepl(x, dfline$ElectrodeConfig)})))){
        lep<-function_aoc_featureExtraction_lfp(ndx, linkFileName = linkFileName, epoch = "raw", ovr)
      }
    }
  }

  if(grepl('I', modality)){
    if(run_aocfeatures){
#      function_aoc_featureExtraction_ica(ndx, linkFileName = linkFileName)
    }
  }

  if(grepl('E', modality)){
    if(run_aocfeatures){
      ##I've set the default epoch to raw because Sam only has llp preprocessing for a few sessions in 2016
      eep<-function_aoc_featureExtraction_eeg(ndx, linkFileName = linkFileName, epoch = "raw", ovr)
    }
  }

  if(grepl('C', modality)){
    if(run_aocfeatures){
      cep<-function_aoc_featureExtraction_csd(ndx, linkFileName = linkFileName, epoch = "raw", ovr)
    }
  }

  if(grepl('M', modality)){
    if(run_aocfeatures){
##putting in the exceptions for which link files we should run the feature extraction on; dual excluded and eeg 
##probably need more
      if(!grepl('dual', linkFileName) & !grepl('EEG', linkFileName)){
        function_aoc_featureExtraction_mua(ndx, linkFileName = linkFileName, ovr)
      }
    }
  }
  
  epl<-c(lep, eep, cep)

  if(grepl('F', modality)){
    if(modality == 'F'){
      modality<-'LMECF'
    }
    fitmod<-gsub('F', '', modality)
    if(run_LRMfit){
      function_aoc_LRMfit(ndx, linkFileName = linkFileName, fitmod, epl)
    }
  }
  
}