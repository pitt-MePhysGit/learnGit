function_aoc_featureExtraction_eeg <- function(ind, linkFileName = 'ampOddclickVProbe_link.txt', epoch = 'raw', ovr = FALSE){
  
  ##The blPup component is specifically excluded; can be turned by flipping:
  blp<-FALSE
  ##code will not error if lacking eye tracking though, the component just won't be calculated
  
  baseDir <- '~/Dropbox/'
  library(R.matlab)
  library(rgenoud)
  library(lattice)
  library(akima)
  library(Rwave)
  library(MASS)
  library(fastICA)

  source( paste(baseDir,'ampOddClick/r-files/functions/prepDFVprobe.R',sep=''))
  source( paste(baseDir,'ampOddClick/r-files/functions/FLA_VProbe_STPSP.R',sep='') )
  
  source( paste(baseDir,'ampOddClick/r-files/functions/functions_ampOddclick.R',sep='') )
  source( paste(baseDir,'ampOddClick/r-files/functions/functions_AOCdual.R',sep='') )
  source( paste(baseDir,'ampOddClick/r-files/functions/functions_VProbe_STPSP.R',sep='') )
  
  source( paste(baseDir,'r-compile/rUtils/R/rUtils.R',sep='') )
  source( paste(baseDir,'r-compile/ePhysLab/R/ePhysLab.R',sep='') )
  
  source( paste(baseDir,'RSkernel/r-files/functions/tapered_ccf.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/lamiarLFPtheory/smooth2DCSD.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/call_fastICA.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/save_fastICA.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/lamiarLFPtheory/showICAcomponents.R',sep='') )
  
  
  source( paste(baseDir,'RSkernel/r-files/functions/eCogFileStr.R',sep='') )
#  source( paste(baseDir,'aMMN_controlled/r-files/functions_eCog.R',sep='') )
  
  
  source( paste(baseDir,'ampOddClick/r-files/functions/getCWTevk.R',sep='') )
  source( paste(baseDir,'ampOddClick/r-files/functions/getCWT.R',sep='') )
  source( paste(baseDir,'ampOddClick/r-files/functions/showCWT.R',sep='') )
  
  source( paste(baseDir,'ampOddClick/r-files/functions/showICx.R',sep='') )
  source( paste(baseDir,'ampOddClick/r-files/functions/getSLix.R',sep='') )
  source( paste(baseDir,'ampOddClick/r-files/functions/analyzeCorrelationStructure.R',sep=''))
  source( paste(baseDir,'ampOddClick/r-files/functions/decodeSVM.R',sep=''))
  
  source( paste(baseDir,'ampOddClick/r-files/figures/showTD.R',sep='') )
  source( paste(baseDir,'ampOddClick/r-files/figures/showTD2.R',sep='') )
  source( paste(baseDir,'ampOddClick/r-files/figures/showBandsByVar.R',sep='') )
  source( paste(baseDir,'ampOddClick/r-files/figures/showMUAByVar.R',sep='') )
  
  source('~/Dropbox/r-compile/signalProcessoR/R/signalProcessoR.R')
  
  ##linkFileName = 'ampOddclickEEG_link.txt'
  # read the curated file
  colClasses <- c('character', rep('numeric', 9), rep('character', 3), 'numeric')
  if(tolower(linkFileName) == 'ampoddclick_link.txt'){
    colClasses <- c('character', rep('numeric', 8), rep('character', 3))
  }
  if(tolower(linkFileName) == 'ampoddclickvprobe_link.txt'){
    colClasses <- c('character', rep('numeric', 2), 'character', rep('numeric', 5), rep('character', 3), rep('numeric', 9))
  }
  #  df         <- read.table(paste(baseDir,'ampOdd_click/ampOddclickVProbe_link.txt',sep=''),header=T,colClasses=colClasses)

  print(paste0(baseDir,'ampOddClick/', linkFileName))
  df         <- read.table(paste0(baseDir,'ampOddClick/', linkFileName), header = T, colClasses = colClasses)
  df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))

  vprobei<-grepl('vprobe', tolower(linkFileName))
  if(vprobei){
    df         <- prepDFVprobe(df)
  }else{
    df         <- prepDFsignalProcessor('ampOddClick', linkFileName)
  }
  dfline<-df[ind,]
  
  ## preliminaries
  #pp   <- function(){paste(baseDir, 'ampOdd_click/word/VProbe/plots/',sep='')}
  #  rda  <- function(){'~/Dropbox/ampOdd_click/Rdata/'}
  rda  <- function(){'/Volumes/EEGLab/EEGLab/ampOddClick'}
##this may need to be adjusted based on what exists in old_rda vs. rda
  trda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
  .GlobalEnv$rda <- rda
  .GlobalEnv$trda <- trda
  
  ##maybe some plotting?
#  .GlobalEnv$pp    <- function(){paste(baseDir, 'ampOddClick/word/VProbe/plots/',df$exp[ind],'/ICA_Npca',Nica,'/',sep='')}
#  pp <- .GlobalEnv$pp
#  if (!exists(pp()))
#    system(paste('mkdir -p ', pp(), sep=''))
  
  
  parMat<-read.table(paste0('~/Dropbox/RSKernel/r-files/parameters/compPar_', dfline$animal, '.txt'), header = TRUE)

  chMat<-paste0('c', 1:(dim(parMat)[2]-4))
  names(parMat)[5:dim(parMat)[2]]<-chMat
  
  chStr<-paste0("ch", 1:ifelse(tolower(dfline$animal)=="rockey", 21, 33))
  tryCatch(eeg  <- loadSignalProcessor(dfline, project = 'ampOddClick', epoch = epoch, chStr = chStr),
           error = function(x){print('This may not be processed at all or with the current signalProcessor; check the old_rda folder')})
  
  odxpp<-dim(eeg$xpp)[2]

  ##erp calc.s one at a time
  ##consider amending for varying weights in the compPar file
  
  ##now we're getting rid of everything after and including frontoCentral
  fceq<-parMat[,1]=='frontoCentral'
  if(sum(fceq)>0){
    parimax<-min(which(fceq))-1
  }else{
    parimax<-dim(parMat)[1]
  }
  
  if(ovr){
    parMatOv<-read.table(paste0('~/Dropbox/RSKernel/r-files/parameters/overridePar_default.txt'), header = TRUE)
    nrOv<-dim(parMatOv)[1]
    ncOv<-dim(parMatOv)[2]
    parMatOv<-cbind(parMatOv, matrix(1, nrOv, 33), matrix(0, nrOv, dim(parMat)[2]-ncOv-33))
    parMatOv<-rbind(parMatOv, parMatOv)
    parMatOv[,1]<-as.character(parMatOv[,1])
    parMatOv[1:nrOv, 1]<-paste0(parMatOv[1:nrOv, 1], "Av")
    
    masklink<-rep(0, nrOv)
    for(i in 1:nrOv){
      masklink[i]<-min(which.min((parMatOv[i,2]-parMat[,2])^2+(parMatOv[i,3]-parMat[,3])^2))
    }
    parMatOv[(nrOv+1):(2*nrOv), (ncOv+1):dim(parMat)[2]]<-parMat[masklink, (ncOv+1):dim(parMat)[2]]

    parMat<-parMatOv
    parimax<-dim(parMat)[1]
    for(i in 1:parimax){
      chInd<-which(parMat[i, 5:dim(parMat)[2]]!=0)
      trng<-as.numeric(parMat[i, 2:3])
      signalStr<-as.character(parMat[i, 1])
      
      eeg$xpp[,(odxpp+i)]<-apply(eeg$lfp[, which(round(eeg$taxis)==trng[1]):which(round(eeg$taxis)==trng[2]), chInd], 1, mean)
      names(eeg$xpp)[odxpp+i]<-signalStr
      
      #We're now calculating average epochs across all eeg electrodes for each override epoch
      #and within the mask that is most temporally similar from that animal's original compPar
#
#
#      trng<-as.numeric(parMat[i, 2:3])
#      tdfo<-as.data.frame(apply(eeg$lfp[, which(round(eeg$taxis)==trng[1]):which(round(eeg$taxis)==trng[2]), 1], 1, mean))
#      names(tdfo)<-paste0(as.character(parMat[i, 1]), "_", 1)
#      for(j in 2:dim(eeg$lfp)[3]){
#        signalStr<-paste0(as.character(parMat[i, 1]), "_", j)
#        tdfplus<-eeg$lfp[, which(round(eeg$taxis)==trng[1]):which(round(eeg$taxis)==trng[2]), j]
#        tdfo<-cbind(tdfo, as.data.frame(apply(tdfplus, 1, mean)))
#        names(tdfo)[j]<-signalStr
#      }
#      eeg$xpp<-cbind(eeg$xpp, tdfo)
      
    }

  }else{
    
    for(i in 1:parimax){
    #for(i in 1:dim(parMat)[1]){
      if(parMat[i, 1]!='blPup'){
        chInd<-which(parMat[i, 5:dim(parMat)[2]]!=0)
        trng<-as.numeric(parMat[i, 2:3])
        signalStr<-as.character(parMat[i, 1])
  
        eeg$xpp[,(odxpp+i)]<-apply(eeg$lfp[, which(round(eeg$taxis)==trng[1]):which(round(eeg$taxis)==trng[2]), chInd], 1, mean)
        names(eeg$xpp)[odxpp+i]<-signalStr
      }
    }
    
    if(blp & sum(parMat[,1]=='blPup')==1){
      blpi<-which(parMat[,1]=='blPup')
      chInd<-which(parMat[blpi, (5:dim(parMat)[2])]!=0)
      tryCatch(eeg  <- loadSignalProcessor(dfline, project = 'ampOddClick', epoch = 'raw', chStr = paste0('ch', chInd)),
               error = function(x){print('This day/animal may not have had eye tracking.')})
  
      trng<-as.numeric(parMat[blpi, 2:3])
      eeg$xpp[,(odxpp+i)]<-apply(eeg$lfp[, which(round(eeg$taxis)==trng[1]):which(round(eeg$taxis)==trng[2]), chInd], 1, mean)
      names(eeg$xpp)[odxpp+i]<-'blPup'
    }
    
  }
  
  
  
  xpp<-eeg$xpp
  
  tdir<-paste0('/Volumes/EEGLab/EEGLab/ampOddClick/', dfline$session)
#  tdir<-paste0('/Volumes/Drobo5D3/EEGLab/ampOddClick/', df$session[ind])
  
#  print(c("dir exists? ", dir.exists(tdir)))
  if (!dir.exists(tdir)){
      system(paste('mkdir -p ', tdir, sep=''))
  }
  
  save(xpp, file = paste0(tdir, '/features_eeg_', epoch, '.rda'))
  print(c("EEG Success: ", as.character(dfline$session)))
  
  return(epoch)
}
  
