function_aoc_featureExtraction_lfp <- function(ind, linkFileName = 'ampOddclickVProbe_link.txt', epoch = "llp", ovr = FALSE){
  
  #ind<-111
  
  newDefined<-FALSE
  
  baseDir <- '~/Dropbox/'
  library(R.matlab)
  library(rgenoud)
  library(lattice)
  library(akima)
  library(Rwave)
  library(MASS)
  library(fastICA)
  library(mvbutils)
  
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
  
  print(paste("Input to function_aoc_fE_lfp, ind: ", ind))
  
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
  #  .GlobalEnv$pp    <- function(){paste(baseDir, 'ampOddClick/word/VProbe/plots/',dfline$exp,'/ICA_Npca',Nica,'/',sep='')}
  #  pp <- .GlobalEnv$pp
  #  if (!exists(pp()))
  #    system(paste('mkdir -p ', pp(), sep=''))

  parMat<-read.table(paste0('~/Dropbox/RSKernel/r-files/parameters/compPar_', dfline$animal, '.txt'), header = TRUE)
  
  if(file.exists(paste0('~/Dropbox/RSKernel/r-files/parameters/muaPar_', dfline$animal, '.txt'))){
    muaPar<-read.table(paste0('~/Dropbox/RSKernel/r-files/parameters/muaPar_', dfline$animal, '.txt'), header = TRUE)
  }else{
    muaPar<-read.table(paste0('~/Dropbox/RSKernel/r-files/parameters/muaPar_default.txt'), header = TRUE)
  }
  
  tname<-names(parMat)
  parMat<-rbdf(unname(cbind(muaPar, matrix(0, dim(muaPar)[1], dim(parMat)[2]-dim(muaPar)[2]))), unname(parMat))
  names(parMat)<-tname

  if(dfline$animal=="Cirque" | dfline$animal=="Jesse"){
    lfpPar<-read.table(paste0('~/Dropbox/RSKernel/r-files/parameters/lfpPar_', dfline$animal, '.txt'), header = TRUE)
    
    tname<-names(parMat)
    parMat<-rbdf(unname(cbind(lfpPar, matrix(0, dim(lfpPar)[1], dim(parMat)[2]-dim(lfpPar)[2]))), unname(parMat))
    names(parMat)<-tname
  }
  
  fceq<-parMat[,1]=='frontoCentral'
  if(sum(fceq)>0){
    parimax<-min(which(fceq))-1
  }else{
    parimax<-dim(parMat)[1]
  }
  
#  parMat<-read.table(paste0('~/Dropbox/RSKernel/r-files/parameters/compPar_Jesse.txt'), header = TRUE)
#  names(parMat)[5:41]<-1:37
###  parMat[,5:41]<-0
#  write.table(parMat, '~/Dropbox/RSKernel/r-files/parameters/compPar_Cirque.txt', col.names = names(parMat), quote = FALSE, sep = "\t", row.names = FALSE)
  
  if(!is.null(dfline$numChan)){
    nelec<-dfline$numChan
    chStr<-paste0('ch', dfline$startVIndex-1+(1:dfline$numChan))
  }else{
    nelec<-dfline$prelimEndMIndex-dfline$prelimStartMIndex+1
    chStr<-paste0('ch', dfline$prelimStartMIndex-1+(1:nelec))
  }
  
  tryCatch(eeg  <- loadSignalProcessor(dfline, project = 'ampOddClick', epoch = epoch, chStr = chStr),
           error = function(x){print(c('This may not be processed at all or with the current signalProcessor; check the old_rda folder for: ', as.character(dfline$session)))})
  
  if(!is.null(dfline$Layer) & newDefined){
    flipi<-dfline$Layer
    
    #define some basic spatial masks used to find local erp analogs; having a bunch is fine for now
    maskplus<-matrix(0, 13, nelec)
    maskminus<-matrix(0, 13, nelec)
    masknames<-c("Av", "Top", "Bot", "TBdiff", "4up", "4down", "4diff", "8up", "8down", "8diff", "up", "down", "diff")
    ##av
    maskplus[1,]<-rep(1, nelec)
    ##top
    maskplus[2, 1]<-1
    ##bottom
    maskplus[3, nelec]<-1
    ##top-bottom
    maskplus[4, 1]<-1
    maskminus[4, nelec]<-1
    ##1-4 above
    if(round(flipi+0.25-4)>=1){
      maskplus[5, round(flipi+0.25-(1:4))]<-1
    }
    ##1-4 below
    if(round(flipi-0.25+4)<=nelec){
      maskplus[6, round(flipi-0.25+(1:4))]<-1
    }
    ##1-4 diff
    if(round(flipi+0.25-4)>=1 & round(flipi-0.25+4)<=nelec){
      maskplus[7,]<-maskplus[5,]
      maskminus[7,]<-maskplus[6,]
    }
    ##5-8 above
    if(round(flipi+0.25-8)>=1){
      maskplus[8, round(flipi+0.25-(5:8))]<-1
    }
    ##5-8 below
    if(round(flipi-0.25+8)<=nelec){
      maskplus[9, round(flipi-0.25+(5:8))]<-1
    }
    ##5-8 diff
    if(round(flipi+0.25-8)>=1 & round(flipi-0.25+8)<=nelec){
      maskplus[10,]<-maskplus[8,]
      maskminus[10,]<-maskplus[9,]
    }
    ##above
    if(flipi>1){
      maskplus[11, 1:round(flipi-0.75)]<-1
    }
    ##below
    if(flipi<nelec){
      maskplus[12, round(flipi+0.75):nelec]<-1
    }
    ##above-below
    if(flipi>1 & flipi<nelec){
      maskplus[13, 1:round(flipi-0.75)]<-1
      maskminus[13, round(flipi+0.75):nelec]<-1
    }
  
    #oeeg<-eeg
    #odxpp<-dim(eeg$xpp)[2]

    ##erp calc.s one at a time
    ##consider amending for varying weights in the compPar file
    for(i in 1:dim(parMat)[1]){
      if(parMat[i, 1]!='blPup'){
        trng<-as.numeric(parMat[i, 2:3])
        
        for(j in 1:dim(maskplus)[1]){
          if(sum(maskplus[j,])>0){
            signalStr<-paste0(as.character(parMat[i, 1]), "_", masknames[j])
            tdfplus<-eeg$lfp[, which(round(eeg$taxis)==trng[1]):which(round(eeg$taxis)==trng[2]), which(maskplus[j,]!=0), drop = FALSE]
            tdfo<-as.data.frame(apply(tdfplus, 1, mean))
            if(sum(maskminus[j,])>0){
              tdfminus<-eeg$lfp[, which(round(eeg$taxis)==trng[1]):which(round(eeg$taxis)==trng[2]), which(maskminus[j,]!=0), drop = FALSE]
              tdfq<-as.data.frame(apply(tdfminus, 1, mean))
              tdfo<-tdfo-tdfq
            }
            names(tdfo)<-signalStr
            eeg$xpp<-cbind(eeg$xpp, tdfo)
          }
        }
      }
    }
  }
  
  if(ovr){
    parMat<-read.table(paste0('~/Dropbox/RSKernel/r-files/parameters/overridePar_default.txt'), header = TRUE)
    parimax<-dim(parMat)[1]
  }
  
  for(i in 1:parimax){
  #for(i in 1:dim(parMat)[1]){
    if(parMat[i, 1]!='blPup'){
      trng<-as.numeric(parMat[i, 2:3])
      tdfo<-as.data.frame(apply(eeg$lfp[, which(round(eeg$taxis)==trng[1]):which(round(eeg$taxis)==trng[2]), 1], 1, mean))
      names(tdfo)<-paste0(as.character(parMat[i, 1]), "_", 1)
      for(j in 2:dim(eeg$lfp)[3]){
        signalStr<-paste0(as.character(parMat[i, 1]), "_", j)
        tdfplus<-eeg$lfp[, which(round(eeg$taxis)==trng[1]):which(round(eeg$taxis)==trng[2]), j]
        tdfo<-cbind(tdfo, as.data.frame(apply(tdfplus, 1, mean)))
        names(tdfo)[j]<-signalStr
      }
      eeg$xpp<-cbind(eeg$xpp, tdfo)
    }
  }
  
  xpp<-eeg$xpp


  tdir<-paste0('/Volumes/EEGLab/EEGLab/ampOddClick/', dfline$session)
  #  tdir<-paste0('/Volumes/Drobo5D3/EEGLab/ampOddClick/', dfline$session)
  
  #  print(c("dir exists? ", dir.exists(tdir)))
  if (!dir.exists(tdir)){
    system(paste('mkdir -p ', tdir, sep=''))
  }
  
  save(xpp, file = paste0(tdir, '/features_lfp_', epoch, '.rda'))
  print(c("LFP success: ", as.character(dfline$session)))

  return(epoch)
}

