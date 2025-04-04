function_aoc_featureExtraction_csd <- function(ind, linkFileName = 'ampOddclickVProbe_link.txt', epoch = "raw", ovr = FALSE){
  
  #ind<-128
  ##FLA_VProbe_STPSP
  
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
  source( paste(baseDir,'RSkernel/r-files/functions/calcCSD.R',sep='') )
  
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
  
  print(paste("Input to function_aoc_fE_csd, ind: ", ind))
  
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
  
  parMat<-read.table(paste0('~/Dropbox/RSKernel/r-files/parameters/compPar_', dfline$animal, '.txt'), header = TRUE)
  
  ## preliminaries
  rda  <- function(){'/Volumes/EEGLab/EEGLab/ampOddClick'}
  trda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
  .GlobalEnv$rda <- rda
  .GlobalEnv$trda <- trda
  
  chStr<-paste0('ch', dfline$startVIndex-1+(1:dfline$numChan))
  tryCatch(lfp  <- loadSignalProcessor(dfline, project = 'ampOddClick', epoch = epoch, chStr = chStr),
           error = function(x){print(c('This may not be processed at all or with the current signalProcessor; check the old_rda folder for: ', as.character(df$session[ind])))})
  
  lag     <- 1
  smooth  <- T
  sd      <- 100
  spacing <- dfline$spacing
  varCor  <- F
  
  lfpnbl <- lfp
  tmp            <- baseLineCorrection(lfp$lfp, range(lfp$taxis), lfp$taxis)
  lfpnbl$lfp <- tmp$res
  csd      <- calcSingleTrialCSD(lfpnbl, spacing = spacing, lag = lag, smooth = T, sd = sd)
  csd$what <- 'csd'
  
  avrec          <- apply(abs(csd$lfp), c(1, 2), sum)
  tmp            <- baseLineCorrection(avrec, range(csd$taxis), csd$taxis)
  avrec          <- tmp$res
  
  
  csd$lfp[,,dim(csd$lfp)[3]] <- avrec
  
  # ----------------
  tmp          <- apply(csd$lfp,c(1,3),max) - apply(csd$lfp,c(1,3),min)
  csd$xpp$p2p  <- apply( tmp, 1, max )
  
  qtmp         <- quantile(csd$xpp$p2p, probs=c(0.5, 0.75))
  mxcsd        <- qtmp[1]+10*diff(qtmp)
  valIndCSD    <- csd$xpp$p2p < mxcsd;
  csd$xpp$valTrial <- valIndCSD
  

  ##In the case that we want to override all components with the sparse common defaults
  if(ovr){
    parMat<-read.table(paste0('~/Dropbox/RSKernel/r-files/parameters/overridePar_default.txt'), header = TRUE)
  }else{
      
    odxpp<-dim(csd$xpp)[2]
    
    for(i in 1:dim(parMat)[1]){
      if(parMat[i, 1]!='blPup'){
        trng<-as.numeric(parMat[i, 2:3])
        signalStr<-as.character(parMat[i, 1])
  
        csd$xpp[,(odxpp+i)]<-apply(avrec[,which(round(csd$taxis)==trng[1]):which(round(csd$taxis)==trng[2])], 1, mean)
        names(csd$xpp)[odxpp+i]<-signalStr
      }
    }
  
    ##this is pulled from the csd scripts - preivous parMats are ignored
    parMat<-matrix(c('Mt', 8, 12, 'Me', 20, 40, 'Mi', 60, 100, 'CSD0', 15, 19, 'CSD1', 25, 45, 'CSD2', 50, 60, 'CSD3', 65, 110, 
                     'CSD4', 100, 150, 'CSD5', 160, 250, 'MUAbl', -50, 0, 'MUAR1', 18, 35, 'MUAR1a', 15, 25, 'MUAR1b', 25, 35, 
                     'MUAR1c', 35, 45, 'MUAR2', 55, 125, 'MUAR3', 160, 250), ncol = 3, byrow = TRUE)
    
  }

  for(i in 1:dim(parMat)[1]){
    twin<-which(round(csd$taxis)==parMat[i, 2]):which(round(csd$taxis)==parMat[i, 3])
    csd$xpp<-cbind(csd$xpp, apply(avrec[,twin], 1, mean))
    names(csd$xpp)[dim(csd$xpp)[2]]<-paste0(parMat[i, 1])
  }
  
  xpp        <- csd$xpp
  
  ## ================================ save the feature map
  tdir<-paste0('/Volumes/EEGLab/EEGLab/ampOddClick/', dfline$session)

  if (!dir.exists(tdir)){
    system(paste('mkdir -p ', tdir, sep=''))
  }
  
  save(xpp, file = paste0(tdir, '/features_csd_', epoch, '.rda'))
  print(c("CSD success: ", as.character(df$session[ind])))
  
  return(epoch)
}

