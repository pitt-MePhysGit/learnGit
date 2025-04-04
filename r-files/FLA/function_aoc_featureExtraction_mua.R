function_aoc_featureExtraction_mua <- function(ind, linkFileName = 'ampOddclickVProbe_link.txt', ovr = FALSE){
  
  #ind<-116
  #ind<-111
  
  calcPC<-TRUE
  nPC<-3
  
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
#  source( paste(baseDir,'aMMN_controlled_temp/r-files/functions_eCog.R',sep='') )
  
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

  print(paste("Input to function_aoc_fE_mua, ind: ", ind))
  
  # read the curated file
  df          <- prepDFsignalProcessor("ampOddclick", linkFile=linkFileName)
  if(tolower(linkFileName)=='ampoddclickvprobe_link.txt'){
    df         <- prepDFVprobe(df)
  }
  dfline<-df[ind,]
  
  ## preliminaries
  rda  <- function(){'/Volumes/EEGLab/EEGLab/ampOddClick'}
  trda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
  .GlobalEnv$rda <- rda
  .GlobalEnv$trda <- trda

  ## ========================= load MUA data
  if(dfline$animal=="Jesse" | dfline$animal=="Cirque"){
    chInd      <- dfline$prelimStartMIndex:dfline$prelimEndMIndex
  }else{
    chInd      <- dfline$startVIndex - 1 + 1:dfline$numChan
  }
  chStr      <- paste('ch',chInd,sep='')
  
#  suppressWarnings(mua        <- loadAOC( df[ind,] , what='mua', chStr=chStr, trda=trda,blcorr=T  ))
  mua  <- loadSignalProcessor(dfline, project = 'ampOddClick', epoch = "mua", chStr = chStr)
  
  ## ======================================================
  tmp          <- apply(mua$lfp, c(1,3), max) - apply(mua$lfp, c(1, 3), min)
  mua$xpp$p2p  <- apply(tmp, 1, max)
  muaquant     <- quantile(mua$xpp$p2p, probs=c(0.25, 0.75))
  mxp2pmua     <- muaquant[2] + 1.5*diff(muaquant)
#  mua$xpp$valTrial <- mua$xpp$p2p < mxp2pmua
  mua$xpp$valTrial <- mua$xpp$p2p < mxp2pmua & 1:dim(mua$xpp)[1]>5
  
  brks <- 0.250 * 2^(0:7)
  brks[1] <- 0.2#; brks[length(brks)]<- 22 

  mua$xpp$isiOct  <- cut(mua$xpp$ISI, breaks = brks)
  mua$xpp$futureISI <- my.lag(mua$xpp$ISI, lag = -1)
  

  ## =================== extract average MUA per trial and time-window
  extractScalars <- function(td = mua, trng = c(8, 12), ext = 't'){
    taxis <- td$taxis
    tndx  <- which(taxis>=trng[1] & taxis<=trng[2])
    res   <- apply(td$lfp[,tndx,], c(1, 3), mean)
    dimnames(res)[[2]] <- paste('M_', ext, '_', 1:dim(td$lfp)[3], sep = '')
    return(res)
  }
  
#  muaPar<-matrix(c('t', 8, 12, 'e', 20, 40, 'i', 60, 100), ncol = 3, byrow = TRUE)
  
  if(file.exists(paste0('~/Dropbox/RSKernel/r-files/parameters/muaPar_', dfline$animal, '.txt'))){
    muaPar<-read.table(paste0('~/Dropbox/RSKernel/r-files/parameters/muaPar_', dfline$animal, '.txt'), header = TRUE)
  }else{
    muaPar<-read.table(paste0('~/Dropbox/RSKernel/r-files/parameters/muaPar_default.txt'), header = TRUE)
  }
  
#  muaPar<-cbind(as.data.frame(muaPar), "mua")
#  names(muaPar)<-c("name", "from", "to", "sgnl")
#  write.table(muaPar, '~/Dropbox/RSKernel/r-files/parameters/muaPar_default.txt', col.names = names(muaPar), quote = FALSE, sep = "\t", row.names = FALSE)

  if(ovr){
    muaPar<-read.table(paste0('~/Dropbox/RSKernel/r-files/parameters/overridePar_default.txt'), header = TRUE)
  }
  
  ##channel-by-channel components
  predMatMua <- c()
  for(i in 1:dim(muaPar)[1]){
    twin<-which(round(mua$taxis)==muaPar[i, 2]):which(round(mua$taxis)==muaPar[i, 3])
    predMatMua <- cbind(predMatMua, extractScalars(td = mua, trng=as.numeric(c(muaPar[i, 2], muaPar[i, 3])), ext = muaPar[i, 1]))
#    assign(paste0('pc', muaPar[i, 1]), prcomp(apply(mua$lfp[,twin,], c(1, 3), mean), scale = TRUE)$rot)
  }
  mua$xpp<-cbind(mua$xpp, predMatMua)
  
  ##PC calcs and resulting spatial PCs over the probe/array
  if(calcPC){
    apcs<-list()
    for(i in 1:dim(muaPar)[1]){
      twin<-which(round(mua$taxis)==muaPar[i, 2]):which(round(mua$taxis)==muaPar[i, 3])
      apcs[[i]]<-prcomp(apply(mua$lfp[,twin,], c(1, 3), mean), scale = TRUE)$rot
      for(j in 1:nPC){
        mua$xpp<-cbind(mua$xpp, apply(mua$lfp[,twin,], c(1, 3), mean)%*%t(t(apcs[[i]][,j])))
        names(mua$xpp)[dim(mua$xpp)[2]]<-paste0("M", muaPar[i, 1], "PC", j)
      }
    }
  }
  
  
  
  
  #old junk
if(1==0){    
  #this needs a loop
  #  apcs<-list(pct, pce, pci)
  #  apcs<-lapply(paste0('pc', muaPar[,1]), get)
  
  
  predMatMua <-                   extractScalars(td=mua, trng=c( 8, 12), ext='t')
  predMatMua <- cbind(predMatMua, extractScalars(td=mua, trng=c(20, 40), ext='e'))
  predMatMua <- cbind(predMatMua, extractScalars(td=mua, trng=c(60,100), ext='i'))

  pct<-prcomp(apply(mua$lfp[,8:12,], c(1, 3), mean), scale = TRUE)$rot
  pce<-prcomp(apply(mua$lfp[,20:40,], c(1, 3), mean), scale = TRUE)$rot
  pci<-prcomp(apply(mua$lfp[,60:100,], c(1, 3), mean), scale = TRUE)$rot
  apcs<-list(pct, pce, pci)

  for(i in 1:3){
    mua$xpp<-cbind(mua$xpp, apply(mua$lfp[,8:12,], c(1, 3), mean)%*%t(t(pct[,i])))
    names(mua$xpp)[dim(mua$xpp)[2]]<-paste0("MtPC", i)

    mua$xpp<-cbind(mua$xpp, apply(mua$lfp[,20:40,], c(1, 3), mean)%*%t(t(pce[,i])))
    names(mua$xpp)[dim(mua$xpp)[2]]<-paste0("MePC", i)
    
    mua$xpp<-cbind(mua$xpp, apply(mua$lfp[,60:100,], c(1, 3), mean)%*%t(t(pci[,i])))
    names(mua$xpp)[dim(mua$xpp)[2]]<-paste0("MiPC", i)
    
#    assign(paste0("mua$xpp$MtPC", i), apply(mua$lfp[,8:12,], c(1, 3), mean)%*%t(t(pct[,i])), inherits = TRUE)
#    assign(paste0("mua$xpp$MePC", i), apply(mua$lfp[,20:40,], c(1, 3), mean)%*%t(t(pce[,i])))
#    assign(paste0("mua$xpp$MiPC", i), apply(mua$lfp[,60:100,], c(1, 3), mean)%*%t(t(pci[,i])))
  }
}
  
  
  
  #library(glmnet)
  #lambdas <- 10^seq(2, -3, by = -.1)
  #stry<-cv.glmnet(apply(mua$lfp[,8:12,1:96], c(1, 3), mean), mua$xpp$ISI, alpha = 1, lambda = lambdas, standardize = TRUE)
  #stry2<-glmnet(apply(mua$lfp[,8:12,1:96], c(1, 3), mean), mua$xpp$ISI, alpha = 1, lambda = stry$lambda.min, standardize = TRUE)
  #plot(stry2)

  #stry<-cv.glmnet(apply(mua$lfp[,8:12,1:96], c(1, 3), mean), mua$xpp$ampInd, alpha = 1, lambda = lambdas, standardize = TRUE)
  #stry2<-glmnet(apply(mua$lfp[,8:12,1:96], c(1, 3), mean), mua$xpp$ampInd, alpha = 1, lambda = stry$lambda.min, standardize = TRUE)
  #plot(stry2)
  
#  browser()
  
  xpp        <- mua$xpp
#  xpp$p2pmua <- mua$xpp$p2p
#  xpp$valTrial <- mua$xpp$p2p<mxp2pmua & 1:dim(xpp)[1]>5
  
  
  ## ================================ save the feature map
  if(is.null(dfline$exp)){
    saveDir  <- paste0(rda(), '/', dfline$session)
  }else{
    saveDir  <- paste0(rda(), '/', dfline$exp)
  }
  
  if(!exists( saveDir )){
    system(paste('mkdir -p ', saveDir, sep=''))
  }
  
  save(xpp, apcs, file = paste0(saveDir, '/features_mua.rda'))
  print(paste0("MUA success: ", as.character(dfline$session)))
  
}

