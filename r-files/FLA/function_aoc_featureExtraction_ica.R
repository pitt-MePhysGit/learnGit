function_aoc_featureExtraction_ica <- function(ind, linkFileName = 'ampOddclickVProbe_link.txt'){
  
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
  source( paste(baseDir,'aMMN_controlled_temp/r-files/functions_eCog.R',sep='') )
  
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
  
  print(paste("Input to function_aoc_fE_ica, ind: ", ind))
  
  # read the curated file
  df          <- prepDFsignalProcessor("ampOddclick", linkFile=linkFileName)
  df         <- prepDFVprobe(df)
  
  ## preliminaries
  rda  <- function(){'/Volumes/EEGLab/EEGLab/ampOddClick'}
  trda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
  .GlobalEnv$rda <- rda
  .GlobalEnv$trda <- trda
  
  chInd      <- 1:6
  chStr      <- paste('ic',chInd,sep='')
  suppressWarnings(llp        <- loadAOC( df[ind,], what='ICA_evk_N6', chStr=chStr, trda=trda, blcorr=T, loadCC=T  ))
  
  ## ======================================================
  tmp          <- apply(llp$lfp,c(1,3),max) - apply(llp$lfp,c(1,3),min)
  llp$xpp$p2p  <- apply( tmp, 1, max )
  llpquant         <- quantile(llp$xpp$p2p, probs=c(.25,.75) )
  ##about 3 s.d. above mean in theory, 2 with practical non-normality
  mxp2plfp         <- llpquant[2] + 1.5*diff(llpquant)
  mx2p2lfp         <- (mxp2plfp%>=%15)%<=%100
  llp$xpp$valTrial <- llp$xpp$p2p < mxp2plfp
  
  brks <- 0.250 * 2^(0:7)
  brks[1] <- 0.2#; brks[length(brks)]<- 22 
  llp$xpp$isiOct  <- cut(llp$xpp$ISI,  breaks=brks)
  llp$xpp$futureISI <- my.lag(llp$xpp$ISI,lag=-1)
  
  llp$srt   <- sort(apply(llp$weight,1,var),decreasing=T, index.return=T)$ix
  llp$var <- apply(llp$weights,1,var)
  
  ## =================== extract average MUA per trial and time-window
  extractScalars <- function(td = mua, trng = c(8,12), ext = 't'){
    taxis <- td$taxis
    tndx  <- which(taxis>=trng[1] & taxis<=trng[2])
    res   <- apply( td$lfp[,tndx,], c(1,3), mean)
    dimnames(res)[[2]] <- paste('M_',ext,'_',1:dim(td$lfp)[3],sep='')
    return(res)
  }
  
  predMatIC <-                  extractScalars(td=llp, trng=c( 8, 12), ext='IC_t')
  predMatIC <- cbind(predMatIC, extractScalars(td=llp, trng=c(20, 40), ext='IC_e'))
  predMatIC <- cbind(predMatIC, extractScalars(td=llp, trng=c(60,100), ext='IC_i'))
  
  xpp        <- llp$xpp
  xpp$p2plfp <- llp$xpp$p2p
  xpp$valTrial <- xpp$p2plfp < mxp2plfp & 1:dim(xpp)[1]>5
  
  xpp <- cbind(xpp, predMatIC)
  
  ## ================================ save the feature map
  saveDir  <- paste0(rda(), '/', df$exp[ind])
  
  if(!exists( saveDir )){
    system(paste('mkdir -p ', saveDir, sep=''))
  }
  
  save(xpp, file = paste0(saveDir, '/features_ica.rda'))
  print(paste0("ICA success: ", as.character(df$session[ind])))
  
}

