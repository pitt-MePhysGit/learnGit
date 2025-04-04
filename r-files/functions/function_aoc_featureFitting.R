function_aoc_featureFitting <- function(ind,Nica=6){
  
  baseDir <- '~/Dropbox/'
  library(R.matlab)
  library(rgenoud)
  library(lattice)
  library(akima)
  library(Rwave)
  library(MASS)
  library(fastICA)
  library(e1071)
  #library(circlize)
  
  
  source(paste(baseDir, 'ampOdd_click/r-files/functions/prepDFVprobe.R',sep=''))
  source( paste(baseDir,'ampOdd_click/r-files/functions/FLA_VProbe_STPSP.R',sep='') )
  
  source( paste(baseDir,'ampOdd_click/r-files/functions/functions_ampOddclick.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/functions/functions_AOCdual.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/functions/functions_VProbe_STPSP.R',sep='') )
  
  source( paste(baseDir,'r-compile/rUtils/R/rUtils.R',sep='') )
  source( paste(baseDir,'r-compile/ePhysLab/R/ePhysLab.R',sep='') )
  
  source( paste(baseDir,'RSkernel/r-files/functions/tapered_ccf.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/lamiarLFPtheory/smooth2DCSD.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/call_fastICA.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/save_fastICA.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/lamiarLFPtheory/showICAcomponents.R',sep='') )
  
  source( paste(baseDir,'RSkernel/r-files/functions/eCogFileStr.R',sep='') )
  source( paste(baseDir,'aMMN_controlled/r-files/functions_eCog.R',sep='') )
  
  source( paste(baseDir,'ampOdd_click/r-files/functions/getCWTevk.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/functions/getCWT.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/functions/showCWT.R',sep='') )
  
  source( paste(baseDir,'ampOdd_click/r-files/functions/showICx.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/functions/getSLix.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/functions/analyzeCorrelationStructure.R',sep=''))
  source( paste(baseDir,'ampOdd_click/r-files/functions/decodeSVM.R',sep=''))
  
  source( paste(baseDir,'ampOdd_click/r-files/figures/showTD.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/figures/showTD2.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/figures/showBandsByVar.R',sep='') )
  source( paste(baseDir,'ampOdd_click/r-files/figures/showMUAByVar.R',sep='') )
  
  # read the curated file
  colClasses <- c('character',rep('numeric',2),'character',rep('numeric',5)  ,rep('character',3),rep('numeric',9) )
  df         <- read.table(paste(baseDir,'ampOdd_click/ampOddclickVProbe_link.txt',sep=''),header=T,colClasses=colClasses)
  df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))
  df         <- prepDFVprobe(df)
  
  ## preliminaries
  #pp   <- function(){paste(baseDir, 'ampOdd_click/word/VProbe/plots/',sep='')}
  rda  <- function(){'~/Dropbox/ampOdd_click/Rdata/'}
  trda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
  .GlobalEnv$rda <- rda
  .GlobalEnv$trda <- trda
  
  # to do:

    ##====================================================== load example ICA component
  .GlobalEnv$pp    <- function(){paste(baseDir, 'ampOdd_click/word/VProbe/plots/',df$exp[ind],'/ICA_Npca',Nica,'/',sep='')}
  pp <- .GlobalEnv$pp
  if (!exists(pp()))
    system(paste('mkdir -p ', pp(), sep=''))
  
  if (!exists(paste( pp(), '/STPSP/', sep='')))
    system(paste('mkdir -p ', pp(),'/STPSP/', sep=''))
  
  
  load( file=paste(rda(),'/',df$exp[ind],'/features.rda',sep='') )
  
  firstLetter <- function(nx,df=xpp){ substr(names(df)[nx],1,1) }
  fLlst <- sapply(1:dim(xpp)[2], firstLetter)
  
  tx     <- xpp[, fLlst %in% c('H','G','B','A','T','D','M')  ]
  tn     <- list(xpp=xpp)
  tn$exp <-  df$exp[ind]
  ## ====================================== run model fits for each predictor
  source( paste(baseDir,'ampOdd_click/r-files/functions/functions_ampOdd_presynapticPlasticity.R',sep='') )
  
  
  # browser()
  # tst    <- callSTPSP(tn,   DVstr='B_e_3', lmOffset=F, gd=T, show=T, showCh=1, who='all', BYstr='exp', FUN=ampFUN,    minISI=.190, maxISI=21, startFromScratch=T, eps=T)
  # tstscl <- callSTPSP(tn,   DVstr='B_e_3', lmOffset=F, gd=T, show=T, showCh=1, who='all', BYstr='exp', FUN=ampFUNscl, minISI=.190, maxISI=21, startFromScratch=T, eps=T)
  # 
  # tstO    <- callSTPSP(tn,   DVstr='B_e_3', lmOffset=T, gd=T, show=T, showCh=1, who='all', BYstr='exp', FUN=ampFUN,    minISI=.190, maxISI=21, startFromScratch=T, eps=T)
  # tstsclO <- callSTPSP(tn,   DVstr='B_e_3', lmOffset=T, gd=T, show=T, showCh=1, who='all', BYstr='exp', FUN=ampFUNscl, minISI=.190, maxISI=21, startFromScratch=T, eps=T)
  # 
  # 
  # anova(tst$res$lm1,tst$res$lm2)
  # 
  # tst$res$modelTest$`Pr(>F)`[2]
  # tstscl$res$modelTest$`Pr(>F)`[2]
  # 
  # tst$res$modelTest
  # tstscl$res$modelTest
  # 
  # tstO$res$modelTest
  # tstsclO$res$modelTest
  # 
  # 
  # tstO$res$lm1
  # tst$res$lm1
  # 
  # tstsclO$res$lm1
  # 

  ampFitList <- lapply( dimnames(tx)[[2]], function(dvs){
    callSTPSP(tn,   DVstr=dvs, gd=T, show=T, showCh=1, who='all', BYstr='exp', FUN=ampFUN, minISI=.190, maxISI=21, startFromScratch=F, eps=T, lmOffset=F)})
  
  ampsclFitList <- lapply( dimnames(tx)[[2]], function(dvs){
    callSTPSP(tn,   DVstr=dvs, gd=T, show=T, showCh=1, who='all', BYstr='exp', FUN=ampFUNscl, minISI=.190, maxISI=21, startFromScratch=F, eps=T, lmOffset=F)})

    
  extractResiduals   <- function(n){ampFitList[[n]]$res$lm1$residuals}
  residMat           <- sapply(1:dim(tx)[2], extractResiduals)
  dimnames(residMat) <- dimnames(tx)
  
  
  extractFitted    <- function(n){ampFitList[[n]]$res$lm1$fitted}
  fitMat           <- sapply(1:dim(tx)[2], extractFitted)
  dimnames(fitMat) <- dimnames(tx)
  
  
  ## ======================================= plot dependence on SOA and AMP for the scalars including the fit
  fLlst <- sapply(1:dim(tx)[2], firstLetter, df=tx)
  Lndx <- which(fLlst %in% c('H','G','B','A','T','D'))
  
  eps.file('Figure5a_Scalars_SOA_AMP_Fit.eps',pdf=T,s=2,mf=c(6,length(Lndx)%/%6),ratio=1)
  par( mfrow=c(length(Lndx)%/%6,6), mar=c(4,4,2,1) )
  disc <- lapply(1:length(Lndx), showTD,  pM=tx[,Lndx], fM=fitMat[,Lndx], xpp=xpp)
  dev.off()
  
  Mndx <- which(fLlst %in% c('M'))
  eps.file('Figure5b_MUA_SOA_AMP_Fit.eps',pdf=T,s=2,mf=c(8,length(Mndx)%/%8),ratio=1)
  par( mfrow=c(length(Mndx)%/%8,8), mar=c(4,4,2,1) )
  disc <- lapply(1:length(Mndx), showTD,  pM=tx[,Mndx], fM=fitMat[,Mndx], xpp=xpp)
  dev.off()
  
  
  
  xpp$valTrial <- xpp$valTrial & xpp$trialType==1
  eps.file('Figure5c_Scalars_Delta_Fit.eps',pdf=T,s=2,mf=c(6,length(Lndx)%/%6),ratio=1)
  par( mfrow=c(length(Lndx)%/%6,6), mar=c(4,4,2,1) )
  disc <- lapply(1:length(Lndx), showTD_delta,  pM=tx[,Lndx], fM=fitMat[,Lndx], xpp=xpp)
  dev.off()
  
  eps.file('Figure5c_Scalars_Fit.eps',pdf=T,s=2,mf=c(6,length(Lndx)%/%6),ratio=1)
  par( mfrow=c(length(Lndx)%/%6,6), mar=c(4,4,2,1) )
  disc <- lapply(1:length(Lndx), showTD_delta,  pM=tx[,Lndx], fM=NULL, xpp=xpp)
  dev.off()
  
  eps.file('Figure5c_MUA_Delta_Fit.eps',pdf=T,s=2,mf=c(8,length(Mndx)%/%8),ratio=1)
  par( mfrow=c(length(Mndx)%/%8,8), mar=c(4,4,2,1) )
  disc <- lapply(1:length(Mndx), showTD_delta,  pM=tx[,Mndx], fM=fitMat[,Mndx], xpp=xpp)
  dev.off()
  
  eps.file('Figure5c_MUA_Delta.eps',pdf=T,s=2,mf=c(8,length(Mndx)%/%8),ratio=1)
  par( mfrow=c(length(Mndx)%/%8,8), mar=c(4,4,2,1) )
  disc <- lapply(1:length(Mndx), showTD_delta,  pM=tx[,Mndx], fM=NULL, xpp=xpp)
  dev.off()
  
  eps.file('Figure5c_MUA_Delta_Fit_only.eps',pdf=T,s=2,mf=c(8,length(Mndx)%/%8),ratio=1)
  par( mfrow=c(length(Mndx)%/%8,8), mar=c(4,4,2,1) )
  disc <- lapply(1:length(Mndx), showTD_delta,  pM=fitMat[,Mndx], fM=NULL, xpp=xpp)
  dev.off()
  
  
  
  #showTD_delta(43, pM=tx, fM=linFitMat,xpp=xpp )
  
  #xpp$ampIndFactor <- as.factor(xpp$ampInd)
  #lm1              <- lm(  tx[,'B_i_1'] ~ log2(xpp$ISI) + xpp$ampIndFactor)
  
  #linFitMat <- fitMat
  #linFitMat[,43]<-lm1$fitted.values
  #showTD(43, pM=tx, fM=linFitMat,xpp=xpp )
  #plot( tapply() )
  
  
  #tst <- xpp[,'B_i_1']
  #image( tapply(lm1$residuals,list(xpp$ampInd[xpp$valTrial], xpp$deltaAmpInd[xpp$valTrial]),mean) )
  #plot(tapply(lm1$residuals,xpp$deltaAmpInd[xpp$valTrial], mean ), type='l' )
}
