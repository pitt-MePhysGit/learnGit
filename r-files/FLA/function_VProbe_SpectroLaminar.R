function_VProbe_SpectroLaminar <- function(j,df=df,baseDir=baseDir,trda=trda,ind=ind){
  
  i <- ind[j]
  
  print(df$session[i])
  
  project <- df$project[1]
  
  library(R.matlab)
  library(rgenoud)
  library(lattice)
  library(akima)
  library(Rwave)
  library(MASS)
  library(fastICA)
  
  #source(paste(baseDir, 'RSkernel/r-files/functions/source_all_functions_parallel.R',sep=''))
  source( paste(baseDir,'ampOddClick/r-files/functions/functions_ampOddclick.R',sep='') )
  source( paste(baseDir,'ampOddClick/r-files/functions/functions_AOCdual.R',sep='') )
  source( paste(baseDir,'ampOddClick/r-files/functions/functions_VProbe_STPSP.R',sep='') )
  #source( paste(baseDir,'RSkernel/r-files/functions/FLA_VProbe_STPSP.R',sep='') )
  
  source( paste(baseDir,'RSkernel/r-files/functions/removeCommonReferenceNoise_ica.R',sep='') )
  
  source( paste(baseDir,'r-compile/rUtils/R/rUtils.R',sep='') )
  source( paste(baseDir,'r-compile/mav.test/R/mav.test.R',sep='') )
  source( paste(baseDir,'r-compile/ePhysLab/R/ePhysLab.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/calcCSD.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/showCSDimage.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/functions/showCSD.R',sep='') )
  source( paste(baseDir,'RSkernel/r-files/multiVplot.R',sep='') )
  
  ## 
  dxt <- ifelse(df$drug[i]=='N','', paste('_',df$drug[i],sep=''))
  pxt <- ifelse(df$probe[i]==1 ,'', paste('_',df$probe[i],sep=''))
  
  .GlobalEnv$pp   <- function(){paste(baseDir, 'ampOddClick/word/VProbe/plots/',df$session[i],dxt,pxt,'/',sep='')}
  print(pp())
  
  if (!file.exists(pp()) )
    system(paste('mkdir -p ', pp(),sep='' ))
  
  ## LOAD RAW DATA ====
  chInd      <- df$prelimStartMIndex[i]:df$prelimEndMIndex[i]
  chStr      <- paste('ch',chInd,sep='')
  
  #mua                      <- loadSignalProcessor(  df[i,],chStr=chStr, project=project, epoch='mua', chBlkN = length(chStr))
  #mua$channelPosition$ypos <- -1 * mua$channelPosition$ypos
  #showSignalProcessor(mua, frml=~iciInd, trng=c(-25,150))
  
  lfp                      <- loadSignalProcessor(  df[i,],chStr=chStr, project=project, epoch='llp', chBlkN = length(chStr))
  #lfp$channelPosition$ypos <- -1 * lfp$channelPosition$ypos
  lfp                      <- reRefSignalProcessor(lfp, reRefChan=1:dim(lfp$lfp)[3])
  #showSignalProcessor( lfp, frml=~ampInd, subsExpr = expression(valTone) )
  
  
  ## Calculate FFT for 500 ms long window encompassing the evoked response ====
  lfp   <- fftSignalProcessor(lfp, trng=c(1,500), taperms=25)
  
  
  ## use spectrum ==== 
  avSpc <- apply(lfp$spc[lfp$xpp$valTone,1:50,], c(2,3), mean)
  #image(lfp$faxis[1:50],-24:-1, avSpc[,seq(24,by=-1,1) ], col=bluered.colors(100), )
  
  for (fx in 1:dim(avSpc)[1])
    avSpc[fx,] <- avSpc[fx,]/max(avSpc[fx,])
  
  eps.file('spectroLaminar.eps',eps=T,pdf=T,s=3)
  image(lfp$faxis[1:50],-24:-1, avSpc[,seq(24,by=-1,1) ], col=bluered.colors(100), )
  dev.off()
  
  ## use log spectrum ==== 
  avSpc <- apply(lfp$logspc[lfp$xpp$valTone,1:50,], c(2,3), mean)
  #image(lfp$faxis[1:50],-24:-1, avSpc[,seq(24,by=-1,1) ], col=bluered.colors(100), )
  for (fx in 1:dim(avSpc)[1])
    avSpc[fx,] <- avSpc[fx,]/max(avSpc[fx,])
  
  eps.file('spectroLaminar_log.eps',eps=T,pdf=T,s=3)
  image(lfp$faxis[1:50],-24:-1, avSpc[,seq(24,by=-1,1) ], col=bluered.colors(100), )
  dev.off()
  
  
}