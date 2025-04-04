function_aoc_featureDecoding <- function(ind,Nica=6){
  
  baseDir <- '~/Dropbox/'
  library(R.matlab)
  library(rgenoud)
  library(lattice)
  library(akima)
  library(Rwave)
  library(MASS)
  library(fastICA)
  library(e1071)
  library(circlize)
  
  
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
#  source( paste(baseDir,'aMMN_controlled/r-files/functions_eCog.R',sep='') )
  
  
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
  # 1) boot-strapped data frame: draw responses for different features from different trials
  # 2) estimate d prime for each feature 
  # 3) decoding of hypothetical recording of all vprobe recordings simulataneously
  # 4) decoding for the few examples of actual two simultaneous Vprobe recordings
  
  ##====================================================== load example ICA component
  .GlobalEnv$pp    <- function(){paste(baseDir, 'ampOdd_click/word/VProbe/plots/',df$exp[ind],'/ICA_Npca',Nica,'/',sep='')}
  pp <- .GlobalEnv$pp
  if (!exists(pp()))
    system(paste('mkdir -p ', pp(), sep=''))
  
  if (!exists(paste( pp(), '/STPSP/', sep='')))
    system(paste('mkdir -p ', pp(),'/STPSP/', sep=''))
  
  
  load( file=paste(rda(),'/',df$exp[ind],'/features.rda',sep='') )
  
  ## ========================================== first pass SVM-based decoding
  firstLetter <- function(nx,df=xpp){ substr(names(df)[nx],1,1) }
  fLlst <- sapply(1:dim(xpp)[2], firstLetter)
  
  tx          <- xpp[, fLlst %in% c('H','G','B','A','T','D','M')  ]
  
  eps.file('Figure7a_Decoding_all.eps',pdf=T,s=2,mf=c(2,1),ratio=1)
  par( mfrow=c(1,2), mar=c(4,4,2,1) )
  
  tx$category <- as.factor(xpp$ampInd)
  svmAMP <- decodeSVM(tx,subs=xpp$ISI<16)
  
  tx$category <- as.factor(xpp$isiOct)
  svmISI <- decodeSVM(tx,subs=xpp$ISI<8)
  dev.off()
  
  ## ==================================== keeping amplitude stable
  Namp <- length(unique(xpp$ampInd))
  Nsoa <- length(unique(xpp$isiOct))
  
  tx          <- xpp[, fLlst %in% c('H','G','B','A','T','D','M')  ]
  
  
  eps.file('Figure7a_Decoding_byAMPandISI_confusion.eps',pdf=T,s=2,mf=c(2,1),ratio=1)
  par( mfrow=c(2,5), mar=c(4,4,2,1) )
  tx$category <- as.factor(xpp$isiOct)
  svmISIlist <- lapply(1:Namp, function(ax){ decodeSVM(tx, subs=xpp$ISI<8 & xpp$ampInd==ax) })
  
  tx$category <- as.factor(xpp$ampInd)
  svmAMPlist <- lapply(1:(Nsoa-2), function(sx){ decodeSVM(tx, subs=xpp$ISI<8 & c(xpp$isiOct)==sx) })
  dev.off()
  
  extractScalar <- function(nx,what='mnErrorR',lst=svmAMPlist){ lst[[nx]][[what]]}
  
  mnErrorRamp <- sapply(1:Namp,  extractScalar, what='mnErrorR',lst=svmAMPlist )
  mnErrorEamp <- sapply(1:Namp,  extractScalar, what='mnErrorE',lst=svmAMPlist )
  
  mnErrorRisi <- sapply(1:(Nsoa-2),  extractScalar, what='mnErrorR',lst=svmISIlist )
  mnErrorEisi <- sapply(1:(Nsoa-2),  extractScalar, what='mnErrorE',lst=svmISIlist )
  
  
  eps.file('Figure7a_Decoding_byAMPandISI.eps',pdf=T,s=2,mf=c(2,1),ratio=1)
  par(mfcol=c(1,2))
  plot(mnErrorRamp, ylim=c(0,2), type='b', pch=20, col='gray' )
  abline( h=svmAMP$mnErrorR, col='gray')
  points(mnErrorEamp, type='b',pch=20, col='blue' )
  abline( h=svmAMP$mnErrorE, col='blue')
  
  plot(mnErrorRisi, ylim=c(0,2), type='b', pch=20, col='gray' )
  abline( h=svmISI$mnErrorR, col='gray')
  points(mnErrorEisi, type='b',pch=20, col='blue' )
  abline(h=svmISI$mnErrorE, col='blue')
  dev.off()
  
  
  
  
  
  ## ================================================  simulation
  simdf          <- xpp[, c('ampInd','isiOct')  ]
  sd             <- 2
  simdf$S1       <-  .5*simdf$ampInd  +  1*c(simdf$isiOct) + rnorm( dim(simdf)[1],m=0,s=sd ) 
  simdf$S2       <-   1*simdf$ampInd  + .5*c(simdf$isiOct) + rnorm( dim(simdf)[1],m=0,s=sd ) 
  simdf$S3       <-  .5*simdf$ampInd  + .5*c(simdf$isiOct) + rnorm( dim(simdf)[1],m=0,s=sd ) 
  simdf$S4       <-  .5*simdf$ampInd  +  1*c(simdf$isiOct) + rnorm( dim(simdf)[1],m=0,s=sd ) 
  simdf$S5       <-   1*simdf$ampInd  + .5*c(simdf$isiOct) + rnorm( dim(simdf)[1],m=0,s=sd ) 
  simdf$S6       <-  .5*simdf$ampInd  + .5*c(simdf$isiOct) + rnorm( dim(simdf)[1],m=0,s=sd ) 
  simdf$S7       <-  .5*simdf$ampInd  +  1*c(simdf$isiOct) + rnorm( dim(simdf)[1],m=0,s=sd ) 
  simdf$S8       <-   1*simdf$ampInd  + .5*c(simdf$isiOct) + rnorm( dim(simdf)[1],m=0,s=sd ) 
  simdf$S9       <-  .5*simdf$ampInd  + .5*c(simdf$isiOct) + rnorm( dim(simdf)[1],m=0,s=sd ) 
  simdf$S10       <-  .5*simdf$ampInd  +  1*c(simdf$isiOct) + rnorm( dim(simdf)[1],m=0,s=sd ) 
  simdf$S11       <-   1*simdf$ampInd  + .5*c(simdf$isiOct) + rnorm( dim(simdf)[1],m=0,s=sd ) 
  simdf$S12       <-  .5*simdf$ampInd  + .5*c(simdf$isiOct) + rnorm( dim(simdf)[1],m=0,s=sd ) 
  
  fLlst <- sapply(1:dim(simdf)[2], firstLetter, df=simdf)
  tx          <- simdf[, fLlst %in% c('S')  ]
  mx <- as.matrix(tx)
  
    
  par( mfrow=c(1,2), mar=c(4,4,2,1) )
  
  tx$category <- as.factor(simdf$ampInd)
  svmAMP <- decodeSVM(tx,subs=xpp$ISI<8)
  
  tx$category <- as.factor(xpp$isiOct)
  svmISI <- decodeSVM(tx,subs=xpp$ISI<8)

  par( mfrow=c(4,3), mar=c(4,4,2,1) )
  disc <- lapply(1:dim(mx)[2], showTD,  pM=mx, fM=NULL, xpp=xpp)
  
  
  
  
  if (1==0){
    ## ==================================== two categories simulataneously
    tx$category <- as.factor( tapply( 1:dim(xpp)[1], list(xpp$isiOct,xpp$ampInd) ))
    decodeSVM(tx,subs=xpp$ISI<8)
    
    ## ========================================== decoding on dimensionally reduced data pca
    tx    <- xpp[, fLlst %in% c('H','G','B','A','T','D','M')  ]
    pcatx <- pca( tx, ncomp=190 )
    
    plot(pcatx$explained_variance,type='l')
    
    px <- data.frame( pcatx$x )
    
    par( mfrow=c(3,4), mar=c(4,4,2,1) )
    disc <- lapply(1:4, showTD,  pM=pcatx$x, fM=NULL, xpp=xpp)
    disc <- lapply(1:4, showTD2,  pM=pcatx$x, fM=NULL, xpp=xpp,IVstr1='isiOct', IVstr2='ampInd')
    disc <- lapply(1:4, showTD2,  pM=pcatx$x, fM=NULL, xpp=xpp,IVstr1='ampInd', IVstr2='isiOct')
    
    eps.file('Figure7b_Decoding_pca.eps',pdf=T,s=2,mf=c(2,1),ratio=1)
    par( mfrow=c(1,2), mar=c(4,4,2,1) )
    
    px$category <- as.factor(xpp$ampInd)
    decodeSVM(px,subs=xpp$ISI<16)
    
    px$category <- as.factor(xpp$isiOct)
    decodeSVM(px,subs=xpp$ISI<8)
    dev.off()
    
    
    
    
    ## ========================================== decoding on dimensionally reduced data ica
    tx    <- xpp[, fLlst %in% c('H','G','B','A','T','D','M')  ]
    icatx <- fastICA( tx, n.comp=4 )
    
    ix <- data.frame( icatx$S )
    
    par( mfrow=c(3,4), mar=c(4,4,2,1) )
    disc <- lapply(1:4, showTD,   pM=icatx$S, fM=NULL, xpp=xpp)
    disc <- lapply(1:4, showTD2,  pM=icatx$S, fM=NULL, xpp=xpp,IVstr1='isiOct', IVstr2='ampInd')
    disc <- lapply(1:4, showTD2,  pM=icatx$S, fM=NULL, xpp=xpp,IVstr1='ampInd', IVstr2='isiOct')
    
    eps.file('Figure7b_Decoding_ica.eps',pdf=T,s=2,mf=c(2,1),ratio=1)
    par( mfrow=c(1,2), mar=c(4,4,2,1) )
    
    ix$category <- as.factor(xpp$ampInd)
    decodeSVM(ix,subs=xpp$ISI<16)
    
    ix$category <- as.factor(xpp$isiOct)
    decodeSVM(ix,subs=xpp$ISI<8)
    dev.off()
  }
  
  
}
