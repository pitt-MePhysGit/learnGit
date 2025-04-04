call_FLA_SUA <- function(j,df=df,baseDir=baseDir,ind=ind, whoList=whoList, dataSet=dataSet, cmpStrList=cmpStrList, suaMethod=suaMethod){
  
  ## preliminaries (old settings)
  .GlobalEnv$pp   <- function(){paste(baseDir, '/ampOdd_click/word/SUA/plots/',sep='')}
  .GlobalEnv$rda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
  .GlobalEnv$trda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
  .GlobalEnv$spd  <- function(){'/Volumes/Drobo5D3/EEG/sortedZJ_ampOddClick/'}
  
  # figure out if this cell has already been run. 
  # this is just a quick and dirty check for one particular file. 
  # typically it is recommended to set startFromScratch to TRUE and rerun the whole thing
  
  n <- ind[j]
  print(n)
  dfline <- df[n,]
  print(dfline)
  
  
  
  simulate         <- FALSE
  startFromScratch <- TRUE
  alreadyRun       <- file.exists( paste( rda(),'/STPSP/STPSPfit_',dfline$session,'_all_FR_ch_',dfline$channel,'_',dfline$cell,'_',dfline$fextn,'_sig1_gd.rda',sep='') )
  #alreadyRun       <- file.exists( paste(pp(),df$session[n],'/ch_',df$channel[n],'_', df$cell[n],'_',df$fextn[n],'/all/',df$session[n],'_ch_',df$channel[n],'_',df$cell[n],'_',df$fextn[n],'_ALL_GA.pdf',sep=''))
  
  if ( (alreadyRun & !startFromScratch)  ){
    print( paste('Skipping ', dfline$session,' ch_',dfline$channel,'_', dfline$cell,'_',dfline$fextn,sep=''))
  }
  
  if (simulate)
    print( paste('Simulate processing of ', rda(),'/STPSP/STPSPfit_',df$session[n],'_all_FR_ch_',df$channel[n],'_',df$cell[n],'_',df$fextn[n],'_sig1_gd.rda',sep=''))
  
  
  if( (startFromScratch | !alreadyRun)){
    
    
    library(R.matlab)
    library(rgenoud)
    library(lattice)
    library(akima)
    library(Rwave)
    library(MASS)
    
    
    source( paste(baseDir,'r-compile/rUtils/R/rUtils.R',sep='') )
    source( paste(baseDir,'r-compile/ePhysLab/R/ePhysLab.R',sep='') )
    source( paste(baseDir,'r-compile/mph.test/R/mph.test.R',sep='') )
    source( paste(baseDir,'r-compile/mav.test/R/mav.test.R',sep='') )
    
    source( paste(baseDir,'RSkernel/r-files/functions/calcCSD.R',sep='') )
    source( paste(baseDir,'RSkernel/r-files/functions/showCSDimage.R',sep='') )
    source( paste(baseDir,'RSkernel/r-files/functions/showCSD.R',sep='') )
    source( paste(baseDir,'RSkernel/r-files/multiVplot.R',sep='') )
    source( paste(baseDir,'RSkernel/r-files/functions/functions_correctForRapitTransitions.R',sep='') )
    
    source( paste(baseDir,'ampOdd_click/r-files/functions/functions_ampOddclick.R',sep='') )
    source( paste(baseDir,'ampOdd_click/r-files/functions/functions_AOCdual.R',sep='') )
    source( paste(baseDir,'ampOdd_click/r-files/functions/functions_VProbe_STPSP.R',sep='') )
    
    source( paste(baseDir,'ampOdd_click/r-files/functions/prepSTPSPVProbeMUA.R',sep='') )
    source( paste(baseDir,'ampOdd_click/r-files/functions/functions_ampOdd_presynapticPlasticity.R',sep='') )
    source( paste(baseDir,'ampOdd_click/r-files/functions/trStats_ampOddClick.R',sep='') )
    
    source( paste(baseDir,'ampOdd_click/r-files/functions/getSessionAverage.R',sep='') )
    source( paste(baseDir,'ampOdd_click/r-files/functions/getValidTrialsAOC.R',sep='') )
    source( paste(baseDir,'ampOdd_click/r-files/functions/analyzeAOCscalar.R',sep='') )
    source( paste(baseDir,'ampOdd_click/r-files/functions/FLA_SUA.R',sep='') )
    
    source( paste(baseDir,'ampOdd_click/r-files/figures/figure_AOC_timeseries.R',sep='') )
    
    options(warn=0)
    
    
    
    
    
    ## load mua channel  
    chStr     <- paste(     'ch', dfline$startVindex-1+dfline$channel,sep='')
    
    #browser()
    mua         <- loadAOC( dfline , what='mua', chStr=chStr, trda=trda, byBlock=F, suaMethod=suaMethod  )
    mua$cell    <- paste(dfline$cell,dfline$fextn,sep='_')
    mua$ch      <- dfline$channel
    mua$animal  <- dfline$animal
    mua$dataSet <- dfline$dataSet
    mua$xpp$ISI[which(mua$xpp$ISI<0)] <- 20
    
    
    ## load LFP channel
    tmpChNr   <- c(-1, 0, 1) + dfline$channel
    tmpChNr   <- ((tmpChNr%>=%1)%<=%(dfline$numChan))+dfline$startVindex-1
    chStr     <- paste(     'ch', tmpChNr,sep='')
    lfp       <- loadAOC( dfline , what='lfp', chStr=chStr, trda=trda,byBlock=F  )
    
    # create a local csd for the three channels bordering the unit (not valid if unit is on channel 1 or 24)
    csd          <- -2*lfp$lfp[,,2]+lfp$lfp[,,1]+lfp$lfp[,,3]
    dim(csd)     <- c( dim(lfp$lfp)[1:2], 1) 
    lfp$lfp      <- zbind(lfp$lfp, csd, dimension=3)
    
    ## ======================================================
    ## preprocess xpp
    tone                <- mua

    ## isi in steps of 1/3 octaves (1/2 for data set 'l')
    #brks <- exp( log(2) * seq(log(0.2)/log(2),by=1/3,to=4.1) )
    #brks[length(brks)] <- 25
    #brks[1] <- 0.190
    #tone$xpp$binISI <- cut( tone$xpp$ISI, breaks=brks )
    
    tone$xpp$binISI <- tone$xpp$isiOct
    
    tone$xpp$binISI2       <- my.lag(tone$xpp$binISI,1)
    tone$xpp$binISI3       <- my.lag(tone$xpp$binISI,2)
    tone$xpp$binISI4       <- my.lag(tone$xpp$binISI,3)
    tone$xpp$binISI5       <- my.lag(tone$xpp$binISI,4)
    
    tone$xpp$ISI2       <- my.lag(tone$xpp$ISI,1)
    tone$xpp$ISI3       <- my.lag(tone$xpp$ISI,2)
    tone$xpp$ISI4       <- my.lag(tone$xpp$ISI,3)
    tone$xpp$ISI5       <- my.lag(tone$xpp$ISI,4)
    
    #
    tone$xpp$deltaISI  <- c(NA, diff(log2(tone$xpp$ISI)))
    tone$xpp$absDeltaISI  <- abs(tone$xpp$deltaISI)
    
    brks    <- seq(-4.5,by=1,4.5)
    brks[1] <- -10; brks[length(brks)]<- 10 
    tone$xpp$deltaISIbin  <- cut(tone$xpp$deltaISI,  breaks=brks)
    
    brks    <- seq(-0.25,by=0.5,4.5)
    brks[length(brks)]<- 10 
    tone$xpp$absDeltaISIbin  <- cut(tone$xpp$absDeltaISI, breaks=brks)
    #
    
    mua$xpp <- tone$xpp
    lfp$xpp <- tone$xpp
    
    ## =============================
    # identify which trial types were presented
    tblTrTy <- c(NA,NA)
    tblTrTy[1] <- length(which(mua$xpp$trialType==1))     # random
    tblTrTy[2] <- length(which(mua$xpp$trialType==2))     # regular
    tmp        <- tblTrTy > (length(mua$xpp$trialType)/4)
    
    doTrTy        <- c( tmp[1] & tmp[2], tmp)
    names(doTrTy) <- c('all','random','regular')
    whoList       <- names(doTrTy)[which(doTrTy) ]
    
    ## prepare spike counts and psth
    # ==========  compute the psths for each trial separately
    width     <- 0.010
    tmp       <- spk2psth(mua,eval.at=0.001*mua$taxis,width=width,rate=TRUE)
    
    psth      <- tmp
    dim(psth) <- c(dim(psth),1) 
    mua$psth  <- psth
    dimnames(mua$psth)[[2]] <- dimnames(tmp)[[2]]
    mua$epsp  <- baseLineCorrection(mua$psth,c(-0.125,-width), as.numeric(dimnames(mua$psth)[[2]]) )$res
    
    
    for (wind in 1:length(whoList)){ 
      mua$who  <- whoList[wind]
      disc     <- lapply(cmpStrList, function(cs){ FLA_SUA(mua, lfp, cmpStr=cs, eps=TRUE ) } )
    }
  }
  
}
