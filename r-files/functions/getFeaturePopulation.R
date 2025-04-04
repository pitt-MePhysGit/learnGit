getFeaturePopulation <- function(nd, df, removeBaseline=T,signal=this_signal,epoch=this_epoch,pval=c(10e-10,10e-5)){
  print(nd)
  load( file=paste(rda(),'/',df$session[nd],'/features.rda',sep='') )
  
  firstLetter <- function(nx,df=xpp,Nl=1){ substr(names(df)[nx],Nl,Nl) }
  signalLst   <- sapply(1:dim(xpp)[2], firstLetter, Nl=1)
  epochLst    <- sapply(1:dim(xpp)[2], firstLetter, Nl=3)
  
  procndx     <- which(signalLst %in% signal & epochLst %in% epoch) 
  
  returnStr <- c('trialNr','seqNr','ampInd','trialType','toneType','ISI','deltaPitch','time','p2p',
                 'powEOG','exp','dB','deltaAmpInd','deltaISI','absDeltaISI','deltaISIbin','absDeltaISIbin',
                 'ampDev','timeDev','repNr','isiOct','isiOctFine','futureISI','blockNr','valTrial')#,paste(signal,epoch,sep='_'))
  
  # ================================= process each potential channel
  getDF_DV <- function(px){
    tx                                <- data.frame(xpp[,returnStr],DV=xpp[,px])
    chNr <- strsplit(names(xpp)[px],'_')[[1]][3]
    DVstr                             <- paste(signal,epoch,sep='_')
    
    names(tx)[which(names(tx)=='DV')] <- DVstr 
    tx$exp                            <- df$session[nd]
    
    tx$firstTrial                     <- 1:length(tx$trialNr)
    offset          <- t.test( tx[[DVstr]] )
    tx$ampIndFactor <- as.factor(tx$ampInd)
    tx$isiFactor    <- as.factor(tx$isiOct)
    
    lmx           <- anova( lm(tx[[DVstr]] ~ tx$ampIndFactor + tx$isiFactor))
    
    tx$offset     <- offset$estimate
    tx$pValoff    <- offset$p.value
    tx$pValAmp    <- lmx$'Pr(>F)'[1]
    tx$pValISI    <- lmx$'Pr(>F)'[2]
    
    tx$session    <- df$session[nd]
    tx$channel    <- names(xpp)[px]
    
    if (removeBaseline){
      if( signal %in% c('H','G','B','A','D')){
        if (epoch %in% c('t','e','i','r') ){
          blDVstr <- paste(signal,'b',chNr,sep='_')
          tx[[DVstr]] <- tx[[DVstr]] - xpp[[blDVstr]]   
        }
      }
    }
    
    mnx <- mean(tx[[DVstr]][tx$valTrial])
    
    DVstrNrm       <- paste(signal,epoch,'nrm',sep='_')
    tx[[DVstrNrm]] <- sign(mnx) * tx[[DVstr]]/mnx
    
    return(tx)
  }
  
  DFcmb <- NULL
  if (length(procndx)>1)
    DFcmb <- getDF_DV(procndx[1])
  if (length(procndx)>2){
    for (i in 2:length(procndx))
      DFcmb <- rbind( DFcmb,  getDF_DV(procndx[i]) )
  }
  
  return(DFcmb)
}