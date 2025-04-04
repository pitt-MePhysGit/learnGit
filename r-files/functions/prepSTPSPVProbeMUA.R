prepSTPSPVProbeMUA <- function(vmua,minISI=0.190, maxISI=25, DVscl=5, sfc=1, who='all',parline=NA, dfline=NA){

  Nchan      <- dim(vmua$lfp)[3]
  tnd        <- which(vmua$taxis >= parline$tmin & vmua$taxis < parline$tmax)
  DVstr      <- parline$cmpStrName
  
  #vmua$xpp$toneIndFactor <- as.factor(vmua$xpp$toneInd)
  
  vmua$xpp$TMP         <- apply(vmua$lfp[,tnd,1], c(1), mean)
  nd                   <- which(names(vmua$xpp)=='TMP')
  names(vmua$xpp)[nd]  <- DVstr

  runChannel <- function(chnd){
    
    vmua$xpp[[DVstr]]     <- apply(vmua$lfp[,tnd,chnd],1,mean)
    
    
    if (identical(parline$what,'muaR')){
      fileName <- paste( rda(), dfline$session, '/MUAbl_brkpoints.rda',sep='' )
      
      if (file.exists(fileName)){# use previously saved break points to remove sharp artifacts
        print(fileName)
        load(file=fileName)# load variable 'brks'
        
        tb   <- unique( c(0,brks,length(vmua$xpp[[DVstr]])))
        x    <- 1:length(vmua$xpp[[DVstr]])
        xbin <- cut(x, tb )
        
        lm0 <- lm(vmua$xpp[[DVstr]] ~ 1)
        lm1 <- lm(vmua$xpp[[DVstr]] ~ 1 + xbin)
        pVal <- anova(lm0,lm1)$P[2]           
        
        if( pVal<0.001 )
          vmua$xpp[[DVstr]] <- lm1$residuals
      }
      
      # correct for slow fluctuations
      if (1==1){ #(dfline$corMUAbl==2){
        print('correcting for slow fluctuations')
        mavDV            <- get.mav(vmua$xpp[[DVstr]], 1:dim(vmua$xpp)[1], eval.at=seq(1,dim(vmua$xpp)[1],by=1 ), show=T, distribution='Gaussian',boot=F, width=100)
        vmua$xpp[[DVstr]] <- vmua$xpp[[DVstr]] - mavDV$mav + mean(mavDV$mav)
      }
    }
    
    
    #if (identical(parline$cmpStrName,'MUAbl')){# try correcting for jumps in the baseline of muaR
    #  tst <- callCorrectForRapidTransitions(vmua$xpp[[DVstr]])
    #  vmua$xpp[[DVstr]] <- tst$corrected.values
    #}
    
    
    thisDir   <- paste('/ch_',chnd,sep='')
    tmpDir <- paste(pp(),thisDir,'/STPSP',sep='' )
    if(!exists(tmpDir))
      system(paste('mkdir -p ',tmpDir) )
    
    #ampFit <- callSTPSP(vmua,                 DVstr=DVstr, gd=T, show=TRUE, showCh=1, who=who, BYstr='exp', 
    #                    FUN=cumsumamp, minISI=minISI, maxISI=maxISI, startFromScratch=F)
    #ampFit$modelFit$par
    
    STPSPfit   <- callSTPSP(vmua, gd=T, startFromScratch=sfc, DVstr=DVstr,who=who, BYstr='ampInd',
                            showCh=chnd,baseName=thisDir,FUN=cumsumamp, 
                            minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=T, eps=T, fileExtn=paste('ch',chnd,sep=''))
    print( c(log(2)/STPSPfit$modelFit$par[1],STPSPfit$model$par[2]) )
    
    #plot( vmua$xpp[[DVstr]][STPSPfit$res$yind], col=vmua$xpp$trialType[STPSPfit$res$yind] )
    #points(STPSPfit$res$lm1$fitted.values)
    
    # seed parameters for STPSPFrqFR
    #sp         <- rep(NA,6)
    #sp[c(1,2)] <- STPSPfit$modelFit$par
    #sp[c(4,5)] <- RFfit$modelFit$par[c(1,2)]
    
    #STPSPFrqFR       <- callFrqSTPSP( vmua, gd=T, startFromScratch=sfc,predFUN=STPSPmodelFun_freq, DVstr=DVstr, who=who,
    #                                  fileExtn=paste('frq3par_ch',chnd,sep=''),baseName=thisDir,sp=sp,    
    #                                  minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=T, eps=T)
    #
    #STPSPFrqFR$model$par
    
    return(STPSPfit)
  }
  resList <- sapply(1:Nchan, runChannel)
  
  extractPar <- function(n,m=1){return(  resList[[ 2*n-1  ]]$par[m])   }
  lambda <- log(2)/sapply(1:Nchan, extractPar, m=1)
  exFrac <- sapply(1:Nchan, extractPar, m=2)    
  
  extractPar <- function(n,m=1){return(  resList[[ 2*n-1  ]]$modelTest$Pr[2])   }
  pval       <- -log(sapply(1:Nchan, extractPar, m=1))/log(10)
  pvalNrm   <- (pval%<=%20)
  pvalNrm   <- pvalNrm%<=a=0%4
  pvalNrm   <- pvalNrm/10
  
  depthcol <- colorRampPalette(c("blue","gray","red"))
  colList <- depthcol(Nchan)
  plot(lambda%<=%5,exFrac,cex=pvalNrm, col=colList,pch=20, xlim=c(0,6),ylim=c(0,1))
  points(rep(6,Nchan),seq(1,by=-.025,length=Nchan), col=colList,pch=20,cex=.5)
  return(resList)
}