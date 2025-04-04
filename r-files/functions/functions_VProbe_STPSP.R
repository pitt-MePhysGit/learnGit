prepSTPSPVProbeMUA <- function(vmua, trng=c(15,50),minISI=0.190,maxISI=25,DVscl=5,sfc=1,who='all',DVstr='mua'){
  
  muaind          <- which( vmua$taxis>=trng[1] & vmua$taxis<=trng[2] )
  Nchan           <- dim(vmua$lfp)[3]
  
  
  runChannel <- function(chnd){
    vmua$xpp[[DVstr]]     <- apply(vmua$lfp[,muaind,chnd],1,mean)
    
    thisDir   <- paste('/ch_',chnd,sep='')
    tmpDir <- paste(pp(),thisDir,'/STPSP',sep='' )
    if(!exists(tmpDir))
      system(paste('mkdir -p ',tmpDir) )
    
    STPSPFrqFR       <- callFrqSTPSP( vmua, gd=T, startFromScratch=sfc,predFUN=STPSPmodelFun_freq, DVstr=DVstr, who=who,
                                      fileExtn=paste('frq3par_ch',chnd,sep=''),baseName=thisDir,    
                                      minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=T, eps=T)
    
    STPSPfit   <- callSTPSP(vmua, gd=T, startFromScratch=sfc, DVstr=DVstr,who=who, minFrq=500, toneBy=1,NOct=3,showCh=chnd,baseName=thisDir,
                            minISI=minISI, maxISI=maxISI, DVsdscale=DVscl, show=T, eps=T, fileExtn=paste('ch',chnd,sep=''))
    print( c(log(2)/STPSPfit$modelFit$par[1],STPSPfit$model$par[2]) )
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
  points(rep(6,Nchan),seq(1,by=-.025,length=24), col=colList,pch=20,cex=.5)
  return(resList)
}