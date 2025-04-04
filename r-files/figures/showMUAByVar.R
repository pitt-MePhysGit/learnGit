showMUAByVar <- function(mua, subs=NULL, varStr='isiOct',logtime=F, fileName='Figure4', freqList=NULL, trngList=NULL){
  
  if (is.null(subs))
    subs <- rep(T,dim(mua$xpp)[1])
  

  Nchan <- dim(mua$lfp)[3]
  
  eps.file(fileName, pdf=T,s=2,mf=c(8,Nchan%/%8),ratio=1)
  par(mfrow=c(Nchan%/%8,8 ), mar=c(2,2,1,1))
  #par(mfcol=c(6,6), mar=c(2,2,1,1))
  
  colFUN <- soa.colors
  if (identical(varStr,'ampInd'))
    colFUN <- bluered.colors
  
  Nvar <- length( unique(mua$xpp[[varStr]]) )
  
  taxis <- mua$taxis
  if(logtime)
    taxis <- mua$logtaxis
  
  for (cx in 1:dim(mua$lfp)[3]){
    
    tlim <- c(3,700)
    if (!logtime)
      tlim <- c(-50,150)
    
    #tndx <- which(mua$taxis>tlim[1] & mua$taxis<tlim[2])
    #tlim <- range( taxis[ mua$taxis[tndx] ] )
    
    tndx <- which(mua$taxis>tlim[1] & mua$taxis<tlim[2])
    tlim <- taxis[ which(mua$taxis%in%tlim) ]
    
    
    vndx <- which(mua$xpp$ampInd==5 & subs)
    mnx  <- apply( mua$lfp[vndx,tndx,cx], 2, mean)
    ylim <- c(-.2,1.2)*max(mnx)
    
    if (logtime){
      plot(NA, xlim=tlim, ylim=ylim, xaxt='n')
      taxat <- sort( unique(unlist(trngList)) )
      taxnd <- which( mua$taxis%in%c(taxat) )
      axis(1, at=mua$logtaxis[taxnd], taxat  )
      abline(v=mua$logtaxis[taxnd],lty=2,lwd=.5, col='gray')
    }
    if (!logtime){
      plot(NA, xlim=tlim, ylim=ylim)
      abline(h=0)
      abline(v=0)
    }
    for (tx in 1:Nvar){
      vndx <- which( c(mua$xpp[[varStr]])==tx  & subs )
      if (length(vndx)>10)
        mean.plot(mua$lfp[vndx,tndx,cx],taxis[tndx], widthFUN=na.se,add=T, xlim=tlim, col=colFUN(Nvar)[tx])
    }
  }
  dev.off()
}