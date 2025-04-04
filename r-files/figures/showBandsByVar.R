showBandsByVar <- function(cwbnd,llp, subs=NULL, varStr='isiOct',logtime=F, fileName='Figure4', freqList=NULL, trngList=NULL){
  
  if (is.null(subs))
    subs <- rep(T,dim(llp$xpp)[1])
  
  #browser()
  eps.file(fileName, pdf=T,s=2,mf=c(dim(llp$lfp)[3],length(cwbnd)),ratio=1)
  par(mfcol=c(length(cwbnd),dim(llp$lfp)[3]), mar=c(2,2,1,1))
  #par(mfcol=c(6,6), mar=c(2,2,1,1))
  
  colFUN <- soa.colors
  if (identical(varStr,'ampInd'))
    colFUN <- bluered.colors
  
  Nvar <- length( unique(llp$xpp[[varStr]]) )
  taxis <- llp$taxis
  if(logtime)
    taxis <- llp$logtaxis
  
  for (cx in 1:dim(llp$lfp)[3]){
    for (bx in seq(from=length(cwbnd),to=1,by=-1)){
      
      tlim <- c(3,700)
      if (!logtime){
        tlim <- c(-50,250)
        if (bx%in%c(5,6))
          tlim = c(-10,70)
        
        if (bx%in%c(1,2))
          tlim = c(-150,500)
      }
      
      tndx <- which(llp$taxis>tlim[1] & llp$taxis<tlim[2])
      tlim <- taxis[ which(llp$taxis%in%tlim) ]
      
      vndx <- which(llp$xpp$ampInd==5 & subs)
      mnx  <- apply(cwbnd[[llp$srt[cx]]]$bands[vndx,tndx,bx], 2, mean)
      ylim <- c(0,1.2)*max(mnx)
      
      if (logtime){
        plot(NA, xlim=tlim, ylim=ylim, xaxt='n')
        taxat <- sort( unique(unlist(trngList)) )
        taxnd <- which( llp$taxis%in%c(taxat) )
        axis(1, at=llp$logtaxis[taxnd], taxat  )
        abline(v=llp$logtaxis[taxnd],lty=2,lwd=.5, col='gray')
      }
      if (!logtime){
        plot(NA, xlim=tlim, ylim=ylim)
        abline(h=0)
        abline(v=0)
      }
      for (tx in 1:Nvar){
        vndx <- which( c(llp$xpp[[varStr]])==tx  & subs )
        if (length(vndx)>10)
          mean.plot(cwbnd[[llp$srt[cx]]]$bands[vndx,tndx,bx],taxis[tndx], widthFUN=na.se,add=T, xlim=tlim, col=colFUN(Nvar)[tx])
      }
    }
  }
  dev.off()
}