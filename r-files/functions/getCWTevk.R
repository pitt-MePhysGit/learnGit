getCWTevk <- function(dat, indx=NA, noctave=8, nvoice=12, w0=pi, chx=1){
  if (is.na(indx[1]))
    indx <- 1:dim(dat$lfp)[1]
  
  sr     <- 1/.001*mean(diff(dat$taxis))
  nF     <- sr/2
  
  taxis <- c(dat$taxis)
  faxis  <- w0/(2*pi) * nF * 2^(seq( to=0, by=1/nvoice, length=nvoice*noctave ))
  faxlog <- seq( to=0, by=1/nvoice, length=nvoice*noctave )
  
  cw     <- list(taxis=taxis, logtaxis=dat$logtaxis, faxis=faxis, faxlog=faxlog, noctave=noctave,nvoice=nvoice, sr=sr, w0=w0, nF=nF)
  cw$mn  <- apply(dat$lfp[indx,,chx],2,na.mean)
  cw$cw  <- cwt(cw$mn, noctave=noctave,  nvoice=nvoice, w0=w0, twoD=TRUE, plot=FALSE)
  
  return(cw)
}