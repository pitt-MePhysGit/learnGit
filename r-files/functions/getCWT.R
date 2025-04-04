getCWT <- function(dat, indx=NA, noctave=8, nvoice=12, w0=pi, chx=1, return.raw=F, freqList=NULL){
  if (is.na(indx[1]))
    indx <- 1:dim(dat$lfp)[1]
  
  sr     <- 1/.001*mean(diff(dat$taxis))
  nF     <- sr/2
  
  taxis <- c(dat$taxis)
  faxis  <- w0/(2*pi) * nF * 2^(seq( to=0, by=1/nvoice, length=nvoice*noctave ))
  faxlog <- seq( to=0, by=1/nvoice, length=nvoice*noctave )
  
  cw     <- list(taxis=taxis, logtaxis=dat$logtaxis, faxis=faxis, faxlog=faxlog, noctave=noctave,nvoice=nvoice, sr=sr, w0=w0, nF=nF)
  cw$mn  <- apply(dat$lfp[indx,,chx],2,na.mean)
  cw$cwe <- cwt(cw$mn, noctave=noctave,  nvoice=nvoice, w0=w0, twoD=TRUE, plot=FALSE)
  
  getCWTsingleTrial   <- function(n){ cwt(dat$lfp[n,,chx], noctave=noctave,  nvoice=nvoice, w0=w0, twoD=TRUE, plot=FALSE)}
  
  tst      <- getCWTsingleTrial(1)
  tmp      <- sapply(indx, getCWTsingleTrial )
  dim(tmp) <- c(dim(tst),length(indx))
  tmp      <- aperm(tmp,c(3,1,2))
  
  if (return.raw)
    cw$raw <- tmp
  
  cw$cw    <-     apply( Mod(tmp),     c(2,3), mean)
  if (!is.null(freqList)){# extract time-series for frequecy ranges
    cw$bands                <- array(c(NA), c(dim(tmp)[1:2],length(freqList)) )
    dimnames(cw$bands)[[3]] <- names(freqList)
    for (fx in 1:length(freqList)){
      fndx           <- which(cw$faxis[seq(from=length(cw$faxis),to=1,by=-1)]>=freqList[[fx]][1] & 
                              cw$faxis[seq(from=length(cw$faxis),to=1,by=-1)]<=freqList[[fx]][2] )
      cw$bands[,,fx] <- apply( Mod(tmp[,,fndx]), c(1,2), mean)
    }  
  }
  
  return(cw)
}