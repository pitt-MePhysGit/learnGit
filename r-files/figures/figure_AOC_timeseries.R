figure_AOC_timeseries <- function(tone, whoStr='lfp', chnd=1, rngT=c(-10,150) ){
  
  if (length(tone$xpp$valInd)==0)
    tone$xpp$valInd <- c(T)
  
  valInd <- tone$xpp$valInd
  
  SOAval <- sort(unique( tapply(tone$xpp$ISI, tone$xpp$isiOct, mean)))
  Nsoa   <- length(SOAval)
  
  AMPval <- sort(unique( tapply(tone$xpp$dB, tone$xpp$ampInd, mean)))
  Namp   <- length(AMPval)
  
  tsdata <- tone[[whoStr]][,,chnd]
  dim(tsdata)
  
  ## ===========================================================================
  ## little helper function that returns valid mean psths even
  ## if the number of repetitions is one or zero
  tmeanPSTH <- function(psth, FUN=mean){
    
    if( dim(psth)[1]==0 )
      res <- rep(NA, dim(tsdata)[2] )  
    
    if( dim(psth)[1]==1 )
      res <- psth[,1]
    
    
    if( dim(psth)[1]>1 )
      res <- apply(psth,c(2),FUN)
    
    return(res)
  }
  
  # ========================================== main effect ISI
  binKerISI       <- tapply( tone$xpp$ISI, tone$xpp$isiOct )  
  NtrISI          <- table(binKerISI)
  mnISI           <- tapply(1:dim(tsdata)[1], list(binKerISI),function(ind){ tmeanPSTH( tsdata[ ind[ind%in%which(valInd)] ,,drop=F]) })
  print(NtrISI)
  
  # ========================================= main effect AMP
  binKerAMP       <- tone$xpp$ampInd 
  NtrAMP          <- table(binKerAMP)
  mnAMP           <- tapply(1:dim(tsdata)[1], list(binKerAMP),function(ind){ tmeanPSTH( tsdata[ ind[ind%in%which(valInd)] ,,drop=F]) })
  print(NtrAMP)
  
  
  # ========================================= interaction ISI and AMP
  NtrISIAMP       <- table(binKerISI,binKerAMP)
  mnISIAMP        <- tapply(1:dim(tsdata)[1], list(binKerISI,binKerAMP),function(ind){ tmeanPSTH( tsdata[ ind[ind%in%which(valInd)] ,,drop=F]) })
  print(NtrISIAMP)
  
  
    
  ## ========================================= make figure  
  #     ISI on the x-axis
  #     AMP on the y-axis
  
  rngY        <- na.range( c(mnISIAMP)  )
  thisSquY    <- function(x){ squeeze(x,to.range=c(0.1,.9), from.range=rngY) }
  
  tmpY        <- pretty( c(0, max(abs(rngY))), n=2 )
  prtY        <- c(0, tmpY[ max(which(tmpY<rngY[2]))]  )
    
  thisSquT    <- function(x){ squeeze(x,to.range=c(0.1,.9), from.range=rngT) }
  tnd         <- which( tone$taxis >= rngT[1] & tone$taxis <= rngT[2])
  
  
  ## ---------------------------------------- set up the plot window
  plot(NULL, xlim=c(0.5, Nsoa+2.75), ylim=c(0.5, Namp+2.75), bty='n', xaxt='n', yaxt='n', xlab='ISI [sec]', ylab='Intensity [dB]' )
  axis(1, at=0.5+(1:Nsoa))
  axis(2, at=0.5+(1:Namp), labels=AMPval)
  
  
  ## ----------------------------------------  adding the interaction plots
  tdf <- expand.grid(ISIbin = sort(unique(binKerISI)), AMPbin = sort(unique(binKerAMP)) )
  
  addPSTH <- function(nx){
    psth <- tsdata[ which(valInd & binKerISI==tdf$ISIbin[nx] & binKerAMP==tdf$AMPbin[nx]),,drop=F]
    points( thisSquT(tone$taxis[tnd])+tdf$ISIbin[nx], thisSquY(0*tnd)+tdf$AMPbin[nx], type='l',lwd=.5, col='gray' )
    points( thisSquT(c(0,0))+tdf$ISIbin[nx], thisSquY(prtY)+tdf$AMPbin[nx], type='l',lwd=.5, col='gray' )
    
    if( dim(psth)[1]>1 )
      mean.plot( thisSquY(psth[,tnd])+tdf$AMPbin[nx], thisSquT(tone$taxis[tnd])+tdf$ISIbin[nx], col='black', widthFUN=na.se)
  }
  disc <- lapply(1:dim(tdf)[1], addPSTH )
  
  ## ----------------------------------------  adding the main effect ISI
  tdf <- data.frame(ISIbin = sort(unique(binKerISI)), AMPbin = c(Namp+1.25) )
  addPSTH <- function(nx){
    psth <- tsdata[ which(valInd & binKerISI==tdf$ISIbin[nx] ),,drop=F]
    points( thisSquT(tone$taxis[tnd])+tdf$ISIbin[nx], thisSquY(0*tnd)+tdf$AMPbin[nx], type='l',lwd=.5, col='gray' )
    points( thisSquT(c(0,0))+tdf$ISIbin[nx], thisSquY(prtY)+tdf$AMPbin[nx], type='l',lwd=.5, col='gray' )
    
    if( dim(psth)[1]>1 )
      mean.plot( thisSquY(psth[,tnd])+tdf$AMPbin[nx], thisSquT(tone$taxis[tnd])+tdf$ISIbin[nx], col=soa.colors(Nsoa+1)[tdf$ISIbin[nx]], widthFUN=na.se, lwd=1.5)
  }
  disc <- lapply(1:dim(tdf)[1], addPSTH )
  
  ## ----------------------------------------  adding the main effect AMP
  tdf <- data.frame(ISIbin = c(Nsoa+1.25), AMPbin = sort(unique(binKerAMP)) )
  addPSTH <- function(nx){
    psth <- tsdata[ which(valInd & binKerAMP==tdf$AMPbin[nx] ),,drop=F]
    points( thisSquT(tone$taxis[tnd])+tdf$ISIbin[nx], thisSquY(0*tnd)+tdf$AMPbin[nx], type='l',lwd=.5, col='gray' )
    points( thisSquT(c(0,0))+tdf$ISIbin[nx], thisSquY(prtY)+tdf$AMPbin[nx], type='l',lwd=.5, col='gray' )
    
    if( dim(psth)[1]>1 )
      mean.plot( thisSquY(psth[,tnd])+tdf$AMPbin[nx], thisSquT(tone$taxis[tnd])+tdf$ISIbin[nx], col=bluered.colors(Namp)[tdf$AMPbin[nx]], widthFUN=na.se,lwd=1.5)
  }
  disc <- lapply(1:dim(tdf)[1], addPSTH )
  
  
  # dummy panel top right --------------------------
  points( thisSquT(tone$taxis[tnd])+ Nsoa+1.25, thisSquY(0*tnd)+Namp+1.25, type='l',lwd=.5, col='black' )
  points( thisSquT(c(0,0))         + Nsoa+1.25, thisSquY(prtY) +Namp+1.25, type='l',lwd=.5, col='black' )
  
  
  x <- thisSquT( 0 )       + Nsoa+1.25
  y <- thisSquY( prtY[2] ) + Namp+1.25
  text( x,y, labels=prtY[2], adj=c(1.5,.5) , cex=.5   )  
  
  x <- thisSquT( rngT[2] ) + Nsoa+1.25
  y <- thisSquY( 0 )       + Namp+1.25
  text( x,y, labels=rngT[2], adj=c(.5,1.5) , cex=.5   )  
  
}
