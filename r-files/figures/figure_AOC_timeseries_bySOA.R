figure_AOC_timeseries_bySOA <- function(tone, whoStr='lfp', chnd=1, rngT=c(-10,150), spreadOut=F, rngY=NA ){
  
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
  tmeanPSTH <- function(psth, FUN=na.mean){
    
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
  
  ## ========================================= make figure  
  #     ISI on the x-axis
  #     AMP on the y-axis
  panelStep <- 0  # spread out plots in different panels?
  to.range    <- c(0,1)
  
  if (spreadOut){
    panelStep <- 1
    to.range  <- c(.1,.9)
  }
  
  if(is.na(rngY)){
    rngY        <- na.range( c(mnISI)  )
    rngY        <- rngY + c(-.1,.1)*diff(rngY)
  }
  thisSquY    <- function(x){ squeeze(x,to.range=to.range, from.range=rngY) }
  
  tmpY        <- pretty( c(0, max(abs(rngY))), n=2 )
  print(tmpY)
  prtY        <- c(0, tmpY[ max(which(tmpY<rngY[2]))]  )
  
  thisSquT    <- function(x){ squeeze(x,to.range=to.range, from.range=rngT) }
  tnd         <- which( tone$taxis >= rngT[1] & tone$taxis <= rngT[2])
  
  ## ---------------------------------------- set up the plot window
    
  plot(NULL, xlim=c(1, 2), ylim=c(1,2), bty='o', xaxt='n', yaxt='n', xlab='ISI [sec]', ylab='Intensity [dB]' )
  axis(1, at=1+thisSquT(pretty(rngT)), pretty(rngT)  )
  axis(2, at=1+thisSquY(pretty(rngY)), pretty(rngY)  )
  
  if (spreadOut){
    plot(NULL, xlim=c(0.5, panelStep*Nsoa+1.75), ylim=c(0.9, 2.1), bty='n', xaxt='n', yaxt='n', xlab='ISI [sec]', ylab='Intensity [dB]' )
    axis(1, at=0.5+(1:Nsoa)  )
    axis(2, at=0.5+(1:Namp), labels=AMPval)
  }
  ## ----------------------------------------  adding the main effect ISI
  tdf <- data.frame(ISIbin = sort(unique(binKerISI)), AMPbin = c(1) )
  addPSTH <- function(nx){
    psth <- tsdata[ which(valInd & binKerISI==tdf$ISIbin[nx] ),,drop=F]
    if (spreadOut){
      points( thisSquT(tone$taxis[tnd])+1+panelStep*tdf$ISIbin[nx], thisSquY(0*tnd)+tdf$AMPbin[nx], type='l',lwd=.5, col='gray' )
      points( thisSquT( c(0,0))        +1+panelStep*tdf$ISIbin[nx], thisSquY(prtY)+tdf$AMPbin[nx], type='l',lwd=.5, col='gray' )
    }
    if (!spreadOut)
      abline( h=1+thisSquY(0) )
      
    if( dim(psth)[1]>1 )
      mean.plot( thisSquY(psth[,tnd])+tdf$AMPbin[nx], thisSquT(tone$taxis[tnd])+1+panelStep*tdf$ISIbin[nx], col=soa.colors(Nsoa+1)[tdf$ISIbin[nx]], widthFUN=na.se, lwd=1.5)
  }
  disc <- lapply(1:dim(tdf)[1], addPSTH )
  
  
}
