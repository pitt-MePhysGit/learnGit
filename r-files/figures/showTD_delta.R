showTD_delta <- function(sx,pM=predMat, fM=NULL, xpp=NULL){
    
  if (is.numeric(fM))
    dim(fM) <- c(length(fM),1)
  
  if (is.numeric(pM))
    dim(pM) <- c(length(pM),1)

  
  txpp <- subset(xpp,xpp$valTrial, drop.levels=T)
  txpp$ampIndFactor <- as.factor(txpp$ampInd)
  txpp$DV <- pM[which(xpp$valTrial),sx]
  
  txpp$logISI      <- log2(txpp$ISI)
  #lm1              <- lm(  txpp$DV ~ txpp$isiOct + txpp$ampIndFactor)
  lm1              <- lm(  txpp$DV ~ txpp$logISI + txpp$ampIndFactor)
  #lm1              <- lm(  pM[which(xpp$valTrial),sx] ~ xpp$isiOct[which(xpp$valTrial)] + xpp$ampIndFactor[which(xpp$valTrial)])
  
  
  if (1==0){
    mnx  <- tapply( txpp$DV, txpp$isiOct, mean)
    mny  <- tapply( txpp$DV, txpp$ampInd, mean)
    
    mnxf  <- tapply( lm1$fitted.values, txpp$isiOct, mean)
    mnyf  <- tapply( lm1$fitted.values, txpp$ampInd, mean)
    
    ylim <- na.range(c(mnx,mny)) + 0.2* c(-1,1)*diff(na.range(c(mnx,mny)))
    
    errorplot(1:length(mnx),mnx, ylim=ylim,lwd=2, type='p', main=dimnames(pM)[[2]][sx])#main=scalarList$cmpStrName[sx])
    errorplot(1:length(mny), mny, ylim=ylim,lwd=2, type='p',add=T, col='red')
    points(1:length(mnx),    mnxf,                 lwd=2,type='l'                 )
    points(1:length(mny),    mnyf,                 lwd=2,type='l', col='red'      )
  }
  
  
  mnx  <- tapply( lm1$residuals, txpp$deltaISIbin, mean)
  sex  <- tapply( lm1$residuals, txpp$deltaISIbin, na.se)
  xax  <- tapply( txpp$deltaISI, txpp$deltaISIbin, mean)
  
  mny  <- tapply( lm1$residuals, txpp$deltaAmpInd, mean)
  sey  <- tapply( lm1$residuals, txpp$deltaAmpInd, na.se)
  yax  <- tapply( txpp$deltaAmpInd, txpp$deltaAmpInd, mean)
  
  
  
  if(!is.null(fM)){
    txpp$DVf <- fM[which(xpp$valTrial),sx]
    #lm2   <- lm(txpp$DVf   ~ txpp$isiOct + txpp$ampIndFactor)
    lm2    <- lm(  txpp$DVf ~ txpp$logISI + txpp$ampIndFactor)
    
    mnxf  <- tapply( lm2$residuals, txpp$deltaISIbin, mean)
    sexf  <- tapply( lm2$residuals, txpp$deltaISIbin, na.se)
    
    mnyf  <- tapply( lm2$residuals, txpp$deltaAmpInd, mean)
    seyf  <- tapply( lm2$residuals, txpp$deltaAmpInd, na.se)
  }
  
  ylim <- na.range(c(mnx,mny)) + 0.2* c(-1,1)*diff(na.range(c(mnx,mny)))
  if (is.null(fM)){
    #errorplot(1:length(mnx),mnx,ey=sex, ylim=ylim, main=scalarList$cmpStrName[sx])
    errorplot(xax,mnx,ey=sex, ylim=ylim,lwd=2,                 main=dimnames(pM)[[2]][sx])
    errorplot(yax,mny,ey=sey, ylim=ylim,lwd=2, add=T, col='red'                          )
  }
  if (!is.null(fM)){
    errorplot(xax, mnx,ey=sex, ylim=ylim,lwd=2, type='p', main=dimnames(pM)[[2]][sx])#main=scalarList$cmpStrName[sx])
    errorplot(yax, mny,ey=sey, ylim=ylim,lwd=2, type='p',add=T, col='red')
    points(   xax, mnxf,                 lwd=2,type='l'                 )
    points(   yax, mnyf,                 lwd=2,type='l', col='red'      )
  }
  abline(h=0)
  abline(v=0)
}