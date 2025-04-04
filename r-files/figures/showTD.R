showTD <- function(sx,pM=predMat, fM=NULL, xpp=NULL){
  
  mnx  <- tapply( pM[xpp$valTrial,sx], xpp$isiOct[xpp$valTrial], mean)
  sex  <- tapply( pM[xpp$valTrial,sx], xpp$isiOct[xpp$valTrial], na.se)
  
  mny  <- tapply( pM[xpp$valTrial,sx], xpp$ampInd[xpp$valTrial], mean)
  sey  <- tapply( pM[xpp$valTrial,sx], xpp$ampInd[xpp$valTrial], na.se)
  
  if(!is.null(fM)){
    mnxf  <- tapply( fM[xpp$valTrial,sx], xpp$isiOct[xpp$valTrial], mean)
    sexf  <- tapply( fM[xpp$valTrial,sx], xpp$isiOct[xpp$valTrial], na.se)
    
    mnyf  <- tapply( fM[xpp$valTrial,sx], xpp$ampInd[xpp$valTrial], mean)
    seyf  <- tapply( fM[xpp$valTrial,sx], xpp$ampInd[xpp$valTrial], na.se)
  }
  
  ylim <- na.range(c(mnx,mny)) + 0.2* c(-1,1)*diff(na.range(c(mnx,mny)))
  if (is.null(fM)){
    #errorplot(1:length(mnx),mnx,ey=sex, ylim=ylim, main=scalarList$cmpStrName[sx])
    errorplot(1:length(mnx),mnx,ey=sex, ylim=ylim,lwd=2,                 main=dimnames(pM)[[2]][sx])
    errorplot(1:length(mny),mny,ey=sey, ylim=ylim,lwd=2, add=T, col='red'                          )
  }
  if (!is.null(fM)){
    errorplot(1:length(mnx), mnx,ey=sex, ylim=ylim,lwd=2, type='p', main=dimnames(pM)[[2]][sx])#main=scalarList$cmpStrName[sx])
    errorplot(1:length(mny), mny,ey=sey, ylim=ylim,lwd=2, type='p',add=T, col='red')
    points(   1:length(mnxf),mnxf,                 lwd=2,type='l'                 )
    points(   1:length(mnyf),mnyf,                 lwd=2,type='l', col='red'      )
  }
}