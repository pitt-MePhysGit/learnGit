showTD2 <- function(sx,pM=predMat, fM=NULL,xpp=NULL, IVstr1='isiOct', IVstr2='ampInd',colFUN=bluered.colors){
  mnx  <- tapply( pM[xpp$valTrial,sx], list(xpp[[IVstr1]][xpp$valTrial],xpp[[IVstr2]][xpp$valTrial]), mean)
  sex  <- tapply( pM[xpp$valTrial,sx], list(xpp[[IVstr1]][xpp$valTrial],xpp[[IVstr2]][xpp$valTrial]), na.se)
  
  mny  <- tapply( pM[xpp$valTrial,sx], list(xpp[[IVstr1]][xpp$valTrial],xpp[[IVstr2]][xpp$valTrial]), mean)
  sey  <- tapply( pM[xpp$valTrial,sx], list(xpp[[IVstr1]][xpp$valTrial],xpp[[IVstr2]][xpp$valTrial]), na.se)
  
  if(!is.null(fM)){
    mnxf  <- tapply( fM[xpp$valTrial,sx], list(xpp[[IVstr1]][xpp$valTrial],xpp[[IVstr2]][xpp$valTrial]), mean)
    sexf  <- tapply( fM[xpp$valTrial,sx], list(xpp[[IVstr1]][xpp$valTrial],xpp[[IVstr2]][xpp$valTrial]), na.se)
    
    mnyf  <- tapply( fM[xpp$valTrial,sx], list(xpp[[IVstr1]][xpp$valTrial],xpp[[IVstr2]][xpp$valTrial]), mean)
    seyf  <- tapply( fM[xpp$valTrial,sx], list(xpp[[IVstr1]][xpp$valTrial],xpp[[IVstr2]][xpp$valTrial]), na.se)
  }
  
  ylim <- na.range(c(mnx,mny)) + 0.2* c(-1,1)*diff(na.range(c(mnx,mny)))
  plot(1:dim(mnx)[2], rep(0,dim(mnx)[2]), ylim=ylim,  main=dimnames(pM)[[2]][sx])
  if (is.null(fM)){
    for (tx in 1:dim(mnx)[1])
      errorplot(1:dim(mnx)[2],mnx[tx,],ey=sex[tx,], add=T, col=colFUN(dim(mnx)[1])[tx] )
  }
  
  
  if (!is.null(fM)){
    errorplot(1:length(mnx), mnx,ey=sex, ylim=ylim, type='p', main=dimnames(pM)[[2]][sx])#main=scalarList$cmpStrName[sx])
    errorplot(1:length(mny), mny,ey=sey, ylim=ylim, type='p',add=T, col='red')
    points(   1:length(mnxf),mnxf,                  type='l'                 )
    points(   1:length(mnyf),mnyf,                  type='l', col='red'      )
  }
}