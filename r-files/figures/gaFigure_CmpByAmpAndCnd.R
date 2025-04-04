gaFigure_CmpByDeltaAmpAndAmp <- function(resList,DVstr='n1.r',norm=T,show=T,trialType=c(1),eps=T,baseName='GAcmpByDeltaAmpAndAmp'){
  
  extractN1 <- function(n){return(resList[[n]]$N1) }
  N1mat     <- array( unlist( lapply(1:length(resList), extractN1 )), c( dim(resList[[1]]$N1),length(resList)) )
  N1        <- apply(N1mat,c(1,2),mean)
  N1se        <- apply(N1mat,c(1,2),na.se)
  
  dbax      <- resList[[1]]$xax
  deltax    <- resList[[1]]$deltax
  
  if (show){
    thisDir <- paste( pp(),'g/FLA/',DVstr,'/',sep='')
    if(!dir.exists(thisDir))
      dir.create(thisDir, recursive=T)    
    
    fileName <- paste('g/FLA/',DVstr,'/',baseName,'_g.eps',sep='')
    if(eps)
      eps.file(fileName, pdf=T, s=2, ratio=3/2)
    
    ylim      <- na.range( c(N1) ) + c(-.1,.5)*diff(na.range(c(N1))) 
    
    colList                <- bluered.colors(dim(N1)[1])
    xnd <- 1:dim(N1)[2]
    for (ix in 1:(dim(N1)[1]))
      errorplot( deltax[xnd], N1[ix,xnd], ey=N1se[ix,xnd],  type='b', lwd=2,pch=20, col=colList[ix],ylim=ylim, add=ifelse(ix==1,F,T) )
    
    #zmx       <- 1#3
    #zlim      <- zmx * c(-1,1)# na.max(abs(N1)) 
    #image(dbax, deltax, (N1%<=%zmx)%>=%-zmx, col=bluered.colors(100), zlim=zlim)
    
    
    
    if(eps)
      dev.off()
  }
  
  res <- list( N1=N1, xax=dbax, deltax=deltax ) 
  return(res)
}