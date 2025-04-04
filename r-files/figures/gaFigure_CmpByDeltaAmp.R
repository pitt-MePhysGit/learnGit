gaFigure_CmpByDeltaAmp <- function(resList,DVstr='n1.r',norm=T,show=T,trialType=c(1,2),eps=T,
                                   col='blue',baseName='GAcmpByDeltaAmp',add=F){
  
  extractN1 <- function(n){return(resList[[n]]$N1) }
  N1mat     <- array( unlist( lapply(1:length(resList), extractN1 )), c( dim(resList[[1]]$N1),length(resList)) )
  N1        <-apply(N1mat,c(1,2),mean)
  
  deltax    <- resList[[1]]$deltax
  
  
  if (show){
    thisDir <- paste( pp(),'g/FLA/',DVstr,'/',sep='')
    if(!dir.exists(thisDir))
      dir.create(thisDir, recursive=T)    
    
    fileName <- paste('g/FLA/',DVstr,'/',baseName,'_g.eps',sep='')
    if(eps)
      eps.file(fileName, pdf=T, s=2, ratio=3/2)
    
    #ylim=c(-1.2,1.2)*na.max(abs(N1[1,]))
    ymx <- 2
    ylim <- ymx * c(-1,1)
    if(!add)
      plot(deltax,rep(-1000,length(deltax)), ylim=ylim, lwd=2)
    
    points(  deltax,N1[1,],col=col,type='l', ylim=ylim, lwd=2)
    addSE <- function(n,...){sebar(N1mat[1,n,], at=deltax[n],...)}
    disc <- sapply(1:length(deltax), addSE, col=col,lwd=2 )
    #points(deltax,N1[2,],col='red', type='l')
    abline(h=0); abline(v=0)
    
    
    if(eps)
      dev.off()
  }
  
  res <- list( N1=N1 ) 
  return(res)
}