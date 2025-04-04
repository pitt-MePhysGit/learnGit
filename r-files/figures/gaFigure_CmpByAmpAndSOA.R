gaFigure_CmpByAmpAndSOA <- function(resList,DVstr='n1.r',norm=T,show=T,trialType=c(1),eps=T,baseName='GAcmpByAmpAndSOA'){
  
  extractN1 <- function(n){return(resList[[n]]$N1) }
  N1mat     <- array( unlist( lapply(1:length(resList), extractN1 )), c( dim(resList[[1]]$N1),length(resList)) )
  N1        <-apply(N1mat,c(1,2),na.mean)
  N1se      <-apply(N1mat,c(1,2),na.se)
  
  dbax      <- resList[[1]]$xax
  isiax      <- resList[[1]]$yax
  vi        <- 1:5
  
  if (show){
    thisDir <- paste( pp(),'g/FLA/',DVstr,'/',sep='')
    if(!dir.exists(thisDir))
      dir.create(thisDir, recursive=T)    
    
    fileName <- paste('g/FLA/',DVstr,'/',baseName,'_g.eps',sep='')
    if(eps)
      eps.file(fileName, pdf=T, s=2, ratio=3/2)
    
    #par(mfcol=c(1,2))
    col <- soa.colors(6)
    errorplot(dbax,N1[,vi[1]], ey=N1se[,vi[1]], col='black', ylim=c(0,2),pch=20)
    for (sx in 1:length(vi))
      errorplot(dbax,N1[ ,vi[sx]], ey=N1se[,vi[sx]], col=col[sx], add=T,pch=20)
    if(eps)
      dev.off()
    
    fileName <- paste('g/FLA/',DVstr,'/',baseName,'2_g.eps',sep='')
    if(eps)
      eps.file(fileName, pdf=T, s=2, ratio=3/2)
    
    col <- bluered.colors(5)
    errorplot(log2(isiax[vi]),N1[1,vi], ey=N1se[1,vi], col='black', ylim=c(0,2),pch=20)
    for (ax in 1:length(dbax))
      errorplot(log2(isiax[vi]),N1[ax,vi], ey=N1se[ax,vi], col=col[ax], add=T,pch=20)
    
    if(eps)
      dev.off()
  }
  
  res <- list( N1=N1, xax=dbax, yax=isiax ) 
  return(res)
}