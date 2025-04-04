subjFigure_CmpByAmpAndCnd <- function(tn,DVstr='n1.r',norm=T,show=T,eps=T,baseName='cmpByAmpAndCnd'){
  
  tn$xpp$DVn <- tn$xpp[[DVstr]]
  if(norm){
    lm0            <- lm(tn$xpp[[DVstr]] ~ as.factor(tn$xpp$exp) +  I( log(tn$xpp$ISI))  )
    tn$xpp$DVn     <- lm0$residuals
  }
  
  ## AMPLITUDE ADAPTATION in regular condition (less extreme amplitudes for both very lound and very soft tones)
  N1 <- tapply(tn$xpp$DVn, list(tn$xpp$trialType, tn$xpp$dB), mean )
  
  dbax      <- tapply(tn$xpp$dB, tn$xpp$ampInd, mean)
  
  
  if (show){
    thisDir <- paste( pp(),tn$exp,'/FLA/',DVstr,'/',sep='')
    if(!dir.exists(thisDir))
      dir.create(thisDir, recursive=T)    
    
    fileName <- paste(tn$exp,'/FLA/',DVstr,'/',baseName,tn$exp,'.eps',sep='')
    if(eps)
      eps.file(fileName, pdf=T, s=3)
    
    plot(  dbax,N1[1,],col='blue',type='b', pch=20, ylim=na.range(c(N1)) )
    points(dbax,N1[2,],col='red', type='b', pch=20)
    abline(h=0)
  
    if(eps)
      dev.off()
  }
  
  res <- list( N1=N1, xax=dbax ) 
  return(res)
}