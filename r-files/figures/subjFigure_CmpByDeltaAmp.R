subjFigure_CmpByDeltaAmp <- function(tn,DVstr='n1.r',norm=T,show=T,trialType=c(1,2),eps=T,baseName='cmpByDeltaAmp'){
  
  tn <- subs.data(tn, tn$xpp$trialType%in%trialType)
  
  tn$xpp$DVn <- tn$xpp[[DVstr]]
  lm0           <- lm(tn$xpp[[DVstr]] ~ as.factor(tn$xpp$exp) + as.factor(tn$xpp$ampInd)*tn$xpp$trialType + I( log(tn$xpp$ISI)) )
  lm1           <- lm(tn$xpp[[DVstr]] ~ as.factor(tn$xpp$exp) + as.factor(tn$xpp$ampInd)*tn$xpp$trialType + I( log(tn$xpp$ISI)) +tn$xpp$deltaAmpInd)
  
  lmresult <- anova(lm0,lm1)
    
  if(norm){
    tn$xpp$DVn    <- lm0$residuals
  }
  
  ## AMPLITUDE ADAPTATION in regular condition (less extreme amplitudes for both very lound and very soft tones)
  N1        <- tapply(tn$xpp$DVn, list(tn$xpp$trialType,tn$xpp$deltaAmpInd), mean )
  
  #dbax      <- tapply(tn$xpp$dB,          tn$xpp$ampInd,      mean)
  deltax    <- tapply(tn$xpp$deltaAmpInd, tn$xpp$deltaAmpInd, mean)
  
  if (show){
    thisDir <- paste( pp(),tn$exp,'/FLA/',DVstr,'/',sep='')
    if(!dir.exists(thisDir))
      dir.create(thisDir, recursive=T)    
    
    fileName <- paste(tn$exp,'/FLA/',DVstr,'/',baseName,tn$exp,'.eps',sep='')
    if(eps)
      eps.file(fileName, pdf=T, s=3)
    
    plot(  deltax,N1[1,],col='blue',type='l', ylim=c(-1,1)*na.max(abs(N1)))
    points(deltax,N1[2,],col='red', type='l')
    abline(h=0); abline(v=0)
    
    
    if(eps)
      dev.off()
  }
  
  res <- list( N1=N1, deltax=deltax,lmresult=lmresult,coeff=lm1$coeff ) 
  return(res)
}