subjFigure_CmpByDeltaAmpAndAmp <- function(tn,DVstr='n1.r',norm=T,show=T,trialType=c(1),eps=T,baseName='cmpByDeltaAmpAndAmp'){
  
  tn <- subs.data(tn, tn$xpp$trialType%in%trialType)
  
  tn$xpp$DVn <- tn$xpp[[DVstr]]
  
  lm0           <- lm(tn$xpp[[DVstr]] ~ as.factor(tn$xpp$exp) + as.factor(tn$xpp$ampInd)*tn$xpp$trialType + I( log(tn$xpp$ISI)) )
  lm1           <- lm(tn$xpp[[DVstr]] ~ as.factor(tn$xpp$exp) + as.factor(tn$xpp$ampInd)*tn$xpp$trialType + I( log(tn$xpp$ISI)) + tn$xpp$deltaAmpInd )
  lmresult <- anova(lm0,lm1)
  
  #if(norm){
  #  tn$xpp$DVn    <- lm0$residuals
  #}
  
  tn$xpp$DVn <- tn$xpp[[DVstr]]
  if(norm){
    
    ## ============ small model
    testLambda <- function(lmd){  
      print(lmd)
      tn$xpp$oKer  <-        exp(  lmd*tn$xpp$ISI )
      lmExp   <- lm(tn$xpp$DV ~ 1 +  tn$xpp$oKer  )
      return( sum(lmExp$residuals^2) )
    }
    
    resExp     <- optimize(f=testLambda, lower=-10, upper=-0.1, maximum = FALSE)
    lmExp      <- lm(tn$xpp$DV ~ 1  +  I( exp(resExp$minimum * tn$xpp$ISI)  ))
    
    mxAmp <- lmExp$coefficients[1]
    tn$xpp$DVn <- tn$xpp$DV/mxAmp
  }
  
  
  N1        <- tapply(tn$xpp$DVn, list(tn$xpp$dB,tn$xpp$deltaAmpInd), mean )
  dbax      <- tapply(tn$xpp$dB,          tn$xpp$ampInd,      mean)
  deltax    <- tapply(tn$xpp$deltaAmpInd, tn$xpp$deltaAmpInd, mean)
  
  
  if (show){
    thisDir <- paste( pp(),tn$exp,'/FLA/',DVstr,'/',sep='')
    if(!dir.exists(thisDir))
      dir.create(thisDir, recursive=T)    
    
    fileName <- paste(tn$exp,'/FLA/',DVstr,'/',baseName,tn$exp,'.eps',sep='')
    if(eps)
      eps.file(fileName, pdf=T, s=3)
    
    zlim      <- c(-1,1) * na.max(abs(N1)) 
    image(dbax, deltax, N1, col=bluered.colors(100), zlim=zlim)
    
    
    if(eps)
      dev.off()
  }
  
  res <- list( N1=N1, xax=dbax, deltax=deltax, lmresult=lmresult  ) 
  return(res)
}