subjFigure_CmpByAmpAndSOA <- function(tn,DVstr='n1.r',norm=T,show=T,eps=T,trialType=c(1),baseName='cmpByAmpAndSOA'){
  
  tn        <- subs.data(tn,   tn$xpp$trialType%in%c(trialType) )
  
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
  
  ## AMPLITUDE by SOA and intensity
  mn     <- tapply(tn$xpp$DVn, list(tn$xpp$dB,tn$xpp$isiOct), na.mean )
  mnse  <- tapply(tn$xpp$DVn, list(tn$xpp$dB,tn$xpp$isiOct), na.se )
  dbax   <- tapply(tn$xpp$dB, list(tn$xpp$dB), mean )
  isiax  <- tapply(tn$xpp$ISI, tn$xpp$isiOct, mean)
  vi <- 1:5
  
  if (show){
    thisDir <- paste( pp(),tn$exp,'/FLA/',DVstr,'/',sep='')
    if(!dir.exists(thisDir))
      dir.create(thisDir, recursive=T)    
    
    fileName <- paste(tn$exp,'/FLA/',DVstr,'/',baseName,tn$exp,'.eps',sep='')
    if(eps)
      eps.file(fileName, pdf=T, s=3)
    
    par(mfcol=c(1,2))
    col <- soa.colors(6)
    errorplot(dbax,mn[,vi[1]], ey=mnse[,vi[1]], col='black', ylim=c(0,2),pch=20)
    for (sx in 1:length(vi))
      errorplot(dbax,mn[ ,vi[sx]], ey=mnse[,vi[sx]], col=col[sx], add=T,pch=20)
    
    
    col <- bluered.colors(5)
    errorplot(log2(isiax[vi]),mn[1,vi], ey=mnse[1,vi], col='black', ylim=c(0,2),pch=20)
    for (ax in 1:length(dbax))
      errorplot(log2(isiax[vi]),mn[ax,vi], ey=mnse[ax,vi], col=col[ax], add=T,pch=20)
      
    
    if(eps)
      dev.off()
  }
  
  res <- list( N1=mn, xax=dbax, yax=isiax ) 
  return(res)
}