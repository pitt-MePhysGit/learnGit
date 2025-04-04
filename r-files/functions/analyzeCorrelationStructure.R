analyzeCorrelationStructure <- function(predMat1, predMat2,fileName='CC.eps',zlimCor=.1, alpha=0.05){
  
  getPmat <- function(predMat1,predMat2){
    N1      <- dim(predMat1)[2]
    N2      <- dim(predMat2)[2]
    corPmat <- array(1, c(N1,N2))
    dimnames(corPmat)[[1]] <- dimnames(predMat1)[[2]]
    dimnames(corPmat)[[2]] <- dimnames(predMat2)[[2]]
    
    for (sx1 in 1:N1){
      for (sx2 in 1:N2){
        corPmat[sx1           ,N2+1-sx2] <- cor.test(predMat1[,sx1],predMat2[,sx2])$p.value
      }
    }
    return(corPmat)
  }
  
  
  N1      <- dim(predMat1)[2]
  N2      <- dim(predMat2)[2]
  
  Ntst  <- (N1-1) *  (N2-1)
  
  alphaB <- alpha/Ntst
  
  str1 <- dimnames(predMat1)[[2]]
  str2 <- dimnames(predMat2)[[2]]
  
  chNr1 <- sapply( str1, function(st){ as.numeric( strsplit(st,'_')[[1]][3]) } )
  chNr2 <- sapply( str2, function(st){ as.numeric( strsplit(st,'_')[[1]][3]) } )
  
  Nchan1 <- length(unique(chNr1))
  Nchan2 <- length(unique(chNr2))
  
  tStr1 <- sapply( str1, function(st){ strsplit(st,'_')[[1]][1] } )
  tStr2 <- sapply( str2, function(st){ strsplit(st,'_')[[1]][1] } )
  
  pStr1 <- sapply( str1, function(st){ strsplit(st,'_')[[1]][2] } )
  pStr2 <- sapply( str2, function(st){ strsplit(st,'_')[[1]][2] } )
  
  cmpStrName1 <- unique( paste(tStr1, pStr1, sep='_'))
  cmpStrName2 <- unique( paste(tStr2, pStr2, sep='_'))
  
  Nstr1 <- length(unique(tStr1))
  Nstr2 <- length(unique(tStr2))
  
  corMat   <- cor(predMat1, predMat2)
  corSqMat <- sign(corMat) * corMat^2
  corPmat  <- getPmat(predMat1, predMat2) 
  corPmatbf <- corPmat
  corPmatbf[ corPmat>alphaB ] <- 1
  
  eps.file(fileName,pdf=T,s=2,mf=c(2,1),ratio=1)
  par( mfrow=c(1,2), mar=c(1,3,3,1) )
  image( 1:N1, 1:N2, -log10(corPmat)%<=%35, zlim=c(-log10(alphaB),50), col=heat.colors(100), yaxt='n', xaxt='n' )
  abline( h=seq(Nchan2+.5,by=Nchan2,to=N2))#, lwd=c(1,2,1,1,2,1,1,2,1,1,2) )
  abline( v=seq(Nchan1+.5,by=Nchan1,to=N1))#, lwd=c(1,2,1,1,2,1,1,2,1,1,2) )
  axis(3, at=seq(Nchan1/2, by=Nchan1, length=length(cmpStrName1)), cmpStrName1,cex.axis=.5  )
  axis(2, at=seq(Nchan2/2, by=Nchan2, length=length(cmpStrName2)), cmpStrName2[seq(length(cmpStrName2),to=1,by=-1)],cex.axis=.5  )
  if(identical(predMat1, predMat2))
    abline(a=N1+1,b= -1)
  
  image( 1:N1, 1:N2, (corSqMat[,seq(dim(corSqMat)[2],to=1,by=-1)]%<=%zlimCor)%>=%-zlimCor, zlim=c(-zlimCor,zlimCor), col=bluered.colors(100), yaxt='n', xaxt='n'  )
  abline( h=seq(Nchan2+.5,by=Nchan2,to=N2))#, lwd=c(1,2,1,1,2,1,1,2,1,1,2) )
  abline( v=seq(Nchan1+.5,by=Nchan1,to=N1))#, lwd=c(1,2,1,1,2,1,1,2,1,1,2) )
  axis(3, at=seq(Nchan1/2, by=Nchan1, length=length(cmpStrName1)), cmpStrName1,cex.axis=.5  )
  axis(2, at=seq(Nchan2/2, by=Nchan2, length=length(cmpStrName2)), cmpStrName2[seq(length(cmpStrName2),to=1,by=-1)],cex.axis=.5  )
  if(identical(predMat1, predMat2))
    abline(a=N1[1]+1,b= -1)
  dev.off()
  
  return(list(corSqMat=corSqMat, corMat=corMat, corPmat=corPmat, corPmatbf=corPmatbf) ) 
}
