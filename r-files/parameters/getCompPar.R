getComponentParameters <- function(cmpStr, animal='default', show=F, add=T, tlim=NULL, ylim=NULL,tFUN=id,yFUN=id, ... ){

  fileStr <- paste(baseDir, '/ampOdd_click/r-files/parameters/compPar_',animal,'.txt',sep='' )
  
  tbl     <- read.table(fileStr,header=T)
  ind     <- which(tbl$name %in% cmpStr)
  
  if(length(ind)==0)
    print(cmpStr)
  
  Nchan   <- dim(tbl)[2] - 4
  chInd   <- which(tbl[ind,5:(4+Nchan)]==1)
  
  cmpList <-  list(trng=c( tbl$from[ind],tbl$to[ind]), chInd=chInd, signalStr=tbl[ind]$sgnl, cmpStr=cmpStr, animal=animal)
  trng    <- cmpList$trng
  
  
  ## ============================================
  if(show){
    if( !add & is.null(tlim) )
      tlim <- ifelse( 'llp'%in%signal | 'raw'%in%signal, c(-50,250), c(-10,50)  ) 
    
    if( !add & is.null(ylim) )
      ylim <- ifelse( 'llp'%in%signal | 'raw'%in%signal, c(-30,30), c(-10,10)  ) 
    
    if(!add)
      plot(c(-1000),c(-1000), ylim=ylim, xlim=tlim)
    
    #polygon( tFUN( c(trng,trng[c(2,1)]) ), yFUN( rep(ylim,each=2) ), col=adjustcolor("white", alpha=0) )
    
    
    negCmp <- substr(cmpStr,1,1)=='n'
    if(nchar(cmpStr)>2 & substr(cmpStr,3,3)=='T')
      negCmp <- T
    
    
    if( negCmp ){
      ylim[2]<-mean(ylim)
      #tcol <- 'blue'
      #tcol <- 'darkgray
    }
    
    if( !negCmp ){
      ylim[1]<-mean(ylim)
      #tcol <- 'red'
      #tcol <- 'darkgray'
    }
    tcol <- componentColor2( cmpStr )
    
    polygon( tFUN( c(trng,trng[c(2,1)])), yFUN( rep(ylim,each=2)),    col=adjustcolor(tcol, alpha=0.25), ... )
  }
  
  ## ===========================================  
  return(cmpList)

}