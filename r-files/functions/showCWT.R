showCWT <- function(cw,flab=c(10,20,30,50,100,200),xlim=c(-50,350),ylim=c(3,250),zlim=c(0,.1),yoff=10, what='cw', logtime=F,freqList=NULL,trngList=trngList,...){
  
  yConv <- function(x){ log(2*pi*x/(cw$w0*cw$nF))/log(2)   }
  cwdata <- cw[[what]]
  
  if (!logtime){
    image(cw$taxis, cw$faxlog, id(Mod(cwdata[,dim(cwdata)[2]:1])),xlim=xlim,ylim=yConv(ylim),zlim=zlim,col=jet.colors(256),useRaster=F,yaxt='n' )
    #image(cw$taxis, cw$faxlog,id(Mod(cw$cw[,dim(cw$cw)[2]:1])),xlim=xlim,ylim=yConv(ylim),zlim=zlim,col=jet.colors(256),useRaster=F,yaxt='n' )
    axis(2,at=yConv(flab),  label=flab)
    abline(h=yConv(flab), col='white',lty=2 )
    abline(v=0,col='white')
    points(cw$taxis,1*cw$mn/max(abs(cw$mn))+yConv(yoff),type='l')
  }
  
  if (logtime){
    tlim <-  cw$logtaxis[ which(cw$taxis%in%xlim) ]
    flab <- sort( unique(unlist(freqList)) )
    taxat <- c(5,10,15,20,30,40,50,100,200,300,500)
    taxat <- sort( unique(unlist(trngList)) )
    
    image(cw$logtaxis, cw$faxlog, id(Mod(cwdata[,dim(cwdata)[2]:1])),xlim=tlim,ylim=yConv(ylim),zlim=zlim,col=jet.colors(256),useRaster=F,yaxt='n',xaxt='n' )
    #image(cw$taxis, cw$faxlog,id(Mod(cw$cw[,dim(cw$cw)[2]:1])),xlim=xlim,ylim=yConv(ylim),zlim=zlim,col=jet.colors(256),useRaster=F,yaxt='n' )
    
    taxnd <- which( cw$taxis%in%c(taxat) )
    axis(1, at=cw$logtaxis[taxnd], taxat  )
    abline(v=cw$logtaxis[taxnd],lty=2,lwd=.5, col='gray')
    
    axis(2,at=yConv(flab),  label=flab)
    abline(h=yConv(flab), col='gray',lty=2 )
    points(cw$logtaxis,1*cw$mn/max(abs(cw$mn))+yConv(yoff),type='l')
  }
}
