showICx <- function(llp,chx, xlim=c(-50,350),ylim=c(1,300),subs=NULL,freqList=NULL,trngList=NULL,...){
  
  if (is.null(subs))
    subs <- rep(T, dim(llp$lfp)[1])
  
  tmp          <- apply(llp$lfp,c(1,3),max) - apply(llp$lfp,c(1,3),min)
  llp$xpp$p2p  <- apply( tmp, 1, max )
  mxp2p        <- median(llp$xpp$p2p) + 4 * diff(quantile(llp$xpp$p2p,probs=c(.25,.75))) 

  cw1 <- getCWTevk(llp, which(llp$xpp$p2p<mxp2p & subs),chx=chx )
  
  flab <- c(1,2,3,10,20,30,50,100,200)
  
  zlim <- c(0,max(Mod(cw1$cw)))
  yoff <- 10
  showCWT(cw1, xlim=xlim, ylim=ylim, zlim=zlim,yoff=yoff,flab=flab,freqList=freqList,trngList=trngList,...)
}