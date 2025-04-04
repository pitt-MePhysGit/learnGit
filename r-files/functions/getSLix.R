getSLix <- function(what='bsp',when='thal',xtn='', freqList=NULL, trngList=NULL){
  
  trng    <- trngList[[ when ]]
  sL      <- data.frame( what=what,   icx=1, cmpStrName=paste(what,when,'1',sep='_'),   tmin=trng[1],   tmax=trng[2], stringsAsFactors=FALSE)
  sL[2,]  <- data.frame( what=what,   icx=2, cmpStrName=paste(what,when,'2',sep='_'),   tmin=trng[1],   tmax=trng[2], stringsAsFactors=FALSE)
  sL[3,]  <- data.frame( what=what,   icx=3, cmpStrName=paste(what,when,'3',sep='_'),   tmin=trng[1],   tmax=trng[2], stringsAsFactors=FALSE)
  sL[4,]  <- data.frame( what=what,   icx=4, cmpStrName=paste(what,when,'4',sep='_'),   tmin=trng[1],   tmax=trng[2], stringsAsFactors=FALSE)
  sL[5,]  <- data.frame( what=what,   icx=5, cmpStrName=paste(what,when,'5',sep='_'),   tmin=trng[1],   tmax=trng[2], stringsAsFactors=FALSE)
  sL[6,]  <- data.frame( what=what,   icx=6, cmpStrName=paste(what,when,'6',sep='_'),   tmin=trng[1],   tmax=trng[2], stringsAsFactors=FALSE)
  return(sL)
}