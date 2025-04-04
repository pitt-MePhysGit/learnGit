## ============================================
## master session-average VProbe
##
rm(list=ls())
defpar <- par()

baseDir <- '~/Dropbox/'
source(paste(baseDir, 'RSkernel/r-files/functions/source_all_functions.R',sep=''))
source(paste(baseDir, 'ampOddClick/r-files/source_all_functions.R',sep=''))

## preliminaries
.GlobalEnv$pp        <- function(){paste(baseDir, '/ampOddClick/word/VProbe/plots/',sep='')}
.GlobalEnv$rda       <- function(){'~/Dropbox/ampOddClick/Rdata/'}
.GlobalEnv$trda      <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
.GlobalEnv$eeglabdir <- function(){'/Volumes/Drobo5D3/EEG/EEGLab/ampOddClick/'}

source( paste(baseDir,'ampOddClick/r-files/functions/prepDFVprobe.R',sep='') )
source( paste(baseDir,'ampOddClick/r-files/functions/prep_df_mua.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/functions_sessionAverageVProbe.R',sep='') )

colClasses <- c('character',rep('numeric',2),'character',rep('numeric',5)  ,rep('character',3),rep('numeric',8) )
df         <- read.table(paste(baseDir,'ampOddClick/ampOddClickVProbe_link.txt',sep=''),header=T,colClasses=colClasses)
df         <- prepDFVprobe(df)
df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))


cleanOnly   <- FALSE
#badDays     <- c(12, 14, 16, 18, 22, 24)


df$use      <- c(1)
if (cleanOnly)
  df$use[ which(df$toneArtefact+df$onsetArtefact>0) ] <- 0
#df$use[badDays] <- 0

dataSet <- 'g'
animal <- c('Walter')
#animal <- c('Walter','Sam')

useInd     <- which( df$dataSet==dataSet & df$animal%in%animal & df$drug=='N' & df$use)# & df$probeString=='d')
#useInd     <- useInd[1:13]
#useInd <- useInd[-25]              # something wrong with Sam_20181206_1215
df[useInd,]

## ======================================== set folder for plots
animalStr <- '_all'
if (length(animal)==1)
  animalStr <- paste('_',animal,sep='')

cleanStr  <- ''
if (cleanOnly)
  cleanStr <- '_clean'

uIdir <- paste('g', animalStr,cleanStr,sep='')
.GlobalEnv$pp   <- function(){paste(baseDir, '/ampOddClick/word/VProbe/plots/',uIdir,'/',sep='')}
if (!file.exists(pp()))
  system( paste('mkdir', pp()) )

## ===================================
muascl <- 1   # 2
lfpscl <- 70  # 13
csdscl <- .4   # 1.5

ylim=c(-1500,1000)

taxis <- seq(-150,749, by=1)
picFile <- paste('GA.eps',sep='')
eps.file(picFile, pdf=T, s=5, ratio=1/2)
disc <- showAll(useInd, extn='mean',xnd=1,minN=6, xlim=c(-50,250), ylim=ylim,stp=50,dskip=4,taxis=taxis,muascl=muascl,lfpscl=lfpscl,csdscl=csdscl)
dev.off()


picFile <- paste('GA_zoom.eps',sep='') 
eps.file(picFile, pdf=T, s=5, ratio=1/2)
disc <- showAll(useInd, extn='mean',xnd=1,minN=6, xlim=c(-10,50), ylim=ylim,stp=50,dskip=2,taxis=taxis,muascl=muascl,lfpscl=lfpscl,csdscl=csdscl)
dev.off()



##
drng <- c(-2000,2000)
dndx <- which(disc$daxis >= drng[1] & disc$daxis <= drng[2])
plot(NA, xlim=c(-50,150), ylim=c(0,1+dim(disc$resMatM)[3]) )
for (dx in 1:dim(disc$resMatM)[3] ){
  points(disc$taxis, dx + apply(disc$resMatM[,dndx,dx],1,na.mean)/max(abs(apply(disc$mnMUA[,dndx],1,na.mean))), type='l' )  
}
abline(v=0)
df[useInd[c(6,8)],]
df[useInd[c(41)],]
##

## =====================================================================================================================
## =================================== load data for logOct
ISIlist <- list()
for (i in 1:6){
  #picFile <- paste('GA_ISIbinOct',i,'.eps',sep='')
  #eps.file(picFile, pdf=T, s=5, ratio=1/2)
  ISIlist[[i]] <- showAll(useInd, extn='mean_ISIbinOct',xnd=i,minN=4, xlim=c(-50,250), ylim=c(-2000,2000),stp=50,taxis=taxis)
  #dev.off()
}

## =================================== load data for ampInd
AMPlist <- list()
for (i in 1:5){
  #picFile <- paste('GA_ISIbinOct',i,'.eps',sep='')
  #eps.file(picFile, pdf=T, s=5, ratio=1/2)
  AMPlist[[i]] <- showAll(useInd, extn='mean_Amp',xnd=i,minN=4, xlim=c(-50,250), ylim=c(-2000,2000),stp=50,taxis=taxis)
  #dev.off()
}


## =============================================================================
## ======================================== grand average across specified depth
runPlots <- function(whatStr,drng=c(-2000,500)){
  
  dndx <- which(ISIlist[[1]]$daxis >= drng[1] & ISIlist[[1]]$daxis <= drng[2])

  ylim <- c(-0.05, .55)
  if (animal=='Walter')
    ylim <- c(-0.15, .95)
  
  csgn <- +1
  if (identical(whatStr,'mnCSD'))
    csgn <- -1
  
  
  if (identical(whatStr,'mnMUA')){
    ylim <- c(-.15, 0.75)
    dylim <- c(-.15, 0.15)
  }
  
  if (identical(whatStr,'mnCSD')){
    ylim <- c(-0.4,0.1)
    dylim <- c(-0.1,0.1)
  }
  
  if (identical(whatStr,'mnLFP')){
    ylim  <- c(-40,40)
    dylim <- c(-10,10)
  }
  
  thisDir   <- paste(pp(),whatStr,sep='')
  if (!file.exists(thisDir))
    system( paste('mkdir', thisDir ))
  
  
  ## ========================================= grand averages
  col <- soa.colors(6)
  picFile <- paste(whatStr,'/GA_SOA_',whatStr,'.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  plot( ISIlist[[1]]$taxis, apply(ISIlist[[1]][[whatStr]][,dndx],1,na.mean), col=col[1], type='l', lwd=2, ylim=ylim, xlim=c(-10,250) )
  for (i in 2:5)
    points( ISIlist[[i]]$taxis, apply(ISIlist[[i]][[whatStr]][,dndx],1,na.mean), col=col[i], lwd=2, type='l' )
  abline(h=0)
  dev.off()
  
  col <- bluered.colors(5)
  picFile <- paste(whatStr,'/GA_AMP_',whatStr,'.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  plot( AMPlist[[1]]$taxis, apply(AMPlist[[1]][[whatStr]][,dndx],1,na.mean), col=col[1], type='l', lwd=2, ylim=ylim, xlim=c(-10,250) )
  for (i in 2:5)
    points( AMPlist[[i]]$taxis, apply(AMPlist[[i]][[whatStr]][,dndx],1,na.mean), col=col[i], lwd=2, type='l' )
  abline(h=0)
  dev.off()
  
  ## ============================================ zoom
  col <- soa.colors(6)
  picFile <- paste(whatStr,'/GA_SOA_',whatStr,'_zoom.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  plot( ISIlist[[1]]$taxis, apply(ISIlist[[1]][[whatStr]][,dndx],1,na.mean), col=col[1], type='l', lwd=2, ylim=ylim, xlim=c(-10,100) )
  for (i in 2:5)
    points( ISIlist[[i]]$taxis, apply(ISIlist[[i]][[whatStr]][,dndx],1,na.mean), col=col[i], lwd=2, type='l' )
  abline(h=0)
  dev.off()
  
  
  col <- bluered.colors(5)
  picFile <- paste(whatStr,'/GA_AMP_',whatStr,'_zoom.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  plot( AMPlist[[1]]$taxis, apply(AMPlist[[1]][[whatStr]][,dndx],1,na.mean), col=col[1], type='l', lwd=2, ylim=ylim, xlim=c(-10,100) )
  for (i in 2:5)
    points( AMPlist[[i]]$taxis, apply(AMPlist[[i]][[whatStr]][,dndx],1,na.mean), col=col[i], lwd=2, type='l' )
  abline(h=0)
  dev.off()
  
  
  ## ============================================ log-time
  taxis                           <- AMPlist[[1]]$taxis
  toff                            <- 0
  logtaxis                        <- log2( taxis + toff)
  logtaxis[ (taxis + toff <= 1)]  <- .01*(taxis[ which(taxis + toff <= 1) ] - 1 + toff) #- log2(taxis[min(which(taxis>-toff))] )# + .01*taxis[ (llp$taxis + toff == 0)]
  logtaxis <- as.vector(logtaxis)
  
  trng <- c(3,550)
  xlim <- range( logtaxis[ taxis>=trng[1] & taxis<trng[2]] )
  
  col <- soa.colors(6)
  picFile <- paste(whatStr,'/GA_SOA_',whatStr,'_log.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  plot( logtaxis, apply(ISIlist[[1]][[whatStr]][,dndx],1,na.mean), col=col[1], type='l', lwd=2, ylim=ylim, xlim=xlim,xaxt='n' )
  for (i in 2:5)
    points( logtaxis, apply(ISIlist[[i]][[whatStr]][,dndx],1,na.mean), col=col[i], lwd=2, type='l' )
  abline(h=0)
  taxat <- c(10,20,30,40,50,100,200,300,400)
  taxnd <- which( round(taxis)%in%c(taxat) )
  axis(1, at=logtaxis[taxnd], taxat  )
  dev.off()
  
  col <- bluered.colors(5)
  picFile <- paste(whatStr,'/GA_AMP_',whatStr,'_log.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  plot( logtaxis, apply(AMPlist[[1]][[whatStr]][,dndx],1,na.mean), col=col[1], type='l', lwd=2, ylim=ylim, xlim=xlim, xaxt='n' )
  for (i in 2:5)
    points( logtaxis, apply(AMPlist[[i]][[whatStr]][,dndx],1,na.mean), col=col[i], lwd=2, type='l' )
  abline(h=0)
  taxat <- c(10,20,30,40,50,100,200,300,400)
  taxnd <- which( round(taxis)%in%c(taxat) )
  axis(1, at=logtaxis[taxnd], taxat  )
  dev.off()
  
  
  ## ============================== Grand average difference between extreme conditions
  ylim <- c(-.25,.25)
  if (animal=='Walter')
    ylim <- c(-0.25, .25)
  
  # ============================================================== normal time-axis
  col <- soa.colors(6)
  #drng <- c(-1500,0)
  picFile <- paste(whatStr,'/GA_Delta_SOA_',whatStr,'.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  plot( ISIlist[[1]]$taxis, apply(ISIlist[[1]][[whatStr]][,dndx],1,na.mean)-apply(ISIlist[[3]][[whatStr]][,dndx],1,na.mean), col=col[1], type='l', lwd=2, ylim=dylim, xlim=c(-10,150) )
  for (i in 2:5)
    points( ISIlist[[i]]$taxis, apply(ISIlist[[i]][[whatStr]][,dndx],1,na.mean)-apply(ISIlist[[3]][[whatStr]][,dndx],1,na.mean), col=col[i], lwd=2, type='l' )
  abline(v=0)
  abline(v=c(8.5,14.5,22),lty=2,lwd=.5)
  dev.off()
  
  
  col <- bluered.colors(5)
  picFile <- paste(whatStr,'/GA_Delta_AMP_',whatStr,'.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  plot( AMPlist[[1]]$taxis, apply(AMPlist[[1]][[whatStr]][,dndx],1,na.mean)-apply(AMPlist[[3]][[whatStr]][,dndx],1,na.mean), col=col[1], type='l', lwd=2, ylim=dylim, xlim=c(-10,150) )
  for (i in 2:5)
    points( AMPlist[[i]]$taxis, apply(AMPlist[[i]][[whatStr]][,dndx],1,na.mean)-apply(AMPlist[[3]][[whatStr]][,dndx],1,na.mean), col=col[i], lwd=2, type='l' )
  abline(v=0)
  abline(v=c(8.5,14.5,22),lty=2,lwd=.5)
  dev.off()
  
  ## ======================================================= zoomed version
  col <- soa.colors(6)
  picFile <- paste(whatStr,'/GA_Delta_SOA_',whatStr,'_zoom.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  plot( ISIlist[[1]]$taxis, apply(ISIlist[[1]][[whatStr]][,dndx],1,na.mean)-apply(ISIlist[[3]][[whatStr]][,dndx],1,na.mean), col=col[1], type='l', lwd=2, ylim=dylim, xlim=c(-10,50) )
  for (i in 2:5)
    points( ISIlist[[i]]$taxis, apply(ISIlist[[i]][[whatStr]][,dndx],1,na.mean)-apply(ISIlist[[3]][[whatStr]][,dndx],1,na.mean), col=col[i], lwd=2, type='l' )
  abline(v=0)
  abline(v=c(8.5,14.5,22),lty=2,lwd=.5)
  dev.off()
  
  
  col <- bluered.colors(5)
  picFile <- paste(whatStr,'/GA_Delta_AMP_',whatStr,'_zoom.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  plot( AMPlist[[1]]$taxis, apply(AMPlist[[1]][[whatStr]][,dndx],1,na.mean)-apply(AMPlist[[3]][[whatStr]][,dndx],1,na.mean), col=col[1], type='l', lwd=2, ylim=dylim, xlim=c(-10,50) )
  for (i in 2:5)
    points( AMPlist[[i]]$taxis, apply(AMPlist[[i]][[whatStr]][,dndx],1,na.mean)-apply(AMPlist[[3]][[whatStr]][,dndx],1,na.mean), col=col[i], lwd=2, type='l' )
  abline(v=0)
  abline(v=c(8.5,14.5,22),lty=2,lwd=.5)
  dev.off()
  
  ## ========================================================== log time-axis
  col <- soa.colors(6)
  picFile <- paste(whatStr,'/GA_Delta_SOA_',whatStr,'_log.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  plot( logtaxis, apply(ISIlist[[1]][[whatStr]][,dndx],1,na.mean)-apply(ISIlist[[3]][[whatStr]][,dndx],1,na.mean), col=col[1], type='l', lwd=2, ylim=dylim, xlim=xlim, xaxt='n' )
  for (i in 2:5)
    points( logtaxis, apply(ISIlist[[i]][[whatStr]][,dndx],1,na.mean)-apply(ISIlist[[3]][[whatStr]][,dndx],1,na.mean), col=col[i], lwd=2, type='l' )
  abline(v=0)
  abline(v=log2(c(8.5,14.5,22)),lty=2,lwd=.5)
  taxat <- c(10,20,30,40,50,100,200,300,400)
  taxnd <- which( round(taxis)%in%c(taxat) )
  axis(1, at=logtaxis[taxnd], taxat  )
  dev.off()
  
  
  col <- bluered.colors(5)
  picFile <- paste(whatStr,'/GA_Delta_AMP_',whatStr,'_log.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  plot( logtaxis, apply(AMPlist[[1]][[whatStr]][,dndx],1,na.mean)-apply(AMPlist[[3]][[whatStr]][,dndx],1,na.mean), col=col[1], type='l', lwd=2, ylim=dylim, xlim=xlim, xaxt='n' )
  for (i in 2:5)
    points( logtaxis, apply(AMPlist[[i]][[whatStr]][,dndx],1,na.mean)-apply(AMPlist[[3]][[whatStr]][,dndx],1,na.mean), col=col[i], lwd=2, type='l' )
  abline(v=0)
  abline(v=log2(c(8.5,14.5,22)),lty=2,lwd=.5)
  taxat <- c(10,20,30,40,50,100,200,300,400)
  taxnd <- which( round(taxis)%in%c(taxat) )
  axis(1, at=logtaxis[taxnd], taxat  )
  dev.off()
  

  ## ============================================= double-delta
  #ylim=c(-.05,.3)
  #if (animal=='Walter')
  #  ylim <- c(-0.15, .55)
  
  picFile <- paste(whatStr,'/GA_Delta_AMP_Delta_SOA_',whatStr,'.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  plot(   AMPlist[[1]]$taxis, csgn*(apply(AMPlist[[5]][[whatStr]][,dndx],1,na.mean)-apply(AMPlist[[1]][[whatStr]][,dndx],1,na.mean)), col='blue', type='l', lwd=2, ylim=c(0.5,2)*dylim, xlim=c(-10,250) )
  points( AMPlist[[1]]$taxis, csgn*(apply(ISIlist[[5]][[whatStr]][,dndx],1,na.mean)-apply(ISIlist[[1]][[whatStr]][,dndx],1,na.mean)), col='red', lwd=2, type='l' )
  abline(v=0)
  abline(h=0)
  dev.off()
  
  
  picFile <- paste(whatStr,'/GA_Delta_AMP_Delta_SOA_',whatStr,'_zoom.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  plot(   AMPlist[[1]]$taxis, csgn*(apply(AMPlist[[5]][[whatStr]][,dndx],1,na.mean)-apply(AMPlist[[1]][[whatStr]][,dndx],1,na.mean)), col='blue', type='l', lwd=2, ylim=c(0.5,2)*dylim, xlim=c(-10,50) )
  points( AMPlist[[1]]$taxis, csgn*(apply(ISIlist[[5]][[whatStr]][,dndx],1,na.mean)-apply(ISIlist[[1]][[whatStr]][,dndx],1,na.mean)), col='red', lwd=2, type='l' )
  abline(v=0)
  abline(h=0)
  #abline(v=c(9,16,22),lty=2,lwd=.5, col='blue')
  #abline(v=c(16,24),lty=2,lwd=.5, col='red')
  dev.off()
}



runPlots('mnCSD', drng=c( -500,   0))
runPlots('mnMUA', drng=c(-2000,1000))
runPlots('mnLFP', drng=c(  250,1000))

## =============================================================================
## ========================================================= modulation by depth
runDepthPlots <- function(whatStr){
  
  csgn <- +1
  if (identical(whatStr,'mnCSD'))
    csgn <- -1
  
  if (identical(whatStr,'mnMUA'))
    dz_lim <- 0.6
  
  if (identical(whatStr,'mnCSD')){
    dz_lim <- 0.4
    if (animal=='Walter')
      dz_lim <- .75
  }
  if (identical(whatStr,'mnLFP'))
    dz_lim <- 40
  
  #ylim=c(-2000,1500)
  #if (animal=='Walter')
  ylim=c(-1500,1000)
  
  mnAMP1 <- apply(AMPlist[[1]][[whatStr]],c(1,2),na.mean)
  mnAMP5 <- apply(AMPlist[[5]][[whatStr]],c(1,2),na.mean)
  
  mnISI1 <- apply(ISIlist[[1]][[whatStr]],c(1,2),na.mean)
  mnISI5 <- apply(ISIlist[[5]][[whatStr]],c(1,2),na.mean)
  
  picFile <- paste(whatStr,'/GA_Delta_AMP_Image_',whatStr,'.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  image(AMPlist[[1]]$taxis, AMPlist[[1]]$daxis,  csgn*(mnAMP5-mnAMP1), col=jet.colors(200), zlim=c(-1,1)*dz_lim, useRaster=T,xlim=c(-10,150), ylim=ylim)
  abline(h=0);abline(v=0)
  dev.off()
  
  picFile <- paste(whatStr,'/GA_Delta_SOA_Image_',whatStr,'.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  image(AMPlist[[1]]$taxis, AMPlist[[1]]$daxis,  csgn*(mnISI5-mnISI1), col=jet.colors(200), zlim=c(-1,1)*dz_lim, useRaster=T,xlim=c(-10,150), ylim=ylim)
  abline(h=0);abline(v=0)
  abline(h=0);abline(v=0)
  dev.off()
  
  picFile <- paste(whatStr,'/GA_Delta_AMP_Image_',whatStr,'_zoom.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  image(AMPlist[[1]]$taxis, AMPlist[[1]]$daxis,  csgn*(mnAMP5-mnAMP1), col=jet.colors(200), zlim=c(-1,1)*dz_lim, useRaster=T,xlim=c(-10,50), ylim=ylim)
  abline(h=0);abline(v=0)
  dev.off()
  
  picFile <- paste(whatStr,'/GA_Delta_SOA_Image_',whatStr,'_zoom.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  image(AMPlist[[1]]$taxis, AMPlist[[1]]$daxis,  csgn*(mnISI5-mnISI1), col=jet.colors(200), zlim=c(-1,1)*dz_lim, useRaster=T,xlim=c(-10,50), ylim=c(-2000,1500))
  abline(h=0);abline(v=0)
  abline(h=0);abline(v=0)
  dev.off()
  
}

runDepthPlots('mnCSD')
runDepthPlots('mnMUA')
runDepthPlots('mnLFP')




#image(AMPlist[[1]]$taxis, AMPlist[[1]]$daxis,  log2( (mnAMP5%>=%.00001)/(mnAMP1%>=%.00001))%<=%3, col=jet.colors(200), zlim=c(-3,3),useRaster=T,xlim=c(5,50))
#image(AMPlist[[1]]$taxis, AMPlist[[1]]$daxis,  log2( (mnISI5%>=%.00001)/(mnISI1%>=%.00001))%<=%3, col=jet.colors(200), zlim=c(-3,3),useRaster=T,xlim=c(5,50))

#image(AMPlist[[1]]$taxis, AMPlist[[1]]$daxis,  mnAMP5-mnAMP1, col=jet.colors(200), zlim=c(-1,1)*dz_lim, useRaster=T,xlim=c(5,50))
#image(AMPlist[[1]]$taxis, AMPlist[[1]]$daxis,  mnISI5-mnISI1, col=jet.colors(200), zlim=c(-1,1)*dz_lim, useRaster=T,xlim=c(5,50))


