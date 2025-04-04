## master channel average BSP
##

## ============================================
## master channel-average VProbe
##  in contrast to master_sessionAverage_VProbe, this script derives channel-based measures and error bars 
# this copied and updated. Previous code has been saved in old_master_channelAverage_VProbe.R

startFromScratch <- FALSE
# Preliminaries     ====
rm(list=ls())
efpar <- par()

baseDir <- '~/Dropbox/'
source('~/Dropbox/r-compile/signalProcessoR/R/signalProcessoR.R')
source('~/Dropbox/r-compile/rUtils/R/rUtils.R')


.GlobalEnv$pp        <- function(){paste(baseDir, '/ampOddClick/word/EEG/plots/',sep='')}

source( paste(baseDir,'ampOddClick/r-files/functions/prepDFVprobe.R',sep='') )
source( paste(baseDir,'ampOddClick/r-files/functions/prep_df_mua.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/functions_sessionAverageVProbe.R',sep='') )

source( paste(baseDir,'ampOddClick/r-files/figures/call_AOC_timeseries_bySOA.R',sep='') )
source( paste(baseDir,'ampOddClick/r-files/figures/figure_AOC_timeseries_bySOA.R',sep='') )

source( paste(baseDir,'ampOddClick/r-files/figures/call_AOC_timeseries_byAMP.R',sep='') )
source( paste(baseDir,'ampOddClick/r-files/figures/figure_AOC_timeseries_byAMP.R',sep='') )

#source( paste(baseDir,'ampOddClick/r-files/functions/loadGASession.R',sep='') )

source( paste(baseDir,'RSkernel/r-files/functions/expand2chdf.R',sep='') )


## LOAD data frame      ====
project               <- 'ampOddClick'
linkFileName          <- 'ampOddclick_link.txt' # mostly sc with some old eeg only
linkFileNameEEG       <- 'ampOddclickEEG_link.txt' # all eeg


df                 <- prepDFsignalProcessor(project,linkFileName=linkFileName)
dfEEG              <- prepDFsignalProcessor(project,linkFileName=linkFileNameEEG)


animal        <- 'Cirque'

## adjust parameters to the animal ====
if (animal == 'Jesse'){
  dataSet       <- 'g'
  year          <- '2016'
  
  tdf[c(8:10,12),]# are bad
  valnd <- c(1,2,3,4,5,6,7,11,13,14) # for jesse
  
  
  chStr      <-  paste('ch', 5,sep='')
  #chStr <- paste('ch',1:33,sep='')
}

if (animal == 'Walter'){
  dataSet       <- 'g'
  year          <- '2016'
  
  #valnd      <- c(1:9, 11, 13:dim(tdf)[1])          # for walter
  valnd <- c(1:8,10,13,15)
  chStr      <-  paste('ch', 5,sep='')
}

if (animal == 'Sam'){
  dataSet       <- 'g'
  year          <- '2016'
  
  valnd <- c(1:3,5:dim(tdfEEG)[1])     # for sam
  valnd <- c(3,5,6,7)
  chStr      <-  paste('ch', 5,sep='')
  
}

if (animal == 'Rockey'){
  dataSet       <- 'g'
  year          <- '2016'
  valnd <- c(1:dim(tdf)[1])     # for rockey
  chStr      <-  paste('ch', 5,sep='')
}

if (animal == 'Ben'){
  dataSet       <- 'g'
  year          <- '2020'
  tdf        <- df[which(df$dataSet==dataSet&df$animal%in%animal&df$year%in%'2020'),]
  valnd <- c(1,3,4,5,6,8,9,11) # for ben, 2 was recorded at 5kS
  #valnd <- c(1,3:dim(tdf)[1])
  chStr      <-  paste('ch', 5,sep='')
}

if (animal == 'Fuzzy'){
  dataSet       <- 'a'
  year          <- '2022'
  
  valnd <- 1:30
  chStr      <-  paste('ch', 5,sep='')
}

if (animal == 'Cirque'){
  dataSet       <- 'a'
  year          <- ''
  
  chStr      <-  paste('ch', 5,sep='')
}



## select sessions ====
tdf           <- df[   which(df$dataSet   ==dataSet & df$animal%in%animal    & df$year%in%year),]
tdfEEG        <- dfEEG[which(dfEEG$dataSet==dataSet & dfEEG$animal%in%animal & dfEEG$year%in%year),]

#identical( tdf$session[valnd], tdfEEG$session[valnd])



# local Functions   ====
getOnset     <- function( ex, thrsh=1e-5, ptn=pval, chx=1 ){
  ons    <- max(251, max(ptn$taxis)+1)
  signdx <- which( ptn$lfp[ex,,chx]<thrsh )
  
  if( length(signdx)>0 )
    ons <- ptn$taxis[ min( signdx )]
  
  return(ons)
}
getOnsetMain <- function(pvres, thrsh=c(1e-50,1e-15,1e-8)){
  pint <- subs.data(pvres, pvres$xpp$effect=='Intercept')
  pamp <- subs.data(pvres, pvres$xpp$effect=='ampInd')
  pici <- subs.data(pvres, pvres$xpp$effect=='iciInd')
  
  tnd <- pvres$taxis<.5 # 
  pint$lfp[,tnd,] <- c(1)
  pamp$lfp[,tnd,] <- c(1)
  pici$lfp[,tnd,] <- c(1)
  
  ## find all channels with a valid multi-unit
  pvres$dfline$OnLat  <- sapply(1:dim(pint$lfp)[1], getOnset, ptn=pint, thrsh=thrsh[1])
  pvres$dfline$AmpLat <- sapply(1:dim(pint$lfp)[1], getOnset, ptn=pamp, thrsh=thrsh[2])
  pvres$dfline$ICILat <- sapply(1:dim(pint$lfp)[1], getOnset, ptn=pici, thrsh=thrsh[3])
  
  return(pvres$dfline)
}
onsetDistributionFigure <- function(pvres, thrsh=c(10e-50,10e-15,10e-8), mxLat=30){
  
  tdf <- getOnsetMain(pvres, thrsh=thrsh)
  
  tdf$use <- T
  tdf$use <- tdf$OnLat < mxLat
  useIndx <- which(tdf$use)
  
  ## make figures       ====
  
  fileName <- paste('/cumLatMain_logp_',-log10(thrsh[1]),'_',-log10(thrsh[2]),'_',-log10(thrsh[3]), '_lat_',mxLat,'.eps',sep='')
  eps.file( fileName,pdf=T,s=2)
  plot(   sort(tdf$OnLat[useIndx]),  seq(from=0,length=length(useIndx),to=1), type='l', col='black',xlim=c(0,40), xlab='Latency [ms]', ylab='Cummulative Probability')
  points( sort(tdf$AmpLat[useIndx]), seq(from=0,length=length(useIndx),to=1), type='l', col='blue' )
  points( sort(tdf$ICILat[useIndx]), seq(from=0,length=length(useIndx),to=1), type='l', col='red' )
  dev.off()
  
  
  
  ## only channels with onset for intersept, amp and ICI
  tdf$use <- tdf$OnLat < 50 & tdf$AmpLat < 50 & tdf$ICILat < 50
  useIndx <- which(tdf$use)
  
  
  # redo cumulative latency with units that have all three effects
  fileName <- paste('/cumLatMainSub_logp_',-log10(thrsh[1]),'_',-log10(thrsh[2]),'_',-log10(thrsh[3]), '_lat_',mxLat,'.eps',sep='')
  eps.file( fileName,pdf=T,s=2)
  plot(   sort(tdf$OnLat[useIndx]),  seq(from=0,length=length(useIndx),to=1), type='l', col='black',xlim=c(0,40), xlab='Latency [ms]', ylab='Cummulative Probability')
  points( sort(tdf$AmpLat[useIndx]), seq(from=0,length=length(useIndx),to=1), type='l', col='blue' )
  points( sort(tdf$ICILat[useIndx]), seq(from=0,length=length(useIndx),to=1), type='l', col='red' )
  dev.off()
  
  
  
  ## distributions for the three main effects
  fileName <- paste('/distLatMain_logp_',-log10(thrsh[1]),'_',-log10(thrsh[2]),'_',-log10(thrsh[3]), '_lat_',mxLat,'.eps',sep='')
  frv <- -5
  tov <- 40
  dI <- density( tdf$OnLat[useIndx], bw=1, from=frv, to=tov )
  dA <- density( tdf$AmpLat[useIndx], bw=1, from=frv, to=tov )
  dC <- density( tdf$ICILat[useIndx], bw=1, from=frv, to=tov )
  
  eps.file( fileName,pdf=T,s=2)
  plot(dI$x, dI$y, type='l', col='black',xlim=c(0,40), xlab='Latency [ms]', ylab='density')
  abline(h=0)
  points(dA$x, dA$y, type='l', col='blue')
  points(dC$x, dC$y, type='l', col='red')
  
  abline(v=dI$x[ which(dI$y==max(dI$y)) ], col='black', lty=2 )
  abline(v=dA$x[ which(dA$y==max(dA$y)) ], col='blue', lty=2 )
  abline(v=dC$x[ which(dC$y==max(dC$y)) ], col='red', lty=2 )
  dev.off()
  
  
  ttst <- t.test(tdf$AmpLat[useIndx],tdf$ICILat[useIndx] )
  fileName <- paste('/deltaLatMain_logp_',-log10(thrsh[1]),'_',-log10(thrsh[2]),'_',-log10(thrsh[3]), '_lat_',mxLat,'.eps',sep='')
  
  eps.file( fileName,pdf=T,s=2,r=2.2/2)
  tdf$sessionID <- tapply(tdf$session, tdf$session)
  plot( jitter(tdf$AmpLat[useIndx],amount=.5), jitter(tdf$ICILat[useIndx],amount=.5), xlim=c(0,40),ylim=c(0,40),pch=20,cex=.3, xlab='Amp latency [ms]', ylab='ICI latency [ms]')#, col=jet.colors(max(tdf$sessionID[useIndx]))[ tdf$sessionID[useIndx] ] )
  abline(a=0,b=1)
  text( 30,1, paste( 'p<10^-',  floor(-log10(ttst$p.value)), sep='') )
  dev.off()
  
  #dCA <- density( tdf$ICILat[useIndx]-tdf$AmpLat[useIndx], bw=1, from=-30, to=30)
  #plot(dCA$x, dCA$y, type='l', col='black',xlim=c(-30,30), xlab='Latency [ms]', ylab='density')
  #abline(h=0);abline(v=0)
  
  
  
  #my.boxplot(   c(table(tdf$OnLat[useIndx])), x=as.numeric(names(table(tdf$OnLat[useIndx]))), col='black',xlim=c(0,40), xlab='Latency [mx]', ylab='Number of instances')
  #axis(1)
  #my.boxplot(   c(table(tdf$AmpLat[useIndx])), x=as.numeric(names(table(tdf$AmpLat[useIndx]))), col='blue', add=T)
  #my.boxplot(   c(table(tdf$ICILat[useIndx])), x=as.numeric(names(table(tdf$ICILat[useIndx]))), col='red', add=T)
  
  
  
}

#chStr <- paste('ch',5,sep='')
#tst      <- loadSignalProcessor(  tdfEEG[24,],chStr=chStr, project='ampOddClick', epoch='bsp')#, analysis='ampInd_iciInd' )
#tst$channelPosition <- getEEGpos('samsym')
#tst      <- loadSignalProcessor(  tdf[valnd[9],],chStr=chStr, project='ampOddClick', epoch='bsp', analysis='ampInd_iciInd' )

#showSignalProcessor( tst, frml=~ampInd, subsExpr = expression(valTone & iciInd<6), trng=c(-10,50), showSE=F)
#showSignalProcessor( tst, frml=~ampInd, subsExpr = expression(valTone & iciInd<6), trng=c(-1,8), showSE=F)

#showSignalProcessor( tst, frml=~sessionID, subsExpr = expression(valTone & iciInd<6 & ampInd>3), trng=c(-1,10), showSE=F)

#valnd <- rep(T, dim(tdfEEG)[1])
mn       <- loadSignalProcessor(  tdfEEG[valnd,],chStr=chStr, project='ampOddClick', epoch='bsp', analysis='ampInd_iciInd', blrng=c(-2,0) )
pvl      <- loadSignalProcessor(  tdfEEG[valnd,],chStr=chStr, project='ampOddClick', epoch='bsp', analysis='trAOV_ampInd_iciInd_trialType_toneType_trialNr' , blrng=NA)
bet      <- loadSignalProcessor(  tdfEEG[valnd,],chStr=chStr, project='ampOddClick', epoch='bsp', analysis='trAOV_parm_ampInd_iciInd_trialType_toneType_trialNr' )


#showSignalProcessor( mn, frml=~sessionID, subsExpr = expression(iciInd<6), trng=c(-2,9) )

ga <- avgSignalProcessor( mn, frml=~sessionID, subsExpr = expression(iciInd<6))
showSignalProcessorXYC( ga, xaxis='sessionID', subsExpr = expression(iciInd<6), trng=c(-2,150)  )

pp <- function(){return(paste('~/Dropbox//ampOddClick/word/VProbe/plots/g_',animal,'/mnBSP/',sep='')) }

#showSignalProcessor(pvl, frml=~effect, subsExpr = expression(effect=='(Intercept)')) 
thrsh1 <- c(1e-50,1e-10,1e-10) # default MUA values
thrsh2 <- c(1e-50,1e-15,1e-8)  # modified MUA values to match significance rate of Amp and ICI
thrsh3 <- c(1e-20,1e-5,1e-5)   # default avRec values
thrsh4 <- c(1e-20,1e-4,1e-4)   # 


lat4 <- getOnsetMain(pvl,thrsh4)
valID <- which( lat4$OnLat < 15)


onsetDistributionFigure(pvl, thrsh=thrsh1)
onsetDistributionFigure(pvl, thrsh=thrsh2)
onsetDistributionFigure(pvl, thrsh=thrsh3)
onsetDistributionFigure(pvl, thrsh=thrsh4)

tcex <- 1.5

eps.file( 'bsp_byAmp_zoom.eps',pdf=T,s=3)
showSignalProcessor( mn, frml=~ampInd, trng=c(-5,25), bty = 'n', scl='numeric',sclfctr = 10,   butterfly=F, showSE=T, subsExpr = expression(iciInd<6 & sessionID%in%valID), tcex=tcex )
dev.off()

eps.file( 'bsp_byICI_zoom.eps',pdf=T,s=3)
showSignalProcessor( mn, frml=~iciInd, trng=c(-5,25), bty = 'l', scl='numeric',sclfctr = 10,   butterfly=F, showSE=T, colFUN=soa.colors2, subsExpr = expression(iciInd<6 & sessionID%in%valID), tcex=tcex )
dev.off()

eps.file( 'bsp_byAmp.eps',pdf=T,s=3)
showSignalProcessor( mn, frml=~ampInd, trng=c(-1,8), bty = 'n', scl='numeric',sclfctr = 3.5,   butterfly=F, showSE=T, subsExpr = expression(iciInd<6 & sessionID%in%valID), tcex=tcex )
dev.off()

eps.file( 'bsp_byICI.eps',pdf=T,s=3)
showSignalProcessor( mn, frml=~iciInd, trng=c(-1,8), bty = 'l', scl='numeric',sclfctr = 3.5,   butterfly=F, showSE=T, colFUN=soa.colors2, subsExpr = expression(iciInd<6 & sessionID%in%valID), tcex=tcex )
dev.off()

eps.file( 'bsp_byAmp_butterfly.eps',pdf=T,s=3)
showSignalProcessor( mn, frml=~ampInd, trng=c(-1,8), bty = 'l', yshft=-0.0, scl='numeric', sclfctr = 2, butterfly=T, showSE=T, subsExpr = expression(iciInd<6 & sessionID%in%valID), tcex=tcex)#, vlineat=c(mean(tst$AmpLat[valNdx])) )
dev.off()

eps.file( 'bsp_byICI_butterfly.eps',pdf=T,s=3)
showSignalProcessor( mn, frml=~iciInd, trng=c(-1,8), bty = 'l', yshft=-0.0, scl='numeric', sclfctr = 2, butterfly=T, showSE=T, colFUN=soa.colors2, subsExpr = expression(iciInd<6 & sessionID%in%valID), tcex=tcex )#, vlineat=c(mean(tst$ICILat[valNdx])) )
dev.off()

eps.file( 'bsp_byAmp_zoom_butterfly.eps',pdf=T,s=3)
showSignalProcessor( mn, frml=~ampInd, trng=c(-5,25), bty = 'l', yshft=-0.0, scl='numeric', sclfctr = 5, butterfly=T, showSE=T, subsExpr = expression(iciInd<6 & sessionID%in%valID), tcex=tcex)#, vlineat=c(mean(tst$AmpLat[valNdx])) )
dev.off()

eps.file( 'bsp_byICI_zoom_butterfly.eps',pdf=T,s=3)
showSignalProcessor( mn, frml=~iciInd, trng=c(-5,25), bty = 'l', yshft=-0.0, scl='numeric', sclfctr = 5, butterfly=T, showSE=T, colFUN=soa.colors2, subsExpr = expression(iciInd<6 & sessionID%in%valID), tcex=tcex )#, vlineat=c(mean(tst$ICILat[valNdx])) )
dev.off()


eps.file( 'beta_Amp_ICI.eps',pdf=T,s=3)
showSignalProcessor( bet, frml=~effect, trng=c(-5,25),scl='numeric', sclfctr=3, yshft=0, subsExpr=expression(effect%in%c('ampInd','iciInd') & sessionID%in%valID), bty='l',showSE=T, tcex=tcex )
dev.off()

