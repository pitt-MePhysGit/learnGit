## playground_RFMmodel
#
baseDir <- '~/Dropbox/'
source(paste(baseDir, 'ampOdd_click/r-files/source_all_functions.R',sep=''))

## preliminaries
pp <- function(){paste(baseDir, 'ampOdd_click/word/plots/',sep='')}
rda <- function(){'/Volumes/rawData/EEG/ampOddClick/rda/'}
rda <- function(){'~/Dropbox/ampOdd_click/Rdata/'}

## load multiple data sets for a first-order analysis

what    <-'raw'
animal  <- 'Walter'
dataSet <- 'k'
chStr   <- 'frontoCentral'
#chStr <- 'p1a'
df      <- df[which(df$dataSet==dataSet),]

#tone        <- loadAmpOddClickAll( tdf , what='raw', chStr=chStr,show=F)
#tone$exp    <- paste(animal,dataSet,sep='_')
#tone        <- subs.data(tone,   tone$xpp$p2p<450 & (!is.na(tone$xpp$isiOct)) )##& abs(tone$xpp$blcor)<120  )
#table(tone$xpp$exp,tone$xpp$trialType)
#table(tone$xpp$toneType,tone$xpp$trialType)

eps      <- F
doseList <- 0.05
xdat     <- showISIbyKet(animal, dataSet=dataSet, speaker='TB', what='raw', eps=eps,doseList=doseList,df=df,chStr=chStr)  
#sdat     <- showISIbyKet('Sam',    dataSet=dataSet, speaker='TB', what='raw', eps=eps,doseList=doseList,df=df)


par(mfcol=c(1,1))
## ====================== split by amp (y-axis) and soa (x-axis)
#xdat <- sdat


mxp2p  <- 2000
lwd    <- 1
NAmp   <- length(unique(xdat[[1]]$xpp$ampInd))
NSOA   <- length(levels(xdat[[1]]$xpp$isiOct))
xscale <- 0.01

eps.file( paste('ampOddClick_',animal,'_',dataSet,'_',chStr,'.eps',sep=''), pdf=T )
aax <- c(1,2,3,4,5,6.25)
sax <- c(1,2,3,4,5,6.25) 
plot(-10,-10, xlim=c(0.75,NSOA+.25), ylim=c(-0.5,NAmp+.5),xaxt='n', yaxt='n' )
axis(1,at=sax+0.25, label=1:(NSOA-1) )
axis(2,at=aax-1, label=1:(NAmp+1) )
valTimes <- which(xdat[[1]]$taxis > -50 & xdat[[1]]$taxis<250)
for(sind in 1:(NSOA-1)){
  for(aind in 1:(NAmp+1)){
    
    lines( sax[sind]+c(0,0.8), rep( aax[aind]-1, 2), col='black' )
    lines( sax[sind]+squeeze( c(0,0),from.range=c(-50,250), to.range=c(0,0.8)), c(aax[aind]-1.25, aax[aind]-.75), col='black' )
    
    valInd <- which( xdat[[1]]$xpp$ampInd == aind & xdat[[1]]$xpp$isiOct == levels(xdat[[1]]$xpp$isiOct)[sind] & xdat[[1]]$xpp$p2p<mxp2p & xdat[[1]]$xpp$trialNr>4)
    
    if (sind == NSOA-1)
      valInd <- which( xdat[[1]]$xpp$ampInd == aind & xdat[[1]]$xpp$ISI < 6.4 & xdat[[1]]$xpp$p2p<mxp2p & xdat[[1]]$xpp$trialNr>4)
    
    if (aind == NAmp+1)
      valInd <- which( xdat[[1]]$xpp$ampInd < aind & xdat[[1]]$xpp$isiOct == levels(xdat[[1]]$xpp$isiOct)[sind] & xdat[[1]]$xpp$p2p<mxp2p & xdat[[1]]$xpp$trialNr>4)
    
    if (sind == NSOA-1 & aind==NAmp+1)
      valInd <- which( xdat[[1]]$xpp$ampInd < aind & xdat[[1]]$xpp$ISI < 6.4 & xdat[[1]]$xpp$p2p<mxp2p & xdat[[1]]$xpp$trialNr>4)
    
    disc   <- mean.plot( (aax[aind]-1)+xscale*xdat[[1]]$lfp[valInd,valTimes,1], sax[sind]+squeeze(xdat[[1]]$taxis[valTimes],to.range=c(0,0.8)), add=T, sd=TRUE, traces=F, 
                         posFUN=na.mean, widthFUN=na.se, range=NULL, col='black',lwd=lwd)
    
    
    ## ====================== MK-801
    valInd <- which( xdat[[2]]$xpp$ampInd == aind & xdat[[2]]$xpp$isiOct == levels(xdat[[2]]$xpp$isiOct)[sind] & xdat[[2]]$xpp$p2p<mxp2p & xdat[[2]]$xpp$trialNr>4)
    
    if (sind == NSOA-1)
      valInd <- which( xdat[[2]]$xpp$ampInd == aind & xdat[[2]]$xpp$ISI < 6.4 & xdat[[2]]$xpp$p2p<mxp2p & xdat[[2]]$xpp$trialNr>4)
    
    if (aind == NAmp+1)
      valInd <- which( xdat[[2]]$xpp$ampInd < aind & xdat[[2]]$xpp$isiOct == levels(xdat[[2]]$xpp$isiOct)[sind] & xdat[[2]]$xpp$p2p<mxp2p & xdat[[2]]$xpp$trialNr>4)
    
    if (sind == NSOA-1 & aind==NAmp+1)
      valInd <- which( xdat[[2]]$xpp$ampInd < aind & xdat[[2]]$xpp$ISI < 6.4 & xdat[[2]]$xpp$p2p<mxp2p & xdat[[2]]$xpp$trialNr>4)
    
    
    disc   <- mean.plot( (aax[aind]-1)+xscale*xdat[[2]]$lfp[valInd,valTimes,1], sax[sind]+squeeze(xdat[[2]]$taxis[valTimes],to.range=c(0,0.8)), add=T, sd=TRUE, traces=F, 
                         posFUN=na.mean, widthFUN=na.se, range=NULL, col='red',lwd=lwd)
  }
}
dev.off()


## effect on LDAEP
eps.file( paste('ampOddClick_LDAEP_',animal,'_',dataSet,'_',chStr,'.eps',sep=''), pdf=T, ratio=3/2 )
plot(-10,-10, xlim=c(0,1), ylim=c(-0.5,NAmp-.5) )
for(aind in 1:(NAmp)){
  
  lines( c(0,0.8), rep( aind-1, 2), col='black' )
  lines( squeeze( c(0,0),from.range=c(-50,250), to.range=c(0,0.8)), c(aind-1.25, aind-.75), col='black' )
  
  valInd <- which( xdat[[1]]$xpp$ampInd == aind & xdat[[1]]$xpp$p2p<mxp2p & xdat[[1]]$xpp$trialNr>4)
  disc   <- mean.plot( (aind-1)+xscale*xdat[[1]]$lfp[valInd,valTimes,1], squeeze(xdat[[1]]$taxis[valTimes],to.range=c(0,0.8)), add=T, sd=TRUE, traces=F, 
                       posFUN=na.mean, widthFUN=na.se, range=NULL, col='black',lwd=lwd)
  
  valInd <- which( xdat[[2]]$xpp$ampInd == aind  & xdat[[2]]$xpp$p2p<mxp2p & xdat[[2]]$xpp$trialNr>4)
  disc   <- mean.plot( (aind-1)+xscale*xdat[[2]]$lfp[valInd,valTimes,1], squeeze(xdat[[2]]$taxis[valTimes],to.range=c(0,0.8)), add=T, sd=TRUE, traces=F, 
                       posFUN=na.mean, widthFUN=na.se, range=NULL, col='red',lwd=lwd)
}
dev.off()

# effect on SOA scaling
eps.file( paste('ampOddClick_TDAEP_',animal,'_',dataSet,'_',chStr,'.eps',sep=''), pdf=T, ratio=3/2 )
plot(-10,-10, xlim=c(0,1), ylim=c(-0.5,NSOA-2.5) )
for(sind in 1:(NSOA-2)){
  
  lines( c(0,0.8), rep( sind-1, 2), col='black' )
  lines( squeeze( c(0,0),from.range=c(-50,250), to.range=c(0,0.8)), c(sind-1.25, sind-.75), col='black' )
  
  valInd <- which( xdat[[1]]$xpp$isiOct == levels(xdat[[1]]$xpp$isiOct)[sind] & xdat[[1]]$xpp$p2p<mxp2p & xdat[[1]]$xpp$trialNr>4)
  disc   <- mean.plot( (sind-1)+xscale*xdat[[1]]$lfp[valInd,valTimes,1], squeeze(xdat[[1]]$taxis[valTimes],to.range=c(0,0.8)), add=T, sd=TRUE, traces=F, 
                       posFUN=na.mean, widthFUN=na.se, range=NULL, col='black',lwd=lwd)
  
  valInd <- which( xdat[[2]]$xpp$isiOct == levels(xdat[[2]]$xpp$isiOct)[sind] & xdat[[2]]$xpp$p2p<mxp2p & xdat[[2]]$xpp$trialNr>4)
  disc   <- mean.plot( (sind-1)+xscale*xdat[[2]]$lfp[valInd,valTimes,1], squeeze(xdat[[2]]$taxis[valTimes],to.range=c(0,0.8)), add=T, sd=TRUE, traces=F, 
                       posFUN=na.mean, widthFUN=na.se, range=NULL, col='red',lwd=lwd)
}
dev.off()

