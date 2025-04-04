project      <- 'ampOddClick'
linkFileName <- 'ampOddClick_link.txt'
index        <- c(109)

## 
baseDir <- '~/Dropbox/'
source('~/Dropbox/r-compile/signalProcessoR/R/signalProcessoR.R')


library(fastICA)

source( paste(baseDir,'RSkernel/r-files/functions/call_fastICA.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/save_fastICA.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/lamiarLFPtheory/showICAcomponents.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/removeICAcomponent.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/tapered_ccf.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/lamiarLFPtheory/smooth2DCSD.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/removeCommonReferenceNoise_ica.R',sep='') )




df     <- prepDFsignalProcessor(project,linkFileName=linkFileName)
dfline <- df[index,]
print(dfline)

pp <- function(){paste(baseDir, project, '/plots/',dfline$session,'/', sep='' )}
if (!file.exists( pp() ))
  system(paste('mkdir ', pp(), sep='') )


# =============== raw
raw                 <- loadSignalProcessor(dfline, project=project, epoch='raw', chStr=paste('ch',c(34:257),sep=''),analysis='ampInd' )
raw$channelPosition <- getSC224pos(sc96ExtnNr = '_Intan')

eps.file( 'ampOddClick_multiplot_raw_overview.pdf', s=6, ratio=2 )
showSignalProcessor(raw, subsExpr = expression(valTone&ampInd<=5), trng=c(-50,250),frml=~ampInd )
dev.off()

mua                 <- loadSignalProcessor(dfline, project=project, epoch='mua', chStr=paste('ch',c(34:257),sep=''),analysis='ampInd' )
mua$channelPosition <- getSC224pos(sc96ExtnNr = '_Intan')

eps.file( 'ampOddClick_multiplot_mua_overview.pdf', s=6, ratio=2 )
showSignalProcessor(mua, subsExpr = expression(valTone&ampInd<=5), trng=c(-50,250),frml=~ampInd )
dev.off()




