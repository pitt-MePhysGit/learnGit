# averages for Siwei
rm(list=ls())
# Preliminaries ====
baseDir <- '~/Dropbox/'
source('~/Dropbox/r-compile/signalProcessoR/R/signalProcessoR.R')

project         <- 'ampOddClick'

linkFileName    <- 'ampOddclick_link.txt'        # includes sessions with SC arrays, and some old EEG only recordings of Ben Jesse Rockey Sam and Walter, wildly varying signal amplitude
linkFileNameEEG <- 'ampOddClickEEG_link.txt'     # the old EEG only recordings 


linkFileNameD   <- 'ampOddClickdual_link.txt'     # all sessions with VProbes
linkFileNameVP  <- 'ampOddclickVProbe_link.txt'   # only good sessions with VProbes


df      <- prepDFsignalProcessor(project,linkFileName=linkFileName)
dfEEG   <- prepDFsignalProcessor(project,linkFileName=linkFileNameEEG)

dfD     <- prepDFsignalProcessor(project,linkFileName=linkFileNameD)
dfVP    <- prepDFsignalProcessor(project,linkFileName=linkFileNameVP)

dfline <- dfD[which(dfD$animal=='Fuzzy' & dfD$date=='20220913'),]

tn <- loadSignalProcessor( dfline, epoch='mua', chStr=c('ch41','ch41'))


pp <- function(){return('~/Dropbox/ampOddClick/plots/forSiwei/')}
eps.file('grandAverage.eps',pdf=T,s=3,ratio=2/3)
showSignalProcessor( tn, frml=~sessionID, trng=c(-50,250), chIndx = 1, scl='numeric',sclfctr=8, bty='n', yshft = -.3)
dev.off()

eps.file('byAmpInd.eps',pdf=T,s=3,ratio=2/3)
showSignalProcessor( tn, frml=~ampInd,    trng=c(-50,250), chIndx = 1, scl='numeric',sclfctr=8, bty='n', yshft = -.3)
dev.off()

eps.file('byICIInd.eps',pdf=T,s=3,ratio=2/3)
showSignalProcessor( tn, frml=~iciInd,    trng=c(-50,250), chIndx = 1, scl='numeric',sclfctr=8, bty='n', yshft = -.3)
dev.off()


av <- avgSignalProcessor(tn, frml=~ampInd+iciInd)
eps.file('byAmp_and_ICI.eps',pdf=T,s=4,ratio=1)
showSignalProcessorXYC( av, xax='ampInd',yax='iciInd',    trng=c(-50,250), chIndx = 1, scl='numeric',sclfctr=8, bty='n', yshft = -.3)
dev.off()
