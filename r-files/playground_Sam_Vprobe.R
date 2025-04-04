## playground SAm VProbe ====

# check pre-processing for Sam's VProbe data

baseDir <- '~/Dropbox/'
source('~/Dropbox/r-compile/signalProcessoR/R/signalProcessoR.R')
project         <- 'ampOddClick'
linkFileNameVP  <- 'ampOddclickVProbe_link.txt'   # only good sessions with VProbes
df <- prepDFsignalProcessor( project, linkFileName = linkFileNameVP)

tdf <- df[df$animal=='Sam',]
table(tdf$ElectrodeConfig, tdf$dataSet)
table(tdf$ElectrodeConfig,tdf$probe)
table(tdf$dataSet,tdf$probe)


## mua exists 
tn <- loadSignalProcessor( tdf, epoch='mua', chStr='ch34', analysis='sessionID'  )
showSignalProcessor( tn, trng=c(-25,250) )
dim(tdf)
dim(tn$dfline)


## bsp is not consistent between sessions
# 1:2 have time=range out to 40ms
# 4                       to 15ms
# 16,17,20,34,36,27 ,48   to 40ms
# 54, 57, 58              to 15ms
tn <- loadSignalProcessor( tdf[60:67,], epoch='bsp', chStr='ch34', analysis='sessionID'  )
showSignalProcessor( tn )
range(tn$taxis)
1/mean(diff(tn$taxis))
dim(tn$dfline)
# rerun with fixed time-range to 15 ms


## hdraw AVG has not been run...
tn <- loadSignalProcessor( tdf[1:67,], epoch='hdraw', chStr='ch34', analysis='sessionID'  )
showSignalProcessor( tn )
range(tn$taxis)
1/mean(diff(tn$taxis))
dim(tn$dfline)
# rerun with fixed time-range to 15 ms




tst <- loadSignalProcessor( df[4,], epoch='mua', chStr='ch34'  )
tst$xpp$iciInd <- as.numeric(tst$xpp$iciInd)
av <- avgSignalProcessor( tst, frml=~ampInd+iciInd )
table(av$xpp$ampInd)
table(av$xpp$iciInd)

