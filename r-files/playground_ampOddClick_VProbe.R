## Preliminaries ====
source('~/Dropbox/r-compile/signalProcessoR/R/signalProcessoR.R')

project         <- 'ampOddClick'
linkFileNameVP  <- 'ampOddclickVProbe_link.txt'   # only good sessions with VProbes
df              <- prepDFsignalProcessor(project,linkFileName=linkFileNameVP)


## Load a random data set ====
ndx <- 177 # fuzzy recording with unequal triggers
ndx <- 69  # sam recording re-processed with new signal processor
ndx <- 140 # walter recording not yet re-processed with new signal processor
ndx <- 198 # Fuzzy recording with wrong ICI index coding

tn <- loadSignalProcessor( df[ndx,], project = project, chStr=paste('ch',34:35,sep=''), epoch='mua' )

tst <- checkISISignalProcessor( tn , show=T)


showSignalProcessor( tn, frml=~ampInd )

table(tn$xpp$ampInd)
table(tn$xpp$iciInd)

plot(diff(tn$xpp$e.MLonset), diff(tn$xpp$e.onset) )


plot(     diff(tn$xpp$e.MLonset[1150:1180]), pch=19)
points(   lm1$coefficients[2]*diff(tn$xpp$e.onset[1150:1180]  ), col='green', pch=20)
