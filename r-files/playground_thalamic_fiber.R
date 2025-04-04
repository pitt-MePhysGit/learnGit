
## ========================================================= playground load sample SUA
dfline    <- chdf[1,]
chStr     <- paste(     'ch', 33+dfline$channel,sep='')
mua       <- loadAOC( dfline , what='mua', chStr=chStr, trda=trda, byBlock=F, suaMethod=suaMethod,suarng=c(-.1,.5)  )

raster.plot(mua$spk, xlim=c(-.05,.15))
plot( (1:length(mua$shp[1,]))/30 , mua$shp[1,], type='l', xlab='time [ms]', ylab='Voltage')


pp <- function(){'/Users/teichert/Dropbox/ampOdd_click/word/VProbe/'}
eps.file('exampleThalamicCell.eps',mf=c(1,2),pdf=T,s=3)
par(mfcol=c(2,1))
raster.plot(mua$spk, xlim=c(-.001,.020), ylab='trialNr', xlab='time [sec]')
#abline(h=1:100)

tax        <- seq(-.001,by=.0001,to=.020)
mua$psth   <- spk2psth(mua,eval.at=tax,width=.0005,rate=TRUE)
mean.plot( mua$psth, 1000*tax, add=F, widthFUN=na.se, ylab='Firing Rate [Hz]', xlab='time [ms]' )
abline(h=0)
dev.off()


mua$xpp$blFR  <- spike.count(mua, range=c(-0.05,-0.001),   rate=TRUE)
mua$xpp$evkFR <- spike.count(mua, range=c(0.0045, 0.0085), rate=FALSE)
mua$xpp$pstFR <- spike.count(mua, range=c(0.01, 0.040),    rate=TRUE)

tapply(mua$xpp$evkFR, mua$xpp$isiOct, mean)
plot(log2(mua$xpp$ISI), mua$xpp$evkFR)


dfline     <- chdf[5,]
chStr      <- paste(     'ch', 33+dfline$channel[n],sep='')
mua2       <- loadAOC( dfline , what='mua', chStr=chStr, trda=trda, byBlock=F, suaMethod=suaMethod,suarng=c(-.1,.5)  )

raster.plot(mua2$spk, xlim=c(-.001,.05))

