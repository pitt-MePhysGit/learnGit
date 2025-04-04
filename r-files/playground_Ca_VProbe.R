defpar <- par()

baseDir <- '~/Dropbox/'
source(paste(baseDir, 'RSkernel/r-files/functions/source_all_functions.R',sep=''))
source( paste(baseDir,'RSkernel/r-files/functions/FLA_VProbe_STPSP.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/prepDFVprobe.R',sep='') )


colClasses <- c('character',rep('numeric',3),'character',rep('numeric',4)  ,rep('character',3),rep('numeric',6) )
df         <- read.table(paste(baseDir,'RSkernel/RSkernelVProbe_link.txt',sep=''),header=T,colClasses=colClasses)
df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))

df <- prepDFVprobe(df)



## preliminaries
options(warn=0)
rda  <- function(){'~/Dropbox/RSkernel/Rdata/'}
trda <- function(){'/Volumes/Drobo5D3/EEG/RSkernel/rda/'}


#
animal  <- 'Sam'
dataSet <- 'g'
speaker <- 'TB'
drug    <- 'loCa'
ind     <- which(df$animal==animal & df$drug==drug,df$dataSet==dataSet & df$ElectrodeConfig=='ev')
tdf     <- df[ ind, ]
tdf$minTR <- NA
tdf

exnd      <- 1
chStr     <- paste('ch',34:57,sep='')
mua       <- loadRS( tdf[exnd,] , what='mua', chStr=chStr, trda=trda , blcorr=T )
lfp       <- loadRS( tdf[exnd,] , what='raw', chStr=chStr, trda=trda , blcorr=F  )

dxt <- ifelse(tdf$drug[exnd]=='N','', paste('_',tdf$drug[exnd],sep=''))
pxt <- ifelse(tdf$probe[exnd]==1 ,'', paste('_',tdf$probe[exnd],sep=''))

.GlobalEnv$pp   <- function(){paste(baseDir, 'RSkernel/word/VProbe/plots/',tdf$session[exnd],dxt,pxt,'/',sep='')}
print(pp())

mua$xpp$ISI[ which( mua$xpp$ISI<0)]  <- 20
lfp$xpp$ISI[ which( lfp$xpp$ISI<0)]  <- 20


## ========================================================= p2p mua
tmp          <- apply(mua$lfp,c(1,3),max) - apply(mua$lfp,c(1,3),min)
mua$xpp$p2p  <- apply( tmp, 1, max )

qtmp         <- quantile(mua$xpp$p2p,probs=c(0.5,0.75))
mxmua        <- qtmp[1]+10*diff(qtmp)
valIndMUA    <- mua$xpp$p2p < mxmua

trng         <- c(20,40)
chrng        <- c( ceiling( tdf$Layer[[exnd]]):24) 
tnd          <- which(mua$taxis >= trng[1] & mua$taxis < trng[2])
mua$xpp$mua  <- apply(mua$lfp[,tnd,chrng], c(1), mean)


## =========================================== calculate single-trial csd
lag     <- 1
smooth  <- T
sd      <- 100
spacing <- tdf[exnd,]$spacing
varCor  <- F

csd      <- calcSingleTrialCSD(lfp, spacing=spacing, lag=lag,smooth=T, sd=sd)
csd$what <- 'csd'
csd$xpp$ISI[ which( csd$xpp$ISI<0)]  <- 20

avrec     <- csd
avrec$lfp <- abs(csd$lfp)
tmp       <- baseLineCorrection(avrec$lfp,c(-50,0), csd$taxis)
avrec$lfp <- tmp$res

# ----------------
tmp          <- apply(avrec$lfp,c(1,3),max) - apply(avrec$lfp,c(1,3),min)
avrec$xpp$p2p  <- apply( tmp, 1, max )

qtmp         <- quantile(avrec$xpp$p2p,probs=c(0.5,0.75))
mxcsd        <- qtmp[1]+10*diff(qtmp)
valIndCSD    <- avrec$xpp$p2p < mxcsd


trng            <- c(25,45)
chrng           <- 1:24#c( (ceiling( tdf$Layer[[exnd]])-4) :  (ceiling( tdf$Layer[[exnd]])+4)) 
tnd             <- which(avrec$taxis >= trng[1] & avrec$taxis < trng[2])
avrec$xpp$avrec <- apply( avrec$lfp[,tnd,chrng], c(1), mean)


#tmp  <- apply( csd$lfp[valIndCSD,,], c(2,3), mean)
#matplot(tmp,tax=csd$taxis, dax=NULL, add=F, type='l', ylim=c(0,24), xlim=c(-60,250), yscl=0, lty=2,lwd=1)
#matplot(tmp,tax=csd$taxis, dax=NULL, add=T, type='l', ylim=c(0,24), xlim=c(-60,250), yscl=2,lwd=2)

#tmp  <- apply( abs(csd$lfp[valIndCSD,,]), c(2), mean)
#plot(csd$taxis, tmp, xlim=c(-50,250), type='l')

tmp  <- apply( avrec$lfp[avrec$xpp$trialNr<4.5&valIndCSD&avrec$xpp$ISI<1,,], c(2), mean)
plot(avrec$taxis, tmp, xlim=c(-50,250), type='l', ylim=c(-.02,.2))

tmp  <- apply( avrec$lfp[avrec$xpp$trialNr<4.5&valIndCSD&avrec$xpp$ISI>2,,], c(2), mean)
points(avrec$taxis, tmp, xlim=c(-50,250), type='l', col='red')


## =================================================================
valInd <- which( valIndCSD & csd$xpp$trialNr<4.5)
soax   <- tapply( log2(avrec$xpp$ISI[valInd]), avrec$xpp$isiOct[valInd], mean)
mnAR1  <- tapply( avrec$xpp$avrec[valInd], avrec$xpp$isiOct[valInd], mean)
seAR1  <- tapply( avrec$xpp$avrec[valInd], avrec$xpp$isiOct[valInd], na.se)

valInd <- which( valIndCSD & csd$xpp$trialNr>4.5)
mnAR2  <- tapply( avrec$xpp$avrec[valInd], avrec$xpp$isiOct[valInd], mean)
seAR2  <- tapply( avrec$xpp$avrec[valInd], avrec$xpp$isiOct[valInd], na.se)

vsoa <- 1:6

ylim <- range(c(mnAR1[vsoa], mnAR2[vsoa]))
errorplot(soax[vsoa], mnAR1[vsoa], ey=seAR1[vsoa], pch=20, lwd=2, col='black', ylim=ylim  )
errorplot(soax[vsoa], mnAR2[vsoa], ey=seAR2[vsoa], pch=20, lwd=2, col='red', add=T  )


## =================================================================
valInd <- which( valIndMUA & mua$xpp$trialNr<4.5)
soax <- tapply( log2(mua$xpp$ISI[valInd]), mua$xpp$isiOct[valInd], mean)
mnMUA1 <- tapply( mua$xpp$mua[valInd], mua$xpp$isiOct[valInd], mean)
seMUA1 <- tapply( mua$xpp$mua[valInd], mua$xpp$isiOct[valInd], na.se)

valInd <- which( valIndMUA & mua$xpp$trialNr>4.5)
mnMUA2 <- tapply( mua$xpp$mua[valInd], mua$xpp$isiOct[valInd], mean)
seMUA2 <- tapply( mua$xpp$mua[valInd], mua$xpp$isiOct[valInd], na.se)

ylim <- range(c(mnMUA1[vsoa], mnMUA2[vsoa]))
errorplot(soax[vsoa], mnMUA1[vsoa], ey=seMUA1[vsoa], pch=20, lwd=2, col='black', ylim=ylim  )
errorplot(soax[vsoa], mnMUA2[vsoa], ey=seMUA2[vsoa], pch=20, lwd=2, col='red', add=T  )







## =========================
lfpscl  <- 4
csdscl  <- 300#500#1300
mxmua <- 40
mxlfp <- 400
interpChanInd <- c();#c(16)

csdGA   <- calcCSD(lfp, mua, lag=lag, mxmua=mxmua,mxlfp=mxlfp,interpChanInd=interpChanInd,smooth=smooth,sd=sd,varCor=varCor,spacing=spacing)
disc <- showCSDImage(csdGA,trng=c(-10,250),lfpscl=lfpscl,csdscl=csdscl,tat=NA,addcsd=T)







