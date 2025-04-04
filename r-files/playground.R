## playground_RFMmodel
#
baseDir <- '~/Dropbox/'
source(paste(baseDir, 'ampOddClick/r-files/source_all_functions.R',sep=''))

## preliminaries
pp <- function(){paste(baseDir, 'ampOddClick/word/plots/',sep='')}
rda <- function(){'/Volumes/rawData/EEG/ampOddClick/rda/'}
rda <- function(){'~/Dropbox/ampOddClick/Rdata/'}

## load multiple data sets for a first-order analysis

what   <-'raw'
animalList <- c('Rockey','Jesse','Walter','Sam')
#animalList <- 'Jesse'
dataSet <- 'g'
chStr='frontoCentral'
#chStr = 'p2'
tdf <- df[which(df$dataSet==dataSet&df$animal%in%animalList),]

cleanOnly <- FALSE

animalStr <- 'all'
if (length(animalList)==1)
  animalStr <- animalList[1]


tone        <- loadAmpOddClickAll( tdf , what='raw', chStr=chStr,show=F)
tone$exp    <- paste(animalStr,dataSet,sep='_')
tone        <- subs.data(tone,   tone$xpp$p2p<450 & (!is.na(tone$xpp$isiOct)) )##& abs(tone$xpp$blcor)<120  )
table(tone$xpp$exp,tone$xpp$trialType)
table(tone$xpp$toneType,tone$xpp$trialType)

tmp            <- baseLineCorrection(tone$lfp,c(-15,0), tone$taxis)
tone$lfp       <- tmp$res

toff                                   <- 0
tone$logtaxis                         <- log2( tone$taxis + toff)
tone$logtaxis[ (tone$taxis + toff <= 1)] <- .01*(tone$taxis[ which(tone$taxis + toff <= 1) ] - 1 + toff) #- log2(taxis[min(which(taxis>-toff))] )# + .01*taxis[ (llp$taxis + toff == 0)]
tone$logtaxis <- as.vector(tone$logtaxis)

## ======================================== set folder for plots

cleanStr  <- ''
if (cleanOnly)
  cleanStr <- '_clean'

uIdir <- paste('g', '_', animalStr,cleanStr,sep='')
.GlobalEnv$pp   <- function(){paste(baseDir, '/ampOddClick/word/plots/',uIdir,'/',sep='')}
if (!file.exists(pp()))
  system( paste('mkdir', pp()) )




## ==================================== ERP by SOA
xlim <- c(-50,450)
ylim <- c(-45,30)
lwd <- 1

Nbin <- length(levels(tone$xpp$isiOct))
colList <- soa.colors(Nbin-1)
eps.file(paste('GA_SOA.eps'), pdf=T, s=2, ratio=3/2)
for(n in 1:(Nbin-2)){
  valInd <- which( tone$xpp$isiOct == levels(tone$xpp$isiOct)[n]  )
  print(tone$xpp$isiOct[valInd[1]])
  print(length(valInd))
  
  disc   <- mean.plot( tone$lfp[valInd,,1], tone$taxis, add=ifelse(n==1,F,T), sd=TRUE, traces=F, 
                       posFUN=na.mean, widthFUN=na.se, range=NULL, col=colList[n],lwd=lwd, ylim=ylim,
                       xlim=xlim)
}
abline(h=0)
dev.off()


Nbin <- length(unique(tone$xpp$ampInd))
colList <- bluered.colors(Nbin)
eps.file(paste('GA_AMP.eps'), pdf=T, s=2, ratio=3/2)
for(n in 1:(Nbin)){
  valInd <- which( tone$xpp$ampInd == n  )
  print(tone$xpp$ampInd[valInd[1]])
  print(length(valInd))
  
  disc   <- mean.plot( tone$lfp[valInd,,1], tone$taxis, add=ifelse(n==1,F,T), sd=TRUE, traces=F, 
                       posFUN=na.mean, widthFUN=na.se, range=NULL, col=colList[n],lwd=lwd, ylim=ylim,
                       xlim=xlim)
}
abline(h=0)
dev.off()


# zoomed version
xlim <- c(-10,150)
Nbin <- length(levels(tone$xpp$isiOct))
colList <- soa.colors(Nbin-1)

eps.file(paste('GA_SOA_zoom.eps'), pdf=T, s=2, ratio=3/2)
for(n in 1:(Nbin-2)){
  valInd <- which( tone$xpp$isiOct == levels(tone$xpp$isiOct)[n]  )
  print(tone$xpp$isiOct[valInd[1]])
  print(length(valInd))
  
  disc   <- mean.plot( tone$lfp[valInd,,1], tone$taxis, add=ifelse(n==1,F,T), sd=TRUE, traces=F, 
                       posFUN=na.mean, widthFUN=na.se, range=NULL, col=colList[n],lwd=lwd, ylim=ylim,
                       xlim=xlim)
}
abline(h=0)
dev.off()


Nbin <- length(unique(tone$xpp$ampInd))
colList <- bluered.colors(Nbin)
eps.file(paste('GA_AMP_zoom.eps'), pdf=T, s=2, ratio=3/2)
for(n in 1:(Nbin)){
  valInd <- which( tone$xpp$ampInd == n  )
  print(tone$xpp$ampInd[valInd[1]])
  print(length(valInd))
  
  disc   <- mean.plot( tone$lfp[valInd,,1], tone$taxis, add=ifelse(n==1,F,T), sd=TRUE, traces=F, 
                       posFUN=na.mean, widthFUN=na.se, range=NULL, col=colList[n],lwd=lwd, ylim=ylim,
                       xlim=xlim)
}
abline(h=0)
dev.off()



## ============================== log version
trng <- c(3,550)
ylim <- c(-45,30)
xlim <- range( tone$logtaxis[ tone$taxis>=trng[1] & tone$taxis<trng[2]] )

Nbin <- length(levels(tone$xpp$isiOct))
colList <- soa.colors(Nbin-1)
eps.file(paste('GA_SOA_log.eps'), pdf=T, s=2, ratio=3/2)
for(n in 1:(Nbin-2)){
  valInd <- which( tone$xpp$isiOct == levels(tone$xpp$isiOct)[n]  )
  print(tone$xpp$isiOct[valInd[1]])
  print(length(valInd))
  
  disc   <- mean.plot( tone$lfp[valInd,,1], tone$logtaxis, add=ifelse(n==1,F,T), sd=TRUE, traces=F, 
                       posFUN=na.mean, widthFUN=na.se, range=NULL, col=colList[n],lwd=lwd, ylim=ylim,
                       xlim=xlim, xaxt='n')
}
abline(h=0)
taxat <- c(10,20,30,40,50,100,200,300,400)
taxnd <- which( round(tone$taxis)%in%c(taxat) )
axis(1, at=tone$logtaxis[taxnd], taxat  )
dev.off()


Nbin <- length(unique(tone$xpp$ampInd))
colList <- bluered.colors(Nbin)
eps.file(paste('GA_AMP_log.eps'), pdf=T, s=2, ratio=3/2)
for(n in 1:(Nbin)){
  valInd <- which( tone$xpp$ampInd == n  )
  print(tone$xpp$ampInd[valInd[1]])
  print(length(valInd))
  
  disc   <- mean.plot( tone$lfp[valInd,,1], tone$logtaxis, add=ifelse(n==1,F,T), sd=TRUE, traces=F, 
                       posFUN=na.mean, widthFUN=na.se, range=NULL, col=colList[n],lwd=lwd, ylim=ylim,
                       xlim=xlim,xaxt='n')
}
abline(h=0)
taxat <- c(10,20,30,40,50,100,200,300,400)
taxnd <- which( round(tone$taxis)%in%c(taxat) )
axis(1, at=tone$logtaxis[taxnd], taxat  )
dev.off()



## ====================== subtract out median amplitude / SOAbin
xlim <- c(-10,150)
ylim <- c(-20,20)

Nbin <- length(levels(tone$xpp$isiOct))
colList <- soa.colors(Nbin-1)
medSOAind <- 3
medInd <- which( tone$xpp$isiOct == levels(tone$xpp$isiOct)[medSOAind]  )
subtract <- apply( tone$lfp[medInd,,1],2,mean  )

eps.file(paste('GA_Delta_SOA.eps'), pdf=T, s=2, ratio=3/2)
for(n in 1:(Nbin-2)){
  valInd <- which( tone$xpp$isiOct == levels(tone$xpp$isiOct)[n]  )
  disc   <- mean.plot( tone$lfp[valInd,,1]-array(rep(subtract,each=length(valInd)),c(length(valInd),length(tone$taxis))), tone$taxis, add=ifelse(n==1,F,T), sd=TRUE, traces=F, 
                       posFUN=na.mean, widthFUN=na.se, range=NULL, col=colList[n],lwd=lwd, ylim=ylim,
                       xlim=xlim)
}
abline(h=0)
abline(v=0)
dev.off()


Nbin <- length(unique(tone$xpp$ampInd))
colList <- bluered.colors(Nbin)
medAMPind <- 3
medInd <- which( tone$xpp$ampInd == medAMPind  )
subtract <- apply( tone$lfp[medInd,,1],2,mean  )

eps.file(paste('GA_Delta_AMP.eps'), pdf=T, s=2, ratio=3/2)
for(n in 1:(Nbin)){
  valInd <- which( tone$xpp$ampInd == n  )
  print(tone$xpp$ampInd[valInd[1]])
  print(length(valInd))
  
  disc   <- mean.plot( tone$lfp[valInd,,1]-array(rep(subtract,each=length(valInd)),c(length(valInd),length(tone$taxis))), tone$taxis, add=ifelse(n==1,F,T), sd=TRUE, traces=F, 
                       posFUN=na.mean, widthFUN=na.se, range=NULL, col=colList[n],lwd=lwd, ylim=ylim,
                       xlim=xlim)
}
abline(h=0)
abline(v=0)
dev.off()


## ============== with logtime
trng <- c(3,550)
ylim <- c(-20,20)
xlim <- range( tone$logtaxis[ tone$taxis>=trng[1] & tone$taxis<trng[2]] )

Nbin <- length(levels(tone$xpp$isiOct))
colList <- soa.colors(Nbin-1)
medSOAind <- 3
medInd <- which( tone$xpp$isiOct == levels(tone$xpp$isiOct)[medSOAind]  )
subtract <- apply( tone$lfp[medInd,,1],2,mean  )

eps.file(paste('GA_Delta_SOA_log.eps'), pdf=T, s=2, ratio=3/2)
for(n in 1:(Nbin-2)){
  valInd <- which( tone$xpp$isiOct == levels(tone$xpp$isiOct)[n]  )
  disc   <- mean.plot( tone$lfp[valInd,,1]-array(rep(subtract,each=length(valInd)),c(length(valInd),length(tone$taxis))), tone$logtaxis, add=ifelse(n==1,F,T), sd=TRUE, traces=F, 
                       posFUN=na.mean, widthFUN=na.se, range=NULL, col=colList[n],lwd=lwd, ylim=ylim,
                       xlim=xlim,xaxt='n')
}
abline(h=0)
taxat <- c(10,20,30,40,50,100,200,300,400)
taxnd <- which( round(tone$taxis)%in%c(taxat) )
axis(1, at=tone$logtaxis[taxnd], taxat  )
dev.off()


Nbin <- length(unique(tone$xpp$ampInd))
colList <- bluered.colors(Nbin)
medAMPind <- 3
medInd <- which( tone$xpp$ampInd == medAMPind  )
subtract <- apply( tone$lfp[medInd,,1],2,mean  )

eps.file(paste('GA_Delta_AMP_log.eps'), pdf=T, s=2, ratio=3/2)
for(n in 1:(Nbin)){
  valInd <- which( tone$xpp$ampInd == n  )
  print(tone$xpp$ampInd[valInd[1]])
  print(length(valInd))
  
  disc   <- mean.plot( tone$lfp[valInd,,1]-array(rep(subtract,each=length(valInd)),c(length(valInd),length(tone$taxis))), tone$logtaxis, add=ifelse(n==1,F,T), sd=TRUE, traces=F, 
                       posFUN=na.mean, widthFUN=na.se, range=NULL, col=colList[n],lwd=lwd, ylim=ylim,
                       xlim=xlim,xaxt='n')
}
abline(h=0)
taxat <- c(10,20,30,40,50,100,200,300,400)
taxnd <- which( round(tone$taxis)%in%c(taxat) )
axis(1, at=tone$logtaxis[taxnd], taxat  )
dev.off()




## ================= min - max condition combine both SOA and AMP in one plot
deltaAMP <- apply( tone$lfp[ which(tone$xpp$ampInd==5),,1],2,mean  ) - apply( tone$lfp[ which(tone$xpp$ampInd==2),,1],2,mean  )
shortInd <- which( tone$xpp$isiOct == levels(tone$xpp$isiOct)[2]  )
longInd <- which( tone$xpp$isiOct == levels(tone$xpp$isiOct)[5]  )
deltaSOA <- apply( tone$lfp[ longInd,,1],2,mean  ) - apply( tone$lfp[ shortInd,,1],2,mean  )

eps.file(paste('GA_Delta_AMPSOA.eps'), pdf=T, s=2, ratio=3/2)
plot(  tone$taxis, deltaAMP, type='l', col='blue', xlim=c(-10,150), ylim=c(-40,15))
points(tone$taxis, deltaSOA, type='l', col='red')
abline(h=0)
abline(v=0)
dev.off()

eps.file(paste('GA_Delta_AMPSOA_zoom.eps'), pdf=T, s=2, ratio=3/2)
plot(  tone$taxis, deltaAMP, type='l', col='blue', xlim=c(-10,50), ylim=c(-5,15))
points(tone$taxis, deltaSOA, type='l', col='red')
abline(h=0)
abline(v=0)
dev.off()


trng <- c(3,550)
ylim <- c(-45,30)
xlim <- range( tone$logtaxis[ tone$taxis>=trng[1] & tone$taxis<trng[2]] )
eps.file(paste('GA_Delta_AMPSOA_log.eps'), pdf=T, s=2, ratio=3/2)
plot(  tone$logtaxis, deltaAMP, type='l', col='blue', xlim=xlim, ylim=c(-40,15),xaxt='n')
points(tone$logtaxis, deltaSOA, type='l', col='red')
abline(h=0)
taxat <- c(10,20,30,40,50,100,200,300,400)
taxnd <- which( round(tone$taxis)%in%c(taxat) )
axis(1, at=tone$logtaxis[taxnd], taxat  )
dev.off()



## ====================== split by amp (y-axis) and soa (x-axis)
NAmp <- length(unique(tone$xpp$ampInd))
NSOA <- length(levels(tone$xpp$isiOct))
xscale <- 0.01
plot(-10,-10, xlim=c(0,NSOA+1), ylim=c(-0.5,NAmp) )
for(sind in 1:(NSOA-1)){
  for(aind in 1:(NAmp)){
    
    valInd <- which( tone$xpp$ampInd == aind & tone$xpp$isiOct == levels(tone$xpp$isiOct)[sind] & tone$xpp$p2p<450)
    valTimes <- which(tone$taxis > -50 & tone$taxis<250)
    
    disc   <- mean.plot( (aind-1)+xscale*tone$lfp[valInd,valTimes,1], sind+squeeze(tone$taxis[valTimes],to.range=c(0,0.8)), add=T, sd=TRUE, traces=F, 
                         posFUN=na.mean, widthFUN=na.se, range=NULL, col='black',lwd=lwd, ylim=ylim,
                         xlim=xlim)
  }
}
abline(h=0)




break

N <- table(tone$xpp$isiOctFine)
vi <- which(N>100)
mn     <- tapply(tone$xpp$n1.r, list(tone$xpp$dB,tone$xpp$isiOctFine), mean )
dbax   <- tapply(tone$xpp$dB, list(tone$xpp$dB), mean )
isiax  <- tapply(tone$xpp$ISI, tone$xpp$isiOctFine, mean)

image(dbax, log(isiax[vi])/log(2),  mn[,vi])
contour(dbax, log(isiax[vi])/log(2),  mn[,vi], add=T)


plot(tapply(tone$xpp$n1.r, tone$xpp$isiOct,mean) ,type='l')
plot(tapply(tone$xpp$n1.r, tone$xpp$ampInd,mean) ,type='l')

## ==================================================================================
## How to normalize different days

# mean by day
lm0 <- lm(tone$xpp$n1.r ~ tone$xpp$exp)
mnByDay <- c(0,lm0$coeff[2:length(lm0$coeff) ])
plot(mnByDay,type='l')

# slope by day
lm1        <- lm(tone$xpp$n1.r ~ tone$xpp$exp + tone$xpp$exp:I(log(tone$xpp$ISI)))
slopeByDay <- lm1$coeff[(length(lm0$coeff)+1):length(lm1$coeff) ]

# robust variance by day
interQuartile <- function(x){   diff(quantile(x,pnorm(c(-.5,.5),m=0,sd=1)))   }
varByDay <- tapply(tone$xpp$n1.r, tone$xpp$exp, interQuartile)

plot(mnByDay, slopeByDay)
cor(mnByDay, slopeByDay)
cor(mnByDay, varByDay)
cor(slopeByDay, varByDay)

plot(mnByDay/abs(mean(mnByDay)),type='l')
points(slopeByDay/abs(mean(slopeByDay)), col='red',type='l')
points(varByDay/mean(varByDay), col='blue',type='l')

## ==================================================================================
## is the EEG amplitude of the current tone dependent on the SPL of previous tone?
DVstr <- 'n1.r'
disc <- subjFigure_CmpByAmpAndCnd(     tone,DVstr=DVstr,norm=F,eps=F)
disc <- subjFigure_CmpByDeltaAmp(      tone,DVstr=DVstr,norm=T,eps=F)
disc <- subjFigure_CmpByDeltaAmpAndAmp(tone,DVstr=DVstr,norm=T,eps=F,trialType=c(1))



lm0       <- lm(tone$xpp[[DVstr]] ~ as.factor(tone$xpp$exp) + as.factor(tone$xpp$ampInd) + I( log(tone$xpp$ISI)) + tone$xpp$trialType)
lm1       <- lm(tone$xpp[[DVstr]] ~ as.factor(tone$xpp$exp) + as.factor(tone$xpp$ampInd) + I( log(tone$xpp$ISI)) + tone$xpp$trialType + tone$xpp$deltaAmpInd)
anova(lm0,lm1)

n1nrm     <- lm0$residuals

tapply(n1nrm, tone$xpp$trialType, mean)
tapply(n1nrm, tone$xpp$dB, mean)

lm1       <- lm(n1nrm ~ 1)
lm2       <- lm(n1nrm ~ 1 + tone$xpp$deltaAmpInd)
anova(lm1,lm2)

dbax      <- tapply(tone$xpp$dB, tone$xpp$ampInd, mean)
deltax    <- tapply(tone$xpp$deltaAmpInd, tone$xpp$deltaAmpInd, mean)

N1        <- tapply(n1nrm, list(tone$xpp$trialType,tone$xpp$deltaAmpInd), mean )
plot(  deltax,N1[1,],col='blue',type='l', ylim=c(-1,1)*na.max(abs(N1)))
points(deltax,N1[2,],col='red', type='l')
abline(h=0); abline(v=0)


trialType <- c(1)
valInd    <- tone$xpp$trialType %in% trialType
tapply(n1nrm[valInd], tone$xpp$dB[valInd], mean)

N1        <- tapply(n1nrm[valInd], list(tone$xpp$dB[valInd],tone$xpp$deltaAmpInd[valInd]), mean )
zlim      <- c(-1,1) * na.max(abs(N1)) 
image(dbax, deltax, N1, col=bluered.colors(100), zlim=zlim)

lm1       <- lm(n1nrm[valInd] ~ 1)
lm2       <- lm(n1nrm[valInd] ~ 1 + tone$xpp$deltaAmpInd[valInd])
anova(lm1,lm2)




## =====================================
## =====================================
Off <- 0
lmbd <- .3
U <- .9
isi <- rep(1,15)
amp <- c(rep(0.2,5),rep(0.8,10))
tamp <- (amp+Off)/(5+Off)
scl  <- (2.5+Off)/(5+Off)
tst1 <- cumsumamp(isi,tamp, lmbd, U)
tst2 <- cumsumscl(isi,tamp/scl, lmbd, (2.5+Off)/(5+Off)*U)
plot(   cumsum(isi), tst1, ylim=c(0,1),type='b' )
points( cumsum(isi), tst2, col='red',type='b')



getSteadyState <- function(soa, db, lmbd=.5,U=.7,Off=0,show=T){
  
  isi <- rep(soa,50)
  amp <- rep(db,50)
  tamp <- (amp+Off)/(5+Off)
  scl  <- (2.5+Off)/(5+Off)
  tst1 <- cumsumamp(isi,tamp, lmbd, U)
  tst2 <- cumsumscl(isi,tamp/scl, lmbd, (2.5+Off)/(5+Off)*U)
  if(show){
    plot(   cumsum(isi), tst1, ylim=c(0,1),type='b' )
    points( cumsum(isi), tst2, col='red',type='b')
  }
  return(c(tst1[50],tst2[50]))
}
dBList  <- seq(1,1,to=5)
isiList <- 2^seq(-2,3,by=1) 

getDB <- function(dbl,soa=1,...){sapply(dbl, function(x){getSteadyState(soa,x,...)})}
getDB(dBList,soa=1,U=.1)

getDBSOA <- function(soal,dbl,...){
  tmp <- sapply(soal, function(x){getDB(dbl,soa=x,...) })
  return( list( amp=tmp[ seq(1,by=2,length=dim(tmp)[1]/2), ],scl=tmp[ seq(2,by=2,length=dim(tmp)[1]/2), ]) )
}
na.max( abs(tst$amp-tst$scl)/max(c(tst$amp,tst$scl)) )
tst <- getDBSOA(isiList,dBList,lmbd=.3,U=.7,show=F)

image(dBList, log(isiList)/log(2), tst$amp-tst$scl, zlim=c(-1.25,1.25)*na.max(abs(tst$amp-tst$scl)),col=bluered.colors(100))
image(dBList, log(isiList)/log(2), tst$amp, zlim=c(0,.6),col=heat.colors(100))
image(dBList, log(isiList)/log(2), tst$scl, zlim=c(0,.6),col=heat.colors(100))


isi <- c(rep(1,10),rep(1.5,9),10,rep(1.5,10))
amp <- c(rep(1.5,10),rep(1.5,10),rep(1.5,10))

lmbd <- .5
U    <- .9
Off  <- 0

tamp <- (amp+Off)/(5+Off)
scl  <- (2.5+Off)/(5+Off)
tst1 <- cumsumamp(isi,tamp, lmbd, U)
tst2 <- cumsumscl(isi,tamp/scl, lmbd, (2.5+Off)/(5+Off)*U)
plot(   cumsum(isi), tst1, ylim=c(0,1),type='b' )
points( cumsum(isi), tst2, col='red',type='b')

isi <- c(rep(.5,10),rep(.5,10),rep(.5,10))
amp <- c(rep(5,10),rep(3,10),rep(1,10))
Off <- 0
tamp <- (amp+Off)/(5+Off)
tst1 <- cumsumamp(isi,tamp, .4, .5)
tst2 <- cumsumscl(isi,tamp, .4, .5)
plot( cumsum(isi), tst1, ylim=c(0,1),type='b' )
points( cumsum(isi), tst2, col='red',type='b')


isi <- c(rep(.25,10),rep(.5,10),rep(.75,10))
amp <- c(rep(3,10),rep(3,10),rep(3,10))
Off <- 0
tamp <- (amp+Off)/(5+Off)
tst1 <- cumsumamp(isi,tamp, .4, .5)
tst2 <- cumsumscl(isi,tamp, .4, .5)
plot( cumsum(isi), tst1, ylim=c(0,1),type='b' )
points( cumsum(isi), tst2, col='red',type='b')


## ===================================================================
## playing with the predictions for a real data set

lmbd <- .3
U    <- .7
Off  <- 0
scl  <- (2.5+Off)/(5+Off)

cumsumV <- function(isi, amp, lmbd, U){
  ampN <- squeeze(amp, from.range=range(amp), to.range=c(1,5) )
  
  getPop <- function(thr){cumsumamp(isi, as.numeric(ampN>thr), lmbd, U )}
  tst    <- sapply(seq(0.5,by=1,4.5),getPop)  
  return( apply(tst,1,mean) )
}

getNonSteady<-function(tn,lmbd=.3,U=.7,Off=0){
  subs <- which(tn$xpp$trialType==1 & tone$xpp$toneType==1)
  isi <- tn$xpp$ISI[subs]
  amp <- (tn$xpp$ampInd[subs]+Off)/(5+Off)
  tst1 <- cumsumamp(isi, amp    , lmbd,                 1*U)
  tst2 <- cumsumscl(isi, amp/scl, lmbd, (2.5+Off)/(5+Off)*U)
  tst3 <- cumsumV(  isi, amp,     lmbd,                 1*U)
  
  print(cor(tst1,tst2))
  print(cor(tst1,tst3))
  #plot( cumsum(isi[1:20]), tst1[1:20], ylim=c(0,1),type='b' )
  #points( cumsum(isi[1:20]), tst2[1:20], col='red',type='b')
  N1amp <- tapply(tst1,list(tn$xpp$dB[subs],tn$xpp$isiOct[subs]),mean)
  N1scl <- tapply(tst2,list(tn$xpp$dB[subs],tn$xpp$isiOct[subs]),mean)
  N1V   <- tapply(tst3,list(tn$xpp$dB[subs],tn$xpp$isiOct[subs]),mean)
  
  return(list(amp=N1amp,scl=N1scl,vee=N1V,tst1=tst1,tst2=tst2,tst3=tst3)) 
}
tst <- getNonSteady(tone,lmbd=.4,U=.6,Off=0)
image(sort(unique(tone$xpp$dB)), seq(-2,by=1,to=4), tst$amp-tst$scl, zlim=c(-1.25,1.25)*na.max(abs(tst$amp-tst$scl)),col=bluered.colors(100))
na.max( abs(tst$amp-tst$scl)/max(c(tst$amp,tst$scl)) )

image(tst$amp, zlim=c(0,.6),col=heat.colors(100))
image(tst$scl, zlim=c(0,.6),col=heat.colors(100))
image(tst$vee, zlim=c(0,.6),col=heat.colors(100))


lm0 <- lm(tst$tst1 ~ tone$xpp$dB[subs]+I(log(tone$xpp$ISI[subs])/log(2)) )
image( tapply(lm0$resid,list(tone$xpp$dB[subs],tone$xpp$deltaAmpInd[subs]),mean), zlim=c(-.1,.1),col=heat.colors(100))
lm0 <- lm(tst$tst2 ~ tone$xpp$dB[subs]+I(log(tone$xpp$ISI[subs])/log(2)))
image( tapply(lm0$resid,list(tone$xpp$dB[subs],tone$xpp$deltaAmpInd[subs]),mean ), zlim=c(-.1,.1),col=heat.colors(100))
lm0 <- lm(tst$tst3 ~ tone$xpp$dB[subs]+I(log(tone$xpp$ISI[subs])/log(2)))
image( tapply(lm0$resid,list(tone$xpp$dB[subs],tone$xpp$deltaAmpInd[subs]),mean ), zlim=c(-.1,.1),col=heat.colors(100))

#plot(tst1,tst2)

plot(   tapply(tst1, tone$xpp$isiOct[subs],mean ),ylim=c(0,1), type='b')
points( tapply(tst2, tone$xpp$isiOct[subs],mean ), col='red', type='b')

plot(   tapply(tst1, tone$xpp$dB[subs],mean ),ylim=c(0,1), type='b')
points( tapply(tst2, tone$xpp$dB[subs],mean ), col='red', type='b')

plot(   rep(1.0,5), tapply(tst1, tone$xpp$dB[subs],mean ),ylim=c(0,1), type='b', xlim=c(0,5))
points( rep(1.2,5), tapply(tst2, tone$xpp$dB[subs],mean ), col='red', type='b')

points( rep(3.0,5), tapply(tst1, tone$xpp$isiOct[subs],mean )[1:5], col='black', type='b')
points( rep(3.2,5), tapply(tst2, tone$xpp$isiOct[subs],mean )[1:5], col='red', type='b')



##===================================================
parm <- c(0.5,.6)
ybin <- tone$xpp$p2p<450 & tone$xpp$powEOG<5000
disc <- runModelSTPSP(parm, tone, ybin, DVstr='z.n1.r', BYstr='exp',shuffel='no', 
                      returnModel=FALSE, show=T, eps=F,chNr=1,fileName='defaultFileName', FUN=cumsumscl)

## no scaling with amplitude at all
disc <- callSTPSP(tone,                 DVstr='n1.r', gd=T, show=TRUE, showCh=1, who='random', BYstr='exp',
                  FUN=cumsum1, minISI=.190, maxISI=12.8, startFromScratch=F)
## scale R/EEGamp with tone amplitude
disc <- callSTPSP(tone,                 DVstr='n1.r', gd=T, show=TRUE, showCh=1, who='random', BYstr='exp', 
                  FUN=cumsumscl, minISI=.190, maxISI=12.8, startFromScratch=F)
## scale expended Fraciton with tone amplitude
disc <- callSTPSP(tone,                 DVstr='p1a.r', gd=T, show=TRUE, showCh=1, who='random', BYstr='exp', 
                  FUN=cumsumamp, minISI=.190, maxISI=12.8, startFromScratch=F)
