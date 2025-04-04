# OVERVIEW all data collected for BRAIN grant ====
rm(list=ls())


# OVERVIEW ANIMALS AND IMPLANTS ====
#
# JESSE   EEG + SC96
# BEN     EEG + SC224
# CIRQUE  EEG + SC224
#
# SAM     EEG + VProbe
# FUZZY   EEG + VProbe
# WALTER  EEG + VProbe
#
# SOLEIL  EEG[+ Mesoscope]



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


colNames <- c('project','session','ElectrodeConfig','animal','year','month','link','date','time','dataSet')
dfEEGOnly <- rbind( df[which(df$ElectrodeConfig=='e'),colNames], dfEEG[,colNames])
table(dfEEGOnly$animal,    dfEEGOnly$ElectrodeConfig)
table(dfEEGOnly$animal,    dfEEGOnly$dataSet)



dfALL <- rbind( df[,colNames], dfEEG[, colNames])#, dfD[,colNames] )
max(table(dfALL$session))
which(table(dfALL$session)==2)
length(which(table(dfALL$session[which(dfALL$animal=='Sam')])==2))
length(which(table(dfALL$session[which(dfALL$animal=='Walter')])==2))


table(df$animal,    df$ElectrodeConfig)
table(dfEEG$animal, dfEEG$ElectrodeConfig)


table(df$animal,    df$dataSet)
table(dfEEG$animal, dfEEG$dataSet)


# includes all sessions regardless of quality
table(dfD$animal,   dfD$ElectrodeConfig)
table(dfD$animal,  dfD$numProbes)
tapply(dfD$numProbes, dfD$animal, sum)


# includes curated sessions (but all of Fuzzy)
# includes one line per probe; thus multiple lines per session are possible
table(dfVP$animal,  dfVP$ElectrodeConfig)
max(table(dfVP$session))





## grand average EEG plots by amplitude and ici ====
animal <- 'Fuzzy'


chStr <- paste('ch',1:33,sep='')
chStrIdx <- 1:33
dataSet <- 'g'
if (animal == 'Fuzzy'){
  chStrIdx <- c(1:19, 32, 33)
  dataSet <- 'a'
}


tdf <- dfEEGOnly[dfEEGOnly$animal==animal & dfEEGOnly$dataSet==dataSet,]


tn <- loadSignalProcessor( tdf, project=project, epoch='raw', chStr=chStr[chStrIdx], analysis='sessionID')
tn$channelPosition <- getEEGpos(animal)[chStrIdx,]

showSignalProcessor( tn, frml=~sessionID, trng=c(-25,250))
showSignalProcessor( tn, trng=c(-25,250), vlineat = c(21,45,85) )

p1ndx <- which( round(tn$taxis) == 21)
p1amp <- tn$lfp[,p1ndx,]
image(1:dim(p1amp)[1],1:dim(p1amp)[2],p1amp, col=jet.colors(100), xlab='sessionID',ylab='chNr')


plot(p1amp[,1], ylim=c(-20,130), type='l')
for (cx in 2:dim(p1amp)[2]){
  points(p1amp[,cx], col=cx, type='l')
}
plot(p1amp[,16])
tn$dfline[which(p1amp[,16]<30),]


tmp <- loadSignalProcessor( tn$dfline[38,], project=project, epoch='raw', chStr=chStr, analysis='sessionID')
tmp$channelPosition <- getEEGpos('Fuzzy')#[chStrIdx,]
showSignalProcessor( tmp )


ici <- loadSignalProcessor( tdf, project=project, epoch='llp', chStr=chStr[chStrIdx], analysis='iciInd')
ici$channelPosition <- getEEGpos(animal)[chStrIdx,]

amp <- loadSignalProcessor( tdf, project=project, epoch='llp', chStr=chStr[chStrIdx], analysis='ampInd')
amp$channelPosition <- getEEGpos(animal)[chStrIdx,]

showSignalProcessor( amp, frml=~ampInd, trng=c(-5,150) )
showSignalProcessor( ici, frml=~iciInd, trng=c(-5,150) )


## GRAND AVERAGE BSP ====
bspNdx <- 1:dim(tdf)[2]
if (animal == 'Fuzzy')
  bspNdx <- 25:dim(tdf)[2] # first 24 sessions were recorded with different sampling rate

if (animal == 'Ben')
  bspNdx <- c(3:6,8:dim(tdf)[2]) # first 2 sessions were recorded with different sampling rate


bsp <- loadSignalProcessor( tdf[bspNdx,], project=project, epoch='bsp', chStr=chStr, analysis='ampInd')
bsp$channelPosition <- getEEGpos(animal)

bspICI <- loadSignalProcessor( tdf[bspNdx,], project=project, epoch='bsp', chStr=chStr, analysis='iciInd')
bspICI$channelPosition <- getEEGpos(animal)

showSignalProcessor( bsp, frml=~ampInd, trng=c(0,15), vlineat = c(4,7,9), scl='numeric',sclfctr=10 )
showSignalProcessor( bspICI, frml=~iciInd, trng=c(0,20), vlineat = c(11), scl='numeric',sclfctr=10  )


showSignalProcessor( bsp, frml=~ampInd, trng=c(0,5), vlineat = c(4,7,9) )
showSignalProcessor( bspICI, frml=~iciInd, trng=c(0,5), vlineat = c(4,7,9) )


## grand average MUA VProbe plots by amplitude and ici ====
animal <- 'Fuzzy'

tdf   <- dfVP[dfVP$animal==animal & dfVP$ElectrodeConfig=='ev' & dfVP$dataSet=='g',]
#tdf <- tdf[-which(c(tdf$session)%in%c("Fuzzy_20230201_0919_zm19000","Fuzzy_20230221_0941_zm21800")),]



chStr       <- paste('ch',34:57,sep='')
tn1          <- loadSignalProcessor( tdf, project=project, epoch='mua', chStr=chStr, analysis='sessionID', chBlkN = 24, chan2var=T)
tn2          <- loadSignalProcessor( tdf, project=project, epoch='mua', chStr=chStr, analysis='sessionID', chBlkN = 24, chan2var=T)

tn$xpp$mx   <- apply(tn$lfp,c(1,3),max)
bltrngndx   <- which(tnn$taxis > -150 & tnn$taxis < 0)
tn$xpp$blsd <- apply(tn$lfp[,bltrngndx,],c(2), sd)

showSignalProcessor(tn, frml=~sessionID )
#tn$channelPosition <- data.frame(xpos=(1:dim(tn$lfp)[3]-1)%/%24,ypos=(1:dim(tn$lfp)[3]-1)%%24, zpos=rep(1,dim(tn$lfp)[3]) )


tnn   <- loadSignalProcessor( tdf, project=project, epoch='mua', chStr=chStr, analysis='ampInd_iciInd', chan2var=T)
tnnrm <- nrmByChanSignalProcessor( tnn )

table(tnnrm$xpp$nrm>2)/30
table(tnn$xpp$sessionID, tnnrm$xpp$nrm>1)/30

showSignalProcessor( tnnrm, frml=~ampInd, subsExpr = expression(nrm>1), trng=c(-10,250), scl='numeric',sclfctr = 1)
showSignalProcessor( tnnrm, frml=~iciInd, subsExpr = expression(nrm>1), trng=c(-10,250), scl='numeric',sclfctr = 1)


mn1 <- avgSignalProcessor( tnn1, frml=~sessionID )
mn2 <- avgSignalProcessor( tnn2, frml=~sessionID )

showSignalProcessor( tnn, frml=~ampInd, trng=c(-10,250) )
showSignalProcessor( tnn, frml=~Group.2, trng=c(-10,250) )



#tnn$channelPosition <- data.frame(xpos=(1:dim(tn$lfp)[3]-1)%/%24,ypos=(1:dim(tn$lfp)[3]-1)%%24, zpos=rep(1,dim(tn$lfp)[3]) )
#tn <- avgSignalProcessor( tnn )


# normalize tnn$lfp by standard deviation on the baseline
bltrngndx <- which(tnn$taxis > -150 & tnn$taxis < 0)
mx   <- apply(tn$lfp,c(1,3),max)
sdbl <- apply(tn$lfp[,bltrngndx,],c(2), sd)

tnnrm <- tnn
for (cx in 1:dim(tn$lfp)[3] ){
  tnnrm$lfp[,,cx] <- tnn$lfp[,,cx]/ sdbl[cx]
}

tnrm      <- avgSignalProcessor( tnnrm )
sdblnrm   <- apply(tnrm$lfp[,bltrngndx,],c(2), sd)
mxNrm     <- apply(tnrm$lfp,c(1,3),max)

plot(sort(mxNrm));abline(h=10)
plot(sort(mx))

showSignalProcessor( tnrm, scl='global',sclfctr=2, trng=c(-5,150) )


tnamp     <- avgSignalProcessor( tnnrm, frml=~ampInd )
#tnamp     <- reRefSignalProcessor( tnamp, dim=1, reRefChan=1:dim(tnamp$lfp)[1] )
tnamp$lfp <- apply(tnamp$lfp[,, which(mxNrm>50&mxNrm<5000) ],c(1,2),mean); dim(tnamp$lfp) <- c(dim(tnamp$lfp),1)

showSignalProcessor( tnamp, frml=~ampInd, trng=c(-5,150), vlineat=c(4,10,16))


tnici     <- avgSignalProcessor( tnnrm, frml=~Group.2 )
tnici     <- reRefSignalProcessor( tnici, dim=1, reRefChan=1:dim(tnici$lfp)[1] )
tnici$lfp <- apply(tnici$lfp[,,mx>1],c(1,2),mean); dim(tnici$lfp) <- c(dim(tnici$lfp),1)

tntn     <- tnn
tntn$lfp <- apply(tntn$lfp[,,mx>1],c(1,2),mean); dim(tntn$lfp) <- c(dim(tntn$lfp),1)


showSignalProcessor( tntn, frml=~Group.2, trng=c(-25,250) )
showSignalProcessor( tntn, frml=~ampInd, trng=c(-25,250) )

showSignalProcessorXYC( tntn, yaxis='ampInd',caxis='Group.2', trng=c(-25,250) )
showSignalProcessorXYC( tntn, caxis='ampInd',yaxis='Group.2', trng=c(-25,250) )



showSignalProcessor( tnamp, frml=~ampInd, trng=c(-5,150), vlineat=c(4,10,16))
showSignalProcessor( tnici, frml=~Group.2, trng=c(-5,150), vlineat=c(4,10,16))



plot(mn$taxis, apply(mn$lfp[,,mx>1,drop=F],2,mean), type='l' )

mntnn          <- tnn
mntnn$lfp      <- apply(tnn$lfp[,,mx>1],c(1,2), mean, drop=F)
dim(mntnn$lfp) <- c(dim(tnn$lfp)[c(1,2)],1)


mnmn <- avgSignalProcessor( tnn, frml=~ampInd+Group.2)

showSignalProcessorXYC( mnmn, yaxis='ampInd',caxis='Group.2', trng=c(-25,250) )
showSignalProcessorXYC( mntnn, caxis='ampInd',yaxis='Group.2', trng=c(-25,250) )


mntnnBF <- mntnn
for (fx in 1:dim(mntnn$xpp)[1] )
  mntnnBF$lfp[fx,,] <- mntnnBF$lfp[fx,,] - mnmn$lfp[1,,]


showSignalProcessorXYC( mntnnBF, yaxis='ampInd',caxis='Group.2', trng=c(-25,250) )
showSignalProcessorXYC( mntnnBF, caxis='ampInd',yaxis='Group.2', trng=c(-25,250) )



