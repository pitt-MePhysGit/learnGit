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


.GlobalEnv$pp        <- function(){paste(baseDir, '/ampOddClick/word/VProbe/plots/',sep='')}

source( paste(baseDir,'ampOddClick/r-files/functions/prepDFVprobe.R',sep='') )
source( paste(baseDir,'ampOddClick/r-files/functions/prep_df_mua.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/functions_sessionAverageVProbe.R',sep='') )

source( paste(baseDir,'ampOddClick/r-files/figures/call_AOC_timeseries_bySOA.R',sep='') )
source( paste(baseDir,'ampOddClick/r-files/figures/figure_AOC_timeseries_bySOA.R',sep='') )

source( paste(baseDir,'ampOddClick/r-files/figures/call_AOC_timeseries_byAMP.R',sep='') )
source( paste(baseDir,'ampOddClick/r-files/figures/figure_AOC_timeseries_byAMP.R',sep='') )

#source( paste(baseDir,'ampOddClick/r-files/functions/loadGASession.R',sep='') )

source( paste(baseDir,'RSkernel/r-files/functions/expand2chdf.R',sep='') )


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

## LOAD data frame      ====
project         <- 'ampOddClick'
linkFileName     <- 'ampOddclickVProbe_link.txt'   # only good sessions with VProbes
linkFileNameS     <- 'ampOddclick_link.txt'   # only good sessions with VProbes

dfV              <- prepDFsignalProcessor(project,linkFileName=linkFileName)
dfS              <- prepDFsignalProcessor(project,linkFileName=linkFileNameS)


dataSet <- 'g'
animal <- c('Fuzzy')#,'Walter')


df <- dfV # default assume that we are deadling with a VProbe animal
if (animal %in% c('Jesse','Ben','Cirque')){
  df <- dfS[dfS$ElectrodeConfig %in% c('escy','essuy'),]
  df$toneArtefact <- F
  df$onsetArtefact <- F
  df$drug <- 'N'
  
  df$numChan <- 96
  if (df$ElectrodeConfig == 'essuy')
    df$numChan <- 224
  
  df$spacing <- 1
  df$Layer   <- 0
  
  
}

cleanOnly   <- FALSE
#badDays     <- c(12, 14, 16, 18, 22, 24)

df$use      <- c(1)
if (cleanOnly)
  df$use[ which(df$toneArtefact+df$onsetArtefact>0) ] <- 0
#df$use[badDays] <- 0

useInd        <- which( df$dataSet==dataSet & df$animal%in%animal & df$drug=='N' & df$use )# & df$probeString=='t')

tdf           <- df[useInd,]
chdf          <- expand2chdf(tdf)

## set folder for plots ====
animalStr <- '_all'
if (length(animal)==1)
  animalStr <- paste('_',animal,sep='')

cleanStr  <- ''
if (cleanOnly)
  cleanStr <- '_clean'

uIdir <- paste('g', animalStr,cleanStr,sep='')
.GlobalEnv$pp   <- function(){paste(baseDir, '/ampOddClick/word/VProbe/plots/',uIdir,'/',sep='')}
if (!file.exists(pp())){
  system( paste('mkdir', pp()) )
  system( paste('mkdir ', pp(), 'mnMUA/', sep=''))
  system( paste('mkdir ', pp(), 'mnAvRec/', sep=''))
  system( paste('mkdir ', pp(), 'mnHdAvRec/', sep=''))
}

## CHECK if DATA exists ====
if (1==0){
  ch34mua <- loadSignalProcessor(tdf, epoch='mua', analysis='ampInd_iciInd', chStr='ch34', blrng=NULL  )
  ch34aov <- loadSignalProcessor(tdf, epoch='mua', analysis='trAOV_ampInd_iciInd_trialType_toneType_trialNr', chStr='ch34', blrng=NULL  )
  
  aov1 <- loadSignalProcessor(tdf, epoch='mua', analysis='trAOV_parm_amp12', chStr='ch34', blrng=NULL  )
  aov3 <- loadSignalProcessor(tdf, epoch='mua', analysis='trAOV_parm_amp23', chStr='ch34', blrng=NULL  )
  aov5 <- loadSignalProcessor(tdf, epoch='mua', analysis='trAOV_parm_amp45', chStr='ch34', blrng=NULL  )
  
  showSignalProcessor(aov1,subsExpr = expression(effect=='ampInd') , showSE=T, scl='numeric',sclfctr = .2)
  showSignalProcessor(aov3,subsExpr = expression(effect=='ampInd') , showSE=T, scl='numeric',sclfctr = .2)
  showSignalProcessor(aov5,subsExpr = expression(effect=='ampInd') , showSE=T, scl='numeric',sclfctr = .2)
}






## LOAD MUA               ====
if (startFromScratch){
  pvl <- loadSignalProcessor(chdf, epoch='mua', analysis='trAOV_ampInd_iciInd_trialType_toneType_trialNr', chStr=NA, blrng=NULL  )
  bet <- loadSignalProcessor(chdf, epoch='mua', analysis='trAOV_parm_ampInd_iciInd_trialType_toneType_trialNr', chStr=NA, blrng=NULL  )
  tn  <- loadSignalProcessor(chdf, epoch='mua', analysis='ampInd_iciInd', chStr=NA  )
  
  # still missing some sessions for the analysis by amp
  bet_amp1 <- loadSignalProcessor(chdf, epoch='mua', analysis='trAOV_parm_amp1', chStr=NA, blrng=NULL  )
  bet_amp2 <- loadSignalProcessor(chdf, epoch='mua', analysis='trAOV_parm_amp2', chStr=NA, blrng=NULL  )
  bet_amp3 <- loadSignalProcessor(chdf, epoch='mua', analysis='trAOV_parm_amp3', chStr=NA, blrng=NULL  )
  bet_amp4 <- loadSignalProcessor(chdf, epoch='mua', analysis='trAOV_parm_amp4', chStr=NA, blrng=NULL  )
  bet_amp5 <- loadSignalProcessor(chdf, epoch='mua', analysis='trAOV_parm_amp5', chStr=NA, blrng=NULL  )
  
  bet_amp12 <- loadSignalProcessor(chdf, epoch='mua', analysis='trAOV_parm_amp12', chStr=NA, blrng=NULL  )
  bet_amp23 <- loadSignalProcessor(chdf, epoch='mua', analysis='trAOV_parm_amp23', chStr=NA, blrng=NULL  )
  bet_amp34 <- loadSignalProcessor(chdf, epoch='mua', analysis='trAOV_parm_amp34', chStr=NA, blrng=NULL  )
  bet_amp45 <- loadSignalProcessor(chdf, epoch='mua', analysis='trAOV_parm_amp45', chStr=NA, blrng=NULL  )
  
  
  
  
  chdf$chIDstr       <- paste(chdf$session,chdf$chNr,sep='_')
  tn$dfline$chIDstr  <- paste(tn$dfline$session,tn$dfline$chNr,sep='_')
  pvl$dfline$chIDstr <- paste(pvl$dfline$session,pvl$dfline$chNr,sep='_')
  bet$dfline$chIDstr <- paste(bet$dfline$session,bet$dfline$chNr,sep='_')
  
  chdf$sessionID       <- 1:dim(chdf)[1]
  tn$dfline$sessionID  <- 1:dim(tn$dfline)[1]
  pvl$dfline$sessionID <- 1:dim(pvl$dfline)[1]
  bet$dfline$sessionID <- 1:dim(bet$dfline)[1]
  
  length(unique(chdf$session))
  length(unique(chdf$chIDstr))
  length(unique(chdf$sessionID))
  
  # check if all sessions have the 5*6 traces
  mn   <- avgSignalProcessor( tn, frml=~sessionID )
  valN <- mn$xpp$Ntrials==30
  
  #mn$dfline[which(mn$xpp$Ntrials<30),]
  #valtn <- mn$xpp$sessionID[which(mn$xpp$Ntrials==30)]
  #table( mn$dfline$session[valtn] )
  #table( mn$dfline$session[-valtn] )
  
  
  pnbool <- chdf$chIDstr %in% pvl$dfline$chIDstr
  tnbool <- chdf$chIDstr %in% tn$dfline$chIDstr
  bebool <- chdf$chIDstr %in% bet$dfline$chIDstr
  
  table(chdf$session[ !pnbool ])
  table(chdf$session[ !tnbool ])
  table(chdf$session[ !bebool ])
  
  if (sum(!valN)>0)
    stop()
  
  #table(chdf$session[ !valN ])
  
  allbool <- pnbool & tnbool & bebool# & valN
  table(allbool)
  
  ## subset tn
  tns        <- tn
  valID      <- tn$dfline$chIDstr %in% chdf$chIDstr[allbool]
  tns$dfline <- tn$dfline[valID,]
  tns <- subs.data(tn, tn$xpp$sessionID %in% tns$dfline$sessionID)
  tns$dfline <- tn$dfline[valID,]
  #min(table(tns$xpp$sessionID))
  
  
  ## subset pn
  pvls       <- pvl
  valID      <- pvl$dfline$chIDstr %in% chdf$chIDstr[allbool]
  pvls$dfline <- pvl$dfline[valID,]
  pvls <- subs.data(pvl, pvl$xpp$sessionID %in% pvls$dfline$sessionID)
  pvls$dfline <- pvl$dfline[valID,]
  
  ## subset bet
  bets       <- bet
  valID      <- bet$dfline$chIDstr %in% chdf$chIDstr[allbool]
  bets$dfline <- bet$dfline[valID,]
  bets <- subs.data(bet, bet$xpp$sessionID %in% bets$dfline$sessionID)
  bets$dfline <- bet$dfline[valID,]
  
  dim(pvls$dfline)
  dim(tns$dfline)
  dim(bets$dfline)
}

## make figures           ====
.GlobalEnv$pp        <- function(){paste(baseDir, '/ampOddClick/word/VProbe/plots/',uIdir,'/mnMUA/',sep='')}

thrsh1 <- c(1e-50,1e-10,1e-10) # default MUA values
thrsh2 <- c(1e-50,1e-15,1e-8)  # modified MUA values to match significance rate of Amp and ICI
thrsh3 <- c(1e-20,1e-5,1e-5)   # default avRec values

onsetDistributionFigure(pvls, thrsh=thrsh1)
onsetDistributionFigure(pvls, thrsh=thrsh2)

tnsnrm <- nrmByChanSignalProcessor(tns, frml=~sessionID)

tst    <- getOnsetMain(pvls, thrsh3)
valNdx <- which(tst$OnLat<20 & tst$AmpLat<20 & tst$ICILat<20) 


showSignalProcessor( tns,    frml=~ampInd, trng=c(-25,150), bty = 'l', yshft=-0.25, scl='numeric',sclfctr = 2, butterfly=F, showSE=T )
showSignalProcessor( tns,    frml=~ampInd, trng=c(-25,150), bty = 'l', yshft=-0.25, scl='numeric',sclfctr = 3, butterfly=F, showSE=T, subsExpr = expression(sessionID%in%valNdx) )
showSignalProcessor( tnsnrm, frml=~ampInd, trng=c(-25,150), bty = 'l', yshft=-0.25, scl='numeric',sclfctr = 2,   butterfly=F, showSE=T, subsExpr = expression(sessionID%in%valNdx) )

#tst <- loadSignalProcessor( tdf[4,],epoch='mua', analysis='trAOV_parm_ampInd_iciInd_trialType_toneType_trialNr', chStr='ch34', blrng=NULL)
#tst$xpp

eps.file( 'muanrm_byAmp.eps',pdf=T,s=3)
showSignalProcessor( tnsnrm, frml=~ampInd, trng=c(-25,150), bty = 'l', yshft=-0.25, scl='numeric',sclfctr = 1.5,   butterfly=F, showSE=T, subsExpr = expression(sessionID%in%valNdx) )
dev.off()

eps.file( 'muanrm_byICI.eps',pdf=T,s=3)
showSignalProcessor( tnsnrm, frml=~iciInd, trng=c(-25,150), bty = 'l', yshft=-0.25, scl='numeric',sclfctr = 1.5,   butterfly=F, showSE=T, subsExpr = expression(sessionID%in%valNdx), colFUN=soa.colors )
dev.off()

eps.file( 'muanrm_byAmp_zoom.eps',pdf=T,s=3)
showSignalProcessor( tnsnrm, frml=~ampInd, trng=c(-5,50), bty = 'l', yshft=-0.25, scl='numeric',sclfctr = 1.5,   butterfly=F, showSE=T, subsExpr = expression(sessionID%in%valNdx) )
dev.off()

eps.file( 'muanrm_byICI_zoom.eps',pdf=T,s=3)
showSignalProcessor( tnsnrm, frml=~iciInd, trng=c(-5,50), bty = 'l', yshft=-0.25, scl='numeric',sclfctr = 1.5,   butterfly=F, showSE=T, subsExpr = expression(sessionID%in%valNdx), colFUN=soa.colors )
dev.off()

eps.file( 'muanrm_byAmp_butterfly.eps',pdf=T,s=3)
showSignalProcessor( tnsnrm, frml=~ampInd, trng=c(-5,50), bty = 'l', yshft=-0.0, scl='numeric', sclfctr = .5, butterfly=T, showSE=T, subsExpr = expression(sessionID%in%valNdx))#, vlineat=c(mean(tst$AmpLat[valNdx])) )
dev.off()

eps.file( 'muanrm_byICI_butterfly.eps',pdf=T,s=3)
showSignalProcessor( tnsnrm, frml=~iciInd, trng=c(-5,50), bty = 'l', yshft=-0.0, scl='numeric', sclfctr = .5, butterfly=T, showSE=T, subsExpr = expression(sessionID%in%valNdx), colFUN=soa.colors)#, vlineat=c(mean(tst$ICILat[valNdx])) )
dev.off()



# betas
Int <- subs.data(bets, bets$xpp$effect=='(Intercept)')
Int$xpp$nrm <- apply( Int$lfp, c(1,3), function(x){max(abs(x))} )
betsnrm <- bets
betsnrm$xpp$nrm <- Int$xpp$nrm[betsnrm$xpp$sessionID]
for (cx in 1:dim(betsnrm$lfp)[1])
  betsnrm$lfp[cx,,]   <- bets$lfp[cx,,] / betsnrm$xpp$nrm[cx]


eps.file( 'beta_Amp_ICI.eps',pdf=T,s=3)
showSignalProcessor( bets, frml=~effect, trng=c(-5,50),scl='numeric', sclfctr=.5, yshft=-.25, subsExpr=expression(effect%in%c('ampInd','iciInd') & sessionID%in%valNdx), bty='l',showSE=T )
dev.off()

eps.file( 'betanrm_Amp_ICI.eps',pdf=T,s=3)
showSignalProcessor( betsnrm, frml=~effect, trng=c(-5,50),scl='numeric', sclfctr=.25, yshft=-.25, subsExpr=expression(effect%in%c('ampInd','iciInd') & sessionID%in%valNdx), bty='l',showSE=T )
dev.off()



# beta ICI for different amplitudes ====
sbet_amp1 <- subs.data(bet_amp1, bet_amp1$xpp$effect=='iciInd')
sbet_amp1$xpp$effect <- 'ici_amp1'

sbet_amp2 <- subs.data(bet_amp2, bet_amp2$xpp$effect=='iciInd')
sbet_amp2$xpp$effect <- 'ici_amp2'

sbet_amp3 <- subs.data(bet_amp3, bet_amp3$xpp$effect=='iciInd')
sbet_amp3$xpp$effect <- 'ici_amp3'

sbet_amp4 <- subs.data(bet_amp4, bet_amp4$xpp$effect=='iciInd')
sbet_amp4$xpp$effect <- 'ici_amp4'

sbet_amp5 <- subs.data(bet_amp5, bet_amp5$xpp$effect=='iciInd')
sbet_amp5$xpp$effect <- 'ici_amp5'


sbet_byAmp <- append.data( append.data( append.data(sbet_amp1, sbet_amp2), append.data(sbet_amp3,sbet_amp4)), sbet_amp5)
eps.file( 'beta_ICIByAmp.eps',pdf=T,s=3)
showSignalProcessor( sbet_byAmp, frml=~effect, trng=c(-5,50),scl='numeric', sclfctr=.5, yshft=-.25, subsExpr=expression( sessionID%in%valNdx), bty='l',showSE=T )
dev.off()
#showSignalProcessor( sbet_byAmp, frml=~effect, trng=c(-5,50),scl='numeric', sclfctr=.5, yshft=-.25, bty='l',showSE=T )

## beta amp for different pairs of amplitude ====
sbet_amp12 <- subs.data(bet_amp12, bet_amp12$xpp$effect=='ampInd')
sbet_amp12$xpp$effect <- 'amp12'

sbet_amp23 <- subs.data(bet_amp23, bet_amp23$xpp$effect=='ampInd')
sbet_amp23$xpp$effect <- 'amp23'

sbet_amp34 <- subs.data(bet_amp34, bet_amp34$xpp$effect=='ampInd')
sbet_amp34$xpp$effect <- 'amp34'

sbet_amp45 <- subs.data(bet_amp45, bet_amp45$xpp$effect=='ampInd')
sbet_amp45$xpp$effect <- 'amp45'

sbet_AmpbyAmp <- append.data( append.data(sbet_amp12, sbet_amp23), append.data(sbet_amp34,sbet_amp45))
eps.file( 'beta_AmpByAmp.eps',pdf=T,s=3)
showSignalProcessor( sbet_AmpbyAmp, frml=~effect, trng=c(-5,50),scl='numeric', sclfctr=.5, yshft=-.25, subsExpr=expression( sessionID%in%valNdx), bty='l',showSE=T )
dev.off()

sbet_ICIamp12 <- subs.data(bet_amp12, bet_amp12$xpp$effect=='iciInd')
sbet_ICIamp12$xpp$effect <- 'ICIamp12'

sbet_ICIamp23 <- subs.data(bet_amp23, bet_amp23$xpp$effect=='iciInd')
sbet_ICIamp23$xpp$effect <- 'ICIamp23'

sbet_ICIamp34 <- subs.data(bet_amp34, bet_amp34$xpp$effect=='iciInd')
sbet_ICIamp34$xpp$effect <- 'ICIamp34'

sbet_ICIamp45 <- subs.data(bet_amp45, bet_amp45$xpp$effect=='iciInd')
sbet_ICIamp45$xpp$effect <- 'ICIamp45'

sbet_ICIvsAmp_12 <- append.data(sbet_amp12, sbet_amp2)
eps.file( 'beta_AmpvsICI_12.eps',pdf=T,s=3)
showSignalProcessor( sbet_ICIvsAmp_12, frml=~effect, trng=c(-5,50),scl='numeric', sclfctr=.5, yshft=-.25, subsExpr=expression( sessionID%in%valNdx), bty='l',showSE=T )
dev.off()

sbet_ICIvsAmp_23 <- append.data(sbet_amp23, sbet_amp3)
eps.file( 'beta_AmpvsICI_23.eps',pdf=T,s=3)
showSignalProcessor( sbet_ICIvsAmp_23, frml=~effect, trng=c(-5,50),scl='numeric', sclfctr=.5, yshft=-.25, subsExpr=expression( sessionID%in%valNdx), bty='l',showSE=T )
dev.off()

sbet_ICIvsAmp_34 <- append.data(sbet_amp34, sbet_amp4)
eps.file( 'beta_AmpvsICI_34.eps',pdf=T,s=3)
showSignalProcessor( sbet_ICIvsAmp_34, frml=~effect, trng=c(-5,50),scl='numeric', sclfctr=.5, yshft=-.25, subsExpr=expression( sessionID%in%valNdx), bty='l',showSE=T )
dev.off()

sbet_ICIvsAmp_45 <- append.data(sbet_amp45, sbet_amp5)
eps.file( 'beta_AmpvsICI_45.eps',pdf=T,s=3)
showSignalProcessor( sbet_ICIvsAmp_45, frml=~effect, trng=c(-5,50),scl='numeric', sclfctr=.5, yshft=-.25, subsExpr=expression( sessionID%in%valNdx), bty='l',showSE=T )
dev.off()



## SAVE GRAND AVERAGE MUA ====
tnsnrm$epoch <- 'mua'
tnsnrm$exp   <- uIdir
gamua <- avgSignalProcessor( tnsnrm, frml=~ampInd+iciInd, subsExpr = expression(sessionID%in%valNdx), saveData=T  )


## LOAD AvRec                ====
.GlobalEnv$pp        <- function(){paste(baseDir, '/ampOddClick/word/VProbe/plots/',uIdir,'/mnAvRec/',sep='')}

## rerunning last fuzzy vprobe recording -A to fix the problem with the iciInd


tdf$chNr <- tdf$startVIndex

tstr  <- loadSignalProcessor(tdf[3,], epoch='mua', chStr='ch34')#, analysis='ampInd_iciInd', chStr=NA  )
tsta  <- loadSignalProcessor(tdf[4,], epoch='avRec', chStr='ch34', analysis='ampInd_iciInd')

table(tstr$xpp$iciInd)
table(tsta$xpp$iciInd)

dsc <- checkISISignalProcessor( tst ,show=T)
table(tst$xpp$ampInd,tst$xpp$iciInd)

#avr  <- loadSignalProcessor(tdf, epoch='avRec', analysis='ampInd_iciInd', chStr=NA  )
pavr <- loadSignalProcessor(tdf, epoch='avRec', analysis='trAOV_ampInd_iciInd_trialType_toneType_trialNr', chStr=NA, blrng=NULL  )
beta <- loadSignalProcessor(tdf, epoch='avRec', analysis='trAOV_parm_ampInd_iciInd_trialType_toneType_trialNr', chStr=NA, blrng=NULL  )
avrx  <- loadSignalProcessor(pavr$dfline, epoch='avRec', analysis='ampInd_iciInd', chStr=NA  )
avrx$xpp$iciInd <- as.numeric(avrx$xpp$Group.2)

avrnrm <- nrmByChanSignalProcessor(avrx, frml=~sessionID)


onsetDistributionFigure(pavr, thrsh=thrsh1)
onsetDistributionFigure(pavr, thrsh=thrsh2)
onsetDistributionFigure(pavr, thrsh=thrsh3)



tst    <- getOnsetMain(pavr, thrsh3)
valNdx <- which(tst$OnLat<50 & tst$AmpLat<50 & tst$ICILat<50) 



eps.file( 'avRec_byAmp.eps',pdf=T,s=3)
showSignalProcessor( avrnrm, frml=~ampInd, trng=c(-5,50), bty = 'l', scl='numeric', yshft=-0.25,sclfctr = 1.5,   butterfly=F, showSE=T , subsExpr = expression(sessionID%in%valNdx & iciInd<6))
dev.off()

eps.file( 'avRec_byICI.eps',pdf=T,s=3)
showSignalProcessor( avrnrm, frml=~iciInd, trng=c(-5,50), bty = 'l', scl='numeric', yshft=-0.25, sclfctr = 1.5,   butterfly=F, showSE=T, colFUN=soa.colors2, subsExpr = expression(sessionID%in%valNdx & iciInd<6) )
dev.off()

eps.file( 'avRec_byAmp_butterfly.eps',pdf=T,s=3)
showSignalProcessor(avrnrm, frml=~ampInd, trng=c(-5,50), bty = 'l', yshft=-0.0, scl='numeric', sclfctr = .45, butterfly=T, showSE=T, subsExpr = expression(sessionID%in%valNdx & iciInd<6))#, vlineat=c(mean(tst$AmpLat[valNdx])) )
dev.off()

eps.file( 'avRec_byICI_butterfly.eps',pdf=T,s=3)
showSignalProcessor( avrnrm, frml=~iciInd, trng=c(-5,50), bty = 'l', yshft=-0.0, scl='numeric', sclfctr = .45, butterfly=T, showSE=T, colFUN=soa.colors2, subsExpr = expression(sessionID%in%valNdx & iciInd<6) )#, vlineat=c(mean(tst$ICILat[valNdx])) )
dev.off()


# betas
bets <- beta
Int <- subs.data(bets, bets$xpp$effect=='(Intercept)')
Int$xpp$nrm <- apply( Int$lfp, c(1,3), function(x){max(abs(x))} )
betsnrm <- bets
betsnrm$xpp$nrm <- Int$xpp$nrm[betsnrm$xpp$sessionID]
for (cx in 1:dim(betsnrm$lfp)[1])
  betsnrm$lfp[cx,,]   <- bets$lfp[cx,,] / betsnrm$xpp$nrm[cx]


eps.file( 'beta_Amp_ICI.eps',pdf=T,s=3)
showSignalProcessor( bets, frml=~effect, trng=c(-5,50),scl='numeric', sclfctr=.5, yshft=-.25, subsExpr=expression(effect%in%c('ampInd','iciInd') & sessionID%in%valNdx), bty='l',showSE=T )
dev.off()

eps.file( 'betanrm_Amp_ICI.eps',pdf=T,s=3)
showSignalProcessor( betsnrm, frml=~effect, trng=c(-5,50),scl='numeric', sclfctr=.25, yshft=-.25, subsExpr=expression(effect%in%c('ampInd','iciInd') & sessionID%in%valNdx), bty='l',showSE=T )
dev.off()



#showSignalProcessor( avrnrm, frml=~ampInd, trng=c(-25,150), bty = 'l', yshft=-0.3, scl='numeric',sclfctr = 1.5, butterfly=F, showSE=T)
#showSignalProcessor( avrnrm, frml=~Group.2, trng=c(-25,150), bty = 'l', yshft=-0.3, scl='numeric',sclfctr = 1.5, butterfly=F, showSE=T)

##showSignalProcessor( avrx, frml=~ampInd, trng=c(-25,150), bty = 'l', yshft=-0.0, scl='numeric', sclfctr = .3, butterfly=T, showSE=T)
#showSignalProcessor( tns, frml=~iciInd, trng=c(-25,150), bty = 'l', yshft=-0.0, scl='numeric', sclfctr = .3, butterfly=T, showSE=T)




## LOAD hdAvRec               ====
.GlobalEnv$pp        <- function(){paste(baseDir, '/ampOddClick/word/VProbe/plots/',uIdir,'/mnHdAvRec/',sep='')}

## rerunning last fuzzy vprobe recording -A to fix the problem with the iciInd
## rerunning fuzzy 2023 03 23, both sessions. One had a wrong link file

wrongSamplingRate <- c()
if (animal == 'Fuzzy')
  wrongSamplingRate <- c(40,41,42,43,45,46,47,48,50,51,52,53,54)
valnd <- which(! 1:dim(tdf)[1] %in% wrongSamplingRate)

avr  <- loadSignalProcessor(tdf[valnd,], epoch='hdAvRec', analysis='ampInd_iciInd', chStr=NA  )
pavr <- loadSignalProcessor(tdf[valnd,], epoch='hdAvRec', analysis='trAOV_ampInd_iciInd_trialType_toneType_trialNr', chStr=NA, blrng=NULL  )
avrnrm <- nrmByChanSignalProcessor(avr, frml=~sessionID)

#thrsh3 <- c(1e-20,1e-5,1e-5)
#onsetDistributionFigure(pavr, thrsh=thrsh1, mxLat=30)
#onsetDistributionFigure(pavr, thrsh=thrsh2, mxLat=30)
onsetDistributionFigure(pavr, thrsh=thrsh3, mxLat=30)



tst    <- getOnsetMain(pavr, thrsh3)
valNdx <- which(tst$OnLat<20 & tst$AmpLat<20 & tst$ICILat<20) 



eps.file( 'hdAvRec_byAmp.eps',pdf=T,s=3)
showSignalProcessor( avrnrm, frml=~ampInd, trng=c(-5,25), bty = 'l', scl='numeric', yshft=-0.25,sclfctr = 1.5,   butterfly=F, showSE=T , subsExpr = expression(sessionID%in%valNdx & iciInd<6))
dev.off()

eps.file( 'hdAvRec_byICI.eps',pdf=T,s=3)
showSignalProcessor( avrnrm, frml=~iciInd, trng=c(-5,25), bty = 'l', scl='numeric', yshft=-0.25, sclfctr = 1.5,   butterfly=F, showSE=T, colFUN=soa.colors2, subsExpr = expression(sessionID%in%valNdx & iciInd<6) )
dev.off()

eps.file( 'hdAvRec_byAmp_butterfly.eps',pdf=T,s=3)
showSignalProcessor(avrnrm, frml=~ampInd, trng=c(-5,25), bty = 'l', yshft=-0.0, scl='numeric', sclfctr = .45, butterfly=T, showSE=T, subsExpr = expression(sessionID%in%valNdx & iciInd<6))#, vlineat=c(mean(tst$AmpLat[valNdx])) )
dev.off()

eps.file( 'hdAvRec_byICI_butterfly.eps',pdf=T,s=3)
showSignalProcessor( avrnrm, frml=~iciInd, trng=c(-5,25), bty = 'l', yshft=-0.0, scl='numeric', sclfctr = .45, butterfly=T, showSE=T, colFUN=soa.colors2, subsExpr = expression(sessionID%in%valNdx & iciInd<6) )#, vlineat=c(mean(tst$ICILat[valNdx])) )
dev.off()




showSignalProcessor( avr,    frml=~ampInd, trng=c(-25,150), bty = 'l', yshft=-0.25, scl='numeric',sclfctr = .25, butterfly=F, showSE=T )
showSignalProcessor( avr,    frml=~ampInd, trng=c(-25,150), bty = 'l', yshft=-0.25, scl='numeric',sclfctr = .25, butterfly=F, showSE=T, subsExpr = expression(sessionID%in%valNdx) )
showSignalProcessor( avrnrm, frml=~ampInd, trng=c(-25,150), bty = 'l', yshft=-0.25, scl='numeric',sclfctr = 2,   butterfly=F, showSE=T )

showSignalProcessor( avrnrm, frml=~ampInd, trng=c(-25,150), bty = 'l', yshft=-0.25, scl='numeric',sclfctr = 2,   butterfly=F, showSE=T, subsExpr = expression(sessionID%in%valNdx) )
showSignalProcessor( avrnrm, frml=~iciInd, trng=c(-25,150), bty = 'l', yshft=-0.25, scl='numeric',sclfctr = 2,   butterfly=F, showSE=T, subsExpr = expression(sessionID%in%valNdx) )


showSignalProcessor( avrnrm, frml=~ampInd, trng=c(-25,150), bty = 'l', yshft=-0.0, scl='numeric', sclfctr = .5, butterfly=T, showSE=T, subsExpr = expression(sessionID%in%valNdx), vlineat=c(mean(tst$AmpLat[valNdx])) )
showSignalProcessor( avrnrm, frml=~iciInd, trng=c(-25,150), bty = 'l', yshft=-0.0, scl='numeric', sclfctr = .5, butterfly=T, showSE=T, subsExpr = expression(sessionID%in%valNdx), vlineat=c(mean(tst$ICILat[valNdx])) )














## =============================== PLOTS FOR BRAIN INITIATIVE MEETING 2019 =============================
rngY <- c(-0.2,1)
rngT <- c(-10,250)
call_AOCtimeseries_bySOA(what='mua', who='all', dpthStr='MUA',  rngT=rngT, rngY=rngY, ganrm=T)
call_AOCtimeseries_byAMP(what='mua', who='all', dpthStr='MUA',  rngT=rngT, rngY=rngY, ganrm=T)

rngY <- c(-0.5,0.5)
call_AOCtimeseries_bySOA(what='csd', who='all', dpthStr='SGSi', rngT=rngT, rngY=rngY, ganrm=T)
call_AOCtimeseries_byAMP(what='csd', who='all', dpthStr='SGSi', rngT=rngT, rngY=rngY, ganrm=T)


## ================================ difference plots
rngY <- c(-0.30,0.30)
rngT <- c(-10,50)
call_AOCtimeseries_bySOA(what='mua', who='all', dpthStr='MUA',  rngT=rngT, rngY=rngY, ganrm=T, dxdn=3)
call_AOCtimeseries_byAMP(what='mua', who='all', dpthStr='MUA',  rngT=rngT, rngY=rngY, ganrm=T, dxdn=3)

call_AOCtimeseries_bySOA(what='csd', who='all', dpthStr='SGSi',  rngT=rngT, rngY=rngY, ganrm=T, dxdn=3)
call_AOCtimeseries_byAMP(what='csd', who='all', dpthStr='SGSi',  rngT=rngT, rngY=rngY, ganrm=T, dxdn=3)




if (1==0){
## =============================================== other stuff
call_AOCtimeseries_bySOA(what='mua',who='all', dpthStr='L2', rngT=c(-10,50),ganrm=F)
call_AOCtimeseries_bySOA(what='mua',who='all', dpthStr='L3', rngT=c(-10,50),ganrm=F)
call_AOCtimeseries_bySOA(what='mua',who='all', dpthStr='L4', rngT=c(-10,50),ganrm=F)
call_AOCtimeseries_bySOA(what='mua',who='all', dpthStr='L5', rngT=c(-10,50),ganrm=F)

call_AOCtimeseries_bySOA(what='raw',who='all', dpthStr='SGSo', rngT=rngT,ganrm=F)


call_AOCtimeseries_bySOA(what='csd',who='all', dpthStr='SGSo', rngT=rngT,ganrm=T)
call_AOCtimeseries_bySOA(what='csd',who='all', dpthStr='GSi' , rngT=rngT,ganrm=T)


# by amplitude
call_AOCtimeseries_byAMP(what='mua',who='all', dpthStr='MUA', rngT=rngT,ganrm=T)

#tst <- figure_AOC_timeseries(tone)













# ------------------------------- test recording to 
#tst   <- loadAllSessions( useInd, extn='mean_ISIAmp',who=who,what=what,xnd=1, ynd=1)
#plot(apply(tst,1,mean), type='l' )

Nchan  <- dim(chdf)[1]   #dim(tst)[2]
SOAndx <- 1:6
AMPndx <- 1:5
dBVal  <- seq(62, by=6, length=5)
Ncnd   <- length(SOAndx) * length(AMPndx) 

for (sx in SOAndx){
  for (ax in AMPndx){
    tmp <- loadAllSessions( useInd, extn='mean_ISIAmp',who='all',what='mua',xnd=sx, ynd=ax)
    if (sx+ax == 2)   
      lfp <- tmp
    if(sx+ax > 2)
      lfp <- cbind(lfp,tmp)
  }
}

# ================================================= select valid channels and create xpp file
depthrng <- c(-1500,0)
chdf$use <- chdf$depth >= depthrng[1] & chdf$depth <= depthrng[2]

tmpxpp     <- expand.grid(isiOct=SOAndx, ampInd=AMPndx)
tmpxpp$ISI <- tmpxpp$isiOct
tmpxpp$dB  <- dBVal[tmpxpp$ampInd]

xpp     <- tmpxpp[rep(1:dim(tmpxpp)[1], each=Nchan), ]
xpp$use <- rep(chdf$use,Ncnd)

# ===================================== create ephys data frame
tone       <- list()
tone$xpp   <- xpp[xpp$use,]
tone$lfp   <- t(lfp[,xpp$use]); dim(tone$lfp) <- c(dim(tone$lfp),1)
tone$taxis <- seq(-150, by=1, length=900)


eps.file(paste('AOC_timeseries_',what,'_',who,'.eps',sep=''), pdf=T,s=4)
disc <- figure_AOC_timeseries(tone, whoStr='lfp', rngT=c(-10,100))
dev.off()

}






