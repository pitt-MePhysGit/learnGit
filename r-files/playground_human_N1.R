#
# model human data
#
library(lattice)

source( paste('~/Dropbox/r-compile/rUtils/R/rUtils.R',sep='') )

source( paste('~/Dropbox/r-compile/mav.test/R/mav.test.R',sep='') )

project   <- 'ampOddClick'

pp <- function(){return(paste('~/Dropbox/',project,'/plots/human/',sep=''))}

fileNameL <- 'Amp_singl_trials_N1.txt'

LLR <- read.table( paste('~/Dropbox/',project,'/RData/',fileNameL,sep=''), header=T,sep = "," )

table(LLR$subjID)
table(LLR$subjID,LLR$sessionNr)
table(LLR$sessionNr, LLR$blockNr,LLR$subjID)

table(LLR$clickIntensity)

# create ISI
LLR$deltaTime                  <- c(20,diff(LLR$absTime))
LLR$deltaTime[LLR$deltaTime<0] <- 20


# 

# ======================= Intensity
plot( tapply(LLR$EEG_N1,      LLR$clickIntensity, na.mean), type='b', pch=20)

## BIN-BASED ANALYSIS ====
brks <- 2^seq(-2,by=1,to=3)
LLR$ISIbin1 <- cut((LLR$ISI%<=%8)%>=%.25,breaks=brks)

brks <- 2^seq(-2,by=.5,to=3)
LLR$ISIbin2 <- cut((LLR$ISI%<=%8)%>=%.25,breaks=brks)

brks <- 2^seq(-2,by=1/3,to=3)
LLR$ISIbin3 <- cut((LLR$ISI%<=%8)%>=%.25,breaks=brks)

brks <- 2^seq(-2,by=1/4,to=3)
LLR$ISIbin4 <- cut((LLR$ISI%<=%8)%>=%.25,breaks=brks)

brks <- 2^seq(-2,by=1/6,to=3)
LLR$ISIbin6 <- cut((LLR$ISI%<=%8)%>=%.25,breaks=brks)


table(LLR$ISIbin1,LLR$clickIntensity)
table(LLR$ISIbin2,LLR$clickIntensity)
table(LLR$ISIbin3,LLR$clickIntensity)
table(LLR$ISIbin4,LLR$clickIntensity)
table(LLR$ISIbin6,LLR$clickIntensity)

plot(tapply(LLR$ISI,      LLR$ISIbin1, na.mean), tapply(LLR$EEG_N1,      LLR$ISIbin1, na.mean), type='b', pch=20)
plot(tapply(LLR$ISI,      LLR$ISIbin2, na.mean), tapply(LLR$EEG_N1,      LLR$ISIbin2, na.mean), type='b', pch=20)
plot(tapply(LLR$ISI,      LLR$ISIbin3, na.mean), tapply(LLR$EEG_N1,      LLR$ISIbin3, na.mean), type='b', pch=20)
plot(tapply(LLR$ISI,      LLR$ISIbin4, na.mean), tapply(LLR$EEG_N1,      LLR$ISIbin4, na.mean), type='b', pch=20)
plot(tapply(LLR$ISI,      LLR$ISIbin6, na.mean), tapply(LLR$EEG_N1,      LLR$ISIbin6, na.mean), type='b', pch=20)

plot(  log2(tapply(LLR$ISI,      LLR$ISIbin1, na.mean)), tapply(LLR$EEG_N1,      LLR$ISIbin1, na.mean), type='b', pch=20, lwd=2,ylim=c(-2,-.5), xlim=c(-2,3))
points(log2(tapply(LLR$ISI,      LLR$ISIbin2, na.mean)), tapply(LLR$EEG_N1,      LLR$ISIbin2, na.mean), type='b', pch=20, lwd=2, col='blue')
points(log2(tapply(LLR$ISI,      LLR$ISIbin3, na.mean)), tapply(LLR$EEG_N1,      LLR$ISIbin3, na.mean), type='b', pch=20, lwd=2,col='green')
points(log2(tapply(LLR$ISI,      LLR$ISIbin4, na.mean)), tapply(LLR$EEG_N1,      LLR$ISIbin4, na.mean), type='b', pch=20, lwd=2, col='orange')
points(log2(tapply(LLR$ISI,      LLR$ISIbin6, na.mean)), tapply(LLR$EEG_N1,      LLR$ISIbin6, na.mean), type='b', pch=20, lwd=2, col='red')
points(log2(tapply(LLR$ISI,      LLR$ISIbin1, na.mean)), tapply(LLR$EEG_N1,      LLR$ISIbin1, na.mean), type='b', pch=20, lwd=2)



showN1facilitation <- function(lr, DVstr='EEG_N1', IVstr='ISIbin1', nrm=F, nrmStr='byAmp',...){
  mn    <- tapply(lr[[DVstr]],      list(lr[[IVstr]], lr$clickIntensity, lr$subjID), na.mean)
  isiax <- tapply(lr$ISI,      list(lr[[IVstr]]), na.mean)
  
  mn1    <- tapply(lr[[DVstr]],      list(lr$ISIbin1, lr$clickIntensity, lr$subjID), na.mean)
  
  
  ylim <- c(-3.25,0)
  if (nrm){
    ylim <- c(-1.25,.125)
    
    if (nrmStr=='byAmp5'){
      for (sx in 1:dim(mn)[3])
        mn[,,sx] <- mn[,,sx]/abs(mn1[dim(mn1)[1],3,sx])
    }
    
    if (nrmStr=='byAmp'){
      for (sx in 1:dim(mn)[3])
        for (ax in 1:dim(mn)[2])
          mn[,ax,sx] <- -1 * mn[,ax,sx]/mn1[dim(mn1)[1],ax,sx]
    }
    
    if (nrmStr=='X'){
      ylim <- c(-2.25,0.25)
      for (sx in 1:dim(mn)[3])
          mn[,,sx] <-  -1 * mn[,,sx]/mean(mn1[,,sx])
    }
  }
  
  errorplot(    log2( isiax), apply(mn,  c(1,2), mean)[,1], ey=apply(mn,  c(1,2), na.se)[,1], type='b', pch=20, lwd=2,ylim=ylim, xlim=c(-2,3), xaxt='n',xlab='ISI', ylab='N1', ... )
  axis(1, at=seq(-2,by=1,to=3), labels = c('1/4','1/2','1','2','4','8') )
  errorplot(    log2( isiax), apply(mn,  c(1,2), mean)[,2], ey=apply(mn,  c(1,2), na.se)[,2], type='b', pch=20, lwd=2, add=T, col='blue')
  errorplot(    log2( isiax), apply(mn,  c(1,2), mean)[,3], ey=apply(mn,  c(1,2), na.se)[,3], type='b', pch=20, lwd=2,add=T, col='red')
  abline(v=seq(-2,by=1,to=3), lwd=.5, col='gray')
}

LR <- subset(LLR, !is.nan(LLR$EEG_N1))
maxN1    <- tapply(LR$EEG_N1,      list(LR$ISIbin1, LR$clickIntensity, LR$subjID), na.mean)[5,3,]

lr <- LR
eps.file('N1_facilitation.eps',pdf=T,ratio = 1/2, s=6)
par( mfcol=c(1,4))
showN1facilitation(lr, IVstr='ISIbin1', main='1 octave bins')
showN1facilitation(lr, IVstr='ISIbin2', main='1/2 octave bins')
showN1facilitation(lr, IVstr='ISIbin3', main='1/3 octave bins')
showN1facilitation(lr, IVstr='ISIbin4', main='1/4 octave bins')
dev.off()

eps.file('N1_facilitation_normByAmp5.eps',pdf=T,ratio = 1/2, s=6)
par( mfcol=c(1,4))
showN1facilitation(lr, IVstr='ISIbin1', nrm=T, nrmStr='byAmp5', main='1 octave bins')
showN1facilitation(lr, IVstr='ISIbin2', nrm=T, nrmStr='byAmp5', main='1/2 octave bins')
showN1facilitation(lr, IVstr='ISIbin3', nrm=T, nrmStr='byAmp5', main='1/3 octave bins')
showN1facilitation(lr, IVstr='ISIbin4', nrm=T, nrmStr='byAmp5', main='1/4 octave bins')
#showN1facilitation(lr, IVstr='ISIbin6', nrm=T, nrmStr='byAmp5', main='1/6 octave bins')
dev.off()

# only subjects with very large N1s
lr     <- subset(LR, LR$subjID %in% names(maxN1)[maxN1 < -3] )
eps.file('N1_facilitation_LargeN1s_normByAmp5.eps',pdf=T,ratio = 1/2, s=6)
par( mfcol=c(1,4))
showN1facilitation(lr, IVstr='ISIbin1', nrm=T, nrmStr='byAmp5', main='1 octave bins')
showN1facilitation(lr, IVstr='ISIbin2', nrm=T, nrmStr='byAmp5', main='1/2 octave bins')
showN1facilitation(lr, IVstr='ISIbin3', nrm=T, nrmStr='byAmp5', main='1/3 octave bins')
showN1facilitation(lr, IVstr='ISIbin4', nrm=T, nrmStr='byAmp5', main='1/4 octave bins')
#showN1facilitation(lr, IVstr='ISIbin6', nrm=T, nrmStr='byAmp5', main='1/6 octave bins')
dev.off()

# only subjects with very small N1s
lr     <- subset(LR, LR$subjID %in% names(maxN1)[maxN1 >= -3] )
eps.file('N1_facilitation_smallN1s_normByAmp5.eps',pdf=T,ratio = 1/2, s=6)
par( mfcol=c(1,4))
showN1facilitation(lr, IVstr='ISIbin1', nrm=T, nrmStr='byAmp5', main='1 octave bins')
showN1facilitation(lr, IVstr='ISIbin2', nrm=T, nrmStr='byAmp5', main='1/2 octave bins')
showN1facilitation(lr, IVstr='ISIbin3', nrm=T, nrmStr='byAmp5', main='1/3 octave bins')
showN1facilitation(lr, IVstr='ISIbin4', nrm=T, nrmStr='byAmp5', main='1/4 octave bins')
#showN1facilitation(lr, IVstr='ISIbin6', nrm=T, nrmStr='byAmp5', main='1/6 octave bins')
dev.off()


## MOVING AVERAGES   ====

lr1  <- subset( LLR, !is.nan(LLR$EEG_N1) & !is.na(LLR$ISI) & LLR$clickIntensity==1 )
lr2  <- subset( LLR, !is.nan(LLR$EEG_N1) & !is.na(LLR$ISI) & LLR$clickIntensity==2 )
lr3  <- subset( LLR, !is.nan(LLR$EEG_N1) & !is.na(LLR$ISI) & LLR$clickIntensity==3 )

doit <- function(lr){res <- get.mav( lr$EEG_N1, log2(lr$ISI), eval.at=seq(-1.9,by=.1,to=2.9), width=.5)} 
mavLst1 <- by(lr1, lr1$subjID, doit)
mavLst2 <- by(lr2, lr2$subjID, doit)
mavLst3 <- by(lr3, lr3$subjID, doit)


extractFcn <- function(sx, lst=mavLst1, whatStr='mav'){ return(lst[[sx]][[whatStr]])  }

mavGA1 <- sapply( 1:length(mavLst1), extractFcn,  whatStr='mav', lst=mavLst1)
mavGA2 <- sapply( 1:length(mavLst2), extractFcn,  whatStr='mav', lst=mavLst2)
mavGA3 <- sapply( 1:length(mavLst3), extractFcn,  whatStr='mav', lst=mavLst3)



eps.file('N1_facilitation_mav_05.eps',pdf=T,ratio = 2, s=4)
par(mfcol=c(1,1))
errorplot( mavLst1[[1]]$eval.at, apply(mavGA1,1,mean), ey=apply(mavGA1,1,na.se), type='l', xlim=c(-2,3), ylim=c(-3.25,0), xaxt='n', xlab='ISI', ylab='N1' )
axis(1, at=seq(-2,by=1,to=3), labels = c('1/4','1/2','1','2','4','8') )

errorplot( mavLst1[[1]]$eval.at, apply(mavGA2,1,mean), ey=apply(mavGA2,1,na.se), type='l', add=T, col='blue' )
errorplot( mavLst1[[1]]$eval.at, apply(mavGA3,1,mean), ey=apply(mavGA3,1,na.se), type='l', add=T, col='red' )

abline( v=mavLst1[[1]]$eval.at[ which( apply(mavGA1,1,mean) == max(apply(mavGA1,1,mean)) )], lty=2 )
abline( v=mavLst1[[1]]$eval.at[ which( apply(mavGA2,1,mean) == max(apply(mavGA2,1,mean)) )]-.01, col='blue',lty=3 )
abline( v=mavLst1[[1]]$eval.at[ which( apply(mavGA3,1,mean) == max(apply(mavGA3,1,mean)) )]+.01, col='red', lty=4 )

dev.off()



## single-subject analysis of facilitation

extractFacPar <- function(sx, lst=mavLst1){
  res    <- rep(NA,4)
  res[1] <- which(lst[[sx]]$mav==max(lst[[sx]]$mav) )
  res[2] <- lst[[sx]]$mav[ res[1] ]
  res[3] <- lst[[sx]]$mav[ res[1] ]
  if (res[1]>1)
    res[3] <- mean(lst[[sx]]$mav[1: min(res[1],which(lst[[sx]]$eval.at< -1.5))])  #[ lst[[sx]]$eval.at< -1.5 ])
  
  res[4] <- mean(lst[[sx]]$mav[ lst[[sx]]$eval.at> 2 ])
  
  return(res)  
  }



facParLst1 <- sapply(1:length(mavLst1), extractFacPar, lst=mavLst1)
facParLst2 <- sapply(1:length(mavLst1), extractFacPar, lst=mavLst2)
facParLst3 <- sapply(1:length(mavLst1), extractFacPar, lst=mavLst3)


N1par3        <- data.frame( t(facParLst3) )
names(N1par3) <- c('latPeakBin','N1peak','N1short','N1long')

N1par3$latPeakLog       <- mavLst3[[1]]$eval.at[ N1par3$latPeakBin ]
N1par3$latPeakMS        <- 2^N1par3$latPeakLog

N1par3$N1facilitation   <- N1par3$N1peak - N1par3$N1short
N1par3$N1refractoriness <- N1par3$N1peak - N1par3$N1long

N1par3$Group <- 3 
N1par3$Group[ which( N1par3$latPeakLog < -1.8 ) ]                                          <- 1
N1par3$Group[ which( N1par3$latPeakLog > -1.8  &  N1par3$N1peak - N1par3$N1short < 0.2)  ] <- 2
N1par3$Group[ which( N1par3$latPeakLog > 0.5 ) ]                                           <- 4

N1par3$SUBJID <- names(mavLst3)



eps.file('N1_facilitation_parameters_amp3_groups.eps',pdf=T,ratio = 1, s=4)
plot(   N1par3$latPeakLog, N1par3$N1facilitation, pch=20, col=N1par3$Group, cex=2, ylim=c(0,2), xlim=c(-2,3), xaxt='n', xlab='Latency of peak N1 amplitude [sec]', ylab='Amount of N1 facilitation [uV]' )
text( N1par3$latPeakLog, N1par3$N1facilitation, N1par3$SUBJID, adj=c(-.3) )
axis(1, at=seq(-2,by=1,to=3), labels = c('1/4','1/2','1','2','4','8') )
dev.off()


eps.file('N1_refractoriness_parameters_amp3_groups.eps',pdf=T,ratio = 1, s=4)
plot(   N1par3$latPeakLog, N1par3$N1refractoriness, pch=20, col=N1par3$Group, cex=2, ylim=c(0,3), xlim=c(-2,3), xaxt='n', xlab='Latency of peak N1 amplitude [sec]', ylab='Amount of N1 facilitation [uV]' )
text( N1par3$latPeakLog, N1par3$N1refractoriness, N1par3$SUBJID, adj=c(-.3) )
axis(1, at=seq(-2,by=1,to=3), labels = c('1/4','1/2','1','2','4','8') )
dev.off()


write.table( t(facParLst3[cNr,]), file= paste(pp(),'test.txt',sep=''), quote=F, col.names = cNms[cNr], row.names=F )

densityplot( N1par3$latPeakLog, bw=.1 )


par(mfrow=c(5,6))
for (sx in 1:length(tstlst))
  mav.show(mavLst1[[ sort(facParLst1[1,], index.return=T)$ix[sx]  ]], col=c('black','gray','gray'))

par(mfrow=c(5,6))
for (sx in 1:length(tstlst))
  mav.show(mavLst2[[ sort(facParLst2[1,], index.return=T)$ix[sx] ]], col=c('black','gray','gray'))


eps.file('N1_facilitation_mav_bySubject_amp3.eps',pdf=T,ratio = 1, s=8)
par(mfrow=c(5,6))
for (sx in 1:length(tstlst))
  mav.show(mavLst3[[ sort(facParLst3[1,], index.return=T)$ix[sx] ]], col=c('black','gray','gray'),xlab='ISI [log2(sec)]', ylab='N1 amplitude', main=names(mavLst3)[sort(facParLst3[1,], index.return=T)$ix[sx] ])
dev.off()




## CORRELATION WITH COGNITIVE TESTS ====

fileNameCog <- 'Cognitive_tests_N1_facilitation.txt'
cog         <- read.table( paste('~/Dropbox/',project,'/RData/',fileNameCog,sep=''), header=T,sep = "\t" )
identical( as.character(cog$SUBJ ), as.character(N1par3$SUBJID) )


cog <- data.frame( N1par3, cog )

N1VARNames  <- c('N1facilitation','latPeakLog','N1short','N1peak', 'N1long','N1refractoriness')
COGVARNames <- c('VOCAB_TS','OVERALLTSCR','MATRIX_TS','FULL2IQ','SPEEDTSCR','ATT_VIGTSCR','WMTSCR','VERBTSCR','VISTSCR','RPSTSCR','SOCCOGTSCR')

cog[,COGVARNames]
tst <- cog[which(!is.na(cog$VERBTSCR)),c(N1VARNames,COGVARNames)]
tst$N1peak <- -1 * tst$N1peak
tst$N1short <- -1 * tst$N1short
tst$N1long <- -1 * tst$N1long

par( mar=c(6,8,1,1) )
image(1:dim(tst)[2], 1:dim(tst)[2], cor(tst), zlim=c(-1,1), col=potential.colors(100), xaxt='n',yaxt='n', ylab='',xlab='' )
axis(1, 1:dim(tst)[2], names(tst), las=2 )
axis(2, 1:dim(tst)[2], names(tst), las=2 )
abline(h=length(N1VARNames)+.5 )
abline(v=length(N1VARNames)+.5 )

#library("psych")
resGrid <- expand.grid( N1VARNames=c(N1VARNames,COGVARNames), COGVARNames=c(N1VARNames,COGVARNames)) 
resGrid$p.value <- sapply( 1:dim(resGrid)[1], function(ix){ cor.test( tst[[ as.character(resGrid$N1VARNames[ix]) ]], tst[[ as.character( resGrid$COGVARNames[ix])]], method='pearson' )$p.value } )

# add markder for significant effecrs
pthrsh  <- c(.1, .05, .01, .001)
sigStr  <- c('','.','*','**','***')

crmat <- cor(tst)
crmat[which(pvalmat>pthrsh[1])] <- 0

xaxm <- rep( 1:dim(crmat)[1],dim(crmat)[1]  )
yaxm <- rep( 1:dim(crmat)[1],each=dim(crmat)[1]  )

pthrshMat <- 1 + (pvalmat<pthrsh[1]) + (pvalmat<pthrsh[2]) + (pvalmat<pthrsh[3]) + (pvalmat<pthrsh[4])

eps.file( 'CorrelationMatrix_N1_COG.eps', pdf=T )
par( mar=c(8,8,1,1) )
image(1:dim(tst)[2], 1:dim(tst)[2], crmat, zlim=c(-1,1), col=potential.colors(100), xaxt='n',yaxt='n', ylab='',xlab='' )
axis(1, 1:dim(tst)[2], names(tst), las=2 )
axis(2, 1:dim(tst)[2], names(tst), las=2 )
abline(h=length(N1VARNames)+.5,lwd=3 )
abline(v=length(N1VARNames)+.5,lwd=3 )
abline(h=2.5 )
abline(v=2.5 )
text( xaxm, yaxm, sigStr[pthrshMat] )
dev.off()


par( mar=c(1,1,1,1) )
plot(tst, pch=20, cex=.5)

#cog$GroupFac <- as.factor(cog$Group)
#cog$N1FF     <- as.factor(cog$Group==3)

cog$smallestN1LatMSlog <- log2( cog$smallestN1LatMS )


# Latency can not reliably be measured in group 2, because they mostly asymptote, thus providing a pretty wide window...
# Dergee of facilitation can be measured in all groups
# 
# create z-scored versions of the predictor variables

#cog$latLogZ <- cog$smallestN1LatMSlog
#cog$latLogZ[ cog$Group==2 ] <- NA

cog$N1f_latLogZ <- (cog$smallestN1LatMSlog - na.mean(cog$smallestN1LatMSlog))/na.sd(cog$smallestN1LatMSlog)
#cog$latLogZ[ is.na(cog$latLogZ) ] <- 0


cog$N1f_ampZ <- (cog$N1facilitationMV - na.mean(cog$N1facilitationMV))/na.sd(cog$N1facilitationMV)

facParLst3[3,]


DVstr <- 'VOCAB_TS'

checkDV <- function(DVstr, subs=cog$Group %in% c(1,2,3,4) ){
  cog$DV <- cog[[DVstr]]
  
  cogS <- subset(cog, !is.na(cog$VOCAB_TS) & subs, drop.levels=T) 
  
  lm0  <- lm( DV ~ 1,                   data=cogS )
  lm1a <- lm( DV ~ 1 + GroupFac,        data=cogS )
  lm1b <- lm( DV ~ 1 + N1FF,            data=cogS )
  lm1c <- lm( DV ~ 1 + latLogZ,         data=cogS )
  lm1d <- lm( DV ~ 1 + ampZ,            data=cogS )
  lm2  <- lm( DV ~ 1 + latLogZ + ampZ,  data=cogS )
  
  print(anova( lm0, lm1a))
  print(anova( lm0, lm1b))
  print(anova( lm0, lm1c))
  print(anova( lm0, lm1d))
  print(anova( lm0, lm2))
  
  par(mfrow=c(2,2))
  boxplot( DV ~ N1FF, data=cogS)
  
  boxplot( DV ~ Group, data=cogS, at=sort(unique(cogS$Group)) )
  points(jitter(cogS$Group), cogS$DV, col=cogS$GroupFac, pch=20)
  
  plot(cogS$latLogZ, cogS$DV, pch=20, col=cogS$Group)
  plot(     cogS$ampZ, cogS$DV, pch=20, col=cogS$Group)
}

checkDV( 'VOCAB_TS')                          # latency significant
checkDV( 'VOCAB_TS', cog$Group %in% c(1,2,3)) # still marginal

checkDV( 'OVERALLTSCR')
checkDV( 'MATRIX_TS')
checkDV( 'FULL2IQ')
checkDV( 'SPEEDTSCR')
checkDV( 'ATT_VIGTSCR')
checkDV( 'WMTSCR')
checkDV( 'VERBTSCR')
checkDV( 'VISTSCR')
checkDV( 'RPSTSCR')
checkDV( 'SOCCOGTSCR') # significant N1FacAmp
checkDV( 'SOCCOGTSCR', cog$Group %in% c(1,2,3))# still significant
checkDV( 'SOCCOGTSCR', cog$Group %in% c(1,3))# still significant


checkDV( 'ROLECURR')
checkDV( 'SOCIALCURR')

#sort( levels(cog$SUBJ), index.return=T)
#sort( c(facParLst3[6,]), index.return=T )
#identical( levels(cog$SUBJ)[c(cog$SUBJ)], c(facParLst3[6,]) )

#cog$smallestN1LatMS     <- facParLst3[7,]
#cog$N1facilitationMV    <- facParLst3[8,]
#cog$N1facilitationGroup <- facParLst3[5,]







## STUFF    ====


par( mfcol=c(1,4))
showN1facilitation(lr, IVstr='ISIbin1', nrm=T, nrmStr='byAmp')
showN1facilitation(lr, IVstr='ISIbin2', nrm=T, nrmStr='byAmp')
showN1facilitation(lr, IVstr='ISIbin3', nrm=T, nrmStr='byAmp')
showN1facilitation(lr, IVstr='ISIbin4', nrm=T, nrmStr='byAmp')

par( mfcol=c(1,4))
showN1facilitation(lr, IVstr='ISIbin1', nrm=T, nrmStr='X')
showN1facilitation(lr, IVstr='ISIbin2', nrm=T, nrmStr='X')
showN1facilitation(lr, IVstr='ISIbin3', nrm=T, nrmStr='X')
showN1facilitation(lr, IVstr='ISIbin4', nrm=T, nrmStr='X')



SOA_AMP_plot <- function(DV, SOA, AMP, ylim=NA, DVstr='EEG_N1',...){
  mnN1SOA <- tapply(DV, SOA, na.mean)
  seN1SOA <- tapply(DV, SOA, na.se)
  soax    <- tapply(SOA,SOA, mean)
  
  mnN1Amp <- tapply(DV, AMP, na.mean)
  seN1Amp <- tapply(DV, AMP, na.se)
  ampx    <- 1+tapply(AMP,AMP, mean)
  
  if (is.na(ylim))
    ylim    <- range( c(mnN1SOA,mnN1Amp) ) + c(-.25,.25)*diff(range( c(mnN1SOA,mnN1Amp) ))
  
  errorplot(1:length(soax),mnN1SOA, ey=seN1SOA, pch=20, col='black', ylim=ylim,...)
  errorplot(ampx,mnN1Amp, ey=seN1Amp, pch=20, add=T, col='red' )
}


SOA_AMP_plot(LLR$EEG_N1, LLR$ISIbin1, LLR$clickIntensity)


## ====================== normalize components

NORM_BY_SUBJ <- function(DV, subjID, type='SUBTRACT'){
  subjFct <- as.factor(subjID)
  lm1 <- lm(DV ~ subjFct)
  subjMean <- tapply(DV, subjID, na.mean)
  
  DVnrm    <- lm1$residuals + na.mean(subjMean) 
  return(DVnrm)
}

getAsymptote <- function(lr, DVstr='EEG_N1', nrm=F){
  ndx <- which(lr$ISI>6)
  mn <- mean(lr[[DVstr]][ndx])
  #print(lr$subjID[1])
  #print(mn)
  if (nrm)
    return( lr[[DVstr]]/mn )
  
  if (!nrm)
    return( lr[[DVstr]] )
  
}# get Asymptote as average activity for ISI>=4


LR$ix <- 1:dim(LR)[1]

tstn1 <- by( LR, list(LR$subjID,LR$clickIntensity), getAsymptote, nrm=T )
tstix <- by( LR, list(LR$subjID,LR$clickIntensity), getAsymptote, nrm=F, DVstr='ix' )

new <- LR[LR$ix[unlist(tstix)],]
new$EEG_N1_asy <- unlist(tstn1)

  
  




LR$EEG_N1_nrm <- NORM_BY_SUBJ(LR$EEG_N1,LR$subjID)

errorplot(    log2(tapply(LLR$ISI,      LLR$ISIbin1, na.mean)), tapply(LLR$EEG_N1,      LLR$ISIbin1, na.mean), ey=tapply(LLR$EEG_N1,      LLR$ISIbin1, na.se), type='b', pch=20, lwd=2,ylim=c(-2,-.5), xlim=c(-2,3))

errorplot(    log2(tapply(LR$ISI,       LR$ISIbin4, na.mean)) , tapply(LR$EEG_N1,      LR$ISIbin4, na.mean), ey=tapply(LR$EEG_N4,      LR$ISIbin1, na.se) , type='b', pch=20, lwd=2,ylim=c(-2.25,-.75), xlim=c(-2,3))
abline(v=seq(-2,by=1,to=3), lwd=.5, lty=2)

errorplot(    log2(tapply(new$ISI,      new$ISIbin1, na.mean)) , tapply(new$EEG_N1_asy,      new$ISIbin1, na.mean), ey=tapply(new$EEG_N1_asy,      new$ISIbin1, na.se) , type='b', pch=20, lwd=2,ylim=c(0,1.25), xlim=c(-2,3))
errorplot(    log2(tapply(new$ISI,      new$ISIbin4, na.mean)) , tapply(new$EEG_N1_asy,      new$ISIbin4, na.mean), ey=tapply(new$EEG_N1_asy,      new$ISIbin4, na.se) , type='b', pch=20, lwd=2,ylim=c(0,1.25), xlim=c(-2,3))
abline(v=seq(-2,by=1,to=3), lwd=.5, lty=2)





errorplot(    log2(tapply(new$ISI,      new$ISIbin4, na.mean)) , tapply(new$EEG_N1_asy,      new$ISIbin4, na.mean), ey=tapply(new$EEG_N1_asy,      new$ISIbin4, na.se) , type='b', pch=20, lwd=2, xlim=c(-2,3))


