## ======================================
## master_compareLinearModels_ketamine

baseDir <- '~/Dropbox/'
source(paste(baseDir, 'ampOdd_click/r-files/source_all_functions.R',sep=''))

## preliminaries
pp <- function(){paste(baseDir, 'ampOdd_click/word/plots/',sep='')}
rda <- function(){'~/Dropbox/ampOdd_click/Rdata/'}

subjList <- c('Sam','Walter')
dataSetList <- c('k')



## ===================================================================
## figure 
who <- 'all'

showCmp <- function(n, eps=T, clr='red', dS='l',DVstr='mnTDAEP'){
  
  loadN1 <- function(sbj, dS='l',who='all', DVstr='mnTDAEP'){ 
    load(file=paste(rda(),sbj,'_',dS,'_lmSummary_',who,'.rda',sep=''))
    return(resList[[n]][[DVstr]])
  }
  
  tst      <- loadN1('Sam',dS=dS, DVstr=DVstr)
  N1matL   <- array( unlist(lapply(subjList, loadN1, dS=dS, who=who, DVstr=DVstr)), c(dim(tst),length(subjList)) )
  mnN1matL <- apply(N1matL, c(1,2), mean) 
  seN1matL <- apply(N1matL, c(1,2), na.se) 
  
  xax <- 1:5
  if(eps)
    eps.file(paste('POPfig_',cmpStrList[n],'_',who,'_',dS,'_',DVstr ,'.eps',sep=''), pdf=T, s=1.65, ratio=1.5)
  
  errorplot(xax, mnN1matL[,1], ey=seN1matL[,1], type='p', pch=20, col=rep('black',length(xax)), ylim=c(-0.5,1.2),wd=0.1, main=cmpStrList[[n]] )
  errorplot(xax, mnN1matL[,2], ey=seN1matL[,2], type='p', pch=20, col=rep(clr,length(xax)),  add=T,wd=0.1 )
  lines(xax, mnN1matL[,1],col='black')
  lines(xax, mnN1matL[,2],col=clr)
  abline(h=c(0,1))
  
  if(eps)
    dev.off()
  
  return(N1matL)
}

TDmat <- lapply( 1:length(cmpStrList), showCmp, dS='k',DVstr='mnTDAEP',clr='red',eps=T )
LDmat <- lapply( 1:length(cmpStrList), showCmp, dS='k',DVstr='mnLDAEP',clr='red',eps=T )




## extracting the single-subject fits using the one-slope per session approach
loadN1 <- function(n, sbj='Sam', dS='k',who='all', DVstr='slopeTD'){ 
  load(file=paste(rda(),sbj,'_',dS,'_lmSummary_',who,'.rda',sep=''))
  return(resList[[n]][[DVstr]]$`Pr(>F)`[1])
}

TS <- unlist( lapply(1:length(cmpStrList), loadN1, DVstr='slopeTDSal',sbj='Sam' ))
TW <- unlist( lapply(1:length(cmpStrList), loadN1, DVstr='slopeTDSal',sbj='Walter' ))

DS <- unlist( lapply(1:length(cmpStrList), loadN1, DVstr='slopeLDSal',sbj='Sam' ))
DW <- unlist( lapply(1:length(cmpStrList), loadN1, DVstr='slopeLDSal',sbj='Walter' ))

## interactions with drug
TDS <- unlist( lapply(1:length(cmpStrList), loadN1, DVstr='slopeTD',sbj='Sam' ))
TDW <- unlist( lapply(1:length(cmpStrList), loadN1, DVstr='slopeTD',sbj='Walter' ))

LDS <- unlist( lapply(1:length(cmpStrList), loadN1, DVstr='slopeLD',sbj='Sam' ))
LDW <- unlist( lapply(1:length(cmpStrList), loadN1, DVstr='slopeLD',sbj='Walter' ))


loadN1 <- function(sbj, n=1, dS='k', who='all', DVstr='resAll' ){ 
  file=paste(rda(),sbj,'_',dS,'_lmSummary_',who,'.rda',sep='')
  print(file)
  load(file=paste(rda(),sbj,'_',dS,'_lmSummary_',who,'.rda',sep=''))
  return(resList[[n]][[DVstr]])
}


pSgn <- c(' ','.','*','**','***')
pVal <- c(0.1,0.05,0.01,0.001)
p2sgn <- function(pv){  pSgn[1+sum(pv<pVal)] }

allPval <- c(TS,TW,DS,DW,TDS,TDW,LDS,LDW)
allPsgn <- unlist(lapply(allPval, p2sgn))
finalPsgn <- ''
for(i in 1:length(allPsgn)){
  finalPsgn <- paste( finalPsgn, allPsgn[i], ifelse( (i)%%length(cmpStrList)==0, '\n','\t'), sep='')
}

cat(finalPsgn)


## ===================================
## quantify the effect size in uV/dB


getEffectSizes <- function(n,whatstr='TD'){
  dfLTlist <- lapply(subjList, loadN1, n=6, DVstr='dfLT')
  
  getLT <- function(n){ tapply(dfLTlist[[n]][[whatstr]], dfLTlist[[n]]$drug, mean)[1] }
  
  res   <- sapply(1:length(dfLTlist), getLT )
  return(  c(mean(res), na.se(res)  )  )
}

getEffectSizes(6, whatstr='TD')
getEffectSizes(6, whatstr='LD')


getReduction <- function(n,whatstr='TD'){
  dfLTlist <- lapply(subjList, loadN1, n=6, DVstr='dfLT')
  
  getLT <- function(n){ 
    tmp <- tapply(dfLTlist[[n]][[whatstr]], dfLTlist[[n]]$drug, mean) 
    return( 1 - tmp[2]/tmp[1] )
    }
  
  res   <- sapply(1:length(dfLTlist), getLT )
  return(  c(mean(res), na.se(res)  )  )
}
getReduction(6, whatstr='TD')
getReduction(6, whatstr='LD')




## try fitting a linear model to population data
cmpInd <- 6 # N1 component
tmp    <- LDmat[[cmpInd]]
tst    <- data.frame( subj=rep(subjList,each=dim(tmp)[1]*dim(tmp)[2]), drug=rep(c(F,T),each=dim(tmp)[1]),IV=c(rep(1:dim(tmp)[1])), DV=c(tmp) ) 

lmTst <- aov( DV~drug+IV+drug:IV + Error(subj),data=tst )
summary(lmTst)

lmTst <- aov( DV~drug+IV+drug:IV+subj,data=tst )
summary(lmTst)

idata <-  data.frame(drug=rep(c(F,T),each=dim(tmp)[1]), IV=c(rep(1:dim(tmp)[1])))
rmdat <- t(array( c(tmp), c( dim(tmp)[1]*dim(tmp)[2], dim(tmp)[3] )  ))

# run the tests: 
lm0    <- lm(  rmdat ~ 1)

anova(lm0, M=~IV+drug, X= ~drug,       idata=idata, test="Spherical")
anova(lm0, M=~IV+drug, X= ~IV,         idata=idata, test="Spherical")
anova(lm0, M=~IV*drug, X= ~IV+drug,    idata=idata, test="Spherical")




## not functional 
getCoefMat <- function(datSet,pThresh=0.05, cmpStr='n1',who='all'){
  
  loadLinearModels <- function(subj, dS='l',who='all', DVstr='lmAllFull',cmpNr=1){ 
    load(file=paste(rda(),subj,'_',dS,'_lmSummary_',who,'.rda',sep=''))
    return(resList[[cmpNr]][[DVstr]])
  }
  
  resListX <- lapply( subjList, loadLinearModels, dS=datSet, cmpNr=which(cmpStrList==cmpStr) )
  coefMatX <- array(  unlist(resListX), c(dim(resListX[[1]]), length(resListX))   )
  dimnames(coefMatX)[[1]] <- dimnames( loadLinearModels(1, dS=datSet) )[[1]]
  dimnames(coefMatX)[[2]] <- cmpStrList
  
  tdfX <- data.frame( subj=subjList, drug=rep(datSet,length(subjList)), t(coefMatX[,cmpStr,]) )
  return(tdfX)
}


drugCol   <- c('red','orange','blue')
pthresh   <- 0.05

coefMatN1 <- getCoefDataFrame('nx', who='all')
coefMatN1[which(coefMatN1$drug=='k'),]

addCmp <- function(n,what='nrmRng', who='all',runttest=F){
  
  coefMatN1 <- getCoefDataFrame(valCmpStrList[n], who=who)
  mnN1      <- tapply(coefMatN1[[what]],    coefMatN1$drug, mean)
  seN1      <- tapply(coefMatN1[[what]],    coefMatN1$drug, na.se)
  
  
  xax <- (n-1) + seq(0,to=0.25, length=length(mnN1))
  errorplot( xax, mnN1, ey=seN1, type='p', add=T, pch=20, lwd=2, col=drugCol )
  
  if(runttest){
    
    offst <- 0
    if (identical(what,'nrmRng'))
      offst <- 1
    
    ttN1      <- tapply(coefMatN1[[what]]-offst,  coefMatN1$drug, function(x){t.test(x)$p.val})
    sigN1 <- which(ttN1<pthresh)
    if(length(sigN1)>0)
      points(xax[sigN1], rep(2.0,length(sigN1)), pch=8, col=drugCol[sigN1], cex=0.75 ) 
  }
}

valCmpStrList <- c('p1a','p1b','px','n1','p2')

## relative attenuation
eps.file('GA_nrmRng_BP_all.eps', pdf=T, s=4.5, ratio=1/2)
plot(  c(-1), c(-1),  ylim=c(-.5, 2),  xlim=c(0, length(valCmpStrList)-.5 ), xaxt='n' )
axis(1, at=0.125+(0:(length(valCmpStrList)-1)), label=valCmpStrList);abline(h=c(0,1,1.9))
disc <- lapply(1:length(valCmpStrList),addCmp, what='nrmRng', who='all', runttest=T)
dev.off()

eps.file('GA_nrmRng_random.eps', pdf=T, s=4.5, ratio=1/2)
plot(  c(-1), c(-1),  ylim=c(-.5, 2),  xlim=c(0, length(valCmpStrList)-.5 ), xaxt='n' )
axis(1, at=0.125+(0:(length(valCmpStrList)-1)), label=valCmpStrList);abline(h=c(0,1,1.9))
disc <- lapply(1:length(valCmpStrList),addCmp, what='nrmRng', who='random',runttest=T)
dev.off()

## offset drug condition
plot(  c(-1), c(-1),  ylim=c(-15, 15),  xlim=c(0, length(valCmpStrList)-0.5 ), xaxt='n' )
axis(1, at=0.125+(0:(length(valCmpStrList)-1)), label=valCmpStrList);abline(h=c(0))
disc <- lapply(1:length(valCmpStrList),addCmp, what='drugTRUE',runttest=T)

## raw amplitude of saline condition
plot(  c(-1), c(-1),  ylim=c(-30, 30),  xlim=c(0, length(valCmpStrList)-0.5 ), xaxt='n' )
axis(1, at=0.125+(0:(length(valCmpStrList)-1)), label=valCmpStrList);abline(h=0)
disc <- lapply(1:length(valCmpStrList),addCmp, what='drugFALSE.isiExpAllNrm')

## raw amplitude of drug condition
plot(  c(-1), c(-1),  ylim=c(-30, 30),  xlim=c(0, length(valCmpStrList)-0.5 ),xaxt='n' )
axis(1, at=0.125+(0:(length(valCmpStrList)-1)), label=valCmpStrList);abline(h=0)
disc <- lapply(1:length(valCmpStrList),addCmp, what='drugTRUE.isiExpAllNrm')



## ===================================================================
## repeated measures ANOVAs

coefMatList <- lapply(cmpStrList, getCoefDataFrame, who='all')
coefMatDF   <-  coefMatList[[1]]
for (i in 2:length(cmpStrList) )
  coefMatDF <- rbind(coefMatDF, coefMatList[[i]])

#valCmp <- c('p1a', 'p1b', 'n1')
valCmp <- c('p1a', 'p1b','px','n1','p2')
#valCmp <- cmpStrList

valDrg <- c('l','k','m')
#valDrg <- dataSetList

y      <- subset(coefMatDF,coefMatDF$cmp%in%valCmp&coefMatDF$drug%in%valDrg, drop.levels=T)

dfDrug <- aggregate(y$nrmRng, list(y$drug), mean)
lm0    <- aov( x ~ -1, data=dfDrug ) 
lm1    <- aov( x ~ 1,  data=dfDrug )
anova(lm0,lm1)

icmp               <- as.ordered(factor(rep(valCmp,each=length(valDrg))  ,levels=valCmp  ))
iDrug              <- as.ordered(factor(rep(valDrg,length(valCmp))  ,levels=valDrg  ))
idata              <- data.frame(cmp=icmp,drug=iDrug) #cmp=gl( length(cmpStrList),   1, length(cmpStrList), labels=cmpStrList))

# creating data matrix with subjects on the rows and repetitions on the columns
rmdat <- array( y$nrmRng, c(length(unique(y$subj)), length(valDrg)*length(valCmp)  ) )


# run the tests: 
lm0    <- lm(  rmdat ~ 1)

anova(lm0, M=~1, X= ~0,                idata=idata, test="Spherical")
anova(lm0, M=~cmp+drug, X= ~cmp,       idata=idata, test="Spherical")
anova(lm0, M=~cmp+drug, X= ~drug,      idata=idata, test="Spherical")
anova(lm0, M=~cmp*drug, X= ~cmp+drug,  idata=idata, test="Spherical")

tapply(y$nrmRng, y$drug, mean)
tapply(y$nrmRng, y$cmp, mean)
tapply(y$nrmRng, list(y$cmp,y$drug), mean)

## =========================================================
postHocPaired(valCmpStrList, valDrg=valDrg ) # full analysis;  main effect component, interaction

vD <- c('l','k','m')
postHocPaired(c('p1a','p1b'),valDrg=vD )  ##x interaction: stronger effect of ket and mk on p1a, stronger effect of mid on p1b
postHocPaired(c('p1a','px'), valDrg=vD )  ### nothing
postHocPaired(c('p1a','n1'), valDrg=vD )  #@# nothing  n1 more affected than p1a
postHocPaired(c('p1a','p2'), valDrg=vD )  #@@ nothing  p2 more affected than p1a;      p2 more affected under midazolam

postHocPaired(c('p1b','px'), valDrg=vD )  #@# nothing  px more affected than p1b
postHocPaired(c('p1b','n1'), valDrg=vD )  #@# nothing  n1 more affected than p1b
postHocPaired(c('p1b','p2'), valDrg=vD )  #x#          stronger suppression for p2

postHocPaired(c('px','n1'), valDrg=vD )   ###  nothing
postHocPaired(c('px','p2'), valDrg=vD )   ##@  nothing                                 p2 less affected by ket and mk, p2 more affected by mid
postHocPaired(c('n1','p2'), valDrg=vD )   ##@  nothing                                 n1 more affected by ket and mk, p2 more affected by mid


## Interaction involving ketamine and MK801
vD <- c('l','k')
postHocPaired(c('p1a','p1b'),valDrg=vD )  #x#           P21>P31
postHocPaired(c('p1a','px'), valDrg=vD )  ### nothing
postHocPaired(c('p1a','n1'), valDrg=vD )  #x#           N85>P21 
postHocPaired(c('p1a','p2'), valDrg=vD )  ###  

postHocPaired(c('p1b','px'), valDrg=vD )  #@#           P55>P31 
postHocPaired(c('p1b','n1'), valDrg=vD )  #x#           N85>P31  
postHocPaired(c('p1b','p2'), valDrg=vD )  ##@                            P135 affected more by ket than mk

postHocPaired(c('px','n1'), valDrg=vD )   ###  nothing
postHocPaired(c('px','p2'), valDrg=vD )   ###  nothing                                
postHocPaired(c('n1','p2'), valDrg=vD )   #@#          N85>P135                              


## Interaction involving ketamine and Midazolam
vD <- c('l','m')
postHocPaired(c('p1a','p1b'),valDrg=vD )  ##x                            P21 more affected by ket, P31 more affected by mid
postHocPaired(c('p1a','px'), valDrg=vD )  ### nothing
postHocPaired(c('p1a','n1'), valDrg=vD )  ### nothing
postHocPaired(c('p1a','p2'), valDrg=vD )  ##x                            P135 more affected by mid, P21 less affected by mid  

postHocPaired(c('p1b','px'), valDrg=vD )  ### nothing     
postHocPaired(c('p1b','n1'), valDrg=vD )  ##@                            P31 less affected by ket, N85 more affected by ket 
postHocPaired(c('p1b','p2'), valDrg=vD )  #@#           P135>P31

postHocPaired(c('px','n1'), valDrg=vD )   ###  nothing
postHocPaired(c('px','p2'), valDrg=vD )   ###  nothing                                
postHocPaired(c('n1','p2'), valDrg=vD )   ##@                            P135 more affected by Mid, N85 more affected by ket


## Interaction involving MK801 and Midazolam
vD <- c('k','m')
postHocPaired(c('p1a','p1b'),valDrg=vD )  ##x                            P21 more affected by mk, P31 more affected by mid
postHocPaired(c('p1a','px'), valDrg=vD )  #@#           P55>P21
postHocPaired(c('p1a','n1'), valDrg=vD )  ### nothing
postHocPaired(c('p1a','p2'), valDrg=vD )  #x@           P135>P21         P135 more affected by mid, P21 less affected by mk  

postHocPaired(c('p1b','px'), valDrg=vD )  ### nothing     
postHocPaired(c('p1b','n1'), valDrg=vD )  ##@                            P31 less affected by mk, N85 more affected by mk 
postHocPaired(c('p1b','p2'), valDrg=vD )  #@#           P135>P31

postHocPaired(c('px','n1'), valDrg=vD )   ###  nothing
postHocPaired(c('px','p2'), valDrg=vD )   ##@                            P135 more affected by Mid, P55 more affected by mk                                 
postHocPaired(c('n1','p2'), valDrg=vD )   ##@                            P135 more affected by Mid, N85 more affected by mk





postHocPaired <- function(valCmp,valDrg=c('l','k','m')){
  y      <- subset(coefMatDF,coefMatDF$cmp%in%valCmp&coefMatDF$drug%in%valDrg, drop.levels=T)
  
  dfDrug <- aggregate(y$nrmRng, list(y$drug), mean)
  lm0    <- aov( x ~ -1, data=dfDrug ) 
  lm1    <- aov( x ~ 1,  data=dfDrug )
  anova(lm0,lm1)
  
  icmp               <- as.ordered(factor(rep(valCmp,each=length(valDrg))  ,levels=valCmp  ))
  iDrug              <- as.ordered(factor(rep(valDrg,length(valCmp))  ,levels=valDrg  ))
  idata              <- data.frame(cmp=icmp,drug=iDrug) #cmp=gl( length(cmpStrList),   1, length(cmpStrList), labels=cmpStrList))
  
  # creating data matrix with subjects on the rows and repetitions on the columns
  rmdat <- array( y$nrmRng, c(length(unique(y$subj)), length(valDrg)*length(valCmp)  ) )
  
  # run the tests: 
  lm0    <- lm(  rmdat ~ 1)
  mainDrug <- anova(lm0, M=~cmp+drug, X= ~cmp,       idata=idata, test="Spherical")
  mainCmp  <- anova(lm0, M=~cmp+drug, X= ~drug,      idata=idata, test="Spherical")
  InterDC  <- anova(lm0, M=~cmp*drug, X= ~cmp+drug,  idata=idata, test="Spherical")
  
  print( tapply(y$nrmRng, y$drug, mean) )
  print( tapply(y$nrmRng, y$cmp, mean)  )
  print( tapply(y$nrmRng, list(y$cmp,y$drug), mean) )
  
  pvals <- c( mainDrug$`H-F Pr`[1],  mainCmp$`H-F Pr`[1], InterDC$`H-F Pr`[1])
  return(pvals)
}










## ==================================================================
## old 

## confirm the correct formating:
apply(rmdat,2,mean)
tapply(y$nrmRng, list(y$cmp,y$drug), mean)

mnbyAnimal <- tapply(y$nrmRng,y$subj,mean)
y$nrmRngDM <- y$nrmRng - mnbyAnimal[y$subj]

lm0 <- aov( I(nrmRng-1) ~ -1 + drug + cmp, data=y)
lm1 <- aov( I(nrmRng-1) ~ +1 + drug + cmp, data=y)
anova(lm0,lm1)

lm0 <- aov( nrmRng ~ 1 + drug + cmp + subj + drug:cmp , data=y)
Anova(lm0)

lm0 <- aov( nrmRngDM ~ 1 + drug + cmp + drug:cmp, data=y)
Anova(lm0)

lm0 <- lm( nrmRng ~ 1 + cmp, data=y)

## this version does not account for violations of sphericity
lm0 <- aov( nrmRng ~ 1 + drug * cmp + Error(subj / (drug*cmp)), data=y)
summary(lm0)

