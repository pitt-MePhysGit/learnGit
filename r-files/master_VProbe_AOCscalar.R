## summary-figures of results of RSscalar analysis

#rm(list=ls())
efpar <- par()

library(VennDiagram)

baseDir <- '~/Dropbox/'
source(paste(baseDir, 'RSkernel/r-files/functions/source_all_functions.R',sep=''))
source( paste(baseDir,'ampOddClick/r-files/functions/functions_AOCdual.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/prepDFVprobe.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/expand2chdf.R',sep='') )
source( paste(baseDir,'ampOddClick/r-files/functions/function_VProbe_AOCscalar.R',sep='') )


df       <- prepDFsignalProcessor('ampOddClick','ampOddclickVProbe_link.txt')
df$exp   <- df$session
df         <- prepDFVprobe(df,getTrTy=F)

dataSetList <- c('g')
animalList <- list()
animalList[[1]]  <- ('Fuzzy')
#animalList[[2]]  <- ('Walter')
#animalList[[1]]  <- c('Sam','Walter')

animals <- c('Sam','Walter')
whoList     <- c('all')
cmpStrList  <- c('MUA1','MUA1a','MUA1b','MUA1c','MUA2','MUA3'
                 ,'CSD1','CSD2','CSD3','CSD4','CSD5'
                 ,'LFP0','LFP1'
                 ,'MUAbl','MUAR1','MUAR1a','MUAR1b','MUAR1c','MUAR2','MUAR3'
)
#cmpStrList <- c('MUA1')#,'MUA2','MUAbl','MUAR1','MUAR2')   
#cmpStrList <- c('MUAR1a','MUAR1b','MUAR1c','MUAR2','MUAR3')
cmpStrList  <- c('CSD1','CSD2','CSD3','CSD4','CSD5')
dpthStrList <- c('All','SG','IG','L2','L3','L4','L5','SGSi','SGSo','GSi','MUA')
#dpthStrList <- c('MUA')


if (1==0){
  cmpStrList  <- c('MUAS')
  
  disc <- function_VProbe_AOCscalar(dataSet='g',animals=c('Fuzzy'),who='all',cmpStr='MUA1',dpthStr='All', pthrAmp = 0.05/1700, pthrISI = 0.05/1700, makeFigs = FALSE)
  #disc <- function_VProbe_AOCscalar(dataSet='g',animals=c('Walter'),who='all',cmpStr='MUAS',dpthStr='All', pthrAmp = 0.05/1700, pthrISI = 0.05/1700, makeFigs = FALSE)
  #disc <- function_VProbe_AOCscalar(dataSet='g',animals=c('Sam'),who='all',cmpStr='MUAS',dpthStr='All', pthrAmp = 0.05/1700, pthrISI = 0.05/1700, makeFigs = FALSE)
  
  #compare 18, 35 ms window to 0, 70
  
  #channelDat<-disc
  #save(channelDat, file = '~/Code/Data/WalterChannelQuality70.rda')
  #save(channelDat, file = '~/Code/Data/SamChannelQuality70.rda')
  
  #load('~/Code/Data/WalterChannelQuality2.rda')
  #load('~/Code/Data/SamChannelQuality.rda')
  
  checkv<-unique(channelDat$session)
  outv<-matrix(0, length(checkv), 2)
  for(i in 1:length(checkv)){
    if(sum(disc$session==checkv[i])>0){
      outv[i, 1]<-min(which(channelDat$session==checkv[i]))
      outv[i, 2]<-min(which(disc$session==checkv[i]))
    }
  }
  outv<-outv[outv[,1]!=0,]
  
  table(channelDat$valBin[unlist(lapply(outv[,1], function(x){x:(x+23)}))], disc$valBin[unlist(lapply(outv[,2], function(x){x:(x+23)}))])
  
  #89% agreement for Walt
  
  
  ##which(channelDat$valBin[which(unlist(lapply(channelDat$session, function(x){grepl('20190716', x)})))])
}

## ====================================== only MUA
cmpStrList  <- c('MUA0','MUA1','MUA1a','MUA1b','MUA1c','MUA2','MUA3')
dpthStrList <- c('MUA')#All','SG','IG','L2','L3','L4','L5','SGSi','SGSo','GSi','MUA')
for (ax in 1:length(animalList)){
  chdfAMM1 <- function_VProbe_AOCscalar(dataSet='g',animals=animalList[[ax]], who='all',    cmpStr='MUA1',  dpthStr='All', pthrAmp = 0.001, pthrISI = 0.001, makeFigs = FALSE)
  for(dx in 1:length(dpthStrList)){
    for(cx in 1:length(cmpStrList)){
      disc <- function_VProbe_AOCscalar(dataSet='g',animals=animalList[[ax]],who='all',cmpStr=cmpStrList[cx],dpthStr=dpthStrList[dx],  valBin=chdfAMM1$valBin, makeFigs = TRUE)      
    }
  }
}


## ====================================== only CSD
cmpStrList  <- c('CSD0','CSD1','CSD2','CSD3','CSD4','CSD5')
dpthStrList <- c('SGSi','SGSo','GSi')

for (ax in 1:length(animalList)){
  chdfAMM1 <- function_VProbe_AOCscalar(dataSet='g',animals=animalList[[ax]], who='all',    cmpStr='CSD1',  dpthStr='All', pthrAmp = 1, pthrISI = 1)
  for(dx in 1:length(dpthStrList)){
    for(cx in 1:length(cmpStrList)){
      disc <- function_VProbe_AOCscalar(dataSet='g',animals=animalList[[ax]],who='all',cmpStr=cmpStrList[cx],dpthStr=dpthStrList[dx],  valBin=chdfAMM1$valBin)      
    }
  }
}





for (cx in 1:length(cmpStrList)){
  for(dx in 1:length(dpthStrList)){
    #disc <- function_VProbe_RSscalar(dataSet='g',animals=animals,who='regular',cmpStr=cmpStrList[cx],dpthStr=dpthStrList[dx], shft=chdfAMM1$shft)
  }
}


#disc <- function_VProbe_RSscalar(dataSet='g',animals=c('Sam'),who='all',cmpStr='MUA0',dpthStr='MUA')
#disc <- function_VProbe_RSscalar(dataSet='g',animals=c('Sam'),who='all',cmpStr='MUA1',dpthStr='MUA')

#disc <- function_VProbe_RSscalar(dataSet='g',animals=c('Sam'),who='all',cmpStr='MUA2',dpthStr='L3')
#disc <- function_VProbe_RSscalar(dataSet='g',animals=c('Sam'),who='all',cmpStr='MUA2',dpthStr='L4')
#disc <- function_VProbe_RSscalar(dataSet='g',animals=c('Sam'),who='all',cmpStr='MUA2',dpthStr='L5')

#disc <- function_VProbe_RSscalar(dataSet='g',animals=c('Sam'),who='all',cmpStr='MUA1a',dpthStr='MUA')
#disc <- function_VProbe_RSscalar(dataSet='g',animals=c('Sam'),who='all',cmpStr='MUA1b',dpthStr='MUA')
#disc <- function_VProbe_RSscalar(dataSet='g',animals=c('Sam'),who='all',cmpStr='MUA1c',dpthStr='MUA')
#disc <- function_VProbe_RSscalar(dataSet='g',animals=c('Sam'),who='all',cmpStr='MUA2',dpthStr='MUA')

#disc <- function_VProbe_RSscalar(dataSet='g',animals=c('Sam'),who='all',cmpStr='CSD1',dpthStr='SGSi')
#disc <- function_VProbe_RSscalar(dataSet='g',animals=c('Sam'),who='all',cmpStr='CSD1',dpthStr='SGSo')

#disc <- function_VProbe_RSscalar(dataSet='g',animals=c('Sam'),who='all',cmpStr='CSD4',dpthStr='SGSi')
#disc <- function_VProbe_RSscalar(dataSet='g',animals=c('Sam'),who='all',cmpStr='CSD4',dpthStr='SGSo')


## Questions
# why is absDeltaPitch >0 across the board? Should be 0 on average?
# why is absDeltaPitch effect so small?



if (1==0){
  
  colClasses <- c('character',rep('numeric',3),'character',rep('numeric',4)  ,rep('character',3),rep('numeric',3) )
  df         <- read.table(paste(baseDir,'RSkernel/RSkernelVProbe_link.txt',sep=''),header=T,colClasses=colClasses)
  df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))
  df$exp     <- df$session
  df         <- prepDFVprobe(df)
  
  # select valid sessions
  dataSet    <- 'g' 
  animals    <- c('Sam')#,'Walter')
  
  who     <- 'all'
  cmpStr  <- 'MUA1c'
  dpthStr <- 'MUA'
  
  
  ##  paths and such
  rda <- function(){'/Volumes/Drobo5D3/EEG/RSkernel/rda/'}
  trda <- function(){'/Volumes/Drobo5D3/EEG/RSkernel/rda/'}
  spd  <- function(){'/Volumes/Drobo5D3/EEG/EEGLab/RS/'}
  options(warn=0)
  
  ## ===============================================
  badDays <- c(9,11,12)+1  # what about 7 and 10?
  df$use  <- c(1)
  
  animalStr <- '_all'
  if (length(animals)==1)
    animalStr <- paste('_',animals,sep='')
  
  cleanOnly <- FALSE
  cleanStr  <- ''
  if (cleanOnly){
    df$use[badDays] <- 0
    cleanStr <- '_clean'
  }
  
  ## ======================================
  useInd <- which(df$drug=='N' &  df$dataSet==dataSet     &   df$animal%in%animals & df$use==1 & df$probe==1)
  #useInd <- useInd[15:20]
  chdf       <- expand2chdf(df[useInd,])
  chdf$depth <- (chdf$Layer - chdf$chIndx) * chdf$spacing
  
  # ================================
  uIdir <- paste('g', animalStr, cleanStr,'/', who,'/',cmpStr, '/', dpthStr,sep='')
  
  pp   <- function(){paste(baseDir, '/RSkernel/word/VProbe/plots/',uIdir,'/',sep='')}
  if (!file.exists(pp()))
    system( paste('mkdir -p ', pp()) )
  
  
  df[useInd,]
  
  ## ================================= find depth
  layers     <- data.frame( cmpStrName='L2',      dmin=  200, dmax= 1000, stringsAsFactors=FALSE) 
  layers[2,] <- data.frame( cmpStrName='L3',      dmin= -700,  dmax=200, stringsAsFactors=FALSE)
  layers[3,] <- data.frame( cmpStrName='L4',      dmin= -1100,  dmax= -700, stringsAsFactors=FALSE)
  layers[4,] <- data.frame( cmpStrName='L5',      dmin= -1500,  dmax= -1100, stringsAsFactors=FALSE)
  layers[5,] <- data.frame( cmpStrName='L6',      dmin= -2000,  dmax= -1500, stringsAsFactors=FALSE)
  layers[6,] <- data.frame( cmpStrName='SGSo',    dmin=     0,  dmax= 500, stringsAsFactors=FALSE)
  layers[7,] <- data.frame( cmpStrName='SGSi',    dmin= -500,  dmax=    0, stringsAsFactors=FALSE)
  layers[8,] <- data.frame( cmpStrName='GSoSi' ,  dmin= -1100,  dmax= -700, stringsAsFactors=FALSE)
  layers[9,] <- data.frame( cmpStrName='MUA' ,    dmin= -1500,  dmax= -300, stringsAsFactors=FALSE)
  layers[10,] <- data.frame( cmpStrName='SG' ,     dmin= -700,  dmax= 1000, stringsAsFactors=FALSE)
  layers[11,] <- data.frame( cmpStrName='IG' ,    dmin= -2000,  dmax= -1100, stringsAsFactors=FALSE)
  
  lnd     <- which( layers$cmpStrName %in% dpthStr )
  dpthRng <- c( unlist(layers[lnd, c('dmin','dmax')] ))
  
  
  ## ==================================================
  readRSscalar <- function(ix,who='all',cmpStr='MUA1'){
    resFileName <- paste(rda(),chdf$exp[ix],'/RSscalar_ch',chdf$chNr[ix],'_',cmpStr,'_',who,'.rda',sep='')
    print(resFileName)
    load(resFileName)  
    return(resList)
  }
  
  resList    <- lapply(1:dim(chdf)[1], readRSscalar, who=who,cmpStr=cmpStr)
  #resListReg <- lapply(1:dim(chdf)[1], readRSscalar, who='regular',cmpStr=cmpStr)
  #resListRnd <- lapply(1:dim(chdf)[1], readRSscalar, who='random',cmpStr=cmpStr)
  
  
  chdf$pValISI   <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$pValISI ) } )
  chdf$pValPitch <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$pValPitch ) } )
  
  
  ## ============================================================== MMN and ND  as a function of SOA
  vlnd           <- which(  chdf$pValISI<10^-2  & chdf$pValPitch<10^-1 & 
                              chdf$depth > dpthRng[1] & chdf$depth < dpthRng[2])
  
  tt        <- lapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVbyToneTypeTrialTypeAndSOA ) }) 
  ttMat     <- array(unlist(tt), c(dim(tt[[1]]),length(tt)) )
  
  xax                 <- as.numeric(unlist(dimnames(resList[[1]]$DVbyToneTypeTrialTypeAndSOA)[3] ))
  mnTT               <- apply(ttMat[,,,vlnd], c(1,2,3), na.mean)
  seTT               <- apply(ttMat[,,,vlnd], c(1,2,3), na.se)
  
  soand <- 1:6
  
  picFile <- paste('MMN_NDbySOA.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  #errorplot( log2(xax[soand]), mnTT[1,1,soand], ey=seTT[1,1,soand],  type='b', lwd=2,pch=20, col='gray',ylim=c(-.6,.6), add=F )# Std Rnd
  errorplot( log2(xax[soand]), mnTT[1,2,soand], ey=seTT[1,2,soand],  type='b', lwd=2,pch=20, col='black',ylim=c(-.6,.6),add=F)                  # Std Reg
  errorplot( log2(xax[soand]), mnTT[2,1,soand], ey=seTT[2,1,soand],  type='b', lwd=2,pch=20, col='blue',add=T)                  # Dev Rnd
  errorplot( log2(xax[soand]), mnTT[2,2,soand], ey=seTT[2,2,soand],  type='b', lwd=2,pch=20, col='red',add=T)                   # Dev Reg 
  dev.off()
  
  
  
  ## =================================================================== compare SOA dependence between regular and random condition
  vlnd           <- which(  chdf$pValISI<10^-2  & chdf$pValPitch<10^-2 & 
                              chdf$depth > dpthRng[1] & chdf$depth < dpthRng[2])
  
  
  rnd        <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVbyTrialTypeAndSOA[1,] ) }) 
  reg        <- sapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVbyTrialTypeAndSOA[2,] ) }) 
  
  xax                 <- as.numeric(unlist(dimnames(resList[[1]]$DVbyTrialTypeAndSOA)[2] ))
  mnReg               <- apply(reg[,vlnd], 1, na.mean)
  seReg               <- apply(reg[,vlnd], 1, na.se)
  
  mnRnd               <- apply(rnd[,vlnd], 1, na.mean)
  seRnd               <- apply(rnd[,vlnd], 1, na.se)
  
  
  picFile <- paste('TrialTypeAndSOA.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  errorplot( log2(xax[soand]), mnReg[soand], ey=seReg[soand],  type='b', lwd=2,pch=20, col='black',ylim=c(-.6,.6), add=F )
  errorplot( log2(xax[soand]), mnRnd[soand], ey=seRnd[soand],  type='b', lwd=2,pch=20, col='red',add=T)
  dev.off()
  
  
  ## ============================================ compare pitch selectivity between random and regular condition
  chdf$bestPitchInd   <- sapply( 1:dim(chdf)[1], function(ix){ return(which( resList[[ix]]$betaPitch==max(resList[[ix]]$betaPitch))[1]) } )
  
  extractFcn <- function(ix){ return(resList[[ix]]$DVbyPitchandTrialType) }
  reg <- array( c(NA), c(21,length(resList)) )
  rnd <- array( c(NA), c(21,length(resList)) )
  for (ix in 1:dim(chdf)[1] ){
    tmp                    <- extractFcn(ix)
    if ( length( which(is.finite(tmp)) )>0){
      
      mn                     <- apply(tmp,1,mean)
      mn[which(is.na(mn))]   <- na.mean(mn)
      tmpmn                  <- mn[ c(rep(1,2),1:length(mn),rep(length(mn), 2))]
      ker                    <- dnorm(-2:2, m=0,sd=.7)
      fmn <- convolve(tmpmn,ker, type='filter')
      
      chdf$bestPitchInd[ix]  <- which(fmn == na.max(fmn))
      shft               <- 11 - chdf$bestPitchInd[ix]
      rnd[shft+1:dim(tmp)[1],ix] <- tmp[,1]
      reg[shft+1:dim(tmp)[1],ix] <- tmp[,2]
    }
  }
  
  vlnd                <- which(  chdf$pValPitch<10^-2 & chdf$bestPitchInd>2 & chdf$bestPitchInd<10 & 
                                   chdf$depth > dpthRng[1] & chdf$depth < dpthRng[2] )
  
  xax                 <- (-10:10)/4
  mnReg               <- apply(reg[,vlnd], 1, na.mean)
  seReg               <- apply(reg[,vlnd], 1, na.se)
  
  mnRnd               <- apply(rnd[,vlnd], 1, na.mean)
  seRnd               <- apply(rnd[,vlnd], 1, na.se)
  
  picFile <- paste('TrialTypeAndPitch.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  errorplot( xax, mnReg, ey=seReg,  type='b', lwd=1.5,pch=20, col='black',ylim=c(-.6,.8), add=F,xlim=c(-2,2) )
  errorplot( xax, mnRnd, ey=seRnd,  type='b', lwd=1.5,pch=20, col='red',add=T)
  dev.off()
  
  dlta                 <- rnd - reg
  mnDlta               <- apply(dlta[,vlnd], 1, na.mean)
  seDlta               <- apply(dlta[,vlnd], 1, na.se)
  picFile <- paste('TrialTypeAndPitchDelta.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  errorplot( xax, mnDlta, ey=seDlta,  type='b', lwd=2,pch=20, col='black',ylim=c(-0.1,.5), add=F )
  abline(h=0);abline(v=0)
  dev.off()
  
  
  ## ============================================================================  DV by RepNr and SOA
  vlnd                <- which(  chdf$pValPitch<10^-2 & chdf$pValISI<10^-2 & 
                                   chdf$depth > dpthRng[1] & chdf$depth < dpthRng[2] )
  
  rep        <- lapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$DVbyRepNrAndSOA ) }) 
  repMat     <- array(unlist(rep), c(dim(rep[[1]]),length(rep)) )
  Nt         <- lapply( 1:dim(chdf)[1], function(ix){ return(resList[[ix]]$NbyRepNrAndSOA ) }) 
  NMat       <- array(unlist(Nt), c(dim(Nt[[1]]),length(Nt)) )
  repMat[which(NMat<20)] <- NA
  
  mnRep <- apply(repMat[,,vlnd], c(1,2), na.mean)
  seRep <- apply(repMat[,,vlnd], c(1,2), na.se)
  xax   <- 0:(dim(repMat)[1]-1)
  
  picFile <- paste('RepNrAndSOA.eps',sep='')
  eps.file(picFile, pdf=T, s=2, ratio=3/2)
  colList                <- soa.colors(length(xax))
  for (ix in 1:(dim(mnRep)[2]-1))
    errorplot( 1:length(xax), mnRep[,ix], ey=seRep[,ix],  type='b', lwd=2,pch=20, col=colList[ix],ylim=c(-.6,.6), add=ifelse(ix==1,F,T) )
  dev.off()
  
}
