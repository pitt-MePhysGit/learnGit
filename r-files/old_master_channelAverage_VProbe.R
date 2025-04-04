## ============================================
## master channel-average VProbe
##  in contrast to master_sessionAverage_VProbe, this script derives channel-based measures and error bars 

rm(list=ls())
efpar <- par()

baseDir <- '~/Dropbox/'
source(paste(baseDir, 'RSkernel/r-files/functions/source_all_functions.R',sep=''))
source(paste(baseDir, 'ampOdd_click/r-files/source_all_functions.R',sep=''))
source( paste(baseDir,'RSkernel/r-files/functions/expand2chdf.R',sep='') )

## preliminaries
.GlobalEnv$pp        <- function(){paste(baseDir, '/ampOdd_click/word/VProbe/plots/',sep='')}
.GlobalEnv$rda       <- function(){'~/Dropbox/ampOdd_click/Rdata/'}
.GlobalEnv$trda      <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
.GlobalEnv$eeglabdir <- function(){'/Volumes/Drobo5D3/EEG/EEGLab/ampOddClick/'}

source( paste(baseDir,'ampOdd_click/r-files/functions/prepDFVprobe.R',sep='') )
source( paste(baseDir,'ampOdd_click/r-files/functions/prep_df_mua.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/functions_sessionAverageVProbe.R',sep='') )

source( paste(baseDir,'ampOdd_click/r-files/figures/call_AOC_timeseries_bySOA.R',sep='') )
source( paste(baseDir,'ampOdd_click/r-files/figures/figure_AOC_timeseries_bySOA.R',sep='') )

source( paste(baseDir,'ampOdd_click/r-files/figures/call_AOC_timeseries_byAMP.R',sep='') )
source( paste(baseDir,'ampOdd_click/r-files/figures/figure_AOC_timeseries_byAMP.R',sep='') )

source( paste(baseDir,'ampOdd_click/r-files/functions/loadGASession.R',sep='') )

colClasses <- c('character',rep('numeric',2),'character',rep('numeric',5)  ,rep('character',3),rep('numeric',8) )
df         <- read.table(paste(baseDir,'ampOdd_click/ampOddClickVProbe_link.txt',sep=''),header=T,colClasses=colClasses)
df         <- prepDFVprobe(df)
df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))


cleanOnly   <- FALSE
#badDays     <- c(12, 14, 16, 18, 22, 24)

df$use      <- c(1)
if (cleanOnly)
  df$use[ which(df$toneArtefact+df$onsetArtefact>0) ] <- 0
#df$use[badDays] <- 0

dataSet <- 'g'
animal <- c('Sam')#,'Walter')

useInd        <- which( df$dataSet==dataSet & df$animal%in%animal & df$drug=='N' & df$use )# & df$probeString=='t')
#useInd        <- useInd[5:6]
#useInd       <- useInd[-c(8)]

tdf           <- df[useInd,]
chdf          <- expand2chdf(tdf)

## ======================================== set folder for plots
animalStr <- '_all'
if (length(animal)==1)
  animalStr <- paste('_',animal,sep='')

cleanStr  <- ''
if (cleanOnly)
  cleanStr <- '_clean'

uIdir <- paste('g', animalStr,cleanStr,sep='')
.GlobalEnv$pp   <- function(){paste(baseDir, '/ampOdd_click/word/VProbe/plots/',uIdir,'/',sep='')}
if (!file.exists(pp()))
  system( paste('mkdir', pp()) )



#tst <- loadGASession(useInd[40], extn='mean',who='all',what='csd',xnd=1, ynd=NA, ganrm=F)
loadAllSessions <- function(useInd, extn='mean',who='all',what='mua',xnd=1, ynd=NA,dxnd=NA,ganrm=F){
  
  resList <- lapply(useInd, loadGASession, extn=extn,who=who,what=what,xnd=xnd, ganrm=ganrm )
  
  Nchan   <- sum( sapply(1:length(resList), function(n){dim(resList[[n]])[2]}))
  resMat  <- array( unlist(resList), c(dim(resList[[1]])[1], Nchan) )
}
#tst <- loadAllSessions(useInd, extn='mean',who='all',what='csd',xnd=1, ynd=NA)

## =============================== PLOTS FOR BRAIN INITIATIVE MEETING 2019 ============================= ##
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


