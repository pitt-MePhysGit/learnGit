#
# compare the old long-latency preprocessing (llp2, higher-order blackman) with new one (nlp,4th-order blackman)
#

baseDir <- '~/Dropbox/'
source(paste(baseDir, 'RSkernel/r-files/functions/source_all_functions.R',sep=''))
source(paste(baseDir, 'ampOddClick/r-files/source_all_functions.R',sep=''))

## preliminaries
.GlobalEnv$pp        <- function(){paste(baseDir, '/ampOddClick/word/VProbe/plots/',sep='')}
.GlobalEnv$rda       <- function(){'~/Dropbox/ampOddClick/Rdata/'}
.GlobalEnv$trda      <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
.GlobalEnv$eeglabdir <- function(){'/Volumes/Drobo5D3/EEG/EEGLab/ampOddClick/'}

source( paste(baseDir,'ampOddClick/r-files/functions/prepDFVprobe.R',sep='') )
source( paste(baseDir,'ampOddClick/r-files/functions/prep_df_mua.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/functions_sessionAverageVProbe.R',sep='') )

source( paste(baseDir,'r-compile/SignalProcessor/R/SignalProcessor.R',sep='') )

#colClasses <- c('character',rep('numeric',2),'character',rep('numeric',5)  ,rep('character',3),rep('numeric',8) )
#df         <- read.table(paste(baseDir,'ampOddClick/ampOddClickVProbe_link.txt',sep=''),header=T,colClasses=colClasses)

colClasses <- c('character',rep('numeric',9),rep('character',3),rep('numeric',1) )
df         <- read.table(paste(baseDir,'ampOddClick/ampOddClickEEG_link.txt',sep=''),header=T,colClasses=colClasses)

#df         <- prepDFVprobe(df)
df$animal  <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))


# ---------------------- load the old preporcessing routine
# load new preprocessing routine
chnd <- 2:33
eeg  <- loadSignalProcessor(df[677,], project='ampOddClick', epoch='llp', chStr=paste('ch',chnd,sep='')) 
raw  <- loadSignalProcessor(df[678,], project='ampOddClick', epoch='raw', chStr=paste('ch',chnd,sep=''),chBlkN=16) 

showSignalProcessor(raw,                 subs=raw$xpp$p2p<Inf, trng=c(-10,150),sclfctr=.5)
showSignalProcessor(raw, IVstr='ampInd', subs=raw$xpp$p2p<Inf, trng=c(-10,150),sclfctr=.5)
showSignalProcessor(raw, IVstr='ICIoct', subs=raw$xpp$p2p<Inf, trng=c(-10,150),sclfctr=.5,colFUN=soa.colors)



mean.plot( eeg$lfp, eeg$taxis, widthFUN=na.se, add=F, ylim=c(-20,25), xlim=c(-20,250))
mean.plot( nrw$lfp, nlp$taxis, widthFUN=na.se, add=T, col='black')
mean.plot( nlp$lfp, nlp$taxis, widthFUN=na.se, add=T, col='red')


plot(  nlp$lfp[1,,1], type='l',col='red')
points(nrw$lfp[1,,1], type='l',col='blue' )

