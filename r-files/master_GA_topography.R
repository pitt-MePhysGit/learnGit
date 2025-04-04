source('~/Dropbox/r-compile/signalProcessoR/R/signalProcessoR.R')

animal      <- 'Jesse'
project     <- 'ampOddClick'
linkFile    <- paste(project,'_link.txt',sep='')
projectDir  <- paste('~/Dropbox/',project,'/',sep='')

pp          <- function(){paste(projectDir, '/plots/GA_topography_',animal,'/',sep='')}
if (!exists(pp()))
  system(paste('mkdir -p ', pp() ))

# ============================ read the link file
df          <- prepDFsignalProcessor(project)

# ============================ select a subset of recordings to use
ind <- which( df$animal==animal & 1:dim(df)[1]%in%grep('20210',df$date))
tdf <- df[ind[c(2,3)],]
print(tdf)

# ================= select epoch, analysis and channels 
epoch    <- 'raw'
analysis <- 'ampInd'
chStr <- paste('ch',34:129,sep='')

# =============== read the data
tn                   <- loadSignalProcessor(tdf, project=project,epoch=epoch,analysis=analysis,chStr=chStr)
tn$channelPosition   <- getSC96pos()


# determine all channels in A1
scx <- tn$channelPosition$x
scy <- tn$channelPosition$y
A1bin <- scx>=8 | (scx==5 & scy <= -9 & scy>= -10)| (scx==6 & scy <= -4) | (scx==7 & scy <= -3)
xxbin <- (scx==7 & scy== -11) | (scx==8 & scy== -11) | (scx==9 & scy== -2) | (scx==10 & scy== -3) | (scx==11 & scy== -11) | (scx==11 & scy== -7)

A1ndx <- which( A1bin & !xxbin)
MPndx <- which(!A1bin & !xxbin)



# average across days
mn <- avgSignalProcessor(tn, frml=~ampInd)


## ------------------------- summary topo-plot
eps.file(paste('GA_MultiPlot_AMP_',epoch,'.eps',sep=''),pdf=T, s=5)
showSignalProcessor(mn, frml=~ampInd,trng=c(-1,10),colFUN=bluered.colors)
dev.off()

eps.file(paste('GA_MultiPlot_',epoch,'.eps',sep=''),pdf=T, s=5)
showSignalProcessor(mn, frml=~sessionID,trng=c(-1,10),colFUN=bluered.colors)
dev.off()



# 
mnGA        <- avgSignalProcessor(mn,frml=~sessionID)

eps.file(paste('GA_tonotopy_onset_',epoch,'.eps',sep=''),pdf=T, s=5)
mn$onsetMap <- mnGA$lfp[1,,]
mn$onsetMap[,54] <- c(0)
showSignalProcessorMap(mn,eval.at=seq(0,by=1,length=20),taxisStr='taxis',dataStr='onsetMap',zlim=c(-2,2))
dev.off()

