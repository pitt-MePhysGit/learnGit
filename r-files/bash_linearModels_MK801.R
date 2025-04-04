## ================================= ##
args <- commandArgs( trailingOnly=TRUE )
defpar <- par()

library(car)
baseDir <- '~/Dropbox/'
source(paste(baseDir, 'ampOdd_click/r-files/source_all_functions.R',sep=''))

pp <- function(){paste(baseDir, 'ampOdd_click/word/plots/',sep='')}
rda <- function(){'~/Dropbox/ampOdd_click/Rdata/'}

animal  <- 'Sam'
dataSet <- 'k'
speaker <- 'TB'
who     <- 'all'

what    <- 'llp2'

if (length(args)>0)
  animal   <- args[1]
if (length(args)>1)
  dataSet  <- args[2]
if (length(args)>2)
  speaker  <- args[3]
if (length(args)>3)
  who  <- args[4]


if (dataSet=='k'){
  #pp <- function(){paste(baseDir, 'RSkernel/word/mk801/plots/',sep='')}
  doseList <- c(0,0.1)
}

trialTypeInd <- c(1,2)
if (identical(who,'random'))
  trialTypeInd <- c(1)
if (identical(who,'regular'))
  trialTypeInd <- c(2)

tdf <- df[which(df$dataSet==dataSet), ]  # n for NMDA

## load the data
chStr <- NA
indS  <- which( (df$animal==animal & df$dataSet==dataSet & df$speaker==speaker & df$ket==doseList[1]) | (df$animal==animal & df$dataSet=='g' & df$speaker==speaker) )
indK  <- which( df$animal==animal & df$dataSet==dataSet & df$speaker==speaker & df$ket==doseList[2])
indA  <- which( df$animal==animal & df$dataSet%in% c(dataSet,'g') & df$speaker==speaker )

#sal        <- loadAmpOddClickAll( df[indS,] , what=what, chStr=chStr, show=F)
#ket        <- loadAmpOddClickAll( df[indK,] , what=what, chStr=chStr, show=F)
all         <- loadAmpOddClickAll( df[indA,], what=what, chStr=chStr, show=F)

#sal$xpp$ISI[sal$xpp$ISI<0]  <- 20
#sal$xpp$ISI[sal$xpp$ISI>20] <- 20

#ket$xpp$ISI[ket$xpp$ISI<0]  <- 20
#ket$xpp$ISI[ket$xpp$ISI>20] <- 20

all$xpp$ISI[all$xpp$ISI<0]  <- 20
all$xpp$ISI[all$xpp$ISI>20] <- 20


mxISI <- 6.4
mxp2p <- 1000
mxpow <- Inf
if(identical(animal,'Sam')){
  mxp2p <- 1000
  mxpow <- Inf # 20000
}

#sal$xpp$drug <- FALSE
#ket$xpp$drug <- TRUE
#all                 <- append.data(sal, ket)

all$xpp$drug <- c(1,0)[ 1 + c(all$xpp$exp) %in% which( df$ket[indA]==0 ) ]
all$xpp$isiLog      <- log2(all$xpp$ISI)
all                 <- subs.data(all,all$xpp$trialNr>4 & all$xpp$p2p<mxp2p  & all$xpp$powEOG<mxpow & all$xpp$trialType%in%trialTypeInd & all$xpp$ISI<mxISI)


all$xpp$DV <- all$xpp$n1

cmpStrList <- c('p1x.r','p1a.r','p1b.r','nx.r','px.r','n1.r','p2.r','n2.r')
cmpStrList <- c('p1x','p1a','p1b','nx','px','n1','p2','n2')
resList    <- lapply(cmpStrList, runAllTests_ampOddClick, thisall=all )

save(resList, file=paste(rda(),animal,'_',dataSet,'_lmSummary_',who,'.rda',sep='') )

