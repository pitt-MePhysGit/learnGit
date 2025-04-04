## preliminaries
args <- commandArgs( trailingOnly=TRUE )

baseDir <- '~/Dropbox/'
source(paste(baseDir, 'ampOdd_click/r-files/source_all_functions.R',sep=''))

## preliminaries
pp <- function(){paste(baseDir, 'ampOdd_click/plots/',sep='')}
rda <- function(){'/Volumes/rawData/EEG/ampOddClick/rda/'}
rda <- function(){'~/Dropbox/ampOdd_click/Rdata/'}

animal <- args[1]
dataSet <- 'g'
speaker <- 'TB'
extn    <- ''

#cmpStrList <- c('n1','p1','n8','nx','px','p2','n2','n8b','p1x','p1a','p1b','p1c','p1Tb')
#cmpStrList <- c('n8b','n8','nx','px','p1x','p1a','p1b')
#cmpStrList <- c('n1','p1','p2','n2','p1c','p1Tb')

#cmpStrList <- c('p1','nx','px','n1','p2','n2','p3','n8','n8b','p1x','p1a','p1Ta','p1b','p1Tb','p1c','pSV','nSN')
#cmpStrList <- c('p2','n2','p3','n8','n8b','p1x','p1a','p1Ta','p1b','p1Tb','p1c','pSV','nSN')
#cmpStrList <- c('p1','n1','p2','n2')
#cmpStrList <- c('px','p3')

#cmpStrList <- c('p1x','p1a','p1b','nx','px','n1','p2','n2',
#                'p1x.s','p1a.s','p1b.s','nx.s','px.s','n1.s','p2.s','n2.s',
#                'p1x.r','p1a.r','p1b.r','nx.r','px.r','n1.r','p2.r','n2.r')

cmpStrList <- c('p1x.r','p1a.r','p1b.r','nx.r','px.r','n1.r','p2.r','n2.r')
cmpStrList <- paste('z.',cmpStrList,sep='')

if (length(args)>=2)
  dataSet <- as.character(args[2])

if (length(args)>=3)
  speaker <- as.character(args[3])

useInd <- which(df$speaker == speaker & df$dataSet==dataSet & df$animal==animal)

if(length(args)>=4)
  useInd <- useInd[ length(useInd) ]

useInd
tdf <- df[useInd, ]

##chStr       <- 'frontoCentral'  
chStr       <- NA
tone        <- loadAmpOddClickAll( tdf , what='raw', chStr=chStr, extn=extn  )
tone$exp    <- paste(animal,'_',dataSet,extn,sep='')

minISI <- 0.190
maxISI <- 12.8
DVscl  <- 5


show <- FALSE
for ( i in 1:length(cmpStrList) ){
  sfc <- TRUE
  
  ## ====================================  Short-term post-synaptic plasticiy  --  limited resource model
  disc <- callSTPSP(         tone, gd=T, startFromScratch=sfc, DVstr=cmpStrList[i],who='all',BYstr='exp',     
                             FUN=cumsum1,minISI=minISI, maxISI=maxISI,DVsdscale=DVscl, show=show )
  disc <- callSTPSP(         tone, gd=T, startFromScratch=sfc, DVstr=cmpStrList[i],who='random',BYstr='exp', 
                             FUN=cumsum1,minISI=minISI, maxISI=maxISI,DVsdscale=DVscl, show=show )
  dsc  <- callSTPSP(         tone, gd=T, startFromScratch=sfc, DVstr=cmpStrList[i],who='regular',BYstr='exp', 
                             FUN=cumsum1,minISI=minISI, maxISI=maxISI,DVsdscale=DVscl, show=show )

  disc <- callSTPSP(         tone, gd=T, startFromScratch=sfc, DVstr=cmpStrList[i],who='all',BYstr='exp',     
                             FUN=cumsumscl,minISI=minISI, maxISI=maxISI,DVsdscale=DVscl, show=show )
  disc <- callSTPSP(         tone, gd=T, startFromScratch=sfc, DVstr=cmpStrList[i],who='random',BYstr='exp', 
                             FUN=cumsumscl,minISI=minISI, maxISI=maxISI,DVsdscale=DVscl, show=show )
  dsc  <- callSTPSP(         tone, gd=T, startFromScratch=sfc, DVstr=cmpStrList[i],who='regular',BYstr='exp', 
                             FUN=cumsumscl,minISI=minISI, maxISI=maxISI,DVsdscale=DVscl, show=show )
  
  
  disc <- callSTPSP(         tone, gd=T, startFromScratch=sfc, DVstr=cmpStrList[i],who='all',BYstr='exp',     
                             FUN=cumsumamp,minISI=minISI, maxISI=maxISI,DVsdscale=DVscl, show=show )
  disc <- callSTPSP(         tone, gd=T, startFromScratch=sfc, DVstr=cmpStrList[i],who='random',BYstr='exp', 
                             FUN=cumsumamp,minISI=minISI, maxISI=maxISI,DVsdscale=DVscl, show=show )
  dsc  <- callSTPSP(         tone, gd=T, startFromScratch=sfc, DVstr=cmpStrList[i],who='regular',BYstr='exp', 
                             FUN=cumsumamp,minISI=minISI, maxISI=maxISI,DVsdscale=DVscl, show=show )
  
  
  #disc <- callSTPSPfixRelease(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='all', show=show )
  #disc <- callSTPSP(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='random', show=show )
  #disc <- callSTPSP(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='regular', show=show )
  
  #disc <- callHistoryOneTone(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='all', show=show )
  #disc <- callHistoryOneTone(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='random', show=show )
  #disc <- callHistoryOneTone(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='regular', show=show )

  
  ## ====================================== ARtau
  #disc <- callARtau(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='all',     show=show, minISI=minISI, maxISI=maxISI,DVsdscale=DVscl )
  #disc <- callARtau(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='random',  show=show, minISI=minISI, maxISI=maxISI,DVsdscale=DVscl  )
  #disc <- callARtau(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='regular', show=show, minISI=minISI, maxISI=maxISI,DVsdscale=DVscl  )
  
  ## ---------- by frequency bin
  #disc <- callARtau(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='regular', frqStr='low',  show=show,minISI=minISI, maxISI=maxISI,DVsdscale=DVscl  )
  #disc <- callARtau(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='regular', frqStr='mid',  show=show, minISI=minISI, maxISI=maxISI,DVsdscale=DVscl  )
  #disc <- callARtau(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='regular', frqStr='high', show=show, minISI=minISI, maxISI=maxISI,DVsdscale=DVscl  )
  
  #disc <- callARtau(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='random', frqStr='low',  show=show,minISI=minISI, maxISI=maxISI,DVsdscale=DVscl  )
  #disc <- callARtau(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='random', frqStr='mid',  show=show, minISI=minISI, maxISI=maxISI,DVsdscale=DVscl  )
  #disc <- callARtau(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='random', frqStr='high', show=show, minISI=minISI, maxISI=maxISI,DVsdscale=DVscl  )
  
  #disc <- callARtau(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='all', frqStr='low',  show=show,minISI=minISI, maxISI=maxISI,DVsdscale=DVscl  )
  #disc <- callARtau(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='all', frqStr='mid',  show=show, minISI=minISI, maxISI=maxISI,DVsdscale=DVscl  )
  #disc <- callARtau(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='all', frqStr='high', show=show, minISI=minISI, maxISI=maxISI,DVsdscale=DVscl  )
  
}


#for ( i in 1:length(cmpStrList) )
#  disc <- callHistoryOneTone(tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='random')

#for ( i in 1:length(cmpStrList) )
#  disc <- callHistoryOneTone(tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='regular')

# for ( i in 1:length(cmpStrList) ){
#   disc <- callHistoryOneTone(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='all' )
#   disc <- callHistoryOneToneOneLambda(tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='all')
# }
# 
# 
# for ( i in 1:length(cmpStrList) ){
#   disc <- callHistoryOneTone(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='random' )
#   disc <- callHistoryOneToneOneLambda(tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='random')
#   disc <- callHistoryOneTone(         tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='regular' )
#   disc <- callHistoryOneToneOneLambda(tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='regular')
# }

#for ( i in 1:length(cmpStrList) )
#  disc <- callHistoryOneToneOneLambda(tone, gd=T, startFromScratch=T, DVstr=cmpStrList[i],who='random')






