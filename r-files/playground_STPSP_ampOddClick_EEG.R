## preliminaries
#args <- commandArgs( trailingOnly=TRUE )

baseDir <- '~/Dropbox/'
source(paste(baseDir, 'ampOdd_click/r-files/source_all_functions.R',sep=''))

## preliminaries
pp  <- function(){paste(baseDir, 'ampOdd_click/plots/g/',sep='')}
rda <- function(){'~/Dropbox/ampOdd_click/Rdata/'}


animalList <- c('Rockey','Jesse','Walter','Sam')
dataSet    <- 'g'
cmpStrList <- c('p1x','p1a','p1b','nx','px','n1','p2','n2')
#cmpStrList <- c('p1x.r','p1a.r','p1b.r','nx.r','px.r','n1.r','p2.r','n2.r')
#cmpStrList <- paste('z.',cmpStrList,sep='')

maxp2p     <- 450
minISI     <- 0.190
maxISI     <- 12.8
DVsdscale  <- 5


chStr       <- NA
getData <- function(subj,DVstr=c('p1a','p1b','nx','px','n1','p2','n2')){
  print(subj)
  tone        <- loadAmpOddClickAll( df[which(df$animal==subj&df$dataSet==dataSet),] , what='raw', chStr=chStr  )
  tone$exp    <- paste(subj,'_',dataSet,sep='')
  
  for (dx in 1:length(DVstr)){
    tone$xpp$DV <- tone$xpp[[DVstr[dx] ]]
    
    testLambda <- function(lmd){  
      print(lmd)
      tone$xpp$oKer  <-        exp(  lmd*tone$xpp$ISI )
      lmExp   <- lm(tone$xpp$DV ~ 1 +  tone$xpp$oKer  )
      return( sum(lmExp$residuals^2) )
    }
    
    resExp     <- optimize(f=testLambda, lower=-20, upper=-0.1, maximum = FALSE)
    lmExp      <- lm(tone$xpp$DV ~ 1  +  I( exp(resExp$minimum * tone$xpp$ISI)  ))
    
    mxAmp <- lmExp$coefficients[1]
    tone$xpp[[ DVstr[dx] ]] <- tone$xpp$DV/mxAmp
  }
  
  return(tone)
}
dataList <- lapply( animalList, getData)

tone <- append.data( append.data( dataList[[1]], dataList[[2]]), append.data(dataList[[3]], dataList[[4]]) )
tone$exp <- 'g'

## no scaling with amplitude at all
cumFit <- callSTPSP(tone, DVstr='p1b', gd=T, show=TRUE, showCh=1, who='all', BYstr='dummy',
                    FUN=cumsum1, minISI=.190, maxISI=12.8, startFromScratch=T)
cumFit$modelFit$par


## scale R/EEGamp with tone amplitude
sclFit <- callSTPSP(tone,                 DVstr='n1', gd=T, show=TRUE, showCh=1, who='random', BYstr='dummy', 
                    FUN=cumsumscl, minISI=.190, maxISI=12.8, startFromScratch=T)
sclFit$modelFit$par


## scale expended Fraction with tone amplitude
ampFit <- callSTPSP(tone,                 DVstr='n1', gd=T, show=TRUE, showCh=1, who='random', BYstr='dummy', 
                    FUN=cumsumamp, minISI=.190, maxISI=12.8, startFromScratch=F, eps=T)
ampFit$modelFit$par


#c(cumFit$modelFit$par,cumFit$modelFit$value)
c(sclFit$modelFit$par,sclFit$modelFit$value)
c(ampFit$modelFit$par,ampFit$modelFit$value)




