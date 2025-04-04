## =======================================
## Example script 
## =======================================

source('~/Dropbox/r-compile/signalProcessoR/R/signalProcessoR.R')

project     <- 'ampOddClick'
linkFile    <- paste(project,'EEG_link.txt',sep='')
linkFile    <- paste(project,'_link.txt',sep='')
projectDir  <- paste('~/Dropbox/',project,'/',sep='')


df          <- prepDFsignalProcessor(project, linkFile=linkFile)

table( df$animal, df$dataSet)
table( df$animal, df$year)


##Ben's first session doesn't load on tinyone, fine on Farnsworth; need to track down why
Bdf <- df[ (df$animal=='Ben'    & df$year=='2020' & df$dataSet=='g'),  ][c(2,3,4,5,6,8,9,10,11),]
Jdf <- df[ (df$animal=='Jesse'  & df$year=='2016' & df$dataSet=='g'),  ][c(1,2,4,5,6,7,11,12,13,14,15),] # 8 is the inverted one
Sdf <- df[ (df$animal=='Sam'    & df$year=='2016' & df$dataSet=='g'),  ]
Rdf <- df[ (df$animal=='Rockey' & df$year=='2016' & df$dataSet=='g'),  ]
Wdf <- df[ (df$animal=='Walter' & df$year=='2016' & df$dataSet=='g'),  ]



#Rscript ~/Dropbox/r-compile/signalProcessoR/R/batchCall_FLAS.R -p ampOddClick -a Walter 20010001 20310001 -dlist 
#paste(Wdf$date, collapse = "+")



Adf <- rbind(Bdf, Jdf, Sdf, Rdf, Wdf)

tdf        <- Jdf
eCogFile   <- 'jessesym'


epoch      <- 'raw'
analysis   <- 'ampInd'
chStr      <- paste('ch',1:33,sep='')
#session      <- loadSignalProcessor(  tdf[2,],chStr='ch4', project=project, epoch=epoch, blrng=c(-25,0) )  



## ======================================== one animal all sessions, all channels 

raw      <- loadSignalProcessor(  tdf,chStr=chStr, project=project, epoch='raw', analysis=analysis, blrng=c(-25,0) )  
bsp      <- loadSignalProcessor(  tdf,chStr=chStr, project=project, epoch='bsp', analysis=analysis, blrng=c(-25,0) )  
llp      <- loadSignalProcessor(  tdf,chStr=chStr, project=project, epoch='llp', analysis=analysis, blrng=c(-25,0) )  

raw$channelPosition <- getEEGpos(eCogFile)
bsp$channelPosition <- getEEGpos(eCogFile)
llp$channelPosition <- getEEGpos(eCogFile)

showSignalProcessor( raw, frml=~sessionID, trng=c(-25,250), sclfctr = 1 )
showSignalProcessor( llp, frml=~sessionID, trng=c(-25,250), sclfctr = .5 )
showSignalProcessor( bsp, frml=~ampInd, trng=c(-1,7), sclfctr = .5 )


epoch <- 'llp'
oneB      <- loadSignalProcessor(  Bdf,chStr='ch4', project=project, epoch=epoch, analysis=analysis, blrng=c(-25,0) )  
oneJ      <- loadSignalProcessor(  Jdf,chStr='ch4', project=project, epoch=epoch, analysis=analysis, blrng=c(-25,0) )  
oneS      <- loadSignalProcessor(  Sdf,chStr='ch4', project=project, epoch=epoch, analysis=analysis, blrng=c(-25,0) )  
oneR      <- loadSignalProcessor(  Rdf,chStr='ch4', project=project, epoch=epoch, analysis=analysis, blrng=c(-25,0) )  
oneW      <- loadSignalProcessor(  Wdf,chStr='ch4', project=project, epoch=epoch, analysis=analysis, blrng=c(-25,0) )  

showSignalProcessorXYC( oneB, xaxis='sessionID',yaxis='dummy',caxis='ampInd', trng=c(-25,150),tsqueezrng=c(0,0.80),colFUN=soa.colors)
showSignalProcessorXYC( oneJ, xaxis='sessionID',yaxis='dummy',caxis='ampInd', trng=c(-25,150),tsqueezrng=c(0,0.80),colFUN=soa.colors)
showSignalProcessorXYC( oneS, xaxis='sessionID',yaxis='dummy',caxis='ampInd', trng=c(-25,150),tsqueezrng=c(0,0.80),colFUN=soa.colors)
showSignalProcessorXYC( oneR, xaxis='sessionID',yaxis='dummy',caxis='ampInd', trng=c(-25,150),tsqueezrng=c(0,0.80),colFUN=soa.colors)
showSignalProcessorXYC( oneW, xaxis='sessionID',yaxis='dummy',caxis='ampInd', trng=c(-25,150),tsqueezrng=c(0,0.80),colFUN=soa.colors)





epoch <- 'llp'
oneA              <- loadSignalProcessor(  Adf,chStr='ch5', project=project, epoch=epoch, analysis='ampInd_iciInd', blrng=c(-25,0) ) 
oneA$xpp$animalID <- c(as.factor(oneA$xpp$animal))  
avgA              <- avgSignalProcessor(oneA, frml=~animalID+ampInd+iciInd)
showSignalProcessorXYC( avgA, xaxis='animalID',yaxis='ampInd',caxis='iciInd', trng=c(-25,150),tsqueezrng=c(0,0.80),colFUN=soa.colors )

showSignalProcessorXYC( avgA, xaxis='ampInd',yaxis='iciInd',caxis='animalID', trng=c(-25,150),tsqueezrng=c(0,0.80) )


avgAA              <- avgSignalProcessor(oneA, frml=~ampInd+iciInd)
showSignalProcessorXYC( avgAA, xaxis='ampInd',yaxis='dummy',caxis='iciInd', trng=c(-25,150),tsqueezrng=c(0,0.80),colFUN=soa.colors )
showSignalProcessorXYC( avgAA, xaxis='iciInd',yaxis='dummy',caxis='ampInd', trng=c(-25,150),tsqueezrng=c(0,0.80),colFUN=bluered.colors )
showSignalProcessorXYC( avgAA, xaxis='ampInd',yaxis='iciInd',caxis='iciInd', trng=c(-25,150),tsqueezrng=c(0,0.80),colFUN=soa.colors )


oneA              <- loadSignalProcessor(  Adf,chStr='ch5', project=project, epoch=epoch, analysis='iciInd', blrng=c(-25,0) ) 
oneA$xpp$animalID <- c(as.factor(oneA$xpp$animal))  
avgA              <- avgSignalProcessor(oneA, frml=~animalID+iciInd)
showSignalProcessorXYC( avgA, xaxis='animalID',yaxis='dummy',caxis='iciInd', trng=c(-25,150),tsqueezrng=c(0,0.80) )







showSignalProcessor( raw, frml=~ampInd, trng=c(-25,250), sclfctr = 1 )
showSignalProcessor( click, frml=~sessionID, trng=c(-25,250), sclfctr = 1 )


showSignalProcessorXYC(onech, xaxis='sessionID', yaxis='ypos', caxis='sessionID', trng=c(-10,150),tsqueezrng=c(0,0.80))









