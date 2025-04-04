.pardefault <- par(no.readonly = T)

source( paste(baseDir,'r-compile/rUtils/R/rUtils.R',sep='') )
source( paste(baseDir,'r-compile/ePhysLab/R/ePhysLab.R',sep='') )
source( paste(baseDir,'r-compile/mph.test/R/mph.test.R',sep='') )

source( paste(baseDir,'aMMN/aMMN_controlled/r-files/functions_RS_kernel.R',sep='') )

#source( paste(baseDir,'aMMN_controlled/r-files/functions_model.R',sep='') )
#source( paste(baseDir,'aMMN_controlled/r-files/functions_RS_modelFit.R',sep='') )
#source( paste(baseDir,'aMMN_controlled/r-files/functions_RS_oneTone.R',sep='') )
#source( paste(baseDir,'aMMN_controlled/r-files/functions_RS_oneTone_oneLambda.R',sep='') )

## some stuff from RSkernel
source( paste(baseDir,'aMMN_controlled_temp/r-files/functions_eCog.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/figures/call_multiplot.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/componentColor.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/componentColor2.R',sep='') )
source( paste(baseDir,'RSkernel/r-files/functions/functions_master_STPSP.R',sep='') )
# stuff not needed right now
#source( paste(baseDir,'RSkernel/r-files/functions/calcCSD.R',sep='') )
#source( paste(baseDir,'RSkernel/r-files/functions/showCSDimage.R',sep='') )
#source( paste(baseDir,'RSkernel/r-files/functions/showCSD.R',sep='') )
#source( paste(baseDir,'RSkernel/r-files/multiVplot.R',sep='') )
#source( paste(baseDir,'RSkernel/r-files/functions/functions_ketamine.R',sep='') )
#source( paste(baseDir,'RSkernel/r-files/functions/functions_model_ARtau.R',sep='') )


sourceDir(paste(baseDir,'ampOddClick/r-files/functions/',sep=''))
sourceDir(paste(baseDir,'ampOddClick/r-files/figures/',sep=''))
source( paste(baseDir,'ampOddClick/r-files/parameters/getCompPar.R',sep='') )


library(rgenoud)
library(R.matlab)
library(lattice)
library(akima)
library(Rwave)
library(MASS)

##rda <- function(){'~/Documents/RSkernel/rda/'}
rda <- function(){'~/Dropbox/ampOddClick/Rdata/'}

colClasses = c('character',rep('numeric',8),rep('character',2))
df   <- read.table(paste(baseDir,'ampOddClick/ampOddclick_link.txt',sep=''),header=T,colClasses=colClasses)
df$animal <- unlist( lapply(1:dim(df)[1], function(n){    animal <- strsplit(df$session[n],'_')[[1]][1] }  ))


