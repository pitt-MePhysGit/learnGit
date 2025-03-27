## pull the latest code from git  ====
system( 'git -C ~/pitt-MePhysDir/learnGit/ pull' )

## define the paths ====
repoDir <- '~/pitt-MePhysGit'
baseDir <- '~/Dropbox'

source(repoDir,'/r-compile/rUtils/R/rUtils.R')
source(repoDir,'/r-compile/ePhysLab/R/ePhysLab.R')
source(repoDir,'/r-compile/mav.test/R/mav.test.R')
source(repoDir,'/r-compile/stereotax/R/stereotax.R')


## load libraries ====
library(R.matlab)
library(lattice)#for eps.file
try(library(car))   # was commented out... not sure why. necessary for averagR; can't be installed on tinyone
library(signal)
library(akima)# for interp
library(Rwave)# for cwt
library(scales)




## plot something ====
plot(1:20, pch=1:20, col=1:20)

## how about a different plot 
plot( 1:10, (1:10)/2 )
