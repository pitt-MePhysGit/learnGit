## define the paths ====
repoDir <- '~/pitt-MePhysGit'
baseDir <- '~/Dropbox'

print(paste("Using ", repoDir ,' repository.', sep=''))

## pull the latest code from git  ====
print(' Local branch: ')
system(paste( 'git -C ',repoDir,'/learnGit/ branch',sep='') )


print(' Pull the latest from the local branch ')
system(paste( 'git -C ',repoDir,'/learnGit/ pull',sep='') )

source(paste(baseDir,'/r-compile/rUtils/R/rUtils.R',sep=''))




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
