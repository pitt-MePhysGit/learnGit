## Rsetup.R                                                               ====
# 
# Rsetup.R is run prior to any R-based analysis.
# it defines the location of the code to be sourced
# options are Dropbox or pitt-MePhysGit
# for code on the git server, you can specify which branch to use
# prior to sourcing the code, it switches to the specified branch
# and performs a pull request to get the latest code.

## define the paths                                                       ====
repoName <- 'pitt-MePhysGit'
repoDir  <- paste('~/',repoName, sep='')
baseDir  <- '~/Dropbox'

print(paste("Using ", repoDir ,' repository.', sep=''))


# list all projects with source-code either on Dropbox or pitt-MePhysGit) ====
projectList <- list('activeTRF','ampOddClick','amRF', 'chait','changeDetection','clickGroup','cyyclicNoise2','detectTone',
                    'DRF','DTRF','durOdd','dynamicChords','eCog','EFR','ephysLab','harmonicity','learnGit','MFF','narrowBandNoise',
                    'noiseTRF','reactivationTRF','resting','rewardAnticipation','RFmap','RSkernel','SMART','SSA','ssrt',
                    'STRF','stroop','sweepRF','tDCS','temporalDiscounting','textures','TRF','vocalizations','whiteNoiseDetection')

# list all projects that have been transferred to pitt-MePhysGit          ====
repoList      <- list('amRF','learnGit')
projectFrame  <- data.frame( name = unlist(projectList), dir = '~/Dropbox/', branch = 'NA')

projectFrame$dir[    which(projectFrame$name %in% repoList) ] <- repoName
projectFrame$branch[ which(projectFrame$name %in% repoList) ] <- 'main'


# adjust branches individually                                            ====

# template:
projectFrame$branch[which(projectFrame$name=='learnGit')] <- 'kifBranch'



## pull the latest code from git   ==== 
rndx <- which(projectFrame$dir == repoName)
for (ix in rndx ){
 
  print(' ========================== ')
  print( paste('Project: ', projectFrame$name[ix], ' - branch: ', projectFrame$branch[ix], sep='') )
  
  print('switch to correct branch')
  system(paste( 'git -C ', repoDir,'/', projectFrame$name[ix],'/ switch ', projectFrame$branch[ix] ,sep='') )
  
  print(' Local branch: ')
  system(paste( 'git -C ',repoDir,'/', projectFrame$name[ix],'/ branch',sep='') )
  
  print(' Pull the latest from the local branch ')
  system(paste( 'git -C ',repoDir,'/',projectFrame$name[ix],'/ pull',sep='') )
  print('>>>>>>>>>>>>>>>>>>>>>>>>>>>')
  
}

source(paste(baseDir,'/r-compile/rUtils/R/rUtils.R',sep=''))




## load libraries                  ====
library(R.matlab)
library(lattice)#for eps.file
try(library(car))   # was commented out... not sure why. necessary for averagR; can't be installed on tinyone
library(signal)
library(akima)# for interp
library(Rwave)# for cwt
library(scales)




## plot something                  ====
plot(1:20, pch=1:20, col=1:20)

## how about a different plot 
plot( 1:10, (1:10)/2 )
