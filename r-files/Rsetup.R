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
projectList <- c('activeTRF','ampOddClick','amRF', 'chait','changeDetection','clickGroup','cyyclicNoise2','detectTone',
                 'DRF','DTRF','durOdd','dynamicChords','eCog','EFR','ephysLab','harmonicity','MFF','narrowBandNoise',
                  'noiseTRF','reactivationTRF','resting','rewardAnticipation','RFmap','RSkernel','SMART','SSA','ssrt',
                  'STRF','stroop','sweepRF','tDCS','temporalDiscounting','textures','TRF','vocalizations','whiteNoiseDetection')

# list all projects that have been transferred to pitt-MePhysGit          ====
repoList      <- c('amRF','learnGit')
projectFrame  <- data.frame( name = unlist(projectList), dir = '~/Dropbox/', branch = 'NA', git=F)

projectFrame$git[    which(projectFrame$name %in% repoList) ] <- T
projectFrame$dir[    which(projectFrame$name %in% repoList) ] <- repoName
projectFrame$branch[ which(projectFrame$name %in% repoList) ] <- 'main'


# adjust branches individually                                            ====

# template:
# projectFrame$branch[which(projectFrame$name=='learnGit')] <- 'kifBranch'


## handle toolboxes                                                       ====
compileList <- list('deadlineModel','ePhysLab','eRnii','fokkerPlank','mav.test','MePhys','mph.test','psyfit','rgenoud',
                    'RPE2PSEShift','rUtils','sdt','SignalProcessor','stereotax','variableDriftrateToolBox',
                    'visualizR','wieneR')

repoList     <- c(unlist(repoList), 'deadlineModel')
compileFrame <- data.frame(name = unlist(compileList), dir = '~/Dropbox/toolbox/', branch = 'NA', git = F)

compileFrame$git[    which(compileFrame$name %in% repoList) ] <- T
compileFrame$dir[    which(compileFrame$name %in% repoList) ] <- paste(repoName,'/toolbox/',sep='')
compileFrame$branch[ which(compileFrame$name %in% repoList) ] <- 'main'


# combine the two frames:
projectFrame <- rbind(projectFrame, compileFrame)

## pull the latest code from git   ==== 
rndx <- which(projectFrame$git == T)
# list the projects in git repository
projectFrame[rndx,]

for (ix in rndx ){
 
  print(' >>>>>>>>>>>>>>>>>>>>>>>>>>> ')
  print( paste('Project: ', projectFrame$name[ix], ' - branch: ', projectFrame$branch[ix], sep='') )
  
  print('switch to correct branch')
  system(paste( 'git -C ', repoDir,'/', projectFrame$name[ix],'/ switch ', projectFrame$branch[ix] ,sep='') )
  
  print(' Local branch: ')
  system(paste( 'git -C ',repoDir,'/', projectFrame$name[ix],'/ branch',sep='') )
  
  print(' Pull the latest from the local branch ')
  system(paste( 'git -C ',repoDir,'/',projectFrame$name[ix],'/ pull',sep='') )
  print(' <<<<<<<<<<<<<<<<<<<<<<<<<<<<  ')
  
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

## source functions                ====
# this is where the list of files to source would go:
# test if I can commit via ssh from kif

