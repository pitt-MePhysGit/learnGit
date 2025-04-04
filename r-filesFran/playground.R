#
# model human data

# /Volumes/EEGCyno/EEG/ampOddClickHuman/
# The path is ~/Dropbox/ampOddClick/r-filesFran.

# Fits should be with the entire subject (all sessions together)

# Read
########################################################################
########################################################################
# Tested that the sequence of sounds within a block is the FULL sequence, 
# only without EEG signal in the rows where data was not analyzed (but the 
# "history of stimulation" still holds). Fixed weird trigger codes 
# (full log of corrections in scripts in brainstorm_pipelines folder in Y drive).

# EEG scalp: In 2022 folder I created new ones with all the subjects

# Sources:
# Y:/HuMon/Events/Single_trials/Sources/Final_Single_trial_source_table_MEG_LLR.mat
# Y:/HuMon/Events/Single_trials/Sources/Final_Single_trial_source_table_MEG_MLR.mat
# These are with the final set of subjects (30 for LLR, 29 for MLR because of one very 
# bad discarded). MLR has the 100Hz low pass filter.

# Before analyzing the sources data set, Tobias preferred to use a fully 
# orientation-constrained (orthogonal to cortex) source rather than the 
# normalized across three orientations. That implies re-computing the entire 
# source solution for HuMon which may take months of computing time alone, 
# so better work with the current data sets for now. To make it simpler, 
# I replaced the entire waveform in those with the average time windows 
# for each peak we used in the EEG analysis. Those files are the ones with 
# scout names on them. 

# There are rare (0.005% of the data) where ISI were above 8s. While that 
# should not happen, I checked at the logs manually and there are indeed 
# instances in which that happened during the recordings, for reasons 
#I don't understand. I decided to replace those by NaNs, as the ones starting 
# a block or a session, which should be interpreted by the model as a "reset" 
# where neuronal resources were replenished. 
# So: any sequence containing NaNs in ISI should be broken down into smaller ones starting at each NaN
########################################################################

# Install packages and define variables
rm(list=ls())
install.packages("gdata")
install.packages("Rcpp")
install.packages("RcppArmadillo")
library(lattice)
baseDir <- '~/Dropbox/'
source( paste(baseDir,'r-compile/rUtils/R/rUtils.R',sep='') )
source( paste(baseDir,'r-compile/ePhysLab/R/ePhysLab.R',sep='') )
source('~/Dropbox/ampOddClick/r-files/functions/functions_ampOdd_presynapticPlasticity.R')
project   <- 'ampOddClick'
fileNameM <- 'Amp_singl_trials_EEG_MLR.txt'
fileNameL <- 'Amp_singl_trials_EEG_LLR.txt'


# /Volumes/EEGCyno/EEG/ampOddClickHuman/
# Load datasets (EEG scalp data)
LLR <- read.table( paste('~/Dropbox/',project,'/RData/',fileNameL,sep=''), header=T,sep = "," )
MLR <- read.table( paste('~/Dropbox/',project,'/RData/',fileNameM,sep=''), header=T,sep = "," )
# /Volumes/EEGCyno/EEG/ampOddClickHuman
# Replace NaN at ISI column by a large number
# One long ISI sequence per subject

# /Volumes/EEGCyno/EEG/ampOddClickHuman/

# Exploring datasets
table(LLR$subjID)
table(LLR$subjID,LLR$sessionNr)
table(LLR$sessionNr, LLR$blockNr,LLR$subjID)
table(MLR$subjID)
table(MLR$subjID,LLR$sessionNr)
table(MLR$sessionNr, LLR$blockNr,LLR$subjID)
names(LLR)
names(MLR)
LLR[15815,]

table(LLR$clickIntensity)
which(LLR$clickIntensity == 8)

# onsets <- c(  1,2,2.5,2.75,3,4,6,9,10,11  )
# cumsumamp is the function for Tobia's model. Inputs are as follows based on example:
# Example: plot(onsets[-1], cumsumamp( diff(onsets), rep(.5,9), 1/.5, .9 ), ylim=c(0,1) )
# Input 1) diff(onsets) is the list of ISIs
# Input 2) rep(.5,9) is the list of intensities for each click: between 0 and 1
# Input 3) 1/.5 is recovery time constant. This example means .5 seconds to recover half of the resources
# Input 4) .9 is cumulative release probability. Between 0 and 1. High means a lot of resources invested in responding to the sound
# The second argument is the list of intensities for each click. 
# The second and third parameter are the key model parameters: 
# recovery time constant and cumulative release probability


abs_time <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
ISI_list <- rep(.02,9)
ci <- 1 # Click intensity: scale must correspond to the ylim?
int_every_click <- rep(ci,length(ISI_list))
rt <- 1/.5 # recovery time Original was 1/.5
plot(abs_time, cumsumamp( ISI_list, int_every_click, rt, 1 ), ylim=c(0,1) )

# differences are squared and sum to get the fit value
# Try the first subject first


