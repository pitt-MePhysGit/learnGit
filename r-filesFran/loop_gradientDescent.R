
rm(list=ls())
##  Notes -------------------------
# playground_gradientDescent
# CALLS load_ampOddClick_Human to define parameters than one in turns calls
# callSTPSP does the gradient descent

# To see results
# parm <- c(ouput of last line of gradient descent)

# All of them are in  Home/Dropbox/ampOdClick/rfiles/functions/pr)
# myEdit(any function in here opens it)
########################################################
## Add packages -----------------------
library(pracma)
library(lattice)
library(ggplot2)
#library(corrplot)
rda <- function(){ return('/Volumes/EEGCyno/RES/ampOddClick/') }
source('~/Dropbox/r-compile/rUtils/R/rUtils.R')
source('~/Dropbox/r-compile/ePhysLab/R/ePhysLab.R')
source('~/Dropbox/ampOddClick/r-files/functions/functions_ampOdd_presynapticPlasticity.R')
source('~/Dropbox/ampOddClick/r-filesFran/load_ampOddClick_Human.R')

## Define variables -------------------------
Wave_list <- c('P0','Na','Pa','Nb','P50','N1','P2')
Type_signal_list <- c('EEG','Source')
Hemisphere_list <- c('L','R') # L or R, only if choosing sources as type of signal
Parcel_list <- c('A1','MBelt','LBelt','PBelt') # Only if choosing sources as type of signal
Subject_list <- c('2193','2199','2235','2238','2250','2320','2259','2262','2274',
             '2275','H001','2214','2212','2285','2337','2372','2389','H002',
             '2299','2406','2426','2494','2427','2456','H005',
             '2005','2355','2387','2288') ## 

## Loop through subjects/Waves/Sources ----------------------------
Header_mega_variable <-c('EEG_P0','EEG_Na','EEG_Pa','EEG_Nb','EEG_P50','EEG_N1','EEG_P2',
                         'A1_L_P0','A1_L_Na','A1_L_Pa','A1_L_Nb','A1_L_P50','A1_L_N1','A1_L_P2',
                         'A1_R_P0','A1_R_Na','A1_R_Pa','A1_R_Nb','A1_R_P50','A1_R_N1','A1_R_P2',
                         'MBelt_L_P0','MBelt_L_Na','MBelt_L_Pa','MBelt_L_Nb','MBelt_L_P50','MBelt_L_N1','MBelt_L_P2',
                         'MBelt_R_P0','MBelt_R_Na','MBelt_R_Pa','MBelt_R_Nb','MBelt_R_P50','MBelt_R_N1','MBelt_R_P2',
                         'LBelt_L_P0','LBelt_L_Na','LBelt_L_Pa','LBelt_L_Nb','LBelt_L_P50','LBelt_L_N1','LBelt_L_P2',
                         'LBelt_R_P0','LBelt_R_Na','LBelt_R_Pa','LBelt_R_Nb','LBelt_R_P50','LBelt_R_N1','LBelt_R_P2',
                         'PBelt_L_P0','PBelt_L_Na','PBelt_L_Pa','PBelt_L_Nb','PBelt_L_P50','PBelt_L_N1','PBelt_L_P2',
                         'PBelt_R_P0','PBelt_R_Na','PBelt_R_Pa','PBelt_R_Nb','PBelt_R_P50','PBelt_R_N1','PBelt_R_P2')
# Load variable to store values
load(paste("~/Dropbox/ampOddClick/r-filesFran/Mega_variable_param",sep=""))
# load(paste("~/Dropbox/ampOddClick/r-filesFran/Mega_variable_fits",sep=""))
load(paste("~/Dropbox/ampOddClick/r-filesFran/Mega_variable_corr",sep=""))

# Mega_variable <- matrix(list(), length(Subject_list), length(Header_mega_variable))
# colnames(Mega_variable) <- Header_mega_variable;
# rownames(Mega_variable) <- Subject_list;

# Define parameters
for (subid in 1:length(Subject_list)) {
  Subject <- Subject_list[subid]
  # Do EEG alone first
  Type_signal <- 'EEG'
  for (w in 1:length(Wave_list)) {
    Wave <- Wave_list[w]
    Hemisphere <- 'NA' # Irrelevant here
    Parcel <- 'NA' # Irrelevant here

    # Run gradient descent if data not in Mega variable (either in Param or fits, as it will be the same)
    if(length(Mega_variable_param[[subid,w]]) == 0) {
      tn <- load_ampOddClick_Human(Wave, Subject, Type_signal, Hemisphere,Parcel)
      tn$xpp$dummy = 1
      print(paste('RUNNING GRADIENT DESCENT FOR ', Subject, Wave, Type_signal, sep = ' '))
      tst <- callSTPSP(tn, DVstr='DV', gd=T, show=TRUE, showCh=1,FUN=ampFUN2q,maxISI=10, minISI=0.2, maxp2p=5000, DVsdscale=5, 
                       maxEOG=Inf, trialType=NA, fitType=NA, manStartPar = NULL, robust = FALSE, startFromScratch = T,
                       fileExtn = paste("p2q",Subject, Parcel, Hemisphere, Wave, sep = '_'),BYstr='dummy',topInt = 1)
      
      tst <- callSTPSP(tn, DVstr='DV', gd=T, show=TRUE, showCh=1,FUN=ampFUN,maxISI=10, minISI=0.2, maxp2p=5000, DVsdscale=5, 
                       maxEOG=Inf, trialType=NA, fitType=NA, manStartPar = NULL, robust = FALSE, startFromScratch = T,
                       fileExtn = paste("p",Subject, Parcel, Hemisphere, Wave, sep = '_'))
      # Store values of best fit in data.frame
      Mega_variable_param[[subid,w]]<- tst[['modelFit']][1]
      Mega_variable_corr[[subid,w]]<- tst[['modelFit']][2]
      Mega_variable_fits[[subid,w]]<- tst[["res"]][["lm0"]][["model"]][,1]
      # Save Mega_variables
      save(Mega_variable_param, file = paste("~/Dropbox/ampOddClick/r-filesFran/Mega_variable_param",sep=""))
      save(Mega_variable_corr, file = paste("~/Dropbox/ampOddClick/r-filesFran/Mega_variable_corr",sep=""))
      save(Mega_variable_fits, file = paste("~/Dropbox/ampOddClick/r-filesFran/Mega_variable_fits",sep=""))
      print(paste('FITTING COMPLETED FOR ', Subject, Wave, Type_signal, sep = ' '))
    } else {
      next
    }

  }
  # Then do Sources
  Type_signal <- 'Source'
  pos <- length(Wave_list) + 1
  for (parc in 1:length(Parcel_list)) {
    Parcel <- Parcel_list[parc]
    for (hem in 1:length(Hemisphere_list)) {
      Hemisphere <- Hemisphere_list[hem]
      for (w in 1:length(Wave_list)) {
        Wave <- Wave_list[w]
        # Run gradient descent if data not in Mega variable (either in Param or fits, as it will be the same)
        if(length(Mega_variable_param[[subid,pos]]) == 0) {
          tn <- load_ampOddClick_Human(Wave, Subject, Type_signal, Hemisphere,Parcel) 
          print(paste('RUNNING GRADIENT DESCENT FOR ', Subject, Parcel, Hemisphere, Wave, sep = ' '))
          tst <- callSTPSP(tn, DVstr='DV', gd=T, show=TRUE, showCh=1,FUN=ampFUN,maxISI=10, minISI=0.2, maxp2p=5000, DVsdscale=5,
                           maxEOG=Inf, trialType=NA, fitType=NA, manStartPar = NULL, robust = FALSE, 
                           fileExtn = paste("p", Subject, Parcel, Hemisphere, Wave, sep = '_'))
          # Store values of best fit in data.frame
          Mega_variable_param[[subid,pos]] <- tst[['modelFit']][1]
          Mega_variable_corr[[subid,pos]] <- tst[['modelFit']][2]
          Mega_variable_fits[[subid,pos]]<- tst[["res"]][["lm0"]][["model"]][,1]
          pos <- pos + 1
          # Save Mega_variables
          save(Mega_variable_param, file = paste("~/Dropbox/ampOddClick/r-filesFran/Mega_variable_param",sep=""))
          save(Mega_variable_corr, file = paste("~/Dropbox/ampOddClick/r-filesFran/Mega_variable_corr",sep=""))
          save(Mega_variable_fits, file = paste("~/Dropbox/ampOddClick/r-filesFran/Mega_variable_fits",sep=""))
          print(paste('FITTING COMPLETED FOR ', Subject, Parcel, Hemisphere, Wave, sep = ' '))
        } else {
          next
        }
      }
    }
  }
}

## For later, run plots ----------------------------
# parm <- c(1.946613e+00,  1.385492e-02,  3.302899e+00,  1.259871e+01,  1.000000e+01)
# resMMN   <- ampFUN( parm, tn, rep(T,dim(tn$xpp)[1]), DVstr='DV',returnStr='Model', show=TRUE, chNr=1)




