
function_parloop_gradientdescent_TT_ampFUN2q <- function(subid){
  
  #  rm(list=ls())
  
  Subject_list <- c('2193','2199','2235','2238','2250','2320','2259','2262','2274',
                    '2275','H001','2214','2212','2285','2337','2372','2389','H002',
                    '2299','2406','2426','2494','2427','2456','H005',
                    '2005','2355','2387','2288') 
  # 2494 failed EEG
  
  Subject <- Subject_list[subid]
  
  ## Row column where values will be stored in Mega_variable with all subjects
  #subid <- c(1:length(Subject_list)) 
  
  
  ## Add packages -----------------------
  library(pracma)
  library(lattice)
  library(ggplot2)
  #library(corrplot)
  rda <<- function(){ return('/Volumes/EEGCyno/RES/ampOddClick/') }
  source('~/Dropbox/r-compile/rUtils/R/rUtils.R')
  source('~/Dropbox/r-compile/ePhysLab/R/ePhysLab.R')
  source('~/Dropbox/ampOddClick/r-files/functions/functions_ampOdd_presynapticPlasticity.R')
  source('~/Dropbox/ampOddClick/r-filesFran/load_ampOddClick_Human.R')
  
  ## Define variables -------------------------
  Wave_list <- c('P0','Na','Pa','Nb','P50','N1','P2')
  Type_signal_list <- c('EEG','Source')
  Hemisphere_list <- c('L','R') # L or R, only if choosing sources as type of signal
  Parcel_list <- c('A1','MBelt','LBelt','PBelt') # Only if choosing sources as type of signal
  
  ## Loop through subjects/Waves/Sources
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
  #load(paste("~/Dropbox/ampOddClick/r-filesFran/Mega_variable_param",sep=""))
  #load(paste("~/Dropbox/ampOddClick/r-filesFran/Mega_variable_corr",sep=""))
  
  # # Do EEG alone first
  # Type_signal <- 'EEG'
  # for (w in 1:length(Wave_list)) {
  #   Wave <- Wave_list[w]
  #   Hemisphere <- 'NA' # Irrelevant here
  #   Parcel <- 'NA' # Irrelevant here
  #   
  #   # Run gradient descent if tst file does not exist
  #   #if (!file.exists(paste("/Volumes/EEGCyno/RES/ampOddClick/STPSP/STPSPFit",Subject, "22223344_5566_zm9999_all_DV_ampFUN_p",Subject,Parcel,Hemisphere,Wave,"gd.rda",sep="_")))
  #   print(paste('RUNNING GRADIENT DESCENT FOR ', Subject, Wave, Type_signal, sep = ' '))
  #   
  #   tn <- load_ampOddClick_Human(Wave, Subject, Type_signal, Hemisphere,Parcel)
  #   tn$xpp$dummy = 1
  #   tst <- callSTPSP(tn, DVstr='DV', gd=T, show=TRUE, showCh=1,FUN=ampFUN2q,maxISI=10, minISI=0.2, maxp2p=5000, DVsdscale=5, 
  #                    maxEOG=Inf, trialType=NA, fitType=NA, manStartPar = NULL, robust = FALSE, startFromScratch = T,
  #                    fileExtn = paste("p2q",Subject, Parcel, Hemisphere, Wave, sep = '_'),BYstr='dummy',topInt = 1)
  #   # Store values of best fit in data.frame
  #   #Mega_variable_param[[subid,w]]<- tst[['modelFit']][1]
  #   #Mega_variable_corr[[subid,w]]<- tst[['modelFit']][2]
  #   ## This next line was deleted as Tobias mentioned with modelFit[1] parameters we can recreate the modeled data easily
  #   ## Mega_variable_fits[[subid,w]]<- tst[["res"]][["lm0"]][["model"]][,1]
  #   ## Save Mega_variables
  #   #save(Mega_variable_param, file = paste("~/Dropbox/ampOddClick/r-filesFran/Mega_variable_param",sep=""))
  #   #save(Mega_variable_corr, file = paste("~/Dropbox/ampOddClick/r-filesFran/Mega_variable_corr",sep=""))
  #   ## save(Mega_variable_fits, file = paste("~/Dropbox/ampOddClick/r-filesFran/Mega_variable_fits",sep=""))
  #   print(paste('FITTING COMPLETED FOR ', Subject, Wave, Type_signal, sep = ' '))
  # }
  # Then do Sources
  Type_signal <- 'Source'
  pos <- length(Wave_list) + 1
  for (parc in 1:length(Parcel_list)) {
    Parcel <- Parcel_list[parc]
    for (hem in 1:length(Hemisphere_list)) {
      Hemisphere <- Hemisphere_list[hem]
      for (w in 1:length(Wave_list)) {
        Wave <- Wave_list[w]
        # Run gradient descent if tst file does not exist
        #if (!file.exists(paste(paste("/Volumes/EEGCyno/RES/ampOddClick/STPSP/STPSPFit",Subject, "22223344_5566_zm9999_all_DV_ampFUN_p",Subject,Parcel,Hemisphere,Wave,"gd.rda",sep="_"))))
        print(paste('RUNNING GRADIENT DESCENT FOR ', Subject, Parcel, Hemisphere, Wave, sep = ' '))
        tn <- load_ampOddClick_Human(Wave, Subject, Type_signal, Hemisphere,Parcel)
        tn$xpp$dummy = 1
        tst <- callSTPSP(tn, DVstr='DV', gd=T, show=TRUE, showCh=1,FUN=ampFUN2q,maxISI=10, minISI=0.2, maxp2p=5000, DVsdscale=5,
                         maxEOG=Inf, trialType=NA, fitType=NA, manStartPar = NULL, robust = FALSE, startFromScratch = T,
                         fileExtn = paste("p2q",Subject, Parcel, Hemisphere, Wave, sep = '_'),BYstr='dummy',topInt = 1)
        ## Store values of best fit in data.frame
        #Mega_variable_param[[subid,pos]] <- tst[['modelFit']][1]
        #Mega_variable_corr[[subid,pos]] <- tst[['modelFit']][2]
        ## This next line was deleted as Tobias mentioned with modelFit[1] parameters we can recreate the modeled data easily
        ## Mega_variable_fits[[subid,pos]]<- tst[["res"]][["lm0"]][["model"]][,1]
        pos <- pos + 1
        ## Save Mega_variables
        #save(Mega_variable_param, file = paste("~/Dropbox/ampOddClick/r-filesFran/Mega_variable_param",sep=""))
        #save(Mega_variable_corr, file = paste("~/Dropbox/ampOddClick/r-filesFran/Mega_variable_corr",sep=""))
        # save(Mega_variable_fits, file = paste("~/Dropbox/ampOddClick/r-filesFran/Mega_variable_fits",sep=""))
        print(paste('FITTING COMPLETED FOR ', Subject, Parcel, Hemisphere, Wave, sep = ' '))

      }
    }
  }
}

