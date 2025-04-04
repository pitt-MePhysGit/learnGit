tic
rm(list=ls())
## Install packages ---------------------
install_and_add_packages <- 'NO' # YES or NO
if (strcmp(install_and_add_packages, 'YES'))  {
  install.packages("gdata")
  install.packages("Rcpp")
  install.packages("RcppArmadillo")
  install.packages('pracma')
  library('pracma')
  library(lattice)
  # install.packages("ggpubr") # did not work
  # library('ggpubr')
  install.packages("ggplot2")
  library(ggplot2)
  install.packages("corrplot")
  library(corrplot)
}


## Load model and define which data to fit and parameters -----------------------------------------------------
# Load Tobias model
source('~/Dropbox/r-compile/rUtils/R/rUtils.R')
source('~/Dropbox/r-compile/ePhysLab/R/ePhysLab.R')
source('~/Dropbox/ampOddClick/r-files/functions/functions_ampOdd_presynapticPlasticity.R')
Wave <- 'P0' # P0, Na, Pa, Nb P50, N1, or P2
if (strcmp(Wave, 'N1') | strcmp(Wave, 'P2'))  {
  sig_latency <- 'LLR'
} else {
  sig_latency <- 'MLR'
}
Type_signal <- 'Source' # EEG or Source
Hemisphere <- 'L' # L or R, only if choosing sources as type of signal
Parcel <- 'A1' # A1, MBelt, LBelt, or PBelt, only if choosing sources as type of signal
isi_breaks <- c(.25, .5, 1, 2, 4, 22)
do_lm_fit <- 'YES' # 'YES' or 'NO' 
model_parameter_selection <- 'find_best' # 'predefined' or 'find_best' 
Subject <- '2005' ## 2355 done 


## Load single trial data and basic QC -----------------------------------------------------------
single_data_files_dir <- '/Volumes/EEGCyno/EEG/ampOddClickHuman/Shared_data_single_trials'
# LLR <- read.table( paste('/Dropbox/',project,'/RData/',fileNameL,sep=''), header=T,sep = "," )
if (strcmp(Type_signal, 'EEG')) {
  ST <- read.table( paste(single_data_files_dir,'/Amp_singl_trials_EEG_',sig_latency,'.txt',sep=''), header=T,sep = "," )
} else {
  ST <- read.table( paste(single_data_files_dir,'/Single_trials_',Parcel,'_',Hemisphere,'_',sig_latency,'.txt',sep=''), header=T,sep = "\t" )
}
# Crop it to contain only the subject at hand
ST <- subset.data.frame(ST,ST$subjID == Subject)
# Basic quality control
table(ST$subjID)
table(ST$subjID,ST$sessionNr)
# table(ST$sessionNr, ST$blockNr,ST$subjID)
table(ST$clickIntensity)
# Ensure no nans in intensity values
for (i in 1:nrow(ST)) {
  if (is.nan(ST$clickIntensity[i]) || ST$clickIntensity[i] > 3 || ST$clickIntensity[i] < 1) {
    stop(paste('Row',i,'has no acceptable intensity value'))
  }
} 
# Fix: replace all nan ISIs by a large number of seconds to 'reset' the resources count
for (i in 1:nrow(ST)) {
  if (is.nan(ST$ISI[i])) {
    ST$ISI[i] <- 20 # 20s
  }
}

## Determine parameters of model or find the best ones -----------------------------------------
dir_subj_corr_tables <- "~/Dropbox/ampOddClick/r-filesFran/Subj_corr_tables"
# Determine name of signal column and extract values
if (strcmp(Type_signal, 'EEG')) {
  Name_column <- paste(Type_signal,'_',Wave,sep='')
} else {
  Name_column <- paste(Wave,sep='')
}
Signal_vector <- eval(parse(text = paste('ST$',Name_column,sep='')))
# Discarded: Normalize real amplitude values to a 0 to 1 scale (to compare later)
# Signal_vector = (Signal_vector-min(na.exclude(Signal_vector)))/(max(na.exclude(Signal_vector))-min(na.exclude(Signal_vector)))
# Find best values iteratively or predefine
if (strcmp(model_parameter_selection, 'predefined')) {
  rec_tim_cons <- 1/.1 # Recovery time constant. For example, 1/.5 means .5 seconds to recover half of the resources
  cum_rel_prob <- 0.5 # Cumulative release probability. Between 0 and 1. High means a lot of resources invested in responding to the sound
  norm_int_values <- c(0.3,0.6,0.9) # values to replace 1, 2 and 3 int values (between 0 and 1)
} else if (strcmp(model_parameter_selection, 'find_best')) {
  
  if (!file.exists(paste(dir_subj_corr_tables,"/S",Subject,"_table_corr",sep=""))) { # Do only if Table corr for this subject does not exist
  # Test different parameters iteratively 
  rtc_to_test <- c(1/.1, 1/.2, 1/.3, 1/.4, 1/.5, 1/.6, 1/.7, 1/.8, 1/.9) # rec_tim_con values to test
  rtc_strings <- c('1/.1', '1/.2', '1/.3', '1/.4', '1/.5', '1/.6', '1/.7', '1/.8', '1/.9')
  crp_to_test <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9) # cum_rel_prob
  norm_int_values_to_test <- c(0.3,0.6,0.9) # normalized int values to replace 1, 2 and 3 int values (between 0 and 1)
  # n65dB_to_test <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  # n75dB_to_test <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  # n85dB_to_test <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  n65dB_to_test <- c(0.5)
  n75dB_to_test <- c(0.6)
  n85dB_to_test <- c(0.7)
  # Create table to store values
  columns_table_corr = c("65dB","75dB","85dB","rtc_string","rtc","crp","corr_value") 
  Table_corr = data.frame(matrix(nrow = 0, ncol = length(columns_table_corr))) 
  colnames(Table_corr) = columns_table_corr
  pos_row <- 1
  for (db1 in 1:length(n65dB_to_test)) {
  for (db2 in 1:length(n75dB_to_test)) {
  for (db3 in 1:length(n85dB_to_test)) {
  for (i in 1:length(rtc_to_test)) {
    for (j in 1:length(crp_to_test)) {
      # Normalize int values first
      Int_vector_norm <- c()
      for (ci in 1:length(ST$clickIntensity)) {
        if (ST$clickIntensity[ci] == 1) { 
          Int_vector_norm[ci] <- n65dB_to_test[db1]
        } else if (ST$clickIntensity[ci] == 2) {
          Int_vector_norm[ci] <- n75dB_to_test[db2]
        } else if  (ST$clickIntensity[ci] == 3) {
          Int_vector_norm[ci] <- n85dB_to_test[db3]
        } else {
          stop(paste('Row',ci,'has no acceptable intensity value'))
        }
      } 
      # Compute prediction with these parameters
      print(paste('Computing predictions (rtc =',rtc_strings[i],'crp =',crp_to_test[j],'int =',n65dB_to_test[db1],n75dB_to_test[db2],n85dB_to_test[db3],')')) 
      predicted_values <- cumsumamp( ST$ISI, Int_vector_norm, rtc_to_test[i], crp_to_test[j] )
      # Adjust with linear model
      if (strcmp(do_lm_fit, 'YES')) {
        lm1 <- lm( Signal_vector ~ predicted_values ) # Should I remove Nans from signal (if so, lm fails for mismatched length)?
        predicted_values <- lm1$fitted
      }
      # Compute correlation
      curr_corr_value <- cor(predicted_values, na.exclude(Signal_vector))
      # Store values in Table
      Table_corr[pos_row,1] <- n65dB_to_test[db1]
      Table_corr[pos_row,2] <- n75dB_to_test[db2]
      Table_corr[pos_row,3] <- n85dB_to_test[db3]
      Table_corr[pos_row,4] <- rtc_strings[i]
      Table_corr[pos_row,5] <- rtc_to_test[i]
      Table_corr[pos_row,6] <- crp_to_test[j]
      Table_corr[pos_row,7] <- curr_corr_value
      pos_row <- pos_row + 1
    } 
  } 
  }
  }
  }
  # Save the table
  print(Table_corr)
  save(Table_corr, file = paste(dir_subj_corr_tables,"/S",Subject,"_table_corr",sep=""))
  } else {
    # Load existing table
    load(paste(dir_subj_corr_tables,"/S",Subject,"_table_corr",sep=""))
  }
  # Now retrieve best  rtc and crp values to use
  idx_high_corr <- which.max(Table_corr[,7])
  norm_int_values <- c(Table_corr[idx_high_corr,1],Table_corr[idx_high_corr,2],Table_corr[idx_high_corr,3])
  rec_tim_cons <- Table_corr[idx_high_corr,5]
  cum_rel_prob <- Table_corr[idx_high_corr,6]
}


## Compute prediction with established parameters  -------------------------------------------------------
# Discarded: Normalize real amplitude values to a 0 to 1 scale (to compare later)
# Signal_vector = (Signal_vector-min(na.exclude(Signal_vector)))/(max(na.exclude(Signal_vector))-min(na.exclude(Signal_vector)))
# Convert intensity values to a determined 0 to 1 scale
Int_vector_norm <- c()
for (ci in 1:length(ST$clickIntensity)) {
  if (ST$clickIntensity[ci] == 1) { 
    Int_vector_norm[ci] <- norm_int_values[1]
  } else if (ST$clickIntensity[ci] == 2) {
    Int_vector_norm[ci] <- norm_int_values[2]
  } else if  (ST$clickIntensity[ci] == 3) {
    Int_vector_norm[ci] <- norm_int_values[3]
  } else {
    stop(paste('Row',ci,'has no acceptable intensity value'))
  }
} 
# Store norm values into data.frame
ST$clickIntensity_norm <- Int_vector_norm
# Compute prediction
predicted_values <- cumsumamp( ST$ISI, ST$clickIntensity_norm, rec_tim_cons, cum_rel_prob ) # Tobias model
if (strcmp(do_lm_fit, 'YES')) {
  lm1 <- lm( Signal_vector ~ predicted_values, na.action=na.exclude) # Should I remove Nans from signal (if so, lm fails for mismatched length)?
  predicted_values <- lm1$fitted
}
corr_value <- cor(predicted_values, na.exclude(Signal_vector)) # Compare time series
print(paste('Correlation value is ',corr_value)) 
toc

## Compute marginal means and plot --------------------------------------------------------------

# Heatmap with Table.corr -------------------------------------
# ggplot(Table_corr, aes(x = rtc_string, y = crp, fill = corr_value)) +
#   geom_tile() +
#  scale_fill_gradient(low="black", high="white")
Table_corr_best_int <- subset.data.frame(Table_corr,Table_corr$`65dB` == norm_int_values[1] 
                                         & Table_corr$`75dB` == norm_int_values[2] 
                                         & Table_corr$`85dB` == norm_int_values[3])
ggplot(Table_corr_best_int, aes(x = rtc_string, y = crp, fill = corr_value)) +
  geom_tile() +
  scale_fill_gradient(low="black", high="white")

## Plot average fits -------------------------------------------
# Store predictions as new column in ST data.frame
sized_predicted_values <- c()
pos_sv <- 1
pos_sp <- 1
for (svl in 1:length(Signal_vector)) {
  if (is.nan(Signal_vector[svl]) ) { 
    sized_predicted_values[pos_sp] <- NaN
    pos_sp <- pos_sp + 1
  } else {
    sized_predicted_values[pos_sp] <- predicted_values[pos_sv]
    pos_sp <- pos_sp + 1
    pos_sv <- pos_sv + 1
  }
}
ST$Pred <- sized_predicted_values
ST$ISI_break <- cut(ST$ISI,breaks=isi_breaks)
# Marginal means
# A <- tapply(ST$ISI,ST$ISI_break, mean)
ISI_means_real <- tapply(Signal_vector,ST$ISI_break, na.mean)
ISI_se_real <- tapply(Signal_vector,ST$ISI_break, na.se)
INT_means_real <- tapply(Signal_vector,ST$clickIntensity, na.mean)
INT_se_real <- tapply(Signal_vector,ST$clickIntensity, na.se)
All_means_real <- tapply(Signal_vector,list(ST$ISI_break,ST$clickIntensity), na.mean)
All_se_real <- tapply(Signal_vector,list(ST$ISI_break,ST$clickIntensity), na.se)
ISI_means_pred <- tapply(ST$Pred,ST$ISI_break, na.mean)
ISI_se_pred <- tapply(ST$Pred,ST$ISI_break, na.se)
INT_means_pred <- tapply(ST$Pred,ST$clickIntensity, na.mean)
INT_se_pred <- tapply(ST$Pred,ST$clickIntensity, na.se)
All_means_pred <- tapply(ST$Pred,list(ST$ISI_break,ST$clickIntensity), na.mean)
All_se_pred <- tapply(ST$Pred,list(ST$ISI_break,ST$clickIntensity), na.se)

dfR<-data.frame(Mean=INT_means_real,
                se=INT_se_real,
                Intensity=c('65dB','75dB','85dB'))
dfP<-data.frame(Mean=c(INT_means_pred),
                se=c(INT_se_pred),
                Intensity=c('65dB','75dB','85dB'))
P <- ggplot()+
  geom_point(data = dfR, aes(x = Intensity, y = Mean, color = "green"), shape = 21)+
  geom_errorbar(data = dfR, aes(x = Intensity, ymin=Mean-se, ymax=Mean+se, color = "green"), width=.1) +
  geom_point(data = dfP, aes(x = Intensity, y = Mean, color = "red"), shape = 21, show.legend = TRUE)+
  geom_errorbar(data = dfP, aes(x = Intensity, ymin=Mean-se, ymax=Mean+se, color = "red"), width=.1) +
  scale_color_manual(values = c("green","red"), labels = c("Real","Predicted"), name = "")
P

dfR<-data.frame(Mean=ISI_means_real,
                se=ISI_se_real,
                ISI_range = c("0.25-0.5s","0.5-1s","1-2s","2-4s","4-8s"))
dfP<-data.frame(Mean=ISI_means_pred,
                se=ISI_se_pred,
                ISI_range = c("0.25-0.5s","0.5-1s","1-2s","2-4s","4-8s"))
P <- ggplot()+
  geom_point(data = dfR, aes(x = ISI_range, y = Mean, color = "green"), shape = 21)+
  geom_errorbar(data = dfR, aes(x = ISI_range, ymin=Mean-se, ymax=Mean+se, color = "green"), width=.1) +
  geom_point(data = dfP, aes(x = ISI_range, y = Mean, color = "red"), shape = 21, show.legend = TRUE)+
  geom_errorbar(data = dfP, aes(x = ISI_range, ymin=Mean-se, ymax=Mean+se, color = "red"), width=.1) +
  scale_color_manual(values = c("green","red"), labels = c("Real","Predicted"), name = "")
P


