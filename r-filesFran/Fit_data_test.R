rm(list=ls())
# Install packages (double check with Tobias if you need to do this every time) ---------------------
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


# Load Tobias model ---------------------------------------------------------------------------------
source('~/Dropbox/r-compile/rUtils/R/rUtils.R')
source('~/Dropbox/r-compile/ePhysLab/R/ePhysLab.R')
source('~/Dropbox/ampOddClick/r-files/functions/functions_ampOdd_presynapticPlasticity.R')


# Define which data to fit and parameters -----------------------------------------------------
Wave <- 'P2' # P0, Na, Pa, Nb P50, N1, or P2
if (strcmp(Wave, 'N1') | strcmp(Wave, 'P2'))  {
  sig_latency <- 'LLR'
} else {
  sig_latency <- 'MLR'
}
Type_signal <- 'EEG' # EEG or Source
Hemisphere <- 'L' # L or R, only if choosing sources as type of signal
Parcel <- 'A1' # A1, MBelt, LBelt, or PBelt, only if choosing sources as type of signal
do_lm_fit <- 'YES' # 'YES' or 'NO' 
model_parameter_selection <- 'find_best' # 'predefined' or 'find_best' 
Subject <- '2355'


# Load single trial data file accordingly -----------------------------------------------------------
single_data_files_dir <- '/Volumes/EEGCyno/EEG/ampOddClickHuman/Shared_data_single_trials'
# LLR <- read.table( paste('/Dropbox/',project,'/RData/',fileNameL,sep=''), header=T,sep = "," )
if (strcmp(Type_signal, 'EEG')) {
  ST <- read.table( paste(single_data_files_dir,'/Amp_singl_trials_EEG_',sig_latency,'.txt',sep=''), header=T,sep = "," )
} else {
  ST <- read.table( paste(single_data_files_dir,'/Single_trials_',Parcel,'_',Hemisphere,'_',sig_latency,'.txt',sep=''), header=T,sep = "\t" )
}


# Basic quality control -----------------------------------------------------------------------------
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

# Create variables to feed Tobias model for predictions ---------------------------------------------
norm_int_values <- c(0.3,0.6,0.9) # values to replace 1, 2 and 3 int values (between 0 and 1)
ISI_vector <- c()
Int_vector <- c()
# Determine column to retrieve amplitude values from
if (strcmp(Type_signal, 'EEG')) {
  Name_column <- paste(Type_signal,'_',Wave,sep='')
} else {
  Name_column <- paste(Wave,sep='')
}
# Retrieve real amplitude values recorded
Signal_vector <- c()
pos_vector <- 1
for (i in 1:nrow(ST)) {
  if (strcmp(ST$subjID[i], Subject)) {
    ISI_vector[pos_vector] <- ST$ISI[i]
    Int_vector[pos_vector] <- ST$clickIntensity[i]
    Signal_vector[pos_vector] <- eval(parse(text = paste('ST$',Name_column,'[i]',sep='')))
    pos_vector <- pos_vector + 1
  }
} 
# Normalize real amplitude values to a 0 to 1 scale (to compare later)
Signal_vector_norm = (Signal_vector-min(na.exclude(Signal_vector)))/(max(na.exclude(Signal_vector))-min(na.exclude(Signal_vector)))
# Convert intensity values to a 0 to 1 scale
for (j in 1:length(Int_vector)) {
  if (Int_vector[j] == 1) { 
    Int_vector[j] <- norm_int_values[1]
  } else if (Int_vector[j] == 2) {
    Int_vector[j] <- norm_int_values[2]
  } else if  (Int_vector[j] == 3) {
    Int_vector[j] <- norm_int_values[3]
  } else {
    stop(paste('Row',j,'has no acceptable intensity value'))
  }
} 
# Also, create x axis for plot later
xaxis <- c()
for (k in 1:length(Int_vector)) {
  xaxis[k] <- k
} 


# Determine parameters of model or find the best ones -----------------------------------------
# USE THE REAL AMPLITUDE VALUES MODIFIED BY INTENSITY AS THE INPUTS FOR THE 0 TO 1 AMPLITUDE SCALE
# INCLUDE A LOOP HERE FOR THAT
if (strcmp(model_parameter_selection, 'predefined')) {
  rec_tim_cons <- 1/.1 # Recovery time constant. For example, 1/.5 means .5 seconds to recover half of the resources
  cum_rel_prob <- 0.5 # Cumulative release probability. Between 0 and 1. High means a lot of resources invested in responding to the sound
} else if (strcmp(model_parameter_selection, 'find_best')) {
  # Test different parameters iteratively 
  rtc_to_test <- c(1/.1, 1/.2, 1/.3, 1/.4, 1/.5, 1/.6, 1/.7, 1/.8, 1/.9) # rec_tim_con values to test
  rtc_strings <- c('1/.1', '1/.2', '1/.3', '1/.4', '1/.5', '1/.6', '1/.7', '1/.8', '1/.9')
  crp_to_test <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9) # cum_rel_prob
  # Create table to store values
  columns_table_corr = c("rtc","crp","corr_value") 
  Table_corr = data.frame(matrix(nrow = 0, ncol = length(columns_table_corr))) 
  colnames(Table_corr) = columns_table_corr
  pos_row <- 1
  for (i in 1:length(rtc_to_test)) {
    for (j in 1:length(crp_to_test)) {
      current_rtc <- rtc_to_test[i]
      current_crp <- crp_to_test[j]
      predicted_values <- cumsumamp( ISI_vector, Int_vector, current_rtc, current_crp )
      if (strcmp(do_lm_fit, 'YES')) {
        lm1 <- lm( Signal_vector_norm ~ predicted_values ) # Should I remove Nans from signal (if so, lm fails for mismatched length)?
        predicted_values <- lm1$fitted
      }
      corr_values <- ccf(predicted_values, na.exclude(Signal_vector_norm)) # Compare time series
      curr_corr_value <- corr_values$acf[which(corr_values$lag == 0)] # Get the one without lag
      Table_corr[pos_row,1] <- rtc_strings[i]
      Table_corr[pos_row,2] <- current_crp
      Table_corr[pos_row,3] <- curr_corr_value
      pos_row <- pos_row + 1
    } 
  } 
  print(Table_corr)
  # Now retrieve best  rtc and crp values to use
  idx_high_corr <- which.max(Table_corr[,3])
  rec_tim_cons <- Table_corr[idx_high_corr,1] # It is a string
  numerador <- str2num(substr(rec_tim_cons, 1, 1))
  denominador <- str2num(substr(rec_tim_cons, 3, 4))
  rec_tim_cons <- numerador/denominador
  cum_rel_prob <- Table_corr[idx_high_corr,2]
}
# PLOT WITH A HEAT MAT THE RTC AND CRP

# Compute prediction with established parameters  -------------------------------------------------------
predicted_values <- cumsumamp( ISI_vector, Int_vector, rec_tim_cons, cum_rel_prob ) # Tobias model
if (strcmp(do_lm_fit, 'YES')) {
  lm1 <- lm( Signal_vector_norm ~ predicted_values ) # Should I remove Nans from signal (if so, lm fails for mismatched length)?
  predicted_values <- lm1$fitted
}
corr_values <- ccf(predicted_values, na.exclude(Signal_vector_norm)) # Compare time series
corr_value <- corr_values$acf[which(corr_values$lag == 0)] # Get the one without lag
print(paste('Correlation value is ',corr_value)) 
# CORR ALONE




# Separate amplitude and predicted values in Int and ISI groups and average for plotting ------
condition_string <- c('11','12','13','31','32','33','51','52','53','71','72','73','91','92','93')
for (i in 1:length(condition_string)) {
    eval(parse(text = paste('PV',condition_string[i],' <- c()',sep='')))
    eval(parse(text = paste('RV',condition_string[i],' <- c()',sep='')))
}
for (l in 1:length(predicted_values)) {
  if (Int_vector[l] == norm_int_values[1] & ISI_vector[l] >= 0.24 & ISI_vector[l] <= 0.5) { 
    PV11[length(PV11)+1] <- predicted_values[l]
    RV11[length(RV11)+1] <- Signal_vector_norm[l]
  } else if (Int_vector[l] == norm_int_values[2] & ISI_vector[l] >= 0.24 & ISI_vector[l] <= 0.5) { 
    PV12[length(PV12)+1] <- predicted_values[l]
    RV12[length(RV12)+1] <- Signal_vector_norm[l]
  } else if (Int_vector[l] == norm_int_values[3] & ISI_vector[l] >= 0.24 & ISI_vector[l] <= 0.5) { 
    PV13[length(PV13)+1] <- predicted_values[l]
    RV13[length(RV13)+1] <- Signal_vector_norm[l]
  } else if (Int_vector[l] == norm_int_values[1] & ISI_vector[l] > 0.5 & ISI_vector[l] <= 1) { 
    PV31[length(PV31)+1] <- predicted_values[l]
    RV31[length(RV31)+1] <- Signal_vector_norm[l]
  } else if (Int_vector[l] == norm_int_values[2] & ISI_vector[l] > 0.5 & ISI_vector[l] <= 1) { 
    PV32[length(PV32)+1] <- predicted_values[l]
    RV32[length(RV32)+1] <- Signal_vector_norm[l]
  } else if (Int_vector[l] == norm_int_values[3] & ISI_vector[l] > 0.5 & ISI_vector[l] <= 1) { 
    PV33[length(PV33)+1] <- predicted_values[l]
    RV33[length(RV33)+1] <- Signal_vector_norm[l]
  } else if (Int_vector[l] == norm_int_values[1] & ISI_vector[l] > 1 & ISI_vector[l] <= 2) { 
    PV51[length(PV51)+1] <- predicted_values[l]
    RV51[length(RV51)+1] <- Signal_vector_norm[l]
  } else if (Int_vector[l] == norm_int_values[2] & ISI_vector[l] > 1 & ISI_vector[l] <= 2) { 
    PV52[length(PV52)+1] <- predicted_values[l]
    RV52[length(RV52)+1] <- Signal_vector_norm[l]
  } else if (Int_vector[l] == norm_int_values[3] & ISI_vector[l] > 1 & ISI_vector[l] <= 2) { 
    PV53[length(PV53)+1] <- predicted_values[l]
    RV53[length(RV53)+1] <- Signal_vector_norm[l]
  } else if (Int_vector[l] == norm_int_values[1] & ISI_vector[l] > 2 & ISI_vector[l] <= 4) { 
    PV71[length(PV71)+1] <- predicted_values[l]
    RV71[length(RV71)+1] <- Signal_vector_norm[l]
  } else if (Int_vector[l] == norm_int_values[2] & ISI_vector[l] > 2 & ISI_vector[l] <= 4) { 
    PV72[length(PV72)+1] <- predicted_values[l]
    RV72[length(RV72)+1] <- Signal_vector_norm[l]
  } else if (Int_vector[l] == norm_int_values[3] & ISI_vector[l] > 2 & ISI_vector[l] <= 4) { 
    PV73[length(PV73)+1] <- predicted_values[l]
    RV73[length(RV73)+1] <- Signal_vector_norm[l]
  } else if (Int_vector[l] == norm_int_values[1] & ISI_vector[l] > 4) { 
    PV91[length(PV91)+1] <- predicted_values[l]
    RV91[length(RV91)+1] <- Signal_vector_norm[l]
  } else if (Int_vector[l] == norm_int_values[2] & ISI_vector[l] > 4) { 
    PV92[length(PV92)+1] <- predicted_values[l]
    RV92[length(RV92)+1] <- Signal_vector_norm[l]
  } else if (Int_vector[l] == norm_int_values[3] & ISI_vector[l] > 4) { 
    PV93[length(PV93)+1] <- predicted_values[l]
    RV93[length(RV93)+1] <- Signal_vector_norm[l]
  } else {
    stop(paste('Row',j,'has no acceptable intensity value'))
  }
} 
# Store average values 
Plot_table_real = data.frame(matrix(nrow = 5, ncol = 3)) 
colnames(Plot_table_real) = c("65dB","75dB","85dB") 
rownames(Plot_table_real) = c("0.25-0.5s","0.5-1s","1-2s","2-4s","4-8s") 
if (1 == 1) { # Silly way to compress many lines
Plot_table_real[1,1] <- mean(RV11,na.rm = TRUE)
Plot_table_real[1,2] <- mean(RV12,na.rm = TRUE)
Plot_table_real[1,3] <- mean(RV13,na.rm = TRUE)
Plot_table_real[2,1] <- mean(RV31,na.rm = TRUE)
Plot_table_real[2,2] <- mean(RV32,na.rm = TRUE)
Plot_table_real[2,3] <- mean(RV33,na.rm = TRUE)
Plot_table_real[3,1] <- mean(RV51,na.rm = TRUE)
Plot_table_real[3,2] <- mean(RV52,na.rm = TRUE)
Plot_table_real[3,3] <- mean(RV53,na.rm = TRUE)
Plot_table_real[4,1] <- mean(RV71,na.rm = TRUE)
Plot_table_real[4,2] <- mean(RV72,na.rm = TRUE)
Plot_table_real[4,3] <- mean(RV73,na.rm = TRUE)
Plot_table_real[5,1] <- mean(RV91,na.rm = TRUE)
Plot_table_real[5,2] <- mean(RV92,na.rm = TRUE)
Plot_table_real[5,3] <- mean(RV93,na.rm = TRUE)
Plot_table_predicted = data.frame(matrix(nrow = 5, ncol = 3)) 
colnames(Plot_table_predicted) = c("65dB","75dB","85dB") 
rownames(Plot_table_predicted) = c("0.25-0.5s","0.5-1s","1-2s","2-4s","4-8s") 
Plot_table_predicted[1,1] <- mean(PV11,na.rm = TRUE)
Plot_table_predicted[1,2] <- mean(PV12,na.rm = TRUE)
Plot_table_predicted[1,3] <- mean(PV13,na.rm = TRUE)
Plot_table_predicted[2,1] <- mean(PV31,na.rm = TRUE)
Plot_table_predicted[2,2] <- mean(PV32,na.rm = TRUE)
Plot_table_predicted[2,3] <- mean(PV33,na.rm = TRUE)
Plot_table_predicted[3,1] <- mean(PV51,na.rm = TRUE)
Plot_table_predicted[3,2] <- mean(PV52,na.rm = TRUE)
Plot_table_predicted[3,3] <- mean(PV53,na.rm = TRUE)
Plot_table_predicted[4,1] <- mean(PV71,na.rm = TRUE)
Plot_table_predicted[4,2] <- mean(PV72,na.rm = TRUE)
Plot_table_predicted[4,3] <- mean(PV73,na.rm = TRUE)
Plot_table_predicted[5,1] <- mean(PV91,na.rm = TRUE)
Plot_table_predicted[5,2] <- mean(PV92,na.rm = TRUE)
Plot_table_predicted[5,3] <- mean(PV93,na.rm = TRUE)
}
# Store standard deviation values 
Stdev_table_real = data.frame(matrix(nrow = 5, ncol = 3)) 
colnames(Stdev_table_real) = c("65dB","75dB","85dB") 
rownames(Stdev_table_real) = c("0.25-0.5s","0.5-1s","1-2s","2-4s","4-8s") 
if (1 == 1) { # Silly way to compress many lines
  Stdev_table_real[1,1] <- sd(RV11,na.rm = TRUE)
  Stdev_table_real[1,2] <- sd(RV12,na.rm = TRUE)
  Stdev_table_real[1,3] <- sd(RV13,na.rm = TRUE)
  Stdev_table_real[2,1] <- sd(RV31,na.rm = TRUE)
  Stdev_table_real[2,2] <- sd(RV32,na.rm = TRUE)
  Stdev_table_real[2,3] <- sd(RV33,na.rm = TRUE)
  Stdev_table_real[3,1] <- sd(RV51,na.rm = TRUE)
  Stdev_table_real[3,2] <- sd(RV52,na.rm = TRUE)
  Stdev_table_real[3,3] <- sd(RV53,na.rm = TRUE)
  Stdev_table_real[4,1] <- sd(RV71,na.rm = TRUE)
  Stdev_table_real[4,2] <- sd(RV72,na.rm = TRUE)
  Stdev_table_real[4,3] <- sd(RV73,na.rm = TRUE)
  Stdev_table_real[5,1] <- sd(RV91,na.rm = TRUE)
  Stdev_table_real[5,2] <- sd(RV92,na.rm = TRUE)
  Stdev_table_real[5,3] <- sd(RV93,na.rm = TRUE)
  Stdev_table_predicted = data.frame(matrix(nrow = 5, ncol = 3)) 
  colnames(Stdev_table_predicted) = c("65dB","75dB","85dB") 
  rownames(Stdev_table_predicted) = c("0.25-0.5s","0.5-1s","1-2s","2-4s","4-8s") 
  Stdev_table_predicted[1,1] <- sd(PV11,na.rm = TRUE)
  Stdev_table_predicted[1,2] <- sd(PV12,na.rm = TRUE)
  Stdev_table_predicted[1,3] <- sd(PV13,na.rm = TRUE)
  Stdev_table_predicted[2,1] <- sd(PV31,na.rm = TRUE)
  Stdev_table_predicted[2,2] <- sd(PV32,na.rm = TRUE)
  Stdev_table_predicted[2,3] <- sd(PV33,na.rm = TRUE)
  Stdev_table_predicted[3,1] <- sd(PV51,na.rm = TRUE)
  Stdev_table_predicted[3,2] <- sd(PV52,na.rm = TRUE)
  Stdev_table_predicted[3,3] <- sd(PV53,na.rm = TRUE)
  Stdev_table_predicted[4,1] <- sd(PV71,na.rm = TRUE)
  Stdev_table_predicted[4,2] <- sd(PV72,na.rm = TRUE)
  Stdev_table_predicted[4,3] <- sd(PV73,na.rm = TRUE)
  Stdev_table_predicted[5,1] <- sd(PV91,na.rm = TRUE)
  Stdev_table_predicted[5,2] <- sd(PV92,na.rm = TRUE)
  Stdev_table_predicted[5,3] <- sd(PV93,na.rm = TRUE)
}
# Compute marginal means
if (1 == 1) { # Silly way to compress many lines
RV_65dB <- mean(as.numeric(Plot_table_real[1:5,1]))
RV_75dB <- mean(as.numeric(Plot_table_real[1:5,2]))
RV_85dB <- mean(as.numeric(Plot_table_real[1:5,3]))
RV_shortest <- mean(as.numeric(Plot_table_real[1,1:3]))
RV_short <- mean(as.numeric(Plot_table_real[2,1:3]))
RV_medium <- mean(as.numeric(Plot_table_real[3,1:3]))
RV_long <- mean(as.numeric(Plot_table_real[4,1:3]))
RV_longest <- mean(as.numeric(Plot_table_real[5,1:3]))
PV_65dB <- mean(as.numeric(Plot_table_predicted[1:5,1]))
PV_75dB <- mean(as.numeric(Plot_table_predicted[1:5,2]))
PV_85dB <- mean(as.numeric(Plot_table_predicted[1:5,3]))
PV_shortest <- mean(as.numeric(Plot_table_predicted[1,1:3]))
PV_short <- mean(as.numeric(Plot_table_predicted[2,1:3]))
PV_medium <- mean(as.numeric(Plot_table_predicted[3,1:3]))
PV_long <- mean(as.numeric(Plot_table_predicted[4,1:3]))
PV_longest <- mean(as.numeric(Plot_table_predicted[5,1:3]))
}
# Compute marginal std dev
if (1 == 1) { # Silly way to compress many lines
STDR_65dB <- sqrt(mean(as.numeric(Stdev_table_real[1:5,1])))
STDR_75dB <- sqrt(mean(as.numeric(Stdev_table_real[1:5,2])))
STDR_85dB <- sqrt(mean(as.numeric(Stdev_table_real[1:5,3])))
STDR_shortest <- sqrt(mean(as.numeric(Stdev_table_real[1,1:3])))
STDR_short <- sqrt(mean(as.numeric(Stdev_table_real[2,1:3])))
STDR_medium <- sqrt(mean(as.numeric(Stdev_table_real[3,1:3])))
STDR_long <- sqrt(mean(as.numeric(Stdev_table_real[4,1:3])))
STDR_longest <- sqrt(mean(as.numeric(Stdev_table_real[5,1:3])))
STDP_65dB <- sqrt(mean(as.numeric(Stdev_table_predicted[1:5,1])))
STDP_75dB <- sqrt(mean(as.numeric(Stdev_table_predicted[1:5,2])))
STDP_85dB <- sqrt(mean(as.numeric(Stdev_table_predicted[1:5,3])))
STDP_shortest <- sqrt(mean(as.numeric(Stdev_table_predicted[1,1:3])))
STDP_short <- sqrt(mean(as.numeric(Stdev_table_predicted[2,1:3])))
STDP_medium <- sqrt(mean(as.numeric(Stdev_table_predicted[3,1:3])))
STDP_long <- sqrt(mean(as.numeric(Stdev_table_predicted[4,1:3])))
STDP_longest <- sqrt(mean(as.numeric(Stdev_table_predicted[5,1:3])))
}


# Plots predicted vs real  -------------------------------------------------------------------------
if (1 == 1) { # Silly way to compress section
#  Plot Intensity
dfR<-data.frame(Mean=c(RV_65dB,RV_75dB,RV_85dB),
               sd=c(STDR_65dB,STDR_75dB,STDR_85dB),
               Intensity=c('65dB','75dB','85dB'))
dfP<-data.frame(Mean=c(PV_65dB,PV_75dB,PV_85dB),
                sd=c(STDP_65dB,STDP_75dB,STDP_85dB),
                Intensity=c('65dB','75dB','85dB'))
P <- ggplot()+
  geom_point(data = dfR, aes(x = Intensity, y = Mean, color = "green"), shape = 21)+
  geom_point(data = dfP, aes(x = Intensity, y = Mean, color = "red"), shape = 21, show.legend = TRUE)+
  scale_color_manual(values = c("green","red"), labels = c("Real","Predicted"), name = "")
P
# Plot ISI
dfR<-data.frame(Mean=c(RV_shortest,RV_short,RV_medium,RV_long,RV_longest),
                sd=c(STDR_shortest,STDR_short,STDR_medium,STDR_long,STDR_longest),
                ISI_range = c("0.25-0.5s","0.5-1s","1-2s","2-4s","4-8s"))
dfP<-data.frame(Mean=c(PV_shortest,PV_short,PV_medium,PV_long,PV_longest),
                sd=c(STDP_shortest,STDP_short,STDP_medium,STDP_long,STDP_longest),
                ISI_range = c("0.25-0.5s","0.5-1s","1-2s","2-4s","4-8s"))
P <- ggplot()+
  geom_point(data = dfR, aes(x = ISI_range, y = Mean, color = "green"), shape = 21)+
  geom_point(data = dfP, aes(x = ISI_range, y = Mean, color = "red"), shape = 21, show.legend = TRUE)+
  scale_color_manual(values = c("green","red"), labels = c("Real","Predicted"), name = "")
P
}

# Play manually with values and plot -----------------------------------------------------
rec_tim_cons <- 1/.1 # Recovery time constant. For example, 1/.5 means .5 seconds to recover half of the resources
cum_rel_prob <- 0.5 # Cumulative release probability. Between 0 and 1. High means a lot of resources invested in responding to the sound
predicted_values_test <- cumsumamp( ISI_vector, Int_vector, rec_tim_cons, cum_rel_prob )
num_in_plot <- 1:1000 # Number of trials in plot and which ones (e.g. 800:2000)
ISI_vector_r <- ISI_vector[num_in_plot]
Int_vector_r <- Int_vector[num_in_plot]
Signal_vector_norm_r <- Signal_vector_norm[num_in_plot]
xaxis_r <- num_in_plot
predicted_values_r <- predicted_values[num_in_plot]
# Plot with dots
plot(xaxis_r,predicted_values_r,ylim=c(0,1),col="red")
points(xaxis_r,Signal_vector_norm_r,ylim=c(0,1),col="green")
legend(0, 1, legend=c("Predicted", "Real"),
       fill = c("red","green"))
# Plot with lines
plot(xaxis_r,predicted_values_r,type="line",ylim=c(0,1),col="red")
lines(xaxis_r,Signal_vector_norm_r,ylim=c(0,1),col="green")
legend(0, 1, legend=c("Predicted", "Real"),
       fill = c("red","green"))

# Plots with all (too crowded)
# plot(xaxis,predicted_values,type="line",ylim=c(0,1),col="red")
# lines(xaxis,Signal_vector_norm,ylim=c(0,1),col="green")

# Classic correlation plot?
plot(predicted_values, Signal_vector_norm)

# corr_values <- ccf(predicted_values_r, na.exclude(Signal_vector_norm_r)) # Compare time series
# corr_value <- corr_values$acf[which(corr_values$lag == 0)] # Get the one without lag
# print(paste('Correlation value is ',corr_value)) 




# do not normalize real data between 0 and 1
# std error of the mean plot

## Old way of computing marginal means and plotting... ------
# condition_string <- c('11','12','13','31','32','33','51','52','53','71','72','73','91','92','93')
# for (i in 1:length(condition_string)) {
#   eval(parse(text = paste('PV',condition_string[i],' <- c()',sep='')))
#   eval(parse(text = paste('RV',condition_string[i],' <- c()',sep='')))
# }
# for (l in 1:length(predicted_values)) {
#   if (Int_vector[l] == norm_int_values[1] & ISI_vector[l] >= 0.24 & ISI_vector[l] <= 0.5) { 
#     PV11[length(PV11)+1] <- predicted_values[l]
#     RV11[length(RV11)+1] <- Signal_vector[l]
#   } else if (Int_vector[l] == norm_int_values[2] & ISI_vector[l] >= 0.24 & ISI_vector[l] <= 0.5) { 
#     PV12[length(PV12)+1] <- predicted_values[l]
#     RV12[length(RV12)+1] <- Signal_vector[l]
#   } else if (Int_vector[l] == norm_int_values[3] & ISI_vector[l] >= 0.24 & ISI_vector[l] <= 0.5) { 
#     PV13[length(PV13)+1] <- predicted_values[l]
#     RV13[length(RV13)+1] <- Signal_vector[l]
#   } else if (Int_vector[l] == norm_int_values[1] & ISI_vector[l] > 0.5 & ISI_vector[l] <= 1) { 
#     PV31[length(PV31)+1] <- predicted_values[l]
#     RV31[length(RV31)+1] <- Signal_vector[l]
#   } else if (Int_vector[l] == norm_int_values[2] & ISI_vector[l] > 0.5 & ISI_vector[l] <= 1) { 
#     PV32[length(PV32)+1] <- predicted_values[l]
#     RV32[length(RV32)+1] <- Signal_vector[l]
#   } else if (Int_vector[l] == norm_int_values[3] & ISI_vector[l] > 0.5 & ISI_vector[l] <= 1) { 
#     PV33[length(PV33)+1] <- predicted_values[l]
#     RV33[length(RV33)+1] <- Signal_vector[l]
#   } else if (Int_vector[l] == norm_int_values[1] & ISI_vector[l] > 1 & ISI_vector[l] <= 2) { 
#     PV51[length(PV51)+1] <- predicted_values[l]
#     RV51[length(RV51)+1] <- Signal_vector[l]
#   } else if (Int_vector[l] == norm_int_values[2] & ISI_vector[l] > 1 & ISI_vector[l] <= 2) { 
#     PV52[length(PV52)+1] <- predicted_values[l]
#     RV52[length(RV52)+1] <- Signal_vector[l]
#   } else if (Int_vector[l] == norm_int_values[3] & ISI_vector[l] > 1 & ISI_vector[l] <= 2) { 
#     PV53[length(PV53)+1] <- predicted_values[l]
#     RV53[length(RV53)+1] <- Signal_vector[l]
#   } else if (Int_vector[l] == norm_int_values[1] & ISI_vector[l] > 2 & ISI_vector[l] <= 4) { 
#     PV71[length(PV71)+1] <- predicted_values[l]
#     RV71[length(RV71)+1] <- Signal_vector[l]
#   } else if (Int_vector[l] == norm_int_values[2] & ISI_vector[l] > 2 & ISI_vector[l] <= 4) { 
#     PV72[length(PV72)+1] <- predicted_values[l]
#     RV72[length(RV72)+1] <- Signal_vector[l]
#   } else if (Int_vector[l] == norm_int_values[3] & ISI_vector[l] > 2 & ISI_vector[l] <= 4) { 
#     PV73[length(PV73)+1] <- predicted_values[l]
#     RV73[length(RV73)+1] <- Signal_vector[l]
#   } else if (Int_vector[l] == norm_int_values[1] & ISI_vector[l] > 4) { 
#     PV91[length(PV91)+1] <- predicted_values[l]
#     RV91[length(RV91)+1] <- Signal_vector[l]
#   } else if (Int_vector[l] == norm_int_values[2] & ISI_vector[l] > 4) { 
#     PV92[length(PV92)+1] <- predicted_values[l]
#     RV92[length(RV92)+1] <- Signal_vector[l]
#   } else if (Int_vector[l] == norm_int_values[3] & ISI_vector[l] > 4) { 
#     PV93[length(PV93)+1] <- predicted_values[l]
#     RV93[length(RV93)+1] <- Signal_vector[l]
#   } else {
#     stop(paste('Row',j,'has no acceptable intensity value'))
#   }
# } 
# # Store average values 
# Plot_table_real = data.frame(matrix(nrow = 5, ncol = 3)) 
# colnames(Plot_table_real) = c("65dB","75dB","85dB") 
# rownames(Plot_table_real) = c("0.25-0.5s","0.5-1s","1-2s","2-4s","4-8s") 
# if (1 == 1) { # Silly way to compress many lines
#   Plot_table_real[1,1] <- mean(RV11,na.rm = TRUE)
#   Plot_table_real[1,2] <- mean(RV12,na.rm = TRUE)
#   Plot_table_real[1,3] <- mean(RV13,na.rm = TRUE)
#   Plot_table_real[2,1] <- mean(RV31,na.rm = TRUE)
#   Plot_table_real[2,2] <- mean(RV32,na.rm = TRUE)
#   Plot_table_real[2,3] <- mean(RV33,na.rm = TRUE)
#   Plot_table_real[3,1] <- mean(RV51,na.rm = TRUE)
#   Plot_table_real[3,2] <- mean(RV52,na.rm = TRUE)
#   Plot_table_real[3,3] <- mean(RV53,na.rm = TRUE)
#   Plot_table_real[4,1] <- mean(RV71,na.rm = TRUE)
#   Plot_table_real[4,2] <- mean(RV72,na.rm = TRUE)
#   Plot_table_real[4,3] <- mean(RV73,na.rm = TRUE)
#   Plot_table_real[5,1] <- mean(RV91,na.rm = TRUE)
#   Plot_table_real[5,2] <- mean(RV92,na.rm = TRUE)
#   Plot_table_real[5,3] <- mean(RV93,na.rm = TRUE)
#   Plot_table_predicted = data.frame(matrix(nrow = 5, ncol = 3)) 
#   colnames(Plot_table_predicted) = c("65dB","75dB","85dB") 
#   rownames(Plot_table_predicted) = c("0.25-0.5s","0.5-1s","1-2s","2-4s","4-8s") 
#   Plot_table_predicted[1,1] <- mean(PV11,na.rm = TRUE)
#   Plot_table_predicted[1,2] <- mean(PV12,na.rm = TRUE)
#   Plot_table_predicted[1,3] <- mean(PV13,na.rm = TRUE)
#   Plot_table_predicted[2,1] <- mean(PV31,na.rm = TRUE)
#   Plot_table_predicted[2,2] <- mean(PV32,na.rm = TRUE)
#   Plot_table_predicted[2,3] <- mean(PV33,na.rm = TRUE)
#   Plot_table_predicted[3,1] <- mean(PV51,na.rm = TRUE)
#   Plot_table_predicted[3,2] <- mean(PV52,na.rm = TRUE)
#   Plot_table_predicted[3,3] <- mean(PV53,na.rm = TRUE)
#   Plot_table_predicted[4,1] <- mean(PV71,na.rm = TRUE)
#   Plot_table_predicted[4,2] <- mean(PV72,na.rm = TRUE)
#   Plot_table_predicted[4,3] <- mean(PV73,na.rm = TRUE)
#   Plot_table_predicted[5,1] <- mean(PV91,na.rm = TRUE)
#   Plot_table_predicted[5,2] <- mean(PV92,na.rm = TRUE)
#   Plot_table_predicted[5,3] <- mean(PV93,na.rm = TRUE)
# }
# # Store standard deviation values 
# Stdev_table_real = data.frame(matrix(nrow = 5, ncol = 3)) 
# colnames(Stdev_table_real) = c("65dB","75dB","85dB") 
# rownames(Stdev_table_real) = c("0.25-0.5s","0.5-1s","1-2s","2-4s","4-8s") 
# if (1 == 1) { # Silly way to compress many lines
#   Stdev_table_real[1,1] <- sd(RV11,na.rm = TRUE)
#   Stdev_table_real[1,2] <- sd(RV12,na.rm = TRUE)
#   Stdev_table_real[1,3] <- sd(RV13,na.rm = TRUE)
#   Stdev_table_real[2,1] <- sd(RV31,na.rm = TRUE)
#   Stdev_table_real[2,2] <- sd(RV32,na.rm = TRUE)
#   Stdev_table_real[2,3] <- sd(RV33,na.rm = TRUE)
#   Stdev_table_real[3,1] <- sd(RV51,na.rm = TRUE)
#   Stdev_table_real[3,2] <- sd(RV52,na.rm = TRUE)
#   Stdev_table_real[3,3] <- sd(RV53,na.rm = TRUE)
#   Stdev_table_real[4,1] <- sd(RV71,na.rm = TRUE)
#   Stdev_table_real[4,2] <- sd(RV72,na.rm = TRUE)
#   Stdev_table_real[4,3] <- sd(RV73,na.rm = TRUE)
#   Stdev_table_real[5,1] <- sd(RV91,na.rm = TRUE)
#   Stdev_table_real[5,2] <- sd(RV92,na.rm = TRUE)
#   Stdev_table_real[5,3] <- sd(RV93,na.rm = TRUE)
#   Stdev_table_predicted = data.frame(matrix(nrow = 5, ncol = 3)) 
#   colnames(Stdev_table_predicted) = c("65dB","75dB","85dB") 
#   rownames(Stdev_table_predicted) = c("0.25-0.5s","0.5-1s","1-2s","2-4s","4-8s") 
#   Stdev_table_predicted[1,1] <- sd(PV11,na.rm = TRUE)
#   Stdev_table_predicted[1,2] <- sd(PV12,na.rm = TRUE)
#   Stdev_table_predicted[1,3] <- sd(PV13,na.rm = TRUE)
#   Stdev_table_predicted[2,1] <- sd(PV31,na.rm = TRUE)
#   Stdev_table_predicted[2,2] <- sd(PV32,na.rm = TRUE)
#   Stdev_table_predicted[2,3] <- sd(PV33,na.rm = TRUE)
#   Stdev_table_predicted[3,1] <- sd(PV51,na.rm = TRUE)
#   Stdev_table_predicted[3,2] <- sd(PV52,na.rm = TRUE)
#   Stdev_table_predicted[3,3] <- sd(PV53,na.rm = TRUE)
#   Stdev_table_predicted[4,1] <- sd(PV71,na.rm = TRUE)
#   Stdev_table_predicted[4,2] <- sd(PV72,na.rm = TRUE)
#   Stdev_table_predicted[4,3] <- sd(PV73,na.rm = TRUE)
#   Stdev_table_predicted[5,1] <- sd(PV91,na.rm = TRUE)
#   Stdev_table_predicted[5,2] <- sd(PV92,na.rm = TRUE)
#   Stdev_table_predicted[5,3] <- sd(PV93,na.rm = TRUE)
# }
# # Compute marginal means
# if (1 == 1) { # Silly way to compress many lines
#   RV_65dB <- mean(as.numeric(Plot_table_real[1:5,1]))
#   RV_75dB <- mean(as.numeric(Plot_table_real[1:5,2]))
#   RV_85dB <- mean(as.numeric(Plot_table_real[1:5,3]))
#   RV_shortest <- mean(as.numeric(Plot_table_real[1,1:3]))
#   RV_short <- mean(as.numeric(Plot_table_real[2,1:3]))
#   RV_medium <- mean(as.numeric(Plot_table_real[3,1:3]))
#   RV_long <- mean(as.numeric(Plot_table_real[4,1:3]))
#   RV_longest <- mean(as.numeric(Plot_table_real[5,1:3]))
#   PV_65dB <- mean(as.numeric(Plot_table_predicted[1:5,1]))
#   PV_75dB <- mean(as.numeric(Plot_table_predicted[1:5,2]))
#   PV_85dB <- mean(as.numeric(Plot_table_predicted[1:5,3]))
#   PV_shortest <- mean(as.numeric(Plot_table_predicted[1,1:3]))
#   PV_short <- mean(as.numeric(Plot_table_predicted[2,1:3]))
#   PV_medium <- mean(as.numeric(Plot_table_predicted[3,1:3]))
#   PV_long <- mean(as.numeric(Plot_table_predicted[4,1:3]))
#   PV_longest <- mean(as.numeric(Plot_table_predicted[5,1:3]))
# }
# # Compute marginal std dev
# if (1 == 1) { # Silly way to compress many lines
#   STDR_65dB <- sqrt(mean(as.numeric(Stdev_table_real[1:5,1])))
#   STDR_75dB <- sqrt(mean(as.numeric(Stdev_table_real[1:5,2])))
#   STDR_85dB <- sqrt(mean(as.numeric(Stdev_table_real[1:5,3])))
#   STDR_shortest <- sqrt(mean(as.numeric(Stdev_table_real[1,1:3])))
#   STDR_short <- sqrt(mean(as.numeric(Stdev_table_real[2,1:3])))
#   STDR_medium <- sqrt(mean(as.numeric(Stdev_table_real[3,1:3])))
#   STDR_long <- sqrt(mean(as.numeric(Stdev_table_real[4,1:3])))
#   STDR_longest <- sqrt(mean(as.numeric(Stdev_table_real[5,1:3])))
#   STDP_65dB <- sqrt(mean(as.numeric(Stdev_table_predicted[1:5,1])))
#   STDP_75dB <- sqrt(mean(as.numeric(Stdev_table_predicted[1:5,2])))
#   STDP_85dB <- sqrt(mean(as.numeric(Stdev_table_predicted[1:5,3])))
#   STDP_shortest <- sqrt(mean(as.numeric(Stdev_table_predicted[1,1:3])))
#   STDP_short <- sqrt(mean(as.numeric(Stdev_table_predicted[2,1:3])))
#   STDP_medium <- sqrt(mean(as.numeric(Stdev_table_predicted[3,1:3])))
#   STDP_long <- sqrt(mean(as.numeric(Stdev_table_predicted[4,1:3])))
#   STDP_longest <- sqrt(mean(as.numeric(Stdev_table_predicted[5,1:3])))
# }
# 
# 
# ## Plots predicted vs real  -------------------------------------------------------------------------
# if (1 == 1) { # Silly way to compress section
#   #  Plot Intensity
#   dfR<-data.frame(Mean=c(RV_65dB,RV_75dB,RV_85dB),
#                   sd=c(STDR_65dB,STDR_75dB,STDR_85dB),
#                   Intensity=c('65dB','75dB','85dB'))
#   dfP<-data.frame(Mean=c(PV_65dB,PV_75dB,PV_85dB),
#                   sd=c(STDP_65dB,STDP_75dB,STDP_85dB),
#                   Intensity=c('65dB','75dB','85dB'))
#   P <- ggplot()+
#     geom_point(data = dfR, aes(x = Intensity, y = Mean, color = "green"), shape = 21)+
#     geom_point(data = dfP, aes(x = Intensity, y = Mean, color = "red"), shape = 21, show.legend = TRUE)+
#     scale_color_manual(values = c("green","red"), labels = c("Real","Predicted"), name = "")
#   P
#   # Plot ISI
#   dfR<-data.frame(Mean=c(RV_shortest,RV_short,RV_medium,RV_long,RV_longest),
#                   sd=c(STDR_shortest,STDR_short,STDR_medium,STDR_long,STDR_longest),
#                   ISI_range = c("0.25-0.5s","0.5-1s","1-2s","2-4s","4-8s"))
#   dfP<-data.frame(Mean=c(PV_shortest,PV_short,PV_medium,PV_long,PV_longest),
#                   sd=c(STDP_shortest,STDP_short,STDP_medium,STDP_long,STDP_longest),
#                   ISI_range = c("0.25-0.5s","0.5-1s","1-2s","2-4s","4-8s"))
#   P <- ggplot()+
#     geom_point(data = dfR, aes(x = ISI_range, y = Mean, color = "green"), shape = 21)+
#     geom_point(data = dfP, aes(x = ISI_range, y = Mean, color = "red"), shape = 21, show.legend = TRUE)+
#     scale_color_manual(values = c("green","red"), labels = c("Real","Predicted"), name = "")
#   P
# }
# 
# ## Play manually with values and plot -----------------------------------------------------
# rec_tim_cons <- 1/.1 # Recovery time constant. For example, 1/.5 means .5 seconds to recover half of the resources
# cum_rel_prob <- 0.5 # Cumulative release probability. Between 0 and 1. High means a lot of resources invested in responding to the sound
# predicted_values_test <- cumsumamp( ISI_vector, Int_vector, rec_tim_cons, cum_rel_prob )
# num_in_plot <- 1:1000 # Number of trials in plot and which ones (e.g. 800:2000)
# ISI_vector_r <- ISI_vector[num_in_plot]
# Int_vector_r <- Int_vector[num_in_plot]
# Signal_vector_r <- Signal_vector[num_in_plot]
# xaxis_r <- num_in_plot
# predicted_values_r <- predicted_values[num_in_plot]
# # Plot with dots
# plot(xaxis_r,predicted_values_r,ylim=c(0,1),col="red")
# points(xaxis_r,Signal_vector_r,ylim=c(0,1),col="green")
# legend(0, 1, legend=c("Predicted", "Real"),
#        fill = c("red","green"))
# # Plot with lines
# plot(xaxis_r,predicted_values_r,type="line",ylim=c(0,1),col="red")
# lines(xaxis_r,Signal_vector_r,ylim=c(0,1),col="green")
# legend(0, 1, legend=c("Predicted", "Real"),
#        fill = c("red","green"))
# 
# # Plots with all (too crowded)
# # plot(xaxis,predicted_values,type="line",ylim=c(0,1),col="red")
# # lines(xaxis,Signal_vector,ylim=c(0,1),col="green")
# 
# # Classic correlation plot?
# plot(predicted_values, Signal_vector)
# 
# # corr_values <- ccf(predicted_values_r, na.exclude(Signal_vector_r)) # Compare time series
# # corr_value <- corr_values$acf[which(corr_values$lag == 0)] # Get the one without lag
# # print(paste('Correlation value is ',corr_value)) 
# 
# 
# # Also, create x axis for plot later
# xaxis <- c()
# for (k in 1:length(Int_vector)) {
#   xaxis[k] <- k
# } 
# 
# # do not normalize real data between 0 and 1
# # std error of the mean plot