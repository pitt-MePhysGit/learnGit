## Plot predicted vs real data gradient descend

rm(list=ls())
## Add packages -----------------------
# install.packages('patchwork')
library(patchwork)
library(pracma)
library(lattice)
library(ggplot2)
#library(corrplot)
rda <- function(){ return('/Volumes/EEGCyno/RES/ampOddClick/') }
source('~/Dropbox/r-compile/rUtils/R/rUtils.R')
source('~/Dropbox/r-compile/ePhysLab/R/ePhysLab.R')
source('~/Dropbox/ampOddClick/r-files/functions/functions_ampOdd_presynapticPlasticity.R')
source('~/Dropbox/ampOddClick/r-filesFran/load_ampOddClick_Human.R')

Generate_EEG_plots <- "YES"
Generate_source_plots <- "NO!" # Average not programmed for this yet
Overwrite_plots <- "YES" # If willing to repeat and overwrite present plots

Subject_list <- c('2193','2199','2235','2238','2250','2320','2259','2262','2274',
                  '2275','H001','2214','2212','2285','2337','2372','2389','H002',
                  '2299','2406','2426','2494','2427','2456','H005',
                  '2005','2355','2387','2288') 
## Define variables -------------------------
Wave_list <- c('P0','Na','Pa','Nb','P50','N1','P2')
# Wave_list <- c('P2')
Type_signal_list <- c('EEG','Source')
Hemisphere_list <- c('L','R') # L or R, only if choosing sources as type of signal
Parcel_list <- c('A1','MBelt','LBelt','PBelt') # Only if choosing sources as type of signal

# EEG first
Type_signal <- 'EEG'
for (w in 1:length(Wave_list)) {
  Wave <- Wave_list[w]
  Hemisphere <- 'NA' # Irrelevant here
  Parcel <- 'NA' # Irrelevant here
  
  ISI_means_real <- matrix(list(), length(Subject_list), 5)
  ISI_se_real <- matrix(list(), length(Subject_list), 5)
  INT_means_real <- matrix(list(), length(Subject_list), 3)
  INT_se_real <- matrix(list(), length(Subject_list), 3)
  All_means_real <- matrix(list(), length(Subject_list), 15)
  All_se_real <- matrix(list(), length(Subject_list), 15)
  ISI_means_pred <- matrix(list(), length(Subject_list), 5)
  ISI_se_pred <- matrix(list(), length(Subject_list), 5)
  INT_means_pred <- matrix(list(), length(Subject_list), 3)
  INT_se_pred <- matrix(list(), length(Subject_list), 3)
  All_means_pred <- matrix(list(), length(Subject_list), 15)
  All_se_pred <- matrix(list(), length(Subject_list), 15)


      for (subid in 1:length(Subject_list)) {
        Subject <- Subject_list[subid]
        # Run gradient descent if tst file does not exist
        if (file.exists(paste("/Volumes/EEGCyno/RES/ampOddClick/STPSP/STPSPFit",Subject, "22223344_5566_zm9999_all_DV_ampFUN_po",Subject,Parcel,Hemisphere,Wave,"gd.rda",sep="_"))){
          # Load original values first (to plot against prediction)
          tn <- load_ampOddClick_Human(Wave, Subject, Type_signal, Hemisphere,Parcel)
          # Get the structure
          ST <- tn[['xpp']]
          # Get the column with real original values
          if (strcmp(Type_signal, 'EEG')) {
            Name_column <- paste(Type_signal,'_',Wave,sep='')
          } else {
            Name_column <- paste(Wave,sep='')
          }
          Signal_vector <- eval(parse(text = paste('ST$',Name_column,sep='')))
          # Load the model fit variable (replaced with the calSTPSP)
          # load(paste("/Volumes/EEGCyno/RES/ampOddClick/STPSP/STPSPFit",Subject, "22223344_5566_zm9999_all_DV_ampFUN_po",Subject,Parcel,Hemisphere,Wave,"gd.rda",sep="_"))
          # Obtain modeled time series with fit values (it won't run gradient descent, but retrieve values of files created before)
          tst <- callSTPSP(tn, DVstr='DV', gd=T, show=TRUE, showCh=1,FUN=ampFUN,maxISI=10, minISI=0.2, maxp2p=5000, DVsdscale=5, 
                           maxEOG=Inf, trialType=NA, fitType=NA, manStartPar = NULL, robust = FALSE, 
                           fileExtn = paste("po",Subject, Parcel, Hemisphere, Wave, sep = '_'))
          # Store predicted time series in structure
          ST$Pred <- tst$res$lm2$fitted.values
          # Extract fit values
          # fit_values <- tst[['modelFit']][1]
          
          # Get marginal means real values 
          ISI_means_real[subid,] <- tapply(Signal_vector,ST$ISI_break, na.mean)
          ISI_se_real[subid,] <- tapply(Signal_vector,ST$ISI_break, na.se)
          INT_means_real[subid,] <- tapply(Signal_vector,ST$clickIntensity, na.mean)
          INT_se_real[subid,] <- tapply(Signal_vector,ST$clickIntensity, na.se)
          All_means_real[subid,] <- tapply(Signal_vector,list(ST$ISI_break,ST$clickIntensity), na.mean)
          All_se_real[subid,] <- tapply(Signal_vector,list(ST$ISI_break,ST$clickIntensity), na.se)
          ISI_means_pred[subid,] <- tapply(ST$Pred,ST$ISI_break, na.mean)
          ISI_se_pred[subid,] <- tapply(ST$Pred,ST$ISI_break, na.se)
          INT_means_pred[subid,] <- tapply(ST$Pred,ST$clickIntensity, na.mean)
          INT_se_pred[subid,] <- tapply(ST$Pred,ST$clickIntensity, na.se)
          All_means_pred[subid,] <- tapply(ST$Pred,list(ST$ISI_break,ST$clickIntensity), na.mean)
          All_se_pred[subid,] <- tapply(ST$Pred,list(ST$ISI_break,ST$clickIntensity), na.se)
        } else {
          print(paste("No EEG files to plot for", Subject,Wave,sep = " "))
        }
      }
      
      # Save matrices for this wave since we have them
      save(list=c("ISI_means_real", "ISI_se_real", "INT_means_real","INT_se_real","All_means_real","All_se_real","ISI_means_pred","ISI_se_pred","INT_se_pred","All_means_pred","All_se_pred"), file = paste("~/Dropbox/ampOddClick/r-filesFran/Summary_table_fits/",Type_signal,"_",Parcel,"_",Hemisphere,"_",Wave,sep=""))
  
      # Average across subjects for this wave
      ISI_means_real <- matrix(unlist(ISI_means_real), ncol = 5, byrow = FALSE)
      ISI_means_real <- colMeans(ISI_means_real)
      ISI_se_real <- matrix(unlist(ISI_se_real), ncol = 5, byrow = FALSE)
      ISI_se_real <- colMeans(ISI_se_real)
      INT_means_real <- matrix(unlist(INT_means_real), ncol = 3, byrow = FALSE)
      INT_means_real <- colMeans(INT_means_real)
      INT_se_real <- matrix(unlist(INT_se_real), ncol = 3, byrow = FALSE)
      INT_se_real <- colMeans(INT_se_real)
      All_means_real <- matrix(unlist(All_means_real), ncol = 15, byrow = FALSE)
      All_means_real <- colMeans(All_means_real)
      All_se_real <- matrix(unlist(All_se_real), ncol = 15, byrow = FALSE)
      All_se_real <- colMeans(All_se_real)
      ISI_means_pred <- matrix(unlist(ISI_means_pred), ncol = 5, byrow = FALSE)
      ISI_means_pred <- colMeans(ISI_means_pred)
      ISI_se_pred <- matrix(unlist(ISI_se_pred), ncol = 5, byrow = FALSE)
      ISI_se_pred <- colMeans(ISI_se_pred)
      INT_means_pred <- matrix(unlist(INT_means_pred), ncol = 3, byrow = FALSE)
      INT_means_pred <- colMeans(INT_means_pred)
      INT_se_pred <- matrix(unlist(INT_se_pred), ncol = 3, byrow = FALSE)
      INT_se_pred <- colMeans(INT_se_pred)
      All_means_pred <- matrix(unlist(All_means_pred), ncol = 15, byrow = FALSE)
      All_means_pred <- colMeans(All_means_pred)
      All_se_pred <- matrix(unlist(All_se_pred), ncol = 15, byrow = FALSE)
      All_se_pred <- colMeans(All_se_pred)

  # Create structure to plot (INTENSITY)
  dfR<-data.frame(Mean=INT_means_real,
                  se=INT_se_real,
                  Intensity=c('65dB','75dB','85dB'))
  dfP<-data.frame(Mean=c(INT_means_pred),
                  se=c(INT_se_pred),
                  Intensity=c('65dB','75dB','85dB'))
  # Plot (INTENSITY)
  P <- ggplot()+
    geom_point(data = dfR, aes(x = Intensity, y = Mean, color = "green"), shape = 21)+
    geom_errorbar(data = dfR, aes(x = Intensity, ymin=Mean-se, ymax=Mean+se, color = "green"), width=.1) +
    geom_point(data = dfP, aes(x = Intensity, y = Mean, color = "red"), shape = 21, show.legend = FALSE)+
    geom_errorbar(data = dfP, aes(x = Intensity, ymin=Mean-se, ymax=Mean+se, color = "red"), width=.1) +
    scale_color_manual(values = c("green","red"), labels = c("Real","Predicted"), name = "") + 
    ggtitle(paste(Wave,sep=""))
  P <- P + theme(legend.position = "none") + theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank())
  nam <- paste("P", w, sep = "")
  assign(nam, P)
  
  # Create structure to plot (ISI)
  dfR<-data.frame(Mean=ISI_means_real,
                  se=ISI_se_real,
                  ISI_range = c("0.25-0.5s","0.5-1s","1-2s","2-4s","4-8s"))
  dfP<-data.frame(Mean=ISI_means_pred,
                  se=ISI_se_pred,
                  ISI_range = c("0.25-0.5s","0.5-1s","1-2s","2-4s","4-8s"))
  # Now plot (ISI)
  P <- ggplot()+
    geom_point(data = dfR, aes(x = ISI_range, y = Mean, color = "green"), shape = 21)+
    geom_errorbar(data = dfR, aes(x = ISI_range, ymin=Mean-se, ymax=Mean+se, color = "green"), width=.1) +
    geom_point(data = dfP, aes(x = ISI_range, y = Mean, color = "red"), shape = 21, show.legend = FALSE)+
    geom_errorbar(data = dfP, aes(x = ISI_range, ymin=Mean-se, ymax=Mean+se, color = "red"), width=.1) +
    scale_color_manual(values = c("green","red"), labels = c("Real","Predicted"), name = "")
  P <- P + theme(legend.position = "none") + theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank()) + 
    theme(axis.text.x=element_blank()) + theme(axis.text.x = element_text(angle = 75, hjust=1))
  nam <- paste("P", w+7, sep = "")
  assign(nam, P)
}

Final_plot_sub <- P1+P2+P3+P4+P5+P6+P7+P8+P9+P10+P11+P12+P13+P14+
  plot_annotation(paste(Type_signal,"Average",sep=" "),theme=theme(plot.title=element_text(hjust=0.5)))+
  plot_layout(ncol = 7) + theme(axis.title=element_text(size=12))
# Save the plot file
ggsave(
  paste(Type_signal,"_Average",".png",sep=""),
  plot = Final_plot_sub,
  device = "png",
  path = "~/Dropbox/ampOddClick/r-filesFran/Figures/Model_fit/EEG/ampFUN/",
  scale = 2,
  width = 15,
  height = 5,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
)

