## Plot parameters gradient descend across subjects and within subjects to compare across ROIs
## Plot correlation values as well (to see differences across subjects)

rm(list=ls())
## Add packages -----------------------
# install.packages('patchwork')
# install.packages('reshape')
library(patchwork)
library(reshape)
library(pracma)
library(lattice)
library(ggplot2)
#library(corrplot)
rda <- function(){ return('/Volumes/EEGCyno/RES/ampOddClick/') }
source('~/Dropbox/r-compile/rUtils/R/rUtils.R')
source('~/Dropbox/r-compile/ePhysLab/R/ePhysLab.R')
source('~/Dropbox/ampOddClick/r-files/functions/functions_ampOdd_presynapticPlasticity.R')
source('~/Dropbox/ampOddClick/r-filesFran/load_ampOddClick_Human.R')

# Define variables
Subject_list <- c('2193','2199','2235','2238','2250','2320','2259','2262','2274',
                  '2275','H001','2214','2212','2285','2337','2372','2389','H002',
                  '2299','2406','2426','2427','2456','H005',
                  '2005','2355','2387','2288')  # '2494' does not have files for some reason
Wave_list <- c('P0','Na','Pa','Nb','P50','N1','P2')
Type_signal_list <- c('EEG','Source')
Hemisphere_list <- c('L','R') # L or R, only if choosing sources as type of signal
Parcel_list <- c('A1','MBelt','LBelt','PBelt') # Only if choosing sources as type of signal

# Load data
load("~/Dropbox/ampOddClick/r-filesFran/Figures/Summary_parameters/ampFUN2q/EEG_Mega_variable_param_v2")
  
# Define mega_variables to plot all subjects together (to see differences across them in fitting/correlation values)
# Header_mega_variable_param_v2 <-c('P0_1','Na_1','Pa_1','Na_1','Pb_1','N1_1','P2_1',
#                                   'P0_2','Na_2','Pa_2','Na_2','Pb_2','N1_2','P2_2',
#                                   'P0_3','Na_3','Pa_3','Na_3','Pb_3','N1_3','P2_3',
#                                   'P0_4','Na_4','Pa_4','Na_4','Pb_4','N1_4','P2_4',
#                                   'P0_5','Na_5','Pa_5','Na_5','Pb_5','N1_5','P2_5',
#                                   'P0_6','Na_6','Pa_6','Na_6','Pb_6','N1_6','P2_6',
#                                   'P0_7','Na_7','Pa_7','Na_7','Pb_7','N1_7','P2_7')
# Mega_variable_param_v2 <- matrix(list(), length(Subject_list), length(Header_mega_variable_param_v2))
# colnames(Mega_variable_param_v2) <- Header_mega_variable_param_v2;
# rownames(Mega_variable_param_v2) <- Subject_list;

Matrix_x <- matrix(unlist(Mega_variable_param_v2), ncol = 49)
averages_sub <- colMeans(Matrix_x)
stdev_sub <- apply(Matrix_x, 2, sd) 
sterr_sub <- stdev_sub/sqrt(length(Subject_list))

# Original
#    parm[1]              time contant of STPSP model
#    parm[2]              release probability of STPSP model
#    parm[3]              mean non-linearity
#    parm[4]              sd non-linearity
#    parm[5]              offset

# New model
#    parm[1]              1/slow time constant of STPSP model
#    parm[2]              slow release probability/depleted fraction of STPSP model
#    parm[3]              1/fast-1/slow time constant of STPSP model
#    parm[4]              fast release probability/depleted fraction of STPSP model
#    parm[5]              mean non-linearity
#    parm[6]              sd non-linearity
#    parm[7]              offset

# Slow time constant of STPSP model
avg_1 <- c(averages_sub[5],averages_sub[6],averages_sub[7])
stderr_1 <- c(sterr_sub[5],sterr_sub[6],sterr_sub[7])

dfR<- data.frame(Mean=avg_1,
                se=stderr_1,
                Intensity=c('P50','N1','P2'),
                order = c(1:3)
                )

Plot_1 <- ggplot()+
  geom_col(data = dfR, aes(x = reorder(Intensity,order), y = Mean))+
  geom_errorbar(data = dfR, aes(x = Intensity, ymin=Mean-se, ymax=Mean+se, color = "green"), width=.1) +
  ggtitle(paste("Slow recovery time constant",sep=""))
Plot_1 <- Plot_1 + theme(legend.position = "none") + theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank())

# Cumulative release probability
avg_1 <- c(averages_sub[12],averages_sub[13],averages_sub[14])
stderr_1 <- c(sterr_sub[12],sterr_sub[13],sterr_sub[14])

dfR<- data.frame(Mean=avg_1,
                 se=stderr_1,
                 Intensity=c('P50','N1','P2'),
                 order = c(1:3)
)

Plot_2 <- ggplot()+
  geom_col(data = dfR, aes(x = reorder(Intensity,order), y = Mean))+
  geom_errorbar(data = dfR, aes(x = Intensity, ymin=Mean-se, ymax=Mean+se, color = "green"), width=.1) +
  ggtitle(paste("Slow cumulative release probability",sep=""))
Plot_2 <- Plot_2 + theme(legend.position = "none") + theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank())

# Fast time constant
avg_1 <- c(averages_sub[19],averages_sub[20],averages_sub[21])
stderr_1 <- c(sterr_sub[19],sterr_sub[20],sterr_sub[21])

dfR<- data.frame(Mean=avg_1,
                 se=stderr_1,
                 Intensity=c('P50','N1','P2'),
                 order = c(1:3)
)

Plot_3 <- ggplot()+
  geom_col(data = dfR, aes(x = reorder(Intensity,order), y = Mean))+
  geom_errorbar(data = dfR, aes(x = Intensity, ymin=Mean-se, ymax=Mean+se, color = "green"), width=.1) +
  ggtitle(paste("Fast time constant",sep=""))
Plot_3 <- Plot_3 + theme(legend.position = "none") + theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank())

# Fast cumulative release probability
avg_1 <- c(averages_sub[26],averages_sub[27],averages_sub[28])
stderr_1 <- c(sterr_sub[26],sterr_sub[27],sterr_sub[28])

dfR<- data.frame(Mean=avg_1,
                 se=stderr_1,
                 Intensity=c('P50','N1','P2'),
                 order = c(1:3)
)

Plot_4 <- ggplot()+
  geom_col(data = dfR, aes(x = reorder(Intensity,order), y = Mean))+
  geom_errorbar(data = dfR, aes(x = Intensity, ymin=Mean-se, ymax=Mean+se, color = "green"), width=.1) +
  ggtitle(paste("Fast cummulative release probability",sep=""))
Plot_4 <- Plot_4 + theme(legend.position = "none") + theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank())

# Mean non-linearity
avg_1 <- c(averages_sub[33],averages_sub[34],averages_sub[35])
stderr_1 <- c(sterr_sub[33],sterr_sub[34],sterr_sub[35])

dfR<- data.frame(Mean=avg_1,
                 se=stderr_1,
                 Intensity=c('P50','N1','P2'),
                 order = c(1:3)
)

Plot_5 <- ggplot()+
  geom_col(data = dfR, aes(x = reorder(Intensity,order), y = Mean))+
  geom_errorbar(data = dfR, aes(x = Intensity, ymin=Mean-se, ymax=Mean+se, color = "green"), width=.1) +
  ggtitle(paste("Mean non-linearity",sep=""))
Plot_5 <- Plot_5 + theme(legend.position = "none") + theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank())

# SD non-linearity
avg_1 <- c(averages_sub[40],averages_sub[41],averages_sub[42])
stderr_1 <- c(sterr_sub[40],sterr_sub[41],sterr_sub[42])

dfR<- data.frame(Mean=avg_1,
                 se=stderr_1,
                 Intensity=c('P50','N1','P2'),
                 order = c(1:3)
)

Plot_6 <- ggplot()+
  geom_col(data = dfR, aes(x = reorder(Intensity,order), y = Mean))+
  geom_errorbar(data = dfR, aes(x = Intensity, ymin=Mean-se, ymax=Mean+se, color = "green"), width=.1) +
  ggtitle(paste("SD non-linearity",sep=""))
Plot_6 <- Plot_6 + theme(legend.position = "none") + theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank())


# Offset
avg_1 <- c(averages_sub[47],averages_sub[48],averages_sub[49])
stderr_1 <- c(sterr_sub[47],sterr_sub[48],sterr_sub[49])

dfR<- data.frame(Mean=avg_1,
                 se=stderr_1,
                 Intensity=c('P50','N1','P2'),
                 order = c(1:3)
)

Plot_7 <- ggplot()+
  geom_col(data = dfR, aes(x = reorder(Intensity,order), y = Mean))+
  geom_errorbar(data = dfR, aes(x = Intensity, ymin=Mean-se, ymax=Mean+se, color = "green"), width=.1) +
  ggtitle(paste("Offset",sep=""))
Plot_7 <- Plot_7 + theme(legend.position = "none") + theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank())




# Merge plots together

Final_plot <- Plot_1 + Plot_2 + Plot_3 + Plot_4 + Plot_5 + Plot_6 + Plot_7 

ggsave(
  paste("GAVR_parameters_EEG.png",sep=""),
  plot = Final_plot,
  device = "png",
  path = "~/Dropbox/ampOddClick/r-filesFran/Figures/Summary_parameters/ampFUN2q/",
  scale = 3,
  width = NA,
  height = NA,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
)



