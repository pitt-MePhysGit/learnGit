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
Generate_files <- "YES" # Already computed

if (strcmp(Generate_files, 'YES')) {
# Define mega_variables to plot all subjects together (to see differences across them in fitting/correlation values)
Header_mega_variable_param_v1 <-c('P0_1','P0_2','P0_3','P0_4','P0_5','Na_1','Na_2','Na_3','Na_4','Na_5',
                         'Pa_1','Pa_2','Pa_3','Pa_4','Pa_5','Nb_1','Nb_2','Nb_3','Nb_4','Nb_5',
                         'Pb_1','Pb_2','Pb_3','Pb_4','Pb_5','N1_1','N1_2','N1_3','N1_4','N1_5',
                         'P2_1','P2_2','P2_3','P2_4','P2_5')
Mega_variable_param_v1 <- matrix(list(), length(Subject_list), length(Header_mega_variable_param_v1))
colnames(Mega_variable_param_v1) <- Header_mega_variable_param_v1;
rownames(Mega_variable_param_v1) <- Subject_list;
Header_mega_variable_param_v2 <-c('P0_1','Na_1','Pa_1','Na_1','Pb_1','N1_1','P2_1',
                                  'P0_2','Na_2','Pa_2','Na_2','Pb_2','N1_2','P2_2',
                                  'P0_3','Na_3','Pa_3','Na_3','Pb_3','N1_3','P2_3',
                                  'P0_4','Na_4','Pa_4','Na_4','Pb_4','N1_4','P2_4',
                                  'P0_5','Na_5','Pa_5','Na_5','Pb_5','N1_5','P2_5')
Mega_variable_param_v2 <- matrix(list(), length(Subject_list), length(Header_mega_variable_param_v2))
colnames(Mega_variable_param_v2) <- Header_mega_variable_param_v2;
rownames(Mega_variable_param_v2) <- Subject_list;

Header_mega_variable_corr <-c('P0','Na','Pa','Nb','Pb','N1','P2')
Mega_variable_corr <- matrix(list(), length(Subject_list), length(Header_mega_variable_corr))
colnames(Mega_variable_corr) <- Header_mega_variable_corr;
rownames(Mega_variable_corr) <- Subject_list;

# EEG first
Type_signal <- 'EEG'
for (subid in 1:length(Subject_list)) {
  Subject <- Subject_list[subid]
  pos_col <- 1
  for (w in 1:length(Wave_list)) {
    Wave <- Wave_list[w]
    Hemisphere <- 'NA' # Irrelevant here
    Parcel <- 'NA' # Irrelevant here
    
    # Run gradient descent if tst file does not exist
    if (file.exists(paste("/Volumes/EEGCyno/RES/ampOddClick/STPSP/STPSPFit",Subject, "22223344_5566_zm9999_all_DV_ampFUN_p",Subject,Parcel,Hemisphere,Wave,"gd.rda",sep="_"))){
      # Load the model fit variable
      load(paste("/Volumes/EEGCyno/RES/ampOddClick/STPSP/STPSPFit",Subject, "22223344_5566_zm9999_all_DV_ampFUN_p",Subject,Parcel,Hemisphere,Wave,"gd.rda",sep="_"))
      # Extract fit values
      Mega_variable_param_v1[[subid,pos_col]] <-  modelFit[['par']][1]
      Mega_variable_param_v1[[subid,pos_col+1]] <-  modelFit[['par']][2]
      Mega_variable_param_v1[[subid,pos_col+2]] <-  modelFit[['par']][3]
      Mega_variable_param_v1[[subid,pos_col+3]] <-  modelFit[['par']][4]
      Mega_variable_param_v1[[subid,pos_col+4]] <-  modelFit[['par']][5]
      Mega_variable_corr[[subid,w]] <-  modelFit[['value']][1]
      Mega_variable_param_v2[[subid,w]] <-  modelFit[['par']][1]
      Mega_variable_param_v2[[subid,w+7]] <-  modelFit[['par']][2]
      Mega_variable_param_v2[[subid,w+14]] <-  modelFit[['par']][3]
      Mega_variable_param_v2[[subid,w+21]] <-  modelFit[['par']][4]
      Mega_variable_param_v2[[subid,w+28]] <-  modelFit[['par']][5]
      # Prepare columns for next iteration
      pos_col <- pos_col + 5
    } else {
      print(paste("No EEG files to plot for", Subject,Wave,sep = " "))
    }
  }
}
# Save variables
save(Mega_variable_param_v1, file = paste("~/Dropbox/ampOddClick/r-filesFran/Figures/Summary_parameters/EEG_Mega_variable_param_v1",sep=""))
save(Mega_variable_param_v2, file = paste("~/Dropbox/ampOddClick/r-filesFran/Figures/Summary_parameters/EEG_Mega_variable_param_v2",sep=""))
save(Mega_variable_corr, file = paste("~/Dropbox/ampOddClick/r-filesFran/Figures/Summary_parameters/EEG_Mega_variable_corr",sep=""))
# Make plots and save plots
# Mega_variable_param_v1
if (1 == 1) { # Silly way to compress a long section
subj_long <- c()
pos_emp <- 0
for (subid in 1:length(Subject_list)) {
  for (i in 1:length(Header_mega_variable_param_v1)) {
  subj_long[pos_emp+i] <- Subject_list[subid]
  }
  pos_emp <- pos_emp + length(Header_mega_variable_param_v1)
}
val_long <- c()
pos_emp <- 0
Matrix_x <- matrix(unlist(Mega_variable_param_v1), ncol = 35)
for (subid in 1:length(Subject_list)) {
  for (i in 1:length(Header_mega_variable_param_v1)) {
    val_long[pos_emp+i] <- Matrix_x[subid,i]
  }
  pos_emp <- pos_emp + length(Header_mega_variable_param_v1)
}

Plot_long_format <- data.frame (Participants  = subj_long,
                    Factors_Waves = Header_mega_variable_param_v1,
                    Factor_value = val_long,
                    order = c(1:980)
  )
plot_heatmap <- ggplot(Plot_long_format, aes(x = reorder(Factors_Waves,order), y = Participants, fill = Factor_value)) +
  geom_tile() +
  scale_fill_gradient(low="black", high="white")

ggsave(
  paste("Parameters_",Type_signal,"_",Parcel,"_",Hemisphere,"_v1.png",sep=""),
  plot = plot_heatmap,
  device = "png",
  path = "~/Dropbox/ampOddClick/r-filesFran/Figures/Summary_parameters/",
  scale = 5,
  width = NA,
  height = NA,
  units = "cm",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
)

}
# Mega_variable_param_v2
if (1 == 1) { # Silly way to compress a long section
  subj_long <- c()
  pos_emp <- 0
  for (subid in 1:length(Subject_list)) {
    for (i in 1:length(Header_mega_variable_param_v2)) {
      subj_long[pos_emp+i] <- Subject_list[subid]
    }
    pos_emp <- pos_emp + length(Header_mega_variable_param_v2)
  }
  val_long <- c()
  pos_emp <- 0
  Matrix_x <- matrix(unlist(Mega_variable_param_v2), ncol = 35)
  for (subid in 1:length(Subject_list)) {
    for (i in 1:length(Header_mega_variable_param_v2)) {
      val_long[pos_emp+i] <- Matrix_x[subid,i]
    }
    pos_emp <- pos_emp + length(Header_mega_variable_param_v2)
  }
  
  Plot_long_format <- data.frame (Participants  = subj_long,
                                  Factors_Waves = Header_mega_variable_param_v2,
                                  Factor_value = val_long,
                                  order = c(1:980)
  )
  plot_heatmap <- ggplot(Plot_long_format, aes(x = reorder(Factors_Waves,order), y = Participants, fill = Factor_value)) +
    geom_tile() +
    scale_fill_gradient(low="black", high="white")
  
  ggsave(
    paste("Parameters_",Type_signal,"_",Parcel,"_",Hemisphere,"_v2.png",sep=""),
    plot = plot_heatmap,
    device = "png",
    path = "~/Dropbox/ampOddClick/r-filesFran/Figures/Summary_parameters/",
    scale = 5,
    width = NA,
    height = NA,
    units = "cm",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
  )
}
# Mega_variable_corr
if (1 == 1) { # Silly way to compress a long section
  subj_long <- c()
  pos_emp <- 0
  for (subid in 1:length(Subject_list)) {
    for (i in 1:length(Header_mega_variable_corr)) {
      subj_long[pos_emp+i] <- Subject_list[subid]
    }
    pos_emp <- pos_emp + length(Header_mega_variable_corr)
  }
  val_long <- c()
  pos_emp <- 0
  Matrix_x <- matrix(unlist(Mega_variable_corr), ncol = 7)
  for (subid in 1:length(Subject_list)) {
    for (i in 1:length(Header_mega_variable_corr)) {
      val_long[pos_emp+i] <- Matrix_x[subid,i]
    }
    pos_emp <- pos_emp + length(Header_mega_variable_corr)
  }
  
  Plot_long_format <- data.frame (Participants  = subj_long,
                                  Waves = Header_mega_variable_corr,
                                  Corr_value = val_long,
                                  order = c(1:196)
  )
  plot_heatmap <- ggplot(Plot_long_format, aes(x = reorder(Waves,order), y = Participants, fill = Corr_value)) +
    geom_tile() +
    scale_fill_gradient(low="black", high="white")
  
  ggsave(
    paste("Parameters_",Type_signal,"_",Parcel,"_",Hemisphere,"_corr.png",sep=""),
    plot = plot_heatmap,
    device = "png",
    path = "~/Dropbox/ampOddClick/r-filesFran/Figures/Summary_parameters/",
    scale = 5,
    width = NA,
    height = NA,
    units = "cm",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
  )
}

# Now same but for sources

Subject_list <- c('2193','2199','2235','2238','2250','2320','2259','2262','2274',
                  '2275','H001','2214','2212','2285','2337','2372','2389','H002',
                  '2299','2406','2426','2456','H005',
                  '2005','2355','2387','2288')  # '2427' AND 2494' does not have source files

# Define mega_variables to plot all subjects together (to see differences across them in fitting/correlation values)
Header_mega_variable_param_v1 <-c('P0_1','P0_2','P0_3','P0_4','P0_5','Na_1','Na_2','Na_3','Na_4','Na_5',
                                  'Pa_1','Pa_2','Pa_3','Pa_4','Pa_5','Nb_1','Nb_2','Nb_3','Nb_4','Nb_5',
                                  'Pb_1','Pb_2','Pb_3','Pb_4','Pb_5','N1_1','N1_2','N1_3','N1_4','N1_5',
                                  'P2_1','P2_2','P2_3','P2_4','P2_5')
Mega_variable_param_v1 <- matrix(list(), length(Subject_list), length(Header_mega_variable_param_v1))
colnames(Mega_variable_param_v1) <- Header_mega_variable_param_v1;
rownames(Mega_variable_param_v1) <- Subject_list;
Header_mega_variable_param_v2 <-c('P0_1','Na_1','Pa_1','Na_1','Pb_1','N1_1','P2_1',
                                  'P0_2','Na_2','Pa_2','Na_2','Pb_2','N1_2','P2_2',
                                  'P0_3','Na_3','Pa_3','Na_3','Pb_3','N1_3','P2_3',
                                  'P0_4','Na_4','Pa_4','Na_4','Pb_4','N1_4','P2_4',
                                  'P0_5','Na_5','Pa_5','Na_5','Pb_5','N1_5','P2_5')
Mega_variable_param_v2 <- matrix(list(), length(Subject_list), length(Header_mega_variable_param_v2))
colnames(Mega_variable_param_v2) <- Header_mega_variable_param_v2;
rownames(Mega_variable_param_v2) <- Subject_list;

Header_mega_variable_corr <-c('P0','Na','Pa','Nb','Pb','N1','P2')
Mega_variable_corr <- matrix(list(), length(Subject_list), length(Header_mega_variable_corr))
colnames(Mega_variable_corr) <- Header_mega_variable_corr;
rownames(Mega_variable_corr) <- Subject_list;

# Then do Sources
Type_signal <- 'Source'
for (parc in 1:length(Parcel_list)) {
  Parcel <- Parcel_list[parc]
  for (hem in 1:length(Hemisphere_list)) {
    Hemisphere <- Hemisphere_list[hem]
    for (subid in 1:length(Subject_list)) {
      Subject <- Subject_list[subid]
      pos_col <- 1
    for (w in 1:length(Wave_list)) {
      Wave <- Wave_list[w]
      if (file.exists(paste(paste("/Volumes/EEGCyno/RES/ampOddClick/STPSP/STPSPFit",Subject, "22223344_5566_zm9999_all_DV_ampFUN_p",Subject,Parcel,Hemisphere,Wave,"gd.rda",sep="_")))) {
        # Load the model fit variable
        load(paste("/Volumes/EEGCyno/RES/ampOddClick/STPSP/STPSPFit",Subject, "22223344_5566_zm9999_all_DV_ampFUN_p",Subject,Parcel,Hemisphere,Wave,"gd.rda",sep="_"))
        # Extract fit values
        Mega_variable_param_v1[[subid,pos_col]] <-  modelFit[['par']][1]
        Mega_variable_param_v1[[subid,pos_col+1]] <-  modelFit[['par']][2]
        Mega_variable_param_v1[[subid,pos_col+2]] <-  modelFit[['par']][3]
        Mega_variable_param_v1[[subid,pos_col+3]] <-  modelFit[['par']][4]
        Mega_variable_param_v1[[subid,pos_col+4]] <-  modelFit[['par']][5]
        Mega_variable_corr[[subid,w]] <-  modelFit[['value']][1]
        Mega_variable_param_v2[[subid,w]] <-  modelFit[['par']][1]
        Mega_variable_param_v2[[subid,w+7]] <-  modelFit[['par']][2]
        Mega_variable_param_v2[[subid,w+14]] <-  modelFit[['par']][3]
        Mega_variable_param_v2[[subid,w+21]] <-  modelFit[['par']][4]
        Mega_variable_param_v2[[subid,w+28]] <-  modelFit[['par']][5]
        # Prepare columns for next iteration
        pos_col <- pos_col + 5
      } else {
        print(paste("No Source files to plot for", Subject,Wave,Hemisphere,Parcel,sep = " "))
      }
    }
    }
    save(Mega_variable_param_v1, file = paste("~/Dropbox/ampOddClick/r-filesFran/Figures/Summary_parameters/",Parcel,"_",Hemisphere,"_Mega_variable_param_v1",sep=""))
    save(Mega_variable_param_v2, file = paste("~/Dropbox/ampOddClick/r-filesFran/Figures/Summary_parameters/",Parcel,"_",Hemisphere,"_Mega_variable_param_v2",sep=""))
    save(Mega_variable_corr, file = paste("~/Dropbox/ampOddClick/r-filesFran/Figures/Summary_parameters/",Parcel,"_",Hemisphere,"_Mega_variable_corr",sep=""))  
    
    # Make plots and save plots
    # Mega_variable_param_v1
    if (1 == 1) { # Silly way to compress a long section
      subj_long <- c()
      pos_emp <- 0
      for (subid in 1:length(Subject_list)) {
        for (i in 1:length(Header_mega_variable_param_v1)) {
          subj_long[pos_emp+i] <- Subject_list[subid]
        }
        pos_emp <- pos_emp + length(Header_mega_variable_param_v1)
      }
      val_long <- c()
      pos_emp <- 0
      Matrix_x <- matrix(unlist(Mega_variable_param_v1), ncol = 35)
      for (subid in 1:length(Subject_list)) {
        for (i in 1:length(Header_mega_variable_param_v1)) {
          val_long[pos_emp+i] <- Matrix_x[subid,i]
        }
        pos_emp <- pos_emp + length(Header_mega_variable_param_v1)
      }
      
      Plot_long_format <- data.frame (Participants  = subj_long,
                                      Factors_Waves = Header_mega_variable_param_v1,
                                      Factor_value = val_long,
                                      order = c(1:945)
      )
      plot_heatmap <- ggplot(Plot_long_format, aes(x = reorder(Factors_Waves,order), y = Participants, fill = Factor_value)) +
        geom_tile() +
        scale_fill_gradient(low="black", high="white")
      
      ggsave(
        paste("Parameters_",Type_signal,"_",Parcel,"_",Hemisphere,"_v1.png",sep=""),
        plot = plot_heatmap,
        device = "png",
        path = "~/Dropbox/ampOddClick/r-filesFran/Figures/Summary_parameters/",
        scale = 5,
        width = NA,
        height = NA,
        units = "cm",
        dpi = 300,
        limitsize = TRUE,
        bg = NULL,
      )
      
    }
    # Mega_variable_param_v2
    if (1 == 1) { # Silly way to compress a long section
      subj_long <- c()
      pos_emp <- 0
      for (subid in 1:length(Subject_list)) {
        for (i in 1:length(Header_mega_variable_param_v2)) {
          subj_long[pos_emp+i] <- Subject_list[subid]
        }
        pos_emp <- pos_emp + length(Header_mega_variable_param_v2)
      }
      val_long <- c()
      pos_emp <- 0
      Matrix_x <- matrix(unlist(Mega_variable_param_v2), ncol = 35)
      for (subid in 1:length(Subject_list)) {
        for (i in 1:length(Header_mega_variable_param_v2)) {
          val_long[pos_emp+i] <- Matrix_x[subid,i]
        }
        pos_emp <- pos_emp + length(Header_mega_variable_param_v2)
      }
      
      Plot_long_format <- data.frame (Participants  = subj_long,
                                      Factors_Waves = Header_mega_variable_param_v2,
                                      Factor_value = val_long,
                                      order = c(1:945)
      )
      plot_heatmap <- ggplot(Plot_long_format, aes(x = reorder(Factors_Waves,order), y = Participants, fill = Factor_value)) +
        geom_tile() +
        scale_fill_gradient(low="black", high="white")
      
      ggsave(
        paste("Parameters_",Type_signal,"_",Parcel,"_",Hemisphere,"_v2.png",sep=""),
        plot = plot_heatmap,
        device = "png",
        path = "~/Dropbox/ampOddClick/r-filesFran/Figures/Summary_parameters/",
        scale = 5,
        width = NA,
        height = NA,
        units = "cm",
        dpi = 300,
        limitsize = TRUE,
        bg = NULL,
      )
    }
    # Mega_variable_corr
    if (1 == 1) { # Silly way to compress a long section
      subj_long <- c()
      pos_emp <- 0
      for (subid in 1:length(Subject_list)) {
        for (i in 1:length(Header_mega_variable_corr)) {
          subj_long[pos_emp+i] <- Subject_list[subid]
        }
        pos_emp <- pos_emp + length(Header_mega_variable_corr)
      }
      val_long <- c()
      pos_emp <- 0
      Matrix_x <- matrix(unlist(Mega_variable_corr), ncol = 7)
      for (subid in 1:length(Subject_list)) {
        for (i in 1:length(Header_mega_variable_corr)) {
          val_long[pos_emp+i] <- Matrix_x[subid,i]
        }
        pos_emp <- pos_emp + length(Header_mega_variable_corr)
      }
      
      Plot_long_format <- data.frame (Participants  = subj_long,
                                      Waves = Header_mega_variable_corr,
                                      Corr_value = val_long,
                                      order = c(1:189)
      )
      plot_heatmap <- ggplot(Plot_long_format, aes(x = reorder(Waves,order), y = Participants, fill = Corr_value)) +
        geom_tile() +
        scale_fill_gradient(low="black", high="white")
      
      ggsave(
        paste("Parameters_",Type_signal,"_",Parcel,"_",Hemisphere,"_corr.png",sep=""),
        plot = plot_heatmap,
        device = "png",
        path = "~/Dropbox/ampOddClick/r-filesFran/Figures/Summary_parameters/",
        scale = 5,
        width = NA,
        height = NA,
        units = "cm",
        dpi = 300,
        limitsize = TRUE,
        bg = NULL,
      )
    }
    
  }
}

Subject_list <- c('2193','2199','2235','2238','2250','2320','2259','2262','2274',
                  '2275','H001','2214','2212','2285','2337','2372','2389','H002',
                  '2299','2406','2426','2456','H005',
                  '2005','2355','2387','2288')  # '2427' AND 2494' does not have source files
## Now create one separate figure per subject to compare parameters and correlations across ROIs
# Only for sources (which have Parcel and Hemisphere dimensions)
Header_mega_variable_ind_sub <-c('P0_1','Na_1','Pa_1','Na_1','Pb_1','N1_1','P2_1',
                                 'P0_2','Na_2','Pa_2','Na_2','Pb_2','N1_2','P2_2',
                                 'P0_3','Na_3','Pa_3','Na_3','Pb_3','N1_3','P2_3',
                                 'P0_4','Na_4','Pa_4','Na_4','Pb_4','N1_4','P2_4',
                                 'P0_5','Na_5','Pa_5','Na_5','Pb_5','N1_5','P2_5')
Rows_mega_variable_ind_sub <- c('A1_L','MBelt_L','LBelt_L','PBelt_L',
                                'A1_R','MBelt_R','LBelt_R','PBelt_R')
Mega_variable_ind_sub <- matrix(list(), length(Rows_mega_variable_ind_sub), length(Header_mega_variable_ind_sub))
colnames(Mega_variable_ind_sub) <- Header_mega_variable_ind_sub;
rownames(Mega_variable_ind_sub) <- Rows_mega_variable_ind_sub;
# Within_subjects
Type_signal <- 'Source'
for (subid in 1:length(Subject_list)) {
Subject <- Subject_list[subid]
pos_h <- 0
for (hem in 1:length(Hemisphere_list)) {
  Hemisphere <- Hemisphere_list[hem]
  for (parc in 1:length(Parcel_list)) {
    Parcel <- Parcel_list[parc]
    for (w in 1:length(Wave_list)) {
      Wave <- Wave_list[w]
      pos_col <- 1
      if (file.exists(paste(paste("/Volumes/EEGCyno/RES/ampOddClick/STPSP/STPSPFit",Subject, "22223344_5566_zm9999_all_DV_ampFUN_p",Subject,Parcel,Hemisphere,Wave,"gd.rda",sep="_")))) {
        # Load the model fit variable
        load(paste("/Volumes/EEGCyno/RES/ampOddClick/STPSP/STPSPFit",Subject, "22223344_5566_zm9999_all_DV_ampFUN_p",Subject,Parcel,Hemisphere,Wave,"gd.rda",sep="_"))
        # Extract fit values
        Mega_variable_ind_sub[[pos_h+parc,w]] <-  modelFit[['par']][1]
        Mega_variable_ind_sub[[pos_h+parc,w+7]] <-  modelFit[['par']][2]
        Mega_variable_ind_sub[[pos_h+parc,w+14]] <-  modelFit[['par']][3]
        Mega_variable_ind_sub[[pos_h+parc,w+21]] <-  modelFit[['par']][4]
        Mega_variable_ind_sub[[pos_h+parc,w+28]] <-  modelFit[['par']][5]
        # Mega_variable_corr[[subid,w]] <-  modelFit[['value']][1]
        # Prepare columns for next iteration
        pos_col <- pos_col + 5
      } else {
        print(paste("No Source files to plot for", Subject,Wave,Hemisphere,Parcel,sep = " "))
      }
    }  
  }
  pos_h <- pos_h + 4
}
save(Mega_variable_ind_sub, file = paste("~/Dropbox/ampOddClick/r-filesFran/Figures/Summary_parameters/Within_subjects/",Subject,"_Mega_variable",sep=""))  
# save(Mega_variable_ind_sub_corr, file = paste("~/Dropbox/ampOddClick/r-filesFran/Figures/Summary_parameters/Within_subjects/",Subject,"_Mega_variable_corr",sep=""))  

# Make plots and save plots
# Mega_variable_ind_sub
if (1 == 1) { # Silly way to compress a long section
  subj_long <- c()
  pos_emp <- 0
  for (subid in 1:length(Rows_mega_variable_ind_sub)) {
    for (i in 1:length(Header_mega_variable_ind_sub)) {
      subj_long[pos_emp+i] <- Rows_mega_variable_ind_sub[subid]
    }
    pos_emp <- pos_emp + length(Header_mega_variable_ind_sub)
  }
  val_long <- c()
  pos_emp <- 0
  Matrix_x <- matrix(unlist(Mega_variable_ind_sub), ncol = 35)
  for (subid in 1:length(Rows_mega_variable_ind_sub)) {
    for (i in 1:length(Header_mega_variable_ind_sub)) {
      val_long[pos_emp+i] <- Matrix_x[subid,i]
    }
    pos_emp <- pos_emp + length(Header_mega_variable_ind_sub)
  }
  
  Plot_long_format <- data.frame (ROIS = subj_long,
                                  Waves = Header_mega_variable_ind_sub,
                                  param_value = val_long,
                                  order = c(1:280)
  )
  plot_heatmap <- ggplot(Plot_long_format, aes(x = reorder(Waves,order), y = ROIS, fill = param_value)) +
    geom_tile() +
    scale_fill_gradient(low="black", high="white")
  
  save(Mega_variable_ind_sub, file = paste("~/Dropbox/ampOddClick/r-filesFran/Figures/Summary_parameters/Within_subjects/",Subject,"_Mega_variable",sep=""))  
  ggsave(
    paste(Type_signal,"_",Subject,"_param.png",sep=""),
    plot = plot_heatmap,
    device = "png",
    path = "~/Dropbox/ampOddClick/r-filesFran/Figures/Summary_parameters/Within_subjects/",
    scale = 5,
    width = NA,
    height = NA,
    units = "cm",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
  )
}

}

}


