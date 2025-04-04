load_ampOddClick_Human <- function(Wave='P2', Subject='2005',  Type_signal='EEG', Hemisphere='L',Parcel='A1'){
  
  if (strcmp(Wave, 'N1') | strcmp(Wave, 'P2'))  {
    sig_latency <- 'LLR'
  } else {
    sig_latency <- 'MLR'
  }
  
  ## Load single trial data and basic QC -----------------------------------------------------------
  single_data_files_dir <- '/Volumes/EEGCyno/EEG/ampOddClickHuman/Shared_data_single_trials'
  if (strcmp(Type_signal, 'EEG')) {
    ST <- read.table( paste(single_data_files_dir,'/Amp_singl_trials_EEG_',sig_latency,'.txt',sep=''), header=T,sep = "," )
  } else {
    ST <- read.table( paste(single_data_files_dir,'/Single_trials_',Parcel,'_',Hemisphere,'_',sig_latency,'.txt',sep=''), header=T,sep = "\t" )
  }
  # Crop it to contain only the subject at hand
  ST <- subset.data.frame(ST,ST$subjID == Subject)
 
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
  
  isi_breaks <- c(.25, .5, 1, 2, 4, 22)
  ST$ISI_break <- cut(ST$ISI,breaks=isi_breaks)
  ST$isiOct    <- ST$ISI_break
  
  
  
  ST$p2p     <- 50
  ST$powEOG  <- 50
  ST$tmp <- ST$ISI > 15
  ST$exp <- cumsum(ST$tmp)
  
  ## create DV 
  ST$DV <- ST[[ paste(Type_signal, '_', Wave, sep='') ]] 
  ST$p2p[ is.na(ST$DV) ]     <- 9999
  ST$DV[ is.na(ST$DV) ]      <- 0
  
  
  ST$trialType <- 1 # random, not predicatble
  
  
  ST$clickIntensityLag1 <- my.lag(ST$clickIntensity,1)
  ST$toneType           <- c(2,1)[1 + (ST$clickIntensity == ST$clickIntensityLag1)]
  
  ST$ampInd <- ST$clickIntensity
  ST$dB     <- c(65,75,85)[ST$ampInd]
  ST$deltaAmpInd <- ST$ampInd - my.lag(ST$ampInd,1)
  
  tn     <- c()
  tn$xpp <- ST
  tn$exp <- paste(Subject,'_22223344_5566_zm9999', sep='')
  
  
  return(tn)
}