function_aoc_LRMfit <- function(ind, linkFileName = 'ampOddclickVProbe_link.txt', fitmod = "", eplist = c('', '', '')){
  
  #ind<-129
  #fitmod<-'LMEC'
  
  baseDir <- '~/Dropbox/'
  library(MASS)
  source('~/Dropbox/r-compile/signalProcessoR/R/signalProcessoR.R')
  source('~/Dropbox/ampOddClick/r-files/functions/prepDFVprobe.R')
  source('~/Dropbox/ampOddClick/r-files/functions/functions_ampOdd_presynapticPlasticity.R')
  
  print(paste("Input to function_aoc_LRMfit, ind: ", ind))

  rda  <- function(){'/Volumes/EEGLab/EEGLab/ampOddClick'}
  trda <- function(){'/Volumes/Drobo5D3/EEG/ampOddClick/rda/'}
  .GlobalEnv$rda <- rda
  .GlobalEnv$trda <- trda
  
  # read the curated file
  df          <- prepDFsignalProcessor("ampOddclick", linkFile = linkFileName)
  if(tolower(linkFileName)=='ampoddclickvprobe_link.txt'){
    df         <- prepDFVprobe(df)
  }
  dfline<-df[ind,]
  
  minISI <- 0.190
  maxISI <- 12.8
  DVscl  <- 5
  show <- FALSE
  DVstr = "active"
  
  tone<-loadSignalProcessor(dfline, project = 'ampOddClick', epoch = "raw", chStr = "ch1")
  
  
  tdir<-paste0('/Volumes/EEGLab/EEGLab/ampOddClick/', dfline$session, '/STPSP/')
  
  if (!dir.exists(tdir)){
    system(paste('mkdir -p ', tdir, sep=''))
  }
  
  
  checkvec<-unlist(strsplit(fitmod, ''))
  for(i in 1:length(checkvec)){
    ldir<-list.files(paste0("/Volumes/EEGLab/EEGLab/ampOddClick/", dfline$session, '/'))
    compexist<-unlist(lapply(ldir, function(x){grepl(paste0('features_', tolower(checkvec[i])), x)}))
    
    if(checkvec[i]%in%c('L', 'E', 'C')){
      cvi<-which(checkvec[i]==c('L', 'E', 'C'))
      ldi<-ldir[which(compexist)] == paste0('features_', c('lfp', 'eeg', 'csd')[cvi], '_', eplist[cvi], '.rda')
      if(any(ldi) & length(ldi)>1){
        compexist[which(compexist)[which(!ldi)]]<-FALSE
      }
    }

    if(sum(compexist)>0){
      for(j in 1:sum(compexist)){
        lrmdat<-tone$xpp[,c('ISI', 'ampInd', 'trialType', 'valTone')]

        rm(xpp)
        load(paste0("/Volumes/EEGLab/EEGLab/ampOddClick/", dfline$session, '/', ldir[which(compexist)[j]]))
      
        valcomp<-(dim(tone$xpp)[2]+1):dim(xpp)[2]
        if(length(valcomp)>0){
          
          for(k in 1:length(valcomp)){
            print(paste0("Fitting: ", checkvec[i], ' - ', ldir[which(compexist)[j]], ' - ', names(xpp)[valcomp[k]]))
      
            fittone<-tone
            fittone$exp<-paste0(dfline$session, '_', strsplit(ldir[which(compexist)[j]], '[.]')[[1]][1], '_', names(xpp)[valcomp[k]])
            fittone$xpp$exp<-fittone$exp

            fittone$xpp$powEOG<-0
            fittone$xpp$toneType<-0

            fittone$xpp$active<-xpp[,valcomp[k]]

            fittone$xpp$p2p<-2*(fittone$xpp$p2p>quantile(xpp$p2p, 0.5)+6*diff(quantile(xpp$p2p, c(0.25, 0.75))) |
            fittone$xpp$active>quantile(fittone$xpp$active, 0.5)+6*diff(quantile(fittone$xpp$active, c(0.25, 0.75))) |
            fittone$xpp$active<quantile(fittone$xpp$active, 0.5)-6*diff(quantile(fittone$xpp$active, c(0.25, 0.75))))
            
            #DVscl for outliers; done in STPSP
            #fittone$xpp$active>quantile(fittone$xpp$active, 0.5)+3*diff(quantile(fittone$xpp$active, c(0.25, 0.75)))
            
            #dumping nested list storage of full results for a few vals
            trlrm<-callSTPSP(fittone, gd = T, startFromScratch = TRUE, DVstr = "active", who = 'all', BYstr = 'exp', FUN = ampFUN, minISI = minISI, maxISI = maxISI, DVsdscale = DVscl, show = F, maxp2p = 1, verbose = TRUE, path = FALSE)

            lrmybin<-trlrm$res$ybin

            lrmdat[[names(xpp)[valcomp[k]]]]<-xpp[,valcomp[k]]
            lrmdat[[paste0(names(xpp)[valcomp[k]], "_R")]]<-NA
            lrmdat[[paste0(names(xpp[valcomp[k]]), "_fit")]]<-NA
            lrmdat[[paste0(names(xpp)[valcomp[k]], "_R")]][lrmybin]<-trlrm$res$lm2$model[['R']]
            lrmdat[[paste0(names(xpp[valcomp[k]]), "_fit")]][lrmybin]<-trlrm$res$lm2$fitted.values

            if(k==1){
              lrmpar<-data.frame(trlrm$modelFit$par)
              names(lrmpar)<-names(xpp)[valcomp[k]]
            }else{
              lrmpar[[names(xpp)[valcomp[k]]]]<-trlrm$modelFit$par
            }
          }

          tpath<-paste0("/Volumes/EEGLab/EEGLab/ampOddClick/", dfline$session, '/STPSP/', strsplit(ldir[which(compexist)[j]], '[.]')[[1]][1], '_fits.rda')
          save(lrmdat, lrmpar, file = tpath)
        
        }
      }
    }
  }

  print(c("LRM success: ", as.character(dfline$session)))
  
  #how readable?
  #t<-proc.time()
  #load('/Volumes/EEGLab/EEGLab/ampOddClick/Walter_20220419_1132_zm5000/LMEC_fits.rda')
  #proc.time()-t
}
