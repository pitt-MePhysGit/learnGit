getValidTrialsAOC <- function(tone,DVstr='DV',maxISI=25,
                                minISI=0.195,
                                maxp2p=Inf,
                                DVsdscale=6,
                                maxPowEOG=Inf,
                                who='all'){
  
  
  trialTypeIndex <- c(1)
  if (who=='all')
    trialTypeIndex <- c(1,2)
  if (who=='regular')
    trialTypeIndex <- c(2)
  
  
  mDV       <- quantile( tone$xpp[[DVstr]], prob=c(0.25,0.5,0.75) )
  valRng    <- mDV[2] + c(-DVsdscale,DVsdscale) * diff(mDV) 
  
  ybin <- ( tone$xpp$p2p<maxp2p                      &   # present in STPSP
              tone$xpp$ISI>=minISI                   &   # present in STPSP
              tone$xpp$ISI<=maxISI                   &   # present in STPSP
              tone$xpp$trialType%in%trialTypeIndex   &   # present in STPSP
              tone$xpp[[DVstr]]>=valRng[1]                 &   # present in STPSP
              tone$xpp[[DVstr]]<=valRng[2]                 &   # present in STPSP
              tone$xpp$powEOG<=maxPowEOG)                # present in STPSP
  
  return(ybin)
}