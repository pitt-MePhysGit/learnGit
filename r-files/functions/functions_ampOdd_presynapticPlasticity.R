## use Rcpp to define a fast routine to calculate the presynaptic plasticity model
library(Rcpp)
library(RcppArmadillo)
library(MASS)

# import the cpp functions from the cpp folder
sourceCpp("~/Dropbox/ampOddClick/r-files/cpp/cumsum1.cpp")
sourceCpp("~/Dropbox/ampOddClick/r-files/cpp/cumsumamp.cpp")
sourceCpp("~/Dropbox/ampOddClick/r-files/cpp/cumsumamp2t.cpp")
##I believe the following function is useless
sourceCpp("~/Dropbox/ampOddClick/r-files/cpp/cumsumamp2q.cpp")
sourceCpp("~/Dropbox/ampOddClick/r-files/cpp/cumsumscl.cpp")


## ============================ work flow
## 1) call_STPSP      : runs the gradient descent
## 2) run_Model_STPSP : calculates predictions for one specific set of parameters, lots of plotting (calls either 3 or 4)
## 3) modelFunctions  : bundle many channels into one global response
## 4) cpp functions   : work-horse building blocks to get responses of individual channels ('cumsum1', 'cumsumamp', etc)


## model functions: intermittent between C++ routines like 'cumsum1' and 'cumsumamp', and runModelSTPSP
## is replacing the direct call to cumsum1 in the old runModel STPSP
## overall function call: STPSPmodelFun_XXX(x,prm)
## x corresponds to tone$xpp, prm corresponds to the set of parameters
## XXX can be any specific function that implements any desired functionality

###10202021 SK added ability to give any subset of parameters to ampFUN (manStartPar is full set)
##
##20211119 SK changed the inputs to lm in many AMPFUNs to take data = y, and var names in formula
##THIS means that using results extracted from STPSPres$res$lmi requires passing DIFFERENT var names



STPSPmodelFun_Amp <- function(x, prm, thrsh=NULL, linear=T, returnStr='df'){
  # based on STPSPmodelFun_freq
  
  # x contains a list of times and amplitudes of tones
  # typically tone amplitude is coded by integers between 1 and 5
  # they represent five tones with 6dB between them
  
  # below is a list of parameters for the fitting routine
  # par[1]: time constant; 
  # par[2]: expended fraction; 
  # par[3]: Off 
  
  
  parList <- list()
  parList$startPar <- c(2,     0.5,        1 )#,      6,   1,        1    )
  parList$lower    <- c(1e-10, 0.01,       0.1)#,    0,   .1,      .001  )
  parList$upper    <- c(10,    1,          10)#,     12,   2,       10    )
  parList$ndeps    <- c(1e-2 , 1e-2,    1e-3)#,    1e-2 , 1e-3,    1e-3  )
  parList$parscale <- c(1e-2,  1e-2,    1e-2)#,    1e-1,  1e-2,    1e-2  )
  
  if (returnStr=='parList' )
    return(parList)
  
  
  if (!returnStr=='parList'){
    
    # by default assume that there are as many channels as there are different tones
    if (is.null(thrsh))
      thrsh <- sort( unique(x$ampInd) )
    
    mxAmp     <- max(x$ampInd)
    Nchan     <- length(thrsh)
    x$ampNorm <- (x$ampInd+prm[3])/(mxAmp+prm[3])
    
    #browser()
    getChanCumsumamp <- function(n){
      # n corresponds to channel number; channel has activation threshold thrsh[n]
      
      relAmp  <- (x$ampNorm - thrsh[n])%>=%0
      if (!linear)
        relAmp <- as.numeric( relAmp > 0)
      
      rsp     <- cumsumamp( x$ISI,  relAmp, prm[1], prm[2]  )
      return(rsp)
    }
    
    #resList       <- lapply(seq(from=-Nchan+1,to=2*Nchan,by=1) ,   getChanCumsumamp)
    resList        <- lapply(seq(from=1,to=Nchan,by=1) ,   getChanCumsumamp)
    resMat         <- array(unlist(resList), c(length(resList[[1]]),length(resList))  )
    
    if (returnStr=='sum' ){
      totalResponse <- apply(resMat,1,na.sum)
      return(totalResponse)
    }
    
    if (returnStr=='df'){
      for (m in 1:dim(resMat)[2])
        resMat[,m] <- resMat[,m]#-mean(resMat[,m])
      
      resDF          <- data.frame(resMat)
      names(resDF)   <- bestFrq 
      return(resDF)
    }
  }
}

ampFUN_U0 <- function(parm, tone, ybin, DVstr='n80', BYstr='exp',shuffel='no', lmOffset=F, returnStr='model', show=FALSE, eps=F,chNr=1,    
                      fileName='defaultFileName'){
  
  # param        parameter vector. at this point, 5 parameters
  #    parm[1]              time contant of STPSP model
  #    parm[2]              release probability of STPSP model
  #    parm[3]              mean non-linearity
  #    parm[4]              sd non-linearity
  #    parm[5]              offset
  
  # DVstr:       which variable to model
  # returnStr   return full model results or not
  # show
  # eps
  # chNr
  
  #par[1]:time constant; par[2]: release probability; par[3]: mu nonlinearlity; par[4]: sd non-linearity par[5] minimium drive
  parList          <- list()
  parList$name     <- 'ampFUN_U0'
  parList$startPar <- c(2,     0.6,   2.5,    1,  0)
  parList$lower    <- c(1e-10, 0.01, -10,  1e-1,  0)
  parList$upper    <- c(10,    1,     10,    20, 10)
  parList$ndeps    <- c(1e-3 , 1e-3, 1e-3, 1e-3, 1e-3)
  parList$parscale <- c(1e-2,  1e-2, 1e-2, 1e-2, 1e-2)
  
  if (returnStr=='parList' ){
    return(parList)
  }
  
  if (!returnStr=='parList'){
    
    lmbd      <- parm[1]
    U         <- .0001
    mu        <- parm[3] # x-shift
    sd        <- parm[4] # width
    off       <- parm[5] # offset D
    
    memnl     <- F   # non-linear memory 
    
    NFrq      <- 1
    deltaFrq  <- NA
    x         <- tone$xpp
    
    x$ISI[ which(x$ISI<0) ] <- 20
    
    ## D stands for drive 
    #x$D     <- (x$ampInd+Off)/(5+Off)
    #umin    <- min(x$D) + range(x$D)*Fu
    #x$Dthrs <- x$D%>=%umin
    x$D      <- pnorm(x$ampInd, m=mu, s=sd) + off
    x$D      <- x$D/max(x$D)
    
    if(sum(is.na(x$D))>0)
      browser()
    
    x$R   <- cumsumscl( x$ISI, x$D, lmbd, U )
    x$R   <- x$R/U 
    x$BY  <- factor(x[[BYstr]])
    
    if(sum(is.na(x$R))>0)
      browser()
    
    yind        <- which(ybin)
    y           <- subset(x, ybin, drop=T)
    
    print(dim(y)[1]/dim(x)[1])
    
    if (length(table(y$BY))>1){
      lm0 <- lm( y[[DVstr]] ~ y$BY       )
      #lm1 <- lm( y[[DVstr]] ~ y$BY + y$R)
      lm1 <- lm( y[[DVstr]] ~ y$BY + y$R:y$BY )
    }
    
    if (length(table(y$BY))==1){
      if (lmOffset){
        lm0     <- lm( y[[DVstr]] ~ 1      )
        lmdummy <- lm( y[[DVstr]] ~ 1 + y$ampFct + y$logISI) 
        lm1     <- lm( y[[DVstr]] ~ 1 + y$R)
      }
      if (!lmOffset){
        lm0 <- lm( y[[DVstr]] ~ -1   )
        lm1 <- lm( y[[DVstr]] ~ -1 + y$R)
      }
      #lm2 <- lm1
    }
    lm2 <- lm1
    print(lm1$coef[ length(lm1$coef) ])
    
    modelTest <- anova(lm0,lm1)
    print( c( parm,  diff(modelTest$RSS), modelTest$Pr[2]) )
    
    if (show)
      showSTPSPFit( list(y=y,DVstr=DVstr,lm0=lm0, lm1=lm1, lm2=lm2, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest, fileName=fileName, eps=eps ) )
    
    if(returnStr=='RSS')
      return( diff(modelTest$RSS) ) 
    
    if(returnStr=='Model')
      return( list(lm0=lm0, lm1=lm1, lm2=lm2, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest ) )#, mnMat=mnMat, mnDV=mnDV, binKer=binKer,
    
  }
}

ampFUN_U1 <- function(parm, tone, ybin, DVstr='n80', BYstr='exp',shuffel='no', lmOffset=F, returnStr='model', show=FALSE, eps=F,chNr=1,    
                      fileName='defaultFileName'){
  
  # param        parameter vector. at this point, 5 parameters
  #    parm[1]              time contant of STPSP model
  #    parm[2]              release probability of STPSP model
  #    parm[3]              mean non-linearity
  #    parm[4]              sd non-linearity
  #    parm[5]              offset
  
  # DVstr:       which variable to model
  # returnStr   return full model results or not
  # show
  # eps
  # chNr
  
  #par[1]:time constant; par[2]: release probability; par[3]: mu nonlinearlity; par[4]: sd non-linearity par[5] minimium drive
  parList          <- list()
  parList$name     <- 'ampFUN_U1'
  parList$startPar <- c(2,     0.6,   2.5,    1,  0)
  parList$lower    <- c(1e-10, 0.01, -10,  1e-1,  0)
  parList$upper    <- c(10,    1,     10,    20, 10)
  parList$ndeps    <- c(1e-3 , 1e-3, 1e-3, 1e-3, 1e-3)
  parList$parscale <- c(1e-2,  1e-2, 1e-2, 1e-2, 1e-2)
  
  if (returnStr=='parList' ){
    return(parList)
  }
  
  if (!returnStr=='parList'){
    
    lmbd      <- parm[1]
    U         <- 1
    mu        <- parm[3] # x-shift
    sd        <- parm[4] # width
    off       <- parm[5] # offset D
    
    memnl     <- F   # non-linear memory 
    
    NFrq      <- 1
    deltaFrq  <- NA
    x         <- tone$xpp
    
    x$ISI[ which(x$ISI<0) ] <- 20
    
    ## D stands for drive 
    #x$D     <- (x$ampInd+Off)/(5+Off)
    #umin    <- min(x$D) + range(x$D)*Fu
    #x$Dthrs <- x$D%>=%umin
    x$D      <- pnorm(x$ampInd, m=mu, s=sd) + off
    x$D      <- x$D/max(x$D)
    
    if(sum(is.na(x$D))>0)
      browser()
    
    x$R   <- cumsumscl( x$ISI, x$D, lmbd, U )
    x$R   <- x$R/U 
    x$BY  <- factor(x[[BYstr]])
    
    if(sum(is.na(x$R))>0)
      browser()
    
    yind        <- which(ybin)
    y           <- subset(x, ybin, drop=T)
    
    print(dim(y)[1]/dim(x)[1])
    
    if (length(table(y$BY))>1){
      lm0 <- lm( y[[DVstr]] ~ y$BY       )
      #lm1 <- lm( y[[DVstr]] ~ y$BY + y$R)
      lm1 <- lm( y[[DVstr]] ~ y$BY + y$R:y$BY )
    }
    
    if (length(table(y$BY))==1){
      if (lmOffset){
        lm0 <- lm( y[[DVstr]] ~ 1   )
        lm1 <- lm( y[[DVstr]] ~ 1 + y$R)
      }
      if (!lmOffset){
        lm0 <- lm( y[[DVstr]] ~ -1   )
        lm1 <- lm( y[[DVstr]] ~ -1 + y$R)
      }
      #lm2 <- lm1
    }
    lm2 <- lm1
    print(lm1$coef[ length(lm1$coef) ])
    
    modelTest <- anova(lm0,lm1)
    print( c( parm,  diff(modelTest$RSS), modelTest$Pr[2]) )
    
    if (show)
      showSTPSPFit( list(y=y,DVstr=DVstr,lm0=lm0, lm1=lm1, lm2=lm2, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest, fileName=fileName, eps=eps ) )
    
    if(returnStr=='RSS')
      return( diff(modelTest$RSS) ) 
    
    if(returnStr=='Model')
      return( list(lm0=lm0, lm1=lm1, lm2=lm2, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest ) )#, mnMat=mnMat, mnDV=mnDV, binKer=binKer,
    
  }
}


ampFUNscl <- function(parm, tone, ybin, DVstr='n80', BYstr='exp',shuffel='no', lmOffset=F, returnStr='model', show=FALSE, eps=F,chNr=1,    
                   fileName='defaultFileName'){
  
  # param        parameter vector. at this point, 5 parameters
  #    parm[1]              time contant of STPSP model
  #    parm[2]              release probability of STPSP model
  #    parm[3]              mean non-linearity
  #    parm[4]              sd non-linearity
  #    parm[5]              offset
  
  # DVstr:       which variable to model
  # returnStr   return full model results or not
  # show
  # eps
  # chNr
  
  #par[1]:time constant; par[2]: release probability; par[3]: mu nonlinearlity; par[4]: sd non-linearity par[5] minimium drive
  parList          <- list()
  parList$name     <- 'ampFUNscl'
  parList$startPar <- c(2,     0.6,   2.5,    1,  0)
  parList$lower    <- c(1e-10, 0.01, -10,  1e-1,  0)
  parList$upper    <- c(10,    1,     10,    20, 10)
  parList$ndeps    <- c(1e-3 , 1e-3, 1e-3, 1e-3, 1e-3)
  parList$parscale <- c(1e-2,  1e-2, 1e-2, 1e-2, 1e-2)
  
  if (returnStr=='parList' ){
    return(parList)
  }
  
  if (!returnStr=='parList'){
    
    lmbd      <- parm[1]
    U         <- parm[2]
    mu        <- parm[3] # x-shift
    sd        <- parm[4] # width
    off       <- parm[5] # offset D
    
    memnl     <- F   # non-linear memory 
    
    NFrq      <- 1
    deltaFrq  <- NA
    x         <- tone$xpp
    
    x$ISI[ which(x$ISI<0) ] <- 20
    
    ## D stands for drive 
    #x$D     <- (x$ampInd+Off)/(5+Off)
    #umin    <- min(x$D) + range(x$D)*Fu
    #x$Dthrs <- x$D%>=%umin
    x$D      <- pnorm(x$ampInd, m=mu, s=sd) + off
    x$D      <- x$D/max(x$D)
    
    if(sum(is.na(x$D))>0)
      browser()
    
    x$R   <- cumsumscl( x$ISI, x$D, lmbd, U )
    x$R   <- x$R/U 
    x$BY  <- factor(x[[BYstr]])
    
    if(sum(is.na(x$R))>0)
      browser()
    
    yind        <- which(ybin)
    y           <- subset(x, ybin, drop=T)
    
    print(dim(y)[1]/dim(x)[1])
    
    if (length(table(y$BY))>1){
      lm0 <- lm( y[[DVstr]] ~ y$BY       )
      #lm1 <- lm( y[[DVstr]] ~ y$BY + y$R)
      lm1 <- lm( y[[DVstr]] ~ y$BY + y$R:y$BY )
    }
    
    if (length(table(y$BY))==1){
      if (lmOffset){
        lm0 <- lm( y[[DVstr]] ~ 1   )
        lm1 <- lm( y[[DVstr]] ~ 1 + y$R)
      }
      if (!lmOffset){
        lm0 <- lm( y[[DVstr]] ~ -1   )
        lm1 <- lm( y[[DVstr]] ~ -1 + y$R)
      }
      #lm2 <- lm1
    }
    lm2 <- lm1
    print(lm1$coef[ length(lm1$coef) ])
    
    modelTest <- anova(lm0,lm1)
    print( c( parm,  diff(modelTest$RSS), modelTest$Pr[2]) )
    
    if (show)
      showSTPSPFit( list(y=y,DVstr=DVstr,lm0=lm0, lm1=lm1, lm2=lm2, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest, fileName=fileName, eps=eps ) )
    
    if(returnStr=='RSS')
      return( diff(modelTest$RSS) ) 
    
    if(returnStr=='Model')
      return( list(lm0=lm0, lm1=lm1, lm2=lm2, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest ) )#, mnMat=mnMat, mnDV=mnDV, binKer=binKer,
    
  }
}







mylmf = function(...) {
  args = list(...)
  
#  browser()
  
  if(names(args)[2]!='data'){
    stop('Please specify data SECOND (after model formula) in the top level model calls.')
  }
  
  isfam<-0
  if(any(names(args)=='family')){
    isfam<-1
    if(length(args)>2){
      isgauss<-unlist(lapply(lapply(args[-1], function(x){all.equal(deparse(x), deparse(gaussian))}), isTRUE))
      islin <- any(isgauss)
      if(islin){
        ilin  <- 1+which(isgauss)
      }
    }else{
      islin <- isTRUE(all.equal(deparse(args[[2]]), deparse(gaussian)))
      if(islin){
        ilin<-2
      }
    }
  }else{
    islin <- 1
  }
  
  isrob <- isTRUE(all.equal(deparse(lmf), deparse(rlm)))
  if(!isrob & islin){
    lmf<-lm
  }
  
  ##CAREFUL is specifying both robust and family
  ##default behavior is for robust to override non-Gaussian glm family
  
  if(isrob & isfam){
    if(length(args)==2){
      y = do.call(lmf, list(args[-ilin]))
    }else{
      y = do.call(lmf, args[-ilin])
    }
  }
  else{
    y = do.call(lmf, args)
  }
  
#  browser()
  
  return(y)
}





ampFUN <- function(parm, tone, ybin, DVstr='n80', BYstr='exp', lmOffset=T, shuffel='no', returnStr='model', show=FALSE, eps=F, chNr=1,    
                   fileName='defaultFileName', robust = FALSE, family = gaussian, testOn = FALSE, parList = NULL){
  
  # param        parameter vector. at this point, 5 parameters
  #    parm[1]              time contant of STPSP model
  #    parm[2]              release probability of STPSP model
  #    parm[3]              mean non-linearity
  #    parm[4]              sd non-linearity
  #    parm[5]              offset
  
  # DVstr:       which variable to model
  # returnStr   return full model results or not
  # show
  # eps
  # chNr
  
#  browser()
  
  #par[1]:time constant; par[2]: release probability; par[3]: mu nonlinearity; par[4]: sd non-linearity par[5] minimum drive
  if(is.null(parList)){
    parList          <- list()
  }
  if(is.null(parList$name)){
    parList$name     <- 'ampFUN'
  }
  
  if(is.null(parList$lower)){
    parList$lower    <- c(1e-10, 0.01, -5,   1e-1, 0.01)
  }
  if(is.null(parList$upper)){
    parList$upper    <- c(10,    0.99,  10,  20,   10)
  }
  if(is.null(parList$startPar)){
#    parList$startPar <- c(8,     0.6,   2.5,    1, 0)
    parList$startPar <- (parList$lower+parList$upper)/2
  }
  if(is.null(parList$ndeps)){
    parList$ndeps    <- c(1e-3 , 1e-3, 1e-3, 1e-3, 1e-3)
  }
  if(is.null(parList$parscale)){
    parList$parscale <- c(1e-2,  1e-2, 1e-2, 1e-2, 1e-2)
  }
  
  if (returnStr=='parList' ){
    return(parList)
  }
  
  if (!returnStr=='parList'){
    
    lmbd      <- parm[1]
    U         <- parm[2]
    mu        <- parm[3] # x-shift
    sd        <- parm[4] # width
    off       <- parm[5] # offset D
    
    memnl     <- F   # non-linear memory 
    
    NFrq      <- 1
    deltaFrq  <- NA
    x         <- tone$xpp
    
    x$ISI[ which(x$ISI<0) ] <- 20
    
    ## D stands for drive 
    #x$D     <- (x$ampInd+Off)/(5+Off)
    #umin    <- min(x$D) + range(x$D)*Fu
    #x$Dthrs <- x$D%>=%umin
    x$D      <- pnorm(x$ampInd, m=mu, s=sd) + off
    x$D      <- x$D/max(x$D)
    
    if(sum(is.na(x$D))>0)
      browser()
    
    x$R   <- cumsumamp( x$ISI, x$D, lmbd, U )
    x$R   <- x$R/U 
    x$BY  <- factor(x[[BYstr]])
    
    if(sum(is.na(x$R))>0)
      browser()
    
    
    yind        <- which(ybin)
    y           <- subset(x, ybin, drop=T)
    
    print(dim(y)[1]/dim(x)[1])

#    browser()
    
    ##should probably make lm default, with option for glm in mylmf
    ##setting global is obviously not the right way to do this, but it works for now
    lmf<<-glm
    if(robust){
      lmf<<-rlm
    }
    
    if (length(table(y$BY))>1){
      lm0 <- mylmf( y[[DVstr]] ~ BY       , data = y, family = family)
      lm1 <- mylmf( y[[DVstr]] ~ BY + R, data = y, family = family)
      lmx <- mylmf( y[[DVstr]] ~ BY + R:BY, data = y, family = family)
    }
    if (length(table(y$BY))==1){
      if (lmOffset){
        lm0 <- mylmf( y[[DVstr]] ~ 1, data = y, family = family)
        lm1 <- mylmf( y[[DVstr]] ~ 1 + R, data = y, family = family)
        lmx <- lm1
        ##note reverse
        if(testOn){
        #  lmp0 <- glm(y$R ~ 1, family = quasibinomial(link = "logit"))
        #  lmp <- glm(y$R ~ 1 + y[[DVstr]], family = quasibinomial(link = "logit"))
          
#          browser()
          lm0 <- mylm(y[[DVstr]] ~ 1, data = y, family = family)
          tx<-pmax(pmin(y$R, 0.99), 0.01)
          lm1 <- mylm(y[[DVstr]] ~ 1 + log(tx/(1-tx)), data = y, family = family)
        }
      }
      if (!lmOffset){
        lm0 <- mylmf( y[[DVstr]] ~ -1, data = y, family = family)
        lm1 <- mylmf( y[[DVstr]] ~ -1 + R, data = y, family = family)
        lmx <- mylmf( y[[DVstr]] ~  1 + R, data = y, family = family)
      }
    }
    lm2 <- lm1

    print(lm1$coef[ length(lm1$coef) ])
    
    modelTest <- anova(lm0, lm1)
#    if(testOn){
#      modelTestp <- anova(lmp, lmp0)
#    }
    
    ##NOTE: I'm replacing the call to RSS within the anova by the call to column 2, 
    #same for lm, but Resid. deviance under glm
    ##SK 20220425
    
    #print( c( parm,  diff(modelTest$RSS), modelTest$Pr[2]) )
    print( c( parm,  diff(modelTest[,2]), modelTest$Pr[2]) )
    
    if (show)
      showSTPSPFit( list(y=y,DVstr=DVstr,lm0=lm0, lm1=lm1, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest, fileName=fileName, eps=eps ) )
    
    if(returnStr=='RSS'){
#      if(testOn){
#        return(modelTestp$Dev[2])
#      }else{
#        return( diff(modelTest$RSS) ) 
      return( diff(modelTest[,2]) ) 
      #      }
    }
    
    if(returnStr=='Model'){
  #    if(testOn){
  #      return( list(lm0=lmp0, lm1=lmp, lm2=lmx, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest ) )#, mnMat=mnMat, mnDV=mnDV, binKer=binKer,
  #    }else{
        return( list(lm0=lm0, lm1=lm1, lm2=lmx, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest ) )#, mnMat=mnMat, mnDV=mnDV, binKer=binKer,
      }
  #  }

  }
}



ampFUN2t <- function(parm, tone, ybin, DVstr='n80', BYstr='exp', lmOffset=T, shuffel='no', returnStr='model', show=FALSE, eps=F, chNr=1,    
                     fileName='defaultFileName', robust = FALSE, family = gaussian, testOn = FALSE, parList = NULL){
  
  # param        parameter vector. at this point, 5 parameters
  #    parm[1]              1/slow time constant of STPSP model
  #    parm[2]              release probability/depleted fraction of STPSP model
  #    parm[3]              1/fast-1/slow time constant of STPSP model
  #    parm[4]              mean non-linearity
  #    parm[5]              sd non-linearity
  #    parm[6]              offset
  
  # DVstr:       which variable to model
  # returnStr   return full model results or not
  # show
  # eps
  # chNr
  
  #par[1]:time constant; par[2]: release probability; par[3]: mu nonlinearity; par[4]: sd non-linearity par[5] minimum drive
  if(is.null(parList)){
    parList          <- list()
  }
  if(is.null(parList$name)){
    parList$name     <- 'ampFUN'
  }
  
  if(is.null(parList$lower)){
    parList$lower    <- c(1e-10, 0.01, 1e-1, -5, 1e-1, 0.01)
  }
  if(is.null(parList$upper)){
    parList$upper    <- c(10,    0.99,   10, 10,   20,   10)
  }
  if(is.null(parList$startPar)){
    parList$startPar <- (parList$lower+parList$upper)/2
  }
  if(is.null(parList$ndeps)){
    parList$ndeps    <- c(1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3)
  }
  if(is.null(parList$parscale)){
    parList$parscale <- c(1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2)
  }
  
  if (returnStr=='parList' ){
    return(parList)
  }
  
  if (!returnStr=='parList'){
    
    lmbd1      <- parm[1]
    U         <- parm[2]
    lmbd2      <- parm[1]+parm[3]
    mu        <- parm[4] # x-shift
    sd        <- parm[5] # width
    off       <- parm[6] # offset D
    
    NFrq      <- 1
    deltaFrq  <- NA
    x         <- tone$xpp
    
    x$ISI[ which(x$ISI<0) ] <- 20
    
    ## D stands for drive 
    x$D      <- pnorm(x$ampInd, m=mu, s=sd) + off
    x$D      <- x$D/max(x$D)
    
    if(sum(is.na(x$D))>0)
      browser()
    
    x$R   <- cumsumamp2t(x$ISI, x$D, lmbd1, U, lmbd2)
    x$R   <- x$R/U 

    x$BY  <- factor(x[[BYstr]])
    
    if(sum(is.na(x$R))>0)
      browser()
    
    
    yind        <- which(ybin)
    y           <- subset(x, ybin, drop=T)
    
    print(dim(y)[1]/dim(x)[1])
    
    #    browser()
    
    lmf<-glm
    if(robust){
      lmf<-rlm
    }
    
    if (length(table(y$BY))>1){
      lm0 <- lmf( y[[DVstr]] ~ y$BY       , data = y, family = family)
      lm1 <- lmf( y[[DVstr]] ~ y$BY + y$R, data = y, family = family)
      lmx <- lmf( y[[DVstr]] ~ y$BY + y$R:y$BY , data = y, family = family)
    }
    if (length(table(y$BY))==1){
      if (lmOffset){
        lm0 <- lmf( y[[DVstr]] ~ 1   , data = y, family = family)
        lm1 <- lmf( y[[DVstr]] ~ 1 + R, data = y, family = family)
        lmx <- lm1
        ##note reverse
        if(testOn){
          lm0 <- lm(y[[DVstr]] ~ 1, data = y, family = family)
          tx<-pmax(pmin(y$R, 0.99), 0.01)
          lm1 <- lm(y[[DVstr]] ~ 1 + log(tx/(1-tx)), data = y, family = family)
        }
      }
      if (!lmOffset){
        lm0 <- lmf( y[[DVstr]] ~ -1   , data = y, family = family)
        lm1 <- lmf( y[[DVstr]] ~ -1 + y$R, data = y, family = family)
        lmx <- lmf( y[[DVstr]] ~  1 + y$R, data = y, family = family)
      }
    }
    lm2 <- lm1
    
    print(lm1$coef[ length(lm1$coef) ])
    
    modelTest <- anova(lm0, lm1)
    #    if(testOn){
    #      modelTestp <- anova(lmp, lmp0)
    #    }
    
    print( c( parm,  diff(modelTest$RSS), modelTest$Pr[2]) )
    
    if (show)
      showSTPSPFit( list(y=y,DVstr=DVstr,lm0=lm0, lm1=lm1, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest, fileName=fileName, eps=eps ) )
    
    if(returnStr=='RSS'){
      #      if(testOn){
      #        return(modelTestp$Dev[2])
      #      }else{
      return( diff(modelTest$RSS) ) 
      #      }
    }
    
    if(returnStr=='Model'){
      #    if(testOn){
      #      return( list(lm0=lmp0, lm1=lmp, lm2=lmx, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest ) )#, mnMat=mnMat, mnDV=mnDV, binKer=binKer,
      #    }else{
      return( list(lm0=lm0, lm1=lm1, lm2=lmx, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest ) )#, mnMat=mnMat, mnDV=mnDV, binKer=binKer,
    }
    #  }
    
  }
}



ampFUN2q <- function(parm, tone, ybin, DVstr='n80', BYstr='exp', lmOffset=T, shuffel='no', returnStr='model', show=FALSE, eps=F, chNr=1,    
                     fileName='defaultFileName', robust = FALSE, family = gaussian, testOn = FALSE, parList = NULL, topInt = 0, Uscale = TRUE){
  
  # param        parameter vector. at this point, 5 parameters
  #    parm[1]              1/slow time constant of STPSP model
  #    parm[2]              slow release probability/depleted fraction of STPSP model
  #    parm[3]              1/fast-1/slow time constant of STPSP model
  #    parm[4]              fast release probability/depleted fraction of STPSP model
  #    parm[5]              mean non-linearity
  #    parm[6]              sd non-linearity
  #    parm[7]              offset
  
  # DVstr:       which variable to model
  # returnStr   return full model results or not
  # show
  # eps
  # chNr
  
  #par[1]:time constant; par[2]: release probability; par[3]: mu nonlinearity; par[4]: sd non-linearity par[5] minimum drive
  if(is.null(parList)){
    parList          <- list()
  }
  if(is.null(parList$name)){
    parList$name     <- 'ampFUN'
  }
  
  if(is.null(parList$lower)){
    parList$lower    <- c(1e-10, 0.01, 1e-1,  0.01,  -5, 1e-1, 0.01)
  }
  if(is.null(parList$upper)){
    parList$upper    <- c(10,    0.99,   10,  0.99,  10,   20,   10)
  }
  if(is.null(parList$startPar)){
    parList$startPar <- (parList$lower+parList$upper)/2
  }
  if(is.null(parList$ndeps)){
    parList$ndeps    <- c(1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3)
  }
  if(is.null(parList$parscale)){
    parList$parscale <- c(1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2)
  }
  
  if (returnStr=='parList' ){
    return(parList)
  }
  
  if (!returnStr=='parList'){
    
    lmbd1      <- parm[1]
    U1         <- parm[2]
    lmbd2      <- parm[1]+parm[3]
    U2         <- parm[4]
    mu        <- parm[5] # x-shift
    sd        <- parm[6] # width
    off       <- parm[7] # offset D
    
    NFrq      <- 1
    deltaFrq  <- NA
    x         <- tone$xpp
    
    x$ISI[ which(x$ISI<0) ] <- 20
    
    ## D stands for drive 
    x$D      <- pnorm(x$ampInd, m=mu, s=sd) + off
    x$D      <- x$D/max(x$D)
      
    if(sum(is.na(x$D))>0)
      browser()
    
    
    x$R1   <- cumsumamp(x$ISI, x$D, lmbd1, U1)
    x$R2   <- cumsumamp(x$ISI, x$D, lmbd2, U2)
    if(Uscale){
      x$R1   <- x$R1/U1 
      x$R2   <- x$R2/U2 
    }
    
    x$BY  <- factor(x[[BYstr]])
    
    if(sum(is.na(x$R))>0)
      browser()
    
    yind        <- which(ybin)
    y           <- subset(x, ybin, drop=T)
    
    print(dim(y)[1]/dim(x)[1])
    
    #    browser()
    
    lmf<-glm
    if(robust){
      lmf<-rlm
    }
    
    if (length(table(y$BY))==1){
      lm0 <- lmf( y[[DVstr]] ~ 1   , data = y, family = family)
      if(topInt == 0){
        lm1 <- lmf( y[[DVstr]] ~ 1 + I(R1 + R2), data = y, family = family)
      }
      if(topInt == 1){
        lm1 <- lmf( y[[DVstr]] ~ 1 + R1 + R2, data = y, family = family)
      }
      if(topInt == 2){
        lm1 <- lmf( y[[DVstr]] ~ 1 + R1*R2, data = y, family = family)
      }
      if(topInt == 3){
        lm1 <- lmf( y[[DVstr]] ~ 1 + R1*R2 + I(R1*R2/D), data = y, family = family)
      }
      if(topInt == 4){
        lm1 <- lmf( y[[DVstr]] ~ 1 + R1 + R2 + pmin(R1, R2), data = y, family = family)
      }
      lmx <- lm1
      # Added to compare results of gradient descent with simple linear model (5/16/23)
      lmdummy <- lm( y[[DVstr]] ~ 1 + y$ampFct + y$logISI) 
    }

    print(lm1$coef[ length(lm1$coef) ])
    
    #browser()
    modelTest <- anova(lm0, lm1)
   
    print( c( parm,  diff(modelTest$RSS), modelTest$Pr[2]) )
    
    if (show)
      showSTPSPFit( list(y=y,DVstr=DVstr,lm0=lm0, lm1=lm1, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest, fileName=fileName, eps=eps ) )
    
    if(returnStr=='RSS'){
      #return( diff(modelTest$RSS) ) 
      return( modelTest$Dev[2] )
    }
    
    if(returnStr=='Model'){
      return( list(lm0=lm0, lm1=lm1, lm2=lmx, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest ) )#, mnMat=mnMat, mnDV=mnDV, binKer=binKer,
    }
  }
}



ampFUN3q <- function(parm, tone, ybin, DVstr='n80', BYstr='exp', lmOffset=T, shuffel='no', returnStr='model', show=FALSE, eps=F, chNr=1,    
                     fileName='defaultFileName', robust = FALSE, family = gaussian, testOn = FALSE, parList = NULL, topInt = 0, Uscale = TRUE){
  
  # param        parameter vector. at this point, 5 parameters
  #    parm[1]              1/slow time constant of STPSP model
  #    parm[2]              slow release probability/depleted fraction of STPSP model
  #    parm[3]              1/mid-1/slow time constant of STPSP model
  #    parm[4]              mid release probability/depleted fraction of STPSP model
  #    parm[5]              1/fast-1/mid-1/slow time constant of STPSP model
  #    parm[6]              fast release probability/depleted fraction of STPSP model
  #    parm[7]              mean non-linearity
  #    parm[8]              sd non-linearity
  #    parm[9]              offset
  
  if(is.null(parList)){
    parList          <- list()
  }
  if(is.null(parList$name)){
    parList$name     <- 'ampFUN'
  }
  
  if(is.null(parList$lower)){
    parList$lower    <- c(1e-10, 0.01, 1e-1,  0.01, 1e-1,  0.01,  -5, 1e-1, 0.01)
  }
  if(is.null(parList$upper)){
    parList$upper    <- c(10,    0.99,   10,  0.99, 10,  0.99, 10,   20,   10)
  }
  if(is.null(parList$startPar)){
    parList$startPar <- (parList$lower+parList$upper)/2
  }
  if(is.null(parList$ndeps)){
    parList$ndeps    <- c(1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3)
  }
  if(is.null(parList$parscale)){
    parList$parscale <- c(1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2)
  }
  
  if (returnStr=='parList' ){
    return(parList)
  }
  
  if (!returnStr=='parList'){
    
    lmbd1      <- parm[1]
    U1         <- parm[2]
    lmbd2      <- parm[1]+parm[3]
    U2         <- parm[4]
    lmbd3      <- parm[1]+parm[3]+parm[5]
    U3         <- parm[6]
    mu        <- parm[7] # x-shift
    sd        <- parm[8] # width
    off       <- parm[9] # offset D
    
    NFrq      <- 1
    deltaFrq  <- NA
    x         <- tone$xpp
    
    x$ISI[ which(x$ISI<0) ] <- 20
    
    ## D stands for drive 
    x$D      <- pnorm(x$ampInd, m=mu, s=sd) + off
    x$D      <- x$D/max(x$D)
    
    if(sum(is.na(x$D))>0)
      browser()
    
    
    x$R1   <- cumsumamp(x$ISI, x$D, lmbd1, U1)
    x$R2   <- cumsumamp(x$ISI, x$D, lmbd2, U2)
    x$R3   <- cumsumamp(x$ISI, x$D, lmbd3, U3)
    if(Uscale){
      x$R1   <- x$R1/U1 
      x$R2   <- x$R2/U2 
      x$R3   <- x$R3/U3 
    }
    
    x$BY  <- factor(x[[BYstr]])
    
    if(sum(is.na(x$R))>0)
      browser()
    
    yind        <- which(ybin)
    y           <- subset(x, ybin, drop=T)
    
    print(dim(y)[1]/dim(x)[1])
    
    #    browser()
    
    lmf<-glm
    if(robust){
      lmf<-rlm
    }
    
    #names(list(1~1, data = y, family = gaussian))[1]!='data'
    #names(args)[1]<-'formula'

    #lmf(args)

    #lmf(1~1, data = y, family = gaussian)
    #mylmf(1~1, data = y, family = gaussian)
    
    
    if (length(table(y$BY))==1){
      lm0 <- lmf( y[[DVstr]] ~ 1, data = y, family = family)
      if(topInt == 0){
        lm1 <- lmf( y[[DVstr]] ~ 1 + I(R1 + R2 + R3), data = y, family = family)
      }
      if(topInt == 1){
        lm1 <- lmf( y[[DVstr]] ~ 1 + R1 + R2 + R3, data = y, family = family)
      }
      if(topInt == 2){
        lm1 <- lmf( y[[DVstr]] ~ 1 + R1*R2*R3, data = y, family = family)
      }
      if(topInt == 3){
#        lm1 <- lmf( y[[DVstr]] ~ 1 + y$R1*y$R2 + I(y$R1*y$R2/y$D))
         lm1 <- lmf( y[[DVstr]] ~ 1 + (R1 + R2 + R3)^2, data = y, family = family)
      }
      if(topInt == 4){
        lm1 <- lmf( y[[DVstr]] ~ 1 + R1 + R2 + R3 + pmin(R1, R2, R3), data = y, family = family)
      }
      if(topInt == 5){
        lm1 <- lmf( y[[DVstr]] ~ 1 + R1 + R2 + R3 + R1:R3 + R2:R3, data = y, family = family)
      }
      lmx <- lm1
    }
    
    print(lm1$coef[ length(lm1$coef) ])
    
    modelTest <- anova(lm0, lm1)
    
    print( c( parm,  diff(modelTest$RSS), modelTest$Pr[2]) )
    
    if (show)
      showSTPSPFit( list(y=y,DVstr=DVstr,lm0=lm0, lm1=lm1, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest, fileName=fileName, eps=eps ) )
    
    if(returnStr=='RSS'){
      return( diff(modelTest$RSS) ) 
    }
    
    if(returnStr=='Model'){
      return( list(lm0=lm0, lm1=lm1, lm2=lmx, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest ) )#, mnMat=mnMat, mnDV=mnDV, binKer=binKer,
    }
  }
}


#no robust option
##do we want separate amplitude functions by resource? probably
##should make this take arbitrarily many response variables and resources
ampFUN2mv3q <- function(parm, tone, ybin, DVstr='active', BYstr='exp', lmOffset=T, shuffel='no', returnStr='model', show=FALSE, eps=F, chNr=1,    
                     fileName='defaultFileName', testOn = FALSE, parList = NULL, topInt = 1, Uscale = TRUE, robust = FALSE, family = gaussian){
  
  # param        parameter vector. at this point, 5 parameters
  #    parm[1]              1/slow time constant of STPSP model
  #    parm[2]              slow release probability/depleted fraction of STPSP model
  #    parm[3]              1/mid-1/slow time constant of STPSP model
  #    parm[4]              mid release probability/depleted fraction of STPSP model
  #    parm[5]              1/fast-1/mid-1/slow time constant of STPSP model
  #    parm[6]              fast release probability/depleted fraction of STPSP model
  #    parm[7]              mean non-linearity
  #    parm[8]              sd non-linearity
  #    parm[9]              offset
  
  if(is.null(parList)){
    parList          <- list()
  }
  if(is.null(parList$name)){
    parList$name     <- 'ampFUN'
  }
  
  if(is.null(parList$lower)){
    parList$lower    <- c(1e-10, 0.01, 1e-1,  0.01, 1e-1,  0.01,  -5, 1e-1, 0.01)
  }
  if(is.null(parList$upper)){
    parList$upper    <- c(10,    0.99,   10,  0.99, 10,  0.99, 10,   20,   10)
  }
  if(is.null(parList$startPar)){
    parList$startPar <- (parList$lower+parList$upper)/2
  }
  if(is.null(parList$ndeps)){
    parList$ndeps    <- c(1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3)
  }
  if(is.null(parList$parscale)){
    parList$parscale <- c(1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2)
  }
  
  if (returnStr=='parList' ){
    return(parList)
  }
  
  if (!returnStr=='parList'){
    
    lmbd1      <- parm[1]
    U1         <- parm[2]
    lmbd2      <- parm[1]+parm[3]
    U2         <- parm[4]
    lmbd3      <- parm[1]+parm[3]+parm[5]
    U3         <- parm[6]
    mu        <- parm[7] # x-shift
    sd        <- parm[8] # width
    off       <- parm[9] # offset D
    
    NFrq      <- 1
    deltaFrq  <- NA
    x         <- tone$xpp
    
    x$ISI[ which(x$ISI<0) ] <- 20
    
    ## D stands for drive 
    x$D      <- pnorm(x$ampInd, m=mu, s=sd) + off
    x$D      <- x$D/max(x$D)
    
    if(sum(is.na(x$D))>0)
      browser()
    
    
    x$R1   <- cumsumamp(x$ISI, x$D, lmbd1, U1)
    x$R2   <- cumsumamp(x$ISI, x$D, lmbd2, U2)
    x$R3   <- cumsumamp(x$ISI, x$D, lmbd3, U3)
    if(Uscale){
      x$R1   <- x$R1/U1 
      x$R2   <- x$R2/U2 
      x$R3   <- x$R3/U3 
    }
    
    x$BY  <- factor(x[[BYstr]])
    
    if(sum(is.na(x$R))>0)
      browser()
    
    yind        <- which(ybin)
    y           <- subset(x, ybin, drop=T)
    
    print(dim(y)[1]/dim(x)[1])
    
#        browser()
    
    lmf<-glm
    if(robust){
      lmf<-rlm
    }
    
    resp<-as.matrix(y[,grepl(DVstr, names(y))])
    if (length(table(y$BY))==1){
      lm0 <- lmf( resp ~ 1  , data = y )
      if(topInt == 0){
        lm1 <- lmf( resp ~ 1 + I(R1 + R2 + R3), data = y)
      }
      if(topInt == 1){
        lm1 <- lmf( resp ~ 1 + R1 + R2 + R3, data = y)
      }
      if(topInt == 2){
        lm1 <- lmf( resp ~ 1 + R1*R2*R3, data = y)
      }
      if(topInt == 3){
        lm1 <- lmf( resp ~ 1 + (R1 + R2 + R3)^2, data = y)
      }
      if(topInt == 4){
        lm1 <- lmf( resp ~ 1 + R1 + R2 + R3 + pmin(R1, R2, R3), data = y)
      }
      if(topInt == 5){
        lm1 <- lmf( resp ~ 1 + R1 + R2 + R3 + R1:R3 + R2:R3, data = y)
      }
      lmx <- lm1
    }
    
    print(lm1$coef[ length(lm1$coef) ])
    
    modelTest <- anova(lm0, lm1)
    
    print( c( parm,  diff(modelTest$RSS), modelTest$Pr[2]) )
    
    if (show)
      showSTPSPFit( list(y=y,DVstr=DVstr,lm0=lm0, lm1=lm1, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest, fileName=fileName, eps=eps ) )
    
    #browser()
    
    if(returnStr=='RSS'){
      return(modelTest$Gen.var.[2]) 
    }
    
    if(returnStr=='Model'){
      return( list(lm0=lm0, lm1=lm1, lm2=lmx, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest ) )#, mnMat=mnMat, mnDV=mnDV, binKer=binKer,
    }
  }
}


ampFUNnp <- function(parm, tone, ybin, DVstr='n80', BYstr='exp',lmOffset=T,shuffel='no', returnStr='model', show=FALSE, eps=F,chNr=1,    
                   fileName='defaultFileName', robust = FALSE, family = gaussian){
  
  # param        parameter vector. at this point, 5 parameters
  #    parm[1]              time contant of STPSP model
  #    parm[2]              release probability of STPSP model
  #    parm[3]              mean non-linearity
  #    parm[4]              sd non-linearity
  #    parm[5]              offset
  
  # DVstr:       which variable to model
  # returnStr   return full model results or not
  # show
  # eps
  # chNr
  
  #par[1]:time constant; par[2]: release probability; par[i] (i>2): unconstrained effect of amplitude i-2
  parList          <- list()
  parList$name     <- 'ampFUN'
  parList$startPar <- c(8,     0.6,  0.5,  0.5,  0.5,  0.5,  0.5)
  parList$lower    <- c(1e-10, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
  parList$upper    <- c(10,    1,    0.99, 0.99, 0.99, 0.99, 0.99)
  parList$ndeps    <- c(1e-3 , 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3)
  parList$parscale <- c(1e-2,  1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2)
  
  if (returnStr=='parList' ){
    return(parList)
  }
  
  if (!returnStr=='parList'){
    
    lmbd      <- parm[1]
    U         <- parm[2]
    a1        <- parm[3]
    a2        <- parm[4]
    a3        <- parm[5]
    a4        <- parm[6]
    a5        <- parm[7]

    NFrq      <- 1
    deltaFrq  <- NA
    x         <- tone$xpp
    
    x$ISI[ which(x$ISI<0) ] <- 20
    x$ampInd[ which(x$ampInd>5) ] <- 5
    
    x$D      <- c(a1, a2, a3, a4, a5)[x$ampInd]

    if(sum(is.na(x$D))>0)
      browser()
    
    x$R   <- cumsumamp( x$ISI, x$D, lmbd, U )
    x$R   <- x$R/U 
    x$BY  <- factor(x[[BYstr]])
    
    if(sum(is.na(x$R))>0)
      browser()
    
    
    yind        <- which(ybin)
    y           <- subset(x, ybin, drop=T)
    
    print(dim(y)[1]/dim(x)[1])
    
    lmf<-glm
    if(robust){
      lmf<-rlm
    }
    
    if (length(table(y$BY))>1){
      lm0 <- lmf( y[[DVstr]] ~ y$BY     , data = y  )
      #lm1 <- lm( y[[DVstr]] ~ y$BY + y$R)
      lm <- lmf( y[[DVstr]] ~ y$BY + y$R:y$BY, data = y )
    }
    
    if (length(table(y$BY))==1){
      if (lmOffset){
        lm0 <- lmf( y[[DVstr]] ~ 1  , data = y )
        lm1 <- lmf( y[[DVstr]] ~ 1 + y$R, data = y)
        lmx <- lm1
      }
      if (!lmOffset){
        lm0 <- lmf( y[[DVstr]] ~ -1   , data = y)
        lm1 <- lmf( y[[DVstr]] ~ -1 + y$R, data = y)
        lmx <- lmf( y[[DVstr]] ~  1 + y$R, data = y)
      }
    }
    lm2 <- lm1
    
    print(lm1$coef[ length(lm1$coef) ])
    
    modelTest <- anova(lm0,lm1)
    print( c( parm,  diff(modelTest$RSS), modelTest$Pr[2]) )
    
    if (show)
      showSTPSPFit( list(y=y,DVstr=DVstr,lm0=lm0, lm1=lm1, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest, fileName=fileName, eps=eps ) )
    
    if(returnStr=='RSS')
      return( diff(modelTest$RSS) ) 
    
    if(returnStr=='Model')
      return( list(lm0=lm0, lm1=lm1, lm2=lmx, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest ) )#, mnMat=mnMat, mnDV=mnDV, binKer=binKer,
    
  }
}


#left as WIP
##trying to add an new fully parameterized ampFUN that captures the diversity of modulation we tend to see
ampFUNfp <- function(parm, tone, ybin, DVstr='n80', BYstr='exp',lmOffset=T,shuffel='no', returnStr='model', show=FALSE, eps=F,chNr=1,    
                     fileName='defaultFileName', robust = FALSE, family = gaussian, t0 = 0.5, t1 = 5.5){
  
  # param        parameter vector. at this point, 5 parameters
  #    parm[1]              time contant of STPSP model
  #    parm[2]              release probability of STPSP model
  #    parm[3]              mean non-linearity
  #    parm[4]              sd non-linearity
  #    parm[5]              offset
  
  # DVstr:       which variable to model
  # returnStr   return full model results or not
  # show
  # eps
  # chNr
  
  #par[1]:time constant; par[2]: release probability; par[i] (i>2): unconstrained effect of amplitude i-2
  parList          <- list()
  parList$name     <- 'ampFUN'
  # parList$startPar <- c(8,     0.6,  0.5,     0.5,    0.5) # Monkey param
  parList$startPar <- c(2,     0.1,  0.5,     0.5,    0.5) # Human param
  parList$lower    <- c(1e-10, 0.01, 0.01,    0.01,  0.01)
  parList$upper    <- c(10,    1,    0.99,    0.99,   0.99)
  parList$ndeps    <- c(1e-3 , 1e-3, 1e-3, 1e-3, 1e-3)
  parList$parscale <- c(1e-2,  1e-2, 1e-2, 1e-2, 1e-2)
  
  if (returnStr=='parList' ){
    return(parList)
  }
  
  if (!returnStr=='parList'){
    
    lmbd      <- parm[1]
    U         <- parm[2]
    
    ##
    ##thresh (in A) prior 0
    ##thresh (in A) post 1
    ##^ these 2 go outside in a profile computation for each combination
    ##a+be^cx
    #####rejected the linear middle
    ##lin beta_0 at 0a
    ##lin beta_1
    
    a1        <- parm[3]
    a2        <- parm[4]
    a3        <- parm[5]

    NFrq      <- 1
    deltaFrq  <- NA
    x         <- tone$xpp
    
    x$ISI[ which(x$ISI<0) ] <- 20
    x$ampInd[ which(x$ampInd>5) ] <- 5

    theta1<-a1
    theta2<-(a2*(1-theta1))+theta1
    theta3<-(a3*(1-theta2+theta1))-theta1
    
    rsig <- 4/(qnorm(theta2)-qnorm(theta1))
    rmu  <- 1-rsig*qnorm(theta1)   
    
    x$D      <- pnorm(x$ampInd, rmu, rsig) + theta3
    x$D[x$ampInd<t0]<-0
    x$D[x$ampInd>t1]<-1
#    x$D[x$D<0]<-0
#    x$D[x$D>0]<-1
#    x$D      <- x$D/max(x$D)

        
    if(sum(is.na(x$D))>0)
      browser()
    
    x$R   <- cumsumamp( x$ISI, x$D, lmbd, U )
    x$R   <- x$R/U 
    x$BY  <- factor(x[[BYstr]])
    
    if(sum(is.na(x$R))>0)
      browser()
    
    
    yind        <- which(ybin)
    y           <- subset(x, ybin, drop=T)
    
    print(dim(y)[1]/dim(x)[1])
    
    lmf<-glm
    if(robust){
      lmf<-rlm
    }
    
    if (length(table(y$BY))>1){
      lm0 <- lmf( y[[DVstr]] ~ y$BY       , data = y)
      #lm1 <- lm( y[[DVstr]] ~ y$BY + y$R)
      lm <- lmf( y[[DVstr]] ~ y$BY + y$R:y$BY , data = y)
    }
    
    if (length(table(y$BY))==1){
      if (lmOffset){
        lm0 <- lmf( y[[DVstr]] ~ 1   , data = y)
        lm1 <- lmf( y[[DVstr]] ~ 1 + y$R, data = y)
        lmx <- lm1
      }
      if (!lmOffset){
        lm0 <- lmf( y[[DVstr]] ~ -1   , data = y)
        lm1 <- lmf( y[[DVstr]] ~ -1 + y$R, data = y)
        lmx <- lmf( y[[DVstr]] ~  1 + y$R, data = y)
      }
    }
    lm2 <- lm1
    
    print(lm1$coef[ length(lm1$coef) ])
    
    modelTest <- anova(lm0,lm1)
    print( c( parm,  diff(modelTest$RSS), modelTest$Pr[2]) )
    
    if (show)
      showSTPSPFit( list(y=y,DVstr=DVstr,lm0=lm0, lm1=lm1, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest, fileName=fileName, eps=eps ) )
    
    if(returnStr=='RSS')
      return( diff(modelTest$RSS) ) 
    
    if(returnStr=='Model')
      return( list(lm0=lm0, lm1=lm1, lm2=lmx, Nx=dim(x),yind=yind,ybin=ybin, Ny=dim(y), modelTest=modelTest ) )#, mnMat=mnMat, mnDV=mnDV, binKer=binKer,
    
  }
}



##Inserting capacity to include whatever desired model here
mylmf = function(...) {
  args = list(...)
  
  #browser()
  
  if(names(args)[2]!='data'){
    stop('Please specify data SECOND (after model formula) in the top level model calls.')
  }
  
  isfam<-0
  if(any(names(args)=='family')){
    isfam<-1
    if(length(args)>2){
      isgauss<-unlist(lapply(lapply(args[-1], function(x){all.equal(deparse(x), deparse(gaussian))}), isTRUE))
      islin <- any(isgauss)
      if(islin){
        ilin  <- 1+which(isgauss)
      }
    }else{
      islin <- isTRUE(all.equal(deparse(args[[2]]), deparse(gaussian)))
      if(islin){
        ilin<-2
      }
    }
  }else{
    islin <- 1
  }
  isrob <- isTRUE(all.equal(deparse(lmf), deparse(rlm)))
  
  #      if(any(lapply(args, function(x){x[[1]]}))){
  #        args$a = bar     # change argument "a"
  #      }
  
  ##CAREFUL is specifying both robust and family
  ##default behavior is for robust to override non-Gaussian glm family
  
  #browser()
  
  if(isrob){
    if(length(args)==2){
      y = do.call(rlm, list(args[-ilin]))
    }else{
      y = do.call(rlm, args[-ilin])
    }
  }
  else{
    y = do.call(glm, args)
  }
  return(y)
}



##added the option to manually pass in starting parameter values; shouldn't affect anything previously written
callSTPSP <- function(tn, DVstr='mmn', gd=T, show=TRUE, showCh=1, fileExtn='', startFromScratch=FALSE,doType=FALSE,who='all',FUN=ampFUN,baseName='',
                      maxISI=25, minISI=0.2,minFrq=500,toneBy=1/2,NOct=3,maxp2p=450,DVsdscale=5, maxEOG=Inf, trialType=NA, fitType=NA, manStartPar = NULL, 
                      robust = FALSE, family = gaussian, parList0 = NULL, verbose = TRUE, path = FALSE,...){
  if(!gd & path){
    path<-FALSE
    warning("Haven't implemented path for genoud yet; can't return the path.")
  }
  
  if(!identical(FUN, ampFUNfp)){
    parList <- FUN( NULL , NA, NA, DVstr=DVstr, returnStr='parList', parList = parList0)
  }else{
    parList <- FUN( NULL , NA, NA, DVstr=DVstr, returnStr='parList')
  }
  
  FUNext <- parList$name
  
  if(!is.null(manStartPar)){
    parList$startPar<-manStartPar
  }
  
  #FUNext = 'cumsumamp'
  #if(identical(FUN,cumsumscl) )
  #  FUNext <- 'cumsumscl'
  #if(identical(FUN,cumsum1) )
  #  FUNext <- 'cumsum1'
  

  fitTypeExtn <- ifelse(gd,'gd','gen')
  fileName    <- paste( baseName, '/STPSP/STPSPFit_',tn$exp,'_',who,'_',DVstr,'_',FUNext, ifelse(nchar(fileExtn)==0,'', '_'), fileExtn,'_',fitTypeExtn,'.rda',sep='' )
  print(fileName)
  
  #fileName    <- paste( baseName,'/STPSP/STPSPFit_',tn$exp,'_',who,'_',DVstr,            ifelse(nchar(fileExtn)==0,'', '_'), fileExtn,'_',fitTypeExtn,'.rda',sep='' )
  print( paste(rda(),fileName,sep=''))
  
  trialType <- c(1,2)
  if(identical(who,'random'))
    trialType <- c(1)
  if(identical(who,'regular'))
    trialType <- c(2)
  
  if(is.na(fitType))
    fitType <- c(0,1)

  animal <- strsplit(tn$exp,'_')[[1]][1]
  
  if(animal=='Rockey') # & what %in% c('raw','mlp','mbp','bsp'))
    maxEOG <- 15000
  
  if(animal=='Walter')# & what %in% c('raw','mlp','mbp','bsp'))
    maxEOG <- 5000
  
  if(animal=='Jesse')# & what %in% c('raw','mlp','mbp','bsp'))
    maxEOG <- 2000
  
  if(animal=='Sam')# & what %in% c('raw','mlp','mbp','bsp'))
    maxEOG <- 15000
  
  frq       <- 1000
  NFrq      <- length(frq)
  deltaFrq  <- NA
  NdeltaFrq <- 0
  
  ## ============================================= do the subsetting before calling the fitting function
  if(is.character(tn$xpp$exp))
    tn$xpp$exp <- factor(tn$xpp$exp)
  
  tn$xpp$firstTrialInDay <- diff( c(0,as.numeric(tn$xpp$exp)) )
  
  x       <- tn$xpp
  
  mDV       <- quantile( x[[DVstr]], prob=c(0.25,0.5,0.75) )
  valRng    <- mDV[2] + c(-DVsdscale,DVsdscale) * diff(mDV) 
  ##just to make spiking work for now
  if(mean(x[[DVstr]]==0)>0.5){
    valRng    <- range(x[[DVstr]])
  }
  propAccpt <- length( which( x[[DVstr]]>valRng[1] &  x[[DVstr]]<valRng[2] ) )/dim(x)[1]
  print(mDV)
  print(valRng)
  print(propAccpt)
  
  if (!exists('p2p',x))
    x$p2p <- maxp2p + c(1,-1)[ 1 + as.numeric( x$valTrial ) ]
  
  ybin <- ( x$p2p<maxp2p                  & 
              x$ISI>minISI                & 
              x$ISI<maxISI                & 
              x$powEOG<maxEOG             &
              x$firstTrialInDay==0        & 
              x[[DVstr]]>=valRng[1]       &  
              x[[DVstr]]<=valRng[2]       & 
              x$toneType%in%fitType       & 
              x$trialType%in%trialType)
  
  print(sum(ybin)/dim(x)[1] )
  
#  browser()
  
  # maxISI:    for fit, exclude trials with an ISI above this value
  # minISI:    for fit, exclude trials with an ISI below this value
  # maxp2p:    for fit, exclude trials with peak-to-peak value above this value
  # DVsdscale: for fit, exclude trials with a dependent variable plus minus DVsdscale time the interquartile distance
  # maxEOG   : for fit, exclude trials with powEOG above this value
  # fitType  : for fit, include only toneType that are in fitType
  # trialType: for fit, include only trialType that are in trialType
  
  ## ========================================================
  doFitNow <- FALSE
  if( !file.exists(paste(rda(),fileName,sep='')) | startFromScratch )
    doFitNow <- TRUE
  
  if (!doFitNow)  
    load( paste(rda(),fileName,sep='') )
  
  parpath<-c()
  
  if(!identical(FUN, ampFUNfp)){
    ## ================================== load start parameters and fitting parameters
    parList <- FUN( NULL , tn, ybin, DVstr=DVstr, returnStr='parList', parList = parList0)
    mnmx    <- array(c(parList$lower, parList$upper) ,c(2,length(parList$upper)))

    ## prepare the output of the fit
    thisModel   <- function(p){FUN( p , tn, ybin, DVstr=DVstr, shuffel='no',returnStr='RSS', robust = robust, family = family, parList = parList,...)}
    
    #browser()
    
    ## run the fit if necessary
    if (doFitNow){
      paste( 'running Fit ...' )
      
      thisOptim <- function(tmdl){
        cl <- makeForkCluster(6)
        res <- genoud(tmdl, length(parList$startPar), max=FALSE, pop.size=1000, max.generations=100, wait.generations=5,Domains=mnmx,
                      hard.generation.limit=TRUE, starting.values=parList$startPar, MemoryMatrix=TRUE, 
                      solution.tolerance=0.0001,gr=NULL, boundary.enforcement=2, lexical=FALSE, gradient.check=TRUE, BFGS=TRUE,
                      BFGSburnin=99,control=list(ndeps=parList$ndeps,parscale=parList$parscale), cluster=cl)
        stopCluster(cl=cl)
        return(res)
      }
    
      print(parList$startPar)
      if(gd & verbose){
        modelFit <- optim( parList$startPar, thisModel, lower=parList$lower, upper=parList$upper, method='L-BFGS-B', control=list(ndeps=parList$ndeps, parscale=parList$parscale))
      }else if(gd & !verbose){
        ##if debugging within a capture.output use:
        ##closeAllConnections()
        invisible(capture.output(modelFit <- optim( parList$startPar, thisModel, lower=parList$lower, upper=parList$upper, method='L-BFGS-B', control=list(ndeps=parList$ndeps,parscale=parList$parscale))))
      }
      parpath<-c()
      if(gd & path){
        set.seed(123)
        modelFit<-optim(parList$startPar, thisModel, lower = parList$lower, upper = parList$upper, method = 'L-BFGS-B', control = list(ndeps = parList$ndeps, parscale = parList$parscale))
        parpath<-matrix(0, length(parList$startPar), modelFit$counts['function'])
        print(modelFit$counts['function'])
        for(i in 1:modelFit$counts['function']){
          set.seed(123)
          print(c(i, (i*(i+1))/(modelFit$counts['function']*(modelFit$counts['function']+1))))
          flush.console()
          invisible(capture.output(parpath[,i]<-optim(parList$startPar, thisModel, lower = parList$lower, upper = parList$upper, method='L-BFGS-B', control = list(ndeps = parList$ndeps, parscale = parList$parscale, maxit = i))$par))
        }
        parpath<-cbind(parList$startPar, parpath)
      }
      if(!gd)
        modelFit <- thisOptim( thisModel  )
      
      modelFit$thisPar <- c(modelFit$par)
    }
  
    #resMMN   <- runModelSTPSP( modelFit$thisPar, tn, ybin, DVstr=DVstr,
    #                           returnModel=TRUE, show=show, chNr=showCh,fileName=fileName,FUN=FUN,...)
    resMMN   <- FUN( modelFit$par, tn, rep(T,dim(tn$xpp)[1]), DVstr=DVstr,
                     returnStr='Model', show=FALSE, chNr=showCh, 
                     fileName=fileName, robust = robust, family = family,...)
    
  }else{
    parpath<-c()
    parList <- FUN( NULL , tn, ybin, DVstr=DVstr, returnStr='parList')

    paste( 'running Fit ...' )
    print(parList$startPar)
    
    tmp1<-cbind(rep(1:(max(tn$xpp$ampInd)+1), (max(tn$xpp$ampInd))+1), rep(1:(max(tn$xpp$ampInd)+1), each = (max(tn$xpp$ampInd))+1))
    tmp2<-tmp1[tmp1[,1]<tmp1[,2],]-0.5
    
    mFl<-list()
    mFp<-c()
    for(i in 1:dim(tmp2)[1]){
      thisModel   <- function(p){FUN(p, tn, ybin, DVstr = DVstr, shuffel = 'no', returnStr = 'RSS', robust = robust, family = family, t0 = tmp2[i,1], t1 = tmp2[i,2], ...)}

##for ampFUNfp, there are no genetic, path or verbose options available
##
##assuming doFitNow
      invisible(capture.output(modelFit <- optim( parList$startPar, thisModel, lower=parList$lower, upper=parList$upper,method='L-BFGS-B' , control=list(ndeps=parList$ndeps,parscale=parList$parscale, maxit = 2))))
      
      resMMN   <- FUN( modelFit$par, tn, rep(T,dim(tn$xpp)[1]), DVstr=DVstr,
                       returnStr='Model', show=FALSE, chNr=showCh, 
                       fileName=fileName, robust = robust, family = family, t0 = tmp2[i, 1], t1 = tmp2[i, 2], ...)

      mFl[[i]]<-modelFit
      mFp[i]<-summary(resMMN$lm1)$sigma
      print(c("Thresholds:", tmp2[i,], "RSE", mFp[i]))      
    }
    
    print(which.min(mFp))
    modelFit<-mFl[[which.min(mFp)]]
    modelFit$thisPar <- c(modelFit$par, tmp2[which.min(mFp),])
    
    resMMN   <- FUN( modelFit$par, tn, rep(T,dim(tn$xpp)[1]), DVstr=DVstr,
                     returnStr='Model', show=FALSE, chNr=showCh, 
                     fileName=fileName, robust = robust, family = family, t0 = tmp2[which.min(mFp), 1], t1 = tmp2[which.min(mFp), 2], ...)
    
  }

  modelFit$modelTest  <- resMMN$modelTest
  modelFit$propAccept <- resMMN$propAccept
  modelFit$coeff      <- c(resMMN$lm1$coeff, summary(resMMN$lm1)$sigma, summary(resMMN$lm1)$adj.r.squared)
  names(modelFit$coeff)[(length(modelFit$coeff)-1):length(modelFit$coeff)]<-c("RSE", "AdjRsquared")
  modelFit$tmpMean    <- resMMN$tmpMean

  fileName    <- paste(         '/STPSP/STPSPFit_',tn$exp,'_',who,'_',DVstr,'_',FUNext, ifelse(nchar(fileExtn)==0,'', '_'), fileExtn,'_',fitTypeExtn,'.rda',sep='' )
  print(  paste(rda(),fileName,sep=''))

  
  
  #browser()
  
  save( file=paste(rda(),fileName,sep=''), modelFit)
  
  resMMN   <- FUN( modelFit$par, tn, rep(T,dim(tn$xpp)[1]), DVstr=DVstr,
                   returnStr='Model', show=show, chNr=showCh, 
                   fileName=fileName, robust = robust, family = family,...)

  return( list( modelFit=modelFit, res=resMMN, ppath = parpath) )
}





# if ( returnModel | show ){
#   
#   ## get some numbers
#   ## changed from y$tpKer to y$oKer
#   Nbin   <- min( c(20, round(dim(y)[1]/500 ))  )%>=%6
#   binX   <- quantile(y$R, prob=seq(1/Nbin,by=1/Nbin,length=Nbin-1) )
#   
#   binKer <- rep(1, dim(y)[1] )
#   for (bnd in 1:length(binX)) {
#     print(bnd)
#     binKer <- binKer + (y$R>binX[bnd])
#   }
#   
#   binKer <- binX[binKer]
#   mnDV   <- tapply(lm0$resid,binKer,mean)
#   mnISI  <- tapply(y$ISI,binKer,mean)
#   mnMat  <- NULL
# }
# 
# # --------------------
# if ( (returnModel | show) & !is.null(tone$lfp)){
#   if(length(dim(tone$lfp))==3 )
#     sdat <- tone$lfp[yind,,chNr]
#   
#   if(length(dim(tone$lfp))==2 )
#     sdat <- tone$lfp[yind,]
#   
#   mnLst <- lapply( sort(unique(binKer)) , function(x){apply( sdat[which(binKer==x),],2,mean) } )
#   mnMat <- array(unlist(mnLst), c(length(mnLst[[1]]),length(mnLst)) );
# }
# 




## make some figures
showSTPSPFit <- function( fitRes ){
  
  #attach(fitRes$)
  lwd <- 2
  
  ## bin according to standard log time
  #brks                <- exp( log(2) * seq(log(0.2)/log(2),by=1/3,to=3.7) )
  #brks                <- c(brks, 25)
  
  #logrng <- round(log(range(y$ISI))/log(2))
  #brks                <- exp( log(2) * seq(logrng[1],by=1/2,to=logrng[2]) )
  
  #binKerISI           <- cut( y$ISI, breaks=brks )
  
  # averages according to log time
  #mnDVbyISI     <- tapply(y[[DVstr]],binKerISI,mean)
  #mnISIbyISI    <- tapply(y$ISI,binKerISI,mean)
  #mnPredbyISI1  <- tapply(lm1$fitted.values,binKerISI,mean)
  
  ## ======================================================================
  ## this is the pretty plot that mirrors the no-model plot
  picFile <- paste(fitRes$fileName,'.eps',sep='')
  if(fitRes$eps)
    eps.file(picFile,pdf=T,ratio=1/.66, s=2 )
  
  #brks <- exp( log(2) * seq(log(0.2)/log(2),by=1/3,to=4.1) )
  #brks[length(brks)] <- 20

  #y$binISI <- cut( y$ISI, breaks=brks )
  fitRes$y$DV <- fitRes$y[[fitRes$DVstr]]
  
  #xax             <- tapply(y$ISI, y$binISI, mean)
  #N1              <- tapply(y$DV, y$binISI, mean)
  #N1p             <- tapply(lm1$fitted.values,y$binISI,mean)

  xax             <- tapply(fitRes$y$ISI, fitRes$y$isiOct, mean)
  N1              <- tapply(fitRes$y$DV             ,fitRes$y$isiOct, mean)
  N1p             <- tapply(fitRes$lm1$fitted.values,fitRes$y$isiOct,mean)
  
  logxax          <- log(xax)/log(2)
  names(N1)       <- xax
  
  ylim    <- na.range(c(N1))+c(-1,1)*diff(na.range(c(N1)))*.1
  par(mfcol=c(1,1))
  plot(logxax, N1 , type='b', pch=20, ylim=ylim,ylab=paste(fitRes$DVstr,' amplitude [uV]',sep=''),xlab='ISI [sec]', xaxt='n',lwd=2)
  addSE <- function(nd){sebar( fitRes$y$DV[nd], at=log(mean(fitRes$y$ISI[nd]))/log(2)) }
  
  disc <- tapply( 1:dim(fitRes$y)[1], fitRes$y$isiOct, addSE ) 
  axis(1, at=log( 0.2*2^seq(0,by=1,to=6))/log(2), label=0.2*2^seq(0,by=1,to=6) )
  
  # add model fit as blue line
  points(logxax, N1p, type='l', col='blue',lwd=2)
  
  if(fitRes$eps)
    dev.off()
  
  
  
  mn     <- tapply(fitRes$y$DV, list(fitRes$y$dB,fitRes$y$isiOct), na.mean )
  mnse   <- tapply(fitRes$y$DV, list(fitRes$y$dB,fitRes$y$isiOct), na.se   )
  dbax   <- tapply(fitRes$y$dB,      fitRes$y$dB, mean              )
  isiax  <- tapply(fitRes$y$ISI,            fitRes$y$isiOct, mean   )
  
  mnp   <- tapply(fitRes$lm1$fitted.values, list(fitRes$y$dB,fitRes$y$isiOct), na.mean )
  
  N <- table(fitRes$y$isiOct)
  vi <- which(N>(length(fitRes$y$isiOct)/20))
  
  ylim <- na.range(c(mn,mnp)) + c(-.1,.25)*diff(na.range(c(mn,mnp)))
  
  picFile <- paste(fitRes$fileName,'_summary.eps',sep='')
  if(fitRes$eps)
    eps.file(picFile,pdf=T, s=2, mf=c(4,1),ratio=3/2 )
  
  par(mfrow=c(1,4))
  col <- bluered.colors(5)
  errorplot(log2(isiax[vi]),mn[1,vi], ey=mnse[1,vi], col='black', ylim=ylim,type='p',pch=20)
  points(   log2(isiax[vi]),mnp[1,vi], type='l',lwd=2)
  for (ax in 1:length(dbax)){
    errorplot(log2(isiax[vi]),mn[ax,vi], ey=mnse[ax,vi], col=col[ax], add=T,type='p',pch=20)
    points(   log2(isiax[vi]),mnp[ax,vi], type='l',lwd=2, col=col[ax])
  }
  
  col <- soa.colors(6)
  errorplot(dbax,mn[,vi[1]], ey=mnse[,vi[1]], col='black', ylim=ylim,type='p',pch=20)
  points(   dbax,mnp[,vi[1]], type='l',lwd=2)
  for (sx in 1:length(vi)){
    errorplot(dbax,mn[ ,vi[sx]], ey=mnse[,vi[sx]], col=col[sx], add=T,type='p',pch=20)
    points(   dbax,mnp[,vi[sx]], type='l',lwd=2, col=col[sx])
  }
  
  ## deltaAmp split by ampInd
  N1        <- tapply(fitRes$y$DV,  list(fitRes$y$dB,fitRes$y$deltaAmpInd), mean )
  N1se      <- tapply(fitRes$y$DV,  list(fitRes$y$dB,fitRes$y$deltaAmpInd), na.se )
  N1p       <- tapply(fitRes$lm1$fitted.values, list(fitRes$y$dB,fitRes$y$deltaAmpInd), mean )
  
  deltax    <- tapply(fitRes$y$deltaAmpInd, fitRes$y$deltaAmpInd, mean)
  col <- bluered.colors(5)
  for (ix in 1:(dim(N1)[1])){
    errorplot( deltax, N1[ix,], ey=N1se[ix,],  type='p', lwd=2,pch=20, col=col[ix],ylim=ylim, 
               add=ifelse(ix==1,F,T), xlab='delta ISI [log time]' )
    points(   deltax, N1p[ix,], type='l',lwd=2, col=col[ix])
  }
  
  ## deltaAmp split by isiOct
  N1        <- tapply(fitRes$y$DV,              list(fitRes$y$isiOct,fitRes$y$deltaAmpInd,fitRes$y$ampInd), mean )[,,3]
  N1se      <- tapply(fitRes$y$DV,              list(fitRes$y$isiOct,fitRes$y$deltaAmpInd,fitRes$y$ampInd), na.se )[,,3]
  N1p       <- tapply(fitRes$lm1$fitted.values, list(fitRes$y$isiOct,fitRes$y$deltaAmpInd,fitRes$y$ampInd), mean )[,,3]
  col       <- soa.colors(6)
  
  for (ix in vi ){
    errorplot( deltax, N1[ix,], ey=N1se[ix,],  type='p', lwd=2,pch=20, col=col[ix],ylim=ylim, 
               add=ifelse(ix==1,F,T), xlab='delta Amp [dB]' )
    points(   deltax, N1p[ix,], type='l',lwd=2, col=col[ix])
  }
  
  if(fitRes$eps)
    dev.off()
  
  
  if(1==0){  
    picFile <- paste(fitRes$fileName,'_TD_delta.eps',sep='')
    if(fitRes$eps)
      eps.file(picFile,pdf=T,ratio=1/.66, s=2 )
    showTD_delta(1,  pM=fitRes$y$DV, fM=fitRes$lm1$fitted.values, xpp=fitRes$y)
    if(fitRes$eps)
      dev.off()
    
    #detach(fitRes)
    
    picFile <- paste(fileName,'_evaluate.eps',sep='')
    if(eps)
      eps.file(picFile,pdf=T)
    
    par(mfrow=c(2,2))
    taxis <- seq(-30,-.200,by=.01)
    plot(    taxis, exp(parm[1]*taxis)/exp(parm[1]*taxis[length(taxis)]), type='l', xlab='time [sec]', ylab='' ,xlim=c(-15,0) )
    
    plot(binX, mnDV, type='l', lwd=2, ylim=range(mnDV)+c(-2,2) )
    addSE <- function(n){ 
      ind <- which(binKer == binX[n])
      sebar(lm0$resid[ind], at=binX[n],width=.5*min(diff(binX)))
    }
    disc <- lapply( 1:length(binX), addSE  )
    
    
    ## difference by log-time
    lmLog          <- lm( y[[DVstr]] ~ y$BY + log(y$ISI) )
    mnPredbyISILog <- tapply(lmLog$fitted.values,binKerISI,mean)
    
    print('Fraction of STPSP to logarithmic model RSS, full model')
    print( 100 * sum(lm1$residuals^2)/sum(lmLog$residuals^2) )
    
    ## by log-time
    plot(  mnISIbyISI, mnDVbyISI, type='b', lwd=2, ylim=na.range(c(mnDVbyISI))*c(-1.1,1.1), log='x',pch=20 )
    points(mnISIbyISI, mnPredbyISI1, type='l',lwd=2,col='red')
    points(mnISIbyISI, mnPredbyISILog, type='l',lwd=2,col='blue')
    
    plot(  mnISIbyISI, mnDVbyISI-mnPredbyISI1, type='l', lwd=2, ylim=na.range(c(mnDVbyISI-mnPredbyISI1))*c(-1.1,1.1), log='x',col='red' )
    points(mnISIbyISI, mnDVbyISI-mnPredbyISILog, type='l',lwd=2, col='blue')
    abline(h=0) 
    
    if(eps)
      dev.off()
  }
  
  if(  1==0 ){  
    picFile <- paste(fileName,'_ERP.eps',sep='')
    
    if(eps)
      eps.file( picFile, pdf=T )
    
    par(mfrow=c(1,1))
    
    cR <- colorRampPalette(c('black','blue','green','orange','red'), bias = 1, space ="Lab",interpolate="spline")
    palette(cR(Nbin))
    xlim=c(-150,600)
    ylim=c(-45,45)
    ylim=c(-100,200)
    
    
    if(tone$what=='mlp'){
      xlim=c(-10,70)
      ylim=c(-15,15)    
    }
    
    if(DVstr=='mua'){
      xlim=c(-50,250)
      ylim=c(-0.1,1)    
    }
    
    
    for(n in 1:(Nbin-1)){
      valInd <- which( binKer == binX[n]  )
      disc   <- mean.plot( sdat[valInd,], tone$taxis, add=ifelse(n==1,F,T), sd=TRUE, traces=F, 
                           posFUN=na.mean, widthFUN=na.se, range=NULL, col=n,lwd=lwd, ylim=ylim,
                           xlim=xlim)
    }
    abline(h=0)
    if(eps)
      dev.off()
  }
  
  ## ===================================================
  #par(mfcol=c(1,2))
  
    
  #image(dbax, deltax, N1, col=bluered.colors(100), zlim=zlim)
  #image(dbax, deltax, N1p, col=bluered.colors(100), zlim=zlim)
  
  ## ================================================== 
  if (1==0){
    isiax           <- log(tapply(y$ISI, y$isiOct, mean))/log(2)
    ampax           <- tapply(y$dB, y$ampInd, mean)
    iv <- which(!is.na(isiax) )
    N1              <- tapply(y$DV             , list(y$isiOct,y$ampInd), mean)
    N1p             <- tapply(lm1$fitted.values, list(y$isiOct,y$ampInd), mean)
    #N1[  which(is.na(N1))] <- 0
    #N1p[which(is.na(N1p))] <- 0
    
    
    zlim <- na.range(c(N1,N1p))
    par(mfcol=c(1,3))
    image(isiax[iv],ampax,N1[iv,], xlab='SOA',ylab='SPL',zlim=zlim, col=heat.colors(100))
    image(isiax[iv],ampax,N1p[iv,], xlab='SOA',ylab='SPL',zlim=zlim, col=heat.colors(100))
    
    image(isiax[iv],ampax, N1[iv,]-N1p[iv,], zlim=c(-1.25,1.25)*na.max(abs(N1[iv,]-N1p[iv,])),col=bluered.colors(100))
    #na.max( abs(tst$amp-tst$scl)/max(c(tst$amp,tst$scl)) )
    
    if (length(table(y$exp))>1){
      lmX     <- lm(y$DV ~ as.factor(y$exp) + as.factor(y$ampInd) + I( log(y$ISI)) )
      n1nrm   <- lmX$residuals
      
      lmY     <- lm( lm1$fitted.values ~ as.factor(y$exp) + as.factor(y$ampInd) + I( log(y$ISI)) )
      n1nrmp   <- lmY$residuals
    }
    
    if (length(table(y$exp))==1){
      lmX     <- lm(y$DV ~  as.factor(y$ampInd) + I( log(y$ISI)) )
      n1nrm   <- lmX$residuals
      
      lmY     <- lm( lm1$fitted.values ~  as.factor(y$ampInd) + I( log(y$ISI)) )
      n1nrmp   <- lmY$residuals
    }
  }
}







## just for debuggin reasonse:
## check if function returns identical lambdas when release probability is fixed to one
callSTPSPfixRelease <- function(tn, DVstr='mmn', gd=T, show=TRUE, showCh=1, fileExtn='', startFromScratch=FALSE,doType=FALSE,who='all',
                                maxISI=25, minISI=0.2,maxp2p=450,DVsdscale=4, maxEOG=Inf, trialType=NA, fitType=NA,...){
  
  
  brks               <- log(500 * 2^seq(0,by=1/2,to=3))/log(2)
  brks[1]            <- log(499)/log(2)
  brks[length(brks)] <- log(4001)/log(2)
  tn$xpp$pitchBin    <- cut(tone$xpp$pitchOct, breaks=brks)
  
  fitTypeExtn <- ifelse(gd,'gd','gen')
  fileName    <- paste( 'STPSP/STPSPFit_Rprob1_',tn$exp,'_',who,'_',DVstr, ifelse(nchar(fileExtn)==0,'', '_'), fileExtn,'_',fitTypeExtn,'.rda',sep='' )
  print( paste(rda(),fileName,sep=''))
  
  trialType <- c(1,2)
  if(identical(who,'random'))
    trialType <- c(1)
  if(identical(who,'regular'))
    trialType <- c(2)
  
  if(is.na(fitType))
    fitType <- c(1,2)
  
  if(animal=='Rockey') # & what %in% c('raw','mlp','mbp','bsp'))
    maxPowEOG <- 5000
  
  if(animal=='Walter')# & what %in% c('raw','mlp','mbp','bsp'))
    maxPowEOG <- 5000
  
  if(animal=='Jesse')# & what %in% c('raw','mlp','mbp','bsp'))
    maxPowEOG <- 2000
  
  if(animal=='Sam')# & what %in% c('raw','mlp','mbp','bsp'))
    maxPowEOG <- 15000
  
  
  frq       <- 1000
  NFrq      <- length(frq)
  deltaFrq  <- NA
  NdeltaFrq <- 0
  
  ## do the subsetting before calling the fitting function
  x       <- tn$xpp
  
  mDV       <- quantile( x[[DVstr]], prob=c(0.25,0.5,0.75) )
  valRng    <- mDV[2] + c(-DVsdscale,DVsdscale) * diff(mDV) 
  propAccpt <- length( which( x[[DVstr]]>valRng[1] &  x[[DVstr]]<valRng[2] ) )/dim(x)[1]
  print(propAccpt)
  
  ybin <- ( x$p2p<maxp2p & x$ISI>minISI & x$ISI<maxISI & x$powEOG<maxEOG &  # at least half the number of tones in kernel as total kernal length, 
              x[[DVstr]]>valRng[1]  &  x[[DVstr]]<valRng[2] & x$toneType%in%fitType & x$trialType%in%trialType)
  
  doFitNow <- FALSE
  if( !file.exists(paste(rda(),fileName,sep='')) | startFromScratch )
    doFitNow <- TRUE
  
  if (!doFitNow)  
    load( paste(rda(),fileName,sep='') )
  
  ## prepare the output of the fit
  thisModel   <- function(p){runModelSTPSP( c(p,1) , tn,ybin, DVstr=DVstr,shuffel='no',...)}
  
  ## run the fit if necessary
  if (doFitNow){
    paste( 'running Fit ...' )
    
    #par[1]:time constant; par[2]: release probability fixed to 1
    startPar <- c(2 )
    lower    <- c(1e-10 )
    upper    <- c(10 )
    #ndeps    <- c(1e-3 , 1e-3 )
    #parscale <- c(1e-2,  1e-2 )
    
    mnmx <- array(c(lower, upper) ,c(2,2))
    
    print(startPar)
    if(gd)
      modelFit <- optimize( f=thisModel, lower=lower, upper=upper) #, method='L-BFGS-B' , control=list(ndeps=ndeps,parscale=parscale ))
    
    if(!gd)
      modelFit <- thisOptim( thisModel  )
    
    modelFit$par <- c(modelFit$minimum,1)
    modelFit$thisPar <- c(modelFit$par)
  }
  
  resMMN   <- runModelSTPSP( modelFit$thisPar, tn, ybin, DVstr=DVstr,
                             returnModel=TRUE, show=show, eps=F, chNr=showCh,fileName=fileName,...)
  
  modelFit$modelTest <- resMMN$modelTest
  modelFit$propAccept <- resMMN$propAccept
  
  save( file=paste(rda(),fileName,sep=''), modelFit)
  
  return( list( modelFit=modelFit, res=resMMN) )
}


