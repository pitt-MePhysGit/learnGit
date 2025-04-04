library(Rcpp)
library(RcppArmadillo)

# import the cpp functions from the cpp folder
sourceCpp("~/Dropbox/ampOddClick/r-files/cpp/recurrent.cpp")
sourceCpp("~/Dropbox/meanFieldSim//cpp/meanFieldSimWorkhorse.cpp")


#rm(list('ALL','Q','RE','RI'))
rm(ALL); rm(RE); rm(RI);rm(Q)

taxis <- seq(0,to=600,by=.0001)
R0    <- 0.0 + as.numeric( (taxis>1 & taxis < 1.005) | (taxis>2 & taxis < 2.00500) )
plot(taxis,R0, type='l')

mode <- 'depressionStabilized'
mode <- 'inhibitionStabilized'

if (mode == 'depressionStabilized'){
  U    <- .05
  lmbd <- 6000
  
  Je0 <- .05
  Ji0 <- 0
  
  Jee <- 2.5
  Jie <- 0
  Jei <- 0
  Jii <- .1
  
  te <- .050
  ti <- .050
}

if (mode == 'inhibitionStabilized'){
  U    <- .1
  lmbd <- 6000
  
  Je0 <- 0.15
  Ji0 <- 0.0
  
  Jee <- .5
  Jie <- .0
  Jei <- .1
  Jii <- .0
  
  te <- .05
  ti <- .05
}

system.time(ALL <- recurrent( R0, lmbd, U, te, ti, Je0, Ji0, Jee, Jie, Jei, Jii))
Q  <- ALL[seq(1,by=1,length=length(taxis))]
RE <- ALL[seq(1+length(taxis),by=1,length=length(taxis))]
RI <- ALL[seq(1+2*length(taxis),by=1,length=length(taxis))]

plot(  taxis, RE, type='l', ylim=c(0,1), xlim=c(.95, 1.5))
points(taxis, R0, type='l', col='green')
points(taxis, RI, type='l', col='red')
points(taxis,  Q, type='l', col='blue')

plot(  RE, Q, type='l' )



## new multi-node function
sourceCpp("~/Dropbox/meanFieldSim/cpp/generalized_meanFieldSim.cpp")





# try to run the workhorse function
# constants, p,  R0
ALL <- s_results( c(lmbd, U, te, ti), c(Je0, Ji0, Jee, Jie, Jei, Jii), R0)




