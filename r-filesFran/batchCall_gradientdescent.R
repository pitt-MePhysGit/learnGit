
## libraries for parfor loop ====
library(foreach)
library(doParallel)

# Define Number of parallel clusters:
nClus <- 29


## COPIED FROM playground_gradientDescent.R

## Add packages -----------------------
library(pracma)
library(lattice)
library(ggplot2)
#library(corrplot)
rda <- function(){ return('/Volumes/EEGCyno/RES/ampOddClick/') }
source('~/Dropbox/r-compile/rUtils/R/rUtils.R')
source('~/Dropbox/r-compile/ePhysLab/R/ePhysLab.R')
source('~/Dropbox/ampOddClick/r-files/functions/functions_ampOdd_presynapticPlasticity.R')
source('~/Dropbox/ampOddClick/r-filesFran/load_ampOddClick_Human.R')
source('~/Dropbox/ampOddClick/r-filesFran/function_parloop_gradientdescent_TT.R')
# source('~/Dropbox/ampOddClick/r-filesFran/function_parloop_gradientdescent_TT_ampFUN2q.R')

Subject_list <- c('2193','2199','2235','2238','2250','2320','2259','2262','2274',
                  '2275','H001','2214','2212','2285','2337','2372','2389','H002',
                  '2299','2406','2426','2494','2427','2456','H005',
                  '2005','2355','2387','2288') 


## RUN TEST CALL OF PARLOOP FUNCTION: ====
# disc <- parloop_gradientdescent( 1 )


## CREATE PARALLEL CLUSTER ====
parallelCluster <- parallel::makeCluster( nClus, outfile = "") 
doParallel::registerDoParallel(parallelCluster)
print(parallelCluster)

## RUN PARALLEL LOOP ====
tryCatch(
  disc <- parallel::parLapply(parallelCluster, 1:length(Subject_list), parloop_gradientdescent), error = function(e) print(e)
)

## SHUT DOWN CLUSTER NEATLY ===
parallel::stopCluster(parallelCluster)
parallelCluster <- c()
