decodeSVM <- function(stx, subs=NULL,fileName='svm.eps',tfrac=.5,...){
  
  print(sum(subs))
  
  stx <- subset(stx,subs,drop=T) 
  stx <- droplevels(stx)
  
  print(dim(stx))
  
  
  tset        <- sample( 1:dim(stx)[1], dim(stx)[1]*tfrac )
  trainingSet <- stx[tset,]
  
  tmpset  <- 1:dim(stx)[1]
  eset    <- setdiff(tmpset,tset)
  evalSet <- stx[eset,]
  
  table(stx$category)
  Ncat <- length(levels(stx$category))
  
  svmfit = svm(category ~ ., data = trainingSet, kernel = "radial", cost = 100, scale = T, cross=10)#, gamma=.00001)
  #svmfit = svm(category ~ ., data = trainingSet, kernel = "linear")#, cost = 100, scale = T, cross=10, gamma=.00001)
  
  print(svmfit)
  trainingSet$prediction <- predict(svmfit, trainingSet)
  trainingSet$correct    <- trainingSet$category == trainingSet$prediction
  mean(trainingSet$correct)
  
  evalSet$prediction <- predict(svmfit,evalSet)
  evalSet$correct    <- evalSet$category == evalSet$prediction
  mean(evalSet$correct)
  
  tbl          <- table(category=evalSet$category, prediction=evalSet$prediction)
  getCnfCat    <- function(cx){ table(category=evalSet$category, prediction=evalSet$prediction)[cx,]/table(evalSet$category)[cx] }
  confusionMat <- sapply(seq(length(levels(stx$category)),by=-1,to=1), getCnfCat)
  
  image(1:Ncat,1:Ncat,confusionMat, col=gray.colors(100), zlim=c(0,1), xlab='decoded cateogry',ylab='actual category', xaxt='n',yaxt='n')
  axis(1,at=1:Ncat)
  axis(2,at=1:Ncat, labels=seq(Ncat,by=-1,to=1) )
  
  mnErrorT <- mean(abs( c(trainingSet$category) -      c(trainingSet$prediction )  ))
  mnErrorE <- mean(abs( c(evalSet$category    ) -      c(evalSet$prediction     )  ))
  mnErrorR <- mean(abs( c(evalSet$category    ) - mean(c(evalSet$category       )) ))
  
  resList <- list(tbl=tbl, confusionMat=confusionMat, pcorT=mean(trainingSet$correct), pcorE=mean(evalSet$correct), mnErrorT=mnErrorT, mnErrorE=mnErrorE,mnErrorR=mnErrorR)
  print('all done')
  return(resList)
}