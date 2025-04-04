combineFeaturePopulation <- function(valnd, df=df, removeBaseline=T, signal='M', epoch='e'){ 
  
  dataList <- lapply( valnd, getFeaturePopulation,   df=df, removeBaseline=removeBaseline, signal=signal, epoch=epoch)

  xpp <- dataList[[1]]
  for (i in 2:(length(dataList))){#                  
    print(i)
    xpp <- rbind( xpp, dataList[[i]])
  }
  return(xpp)
}