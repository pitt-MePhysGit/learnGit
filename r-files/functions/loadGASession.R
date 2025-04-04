loadGASession <- function(n, who='all',extn='mean',what='csd',fextn='',xnd=1, ynd=NA, ganrm=T){
  
  if(ganrm){## normalize results by the max(abs( of the grand average ))
    mnGA       <- loadGASession(n, who=who,extn='mean',what=what,fextn=fextn,xnd=xnd, ynd=ynd,ganrm=F)
    mnGAnrm    <- apply(mnGA,2,function(x){max(abs(x))}) 
    mnGAnrmMat <- array( mnGAnrm[rep( 1:dim(mnGA)[2],each=dim(mnGA)[1] )], dim(mnGA) )
  }
  
  dfline <- df[n,]
  exp    <- df$session[n]
  print(exp)
  #print(dfline)
  
  if(df$probe[n]==2)
    exp <- df$exp[n]
  
  resFolder <- paste(trda(), exp,'/', sep='')
  fileName  <- paste(resFolder,who,'_',extn,'_', what, fextn,'.rda',sep='')
  load(fileName)
  
  if (length(mn)>100)# just one time-series
    tmp <- mn
  
  if (length(mn)<100){# array of time-series
    
    if (is.na(ynd))
      tmp <- mn[[xnd]]
    
    if (!is.na(ynd))
      tmp <- mn[[xnd,ynd]]
  }
  if(ganrm)
    tmp <- tmp / mnGAnrmMat
  
  #print(dim(tmp))
  if(dim(tmp) )
    
    
    #tmp <- tmp[,dim(tmp)[2]:1]
    return(tmp)
}
