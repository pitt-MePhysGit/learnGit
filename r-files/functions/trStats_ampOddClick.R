## ==================================================================
##  this is the latest routine that uses the Anova from the car library
##  it uses Type 2 sums of square and includes all possible interactions
##  should be doing the same thing that decide.trStats
library(car)
trStats.ampOddClick <- function(input,width=.050,what='spk',chNr=1,
                              eval.at=seq(-.5,1,by=.02),alpha=0.05,
                              variables=c('sp','isiLog','ampInd','deltaAmpInd', 'absDeltaAmpInd','repNr','trialType','trialNr'),
                              subs=NULL,extn='',show=T, fix=FALSE){
  
  if (is.null(subs))
    subs <- rep(TRUE, dim(input$xpp)[1])
  
  if(identical(what,'spk'))
    filename <- paste(rda(),'glm_work/ampOddClick.trStatsCar_',strsplit(input$exp,'/')[[1]][1],
                      '_ch',input$ch,'_cell',input$cell,'_',extn,'width',1000*width,'.rda',sep='')
  
  
  
  if(!identical(what,'spk'))
    filename <- paste(rda(),'glm_work/ampOddClick.trStatsCar_',strsplit(input$exp,'/')[[1]][1],
                      '_ch',input$ch,extn,'.rda',sep='')
  
  
  ## normalize numeric predictor variables
  input$xpp$isiLog           <- input$xpp$isiLog           - na.mean(input$xpp$isiLog[which(subs)]          ) 
  input$xpp$deltaAmpInd      <- input$xpp$deltaAmpInd      - na.mean(input$xpp$deltaAmpInd[which(subs)]   ) 
  input$xpp$absDeltaAmpInd   <- input$xpp$absDeltaAmpInd   - na.mean(input$xpp$absDeltaAmpInd[which(subs)]) 
  input$xpp$trialNr          <- input$xpp$trialNr          - na.mean(input$xpp$trialNr[which(subs)]) 
  input$xpp$repNr            <- input$xpp$repNr            - na.mean(input$xpp$repNr[which(subs)]) 
  
  
  all <- list()
  if(length(variables)>0){
    ##============================== 
    this.test <- function(subinput){
      
      if(max(subinput$xpp$sp)==0 && min(subinput$xpp$sp)==0) ## not so nice trick to avoid a singular system
        input$xpp$sp[1] <- 1
      
      res         <- list()
      res$aovFull <- aov(sp ~ 1+., data=subinput$xpp[,variables])
      res$car     <- Anova(res$aovFull,type=2)
      
      lm0           <- aov(sp ~ -1., data=subinput$xpp[,variables])
      res$lm1       <- aov(sp ~  1., data=subinput$xpp[,variables])
      res$blcar     <- anova(lm0,res$lm1)
      return( res )
    }## ============================================= this is the test that is handed to tres.test2
    res <- tres.test2(input,what=what,main='trial number',eval.at=eval.at,test=this.test,width=width,subs=subs,show=show,alpha=alpha,chNr=chNr)
  }
  
  save(res, file=filename)
  return(res)
}



