LRT <-
function(object,object2,boot.repl=0) { ## compare two HM objects
  info <- compare.model.structures(object,object2)
  nullm <- info$nullm; fullm <- info$fullm; testlik <- info$testlik;df <- info$df
  LRTori <- 2*(fullm$APHLs[[testlik]]-nullm$APHLs[[testlik]])
  pvalue <- 1-pchisq(LRTori,df=df) ## but not valid for testing null components of variance
  resu <- list(nullfit=nullm,fullfit=fullm,basicLRT = data.frame(LR2=LRTori,df=df,pvalue=pvalue)) ## format appropriate for more tests  
  if (boot.repl>0) {
    simbData <- nullm$data
    if (tolower(nullm$family$family)=="binomial") {
      form <- attr(nullm$predictor,"oriFormula") ## this must exists...  
      if (is.null(form)) {
        mess <- pastefrom("a 'predictor' object must have an 'oriFormula' member.",prefix="(!) From ")
        stop(mess)
      }
      exprL <- as.character(form[[2]][[2]]) 
      exprR <- as.character(form[[2]][[3]]) 
    }
    if (boot.repl<100) print("It is recommended to set boot.repl>=100 for Bartlett correction",quote=FALSE)
    aslistfull <- as.list(getCallHL(fullm)) 
    ## problem is for corrHLfit etc this is the call of the final HLfit call with $processed and a lot of missing original arguments  
    aslistfull$processed <- NULL ## may capture bugs 
    aslistnull <- as.list(getCallHL(nullm))
    aslistnull$processed <- NULL ## may capture bugs
    computeBootRepl <- function() {
      ## draw sample
      newy <- simulate(nullm)  ## only a vector of response value ## cannot simulate all samples in one block since some may not be analyzable  
      if (tolower(nullm$family$family)=="binomial") {
        ## c'est bouseux: soit j'ai (pos, neg) et le remplacement est possible
        ##    soit j'ai (pos,ntot -pos) et le 2e remplacment n'est pas poss (et pas necess)
        ##    aussi (ntot - pos, pos) ...
        ## would be simple if always ntot-pos, but how to control this ? 
        if (length(exprL)==1) simbData[[exprL]] <- newy 
        if (length(exprR)==1) simbData[[exprR]] <- nullm$weights - newy                    
        ## if (length(exprR)! =1) exprRdoes not correspond to a column in the data;frmae so there is no column to replace                     
      } else {simbData[[as.character(nullm$predictor[[2]])]] <- newy}
      ## analyze under both models
      aslistfull$data <- simbData
      aslistnull$data <- simbData
      fullfit <- (eval(as.call(aslistfull)))
      if (inherits(fullfit,"try-error")) return(fullfit) ## eg separation in binomial models
      nullfit <- try(eval(as.call(aslistnull)))
      if (inherits(nullfit,"try-error")) return(nullfit) 
      ## return pair of likelihoods
      if (inherits(fullfit,"HLfitlist")) {
        return(c(attr(fullfit,"APHLs")[[testlik]],attr(nullfit,"APHLs")[[testlik]]))
      } else return(c(fullfit$APHLs[[testlik]],nullfit$APHLs[[testlik]]))
    }
    bootLs<-matrix(,nrow=boot.repl,ncol=2) 
    colnames(bootLs) <- paste(c("full.","null."),testlik,sep="")
    msg <- "bootstrap replicates: "
    msglength <- nchar(msg)
    cat(msg)
    t0 <- proc.time()["user.self"]
    for (ii in 1:boot.repl) {
      locitError <- 0
      repeat { ## for each ii!
        bootrep <- (computeBootRepl())
        if ( ! inherits(bootrep,"try-error")) { 
          bootLs[ii,] <- bootrep
          break ## replicate performed, breaks the repeat
        } else { ## there was one error
          locitError <- locitError + 1
          if (locitError>10) { ## to avoid an infinite loop
            stop("Analysis of bootstrap samples fails repeatedly. Maybe no statistical information in them ?")
          } ## otherwise repeat!
        }
      } 
      tused <- proc.time()["user.self"]-t0
      ttotal <- tused* boot.repl/ii
      if (interactive()) {
        for (bidon in 1:msglength) cat("\b")
        msg <- paste("Estimated time remaining for bootstrap: ",signif(ttotal-tused,2)," s.",sep="")
        msglength <- nchar(msg)
        cat(msg)
      } else {
        cat(ii);cat(" ")
        if ((ii %% 40)==0) cat("\n")
      }
    }
    cat("\n")
    bootdL <- bootLs[,1]-bootLs[,2]
    meanbootLRT <- 2*mean(bootdL)
    resu <- c(resu,list(rawBootLRT = data.frame(LR2=LRTori,df=df,pvalue=(1+sum(bootdL>=LRTori/2))/(boot.repl+1)))) ## format appropriate for more tests  
    LRTcorr <- LRTori*df/meanbootLRT
    resu <- c(resu,list(BartBootLRT = data.frame(LR2=LRTcorr,df=df,pvalue=1-pchisq(LRTcorr,df=df)))) ## format appropriate for more tests  
    bootInfo <- list(meanbootLRT = meanbootLRT,bootreps = bootLs)
    resu <- c(resu,list(bootInfo=bootInfo)) ## keeps the sublist structure, which is not compatible with hglmjob.R...  
  }
  class(resu) <- c("fixedLRT",class(resu)) 
  return(resu)
}
