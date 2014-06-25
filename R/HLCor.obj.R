HLCor.obj <-
function(ranefParsVec,skeleton,HLCor.obj.value="p_bv",trace=NULL,family=gaussian(),...) { ## name of first arg MUST differ from names in dotlist...
  mc <- match.call(expand.dots=TRUE) ## (1) expand.dots added 11/04/2014 for the mutlinomial... eval 
  dotlist <- list(...)
  ## (2) (written long before (1)and (3) )
  ## potentially used by getCall(object) in update.HL... if HLfit was called by HLCor through a do.call() this contains the body of the function 
  ## Pour resoudre le probleme de memoire (mais pas du programmeur): 
  ## In that case HLCor removes this from the HLfit object and gives its own call. Otherwise we can improve a bit by 
  ## mc[[1]] <-  call("HLfit")[[1]] ## replace the body with the call; eval(mc) will still work
  ## but all other arguments are still evaluated... cf HLCor
  ## (3) quand la fonction est appelee par optim there is no appropriate info about the called objective function in mc... samefor function argument
  if ( inherits(mc$data,"list")) {
    ## then processed should already be a list
    family <- mc$family
    fitlist <- lapply(seq_len(length(mc$data)),function(it){
      locmc <- mc
      locmc[[1L]] <- as.name("HLCor.obj") ## replaces "f" !
      locmc$ranefParsVec <- ranefParsVec ## replaces "arg" !
      if (family$family=="multi") locmc$family <- family$binfamily
      locmc$distMatrix <- mc$distMatrix[[it]]
      locmc$uniqueGeo <- mc$uniqueGeo[[it]]
      locmc$data <- mc$data[[it]]
      locmc$processed <- mc$processed[[it]]
      eval(locmc) ## this will execute all the code below starting from dotlist <- list(...) 
    })
    resu <- sum(unlist(fitlist))
    if (is.character(trace)) {
      verif <- paste("#global:",ranefParsVec,resu) 
      write(verif,file=trace,append=T) ## the file is unlink'ed in corrHLfit()  
    }
    return(resu)
  }
  HLCor.formals <- names(formals(HLCor))
  HLfit.formals <- names(formals(HLfit))
  designL.formals <- names(formals(designL.from.Corr))
  makescaled.formals <- names(formals(make.scaled.dist))
  HLnames <- (c(HLCor.formals,HLfit.formals,designL.formals,makescaled.formals))  ## cf parallel code in corrHLfit
  HLCor.args <- dotlist[intersect(names(dotlist),HLnames)]
  forGiven <- relist(ranefParsVec,skeleton) ## given values of the optimized variables
  HLCor.args$ranPars[names(forGiven)] <- forGiven ## do not wipe out other fixed, non optimized variables
  HLCor.args$family <- family 
  if (is.character(trace)) {
    if(.SpaMM$TRACE.UNLINK) unlink("HLCor.args.*.RData")
    zut <- paste(ranefParsVec,collapse="")  
    save(HLCor.args,file=paste("HLCor.args.",zut,".RData",sep="")) ## for replicating the problem
  }
  hlfit <- do.call(HLCor,HLCor.args)
  aphls <- hlfit$APHLs
  resu <- aphls[[HLCor.obj.value]]
  verif <- c(unlist(aphls),hlfit$lambda,hlfit$phi,ranefParsVec) ## hlfit$phi may be NULL
  if (is.character(trace)) {
    write(verif,file=trace,ncolumns=length(verif),append=T) ## the file is unlink'ed in corrHLfit()  
  }
  return(resu) #
}
