HLCor.obj <-
function(ranefParsVec,skeleton,HLCor.obj.value="p_bv",trace=NULL,family=gaussian(),...) { ## name of first arg MUST differ from names in dotlist...
  mc <- match.call(expand.dots=TRUE) ## (1) expand.dots added 11/04/2014 for the mutlinomial... eval 
  dotlist <- list(...)
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
  hlfit <- do.call("HLCor",HLCor.args)
  aphls <- hlfit$APHLs
  resu <- aphls[[HLCor.obj.value]]
  readable <- unlist(toCanonical(ranPars=forGiven,corr.model=dotlist$`corr.model`,checkComplete=FALSE)$ranPars) ## FR->FR use of dotlist...
  verif <- c(unlist(aphls),hlfit$lambda,hlfit$phi,readable,ranefParsVec) ## hlfit$phi may be NULL
  if (is.character(trace)) {
    write(verif,file=trace,ncolumns=length(verif),append=T) ## the file is unlink'ed in corrHLfit()  
  }
  return(resu) #
}
