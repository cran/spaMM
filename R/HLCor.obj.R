HLCor.obj <-
function(ranefParsVec,skeleton,HLCor.obj.value="p_bv",trace=NULL,...) { ## name of first arg MUST differ from names in dotlist...
  dotlist <- list(...)
  HLCor.formals <- names(formals(HLCor))
  HLfit.formals <- names(formals(HLfit))
  designL.formals <- names(formals(designL.from.Corr))
  HLnames <- (c(HLCor.formals,HLfit.formals,designL.formals))  ## cf parallel code in corrHLfit
  HLCor.args <- dotlist[intersect(names(dotlist),HLnames)]
  forGiven <- relist(ranefParsVec,skeleton) ## given values of the optimized variables
  HLCor.args$ranPars[names(forGiven)] <- forGiven ## do not wipe out other fixed, non optimized variables
  if (is.character(trace)) {
    if(.SpaMM$TRACE.UNLINK) unlink("HLCor.args.*.RData")
    zut <- paste(ranefParsVec,collapse="")  
    save(HLCor.args,file=paste("HLCor.args.",zut,".RData",sep="")) ## for replicating the problem
  }
  hlfit <- do.call(HLCor,HLCor.args)$hlfit
  aphls <- hlfit$APHLs
  resu <- aphls[[HLCor.obj.value]]
  verif <- c(unlist(aphls),hlfit$lambda,hlfit$phi,ranefParsVec) ## hlfit$phi may be NULL
  if (is.character(trace)) {
    write(verif,file=trace,ncolumns=length(verif),append=T) ## the file is unlink'ed in corrHLfit()  
  }
  return(resu) #
}
