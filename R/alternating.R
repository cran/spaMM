alternating <-
function(init.optim,LowUp,anyOptim.args,maxIter,ranPars,HLCor.args,trace,Optimizer,optimizers.args,corners) {
  nam <- names(init.optim)
  if (any(c("trPhi","trLambda") %in% nam )) {
    mess <- pastefrom("Dispersion parameters non allowed in 'init.corrHLfit' with alternating algorithm.",prefix="(!) From ")
    stop(mess)
  }
  initcorr <- init.optim[nam %in% c("trRho","trNu","Nugget","ARphi")]
  HLfitLowUp <- LowUp
  HLfitLowUp$lower[c("trRho","trNu","Nugget","ARphi")] <- NULL
  HLfitLowUp$upper[c("trRho","trNu","Nugget","ARphi")] <- NULL
  corrLowUp <- LowUp
  corrLowUp$lower[c("trPhi","trLambda")] <- NULL
  corrLowUp$upper[c("trPhi","trLambda")] <- NULL
  anycorrOptim.args <- anyOptim.args
  iter <- 0
  conv <- 1
  currentLik <- -Inf
  while (iter < maxIter && conv > 1e-5 ) { ## if alternating: alternate HLCor and locoptim
    ranPars[names(initcorr)] <- initcorr
    attr(ranPars,"type")[names(initcorr)] <- "var"
    HLCor.args$ranPars <- ranPars
    oldLik <- currentLik
    if (is.character(trace$file)) {
      if(.SpaMM$TRACE.UNLINK) unlink("HLCor.args.*.RData")
      zut <- paste(unlist(initcorr),collapse="")  
      save(HLCor.args,file=paste("HLCor.args.",zut,".RData",sep="")) ## for replicating the problem
    }
    givencorr <- do.call(HLCor,HLCor.args) ## optim disp and beta given corr param
    currentLik <- givencorr$APHLs$p_v ## iterations maximize p_v
    conv <- currentLik-oldLik
    anycorrOptim.args$ranPars$lambda <- givencorr$lambda
    anycorrOptim.args$ranPars$phi <- givencorr$phi
    #### anycorrOptim.args$etaFix <- list(beta=givencorr$fixef,v_h=givencorr$v_h) ## that's what LeeN01sm say, but this does not work
    anycorrOptim.args$etaFix <- list(beta=givencorr$fixef) 
    initcorr <- locoptim(initcorr,corrLowUp,anyOptim.args=anycorrOptim.args,trace,Optimizer=Optimizer,optimizers.args=optimizers.args,corners=corners) 
    iter <- iter+1
  }
  optPars <- c(initcorr,givencorr$lambda,givencorr$phi)
  optPars
}
