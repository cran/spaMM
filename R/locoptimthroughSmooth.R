locoptimthroughSmooth <-
function(pargrid,anyHLCor.args,prevPtls=NULL,control.smooth=list()) {
  anyHLCor.args$processed$SEMseed <- NULL ## SEMseed bien pour controler individuellement un SEM mais pas une série 
  ranges <- apply(pargrid,2,range)
  LowUp <- apply(ranges,1,as.list)
  names(LowUp) <- c("lower","upper")
  lower <- LowUp$lower
  upper <- LowUp$upper
  grid.obj <- apply(pargrid,1,function(v) {
    anyHLCor.args$ranPars[names(lower)] <- relist(v,lower)
    hlcor <- do.call(HLCor,anyHLCor.args) ## this reconstructs a list of the form of initvec, then adds it to other anyHLCor$ranPars information
    c(logL=hlcor$APHLs$p_v,lambda=hlcor$lambda)
  })
  nrepl <- control.smooth$nrepl
  if (anyHLCor.args$HLmethod == "SEM") {
    if (is.null(nrepl)) {
      nrepl <- Inf ## default is to duplicate all simuls !
    }
    if (nrepl == 0 && (length(control.smooth$ranFix) != length(lower))) {
      message("(!) From locoptimthroughSmooth: suspect nrepl==0 for stochastic simulation with correlation parameters to be estimated.")
    }
  } else {
    if (is.null(nrepl)) {
      nrepl <- 0
    }
    if (nrepl>0) message("(!) From locoptimthroughSmooth: suspect nrepl>0 for deterministic simulation.")
  }  
  if (nrepl>0) {
    subpargrid <- pargrid[sample(nrow(pargrid),min(nrepl,nrow(pargrid))),,drop=FALSE]
    subgrid.obj <- apply(subpargrid,1,function(v) {
      anyHLCor.args$ranPars[names(lower)] <- relist(v,lower)
      hlcor <- do.call(HLCor,anyHLCor.args) ## this reconstructs a list of the form of initvec, then adds it to other anyHLCor$ranPars information
      c(logL=hlcor$APHLs$p_v,lambda=hlcor$lambda)
    })
    forSmooth <- cbind(rbind(pargrid,subpargrid),rbind(t(grid.obj),t(subgrid.obj)))      
  } else {
    forSmooth <- cbind(pargrid,t(grid.obj))      
  }
  forSmooth <- data.frame(forSmooth)
  if ( ! is.null(prevPtls)) forSmooth <- rbind(prevPtls,forSmooth) ## older points first!
  RHS <- paste(names(lower),collapse="+")
  form <- as.formula(paste("logL ~ 1+Matern(1|",RHS,")"))
  ## ****** predict a (profile) likelihood surface for the correlation parameters ******
  ## phi here represent varSEMandMCint the variance of MC integration and that of the SEM algo
  ##  phi = a [SEM]+ b lambda [MCint]: 
  ## i e log(phi) = log(a+ b lambda), not a +b log(lambda) => code pas ideal
  ## fortunately lambda may be roughly indep of the corrpars
  ## Further the IRLS in glm.fit is fussy... 
  ## et le fit des dev resids peut etre surprenant, cf plot(var explicative,log(resp)) dans dispGammaGLM -> glm(resp ~...)
  #Krigobj <- corrHLfit(form,data=forSmooth,resid.formula= ~log(lambda),init.corrHLfit=list(rho=rep(1,length(lower))),ranFix=list(nu=4))
  ranFix <- list(nu=4)
  if (anyHLCor.args$HLmethod != "SEM") ranFix$phi <- 1e-06 ## interpolation
  ranFix[names(control.smooth$ranFix)] <- control.smooth$ranFix
  init.corrHLfit <- list(rho=rep(1,length(lower))) ## important as it gives the length of rho to corrHLfit
  init.corrHLfit[names(ranFix)] <- NULL
  Krigobj <- corrHLfit(form,data=forSmooth,resid.formula= ~1,init.corrHLfit=init.corrHLfit,ranFix=ranFix)
#  browser()
  ### quick check of variance of logL estimation:
  if (FALSE) {
   essai <- forSmooth[do.call(order,forSmooth[,seq_len(length(lower))]),,drop=FALSE]
   diffs <- apply(essai,2,diff)
   ## next line selects diff lines for identical coordinates and takes from them the diff value for logL 
   logLdiffs <- diffs[apply(diffs[,seq_len(length(lower))]==rep(0,length(lower)),1,all),length(lower)+1]
   var(logLdiffs)/2 ## since diff= sum(eps_1+eps_2)
  }
  ###
  ## ****** optimize in the predicted likelihood surface ******
  predictions <- predict(Krigobj)
  predictions <- predictions[order(predictions$fitted,decreasing=TRUE),]
  initvec <- predictions[1,names(lower)]
  ## redefines lower, upper, for maximization
  ranges <- apply(forSmooth[,names(lower),drop=FALSE],2,range)
  LowUp <- apply(ranges,1,as.list)
  names(LowUp) <- c("lower","upper")
  lower <- LowUp$lower
  upper <- LowUp$upper
  ##
  if (length(initvec)==1) {
    optr <- optimize(function(v) {predict(Krigobj,v)$fitted},maximum=TRUE,lower=unlist(lower),upper=unlist(upper)) 
    optr$par <- optr$maximum
    optr$value <- optr$objective
    optr$maximum <- NULL
    optr$objective <- NULL
  } else { 
    optr <- optim(initvec,function(v) {predict(Krigobj,v)$fitted},method="L-BFGS-B",
                  control=list(fnscale=-1),lower=unlist(lower),upper=unlist(upper))
  }
  names(optr$par) <- names(lower)
  attr(predictions,"fittedPars") <- names(lower) 
  attr(predictions,"respName") <- "fitted" 
  attr(predictions,"MSy") <- Krigobj$phi ## FR->FR ecrire un extracteur pour phi... 
  optr$predictions <- predictions 
#  ranges <- apply(uppersurf,2,range)
#  nextLowUp <- apply(ranges,1,as.list)
#  names(nextLowUp) <- c("lower","upper")
  ######
#  optr$nextLowUp <- nextLowUp
  optr$Krigobj <- Krigobj ## 
  optr$forSmooth <- forSmooth
  return(optr)
}
