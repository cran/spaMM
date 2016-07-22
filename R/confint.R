confint.HLfit <- function(object,parm,level=0.95,verbose=TRUE,...) {
  dlogL <- qchisq(level,df=1)/2
  znorm <- qnorm((1+level)/2)
  if (is.character(parm)) {
    whichcol <- which(names(object$fixef)==parm)
    if (length(whichcol)==0L) stop("Parameter not in the model")
    attr(parm,"col") <- whichcol 
  } else {
    parmcol <- parm
    if (parm > length(object$fixef)) stop("'parm' not compatible with # of fixed-effects coefficients")
    parm <- names(object$fixef)[parmcol]
    attr(parm,"col") <- parmcol
  }
  llc <- as.list(getCall(object))
  if (paste(llc[[1]]) =="corrHLfit") {
    lc <- attr(object,"HLCorcall") ## optimInfo+HLCorcall always present
  } else if (paste(llc[[1]]) =="fitme") { #optimInfo always present
    lc <- attr(object,"HLCorcall") ## may or may not be present
    if (is.null(lc)) {
      lc <- object$call ## HLfit call
    } else lc <- as.list(lc) ## optimInfo+HLCor
  } else lc <- llc ## no optimInfo
  HL <- object$HL
  lik <- switch(paste(HL[1]),
                "0"="hlik",
                "1"="p_v",
                stop(paste("confint does not yet handle HLmethod",paste(HL,collapse=" "),"(or ",llc$HLmethod,").",sep=" ")))
  lc$control.HLfit$intervalInfo$fitlik <- object$APHLs[[lik]]
  lc$control.HLfit$intervalInfo$targetlik <- object$APHLs[[lik]]-dlogL
  lc$control.HLfit$intervalInfo$MLparm <- object$fixef[parm]
  lc$control.HLfit$intervalInfo$parm <- parm
  lc$control.HLfit$LevenbergM <- FALSE ## simple... but read only in preprocess which will usually not be run...
  if (! is.null(lc$processed)) lc$processed <- setProcessed(lc$processed,"LevenbergM","FALSE") ## same idea... (10/2015) but with return value 07/2016 ! 
  beta_se <- sqrt(diag(object$beta_cov))
  lc$control.HLfit$intervalInfo$asympto_abs_Dparm <- asympto_abs_Dparm <- znorm* beta_se
  optimInfo <- attr(object,"optimInfo")
  if ( ! is.null(optimInfo)) {
    olc <- lc ## olc is working copy
    LUarglist <- optimInfo$LUarglist
    trTemplate <- optimInfo$`optim.pars` ## optr$value in transformed scale 
    attr(trTemplate,"method") <-NULL
    if (paste(lc[[1]])=="HLCor") {
      ## good starting values are important... important to use canonizeRanPars as in HLCor
      attr(trTemplate,"RHOMAX") <- LUarglist$RHOMAX ## an ugly bit of coding, but canonizeRanpars()  expects them 
      attr(trTemplate,"NUMAX") <- LUarglist$NUMAX
      objfn <- function(ranefParsVec) { 
        ## bc locoptim expects a fn with first arg ranefParsVec
        ## ca serait mieux de pas avoir de contrainte la dessus et de pvr nommer l'arg trParsVec
        ## bc HLCor call uses transformed scale for ranPars
        ranefParsList <-relist(ranefParsVec,trTemplate)
        olc$ranPars[names(trTemplate)] <- ranefParsList 
        locfit <- eval(as.call(olc)) ## HLCor call with given ranefParsVec
        resu <- locfit$fixef[parm]
        attr(resu,"info") <- locfit$APHLs$p_v 
        ## attribute lost by optim but otherwise useful for debugging 
        return(resu) ## return value to be optimized is a parameter value, not a likelihood
      }
    } else if (paste(lc[[1]])=="HLfit") {
      objfn <- function(ranefParsVec) { 
        ## bc locoptim expects a fn with first arg ranefParsVec
        ranefParsList <-relist(ranefParsVec,trTemplate)
        olc$ranFix[names(trTemplate)] <- ranefParsList 
        locfit <- eval(as.call(olc)) ## HLfit call with given ranefParsVec
        resu <- locfit$fixef[parm]
        attr(resu,"info") <- locfit$APHLs$p_v 
        ## attribute lost by optim but otherwise useful for debugging 
        return(resu) ## return value to be optimized is a parameter value, not a likelihood
      }
    }
    canonTemplate <- canonizeRanPars(ranPars=trTemplate,
                                     corr.model=LUarglist$corr.model, #corr.model=lc$`corr.model`,
                                     checkComplete=FALSE)
    canonTemplate <- canonTemplate$ranPars
    LUarglist$canon.init <- canonTemplate
    LowUp <- do.call(makeLowerUpper,LUarglist)
    loclist <- list(init.optim=trTemplate,LowUp=LowUp,objfn=objfn,anyObjfnCall.args=list(),optimizers.args=list())
  }
  ## lowerfit
  fac <- 1L 
  warnori <- options(warn=-1)
  prevmsglength <- 0L
  while(fac < 1e6) {
    lc$control.HLfit$intervalInfo$init <- (object$fixef-asympto_abs_Dparm/fac)[parm]
    if (! is.null(optimInfo)) {
      olc <- lc
      # The objective function 'objfn' returns the confint bound given the corr pars. Thus locptim maximizes the confint bound over the the corr pars
      optr <- do.call(locoptim,loclist) 
      ## recover fit for optimized params (must use call with intervalInfo and LevenbergM=FALSE)
      if (paste(lc[[1]])=="HLCor") {
        olc$ranPars[names(trTemplate)] <- optr[names(trTemplate)] ## list <- list
      } else {
        olc$ranFix[names(trTemplate)] <- optr[names(trTemplate)] 
      }
      lowerfit <- eval(as.call(olc))
      attr(lowerfit,"optimInfo") <- optimInfo ## expected by summary.HLfit
      ##
    } else lowerfit <- eval(as.call(lc))
    if (is.null(lowerfit$warnings$innerNotConv)) {
      if (fac > 1.1 && verbose) overcat(" ...converged                                                          \n",prevmsglength) 
      break
    } else {
      if (verbose) prevmsglength <- overcat("Convergence problem, trying another starting value for lower bound...",prevmsglength)
      fac <- 2L*fac
    }
  }
  
  options(warnori)
  ## upperfit:
  fac <- 1L
  warnori <- options(warn=-1)
  prevmsglength <- 0L
  while(fac < 1e6) {
    lc$control.HLfit$intervalInfo$init <- (object$fixef+asympto_abs_Dparm/fac)[parm]
    if (! is.null(optimInfo)) {
      olc <- lc
      loclist$maximize <- TRUE ## 
      # The objective function returns the confint bound given the corr pars. Thus locptim maximizes the confint bound over the the corr pars
      optr <- do.call(locoptim,loclist) 
      ## recover HLCor fit for optimized params
      if (paste(lc[[1]])=="HLCor") {
        olc$ranPars[names(trTemplate)] <- optr[names(trTemplate)] ## list <- list
      } else {
        olc$ranFix[names(trTemplate)] <- optr[names(trTemplate)] 
      }
      upperfit <- eval(as.call(olc))
      attr(upperfit,"optimInfo") <- optimInfo ## expected by summary.HLfit
      ##
    } else upperfit <- eval(as.call(lc))
    if (is.null(upperfit$warnings$innerNotConv)) {
      if (fac > 1.1 && verbose) overcat(" ...converged                                                          \n",prevmsglength) 
      break
    } else {
      if (verbose) prevmsglength <- overcat("Convergence problem, trying another starting value for upper bound...",prevmsglength)
      fac <- 2L*fac
    }
  }
  options(warnori)
  interval <- c(lowerfit$fixef[parm],upperfit$fixef[parm])
  names(interval) <- paste(c("lower","upper"),parm)
  if (verbose) print(interval)
  invisible(list(lowerfit=lowerfit,upperfit=upperfit,interval=interval))
}
