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
  if (paste(llc[[1]]) =="corrHLfit") { ## optimInfo always present
    lc <- get_HLCorcall(object,fixed=getCall(object)$ranFix) ## explicit ranFix: cf get_HLCorcall source; 
                                          ## ranPars are optimized again as controlled by code below
  } else if (paste(llc[[1]]) =="fitme") { ## optimInfo always present
    if (inherits(object,"HLCor")) {
      lc <- get_HLCorcall(object,fixed=getCall(object)$fixed)
    } else {
      lc <- object$call ## HLfit call
    }
  } else lc <- llc 
  HL <- object$HL
  lik <- switch(paste(HL[1]),
                "0"="hlik",
                "1"="p_v",
                stop(paste("confint does not yet handle HLmethod",paste(HL,collapse=" "),"(or ",llc$HLmethod,").",sep=" ")))
  lc$control.HLfit$intervalInfo$fitlik <- object$APHLs[[lik]]
  lc$control.HLfit$intervalInfo$targetlik <- object$APHLs[[lik]]-dlogL
  lc$control.HLfit$intervalInfo$MLparm <- object$fixef[parm]
  lc$control.HLfit$intervalInfo$parm <- parm
  if (! is.null(lc$processed)) .assignWrapper(lc$processed,"LevenbergM['force'] <- FALSE") ## inhibits LevM for confint 
  beta_cov <- .get_beta_cov_any_version(object)
  beta_se <- sqrt(diag(beta_cov))
  lc$control.HLfit$intervalInfo$asympto_abs_Dparm <- asympto_abs_Dparm <- znorm* beta_se
  optimInfo <- attr(object,"optimInfo") ## may be NULL
  trTemplate <- optimInfo$`optim.pars` ## may be NULL if optimInfo is NULL or if  optimInfo is not NULL but no par as outer optimized
  if ( ! is.null(trTemplate)) {
    olc <- lc ## olc is working copy
    LUarglist <- optimInfo$LUarglist
    attr(trTemplate,"method") <- NULL
    if (paste(lc[[1]])=="HLCor") {
      ## good starting values are important... important to use canonizeRanPars as in HLCor
      attr(trTemplate,"RHOMAX") <- LUarglist$RHOMAX ## an ugly bit of coding, but canonizeRanpars()  expects them 
      attr(trTemplate,"NUMAX") <- LUarglist$NUMAX
      objfn <- function(ranefParsVec, anyHLCor_obj_args=NULL, HLcallfn.obj=NULL) { 
        ## ca serait mieux de pas avoir de contrainte la dessus et de pvr nommer l'arg trParsVec
        ## bc HLCor call uses transformed scale for ranPars
        ranefParsList <- relist(ranefParsVec,trTemplate)
        olc$ranPars[names(trTemplate)] <- ranefParsList 
        locfit <- eval(as.call(olc)) ## HLCor call with given ranefParsVec
        resu <- (posforminimiz)* locfit$fixef[parm]
        attr(resu,"info") <- locfit$APHLs$p_v 
        ## attribute lost by optim but otherwise useful for debugging 
        return(resu) ## return value to be optimized is a parameter value, not a likelihood
      }
    } else if (paste(lc[[1]])=="HLfit") {
      objfn <- function(ranefParsVec, anyHLCor_obj_args=NULL, HLcallfn.obj=NULL) { 
        ## bc locoptim expects a fn with first arg ranefParsVec
        ranefParsList <- relist(ranefParsVec,trTemplate)
        olc$ranFix[names(trTemplate)] <- ranefParsList 
        locfit <- eval(as.call(olc)) ## HLfit call with given ranefParsVec
        resu <- (posforminimiz)*locfit$fixef[parm]
        attr(resu,"info") <- locfit$APHLs$p_v 
        ## attribute lost by optim but otherwise useful for debugging 
        return(resu) ## return value to be optimized is a parameter value, not a likelihood
      }
    }
    ad_hoc_fn <- function(posforminimiz) {
      if (length(unlist(trTemplate))==1L) { ## .new_locoptim cannot be used in 1D case bc -> new_locoptim -> optimize does not use objfn_locoptim!
        optr <- optimize(f=objfn,interval=unlist(LowUp)) # objfn uses posforminimiz in its definition 
        return(optr$m) ## minimum or maximum...
      } else { # maximization (profiled out pars) may be multidimensional, 
        optr <- .new_locoptim(init.optim=trTemplate,LowUp=LowUp,
                              objfn_locoptim=objfn, # uses posforminimiz in its definition 
                              # default HLcallfn.obj="HLCor.obj"
                              anyHLCor_obj_args=anyObjfnCall.args,
                              control=list())  
        return(optr[names(trTemplate)])
      }
    }
    canonTemplate <- .canonizeRanPars(ranPars=trTemplate,
                                      corr_types=LUarglist$corr_types, #corr.model=lc$`corr.model`,
                                     checkComplete=FALSE)
    canonTemplate <- canonTemplate$ranPars
    LUarglist$canon.init <- canonTemplate
    LowUp <- do.call(.makeLowerUpper,LUarglist)
  }
  ## lowerfit
  fac <- 1L 
  warnori <- options(warn=-1)
  prevmsglength <- 0L
  while(fac < 1e6) {
    lc$control.HLfit$intervalInfo$init <- (object$fixef-asympto_abs_Dparm/fac)[parm]
    if (! is.null(trTemplate)) {
      anyObjfnCall.args <- as.list(lc[-1L]) ## includes processed, ranPars, controlS.dist, control.HLfit...
      anyObjfnCall.args$skeleton <- trTemplate
      # The objective function 'objfn' returns the confint bound given the corr pars. Thus locoptim maximizes the confint bound over the the corr pars
      olc <- lc ## that's olc that is used in the objective fn !
      posforminimiz <- 1 ## defined in the envir where objfn is defined... (bad style)
      if (paste(lc[[1]])=="HLCor") {
        olc$ranPars[names(trTemplate)] <- ad_hoc_fn() 
      } else {
        olc$ranFix[names(trTemplate)] <- ad_hoc_fn()
      }
      ## recover fit for optimized params (must use call witSh intervalInfo and LevenbergM=FALSE)
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
    if (! is.null(trTemplate)) {
      olc <- lc
      posforminimiz <- -1 ## maximization
      if (paste(lc[[1]])=="HLCor") {
        olc$ranPars[names(trTemplate)] <- ad_hoc_fn() 
      } else {
        olc$ranFix[names(trTemplate)] <- ad_hoc_fn()
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
