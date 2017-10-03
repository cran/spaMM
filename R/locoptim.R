.corners_objfn <- function(v,objfn,notoptimObjfnCall.args) {
  notoptimObjfnCall.args$ranefParsVec <- v
  do.call(objfn,notoptimObjfnCall.args) ## this reconstructs a list of the form of initvec, then adds it to other anyHLCor$ranPars information
}

## this function wraps optim, optimize and NLOPT_LN_BOBYQA
# It by default minimizes its function; to maximize, use maximize=TRUE (minimization used only in confint 29/08/2014)
## the objective function 'objfn' is by default the objective function HLCor.obj (currently not by name bc do.call(optimize,...name)) fails)
##   the first arg of this fn must be ranefParsVec
.locoptim <- function(init.optim,LowUp,anyObjfnCall.args,
                     control, 
                     objfn=HLCor.obj, #FR->FR do.call("optim", c(<list>, list(fn = objfn))) does not work with a char string
                     maximize=FALSE) {
  Optimizer <- control$Optimizer
  if (is.null(Optimizer)) Optimizer <- "NLOPT_LN_BOBYQA"
  processedHL1 <- .getProcessed(anyObjfnCall.args$processed,"HL[1]",from=1L) 
  initvec <- unlist(init.optim)
  if ( ! is.null(anyObjfnCall.args$ranefParsVec) ) stop("! is.null(anyObjfnCall.args$ranefParsVec)") ## catch programming errors
  optimizerObjfnCall.args <- notoptimObjfnCall.args <- anyObjfnCall.args
  notoptimObjfnCall.args$ranefParsVec <- initvec ## logscale, inherits names from init.optim
  lower <- unlist(LowUp$lower); upper <- unlist(LowUp$upper)
  if (length(initvec)==0L) {
    optPars <- NULL
  } else if (length(initvec)==1L) { ## one-dimensional optimization; perhaps easily trapped in local maxima ?
    ## corrHLfit adjacency by default reaches this point
    anyOptimize.args <- optimizerObjfnCall.args
    anyOptimize.args$interval <- c(lower,upper)
    anyOptimize.args$tol <- control$optimize$tol
    if (maximize) anyOptimize.args$maximum <- TRUE
    optr <- do.call(optimize,c(anyOptimize.args,list(f=objfn)))  ## optimize p_bv
    if(maximize) {
      optPars <- relist(optr$maximum,init.optim)
    } else optPars <- relist(optr$minimum,init.optim)
    attr(optPars,"method") <- "optimize"
    ## HLCor.args$ranPars[names(optPars)] <- optPars 
  } else if (Optimizer=="NLOPT_LN_BOBYQA") {
    objfn_nloptr <- function(x,anyObjfnCall.args) { ## all functions should have the same args.
      arglist <- c(list(ranefParsVec=x),anyObjfnCall.args)
      return( - do.call(objfn,arglist))
    }
    nloptr_controls <- spaMM.getOption("nloptr")
    nloptr_controls[names(control$nloptr)] <- control$nloptr ## Overwrite defaults with any element of $nloptr
    lowerb <- unlist(lower)
    upperb <- unlist(upper)
    optr <- nloptr::nloptr(x0=initvec,eval_f=objfn_nloptr,lb=lowerb,ub=upperb,
                   opts=nloptr_controls,anyObjfnCall.args=anyObjfnCall.args)
    while (optr$status==5L) { ## 5 => termination bc maxeval has been reached
      prevlik <- optr$objective
      reinit <- pmax(lowerb,pmin(upperb,optr$solution))
      optr <- nloptr::nloptr(x0=reinit,eval_f=objfn_nloptr,lb=lowerb,ub=upperb,
                             opts=nloptr_controls,anyObjfnCall.args=anyObjfnCall.args)
      loc_ftol <- max(1e-8, optr$options$ftol_abs)
      if (- optr$objective < - prevlik+loc_ftol) break ## no progress in <= maxeval iterations
    }
    optPars <- relist(optr$solution,init.optim)
    ## full optr is big. We take out the two items that contribute much to saveSize:
    optr$eval_f <- NULL
    optr$nloptr_environment <- NULL
    optPars <- structure(optPars,method="nloptr",optr=optr) 
  } else {
    parscale <- (upper-lower) ## unlist because some list elements may have length >1 (eg if rho has several     ###### basic unrefined fit for initvec...
    # nlminb code removed 11/2016 ## maxcorners code removed in v2.1.94
    control_optim <- list(parscale=parscale,factr=1e9) ## factr was the stricter 1e8 up to 23/01/13
    if (maximize) control_optim$fnscale <- -1     
    control_optim[names(control$optim$control)] <- control$optim$control ## ...which may be overwritten 
    ## lower and upper were missing 7/11/2014!
    anyBaseOptim.args <- c(optimizerObjfnCall.args,list(par=initvec,lower=lower,upper=upper,control=control_optim,method=Optimizer)) 
    optr <- do.call("optim",c(anyBaseOptim.args,list(fn=objfn))) ## optimize HLCor.obj()'s 'objective'
    optPars <- relist(optr$par,init.optim)
    ## full optr is big. We take out the two items that contribute much to saveSize:
    optr$eval_f <- NULL
    optr$nloptr_environment <- NULL
    attr(optPars,"optr") <- optr  
    attr(optPars,"method") <- "optim"  
    ## HLCor.args$ranPars[names(optPars)] <- optPars 
  }  
  return(optPars)
} ## end def locoptim
