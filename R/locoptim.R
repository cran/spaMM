## .new_locoptim() function wraps optim, optimize and NLOPT_LN_BOBYQA
# It uses optimize(HLcallfn.obj .. maximum=TRUE) or OTHERWISE the .objfn_locoptim() wrapper to maximize likelihood:
#  this is messy to control hence the ad_hoc_fn wrapper in confint.HLfit
# The first arg of the objective functions must be ranefParsVec
.objfn_locoptim <- function(x, anyHLCor_obj_args, HLcallfn.obj) { ## the (more or less) default value of .new_locoptim <- function( .. objfn_locoptim .. ) 
  anyHLCor_obj_args$ranefParsVec <- x 
  return( - do.call(HLcallfn.obj, anyHLCor_obj_args))
}

.optim_by_nloptr <- function(lowerb, upperb, initvec, objfn_locoptim, anyHLCor_obj_args, HLcallfn.obj, init.optim, 
                             nloptr_controls, grad_locoptim=NULL) {
  ## this is also called if length(lower)=0 by  (SEM or not) and optPars is then null 
  optr <- nloptr::nloptr(x0=initvec,eval_f=objfn_locoptim,
                         eval_grad_f=grad_locoptim, # ignored with NLOPT_LN_BOBYQA
                         lb=lowerb,ub=upperb,
                         opts=nloptr_controls, anyHLCor_obj_args=anyHLCor_obj_args, HLcallfn.obj=HLcallfn.obj)
  while (optr$status==5L) { ## optr$status=5 => termination bc maxeval has been reached 
    # met status=4: nloptr message in normal termination due toxtol_rel, but is this true ?
    message("maxeval reached in nloptr(); nloptr() called again until apparent convergence of objective.") 
    prevlik <- optr$objective
    reinit <- pmax(lowerb,pmin(upperb,optr$solution))
    optr <- nloptr::nloptr(x0=reinit,eval_f=objfn_locoptim, eval_grad_f=grad_locoptim, lb=lowerb,ub=upperb,
                           opts=nloptr_controls, anyHLCor_obj_args=anyHLCor_obj_args, HLcallfn.obj=HLcallfn.obj)
    loc_ftol <- max(1e-8, optr$options$ftol_abs)
    if (- optr$objective < - prevlik+loc_ftol) break ## no progress in <= maxeval iterations
  }
  return(optr)
}

.optim_by_bobyqa <- function(lowerb, upperb, initvec, objfn_locoptim, anyHLCor_obj_args, HLcallfn.obj, init.optim, 
                             bobyqa_controls=list()) {
  if ( ! requireNamespace("minqa",quietly=TRUE) ) {stop("Package minqa not installed.")}
  optr <- minqa::bobyqa(par=initvec,fn=objfn_locoptim,lower=lowerb,upper=upperb,control=bobyqa_controls,
                        anyHLCor_obj_args=anyHLCor_obj_args, HLcallfn.obj=HLcallfn.obj)
  while(optr$ierr==1L) { #maximum number of function evaluations exceeded
    message("maxeval reached in bobyqa(); bobyqa() called again until apparent convergence of objective.") 
    prevmlik <- optr$fval
    reinit <- pmax(lowerb,pmin(upperb,optr$par))
    optr <- minqa::bobyqa(par=initvec,fn=objfn_locoptim,lower=lowerb,upper=upperb,control=bobyqa_controls,
                          anyHLCor_obj_args=anyHLCor_obj_args, HLcallfn.obj=HLcallfn.obj)
    if (optr$fval > prevmlik-1e-8) break ## not enough progress in <= maxeval iterations
  }
  optr$value <- - optr$fval
  return(optr)
}

.xtol_abs_fn <- function(LowUp,factors=.spaMM.data$options$xtol_abs_factors, rC_transf=.spaMM.data$options$rC_transf) {
  xtol_abs <- vector("list",length(LowUp$lower))
  parnames <- names(LowUp$lower)
  if ("trRanCoefs" %in% parnames) {
    for (st in names(LowUp$lower)) if(st=="trRanCoefs") {
      len <- length(unlist(LowUp$lower[[st]]))
      Xi_ncol <- floor(sqrt(len*2))
      if (rC_transf=="chol") {
        xtol_abs[[st]] <- rep(1e-12,len) # note that order of elements is that of upper.tri 
      } else xtol_abs[[st]] <- c(rep(factors["rcLam"],Xi_ncol),rep(factors["rcCor"],len-Xi_ncol)) # "sph" etc
    } else {xtol_abs[[st]] <- rep(factors["others"],length(unlist(LowUp$lower[[st]])))}
    rng <- unlist(LowUp$upper)-unlist(LowUp$lower)
    xtol_abs <- unlist(xtol_abs)# not clear why
    rng_finite <- is.finite(rng)
    xtol_abs[rng_finite] <- xtol_abs[rng_finite] * rng[rng_finite]
  } else xtol_abs <- factors["abs"]
  return(xtol_abs)
}

# returns optPars which is a list given by relist(.,init.optim), with attributes the optimMethod and (+:- raw) optr 
.new_locoptim <- function(init.optim, LowUp, control, objfn_locoptim, 
                          anyHLCor_obj_args, HLcallfn.obj="HLCor.obj", 
                          user_init_optim, ## only purpose is to check whether (some of the) init.optim comes from explicit user info.
                          grad_locoptim=NULL
) {
  initvec <- unlist(init.optim) 
  if ( ! length(initvec)) return(NULL)
  refit_info <- control[["refit"]]
  lowerb <- unlist(LowUp$lower)
  upperb <- unlist(LowUp$upper) 
  Optimizer <- control[["optimizer"]] ## consistent with control.corrHLfit
  if (is.null(Optimizer)) {
    if (use_optimizer1D <- (length(initvec)==1L)) {
      uuinit <- unlist(user_init_optim)
      uuinit_not_nan <- uuinit[ ! is.nan(uuinit)]
      use_optimizer1D <- (! length(uuinit_not_nan))
    }
    if (use_optimizer1D) { 
      Optimizer <- spaMM.getOption("optimizer1D")
      if (Optimizer=="default") Optimizer <- "optimize" ## no control of initial value (but it _is_ faster that the other optimizers)
    } else {
      Optimizer <- spaMM.getOption("optimizer")
      if (Optimizer=="default") Optimizer <- "nloptr" ## old default; but new spaMM.getOption() default is ...?
    }
  }
  if (Optimizer=="optimize") {
    if (length(unlist(user_init_optim))) {
      message(paste("1D optimization by optimize(): user-provided initial value is ignored.\n",
                    "Change spaMM option 'optimizer1D' to override this.")
              )
    }
    if (is.character(HLcallfn.obj)) HLcallfn.obj <- eval(as.name(HLcallfn.obj)) # ## do.call("optimize", c(<list>, list(fn = objfn))) does not work with a char string
    locarglist <- c(anyHLCor_obj_args,list(f=HLcallfn.obj, interval=c(lowerb,upperb), maximum=TRUE))
    tol <- control[["optimize"]]$tol
    if (is.null(tol)) tol <- spaMM.getOption("optimize_tol")
    locarglist$tol <- tol
    optr <- do.call("optimize",locarglist) ## MAXimization of +logL <- HLcallfn.obj(...)
    optPars <- relist(optr$m,init.optim)
  } else if (Optimizer=="bobyqa") { ## May more narrowly approach lowerb and upperb, ~> longer computation times
    bobyqa_controls <- spaMM.getOption("bobyqa")
    bobyqa_controls$npt <- 2*length(initvec)+1
    bobyqa_controls$rhobeg <- min(0.95, 0.2*min(upperb-lowerb))
    bobyqa_controls$rhoend <- min(bobyqa_controls$rhobeg/1e6, 1e-6)
    bobyqa_controls[names(control$bobyqa)] <- control$bobyqa ## Overwrite defaults with any element of $bobyqa
    if (is.null(bobyqa_controls$maxfun)) bobyqa_controls$maxfun <- max(2*( eval(.spaMM.data$options$maxeval)),1+10*length(initvec)^2) # bobyqa will complain if not > second value
    optr <- .optim_by_bobyqa(lowerb, upperb, initvec, objfn_locoptim, anyHLCor_obj_args, HLcallfn.obj, init.optim, 
                             bobyqa_controls) ## does not use gradients
    optPars <- relist(optr$par,init.optim)
  } else if (Optimizer=="L-BFGS-B") {
    parscale <- (upperb-lowerb) 
    # nlminb code removed 11/2016 ## maxcorners code removed in v2.1.94
    control_optim <- list(parscale=parscale,factr=1e9) ## factr was the stricter 1e8 up to 23/01/13
    control_optim[names(control[["optim"]]$control)] <- control[["optim"]]$control ## ...which may be overwritten 
    optr <- optim(par=initvec,fn=objfn_locoptim,lower=lowerb,upper=upperb,control=control_optim,method="L-BFGS-B",
                  gr=grad_locoptim, 
                  anyHLCor_obj_args=anyHLCor_obj_args, HLcallfn.obj=HLcallfn.obj) ## optimize HLCor.obj()'s 'objective'
    optPars <- relist(optr$par,init.optim)
    ## full optr is big. We take out the two items that contribute much to saveSize:
    optr$eval_f <- NULL
    optr$nloptr_environment <- NULL
    attr(optPars,"optr") <- optr  
    attr(optPars,"method") <- "optim"  
  } else { ##"nloptr"
    nloptr_controls <- spaMM.getOption("nloptr")
    #if (spaMM.getOption("TRY_R_new")) nloptr_controls$xtol_rel <- nloptr_controls$xtol_rel * length(lowerb) # weird; why ?
    nloptr_controls[names(control$nloptr)] <- control$nloptr ## Overwrite defaults with any element of $nloptr
    if (is.null(nloptr_controls$maxeval)) nloptr_controls$maxeval <- eval(.spaMM.data$options$maxeval)
    if (is.null(nloptr_controls$xtol_abs)) nloptr_controls$xtol_abs <- eval(.spaMM.data$options$xtol_abs, 
                                                                            list(LowUp=LowUp, rC_transf=.spaMM.data$options$rC_transf))
    optr <- .optim_by_nloptr(lowerb=lowerb, upperb=upperb, initvec=initvec, objfn_locoptim=objfn_locoptim, 
                             anyHLCor_obj_args=anyHLCor_obj_args, HLcallfn.obj=HLcallfn.obj, init.optim=init.optim, 
                             nloptr_controls=nloptr_controls, grad_locoptim = grad_locoptim) 
    optPars <- relist(optr$solution,init.optim)
    if (anyNA(refit_info)) refit_info <- (nloptr_controls$xtol_rel > (5e-6 + 1e-8)) ## FIXME not documented (& anyNA to handle NULL)
    ## full optr is big. We take out the two items that contribute much to saveSize:
    optr$eval_f <- NULL
    optr$nloptr_environment <- NULL
  }
  optPars <- structure(optPars,method=Optimizer,optr=optr,
                       refit_info=refit_info) ## refit_info is control[["refit"]] if code follows the doc (but there is an undocumented 'FIXME')
  return(optPars)
}



