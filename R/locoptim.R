## .new_locoptim() function wraps optim, optimize and NLOPT_LN_BOBYQA
# It uses optimize(HLcallfn.obj .. maximum=TRUE) or OTHERWISE the .objfn_locoptim() wrapper to maximize likelihood:
#  this is messy to control hence the ad_hoc_fn wrapper in confint.HLfit
# The first arg of the objective functions must be ranefParsVec
.objfn_locoptim <- function(x, anyHLCor_obj_args, HLcallfn.obj) { ## the (more or less) default value of .new_locoptim <- function( .. objfn_locoptim .. ) 
  anyHLCor_obj_args$ranefParsVec <- x 
  return( - do.call(HLcallfn.obj, anyHLCor_obj_args))
}

.optim_by_nloptr <- function(lowerb, upperb, initvec, objfn_locoptim, local_control, grad_locoptim=NULL, LowUp, ...) {
  nloptr_controls <- .get_nloptr_controls(init=initvec, LowUp=LowUp)
  nloptr_controls[names(local_control)] <- local_control ## Overwrite defaults with any element of $nloptr
  ## this is also called if length(lower)=0 by  (SEM or not) and optPars is then null 
  optr <- nloptr::nloptr(x0=initvec, eval_f=objfn_locoptim,
                         eval_grad_f=grad_locoptim, # ignored with NLOPT_LN_BOBYQA
                         lb=lowerb,ub=upperb, opts=nloptr_controls, ...)
  while (optr$status==5L) { ## optr$status=5 => termination bc maxeval has been reached 
    # met status=4: nloptr message in normal termination due toxtol_rel, but is this true ?
    message("maxeval reached in nloptr(); nloptr() called again until apparent convergence of objective.") 
    prevlik <- optr$objective
    reinit <- pmax(lowerb,pmin(upperb,optr$solution))
    optr <- nloptr::nloptr(x0=reinit, eval_f=objfn_locoptim,
                           eval_grad_f=grad_locoptim, # ignored with NLOPT_LN_BOBYQA
                           lb=lowerb,ub=upperb, opts=nloptr_controls, ...)
    loc_ftol <- max(1e-8, optr$options$ftol_abs)
    if (- optr$objective < - prevlik+loc_ftol) break ## no progress in <= maxeval iterations
  }
  return(optr)
}

#
.optim_by_bobyqa <- function(lowerb, upperb, initvec, objfn_locoptim, local_control, adjust_init=list(), ...) {
  bobyqa_controls <- .get_bobyqa_controls(init=initvec, upper=upperb, lower=lowerb)
  bobyqa_controls[names(local_control)] <- local_control ## Overwrite defaults with any element of $bobyqa
  bobyqa_margin <- .spaMM.data$options$bobyqa_margin
  margin <- (upperb-lowerb)*bobyqa_margin 
  margin <- pmin(bobyqa_margin,margin) # handles infinite ranges (but not only)
  init <- pmax(lowerb+margin,pmin(upperb-margin,initvec))
  # And this is a more substantial adjustment at the margin
  if ( ! is.null(adjust_init$lower)) init <- pmax(adjust_init$lower, init)
  if ( ! is.null(adjust_init$upper)) init <- pmin(adjust_init$upper, init)
  optr <- bobyqa(par=init,fn=objfn_locoptim,lower=lowerb,upper=upperb,control=bobyqa_controls, ...)
  while(optr$ierr==1L) { #maximum number of function evaluations exceeded
    message("maxeval reached in bobyqa(); bobyqa() called again until apparent convergence of objective.") 
    prevmlik <- optr$fval
    reinit <- pmax(lowerb,pmin(upperb,optr$par))
    optr <- bobyqa(par=reinit,fn=objfn_locoptim,lower=lowerb,upper=upperb,control=bobyqa_controls, ...)
    if (optr$fval > prevmlik-1e-8) break ## not enough progress in <= maxeval iterations
  }
  optr$value <- - optr$fval
  return(optr)
}

.xtol_abs_fn <- function(LowUp, # must be a structured list, not simply a list of two vectors, for it to have some effect.
                         factors=.spaMM.data$options$xtol_abs_factors, rC_transf=.spaMM.data$options$rC_transf) {
  parnames <- names(LowUp$lower)
  rng <- unlist(LowUp$upper, use.names = FALSE)-unlist(LowUp$lower, use.names = FALSE)
  rng_finite <- is.finite(rng)
  if ("trRanCoefs" %in% parnames) {
    xtol_abs <- .relist_rep(NA,LowUp$lower)
    for (st in parnames) {
      if (st=="trRanCoefs") {
        trRanCoefs <- LowUp$lower$trRanCoefs
        for (rc in names(trRanCoefs)) {
          len <- length(trRanCoefs[[rc]])
          Xi_ncol <- floor(sqrt(len*2))
          # if (rC_transf=="chol") {
          #   xtol_abs[[st]] <- rep(1e-12,len) # note that order of elements is that of upper.tri 
          # } else 
          xtol_abs$trRanCoefs[[rc]] <- c(rep(factors["rcLam"],Xi_ncol),rep(factors["rcCor"],len-Xi_ncol)) # "sph" etc
        }
      } else {xtol_abs[[st]] <- rep(factors["others"],length(unlist(LowUp$lower[[st]], use.names = FALSE)))}
    }
    xtol_abs <- unlist(xtol_abs, use.names = FALSE)
  } else {
    # For a long time factors$abs has been 1e-7 and there was no scaling by finite range
    # With sclaing, 1e-7 is too large.
    xtol_abs <- rep(factors[["abs"]],length(unlist(LowUp$lower, recursive=TRUE, use.names = FALSE)))
  }
  xtol_abs[rng_finite] <- xtol_abs[rng_finite] * rng[rng_finite]
  return(xtol_abs)
}

# This function wraps TMB::MakeADFun():  It converts a 'processed' to input for this function, and return the 
# AD structure suitable for inputs to nlminb().
# The spaMM_ADFun DLL must be loaded as in devel/TMB/generic: this DLL may have very limited functionality,
# being first developed for the diagnosis [see devel/TMB/diagnosis_poly6.R] 
# of [discrepancies between packages on poly(cbind(age, parity), 6) fit in test-back-compat.R]
.wrap_MakeADFun <- function(processed, init.optim, DLL) {
  y <- processed$y
  if (is.null(scX <- environment(processed$X_off_fn)$X_off)) { # this scX is available if an init beta was given; then $AUGI0_ZX$X.pv has zero cols
    stop("won't work bc processed init, upper and lower do not have beta.")
    scX <- processed$AUGI0_ZX$X.pv
  }
  ZAlist <- processed$ZAlist
  nrand <- length(ZAlist)
  
  parlist <- list(b=.scale(scX, 0))
  
  Z <- vector("list", nrand)
  if (nrand) {
    for (rd in seq_len(nrand)) {
      Zrd <- ZAlist[[rd]]
      attr(Zrd,"is_incid") <- NULL
      attr(Zrd,"namesTerm") <- NULL
      Z[[rd]] <- Zrd
    }
    parlist$trLambda <- c(init.optim$trLambda) # c() dropping attributes 
    cum_n_u_h <- processed$cum_n_u_h 
    parlist$u <- rep(0,tail(cum_n_u_h,1L))
  }
  
  famfam <- processed$family$family
  if (famfam %in% c("negbin1","negbin2")) {
    parlist$trShape <- init.optim$trNB_shape
  } else stop("family not yet handled.")
  
  template_order <- c("b","trLambda","trShape","u")
  parlist <- parlist[intersect(template_order,names(parlist))]
  MakeADFun <- get("MakeADFun", envir = asNamespace("TMB") ) 
  adfun <- MakeADFun(data = list(n=length(y), y=y, scX=scX, ZAlist=Z, nrand=nrand),
                     parameters=parlist,
                     DLL = DLL,
                     random = "u")
  adfun
}


# returns optPars which is a list given by relist(.,init.optim), with attributes the optimMethod and (+:- raw) optr 
.new_locoptim <- function(init.optim, LowUp, control, objfn_locoptim, 
                          anyHLCor_obj_args, HLcallfn.obj="HLCor.obj", 
                          user_init_optim, ## only purpose is to make sure that if the user provides an explicit init in 1D, optimize() is not used.
                          grad_locoptim=NULL,
                          verbose,
                          ADFun=anyHLCor_obj_args$processed$ADFun # 
) {
  initvec <- unlist(init.optim) 
  if ( ! length(initvec)) return(NULL)
  refit_info <- control[["refit"]]
  lowerb <- unlist(LowUp$lower)
  upperb <- unlist(LowUp$upper) 
  Optimizer <- control[["optimizer"]] ## consistent with control.corrHLfit
  # If user provides an explicit init in 1D, optimize() is not used:
  if (is.null(Optimizer)) {
    if ( ! is.null(ADFun)) {
      Optimizer <- "nlminb"
      use_optimizer1D <- FALSE
      if ( ! is.list(ADFun)) { # presumably user input ADFun=TRUE or =<a DLL name>; TMB::makeADFun returns a list
        if ( ! is.character(ADFun)) ADFun <- "spaMM_ADFun"
        ADFun <- .wrap_MakeADFun(anyHLCor_obj_args$processed, init.optim=init.optim, DLL=ADFun)
      }
    } else {
      if (use_optimizer1D <- (length(initvec)==1L)) {
        uuinit <- unlist(user_init_optim)
        uuinit_not_nan <- uuinit[ ! is.nan(uuinit)]
        use_optimizer1D <- (! length(uuinit_not_nan))
      }
      if (use_optimizer1D) { 
        Optimizer <- spaMM.getOption("optimizer1D")
        ## if (Optimizer=="default") Optimizer <- "optimize" ## no control of initial value (but it _is_ faster that the other optimizers)
      } else {
        Optimizer <- spaMM.getOption("optimizer")
        ## if (Optimizer=="default") Optimizer <- ".safe_opt" # in case "default" would be used; but the latter case looks obsolete.
      }
    }
  }
  
  if (is.function(Optimizer)) { # user provided optimizer; private for devel purposes: cf # cf devel/TMB/spatial_poisson/spatial_TMB_python_devel.R
    # Then this function should have the same interface as .optim_by_nloptr()
    user_def_optimizer <- Optimizer
    Optimizer <- "user-defined" # used to build returned structure, whose API assumes it compares to a character string.
    optr <- user_def_optimizer(lowerb=lowerb, upperb=upperb, initvec=initvec, objfn_locoptim=objfn_locoptim, 
                                      local_control=control[["nloptr"]], grad_locoptim = grad_locoptim, LowUp=LowUp,
                                      anyHLCor_obj_args=anyHLCor_obj_args, HLcallfn.obj=HLcallfn.obj) 
    optPars <- relist(optr$solution,init.optim)
    attr(optPars,"optr") <- optr  
  } else if (Optimizer=="optimize") {
    # since explicit init by user is heeded, the following message is only helpful to me in a tracing session...
    if (verbose) message(paste("1D optimization by optimize(): spaMM's *default* initial value is ignored.\n",
                  "Provide explicit initial value, or change spaMM option 'optimizer1D' for initial value to be taken into account."))
    if (is.character(HLcallfn.obj)) HLcallfn.obj <- eval(as.name(HLcallfn.obj)) # ## do.call("optimize", c(<list>, list(fn = objfn))) does not work with a char string
    locarglist <- c(anyHLCor_obj_args,list(f=HLcallfn.obj, interval=c(lowerb,upperb), maximum=TRUE))
    tol <- control[["optimize"]]$tol
    if (is.null(tol)) tol <- spaMM.getOption("optimize_tol")
    locarglist$tol <- tol
    optr <- do.call("optimize",locarglist) ## MAXimization of +logL <- HLcallfn.obj(...)
    optPars <- relist(optr$maximum,init.optim)
  } else if (Optimizer=="nloptr") { 
    optr <- .optim_by_nloptr(lowerb=lowerb, upperb=upperb, initvec=initvec, objfn_locoptim=objfn_locoptim, 
                             local_control=control[["nloptr"]], grad_locoptim = grad_locoptim, LowUp=LowUp,
                             anyHLCor_obj_args=anyHLCor_obj_args, HLcallfn.obj=HLcallfn.obj) 
    optPars <- relist(optr$solution,init.optim)
    if (anyNA(refit_info)) refit_info <- (optr$options$xtol_rel > (5e-6 + 1e-8)) ## FIXME not documented (& anyNA to handle NULL)
    ## full optr is big. We take out the two items that contribute much to saveSize:
    optr$eval_f <- NULL
    optr$nloptr_environment <- NULL
  } else if (Optimizer==".safe_opt") { ## May more narrowly approach lowerb and upperb, ~> longer computation times
    optr <- .safe_opt(init=initvec, lower=lowerb, upper=upperb, 
                      objfn=objfn_locoptim, # minimization of -logL
                      verbose=max(0L,verbose-1L), 
                      anyHLCor_obj_args=anyHLCor_obj_args, HLcallfn.obj=HLcallfn.obj , LowUp=LowUp
                      ) ## does not use gradients
    optPars <- relist(optr$solution,init.optim)
  } else if (Optimizer=="bobyqa") { ## May more narrowly approach lowerb and upperb, ~> longer computation times
    optr <- .optim_by_bobyqa(lowerb, upperb, initvec, objfn_locoptim,
                             local_control=control[["bobyqa"]], anyHLCor_obj_args=anyHLCor_obj_args, 
                             HLcallfn.obj=HLcallfn.obj) ## does not use gradients
    optPars <- relist(optr$par,init.optim)
    optr$objective <- optr$fval # for easy tests on the results, e.g. test-ranCoefs.R
  } else if (Optimizer=="nlminb") { 
    nlminb_controls <- .get_nlminb_controls(init=initvec, upper=upperb, lower=lowerb)
    local_control <- control[["nlminb"]]
    nlminb_controls[names(local_control)] <- local_control ## Overwrite defaults with any element of $nlminb
    if ( ! is.null(ADFun)) { 
      # cf devel/TMB/spatial_poisson/spatial_TMB_python_devel.R 
      #    devel/TMB/generic.R
      #    devel/TMB/diagnosis_poly6.R
      objfn <- function(x, ...) ADFun$fn(x) # -logL
      gr <- function(x, ...) ADFun$gr(x) # gradient of -logL
      optr <- stats::nlminb(initvec, objfn, lower=lowerb,upper=upperb, gradient=gr,
                            anyHLCor_obj_args=anyHLCor_obj_args, HLcallfn.obj=HLcallfn.obj, 
                            control=nlminb_controls)
    } else optr <- stats::nlminb(initvec, objfn_locoptim, lower=lowerb,upper=upperb, 
                   anyHLCor_obj_args=anyHLCor_obj_args, HLcallfn.obj=HLcallfn.obj, 
                   control=nlminb_controls)
    if ( ! optr$iterations < nlminb_controls$iter.max) warning(paste0("nlminb() reached control$iter.max=", 
                                                                      nlminb_controls$iter.max))
    if ( ! sum(optr$evaluations) < nlminb_controls$eval.max) warning(paste0("nlminb() reached control$eval.max=", 
                                                                            nlminb_controls$eval.max))
    ## with the same arguments one can try BB::spg() and dfoptim:hjkb, both unconvincing
    optPars <- relist(optr$par,init.optim)
  } else if (Optimizer=="L-BFGS-B") { # legacy
    parscale <- (upperb-lowerb) 
    parscale[is.infinite(parscale)] <- 2000 # ad hoc patch. Inf occurs for ranCoefs and 2000 is of the order of parscale for 'simple lambdas'
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
  } else stop("Unhandled optimizer")
  optPars <- structure(optPars,method=Optimizer,optr=optr,
                       refit_info=refit_info) ## refit_info is control[["refit"]] if code follows the doc (but there is an undocumented 'FIXME')
  return(optPars)
}



