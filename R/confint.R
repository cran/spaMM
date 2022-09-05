.boot_single_par <- function(boot.ci_args, t0, ts, verbose) {
  boot.ci_args$t0 <- t0
  boot.ci_args$t <- ts
  boot.ci_args$boot.out <- list(R = length(ts), sim="parametric")
  resu <- do.call("boot.ci", boot.ci_args)
  ## resu$call is shown and this may be ugly bc of the long t vector. print.bootci() uses dput(), which has no generic for that. We will wrap the 
  ## print.bootci() call in a print.bootci4call() that locally alterns the $call for nicer printing.
  class(resu) <- c("bootci4print", class(resu))
  if (verbose) print(resu)
  resu
}


.confint_boot <- function(boot_args, object, expr_t, t_fn=NULL, parm, boot.ci, level, verbose, ...) {
  spaMM_boot_args <- intersect(names(boot_args),names(formals(spaMM_boot))) # possible conflict between boot.ci 'type' arg and spaMM_boot 'type' arg
  spaMM_boot_args <- boot_args[spaMM_boot_args]
  spaMM_boot_args$object <- object
  if (is.null(spaMM_boot_args$type)) {
    spaMM_boot_args$type <- "residual"
  } else if (length(intersect(spaMM_boot_args$type,c("basic","perc","norm")))) {
    warning(
      paste0("Hmmm. It looks like you are using 'boot_args$type' to pass\n", 
             "the boot.ci() 'type' argument. Use 'boot_args$ci_type' for that purpose.\n",
             "'boot_args$type' is for passing the spaMM_boot() 'type' argument.\n",
             "[see help(\"confint.HLfit\") for further details].")
      ,immediate. = TRUE)
  }
  #
  if (is.null(boot_args$ci_type)) {
    boot_args$type <- c("basic","perc","norm")
  } else {
    if (length(setdiff(boot_args$ci_type,c("basic","perc","norm")))) 
      stop("Check the 'boot_args$ci_type' argument: only \"basic\", \"perc\", and \"norm\" CI types can be computed.")
    boot_args$type <- boot_args$ci_type
    boot_args$ci_type <- NULL
  }
  boot.ci_args <- intersect(names(boot_args),names(formals(boot.ci)))
  boot.ci_args <- boot_args[boot.ci_args]
  boot.ci_args$conf <- level
  if (is.null(t_fn)) {
    spaMM_boot_args$simuland <- function(y, ...) {
      upd <- update_resp(object, newresp=y)
      eval(expr_t, list(hlfit=upd))
    }
    t0 <- eval(expr_t, list(hlfit=object))
  } else {
    spaMM_boot_args$simuland <- function(y, ...) {
      upd <- update_resp(object, newresp=y)
      t_fn(upd)
    }
    t0 <- t_fn(object, ...)
  }
  ts <- drop(do.call(spaMM_boot,spaMM_boot_args, ...)[["bootreps"]])
  hasnorm <- "norm" %in% boot.ci_args$type
  hasperc <- "perc" %in% boot.ci_args$type
  hasbasic <- "basic" %in% boot.ci_args$type
  np <- NCOL(ts)
  template <- matrix(NA,ncol=2,nrow=np)
  marg <- (1-level)*50
  colnames(template) <- paste(c(marg, 100-marg),"%")
  if (is.character(parm)) rownames(template) <- parm #colnames(ts)   # parm is not correct if parm is a function etc. # in which case the names may remain NULL
  
  tl <- list()
  if (hasnorm) tl$normal <- template
  if (hasperc) tl$percent <- template
  if (hasbasic) tl$basic <- template
  if (np>1L) {
    resu <- vector("list",np)
    for (colit in seq_len(np)) {
      resu[[colit]] <- .boot_single_par(boot.ci_args, t0=t0[colit], ts=ts[,colit], verbose=verbose)
      if (hasnorm) tl$normal[colit,] <- tail(resu[[colit]]$normal[1,],n=2L)
      if (hasperc) tl$percent[colit,] <- tail(resu[[colit]]$percent[1,],n=2L)
      if (hasbasic) tl$basic[colit,] <- tail(resu[[colit]]$basic[1,],n=2L)
    }
  } else {
    resu <- .boot_single_par(boot.ci_args, t0=t0, ts=ts, verbose=verbose)
    if (hasnorm) tl$normal[1,] <- tail(resu$normal[1,],n=2L)
    if (hasperc) tl$percent[1,] <- tail(resu$percent[1,],n=2L)
    if (hasbasic) tl$basic[1,] <- tail(resu$basic[1,],n=2L)
  }
  attr(resu,"table") <- tl # a list with elements 'normal', 'percent' 'basic', each 
  return(resu)
}

.confint_LRT <- function(level, parm, object, verbose) {
  if ((np <- length(parm))>1L) {
    resu <- vector("list",np)
    lower <- upper <- numeric(np)
    for (colit in seq_len(np)) {
      resu[[colit]] <- .confint_LRT_single_par(level=level, parm=parm[colit], object=object, verbose=verbose)
      lower[colit] <- resu[[colit]]$interval[[1]]
      upper[colit] <- resu[[colit]]$interval[[2]]
    }
    names(resu) <- parm
    table. <- cbind(lower,upper)
  } else {
    resu <- .confint_LRT_single_par(level=level, parm=parm, object=object, verbose=verbose)
    table. <- resu$interval
    dim(table.) <- c(1L,2L)
  }
  rownames(table.) <- parm
  marg <- (1-level)*50
  colnames(table.) <- paste(c(marg, 100-marg),"%")
  attr(resu,"table") <- table.
  return(resu)
}

.confint_LRT_single_par <- function(level, parm, object, verbose) {
  dlogL <- qchisq(level,df=1)/2
  znorm <- qnorm((1+level)/2)
  if (is.character(parm)) {
    parmcol <- which(names(object$fixef)==parm)
    if (length(parmcol)==0L) stop("Parameter not in the model")
    attr(parm,"col") <- parmcol 
  } else {
    parmcol <- parm
    if (parm > length(object$fixef)) stop("'parm' not compatible with # of fixed-effects coefficients")
    parm <- names(object$fixef)[parmcol]
    attr(parm,"col") <- parmcol
  }
  llc <- getCall(object)
  fnname <-.get_bare_fnname.HLfit(object, call.=llc)
  lc <- switch(fnname,
               "corrHLfit" = get_HLCorcall(object,fixed=llc$ranFix),
               "fitme" = get_HLCorcall(object,fixed=llc$fixed), # HLfit or HLCor call
               "HLCor" = get_HLCorcall(object,fixed=llc$ranPars),
               "HLfit" = get_HLCorcall(object,fixed=llc$ranFix),
               "fitmv" = get_HLCorcall(object,fixed=llc$fixed), 
               stop("Unhandled getCall(object) in confint.HLfit()")
  ) # uniformly calls get_HLCorcall() -> .preprocess() so that $processed is always available, 
    # .makeLowerUpper() is called. THe optimInfo$LUarglist is used only later 
    # For fitmv, .makeLowerUpper() is called several times.
  HL <- object$HL
  fixeflik <- switch(paste(HL[1L]),
                "0"=if (object$models[["eta"]]=="etaGLM") "p_v" else "hlik", # if the user fitted a GLM by PQL/L there is no hlik 
                "1"="p_v",
                stop(paste("confint does not yet handle HLmethod",paste(HL,collapse=" "),
                           "(or ",c(llc$method,llc$HLmethod),").",sep=" ")))
  beta_cov <- .get_beta_cov_any_version(object)
  beta_se <- sqrt(diag(x=beta_cov))[parm]
  asympto_abs_Dparm <- znorm* beta_se
  X.pv <- lc$processed$AUGI0_ZX$X.pv
  X_is_scaled <- ( ! is.null(attr(X.pv,"scaled:scale")))
  #
  # FIXME test processed nature to prevent multinom stuff here
  warnlik <- "p_v" # bc calling on an REML fit is poor anyway.
  likfns <- unique(c(fixeflik,warnlik))
  fixeflik <- unlist(object$APHLs[fixeflik]) # named
  warnlik <- unlist(object$APHLs[warnlik])
  intervalinfo <- list(fixeflik=fixeflik,
                       warnlik=warnlik, 
                       likfns=likfns,
                       targetlik=fixeflik-dlogL,
                       parm=parm, # name, vs $MLparm: ML value
                       parmcol_X=parmcol,
                       parmcol_ZX=length(object$lambda.object$lambda_est)+parmcol, 
                       no_phi_pred= ! ("phiHGLM" %in% object$models[["phi"]]), # %in% for mv
                       asympto_abs_Dparm=asympto_abs_Dparm)
  oldopt <- .options.processed(lc$processed, intervalInfo=intervalinfo, augZXy_cond=FALSE)
  if (X_is_scaled) {
    lc$processed$intervalInfo$MLparm <- .scale(beta=object$fixef,X=X.pv)[parm]
  } else lc$processed$intervalInfo$MLparm <- object$fixef[parm]   
  #
  .assignWrapper(lc$processed,"LevenbergM['force'] <- FALSE") ## inhibits LevM for confint 
  optimInfo <- attr(object,"optimInfo") ## may be NULL
  trTemplate <- optimInfo$`optim.pars` ## may be NULL if optimInfo is NULL or if  optimInfo is not NULL but no par as outer optimized
  attr(trTemplate,"optr") <- NULL 
  attr(trTemplate,"method") <- NULL 
  if ( ! is.null(augZXy_phi_est <- optimInfo$augZXy_phi_est)) { ## augZXy not used for confint but may have been used in the original fit.
    lambda <- .dispInv(trTemplate$trLambda) ## assuming $trLambda exists...
    lambda <- lambda * augZXy_phi_est 
    trTemplate$trLambda <- .dispFn(lambda) 
  }
  if ( ! is.null(trTemplate)) {
    olc <- lc ## olc is working copy
    LUarglist <- optimInfo$LUarglist
    attr(trTemplate,"method") <- NULL
    if (paste(lc[[1]])=="HLCor") {
      ## good starting values are important... important to use canonizeRanPars as in HLCor
      attr(trTemplate,"moreargs") <- LUarglist$moreargs
      ## locoptim expects a fn with first arg ranefParsVec
      ## ca serait mieux de pas avoir de contrainte la dessus et de pvr nommer l'arg trParsVec
      ## bc HLCor call uses transformed scale for ranPars
      objfn <- function(ranefParsVec, anyHLCor_obj_args=NULL, HLcallfn.obj=NULL) { 
        ranefParsList <- relist(ranefParsVec,trTemplate)
        olc$ranPars <- structure(.modify_list(olc$ranPars,ranefParsList)) ## replaces ! some elements and keeps the "type" !
        locfit <- eval(as.call(olc)) ## HLCor call with given ranefParsVec
        resu <- (posforminimiz)* locfit$fixef[parm]
        attr(resu,"info") <- locfit$APHLs$p_v 
        ## attribute lost by optim but otherwise useful for debugging 
        return(resu) ## return value to be optimized is a parameter value, not a likelihood
      }
    } else if (paste(lc[[1]])=="HLfit") {
      objfn <- function(ranefParsVec, anyHLCor_obj_args=NULL, HLcallfn.obj=NULL) { 
        ranefParsList <- relist(ranefParsVec,trTemplate)
        olc$ranFix <- structure(.modify_list(olc$ranFix,ranefParsList)) ## replaces ! some elements and keeps the "type" !
        locfit <- eval(as.call(olc)) ## HLfit call with given ranefParsVec
        resu <- (posforminimiz)*locfit$fixef[parm]
        attr(resu,"info") <- locfit$APHLs$p_v 
        ## attribute lost by optim but otherwise useful for debugging 
        return(resu) ## return value to be optimized is a parameter value, not a likelihood
      }
    }
    rC_transf <- .spaMM.data$options$rC_transf
    LUarglist$canon.init <- .canonizeRanPars(ranPars=trTemplate,
                                             corr_info=.get_from_ranef_info(object), 
                                             checkComplete=FALSE, rC_transf=rC_transf)
    LowUp <- do.call(.makeLowerUpper,LUarglist)
    optim_bound_over_nuisance_pars <- function(posforminimiz) { ## optimize the CI bound and returns the parameters that optimize this bound
      user_init_optim <- switch(fnname,
                                "corrHLfit" = llc[["init.corrHLfit"]],
                                "fitme" = llc[["init"]], 
                                NULL)
      init <- unlist(trTemplate)
      if (paste(lc[[1]])=="HLCor") { HLcallfn_obj <- "HLCor.obj" } else HLcallfn_obj <- "HLfit.obj"
      .assignWrapper(anyObjfnCall.args$processed,
                     paste0("return_only <- \"confint_bound\""))
      optr <- .new_locoptim(init.optim=trTemplate,LowUp=LowUp,
                            objfn_locoptim=objfn, # uses posforminimiz in its definition 
                            HLcallfn.obj=HLcallfn_obj,
                            user_init_optim=user_init_optim,
                            anyHLCor_obj_args=anyObjfnCall.args,
                            control=list(optimizer=spaMM.getOption("optimizer")), # important to avoid use of optimize()
                            verbose=FALSE) 
      .assignWrapper(anyObjfnCall.args$processed,
                     paste0("return_only <- NULL"))
      # We need an optimizer with control of the initial value (hence not optimize());
      # otherwise the optimizer may never find a value of the nuisance pars that results in a focal parameter value 
      #  with high enough lik. In that case the returned value of the focal parameter is the ML estimate and the interval reduces to the ML estimate.
      return(optr) ## the bound, relist()'ed according to trTemplate
    }
  }
  ## lowerfit
  fac <- 1L 
  warnori <- options(warn=-1)
  prevmsglength <- 0L
  while(fac < 1e6) {
    init_beta <- object$fixef-asympto_abs_Dparm/fac
    if (X_is_scaled) init_beta <- .scale(beta=init_beta,X=X.pv) ## using locally saved X.pv
    lc$processed$intervalInfo$init <- init_beta[parm]
    lc$processed$intervalInfo$init_v_h <- object$v_h
    if (! is.null(trTemplate)) {
      anyObjfnCall.args <- as.list(lc[-1L]) ## includes processed, ranPars, controlS.dist, control.HLfit...
      anyObjfnCall.args$skeleton <- trTemplate
      # The objective function 'objfn' returns the confint bound given the corr pars. Thus locoptim maximizes the confint bound over the the corr pars
      olc <- lc ## that's olc that is used in the objective fn !
      posforminimiz <- 1 ## defined in the envir where objfn is defined... (bad style)
      bound <- optim_bound_over_nuisance_pars()
      if (paste(lc[[1]])=="HLCor") {
        olc$ranPars <- structure(.modify_list(olc$ranPars,bound)) ## replaces ! some elements and keeps the "type" (lazyness)!
      } else {
        olc$ranFix <- structure(.modify_list(olc$ranFix,bound)) ## replaces ! some elements and keeps the "type" !
      }
      ## recover fit for optimized params (must use call with intervalInfo and LevenbergM=FALSE)
      lowerfit <- eval(as.call(olc)) ## full HLfit objectobject
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
    init_beta <- object$fixef+asympto_abs_Dparm/fac
    if (X_is_scaled) {
      init_beta <- .scale(beta=init_beta,X=X.pv)
    }
    lc$processed$intervalInfo$init <- init_beta[parm]
    if (! is.null(trTemplate)) {
      olc <- lc
      posforminimiz <- -1 ## maximization
      bound <- optim_bound_over_nuisance_pars()
      if (paste(lc[[1]])=="HLCor") {
        olc$ranPars <- structure(.modify_list(olc$ranPars,bound)) ## replaces ! some elements and keeps the "type" !
      } else {
        olc$ranFix <- structure(.modify_list(olc$ranFix,bound)) ## replaces ! some elements and keeps the "type" !
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
  if (verbose && prevmsglength) cat("\n")
  options(warnori)
  .options.processed(lc$processed, oldopt)
  interval <- c(lowerfit$fixef[parm],upperfit$fixef[parm])
  names(interval) <- paste(c("lower","upper"),parm)
  if (verbose) print(interval)
  resu <- list(lowerfit=lowerfit,upperfit=upperfit,interval=interval)
  if (! is.null(resu$confint_best_fit <- lc$processed$envir$confint_best)) {
    locmess <- paste("Element 'confint_best_fit' of the return object contain information about a possible better fit to the data",
                     "\nYou can for example refit the model using the parameter values shown in that element as initial values.")
    message(locmess)
  }
  return(resu)
}


confint.HLfit <- function(object, parm, level=0.95, verbose=TRUE, 
                          boot_args=NULL, format="default",...) {
   
  if (is.character(parm)) {
    if (is.list(boot_args)) expr_t <- substitute(fixef(hlfit)[parm], list(parm=parm))
    t_fn <- NULL
  } else if (is.function(parm)) {
    t_fn <- parm
  } else {
    expr_t <- parm
    t_fn <- NULL
    #if (is.null(boot_args)) boot_args <- list(nsim=999L) fail is the example from the doc...
  }
  if (  is.list(boot_args)) {
    boot_res <- .confint_boot(boot_args, object, expr_t, t_fn, parm, boot.ci, level, verbose = verbose) 
    if (format=="stats") {
      attr(boot_res,"table")
    } else invisible(boot_res)
  } else {
    if ( .REMLmess(object,return_message=FALSE)) {
      warning("REML fits are not quite suitable for computing intervals for fixed effects.")
    }
    resu <- .confint_LRT(level, parm, object, verbose)
    if (format=="stats") {
      attr(resu,"table")
    } else invisible(resu)
  }
}
