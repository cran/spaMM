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


.confint_boot <- function(boot_args, object, expr_t, t_fn=NULL, parm, boot.ci, level, 
                          verbose, # will be ignored by update_resp(), but used by .boot_single_par()
                          ...) {
  spaMM_boot_args <- intersect(names(boot_args),names(formals(spaMM_boot))) # possible conflict between boot.ci 'type' arg and spaMM_boot 'type' arg
  spaMM_boot_args <- boot_args[spaMM_boot_args]
  spaMM_boot_args$object <- object
  if (is.null(spaMM_boot_args$type)) {
    spaMM_boot_args$type <- "marginal"
    message("Missing 'type' in 'boot_args' is set to '",
            spaMM_boot_args$type,"'\n (this default type has been changed in version 4.4.23).")
  } else if (length(intersect(spaMM_boot_args$type,c("basic","perc","norm")))) {
    warning(
      cli::format_warning(paste0("Hmmm. It looks like you are using 'boot_args$type' to pass\n", 
             "the boot.ci() 'type' argument. Use 'boot_args$ci_type' for that purpose.\n",
             "'boot_args$type' is for passing the spaMM_boot() 'type' argument.\n",
             "[see {.help [{.fun confint.HLfit}](spaMM::confint.HLfit)} for further details]."))
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
      upd <- update_resp(object, newresp=y, ...)
      eval(expr_t, list(hlfit=upd))
    }
    t0 <- eval(expr_t, list(hlfit=object))
  } else {
    spaMM_boot_args$simuland <- function(y, ...) {
      upd <- update_resp(object, newresp=y, ...)
      t_fn(upd)
    }
    t0 <- t_fn(object)
  }
  ts <- drop(do.call(spaMM_boot,spaMM_boot_args, ...)[["bootreps"]])
  hasnorm <- "norm" %in% boot.ci_args$type
  hasperc <- "perc" %in% boot.ci_args$type
  hasbasic <- "basic" %in% boot.ci_args$type
  np <- .old_NCOL(ts)
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

# This rests on .confint_LRT_single_par() which hacks the 'inner' fitting algo for fixed effects to optimize the CI bound over all other parameters 
# and returns the parameters that optimize this bound, so that no numerical profiling of the likelihood has to be performed.
.confint_LRT <- function(level, parm, object, verbose, ..., 
                         is_p4m=inherits(object,"pois4mlogit")) {
  if ((np <- length(parm))>1L) {
    resu <- vector("list",np)
    lower <- upper <- numeric(np)
    for (colit in seq_len(np)) {
      if (is_p4m) {
        resu[[colit]] <- .confint_LRT_single_par_p4m(level=level, parm=parm[colit], 
                                                     object=object, verbose=verbose, ...)
      } else resu[[colit]] <- .confint_LRT_single_par(level=level, parm=parm[colit], object=object, verbose=verbose)
      lower[colit] <- resu[[colit]]$interval[[1]]
      upper[colit] <- resu[[colit]]$interval[[2]]
    }
    names(resu) <- parm
    table. <- cbind(lower,upper)
  } else {
    if (is_p4m) {
      resu <- .confint_LRT_single_par_p4m(level=level, parm=parm, object=object, verbose=verbose, ...)
    } else resu <- .confint_LRT_single_par(level=level, parm=parm, object=object, verbose=verbose)
    table. <- resu$interval
    dim(table.) <- c(1L,2L)
  }
  rownames(table.) <- parm
  marg <- (1-level)*50
  oldopt <- options(OutDec = ".")
  colnames(table.) <- paste(c(marg, 100-marg),"%")
  options(oldopt)
  attr(resu,"table") <- table.
  return(resu)
}

.old_hack_trTemplate <- function(optimInfo, object) {
  trTemplate <- optimInfo$`optim.pars` ## may be NULL if optimInfo is NULL or if  optimInfo is not NULL but no par as outer optimized
  # For partial ranCoefs, optimInfo$`optim.pars` contains a full vector 
  # when HLfit_body() -> .canonizeRanPars() is reached, there must be ranCoefs with the constraints and 
  # trRancoefs with a fully varying vector
  if ( ! is.null(trTemplate)) { # ... not from HLfit()...
    attr(trTemplate,"optr") <- NULL 
    attr(trTemplate,"method") <- NULL 
    if ( ! is.null(augZXy_phi_est <- optimInfo$augZXy_phi_est )) {## augZXy not used for confint but may have been used in the original fit.
      warning("confint() called on this LMM fit obtained with spaMM < 4.1.42 may be unreliable.\n It would be safer to refit the model.",
              immediate. = TRUE)
      fittedpars <- .get_fittedPars(object, partial_rC="rm", phiPars=FALSE, phifits=FALSE, verbose=FALSE)
      if ( ! is.null(trTemplate$trLambda)) { ## could be NULL for random-coef model
        trTemplate$trLambda <- .dispFn(fittedpars$lambda) 
      } 
      if ( ! is.null(rC <- fittedpars$ranCoefs)) { # partially fixed ranCoefs are not fitted by augZXy...
        for (char_rd in names(trTemplate$trRanCoefs)) {
          rC[[char_rd]] <- .ranCoefsFn(rC[[char_rd]], rC_transf = .spaMM.data$options$rC_transf)
        }
        trTemplate$trRanCoefs <- rC 
      } 
    }
  }
  trTemplate
}


.confint_LRT_single_par <- function(level, parm, object, verbose) {
  if (length(.unlist(.get_rC_inits_from_hlfit(object, type="inner")))) {
    # inner-estimation of ranCoefs is incompatible with the confint hack.
    stop("confint() attempted on fit with inner-estimated random-coefficient model.")
  }
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
  HL <- object$HL
  fixeflik <- switch(paste(HL[1L]),
                "0"=if (object$models[["eta"]]=="etaGLM") "p_v" else "hlik", # if the user fitted a GLM by PQL/L there is no hlik 
                "1"="p_v",
                stop(paste("confint does not yet handle HLmethod",paste(HL,collapse=" "),
                           "(or ",c(llc$method,llc$HLmethod),").",sep=" ")))
  beta_cov <- .get_beta_cov_any_version(object)
  beta_se <- sqrt(diag(x=beta_cov))[parm]
  asympto_abs_Dparm <- znorm* beta_se
  #
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
  ### In previous versions of this fn I hacked 'processed' after creation but nw I can avoid that
  ###  + get suitable optimInfo's trTemplate rather than need to hack the one of the original fit 
  ###         (=> was causing pbs if original augZXy phi scaling) 
  ## modif control.HLfit for intervalInfo to be taken into account eg to inhibit augZXy
  control.HLfit <- llc$control.HLfit
  control.HLfit$intervalInfo <- intervalinfo
  control.HLfit$LevenbergM <- FALSE # maybe tat could be automatic too
  # 
  fittingFunction <- .get_bare_fnname.HLfit(object) 
  lc <- .get_HLCorcall_W_init(object, llc=llc, fittingFunction=fittingFunction, control.HLfit)
  processed <- lc$processed
  X.pv <- processed$AUGI0_ZX$X.pv
  X_is_scaled <- ( ! is.null(attr(X.pv,"scaled:scale")))
  if (X_is_scaled) {
    processed$intervalInfo$MLparm <- .scale(beta=object$fixef,X=X.pv)[parm]
  } else processed$intervalInfo$MLparm <- object$fixef[parm]   
  if (object$spaMM.version < "4.1.42") {
    optimInfo <- attr(object,"optimInfo") ## may be NULL
    trTemplate <- .old_hack_trTemplate(optimInfo, object)
  } else {
    # attr(object,"optimInfo") still exists but is much ess suitable that the new one from 'lc'.
    optimInfo <- attr(lc,"optimInfo") ## may be NULL (from HLCor, HLfit, or length(initvec)=0 in other fitting fns)
    trTemplate <- optimInfo$`init.optim` ## may be NULL if optimInfo is NULL or if  optimInfo is not NULL but no par as outer optimized
  }  
  # For partial ranCoefs, it should contains a full vector 
  # when HLfit_body() -> .canonizeRanPars() is reached, there must be ranCoefs with the constraints and 
  # trRancoefs with a fully varying vector
  if ( ! is.null(trTemplate)) { # ... not from HLfit()...
    olc <- lc ## olc is working copy
    LUarglist <- optimInfo$LUarglist
    if (paste(lc[[1]])=="HLCor") attr(trTemplate,"moreargs") <- .get_moreargs(object)
    ## locoptim expects a fn with first arg ranefParsVec
    objfn <- function(ranefParsVec, anyHLCor_obj_args=NULL, HLcallfn.obj=NULL,
                      objfn.extras) { 
      ranefParsList <- relist(ranefParsVec,trTemplate)
      if (length(.unlist(trTemplate$trRanCoefs)) &&
          (length(objfn.extras[["user.lower"]]$ranCoefs) || length(objfn.extras[["user.upper"]]$ranCoefs))
      ) ranefParsList  <- .apply_transformed_box_constr(fix=ranefParsList, skeleton=trTemplate, 
                                            user.lower=objfn.extras[["user.lower"]], 
                                            user.upper=objfn.extras[["user.upper"]], transf=TRUE)
      olc$fixed <- structure(.modify_list(olc$fixed,ranefParsList)) ## replaces ! some elements and keeps the "type" !
      locfit <- eval(as.call(olc)) ## HLfit call with given ranefParsVec
      resu <- (posforminimiz)*locfit$fixef[parm]
      attr(resu,"info") <- locfit$APHLs$p_v 
      ## attribute lost by optim but otherwise useful for debugging 
      #print(olc$fixed)
      #print(resu)
      return(resu) ## return value to be optimized is a parameter value, not a likelihood
    }
    rC_transf <- .spaMM.data$options$rC_transf
    LUarglist$canon.init <- .canonizeRanPars(ranPars=trTemplate,
                                             corr_info=.get_from_ranef_info(object), 
                                             checkComplete=FALSE, rC_transf=rC_transf)
    LowUp <- do.call(.makeLowerUpper,LUarglist)
    optim_bound_over_nuisance_pars <- function(posforminimiz) { ## optimize the CI bound and returns the parameters that optimize this bound
      user_init_optim <- switch(fittingFunction,
                                "corrHLfit" = llc[["init.corrHLfit"]],
                                "fitme" = llc[["init"]], 
                                NULL)
      init <- unlist(trTemplate)
      if (paste(lc[[1]])=="HLCor") { HLcallfn_obj <- "HLCor.obj" } else HLcallfn_obj <- "HLfit.obj"
      .assignWrapper(anyObjfnCall.args$processed,
                     paste0("return_only <- \"confint_bound\""))
      optr <- .new_locoptim(init.optim=trTemplate,LowUp=LowUp, objfn.extras=LUarglist,
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
    processed$intervalInfo$init <- init_beta[parm]
    processed$intervalInfo$init_v_h <- object$v_h
    if (! is.null(trTemplate)) {
      anyObjfnCall.args <- as.list(lc[-1L]) ## includes processed, ranPars, controlS.dist, control.HLfit...
      anyObjfnCall.args$skeleton <- trTemplate
      # The objective function 'objfn' returns the confint bound given the corr pars. Thus locoptim maximizes the confint bound over the the corr pars
      olc <- lc ## that's olc that is used in the objective fn !
      posforminimiz <- 1 ## defined in the envir where objfn is defined... (bad style)
      bound <- optim_bound_over_nuisance_pars()
      # if (paste(lc[[1]])=="HLCor") {
      #   olc$fixed <- structure(.modify_list(olc$fixed,bound)) ## replaces ! some elements and keeps the "type" (lazyness)!
      # } else {
      #   olc$ranFix <- structure(.modify_list(olc$ranFix,bound)) ## replaces ! some elements and keeps the "type" !
      # }
      olc$fixed <- structure(.modify_list(olc$fixed,bound)) ## replaces ! some elements and keeps the "type" !
      ## recover fit for optimized params (must use call with intervalInfo and LevenbergM=FALSE)
      lowerfit <- eval(as.call(olc)) ## full HLfit objectobject
      attr(lowerfit,"optimInfo") <- optimInfo ## expected by summary.HLfit
      # lowerfit <- .update_ranef_info(lowerfit, moreargs=LUarglist$moreargs)
      ##
    } else lowerfit <- eval(as.call(lc))
    notconv <- logLik(lowerfit)<  processed$intervalInfo$targetlik-0.001 || 
      ! is.null(lowerfit$warnings$innerNotConv)
    if (notconv) {
      fac <- 2L*fac
    } else break
  }
  if (notconv) {
    cat("possible convergence problem for lower bound.\n")
  } else if (verbose) prevmsglength <- overcat(" ...lower bound converged                               \n",prevmsglength) 
  
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
    processed$intervalInfo$init <- init_beta[parm]
    if (! is.null(trTemplate)) {
      olc <- lc
      posforminimiz <- -1 ## maximization
      bound <- optim_bound_over_nuisance_pars()
      # if (paste(lc[[1]])=="HLCor") {
      #   olc$fixed <- structure(.modify_list(olc$fixed,bound)) ## replaces ! some elements and keeps the "type" !
      # } else {
      #   olc$ranFix <- structure(.modify_list(olc$ranFix,bound)) ## replaces ! some elements and keeps the "type" !
      # }
      olc$fixed <- structure(.modify_list(olc$fixed,bound)) ## replaces ! some elements and keeps the "type" !
      upperfit <- eval(as.call(olc))
      attr(upperfit,"optimInfo") <- optimInfo ## expected by summary.HLfit
      # upperfit <- .update_ranef_info(upperfit, moreargs=LUarglist$moreargs)
      ##
    } else upperfit <- eval(as.call(lc))
    notconv <- logLik(upperfit)<  processed$intervalInfo$targetlik-0.001 || 
      ! is.null(upperfit$warnings$innerNotConv)
    if (notconv) {
      fac <- 2L*fac
    } else break
  }
  if (notconv) {
    cat("possible convergence problem for lower bound.\n")
  } else if (verbose) prevmsglength <- overcat(" ...upper bound converged                               ",prevmsglength) 
  
  if (verbose && prevmsglength) cat("\n")
  options(warnori)
  # .options.processed(lc$processed, oldopt)
  interval <- c(lowerfit$fixef[parm],upperfit$fixef[parm])
  names(interval) <- paste(c("lower","upper"),parm)
  if (verbose) {
    if (any(object$models[["phi"]]=="phiHGLM") && processed$verbose["phifit"]) { 
      cat("\n")
      processed$fitenv$prevmsglength <- 0L
    } # newline after the phifit progress mess before printing CIs
    print(interval)
  }
  resu <- list(lowerfit=lowerfit,upperfit=upperfit,interval=interval)
  if (! is.null(resu$confint_best_fit <- processed$envir$confint_best)) {
    locmess <- paste("Element 'confint_best_fit' of the return object contain information about a possible better fit to the data",
                     "\nYou can for example refit the model using the parameter values shown in that element as initial values.")
    message(locmess)
  }
  return(resu)
}

# With preprocessed p4m:
# (1) Get 
# update(object, 
#   etaFix=list(beta=setNames(0,parm)), 
#   control=list(p4m="H", get_proc_call_4_p4m=TRUE))
# -> The get_proc_call_4_p4m flag to retrieve a (HL...body) processed call
# from which we use the $processed envir (in a .p4m_by_iters call).
# -> etaFix, to set up in this case a processed$X_off_Xb_fn function
# that will then be called in HLfit_body() -> possibly inner fns -> .get_off() for any value of etaFix.
# The objfn of a function outer optimizing a beta parm must then provide 
# each new value of parm through etaFix to the .p4m_by_iters call.
#
# X_off_Xb_fn() updates processed$off (and returns it) for new
# 'offsets' in the standard sense (offset() terms 
#   whose value is given by model.offset())
# and for new etaFix values (always taken into account through $off). 
# X_off_Xb_fn must be called in two contexts: 
# -> in HLfit_body() and inner fns -> .get_off(), to update the processed$off 
#   as fn of new etaFix and retrieve this off; (confint, numinfo, & 
#    outer-beta procedures by fitme/fitmv)
# -> .p4m_by_iters, each time the .dynoffset is updated. 

# Overall the idea is to optimize the fixed 'parm' coef by using 
# (1) an objective function that takes it and ranPars as arguments and *returns the parm value*, 
# (2) a constraint function that checks the logLik.
# In addition, an ad-hoc opt_env keep info as side effect of the objfn, to 
# (1) avoid fitting the same model both in the objfn and the constrfn, 
# (2) compensate for some possible deficiencies of nloptr use or of nloptr's constrained optim procedure.
.confint_LRT_single_par_p4m <- function(level, parm, object, verbose, ...) {
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
  llc$data <- object$data # llc$data may not even have .dynoffset 
  control <- llc$control
  control$p4m <- "H" # ranPars are fixed in all fits, so this should be OK.
  # hack for the progress bar (vs control$port_env which has distinct role signalling 2nd step):
  control[["meta_port_env"]] <- list2env(list(prevmsglength=0L,IT=0L), parent = emptyenv())
  control$"p4m_reactvt_warn" <- FALSE
  dataWdyn <- object$data # these ones include a .dynoffset (and the "best" one)
  
  if (use_proc_call <- identical(control[["use_proc_call"]],TRUE) ) { 
    if (...length() &&
        length(badmatches <- intersect(...names(), c("control","etaFix")))) 
      stop(paste(paste(badmatches,collapse=","),"are not allowed in this context."))
    HL_body_call <- update(object, 
                           etaFix=list(beta=setNames(0,parm)), # cf comments above
                           control=list(p4m="H", get_proc_call_4_p4m=TRUE),
                           ...)
    # .p4m_by_iters call rather than update fit:
    processed_p4m_call <- llc
    processed_p4m_call$processed <- HL_body_call$processed
    processed_p4m_call[[1L]] <- get(".p4m_by_iters", asNamespace("spaMM"), inherits=FALSE)  
    processed_p4m_call$control <- .reformat_p4m_controls(control, has_bar=TRUE)
    processed_p4m_call$submodels <- llc$submodels
  } else {
    update_args <- list(object=object, control=control, data=dataWdyn)
    if (length(setdiff(names(verbose),"confint"))) {
      fitverbose <- llc$verbose
      fitverbose <- .reformat_verbose(fitverbose)
      fitverbose <- .modify_list(fitverbose, verbose)
      update_args$verbose <- fitverbose
    } 
    if (...length()) update_args <- .modify_list(update_args, list(...))
  }
  targetlik <- logLik(object)-dlogL
  HL <- object$HL
  fixeflik <- switch(paste(HL[1L]),
                     "0"=if (object$models[["eta"]]=="etaGLM") "p_v" else "hlik", # if the user fitted a GLM by PQL/L there is no hlik 
                     "1"="p_v",
                     stop(paste("confint does not yet handle HLmethod",paste(HL,collapse=" "),
                                "(or ",c(llc$method,llc$HLmethod),").",sep=" ")))
  beta_cov <- .get_beta_cov_any_version(object)
  beta_se <- sqrt(diag(x=beta_cov))[parm]
  asympto_abs_Dparm <- znorm* beta_se
  posforminimiz <- NULL
  #
  X.pv <- object$X.pv
  X_is_scaled <- ( ! is.null(attr(X.pv,"scaled:scale")))
  
  # Maybe not general code:
  ranPars <- get_ranPars(object, lambda_names=TRUE) 
  fittedPars <- get_fittedPars(object)
  fixedPars <- .remove_from_cP(ranPars, u_names=names(unlist(fittedPars)))
  
  # attr(object,"optimInfo") still exists ?? less suitable that the new one from some internal processed call ??.
  optimInfo <- attr(object,"optimInfo") ## may be NULL (from HLCor, HLfit, or length(initvec)=0 in other fitting fns)
  trTemplate <- optimInfo$`optim.pars` # (solution), NOT optimInfo$init.optim 
  # => trTemplate may be NULL if optimInfo is NULL or if  optimInfo is not NULL but no par as outer optimized
  # The objfn calls only .p4m_by_iters, not .p4m_by_outer_optim. Hence all ranPars must be in trTemplate,
  # otherwise the CI wouldn't be based on the profile lik bc ranPars would not maximize logL.
  # For partial ranCoefs, it should contains a full vector 
  # when HLfit_body() -> .canonizeRanPars() is reached, there must be ranCoefs with the constraints and 
  # trRancoefs with a fully varying vector
  if ( ! is.null(trTemplate)) { # 
    LUarglist <- optimInfo$LUarglist
    ## minimize focal parameter subject to constraint on logL.
    ## both the objfn and the constrfn require a model fit.
    ## We create an environment shared by both fns. 
    ## The first fn called (objfn in most case, but not always) 
    ## provides the fit and stores it in the envir
    opt_env <- new.env(parent = emptyenv())
    opt_env$bestfit <- object
    
    .nuispars2fullfix_with_constr <- function(ranefParsVec) {
      ranefParsList <- relist(ranefParsVec,trTemplate)
      if (length(.unlist(trTemplate$trRanCoefs)) &&
          (length(LUarglist[["user.lower"]]$ranCoefs) || length(LUarglist[["user.upper"]]$ranCoefs))
      ) ranefParsList  <- .apply_transformed_box_constr(fix=ranefParsList, skeleton=trTemplate, 
                                                        user.lower=LUarglist[["user.lower"]], 
                                                        user.upper=LUarglist[["user.upper"]], transf=TRUE)
      
      fixed <- .modify_list(fixedPars, ranefParsList)
      fixed <- .canonizeRanPars(fixed, corr_info=object$ranef_info$sub_corr_info, 
                                checkComplete=FALSE, rC_transf=.spaMM.data$options$rC_transf)
      fixed
    }
    
    data <- object$data
    
    if (use_proc_call) {
      .eval_fit <- function(parvec, init.HLfit, init.dynoffset) {
        if ( ! is.null(init.dynoffset)) 
          processed_p4m_call$processed$data$".dynoffset" <- init.dynoffset
        if (length(init.HLfit)) processed_p4m_call$init.HLfit <- init.HLfit
        processed_p4m_call$etaFix <- list(beta=setNames(parvec[1], parm))
        processed_p4m_call$fixed <- .nuispars2fullfix_with_constr(parvec[-1])
        locfit <- eval(processed_p4m_call, parent.frame())
        if (targetlik <= logLik(locfit) && # 2nd crit is TRUE if parm is more extreme than previous 'best' parm 
            posforminimiz*(locfit$fixef[parm]-opt_env$bestfit$fixef[parm])<0) {
          # cat(cli::col_green("bestparvec updated"))
          opt_env$bestfit <- locfit
          opt_env$bestparvec <- parvec
          # cat(parvec[1])
        }
        locfit
      }
    } else {
      .eval_fit <- function(parvec, init.HLfit, init.dynoffset) {
        if ( ! is.null(init.dynoffset)) {
          # cat("\n+")
          data$".dynoffset" <- init.dynoffset
          update_args$data <- data
        }
        update_args$etaFix <- list(beta=setNames(parvec[1], parm))
        update_args$fixed <- .nuispars2fullfix_with_constr(parvec[-1])
        if (length(init.HLfit)) update_args$init.HLfit <- init.HLfit
        locfit <- do.call(update, update_args) 
        if (targetlik <= logLik(locfit) && # 2nd crit is TRUE if parm is more extreme than previous 'best' parm 
            posforminimiz*(locfit$fixef[parm]-opt_env$bestfit$fixef[parm])<0) {
          # cat(cli::col_green("bestparvec updated"))
          opt_env$bestfit <- locfit
          opt_env$bestparvec <- parvec
          # cat(parvec[1])
        }
        locfit
      }
    }
    
    thresh_old <- c(beta_se/10, 
                    rep(0.01, length(unlist(trTemplate,use.names = FALSE))))
    
    get_init.dynoffset <- function(parvec, oldparvec, verbose=FALSE) {
      if ( (! is.null(oldparvec)) &&
          all(abs(oldparvec-parvec)< thresh_old)) {
        opt_env$locfit$data$.dynoffset
      } else if (( ! is.null(bestparvec <- opt_env$bestparvec)) &&
                  all(abs(bestparvec-parvec)< thresh_old)) {
        opt_env$bestfit$data$.dynoffset
      } else NULL # presumably equivalent to object$data$.dynoffset
    }
    
    
    get_init.HLfit <- function(parvec, oldparvec, bestparvec=opt_env$bestparvec) { 
      # return(list())  # cf comment 'weird' below => this fn is kept for future attempts.
      init.HLfit <- list()
      if ((! is.null(oldparvec)) && 
          abs(oldparvec-parvec)[1] < thresh_old[1]) {
        beta <- na.omit(fixef(opt_env$locfit))
        beta <- beta[names(beta) !=parm]
        init.HLfit$fixef <- beta
      } else if ((! is.null(bestparvec)) &&
                 abs(bestparvec-parvec)[1] < thresh_old[1]) {          
          beta <- na.omit(fixef(opt_env$bestfit))
          beta <- beta[names(beta) !=parm]
          init.HLfit$fixef <- beta
      } # else no init beta

      if (( ! is.null(oldparvec)) && 
          all(abs(oldparvec-parvec)[-1] < thresh_old[-1])) {
        init.HLfit$v_h <- ranef(opt_env$locfit, type="bare.init")
      } else if ((! is.null(bestparvec)) &&
                 abs(bestparvec-parvec)[-1] < thresh_old[-1]) {          
        init.HLfit$v_h <- ranef(opt_env$bestfit, type="bare.init")
      } # else no init v_h
      init.HLfit
    }
    
    constrfn <- function(parvec) {
      if (identical(opt_env$parvec,parvec)) {
        locfit <- opt_env$locfit
        # print("locfit from env in constrfn")
      } else {
        locfit <- .eval_fit(parvec, 
                            init.HLfit = get_init.HLfit(parvec, oldparvec=opt_env$parvec), 
                            init.dynoffset = get_init.dynoffset(parvec, oldparvec=opt_env$parvec))
        opt_env$parvec <- parvec
        opt_env$locfit <- locfit
        # print("locfit refit in constrfn")
      }
      targetlik-logLik(locfit) # must be negative, ie logLik > targetlik
    }
    
    #  objfn -> .eval_fit optimizes logL over the other fixed-effect coefficients
    # It returns *the input focal parm* (with p_v as attribute).
    objfn <- function(parvec) {
      if (identical(opt_env$parvec,parvec)) {
        # print("locfit from env in objfn")
        locfit <- opt_env$locfit
      } else {
        locfit <- .eval_fit(parvec, 
                            init.HLfit = get_init.HLfit(parvec, oldparvec=opt_env$parvec), 
                            init.dynoffset = get_init.dynoffset(parvec, oldparvec=opt_env$parvec))
        opt_env$parvec <- parvec
        opt_env$locfit <- locfit
        # print("locfit from fit in objfn")
      }
      resu <-  (posforminimiz)*parvec[1] # (posforminimiz)*locfit$fixef[parm]
      attr(resu,"info") <- locfit$APHLs$p_v 
      return(resu) ## return value to be optimized is a parameter value, not a likelihood
    }
    
    # I need to convert back the trTemplate to get LowUp, but this LowUp
    # is itself in transformed space, as is fitting for nloptr input.
    rC_transf <- .spaMM.data$options$rC_transf
    LUarglist$canon.init <- .canonizeRanPars(ranPars=trTemplate,
                                             corr_info=.get_from_ranef_info(object), 
                                             checkComplete=FALSE, rC_transf=rC_transf)
    LowUp <- do.call(.makeLowerUpper,LUarglist) # _____F I X M E_____ high lambda values make everything difficult

    optim_bound_over_nuisance_pars <- function(posforminimiz) { ## optimize the CI bound and returns the parameters that optimize this bound
      init <- unlist(trTemplate)
      if (posforminimiz>0) {
        init <- c(fixef(object)[parm]-asympto_abs_Dparm,init)
        fac <- 2
        upper <- c(fixef(object)[parm],unlist(LowUp$upper))
        lower <- c(fixef(object)[parm]-asympto_abs_Dparm*fac,unlist(LowUp$lower)) 
      } else {
        init <- c(fixef(object)[parm]+asympto_abs_Dparm,init)
        lower <- c(fixef(object)[parm],unlist(LowUp$lower))
        fac <- 2
        upper <- c(fixef(object)[parm]+asympto_abs_Dparm*fac,unlist(LowUp$upper))
      }
      # See alternative experiments in ...doc_code/constrOptim_auglag.R
      optr <- .constrOptim(init=init,
                           lower=lower, 
                           upper=upper,
                           neg_ineq_constrfn=constrfn,
                           eq_constrfn=NULL,
                           objfn=objfn, # uses posforminimiz in its definition 
                           control=list(), 
                           verbose=FALSE) 
      return(optr) # with $objective= focal parameter bound and $solution = a parvec including it and a ranefParsVec
    }
    
    ## lower fit (objfn using posforminimiz <- 1); then upper fit (... -1)
    posforminimiz <- 1 
    abyss <- optim_bound_over_nuisance_pars(posforminimiz=posforminimiz) 
    lowerfit <- opt_env$bestfit
    
    opt_env <- new.env(parent = emptyenv())
    opt_env$bestfit <- object  # but no $bestparvec yet.
    posforminimiz <- -1 
    abyss <- optim_bound_over_nuisance_pars(posforminimiz=posforminimiz)
    upperfit <- opt_env$bestfit
    
    interval <- c(lowerfit$fixef[parm],upperfit$fixef[parm])
    
  } else { # no fitted ranPars
    ## uniroot: define range
    fac <- 1L 
    while(fac < 1e6) {
      init_beta <- object$fixef-asympto_abs_Dparm*fac
      # if (X_is_scaled) init_beta <- .scale(beta=init_beta,X=X.pv) ## using locally saved X.pv
      # : do not scale a value that will be passed as etaFix argument!
      etaFix <- .modify_list(llc$etaFix, 
                             list(beta=setNames(init_beta[parm], parm)))
      lowerfit <- update(object, etaFix=etaFix, data=dataWdyn, control=control) 
      if (logLik(lowerfit)<targetlik) break
      fac <- fac*2L
    }
    lo <- init_beta[parm]
    fac <- 1L 
    while(fac < 1e6) {
      init_beta <- object$fixef+asympto_abs_Dparm*fac
      # if (X_is_scaled) init_beta <- .scale(beta=init_beta,X=X.pv) ## using locally saved X.pv
      etaFix <- .modify_list(llc$etaFix, 
                             list(beta=setNames(init_beta[parm], parm)))
      upperfit <- update(object, etaFix=etaFix, data=dataWdyn, control=control) 
      if (logLik(upperfit)<targetlik) break
      fac <- fac*2L
    }
    hi <- init_beta[parm]
    
    obj_fn <- function(v) {
      etaFix <- .modify_list(llc$etaFix, 
                             list(beta=setNames(v, parm)))
      lowerfit <- update(object, etaFix=etaFix, data=dataWdyn, control=control) 
      logLik(lowerfit)-targetlik
    }
    
    interval <- c(
      uniroot(obj_fn,lower=lo,upper=object$fixef[parm])$root,
      uniroot(obj_fn,lower=object$fixef[parm],upper=hi)$root
    )
    etaFix <- .modify_list(llc$etaFix, 
                           list(beta=setNames(interval[1], parm)))
    lowerfit <- update(object, etaFix=etaFix, data=dataWdyn, control=control) 
    if (lowerfit$warnings$succInnerNotConv) { 
      lowerfit$warnings$succInnerNotConv <- "Non-convergence issue in this fit." 
    } else lowerfit$warnings$succInnerNotConv <- NULL
    etaFix <- .modify_list(llc$etaFix, 
                           list(beta=setNames(interval[2], parm)))
    upperfit <- update(object, etaFix=etaFix, data=dataWdyn, control=control) 
    if (upperfit$warnings$succInnerNotConv) { 
      upperfit$warnings$succInnerNotConv <- "Non-convergence issue in this fit." 
    } else upperfit$warnings$succInnerNotConv <- NULL
  }
  if (control[["meta_port_env"]]$IT>0L) cat("\n")
  names(interval) <- paste(c("lower","upper"),parm)
  if (verbose) print(interval)
  resu <- list(lowerfit=lowerfit,upperfit=upperfit,interval=interval)
  return(resu)
}

confint.HLfit <- function(object, parm, level=0.95, verbose=TRUE, 
                          boot_args=NULL, format="default", ...) {
   
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
    boot_res <- .confint_boot(boot_args, object, expr_t, t_fn, parm, boot.ci, 
                              level, verbose = verbose, ...) 
    if (format=="stats") {
      attr(boot_res,"table")
    } else invisible(boot_res)
  } else {
    if ( .REMLmess(object,return_message=FALSE)) {
      warning("REML fits are not quite suitable for computing intervals for fixed effects.")
    }
    resu <- .confint_LRT(level, parm, object, verbose, ...)
    if (format=="stats") {
      attr(resu,"table")
    } else invisible(resu)
  }
}

# for a respint() concept, see LawlessF05
# requires the user to provide a function that simulates the new obs (!= data) to be predicted
# and a function that evaluates the pdf of such new obs for given parameters (argh)

