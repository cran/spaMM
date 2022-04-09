fitme_body <- function(processed,
                       init=list(),
                       #                       init.HLfit=list(),
                       fixed=list(), ## NULL not valid (should be handled in preprocess?)
                       lower=list(),upper=list(),
                       # control.dist=list(), ## info in processed
                       control=list(), ## optimizer, <optimizer controls>, precision
                       nb_cores=NULL,
                       ... ## cf dotnames processing below 
) {
  dotlist <- list(...) ## forces evaluations, which makes programming easier...
  if (is.list(processed)) {
    proc1 <- processed[[1L]]
  } else proc1 <- processed
  verbose <-  proc1$verbose
  HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                names(formals(mat_sqrt)),names(formals(make_scaled_dist))))  ## cf parallel code in HLCor.obj
  ## fill HLCor.args
  good_dotnames <- intersect(names(dotlist),HLnames) ## those specifically for the called fns as def'd by HLnames
  if (length(good_dotnames)) {
    HLCor.args <- dotlist[good_dotnames]
  } else HLCor.args <- list() 
  ## replace some HLCor.args members  
  if (  is.list(processed) ) {
    pnames <- names(processed[[1]])
  } else pnames <- names(processed)
  for (st in pnames) HLCor.args[st] <- NULL 
  # 'processed' may be modified below, then will be copied in HLCor.args (and then removed from this envir for safety)
  optim.scale <- control[["optim.scale"]] 
  if (is.null(optim.scale)) optim.scale="transformed" ## currently no public alternative
  sparse_precision <- proc1$is_spprec
  #
  user_init_optim <- init #
  optim_blob <- .calc_optim_args(proc_it=proc1, processed=processed,
                                 user_init_optim=user_init_optim, fixed=fixed, lower=lower, upper=upper, 
                                 verbose=verbose, optim.scale=optim.scale, For="fitme") 
  # modify HLCor.args and <>bounds;   ## distMatrix or uniqueGeo potentially added to HLCor.args:
  # init <- optim_blob$inits$`init` ## list; keeps all init values, all in untransformed scale
  init.optim <- optim_blob$inits$`init.optim` ## list; subset of all estimands, as name implies, and in transformed scale
  init.HLfit <- optim_blob$inits$`init.HLfit` ## list; subset as name implies 
  fixed <- optim_blob$fixed
  corr_types <- optim_blob$corr_types
  LUarglist <- optim_blob$LUarglist
  moreargs <- LUarglist$moreargs
  LowUp <- optim_blob$LowUp
  lower <- LowUp$lower ## list ! which elements may have length >1 !
  upper <- LowUp$upper ## list !
  #

  #
  if ( ! is.null(residproc1 <- proc1$residProcessed) && identical(spaMM.getOption("outer_optim_resid"),TRUE)) {
    ## Problem is that outer optim at the mean model level is useful if we can avoid computation of the leverages 
    ## But here anyway we need the leverages of the 'mean' model to define the resid model response.
    resid_optim_blob <- .calc_optim_args(proc_it=residproc1, processed=proc1,
                                         user_init_optim=proc1$residModel$init, fixed=proc1$residModel$fixed, ## all user input must be in proc1$residModel
                                         lower=proc1$residModel$lower, upper=proc1$residModel$upper, ## all user input must be in proc1$residModel
                                         verbose=c(SEM=FALSE), optim.scale=optim.scale, For="fitme") 
    resid_init.optim <- resid_optim_blob$inits$`init.optim` ## list; subset of all estimands, as name implies, and in transformed scale
    proc1$residModel$`init.HLfit` <- resid_optim_blob$inits$`init.HLfit` ## list; subset as name implies 
    proc1$residModel$fixed <- resid_optim_blob$fixed
    ##### residproc1$corr_types <- resid_optim_blob$corr_types
    # resid_user.lower <- resid_optim_blob$user.lower
    # resid_user.upper <- resid_optim_blob$user.upper
    init.optim <- c(init.optim,list(resid=resid_init.optim))
    lower <- c(lower,list(resid=resid_optim_blob$LowUp$lower))
    upper <- c(upper,list(resid=resid_optim_blob$LowUp$upper))
    LowUp <- list(lower=lower,upper=upper)
    moreargs <- c(moreargs,list(resid=resid_optim_blob$LUarglist$moreargs))
  }
  
  
  processedHL1 <- proc1$HL[1] 
  needHLCor_specific_args <- (length(unlist(lower$corrPars)) || 
                                length(intersect(corr_types,c("Matern","Cauchy","adjacency","AR1","corrMatrix", "IMRF","corrFamily"))))
  if (needHLCor_specific_args) {
    HLcallfn.obj <- "HLCor.obj" 
    HLcallfn <- "HLCor"
    control.dist <- vector("list",length(moreargs))
    for (nam in names(moreargs)) control.dist[[nam]] <- moreargs[[nam]]$control.dist 
    HLCor.args[["control.dist"]] <- control.dist ## always reconstructed locally, not in the fitme_body call
    HLCor.args$ranPars <- fixed ## to be modified by objective function
  # } else if (processed$augZXy_cond) { # an interesting experiment, but not consistently used as replacement for HLfit.obj throughout spaMM)
  #   HLcallfn.obj <- ".augZXy_obj"
  #   HLcallfn <- "HLfit"
  #   HLCor.args$ranFix <- fixed  
  } else {
    HLcallfn.obj <- "HLfit.obj"
    HLcallfn <- "HLfit"
    HLCor.args$ranFix <- fixed  
  }
  HLCor.args$init.HLfit <- init.HLfit
  HLCor.args$processed <- processed ## for the <...>.obj and <...>_body functions  
  ## 
  anyHLCor_obj_args <- HLCor.args
  ## HLCor.obj uses a vector + skeleton
  initvec <- unlist(init.optim)
  if (needHLCor_specific_args) {
    anyHLCor_obj_args$skeleton <- structure(init.optim, moreargs=moreargs, ## moreargs is a list over ranefs 
                                            type=relist(rep("fix",length(initvec)),init.optim) )
  } else {
    anyHLCor_obj_args$skeleton <- structure(init.optim,
                                            type=relist(rep("fix",length(initvec)),init.optim))
  }
  .assignWrapper(anyHLCor_obj_args$processed,
                   paste0("return_only <- \"",proc1$objective,"APHLs\""))
  if (length(initvec)) {
    augZXy_phi_est <- NULL
    if (identical(verbose["getCall"][[1L]],TRUE)) { ## to get an optim call with its initial value. Then HLcallfn is called and its call returned.
      ## confint -> get_HLCorcall needs an HLCor call with the following ranFix
      ranPars_in_refit <- structure(.modify_list(fixed,init.optim),
                                    # I label this "F I X" as a TAG for this modif type attribute:
                                    type=.modify_list(relist(rep("fix",length(unlist(fixed))),fixed), #attr(fixed,"type"), 
                                                     relist(rep("outer",length(initvec)),init.optim)) )
    } else {
      use_SEM <- (!is.null(processedHL1) && processedHL1=="SEM")
      time2 <- Sys.time()
      if (use_SEM) {
        optimMethod <- "iterateSEMSmooth"
        if (is.null(proc1$SEMargs$control_pmvnorm$maxpts)) {
          .assignWrapper(processed,"SEMargs$control_pmvnorm$maxpts <- quote(250L*nobs)") 
        } ## else default visible in SEMbetalambda
        ## its names should match the colnames of the data in Krigobj = the  parameters of the likelihood surface. Current code maybe not general.
        loclist <- list(anyHLCor_obj_args=anyHLCor_obj_args,  ## contains $processed
                        LowUp=LowUp,init.corrHLfit=init.optim, ## F I X M E usage of user_init_optim probably not definitive
                        control.corrHLfit=control,
                        verbose=verbose[["iterateSEM"]],
                        nb_cores=nb_cores)
        optr <- .probitgemWrap("iterateSEMSmooth",arglist=loclist, pack="probitgem") # do.call("iterateSEMSmooth",loclist) 
        optPars <- relist(optr$par,init.optim)
        if (!is.null(optPars)) attr(optPars,"method") <-"optimthroughSmooth"
        refit_info <- control[["refit"]] ## may be a NULL/ NA/ boolean or a list of booleans 
        if ( is.null(refit_info) || is.na(refit_info)) refit_info <- FALSE ## alternatives are TRUE or an explicit list or NULL
      } else {
        optPars <- .new_locoptim(init.optim, ## try to use gradient? But neither minqa nor _LN_BOBYQA use gradients. optim() can
                                 LowUp, 
                                 control, objfn_locoptim=.objfn_locoptim, 
                                 HLcallfn.obj=HLcallfn.obj, anyHLCor_obj_args=anyHLCor_obj_args, 
                                 user_init_optim=user_init_optim,
                                 grad_locoptim=NULL, verbose=verbose[["TRACE"]])
        refit_info <- attr(optPars,"refit_info") ## 'derives' from control[["refit"]] with modification ## may be NULL but not NA
      }
      optim_time <- .timerraw(time2)
      augZXy_phi_est <- proc1$augZXy_env$phi_est ## may be NULL: if phi was not estimated by augZXy
      refit_args <- .get_refit_args(fixed, optPars, processed, moreargs, proc1, refit_info, HLCor.args, augZXy_phi_est)
      HLCor.args <- refit_args$HLCor.args
      ranPars_in_refit <- refit_args$ranPars_in_refit
      if ( ! is.null(trBeta <- ranPars_in_refit$trBeta)) { # outer beta
        processed$off <- environment(processed$X_off_fn)$ori_off
        processed$X.pv <- environment(processed$X_off_fn)$X_off
        processed$AUGI0_ZX <- .init_AUGI0_ZX(processed$X.pv, processed$AUGI0_ZX$vec_normIMRF, processed$ZAlist, nrand=length(processed$ZAlist), n_u_h=nrow(processed$AUGI0_ZX$ZeroBlock), 
                                             sparse_precision=processed$is_spprec, 
                                             as_mat=.eval_as_mat_arg(processed))
        HLCor.args$init.HLfit$fixef <- .betaInv(trBeta)
        ranPars_in_refit$trBeta <- NULL
        processed$port_env$port_fit_values$fixef <- NULL
      }
    } ## end if ...getCall... else
    #
    # refit_info is list if so provided by user, else typically boolean. An input NA should have been converted to something else (not documented).
    if (needHLCor_specific_args) {
      attr(ranPars_in_refit,"moreargs") <- moreargs 
      HLCor.args$ranPars <- ranPars_in_refit 
    } else HLCor.args$ranFix <- ranPars_in_refit 
  } else if (len_ranPars <- length(unlist(HLCor.args$ranPars))) { ## Set attribute
    HLCor.args$ranPars <- structure(HLCor.args$ranPars,
                                  type = relist(rep("fix", len_ranPars), HLCor.args$ranPars),
                                  moreargs=moreargs) ## moreargs needed if user handles fixed(<transformed params>) ('hyper' tests)
  }
  #
  # not local to anyHLCor_obj_args$processed: change processed globally
  .assignWrapper(HLCor.args$processed,"return_only <- NULL") 
  .assignWrapper(HLCor.args$processed,"verbose['warn'] <- TRUE") ## important!
  # Run in all cases to produce the full object (rather than only the optimization result):
  hlcor <- do.call(HLcallfn,HLCor.args) ## recomputation post optimization, or only computation if length(initvec)=0
  if (is.call(hlcor)) {
    ## then do.call(HLcallfn,HLCor.args) has retuned the call, not the fit. 
    ## see def of get_HLCorcall() for further explanation
    return(hlcor) ## HLCorcall
  }
  if (length(initvec)) {
    attr(hlcor,"optimInfo") <- list(LUarglist=LUarglist, optim.pars=optPars, 
                                    objective=proc1$objective,
                                    augZXy_phi_est=augZXy_phi_est, ## gives info to interpret optim.pars in confint.HLfit()
                                    optim_time=optim_time) ## processed was erased for safety
    locoptr <- attr(optPars,"optr")
    if (attr(optPars,"method")=="nloptr") {
      if (locoptr$status<0L) hlcor$warnings$optimMessage <- paste0("nloptr() message: ",
                                                                   locoptr$message," (status=",locoptr$status,")")
    } else if ( attr(optPars,"method")=="optim" ) {
      if (locoptr$convergence) hlcor$warnings$optimMessage <- paste0("optim() message: ",locoptr$message,
                                                                     " (convergence=",locoptr$convergence,")")
    } else if ( attr(optPars,"method")== "optimthroughSmooth") {
      # provide logL estimate from the smoothing, to be used rather than the hlcor logL :
      logLapp <- optr$value
      attr(logLapp,"method") <- "  logL (smoothed)" 
      hlcor$APHLs$logLapp <- logLapp
    }
    hlcor$warnings$suspectRho <- .check_suspect_rho(corr_types, ranPars_in_refit, LowUp)
    if ( ! is.null(PQLdivinfo <- processed$envir$PQLdivinfo)) {
      hlcor$divinfo <- PQLdivinfo
      hlcor$warnings$divinfo <- "Numerical issue detected; see div_info(<fit object>) for more information."
      warning(hlcor$warnings$divinfo)
    }
  }
  ## substantial effect on object size! :
  lsv <- c("lsv",ls())
  rm(list=setdiff(lsv,"hlcor")) 
  ##
  return(hlcor)
}