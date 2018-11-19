## next steps are COMPoisson optim, detection des modeles lambda, gestion des lambda multiples
#  et meilleurs default phi, lambda 
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
                names(formals(.spaMM.data$options$mat_sqrt_fn)),names(formals(make_scaled_dist))))  ## cf parallel code in HLCor.obj
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
  sparse_precision <- proc1$sparsePrecisionBOOL
  #
  user_init_optim <- init ## user_init_optim only for a check in new_locoptim, true initial value init.optim is modified below
  optim_blob <- .calc_optim_args(proc1=proc1, processed=processed,
                                 init=init, fixed=fixed, lower=lower, upper=upper, 
                                 verbose=verbose, optim.scale=optim.scale, For="fitme") 
  # modify HLCor.args and <>bounds;   ## distMatrix or uniqueGeo potentially added to HLCor.args:
  # init <- optim_blob$inits$`init` ## list; keeps all init values, all in untransformed scale
  init.optim <- optim_blob$inits$`init.optim` ## list; subset of all estimands, as name implies, and in transformed scale
  init.HLfit <- optim_blob$inits$`init.HLfit` ## list; subset as name implies 
  fixed <- optim_blob$fixed
  corr_types <- optim_blob$corr_types
  moreargs <- optim_blob$moreargs
  LUarglist <- optim_blob$LUarglist
  LowUp <- optim_blob$LowUp
  lower <- LowUp$lower ## list ! which elements may have length >1 !
  upper <- LowUp$upper ## list !
  #

  #
  if ( ! is.null(residproc1 <- proc1$residProcessed) && identical(spaMM.getOption("outer_optim_resid"),TRUE)) {
    ## problem is that this is useful if we can avoid computation of the leverages of the resid Model
    ## (we already need the leverages of the 'mean' model to define the resid model response)
    resid_optim_blob <- .calc_optim_args(proc1=residproc1, processed=proc1,
                                         init=proc1$residModel$init, fixed=proc1$residModel$fixed, ## all user input must be in proc1$residModel
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
    moreargs <- c(moreargs,list(resid=resid_optim_blob$moreargs))
  }
  
  
  processedHL1 <- proc1$HL[1] ## there's also HLmethod in processed<[[]]>$callargs
  needHLCor_specific_args <- (length(unlist(lower$corrPars)) || 
                                length(intersect(corr_types,c("Matern","Cauchy","adjacency","AR1","corrMatrix"))))
  if (needHLCor_specific_args) {
    HLcallfn.obj <- "HLCor.obj" 
    HLcallfn <- "HLCor"
    control.dist <- vector("list",length(moreargs))
    for (nam in names(moreargs)) control.dist[[nam]] <- moreargs[[nam]]$control.dist
    HLCor.args$control.dist <- control.dist # lapply(moreargs,getElement,name="control.dist")
    HLCor.args$ranPars <- fixed ## to be modified by objective function
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
    anyHLCor_obj_args$skeleton <- structure(init.optim, moreargs=moreargs,
                                            type=relist(rep("fix",length(initvec)),init.optim),
                                            moreargs=moreargs)
  } else {
    anyHLCor_obj_args$skeleton <- structure(init.optim,
                                            type=relist(rep("fix",length(initvec)),init.optim))
  }
  .assignWrapper(anyHLCor_obj_args$processed,
                   paste0("return_only <- \"",proc1$objective,"APHLs\""))
  if (length(initvec)) {
    
    if (identical(verbose["getCall"][[1L]],TRUE)) { ## toget an optim call with its initial value. Then HLcallfn is called and its call returned.
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
                        #preprocess.formal.args=preprocess.formal.args, 
                        control.corrHLfit=control,
                        verbose=verbose[["iterateSEM"]],
                        nb_cores=nb_cores)
        optr <- .probitgemWrap("iterateSEMSmooth",arglist=loclist, pack="probitgem") # do.call("iterateSEMSmooth",loclist) 
        optPars <- relist(optr$par,init.optim)
        if (!is.null(optPars)) attr(optPars,"method") <-"optimthroughSmooth"
        refit_info <- control[["refit"]] ## may be a NULL/ NA/ boolean or a list of booleans 
        if ( is.null(refit_info) || is.na(refit_info)) refit_info <- FALSE ## alternatives are TRUE or an explicit list or NULL
      } else {
        optPars <- .new_locoptim(init.optim, ## F I X M E try to use gradient? But neither minnqa nor _LN_BOBYQA use gradients. optim() can
                                 LowUp, 
                                 control, objfn_locoptim=.objfn_locoptim, 
                                 HLcallfn.obj=HLcallfn.obj, anyHLCor_obj_args=anyHLCor_obj_args, 
                                 user_init_optim=user_init_optim,
                                 grad_locoptim=NULL) # .grad_locoptim)
        refit_info <- attr(optPars,"refit_info") ## 'derives' from control[["refit"]] with modification ## may be NULL but not NA
      }
      optim_time <- .timerraw(time2)
      ranPars_in_refit <- structure(.modify_list(fixed,optPars),
                                   type=.modify_list(relist(rep("fix",length(unlist(fixed))),fixed), #attr(fixed,"type"),
                                                     relist(rep("outer",length(unlist(optPars))),optPars)))
      augZXy_phi_est <- proc1$augZXy_env$phi_est ## may be NULL: if phi was not estimated by augZXy
      
      any_nearly_singular_covmat <- FALSE
      if (! is.null(optPars$trRanCoefs)) {
        ranCoefs <- optPars$trRanCoefs # copy names...
        for (char_rd in names(optPars$trRanCoefs)) {
          ranCoefs[[char_rd]] <- .ranCoefsInv(optPars$trRanCoefs[[char_rd]])
          covmat <- .calc_cov_from_ranCoef(ranCoefs[[char_rd]])
          if (kappa(covmat)>1e14 || min(eigen(covmat,only.values = TRUE)$values)<1e-6) any_nearly_singular_covmat <- TRUE
        }
      }
      init_refit <- list()
      if ( is.null(refit_info) ) { ## not the case with SEM
        # if (any_nearly_singular_covmat) message("Nearly-singular cov matrix => refitting random-coefficient parameters by default.")
        if (any_nearly_singular_covmat) message("Nearly-singular cov matrix => Consider refitting with refit_info=TRUE.")
        refit_info <- list(phi=FALSE,lambda=(! is.null(optPars$trRanCoefs)),ranCoefs=(FALSE && any_nearly_singular_covmat))
      } else if ( is.list(refit_info) ) { ## never the default
        if (any_nearly_singular_covmat) {
          if (is.null(refit_info$ranCoefs)) {
            # message("Nearly-singular cov matrix => refitting random-coefficient parameters by default.")
          } else if ( ! refit_info$ranCoefs ) message("Nearly-singular cov matrix => Consider refitting with refit_info$ranCoefs=FALSE.")
        } 
        refit_info <- .modify_list( list(phi=FALSE,lambda=(! is.null(optPars$trRanCoefs)),
                                         ranCoefs=(FALSE && any_nearly_singular_covmat)), refit_info) 
        ## lambda=TRUE becomes the dafault (test random-slope 'ares' shows the effect)
      } else {
        if (any_nearly_singular_covmat) {
          if ( ! refit_info) message("Nearly-singular cov matrix => Consider refitting with refit_info=TRUE.")
        } 
        refit_info <- list(phi=refit_info,lambda=refit_info,ranCoefs=refit_info) # default with SEM is all FALSE
      }
      # refit:
      if (identical(refit_info$phi,TRUE)) {
        if ( ! is.null(optPars$trPhi)) {
          init_refit$phi <- .dispInv(optPars$trPhi)
          ranPars_in_refit$trPhi <- NULL
          attr(ranPars_in_refit,"type")$trPhi <- NULL 
        } else if (proc1$augZXy_cond ) {
          init_refit$phi <- augZXy_phi_est # might still be NULL...
        }
      }
      # refit, or rescale by augZXy_phi_est even if no refit:
      if ( ! is.null(augZXy_phi_est) || identical(refit_info$lambda,TRUE)) {
        if (! is.null(optPars$trLambda)) { # does NOT include the lambdas of the ranCoefs
          lambdapos <- as.integer(names(optPars$trLambda))
          lambda <- rep(NA,max(lambdapos))
          names(lambda) <- paste(seq_len(length(lambda)))
          lambda[lambdapos] <- .dispInv(optPars$trLambda)
          if ( ! is.null(augZXy_phi_est)) { lambda <- lambda * augZXy_phi_est }
          if (identical(refit_info$lambda,TRUE)) {
            init_refit$lambda <- lambda
            ranPars_in_refit$trLambda <- NULL
            attr(ranPars_in_refit,"type")$trLambda <- NULL
          } else ranPars_in_refit$trLambda <- .dispFn(lambda) ## rescale, but no refit
        }
      }
      ## refit, or rescale by augZXy_phi_est even if no refit:
      if ( ! is.null(augZXy_phi_est)  || identical(refit_info$ranCoefs,TRUE)) {
        if (! is.null(optPars$trRanCoefs)) {
          for (char_rd in names(ranCoefs)) {
            rancoef <- ranCoefs[[char_rd]]
            if ( ! is.null(augZXy_phi_est)) {
              Xi_ncol <- attr(rancoef,"Xi_ncol")
              lampos <- rev(length(rancoef) -cumsum(seq(Xi_ncol))+1L)  ## NOT cumsum(seq(Xi_cols))
              rancoef[lampos] <- rancoef[lampos] *augZXy_phi_est
              rancoef <- as.vector(rancoef) ## as.vector drops attributes
            } 
            if (identical(refit_info$ranCoefs,TRUE)) { 
              init_refit$ranCoefs[[char_rd]] <- rancoef
              ranPars_in_refit$trRanCoefs[char_rd] <- NULL
              attr(ranPars_in_refit,"type")$trRanCoefs[char_rd] <- NULL
            } else ranPars_in_refit$trRanCoefs[[char_rd]] <- .ranCoefsFn(rancoef)  
          }
        }      
      }
      if (length(init_refit)) HLCor.args$init.HLfit <- .modify_list(HLCor.args$init.HLfit, init_refit)  
    } ## end if ...getCall... else
    #
    # refit_info is list if so provided by user, else typically boolean. An input NA should have been converted to something else (not documented).
    if (needHLCor_specific_args) {
      attr(ranPars_in_refit,"moreargs") <- moreargs 
      HLCor.args$ranPars <- ranPars_in_refit 
    } else HLCor.args$ranFix <- ranPars_in_refit 
  } else if (len_ranPars <- length(unlist(HLCor.args$ranPars))){ ## Set attribute
    HLCor.args$ranPars <- structure(HLCor.args$ranPars,
                                  type = relist(rep("fix", len_ranPars), HLCor.args$ranPars))
  }
  #
  # if processed is an envir, the following is not local to anyHLCor_obj_args$processed but change processed globally
  .assignWrapper(HLCor.args$processed,"return_only <- NULL") 
  .assignWrapper(HLCor.args$processed,"verbose['warn'] <- TRUE") ## important!
  #browser()
  hlcor <- do.call(HLcallfn,HLCor.args) ## recomputation post optimization, or only computation if length(initvec)=0
  if (is.call(hlcor)) {
    ## then do.call(HLcallfn,HLCor.args) has retuned the call, not the fit. 
    ## see def of get_HLCorcall() for further explanation
    return(hlcor) ## HLCorcall
  }
  # hlcor may have received attr(.,"info.uniqueGeo") from HLCor_body.
  if (length(initvec)) {
    attr(hlcor,"optimInfo") <- list(LUarglist=LUarglist, optim.pars=optPars, 
                                    objective=proc1$objective,
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
  }
  ## substantial effect on object size! :
  lsv <- c("lsv",ls())
  rm(list=setdiff(lsv,"hlcor")) 
  ##
  return(hlcor)
}