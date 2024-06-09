corrHLfit_body <- function(processed, ## possibly a list of environments
                           init.corrHLfit=list(),
                      fixed=list(),
                      lower=list(),upper=list(),
                      #control.dist=list(), ## info in processed
                      control.corrHLfit=list(), ## optim.scale, optimizer, <optimizer controls>
                      nb_cores=NULL,
                      ... ## cf dotnames processing below
) { 
  dotlist <- list(...) ## forces evaluations, which makes programming easier... contains the user's init.HLfit
  if (is.list(processed)) {
    proc1 <- processed[[1L]]
  } else proc1 <- processed
  verbose <-  proc1$verbose
  HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                names(formals(mat_sqrt)),names(formals(make_scaled_dist))))  ## cf parallel code in HLCor.obj
  ## fill HLCor.args
  good_dotnames <- intersect(names(dotlist),HLnames) ## those specifically for the called fns as def'd by HLnames # potentially contains "init.HLfit"
  if (length(good_dotnames)) {
    HLCor.args <- dotlist[good_dotnames] # may contain the user's init.HLfit
  } else HLCor.args <- list() 
  ## replace some HLCor.args members  
  if ( is.list(processed)) { ## list of environments
    pnames <- names(processed[[1]])
  } else pnames <- names(processed) # there's an "init_HLfit" but no "init.HLfit"
  for (st in pnames) HLCor.args[st] <- NULL 
  # 'processed' may be modified below, then will be copied in HLCor.args (and then removed from this envir for safety)
  optim.scale <- control.corrHLfit$optim.scale 
  if (is.null(optim.scale)) optim.scale="transformed" ## currently no public alternative
  sparse_precision <- proc1$is_spprec
  #
  user_init_optim <- init.corrHLfit 
  
  optim_blob <- .calc_optim_args(proc_it=proc1, processed=processed,
                                 user_init_optim=user_init_optim, fixed=fixed, lower=lower, upper=upper, 
                                 verbose=verbose, optim.scale=optim.scale, For="corrHLfit") 
  init.corrHLfit <- NaN ## make clear it's not to be used
  # #  code correct mais opaque pour l'utilisateur: init.optim va être ignore en cas de 1D optimization... (OK in fitme)
  # ## maximization VIA HLCor.obj
  # ## by maxim over (corrpars,phi,lambda, (beta)_[corrpars,phi,lambda])
  # ##     if trPhi,trLambda are in the init.optims
  # ## or by maxim over (corrpars,(beta,phi,lambda)_corrpars)
  # ##     otherwise.
  # ################ construct intervals for this maximization
  init.optim <- optim_blob$inits$`init.optim` ## list; subset of all estimands, as name implies, and in transformed scale
  init.HLfit <- optim_blob$inits$`init.HLfit` ## list; subset as name implies 
  fixed <- optim_blob$fixed
  corr_types <- optim_blob$corr_types
  LUarglist <- optim_blob$LUarglist
  LowUp <- optim_blob$LowUp
  lower <- LowUp$lower ## list ! which elements may have length >1 !
  upper <- LowUp$upper ## list !
  moreargs <- LUarglist$moreargs
  
  ## there is a block for residual disp models here in fitme_body(), where moreargs lower upper LowUp are modified
  
  control.dist <- vector("list",length(moreargs))
  for (nam in names(moreargs)) control.dist[[nam]] <- moreargs[[nam]]$control.dist
  HLCor.args$control.dist <- control.dist # lapply(moreargs,getElement,name="control.dist")
  processedHL1 <- proc1$HL[1] 
  HLCor.args$fixed <- fixed
  HLCor.args$processed <- processed
  ## 
  anyHLCor_obj_args <- HLCor.args
  ## HLCor.obj uses a vector + skeleton
  initvec <- unlist(init.optim)
  anyHLCor_obj_args$skeleton <- structure(init.optim, 
                                          moreargs=moreargs, ## moreargs is a list over ranefs 
                                          type=relist(rep("fix",length(initvec)),init.optim) )
  .assignWrapper(anyHLCor_obj_args$processed,
                   paste0("return_only <- \"",proc1$objective,"APHLs\""))
  use_SEM <- (!is.null(processedHL1) && processedHL1=="SEM"  && length(lower))
  if (use_SEM) {
    if (is.null(proc1$SEMargs$control_pmvnorm$maxpts)) {
      .assignWrapper(processed,"SEMargs$control_pmvnorm$maxpts <- quote(250L*nobs)") 
    } ## else default visible in SEMbetalambda
    ## its names should match the colnames of the data in Krigobj = the  parameters of the likelihood surface. Current code maybe not general.
    iterateSEMSmooth <- get("iterateSEMSmooth",envir = asNamespace("probitgem"), inherits=FALSE)
    optr <- iterateSEMSmooth(anyHLCor_obj_args=anyHLCor_obj_args,  ## contains $processed
                             LowUp=LowUp,init.corrHLfit=user_init_optim, ## FIXME usage of user_init_optim probably not definitive
                             control.corrHLfit=control.corrHLfit,
                             verbose=verbose[["iterateSEM"]],
                             nb_cores=nb_cores) 
    # optr <- .probitgemWrap("iterateSEMSmooth",arglist=loclist, pack="probitgem") # do.call("iterateSEMSmooth",loclist) 
    optPars <- relist(optr$par,init.optim)
    if (!is.null(optPars)) attr(optPars,"method") <-"optimthroughSmooth"
  } else { ## this is also called if length(lower)=0 by  (SEM or not) and optPars is then null
    if (.safe_true(verbose["getCall"][[1L]])) {
      optPars <- init.optim
    } else {
      #
      ## maximization VIA HLCor.obj
      ## by maxim over (corrpars,phi,lambda, (beta)_[corrpars,phi,lambda])
      ##     if trPhi,trLambda are in the init.optims
      ## or by maxim over (corrpars,(beta,phi,lambda)_corrpars)
      ##     otherwise.
      #
      optPars <- .new_locoptim(init.optim, LowUp=LowUp,anyHLCor_obj_args=anyHLCor_obj_args,
                               objfn_locoptim=.objfn_locoptim,
                               control=control.corrHLfit, 
                               user_init_optim=user_init_optim, verbose=verbose[["TRACE"]])
    }
  }
  if (!is.null(optPars)) {
    ranPars_in_refit <- structure(.modify_list(HLCor.args$fixed,optPars),
                                    # I write "F I X" as a TAG for this modif type attribute:
                                    type=.modify_list(.relist_rep("fix",HLCor.args$fixed), #attr(HLCor.args$fixed,"type"),
                                                      .relist_rep("outer",optPars)),
                                    moreargs=moreargs)
  } else {
    ranPars_in_refit <- structure(HLCor.args$fixed,
                                  type = .relist_rep("fix", HLCor.args$fixed),
                                  moreargs=moreargs)
  }
  ranPars_in_refit <- .expand_hyper(ranPars_in_refit, processed$hyper_info, moreargs=moreargs)
  HLCor.args$fixed <- ranPars_in_refit
  # not local to anyHLCor_obj_args$processed: change processed globally
  .assignWrapper(HLCor.args$processed,"return_only <- NULL") 
  .assignWrapper(HLCor.args$processed,"verbose['warn'] <- TRUE") ## important!
  hlcor <- do.call("HLCor",HLCor.args) ## recomputation post optimization (or only computation, if length(lower)=0)
  #
  if (is.call(hlcor)) { 
    if (length(initvec)) { # always TRUE ?
      attr(hlcor,"optimInfo") <- list(LUarglist=LUarglist, init.optim=init.optim,
                                      objective=proc1$objective,
                                      rC_transf=.spaMM.data$options$rC_transf)
    }
    return(hlcor) ## HLCorcall
  } else {
    if ( proc1$fitenv$prevmsglength) { # there was output for a phi-resid.model. The fit object may then be printed...
      cat("\n")
      proc1$fitenv$prevmsglength <- 0L
    }
  }
  #
  if (length(initvec)) { # always TRUE ?
    attr(hlcor,"optimInfo") <- list(LUarglist=LUarglist, optim.pars=optPars,
                                    objective=proc1$objective,
                                    rC_transf=.spaMM.data$options$rC_transf
                                    #augZXy_phi_est=NULL, # implicit NULL bc augZXy in prevented for outer optim by corrHLfit
                                    )
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