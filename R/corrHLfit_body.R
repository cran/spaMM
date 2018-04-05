corrHLfit_body <- function(processed, ## possibly a list of environments
                           init.corrHLfit=list(),
#                      init.HLfit=list(),
                      ranFix=list(), 
                      lower=list(),upper=list(),
                      #control.dist=list(), ## info in processed
                      control.corrHLfit=list(), ## optim.scale, optimizer, <optimizer controls>
                      nb_cores=NULL,
                      ... ## cf dotnames processing below
) { 

  

  #########################################################
  #
  dotlist <- list(...) ## forces evaluations, which makes programming easier...
  if (is.list(processed)) {
    proc1 <- processed[[1L]]
  } else proc1 <- processed
  verbose <-  proc1$verbose
  HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                names(formals(designL.from.Corr)),names(formals(make_scaled_dist))))  ## cf parallel code in HLCor.obj
  ## fill HLCor.args
  good_dotnames <- intersect(names(dotlist),HLnames) ## those specifically for the called fns as def'd by HLnames
  if (length(good_dotnames)) {
    HLCor.args <- dotlist[good_dotnames]
  } else HLCor.args <- list() 
  ## replace some HLCor.args members  
  if ( is.list(processed)) { ## list of environments
    pnames <- names(processed[[1]])
  } else pnames <- names(processed)
  for (st in pnames) HLCor.args[st] <- NULL 
  # 'processed' may be modified below, then will be copied in HLCor.args (and then removed from this envir for safety) => NO NEED to copy any processed element in HLCor.args now.
  optim.scale <- control.corrHLfit$optim.scale 
  if (is.null(optim.scale)) optim.scale="transformed" ## currently no public alternative
  sparse_precision <- processed$sparsePrecisionBOOL
  # test-Nugget -> different p_v / p_bv compromise (!)
  #
  corr_info <- proc1$corr_info 
  #
  # modify HLCor.args and <>bounds;   
  user_init_optim <- init.corrHLfit 
  corr_types <- corr_info$corr_types
  fixed <- .post_process_fixed(ranFix,corr_types=corr_types)
  init.optim <- .post_process_parlist(init.corrHLfit,corr_types=corr_types)   
  init.HLfit <- processed$init_HLfit #to be modified below ## dhglm uses fitme_body not (fitme-> .preprocess) => dhglm code modifies processed$init_HLfit
  init.corrHLfit <- NaN ## make clear it's not to be used
  spatial_terms <- attr(proc1$ZAlist,'exp_spatial_terms')
  family <- proc1$family
  if (family$family=="COMPoisson") {
    checknu <- suppressWarnings(try(environment(family$aic)$nu,silent=TRUE))
    if (inherits(checknu,"try-error") && is.null(init.optim$COMP_nu)) init.optim$COMP_nu <- 1
  } else if (family$family == "negbin") {
    checktheta <- suppressWarnings(try(environment(family$aic)$shape,silent=TRUE))
    if (inherits(checktheta,"try-error") && is.null(init.optim$NB_shape)) init.optim$NB_shape <- 1
  }
  #
  user.lower <- .post_process_parlist(lower,corr_types)
  user.upper <- .post_process_parlist(upper,corr_types) ## keep user input 
  .check_conflict_init_fixed(fixed,init.optim, "given as element of both 'fixed' and 'init'. Check call.")
  .check_conflict_init_fixed(init.HLfit,init.optim, "given as element of both 'init.HLfit' and 'init'. Check call.") ## has quite poor effect on fits
  moreargs <- .calc_moreargs(processed=processed, # possibly a list of environments -> .calc_range_info -> scans then to compute a mean(nbUnique) 
                             corr_types=corr_types, fixed=fixed, init.optim=init.optim, 
                             control_dist=processed$control_dist, 
                             init.HLfit=init.HLfit, corr_info=corr_info, verbose=verbose, lower=lower, upper=upper)
  inits <- .calc_inits(init.optim=init.optim, init.HLfit=init.HLfit,
                       ranFix=fixed, corr_types=corr_types,
                       moreargs=moreargs,
                       user.lower=user.lower, user.upper=user.upper,
                       optim.scale=optim.scale, 
                       For="corrHLfit"
                       )
  #
  init <- inits$`init` ## keeps all init values, all in untransformed scale
  #  code correct mais opaque pour l'utilisateur: init.optim va être ignore en cas de 1D optimization... (OK in fitme)
  init.optim <- inits$`init.optim` ## subset of all optim estimands, as name implies, and in transformed scale
  init.HLfit <- inits$`init.HLfit` ## subset as name implies 
  #
  ## maximization VIA HLCor.obj
  ## by maxim over (corrpars,phi,lambda, (beta)_[corrpars,phi,lambda])
  ##     if trPhi,trLambda are in the init.optims
  ## or by maxim over (corrpars,(beta,phi,lambda)_corrpars)
  ##     otherwise.
  ################ construct intervals for this maximization
  ## construct default upper and lower values ; on transformed scale by default
  user.lower <- lower; user.upper <- upper ## keep user input 
  if ("lambda" %in% c(names(user.lower),names(user.lower)) 
      && is.null(init$lambda)) {
    stop("'lambda' in 'lower' or 'upper' has no effect if absent from 'init.corrHLfit'.")   
    ## corrHLfit-specific logic: if lambda absent from init.corrHLfit it is not outer optimized and then 'as the message says'.
  }
  ################
  LUarglist <- list(canon.init=init, ## canonical scale
                    init.optim=init.optim, ## transformed scale, used only to initialize lower and upper 
                    user.lower=user.lower,user.upper=user.upper,
                    corr_types=corr_types,
                    ranFix=fixed,
                    optim.scale=optim.scale, 
                    moreargs=moreargs)
  LowUp <- do.call(".makeLowerUpper",LUarglist)
  ## LowUp: a list with elements lower and upper that inherits names from init.optim, must be optim.scale as init.optim is by construction
  lower <- LowUp$lower ## list ! which elements may have length >1 !
  upper <- LowUp$upper ## list !
  HLCor.args$ranPars <- fixed  
  #
  control.dist <- vector("list",length(moreargs))
  for (nam in names(moreargs)) control.dist[[nam]] <- moreargs[[nam]]$control.dist
  HLCor.args$control.dist <- control.dist # lapply(moreargs,getElement,name="control.dist")
  processedHL1 <- proc1$HL[1] ## there's also HLmethod in processed<[[]]>$callargs
  if (!is.null(processedHL1) && processedHL1=="SEM" && length(lower)) {
    optimMethod <- "iterateSEMSmooth"
    if (is.null(proc1$SEMargs$control_pmvnorm$maxpts)) {
      if (length(LowUp$lower)) {
        .assignWrapper(processed,"SEMargs$control_pmvnorm$maxpts <- quote(250L*nobs)") 
      } ## else default visible in SEMbetalambda
    }
  } else optimMethod <- ".new_locoptim"
  HLCor.args$processed <- processed
  ## 
  anyHLCor_obj_args <- HLCor.args
  ## HLCor.obj uses a vector + skeleton
  initvec <- unlist(init.optim)
  anyHLCor_obj_args$skeleton <- structure(init.optim, moreargs=moreargs, 
                                          type=relist(rep("fix",length(initvec)),init.optim),
                                          moreargs=moreargs)
  .assignWrapper(anyHLCor_obj_args$processed,
                   paste("return_only <- \"",anyHLCor_obj_args$processed$objective,"APHLs\"",sep=""))
  if (optimMethod=="iterateSEMSmooth") {   
    ## its names should match the colnames of the data in Krigobj = the  parameters of the likelihood surface. Current code maybe not general.
    loclist <- list(anyHLCor_obj_args=anyHLCor_obj_args,  ## contains $processed
                    LowUp=LowUp,init.corrHLfit=user_init_optim, ## FIXME usage of user_init_optim probably not definitive
                    #preprocess.formal.args=preprocess.formal.args, 
                    control.corrHLfit=control.corrHLfit,
                    verbose=verbose[["iterateSEM"]],
                    nb_cores=nb_cores)
    optr <- .probitgemWrap("iterateSEMSmooth",arglist=loclist, pack="probitgem") # do.call("iterateSEMSmooth",loclist) 
    optPars <- relist(optr$par,init.optim)
    if (!is.null(optPars)) attr(optPars,"method") <-"optimthroughSmooth"
  } else { ## this is also called if length(lower)=0 by  (SEM or not) and optPars is then null
    if (identical(verbose["getCall"][[1L]],TRUE)) {
      optPars <- init.optim
    } else {
      optPars <- .new_locoptim(init.optim, LowUp=LowUp,anyHLCor_obj_args=anyHLCor_obj_args,
                               objfn_locoptim=.objfn_locoptim,
                               control=control.corrHLfit, 
                               user_init_optim=user_init_optim)
    }
  }
  if (!is.null(optPars)) {
    HLCor.args$ranPars <- structure(.modify_list(HLCor.args$ranPars,optPars),
                                    type=.modify_list(attr(HLCor.args$ranPars,"type"),
                                                      relist(rep("outer",length(unlist(optPars))),optPars)),
                                    moreargs=moreargs)
  }
  # if processed is an envir, the following is not local to anyHLCor_obj_args$processed but change processed globally
  .assignWrapper(HLCor.args$processed,"return_only <- NULL") 
  .assignWrapper(HLCor.args$processed,"verbose['warn'] <- TRUE") ## important!
  hlcor <- do.call("HLCor",HLCor.args) ## recomputation post optimization (or only computation, if length(lower)=0)
  #
  if (is.call(hlcor)) { return(hlcor[]) } ## HLCorcall
  #
  # hlcor should have received attr(.,"info.uniqueGeo") from HLCor_body.
  attr(hlcor,"optimInfo") <- list(LUarglist=LUarglist, optim.pars=optPars, objective=HLCor.args$processed$objective)
  if ( ! is.null(optPars)) {
    locoptr <- attr(optPars,"optr")
    if (attr(optPars,"method")=="nloptr") {
      if (locoptr$status<0L) hlcor$warnings$optimMessage <- paste("nloptr() message: ",
                                                                  locoptr$message," (status=",locoptr$status,")",sep="")
    } else if ( attr(optPars,"method")=="optim" ) {
      if (locoptr$convergence) hlcor$warnings$optimMessage <- paste("optim() message: ",locoptr$message,
                                           " (convergence=",locoptr$convergence,")",sep="")
    } else if ( attr(optPars,"method")== "optimthroughSmooth") {
      # provide logL estimate from the smoothing, to be used rather than the hlcor logL :
      logLapp <- optr$value
      attr(logLapp,"method") <- "  logL (smoothed)" 
      hlcor$APHLs$logLapp <- logLapp
    }
  }
  lsv <- c("lsv",ls())
  rm(list=setdiff(lsv,"hlcor")) 
  return(hlcor) ## it's the call which says it was returned by corrHLfit
}

    
