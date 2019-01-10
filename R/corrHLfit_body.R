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
                names(formals(.spaMM.data$options$mat_sqrt_fn)),names(formals(make_scaled_dist))))  ## cf parallel code in HLCor.obj
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
  sparse_precision <- proc1$sparsePrecisionBOOL
  # test-Nugget -> different p_v / p_bv compromise (!)
  #
  # modify HLCor.args and <>bounds;   
  user_init_optim <- init.corrHLfit 
  
  optim_blob <- .calc_optim_args(proc1=proc1, processed=processed,
                                 init=init.corrHLfit, fixed=ranFix, lower=lower, upper=upper, 
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
  moreargs <- optim_blob$moreargs
  LUarglist <- optim_blob$LUarglist
  LowUp <- optim_blob$LowUp
  lower <- LowUp$lower ## list ! which elements may have length >1 !
  upper <- LowUp$upper ## list !
  HLCor.args$ranPars <- fixed
  #
  ## maximization VIA HLCor.obj
  ## by maxim over (corrpars,phi,lambda, (beta)_[corrpars,phi,lambda])
  ##     if trPhi,trLambda are in the init.optims
  ## or by maxim over (corrpars,(beta,phi,lambda)_corrpars)
  ##     otherwise.
  #
  control.dist <- vector("list",length(moreargs))
  for (nam in names(moreargs)) control.dist[[nam]] <- moreargs[[nam]]$control.dist
  HLCor.args$control.dist <- control.dist # lapply(moreargs,getElement,name="control.dist")
  processedHL1 <- proc1$HL[1] ## there's also HLmethod in processed<[[]]>$callargs
  if (!is.null(processedHL1) && processedHL1=="SEM" && length(lower)) {
    optimMethod <- "iterateSEMSmooth"
    if (is.null(proc1$SEMargs$control_pmvnorm$maxpts)) {
      .assignWrapper(processed,"SEMargs$control_pmvnorm$maxpts <- quote(250L*nobs)") 
    } ## else default visible in SEMbetalambda
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
                   paste0("return_only <- \"",proc1$objective,"APHLs\""))
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
                                    # I write "F I X" as a TAG for this modif type attribute:
                                    type=.modify_list(relist(rep("fix",length(unlist(HLCor.args$ranPars))),HLCor.args$ranPars), #attr(HLCor.args$ranPars,"type"),
                                                      relist(rep("outer",length(unlist(optPars))),optPars)),
                                    moreargs=moreargs)
  } else {
      HLCor.args$ranPars <- structure(HLCor.args$ranPars,
                                      type = relist(rep("fix", length(unlist(HLCor.args$ranPars))), HLCor.args$ranPars))
  }
  # if processed is an envir, the following is not local to anyHLCor_obj_args$processed but change processed globally
  .assignWrapper(HLCor.args$processed,"return_only <- NULL") 
  .assignWrapper(HLCor.args$processed,"verbose['warn'] <- TRUE") ## important!
  hlcor <- do.call("HLCor",HLCor.args) ## recomputation post optimization (or only computation, if length(lower)=0)
  #
  if (is.call(hlcor)) { return(hlcor[]) } ## HLCorcall
  #
  # hlcor should have received attr(.,"info.uniqueGeo") from HLCor_body.
  attr(hlcor,"optimInfo") <- list(LUarglist=LUarglist, optim.pars=optPars, 
                                  #augZXy_phi_est=NULL, # implicit NULL
                                  objective=proc1$objective)
  if ( ! is.null(optPars)) {
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
  lsv <- c("lsv",ls())
  rm(list=setdiff(lsv,"hlcor")) 
  return(hlcor) ## it's the call which says it was returned by corrHLfit
}

    
