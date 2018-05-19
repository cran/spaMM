## next steps are COMPoisson optim, detection des modeles lambda, gestion des lambda multiples
#  et meilleurs default phi, lambda 
fitme_body <- function(processed,
                       init=list(),
#                       init.HLfit=list(),
                       fixed=list(), ## NULL not valid (should be handled in preprocess?)
                       lower=list(),upper=list(),
                       # control.dist=list(), ## info in processed
                       control=list(), ## optimizer, <optimizer controls>, precision
                       ... ## cf dotnames processing below 
) {
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
  needHLCor_specific_args <- (length(unlist(lower$corrPars)) || length(intersect(corr_types,c("Matern","adjacency","AR1","corrMatrix"))))
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
                                    # I write "F I X" as a TAG for this modif type attribute:
                                    type=.modify_list(relist(rep("fix",length(unlist(fixed))),fixed), #attr(fixed,"type"), 
                                                     relist(rep("outer",length(initvec)),init.optim)) )
    } else {
      use_SEM <- (!is.null(processedHL1) && processedHL1=="SEM")
      if (use_SEM) {
        if (is.null(proc1$SEMargs$control_pmvnorm$maxpts)) {
          if (length(LowUp$lower)) {
            .assignWrapper(processed,"SEMargs$control_pmvnorm$maxpts <- quote(250L*nobs)") 
          } ## else default visible in SEMbetalambda
        }
        stop("reimplement iterateSEMSmooth() in fitme() later")
        ## and then we need to put back the code for logL from smoothing at the end of this function
        ## beware of SEM case with length(lower)=0
        refit_info <- control[["refit"]] ## may be a NULL/ NA/ boolean or a list of booleans 
        if (is.null(refit_info) || is.na(refit_info)) refit_info <- FALSE ## alternatives are TRUE or an explicit list, but no longer NULL
      } else {
        optPars <- .new_locoptim(init.optim, 
                                 LowUp, 
                                 control, objfn_locoptim=.objfn_locoptim, 
                                 HLcallfn.obj=HLcallfn.obj, anyHLCor_obj_args=anyHLCor_obj_args, 
                                 user_init_optim=user_init_optim)
        refit_info <- attr(optPars,"refit_info") ## 'derives' from control[["refit"]]
      }
      ranPars_in_refit <- structure(.modify_list(fixed,optPars),
                                   type=.modify_list(relist(rep("fix",length(unlist(fixed))),fixed), #attr(fixed,"type"),
                                                     relist(rep("outer",length(unlist(optPars))),optPars)))
      if ( is.list(refit_info)) {refit_phi <- refit_info$phi} else refit_phi <- refit_info ## result may be NULL in list case
      init_refit <- list()
      if (identical(refit_phi,TRUE) && ! is.null(optPars$trPhi)) {
        init_refit$phi <- .dispInv(optPars$trPhi)
        ranPars_in_refit$trPhi <- NULL
        attr(ranPars_in_refit,"type")$trPhi <- NULL 
      }
      if ( is.list(refit_info)) {refit_lambda <- refit_info$lambda} else refit_lambda <- refit_info
      if (identical(refit_lambda,TRUE) && ! is.null(optPars$trLambda)) {
        lambdapos <- as.integer(names(optPars$trLambda))
        lambda <- rep(NA,max(lambdapos))
        lambda[lambdapos] <- .dispInv(optPars$trLambda)
        init_refit$lambda <- lambda
        ranPars_in_refit$trLambda <- NULL
        attr(ranPars_in_refit,"type")$trLambda <- NULL 
      }
      if (length(init_refit)) HLCor.args$init.HLfit <- .modify_list(HLCor.args$init.HLfit, init_refit)  
    } ## end if ...getCall... else
    #
    # refit_info is list if so provided by user, else typically boolean. An input NA should have been converted to something else (not documented).
    if (needHLCor_specific_args) {
      attr(ranPars_in_refit,"moreargs") <- moreargs 
      HLCor.args$ranPars <- ranPars_in_refit ## variable locally
    } else HLCor.args$ranFix <- ranPars_in_refit ## variable locally  
    
  }
  #
  # if processed is an envir, the following is not local to anyHLCor_obj_args$processed but change processed globally
  .assignWrapper(HLCor.args$processed,"return_only <- NULL") 
  .assignWrapper(HLCor.args$processed,"verbose['warn'] <- TRUE") ## important!
  hlcor <- do.call(HLcallfn,HLCor.args) ## recomputation post optimization, or only computation if length(initvec)=0
  if (is.call(hlcor)) {
    ## then do.call(HLcallfn,HLCor.args) has retuned the call, not the fit. 
    ## see def of get_HLCorcall() for further explanation
    return(hlcor) ## HLCorcall
  }
  # hlcor may have received attr(.,"info.uniqueGeo") from HLCor_body.
  if (length(initvec)) {
    attr(hlcor,"optimInfo") <- list(LUarglist=LUarglist, optim.pars=optPars, 
                                    objective=proc1$objective) ## processed was erased for safety
  }
  ## substantial effect on object size! :
  lsv <- c("lsv",ls())
  rm(list=setdiff(lsv,"hlcor")) 
  ##
  return(hlcor)
}