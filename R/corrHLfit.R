.check_args_corrHLfit <- function(...,ranFix=list(),ranPars=NULL,fixed=NULL, init.corrHLfit=list(), lower=list(), upper=list()) {
  mc <- match.call(expand.dots = TRUE)
  if (is.null(ranFix)) mc$ranFix <- list() ## deep reason is that relist(., HLCor$ranPars) will need a list ## corrHLfit-spacific
  ## Preventing confusions
  if (!is.null(mc$ranPars)) {
    stop("incorrect 'ranPars' argument in corrHLfit call. Use ranFix (ranPars is for HLCor only)")
  }
  if (!is.null(mc$fixed)) {
    stop("incorrect 'fixed' argument in corrHLfit call: fixed is for fitme only.")
  }
  if ( ! (is.list(lower) && is.list(upper))) {
    wrongclass <- setdiff(unique(c(class(lower),class(upper))),"list")
    stop(paste("'lower' and 'upper' must be of class list, not",paste(wrongclass,collapse=" or ")))
    ## as.list() would flatten rho vectors
  } # else if ((length(lower) || length(upper)) && ! length(init.corrHLfit)) {
  #   warning("'lower' or 'upper' specifications without matching 'init.corrHLfit' have no effect",immediate. = TRUE)
  # }
  HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                names(formals(mat_sqrt)),names(formals(make_scaled_dist))))  
  dotnames <- setdiff(names(mc)[-1],names(formals(corrHLfit)))
  argcheck <- setdiff(dotnames,HLnames)
  if (length(argcheck)) {
    warning(paste("suspect argument(s) ",paste(argcheck, collapse=",")," in corrHLfit call."))
    if ("offset" %in% argcheck) {
      stop("the offset should be a formula term, not a distinct argument.")
    }
  }
  # 
  return(mc)
}

.preprocess_corrHLfit <- function(formula,data, ## matches minimal call of HLfit
                                       family=gaussian(),
                                       init=list(),
                                       ranFix=list(), ## replaces ranFix
                                       lower=list(),upper=list(),
                                       resid.model=~1,
                                       init.HLfit=list(),
                                       control=list(), ## optim.scale (private), nloptr, refit
                                       control.dist=list(),
                                       method="REML", 
                                       HLmethod=method, ## LRT fns assume HLmethod when they are called and when calling
                                       processed=NULL, 
                                       nb_cores = NULL, # to be used by SEM...
                                       objective=NULL,
                                       ... ) {
  mc <- match.call(expand.dots=TRUE) ## mc including dotlist
  if (is.null(processed)) {
    family <- .checkRespFam(family)
    preprocess_args <- .get_inits_preprocess_args(For="corrHLfit")
    names_nondefault  <- intersect(names(mc),names(preprocess_args)) ## mc including dotlist
    preprocess_args[names_nondefault] <- mc[names_nondefault] 
    preprocess_args$family <- family ## already checked 
    if ( ! is.null(mc$rand.family)) preprocess_args$rand.families <- mc$rand.family ## because preprocess expects $rand.families 
    preprocess_args$predictor <- mc$formula ## because preprocess stll expects $predictor 
    preprocess_args$init <- mc$init.corrHLfit ## because preprocess init
    if ( ! missing(method)) preprocess_args$HLmethod <- method
    if ( identical(family$family,"multi")) {
      ## then data are reformatted as a list. Both HLCor and HLfit can analyse such lists for given corrPars and return the joint likelihood
      ## By contrast HLCor should not fit different corrPars to each data, so it does not lapply("corrHLfit",...)
      ## Rather, it calls preprocess which will construct a list of processed objects, to be used conjointly with the data list.
      ## But then we should not attempt to modify an element of 'pocessed' as if it was a single processed object
      ## We must use setProcessed / getProcessed to access list elements.
      if ( ! inherits(data,"list")) {
        familyargs <- family
        familyargs$family <- NULL
        familyargs$binfamily <- NULL
        ## we need the data list in the corrHLfit envir for the call to .makeCheckGeoMatrices
        preprocess_args$data <- do.call(binomialize,c(list(data=data),familyargs)) ## if data not already binomialized
      }     
    }
    mc$processed <- do.call(.preprocess,preprocess_args,envir=parent.frame(1L))
    pnames <- c("data","family","formula","prior.weights","HLmethod","method","rand.family","control.glm","REMLformula",
                "resid.model", "verbose","distMatrix","uniqueGeo","adjMatrix") 
    for (st in pnames) mc[st] <- NULL 
  } 
  return(mc)  
} 

## wrapper for optimization of HLCor.obj OR (iterateSEMSmooth -> HLCor directly)
corrHLfit <- function(formula,data, ## matches minimal call of HLfit
                       init.corrHLfit=list(),
                       init.HLfit=list(),
                       ranFix=list(), 
                       lower=list(),upper=list(),
                       objective=NULL, ## return value of HLCor.obj for optim calls... FR->FR meaningless for full SEM
                       resid.model=~1, 
                       control.dist=list(),
                      control.corrHLfit=list(), ## optim.scale, optimizer, <Optimizer controls>
                      processed=NULL, ## added 2014/02 for programming purposes
                       family=gaussian(),
                      method="REML", 
                      nb_cores=NULL,
                       ... ## pb est risque de passer des args mvs genre HL.method et non HLmethod...
) {
  assign("spaMM_glm_conv_crit",list(max=-Inf) , envir=environment(spaMM_glm.fit))
  time1 <- Sys.time()
  oricall <- match.call(expand.dots=TRUE) ## mc including dotlist
  #oricall$formula <- .preprocess_formula(formula, env=control.corrHLfit$formula_env) ## Cf comment in .GetValidData_info()
  oricall$control.HLfit <- eval(oricall$control.HLfit, parent.frame()) # to evaluate variables in the formula_env, otherwise there are bugs in waiting 
  mc <- oricall
  mc[[1L]] <- get(".preprocess_formula", asNamespace("spaMM"), inherits=FALSE)  ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  oricall$formula <- mc$formula <- eval(mc,parent.frame()) # 
  mc[[1L]] <- get(".check_args_corrHLfit", asNamespace("spaMM"), inherits=FALSE) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  mc <- eval(mc,parent.frame()) # returns modified call 
  mc[[1L]] <- get(".preprocess_corrHLfit", asNamespace("spaMM"), inherits=FALSE) 
  mc <- eval(mc,parent.frame()) # returns modified call including an element 'processed'
  mc[[1L]] <- get("corrHLfit_body", asNamespace("spaMM"), inherits=FALSE) 
  hlcor <- eval(mc,parent.frame()) 
  .check_conv_glm_reinit()
  if (inherits(hlcor,"HLfitlist")) {
    attr(hlcor,"call") <- oricall
  } else {
    oricall$control.dist <- mc$processed$control_dist ## but never in the corrHLfit_body() call
    hlcor$call <- oricall ## this is a call to corrHLfit()
  }
  lsv <- c("lsv",ls())
  if ( ! is.call(hlcor) ) {
    if ( inherits(hlcor,"HLfitlist") ) {
      attr(hlcor,"how") <- list(fit_time=.timerraw(time1), fnname="corrHLfit", spaMM.version=hlcor[[1L]]$how$spaMM.version)
    } else {
      hlcor$how$fit_time <- .timerraw(time1)
      hlcor$how$fnname <- "corrHLfit"
      hlcor$fit_time <- structure(hlcor$how$fit_time,
                                  message="Please use how(<fit object>)[['fit_time']] to extract this information cleanly.")
    }
  }
  rm(list=setdiff(lsv,"hlcor")) 
  return(hlcor)
}

