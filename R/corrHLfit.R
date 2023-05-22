.check_args_corrHLfit <- function(...,ranFix=list(),ranPars=NULL,fixed=NULL, init.corrHLfit=list(), lower=list(), upper=list()) {
  mc <- match.call(expand.dots = TRUE)
  if (is.null(ranFix)) mc$ranFix <- list() ## deep reason is that relist(., HLCor$ranPars) will need a list ## corrHLfit-spacific
  ## Preventing confusions
  if (!is.null(mc$ranPars)) {
    stop("incorrect 'ranPars' argument in corrHLfit call. Use ranFix (ranPars is for HLCor only)")
  }
  # if (!is.null(mc$fixed)) {
  #   stop("incorrect 'fixed' argument in corrHLfit call: fixed is for fitme only.")
  # }
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
                                       fixed=list(), ## currentl ignored
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
    preprocess_args$ranFix <- mc$fixed ## because preprocess expects ranFix
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
    ## removing all elements that are matched in processed:
    # We should remove all processed arguments, in particular those that go into the 'dotlist", otherwise their promises are evaluated again
    ## which is a waste of time (cf corrMatrix=as_precision(...))
    pnames <- c("data","family","formula","prior.weights", "weights.form","HLmethod","method","rand.family","control.glm","REMLformula",
                "resid.model", "verbose","distMatrix","adjMatrix", "control.dist", "corrMatrix","covStruct") 
    # control.HLfit" "init.HLfit"    "etaFix"  remain.
    for (st in pnames) mc[st] <- NULL 
  } 
  return(mc)  
} 

## wrapper for optimization of HLCor.obj OR (iterateSEMSmooth -> HLCor directly)
corrHLfit <- function(formula,data, ## matches minimal call of HLfit
                      init.corrHLfit=list(),
                      init.HLfit=list(),
                      ranFix, 
                      fixed=list(),
                      lower=list(),upper=list(),
                      objective=NULL, ## return value of HLCor.obj for optim calls... FR->FR meaningless for full SEM
                      resid.model=~1, 
                      control.dist=list(),
                      control.corrHLfit=list(), ## optim.scale, optimizer, <Optimizer controls>
                      processed=NULL, ## added 2014/02 for programming purposes
                      family=gaussian(),
                      method="REML", 
                      nb_cores=NULL,
                      weights.form=NULL,
                      ... ## pb est risque de passer des args mvs genre HL.method et non HLmethod...
) {
  .spaMM.data$options$xLM_conv_crit <- list(max=-Inf)
  time1 <- Sys.time()
  oricall <- match.call(expand.dots=TRUE) ## mc including dotlist
  oricall$control.HLfit <- eval(oricall$control.HLfit, parent.frame()) # to evaluate variables in the formula_env, otherwise there are bugs in waiting
  if ( ! missing(ranFix)) {
    oricall$fixed <- ranFix
    oricall$ranFix <- NULL
  }
  oricall$fixed <- eval(oricall$fixed, parent.frame()) # allows modif in post-fit code (cf get_HLCorcall) 
  mc <- oricall
  #
  if ( ! is.null(weights.form)) {
    mc[["prior.weights"]] <-  weights.form[[2]]
  } else if ("prior.weights" %in%  evalq(names(substitute(...())))) {
    p_weights <- substitute(alist(...))$prior.weights # necessary when prior weights has been passed to fitme 
    # through the '...' of another function. In that case we reconstruct the call argument as if they had not been passed in this way.
    # is user quoted the pw, the str() of the result of the substitute() calls is language quote(...)  ~  doubly quoted stuff... => eval 
    if ( (inherits(p_weights,"call") && p_weights[[1L]] == "quote") ) p_weights <- eval(p_weights)
    mc[["prior.weights"]] <- p_weights
  }
  #
  mc[[1L]] <- get(".preprocess_formula", asNamespace("spaMM"), inherits=FALSE)  ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  oricall$formula <- mc$formula <- eval(mc,parent.frame()) # 
  mc[[1L]] <- get(".check_args_corrHLfit", asNamespace("spaMM"), inherits=FALSE) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  mc <- eval(mc,parent.frame()) # returns modified call 
  # mc[["ranFix"]] <- .preprocess_fixed(ranFix) # useless bc 'Partial constraints on ranCoefs are not handled by this function. Use fitme() instead.'
  mc[[1L]] <- get(".preprocess_corrHLfit", asNamespace("spaMM"), inherits=FALSE) 
  mc <- eval(mc,parent.frame()) # returns modified call including an element 'processed'
  mc[[1L]] <- get("corrHLfit_body", asNamespace("spaMM"), inherits=FALSE) 
  hlcor <- eval(mc,parent.frame()) 
  .check_conv_dispGammaGLM_reinit()
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
    if ( ! is.null(mc$control.HLfit$NbThreads)) .setNbThreads(thr=.spaMM.data$options$NbThreads)
  }
  rm(list=setdiff(lsv,"hlcor")) 
  return(hlcor)
}

