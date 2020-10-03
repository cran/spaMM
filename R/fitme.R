.check_args_fitme <- function(...,HLmethod, ranPars=NULL,ranFix=NULL, fixed=list(),lower=list(),upper=list()) {
  mc <- match.call(expand.dots = TRUE)
  if ( missing(HLmethod)) {
    mc$HLmethod <- mc$method
  } else message("'HLmethod' argument for fitme() may become obsolete: use 'method' instead. ")  
  if (is.null(fixed)) mc$fixed <- list() ## deep reason is that relist(., fixed) will need a list ## fitme-specific
  ## Preventing confusions
  if (!is.null(mc$ranPars)) {
    stop("incorrect 'ranPars' argument in fitme() call. Use 'fixed' (ranPars is for HLCor() only)")
  }
  if (!is.null(mc$ranFix)) {
    stop("incorrect 'ranFix' argument in fitme() call. Use 'fixed' (ranFix is for HLfit() and corrHLfit() only)")
  }
  if ( ! (is.list(lower) && is.list(upper))) {
    wrongclass <- setdiff(unique(c(class(lower),class(upper))),"list")
    stop(paste("'lower' and 'upper' must be of class list, not",paste(wrongclass,collapse=" or ")))
  } ## as.list() would flatten rho vectors
  #
  HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                names(formals(mat_sqrt)),names(formals(make_scaled_dist))))  
  dotnames <- setdiff(names(mc)[-1],names(formals(fitme)))
  argcheck <- setdiff(dotnames,HLnames)
  if (length(argcheck)) {
    warning(paste("suspect argument(s) ",paste(argcheck, collapse=",")," in fitme call."))
    if ("offset" %in% argcheck) {
      stop("the offset should be a formula term, not a distinct argument.")
    }
  }
  # 
  return(mc)
}

.preprocess_fitme <- function(formula,data,
                             family=gaussian(),
                             init=list(),
                             fixed=list(), ## replaces ranFix
                             lower=list(),upper=list(),
                             resid.model=~1,
                             init.HLfit=list(),
                             control=list(), ## optim.scale (private), nloptr, refit
                             control.dist=list(),
                             method="ML", 
                             HLmethod=method, ## LRT fns assume HLmethod when they are called and when calling
                             processed=NULL, 
                             nb_cores = NULL, # to be used by SEM...
                             objective=NULL,
                             ... 
) {
  # Here if e.g. 'data' is a promise, str(data) is OK but eval(mc$data) fails. We can manipulate mc elements 
  # but cannot evaluate (hence, test) them easily in the current envir (nor in a fn called from here) 
  # Hence args that must be tested in a fn called from here cannot be passed through mc 
  # A fortiori, we must call .preprocess() in the current envir (where we created mc), not in a fn called from here
  mc <- match.call(expand.dots = TRUE)
  if (is.null(processed)) {
    family <- .checkRespFam(family) ## beware negbin not shadowed by mgcv::negbin()
    preprocess_args <- .get_inits_preprocess_args(For="fitme")
    names_nondefault  <- intersect(names(mc),names(preprocess_args)) ## mc including dotlist
    preprocess_args[names_nondefault] <- mc[names_nondefault] 
    preprocess_args$family <- family ## checked version of 'family'
    if ( ! is.null(mc$rand.family)) preprocess_args$rand.families <- mc$rand.family ## because preprocess expects $rand.families 
    preprocess_args$predictor <- mc$formula ## because preprocess still expects $predictor 
    preprocess_args$ranFix <- fixed ## because preprocess expects ranFix
    preprocess_args$HLmethod <- HLmethod ## forces evaluation
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
    mc$processed <- do.call(.preprocess, preprocess_args, envir=parent.frame(1L))
    ## removing all elements that are matched in processed:
    pnames <- c("data","family","formula","prior.weights","HLmethod","method","rand.family","control.glm","REMLformula",
                "resid.model", "verbose","distMatrix","uniqueGeo","adjMatrix", "control.dist") 
    for (st in pnames) mc[st] <- NULL 
  }  
  return(mc)
}



fitme <- function(formula,data, ## matches minimal call of HLfit
                  family=gaussian(),
                  init=list(),
                  fixed=list(), ## replaces ranFix
                  lower=list(),upper=list(),
                  resid.model=~1,
                  init.HLfit=list(),
                  control=list(), ## optim.scale (private), nloptr, refit
                  control.dist=list(),
                  method="ML", 
                  HLmethod=method, ## LRT fns assume HLmethod when they are called and when calling
                  processed=NULL, 
                  nb_cores = NULL, # to be used by SEM...
                  objective=NULL,
                  ... 
) {
  assign("spaMM_glm_conv_crit",list(max=-Inf) , envir=environment(spaMM_glm.fit))
  time1 <- Sys.time()
  oricall <- match.call(expand.dots=TRUE) ## mc including dotlist
  mc <- oricall
  mc[[1L]] <- get(".preprocess_formula", asNamespace("spaMM"))  ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  oricall$formula <- mc$formula <- eval(mc,parent.frame()) # 
  ## : among other effects, forces eval of promise for formula, so re-evaluating the call later will work 
  ## [cf probitgem re-evaluating fitme(form,.....) in eval_smoothtest()]
  mc[[1L]] <- get(".check_args_fitme", asNamespace("spaMM")) 
  mc <- eval(mc,parent.frame()) # 
  mc[[1L]] <- get(".preprocess_fitme", asNamespace("spaMM"))
  mc <- eval(mc,parent.frame()) # returns modified call including an element 'processed'
  mc[[1L]] <- get("fitme_body", asNamespace("spaMM")) 
  hlcor <- eval(mc,parent.frame()) 
  .check_conv_glm_reinit()
  if (inherits(hlcor,"HLfitlist")) {
    attr(hlcor,"call") <- oricall
  } else {
    oricall$control.dist <- mc$processed$control_dist ## but never in the fitme_body() call
    hlcor$call <- oricall ## this is a call to fitme()
  }
#  attr(hlcor,"HLCorcall") <- NULL # presumably no more needed
  lsv <- c("lsv",ls())
  if ( ! .is.multi(family) && ! is.call(hlcor) ) {
    hlcor$how$fit_time <- .timerraw(time1)
    hlcor$how$fnname <- "fitme"
    hlcor$fit_time <- structure(hlcor$how$fit_time,
                                message="Please use how(<fit object>)[['fit_time']] to extract this information cleanly.")
  }
  rm(list=setdiff(lsv,"hlcor")) ## empties the whole local envir except the return value
  return(hlcor)
}

