HLCor <- function(formula,
                  data,family=gaussian(),
                  fixed=NULL, ## all dispersion and correlation params 
                  ranPars,
                  distMatrix, adjMatrix, corrMatrix, covStruct=NULL,
                  method="REML",
                  verbose=c(inner=FALSE), 
                  control.dist=list(), ## provided by <corrfitme>_body if called through this function. Otherwise processed in not available and control.dist will be preprocessed.
                  weights.form=NULL,
                  ...) { # may contain processed
  .spaMM.data$options$xLM_conv_crit <- list(max=-Inf)
  time1 <- Sys.time()
  oricall <- match.call(expand.dots = TRUE) 
  oricall <- ..n_names2expr(oricall) 
  if ( ! is.null(oricall$ranFix)) { ## avoiding user's confusion
    stop("!From HLCor: 'ranFix' found in '...'. Make sure to use 'fixed' or 'ranPars' only")
  }
  if (!is.null(oricall$LamFix)) {
    stop("argument 'LamFix' of HLCor is obsolete")
  }
  if (!is.null(oricall$PhiFix)) {
    stop("argument 'PhiFix' of HLCor is obsolete")
  }  
  oricall$control.HLfit <- eval(oricall$control.HLfit, parent.frame()) # to evaluate variables in the formula_env, otherwise there are bugs in waiting 
  # frst steps as in HLFit: (no need to test missing(data) in several functions)
  if ( ! missing(ranPars)) { 
    oricall$fixed <- ranPars
    oricall$ranPars <- NULL
  }
  oricall$fixed <- eval(oricall$fixed, parent.frame()) # allows modif in post-fit code (cf get_HLCorcall) 
  mc <- oricall
  if (is.null(processed <- oricall$processed)) { ## no 'processed'
    ## FR->FR suggests we should add processed as argument of HLCor...
    family <- .checkRespFam(family)
    if ( identical(family$family,"multi")) {
      if ( ! inherits(data,"list")) {
        if(family$binfamily$family=="binomial") {
          familyargs <- family
          familyargs$family <- NULL
          familyargs$binfamily <- NULL
          data <- do.call(binomialize,c(list(data=data),familyargs)) ## if data not already binomialized
        }
      }
    }
    #
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
    if ( inherits(data,"list")) {
      ## RUN THIS LOOP and return
      fitlist <- lapply(seq_len(length(data)), function(data_it){
        if (identical(family$family,"multi")) mc$family <- family$binfamily
        mc$data <- data[[data_it]]
        mc$distMatrix <- oricall$distMatrix[[data_it]]
        eval(mc)
      }) ## a pure list of HLCor objects
      liks <- sapply(fitlist, function(v) {unlist(v$APHLs)})
      liks <- apply(liks,1,sum)
      attr(fitlist,"APHLs") <- as.list(liks)
      class(fitlist) <- c("HLfitlist",class(fitlist)) 
      return(fitlist) ## list of HLfit object + one attribute
    } else {## there is a single data set, still without processed
      #mc$formula <- .preprocess_formula(formula, env=eval(oricall$control.HLfit)$formula_env)
      mc[[1L]] <- get(".preprocess_formula", asNamespace("spaMM"), inherits=FALSE)  ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
      oricall$formula <- mc$formula <- eval(mc,parent.frame()) # 
      preprocess_args <- .get_inits_preprocess_args(For="HLCor")
      names_nondefault  <- intersect(names(mc),names(preprocess_args)) ## mc including dotlist
      preprocess_args[names_nondefault] <- mc[names_nondefault] 
      preprocess_args$family <- family ## already checked 
      if ( ! is.null(mc$rand.family)) preprocess_args$rand.families <- mc$rand.family ## because preprocess expects $rand.families 
      preprocess_args$predictor <- mc$formula ## because preprocess stll expects $predictor 
      preprocess_args$init <- mc$init.corrHLfit ## because preprocess init
      if ( ! missing(method)) preprocess_args$HLmethod <- method
      preprocess_args$ranFix <- oricall$fixed ## because preprocess expects ranFix
      mc$processed <- do.call(.preprocess, preprocess_args, envir=parent.frame(1L))
      mc$processed$fitenv <- list2env(list(prevmsglength=0L))
      # HLCor() DOES expect a DISTINCT control.dist argument in a call with a 'processed' argument so we extract it:
      oricall$control.dist <- mc$processed$control_dist ## fix bug 26/12/2018 (wrong name) => v2.5.32.
    }
  } else { ## 'processed' is available
    if (  is.list(processed) )  { ## "multiple" processed list 
      ## RUN THIS LOOP and return
      fitlist <- lapply(seq_len(length(processed)), function(it){
        mc$processed <- processed[[it]] ## The data are in processed !
        eval(mc)
      }) ## a pure list of HLCor objects
      liks <- sapply(fitlist, function(v) {unlist(v$APHLs)})
      liks <- apply(liks,1,sum)
      attr(fitlist,"APHLs") <- as.list(liks)
      class(fitlist) <- c("HLfitlist",class(fitlist)) 
      return(fitlist) ## list of HLfit object + one attribute
    } else { ## there is one processed for a single data set 
      # HLCor_body() called below
    }
  }
  ################# single processed, single data analysis: 
  if (.safe_true(mc$processed[["verbose"]]["getCall"][[1L]])) return(mc) ## returns a call is verbose["getCall"'"] is TRUE or 1
  #
  # In a call from fitme -> fitme_body -> .new_locoptim -> .safe_opt -> nloptr::nloptr -> HLcallfn.obj -> HLCor.obj -> here,
  # the mc elements are "ranPars"      "processed" and possibly "control.dist" ... which correspond to the formals of HLCor_body:
  # processed, ranPars, 
  #  control.dist=list(), "# possibly distinct from processed info bc corrHLfit_body/fitme_body may have modified it" 
  #  and ... for HLfit.
  # Removing control.dist > bug in fixedLRT routine test, which has a rho mapping.
  pnames <- c("data","family","formula","prior.weights", "weights.form","HLmethod","method","rand.family","control.glm","REMLformula",
              "resid.model", "verbose","distMatrix","adjMatrix", "corrMatrix","covStruct", "ranPars") 
  for (st in pnames) mc[st] <- NULL 
  # mc[[1L]] <- get("HLCor_body", asNamespace("spaMM"), inherits=FALSE) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  mc[[1L]] <- mc$processed$HLCor_body
  hlcor <- eval(mc,parent.frame())
  .check_conv_dispGammaGLM_reinit()
  if ( ! is.null(processed$return_only)) {
    return(hlcor)    ########################   R E T U R N   a list with $APHLs
  }
  hlcor$call <- oricall ## potentially used by getCall(object) in update.HL ./., NOT processed ; potentially overwritten by calling function
  # ./. and more directly by confint (very convenient)
  if ( inherits(hlcor,"HLfitlist") ) {
    attr(hlcor,"how") <- list(fit_time=.timerraw(time1),fnname="HLCor", spaMM.version=hlcor[[1L]]$how$spaMM.version)
  } else {
    hlcor$how$fit_time <- .timerraw(time1)
    hlcor$how$fnname <- "HLCor"
    hlcor$fit_time <- structure(hlcor$how$fit_time,
                                message="Please use how(<fit object>)[['fit_time']] to extract this information cleanly.")
  }
  if ( ! is.null(mc$control.HLfit$NbThreads)) .setNbThreads(thr=.spaMM.data$options$NbThreads)
  return(hlcor)
}

