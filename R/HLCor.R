HLCor <- function(formula,
                  data,family=gaussian(),
                  ranPars=NULL, ## all dispersion and correlation params ideally provided through ranPars
                  distMatrix, adjMatrix, corrMatrix, covStruct=NULL,
                  method="REML",
                  verbose=c(trace=FALSE), 
                  control.dist=list(), ## provided by <corrfitme>_body if called through this function. Otherwise processed in not available and control.dist will be preprocessed.
                  ...) { # may contain processed
  assign("spaMM_glm_conv_crit",list(max=-Inf) , envir=environment(spaMM_glm.fit))
  time1 <- Sys.time()
  oricall <- match.call(expand.dots = TRUE) 
  if ( ! is.null(oricall$ranFix)) { ## avoiding user's confusion
    stop("!From HLCor: ranFix found in '...'. Make sure to use ranPars only")
  }
  if (!is.null(oricall$LamFix)) {
    stop("argument 'LamFix' of HLCor is obsolete")
  }
  if (!is.null(oricall$PhiFix)) {
    stop("argument 'PhiFix' of HLCor is obsolete")
  }  
  oricall$control.HLfit <- eval(oricall$control.HLfit, parent.frame()) # to evaluate variables in the formula_env, otherwise there are bugs in waiting 
  # frst steps as in HLFit: (no need to test missing(data) in several functions)
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
    if ( inherits(data,"list")) {
      ## RUN THIS LOOP and return
      fitlist <- lapply(seq_len(length(data)), function(data_it){
        locmc <- oricall
        if (identical(family$family,"multi")) locmc$family <- family$binfamily
        locmc$data <- data[[data_it]]
        locmc$distMatrix <- oricall$distMatrix[[data_it]]
        # locmc$uniqueGeo <- oricall$uniqueGeo[[data_it]] # formal removal of uniqueGeo argument 2021/08/03
        eval(locmc)
      }) ## a pure list of HLCor objects
      liks <- sapply(fitlist, function(v) {unlist(v$APHLs)})
      liks <- apply(liks,1,sum)
      attr(fitlist,"APHLs") <- as.list(liks)
      class(fitlist) <- c("HLfitlist",class(fitlist)) 
      return(fitlist) ## list of HLfit object + one attribute
    } else {## there is a single data set, still without processed
      mc <- oricall
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
      preprocess_args$ranFix <- ranPars ## because preprocess expects ranFix
      mc$processed <- do.call(.preprocess, preprocess_args, envir=parent.frame(1L))
      # HLCor() DOES expect a DISTINCT control.dist argument in a call with a 'processed' argument so we extract it:
      oricall$control.dist <- mc$processed$control_dist ## fix bug 26/12/2018 (wrong name) => v2.5.32.
    }
  } else { ## 'processed' is available
    if (  is.list(processed) )  { ## "multiple" processed list 
      ## RUN THIS LOOP and return
      fitlist <- lapply(seq_len(length(processed)), function(it){
        locmc <- oricall
        locmc$processed <- processed[[it]] ## The data are in processed !
        eval(locmc)
      }) ## a pure list of HLCor objects
      liks <- sapply(fitlist, function(v) {unlist(v$APHLs)})
      liks <- apply(liks,1,sum)
      attr(fitlist,"APHLs") <- as.list(liks)
      class(fitlist) <- c("HLfitlist",class(fitlist)) 
      return(fitlist) ## list of HLfit object + one attribute
    } else { ## there is one processed for a single data set 
      mc <- oricall
      # HLCor_body() called below
    }
  }
  ################# single processed, single data analysis: 
  if (identical(mc$processed[["verbose"]]["getCall"][[1L]],TRUE)) return(oricall) ## returns a call is verbose["getCall"'"] is TRUE
  #
  pnames <- c("data","family","formula","prior.weights","HLmethod","method","rand.family","control.glm","REMLformula",
              "resid.model", "verbose","distMatrix","uniqueGeo","adjMatrix") ## try covStruct too...
  for (st in pnames) mc[st] <- NULL 
  mc[[1L]] <- get("HLCor_body", asNamespace("spaMM"), inherits=FALSE) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  hlcor <- eval(mc,parent.frame())
  .check_conv_glm_reinit()
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
  return(hlcor)
}

