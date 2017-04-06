## wrapper for optimization of HLCor.obj OR (iterateSEMSmooth -> HLCor directly)
.def_call_corrHLfit_body <- function(formula,data, ## matches minimal call of HLfit
                      init.corrHLfit=list(),
                      init.HLfit=list(),
                      ranFix=list(), 
                      lower=list(),upper=list(),
                      objective=NULL, ## return value of HLCor.obj for optim calls... FR->FR meaningless for full SEM
                      resid.model=~1, resid.formula,
                      control.dist=list(),
                      control.corrHLfit=list(), ## optim.scale, Optimizer, optimizer.args, maxIter, maxcorners, precision
                      processed=NULL, ## added 2014/02 for programming purposes
                      family=gaussian(),
                      nb_cores=NULL,
                      ... ## pb est risque de passer des args mvs genre HL.method et non HLmethod...
) {
  mc <- match.call(expand.dots=TRUE) ## mc including dotlist
  if (!is.null(mc$ranPars)) {
    stop("incorrect 'ranPars' argument in corrHLfit call. Use ranFix (ranPars is for HLCor only)")
  }
  if (!is.null(mc$fixed)) {
    stop("incorrect 'fixed' argument in corrHLfit call: fixed is for fitme only.")
  }
  if ( ! (is.list(lower) && is.list(upper))) {
    wrongclass <- setdiff(unique(c(class(lower),class(upper))),"list")
    stop(paste("'lower' and 'upper' must be of class list, not",paste(wrongclass,collapse=" or ")))
  } ## as.list() would flatten rho vectors
  #
  HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                names(formals(designL.from.Corr)),names(formals(make_scaled_dist))))  
  dotnames <- setdiff(names(mc)[-1],names(formals(corrHLfit)))
  argcheck <- setdiff(dotnames,HLnames)
  if (length(argcheck)>0) warning(paste("suspect argument(s) ",paste(argcheck, collapse=",")," in corrHLfit call."))
  # 
  if (is.null(processed)) {
    family <- checkRespFam(family)
    #family <- as_call_family(family)
    FHF <- formals(HLfit) ## makes sure about default values 
    names_FHF <- names(FHF)
    if ( ! is.null(mc$resid.formula)) mc$resid.model <- mc$resid.formula
    names_nondefault  <- intersect(names(mc),names_FHF) ## mc including dotlist
    FHF[names_nondefault] <- mc[names_nondefault] ##  full HLfit args
    preprocess.formal.args <- FHF[which(names_FHF %in% names(formals(preprocess)))] 
    preprocess.formal.args$family <- family ## already checked 
    preprocess.formal.args$rand.families <- FHF$rand.family ## because preprocess expects $rand.families 
    preprocess.formal.args$predictor <- FHF$formula ## because preprocess stll expects $predictor 
    preprocess.formal.args$objective <- objective ## not useful except perhaps to reproduce results from ori paper. 
    preprocess.formal.args$adjMatrix <- mc$adjMatrix ## because adjMatrix not in formals(HLfit    #
    famfam <- family$family
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
        ## we need the data list in the corrHLfit envir for the call to makeCheckGeoMatrices
        preprocess.formal.args$data <- do.call(binomialize,c(list(data=data),familyargs)) ## if data not already binomialized
      }     
    }
    mc$processed <- do.call(preprocess,preprocess.formal.args,envir=parent.frame(1L))
    mc$verbose <- .reformat_verbose(eval(mc$verbose),For="corrHLfit")
    ## removing all elements that are matched in processed:
    mc$data <- NULL
    mc$family <- NULL
    mc$formula <- NULL ## processed
    mc$prior.weights <- NULL ## processed$prior.weights  
    mc$HLmethod <- NULL ## processed$HL  
    mc$rand.family <- NULL ## processed$rand.families  
    mc$control.glm <- NULL ## processed$control.glm  
    mc$resid.formula <- NULL ## mc$resid.model  
    mc$REMLformula <- NULL ## processed$REMLformula
  }  
  
  mc[[1L]] <- quote(spaMM::corrHLfit_body) 
  return(mc)
}

## wrapper for optimization of HLCor.obj OR (iterateSEMSmooth -> HLCor directly)
.corrHLfit_future <- function(formula,data, ## matches minimal call of HLfit
                      init.corrHLfit=list(),
                      init.HLfit=list(),
                      ranFix=list(), 
                      lower=list(),upper=list(),
                      objective=NULL, ## return value of HLCor.obj for optim calls... FR->FR meaningless for full SEM
                      resid.model=~1, resid.formula,
                      control.dist=list(),
                      control.corrHLfit=list(), ## optim.scale, Optimizer, optimizer.args, maxIter, maxcorners, precision
                      processed=NULL, ## added 2014/02 for programming purposes
                      family=gaussian(),
                      nb_cores=NULL,
                      ... ## pb est risque de passer des args mvs genre HL.method et non HLmethod...
) {
  spaMM.options(spaMM_glm_conv_crit=list(max=-Inf),COMP_maxn_warned=FALSE,COMP_geom_approx_warned=FALSE)
  time1 <- Sys.time()
  oricall <- mc <- match.call(expand.dots=TRUE) ## mc including dotlist
  oricall$formula <- .stripFormula(formula)  
  if ( ! missing(resid.formula)) oricall$resid.formula <- .stripFormula(resid.formula)
  if ( ! is.null(form <- resid.model$formula)) {
    ## oricall$resid.model$formula <- ... if wrong if oricall$resid.model is a promise
    ## we can get info about the promise by using pryr::uneval(resid.model)
    ## class(uneval(resid.model)) may be "bytecode" if the corrHLfit argument was already an evaluated promise
    ## otherwise class(uneval(resid.model)) may be a "call".
    ## if needed, uneval is minimal R code  + the C++ code in promise.cpp: promise_code(resid_model) might suffice.
    resid.model$formula <- .stripFormula(form)
    oricall$resid.model <- resid.model
  }
  mc[[1L]] <- .def_call_corrHLfit_body
  mc <- eval(mc,parent.frame())  
  hlcor <- eval(mc,parent.frame()) 
  .check_conv_glm_reinit()
  attr(hlcor,"corrHLfitcall") <- oricall ## this says the hlcor was returned by corrHLfit
  attr(hlcor,"HLCorcall") <- NULL
  lsv <- c("lsv",ls())
  hlcor$fit_time <- .timerraw(time1)
  rm(list=setdiff(lsv,"hlcor")) 
  #class(hlcor) <- c(class(hlcor),"corrHLfit")
  return(hlcor)
}


## wrapper for optimization of HLCor.obj OR (iterateSEMSmooth -> HLCor directly)
corrHLfit <- function(formula,data, ## matches minimal call of HLfit
                       init.corrHLfit=list(),
                       init.HLfit=list(),
                       ranFix=list(), 
                       lower=list(),upper=list(),
                       objective=NULL, ## return value of HLCor.obj for optim calls... FR->FR meaningless for full SEM
                       resid.model=~1, resid.formula,
                       control.dist=list(),
                       control.corrHLfit=list(), ## optim.scale, Optimizer, optimizer.args, maxIter, maxcorners, precision
                       processed=NULL, ## added 2014/02 for programming purposes
                       family=gaussian(),
                       nb_cores=NULL,
                       ... ## pb est risque de passer des args mvs genre HL.method et non HLmethod...
) {
  spaMM.options(spaMM_glm_conv_crit=list(max=-Inf),COMP_maxn_warned=FALSE,COMP_geom_approx_warned=FALSE)
  time1 <- Sys.time()
  oricall <- mc <- match.call(expand.dots=TRUE) ## mc including dotlist
  oricall$formula <- .stripFormula(formula)  
  if ( ! missing(resid.formula)) oricall$resid.formula <- .stripFormula(resid.formula)
  if ( ! is.null(form <- resid.model$formula)) {
    ## oricall$resid.model$formula <- ... if wrong if oricall$resid.model is a promise
    ## we can get info about the promise by using pryr::uneval(resid.model)
    ## class(uneval(resid.model)) may be "bytecode" if the corrHLfit argument was already an evaluated promise
    ## otherwise class(uneval(resid.model)) may be a "call".
    ## if needed, uneval is minimal R code  + the C++ code in promise.cpp: promise_code(resid_model) might suffice.
    resid.model$formula <- .stripFormula(form)
    oricall$resid.model <- resid.model
  }
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
  } ## as.list() would flatten rho vectors
  #
  HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                names(formals(designL.from.Corr)),names(formals(make_scaled_dist))))  
  dotnames <- setdiff(names(mc)[-1],names(formals(corrHLfit)))
  argcheck <- setdiff(dotnames,HLnames)
  if (length(argcheck)>0) warning(paste("suspect argument(s) ",paste(argcheck, collapse=",")," in corrHLfit call."))
  # 
  if (is.null(processed)) {
    family <- checkRespFam(family)
    #family <- as_call_family(family)
    FHF <- formals(HLfit) ## makes sure about default values 
    names_FHF <- names(FHF)
    if ( ! is.null(mc$resid.formula)) mc$resid.model <- mc$resid.formula
    names_nondefault  <- intersect(names(mc),names_FHF) ## mc including dotlist
    FHF[names_nondefault] <- mc[names_nondefault] ##  full HLfit args
    preprocess.formal.args <- FHF[which(names_FHF %in% names(formals(preprocess)))] 
    preprocess.formal.args$family <- family ## already checked 
    preprocess.formal.args$rand.families <- FHF$rand.family ## because preprocess expects $rand.families 
    preprocess.formal.args$predictor <- FHF$formula ## because preprocess stll expects $predictor 
    preprocess.formal.args$objective <- objective ## not useful except perhaps to reproduce results from ori paper. 
    preprocess.formal.args$adjMatrix <- mc$adjMatrix ## because adjMatrix not in formals(HLfit    #
    famfam <- family$family
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
        ## we need the data list in the corrHLfit envir for the call to makeCheckGeoMatrices
        preprocess.formal.args$data <- do.call(binomialize,c(list(data=data),familyargs)) ## if data not already binomialized
      }     
    }
    mc$processed <- do.call(preprocess,preprocess.formal.args,envir=parent.frame(1L))
    mc$verbose <- .reformat_verbose(eval(mc$verbose),For="corrHLfit")
    ## removing all elements that are matched in processed:
    mc$data <- NULL
    mc$family <- NULL
    mc$formula <- NULL ## processed
    mc$prior.weights <- NULL ## processed$prior.weights  
    mc$HLmethod <- NULL ## processed$HL  
    mc$rand.family <- NULL ## processed$rand.families  
    mc$control.glm <- NULL ## processed$control.glm  
    mc$resid.formula <- NULL ## mc$resid.model  
    mc$REMLformula <- NULL ## processed$REMLformula
  }  
  
  mc[[1L]] <- quote(spaMM::corrHLfit_body) 
  hlcor <- eval(mc,parent.frame()) 
  .check_conv_glm_reinit()
  attr(hlcor,"corrHLfitcall") <- oricall ## this says the hlcor was returned by corrHLfit
  attr(hlcor,"HLCorcall") <- NULL
  lsv <- c("lsv",ls())
  hlcor$fit_time <- .timerraw(time1)
  rm(list=setdiff(lsv,"hlcor")) 
  #class(hlcor) <- c(class(hlcor),"corrHLfit")
  return(hlcor)
}

