.def_call_fitme_body <- function(formula,data, ## matches minimal call of HLfit
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
                  ... 
) {
  mc <- match.call(expand.dots=TRUE) ## mc including dotlist
  if ( missing(HLmethod)) {
    HLmethod <- method
  } else message("'HLmethod' argument for fitme() may become obsolete: use 'method' instead. ")
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
                names(formals(designL.from.Corr)),names(formals(make_scaled_dist))))  
  dotnames <- setdiff(names(mc)[-1],names(formals(fitme)))
  argcheck <- setdiff(dotnames,HLnames)
  if (length(argcheck)>0) warning(paste("suspect argument(s) ",paste(argcheck, collapse=",")," in fitme call."))
  # 
  if (is.null(processed)) {
    family <- checkRespFam(family)
    FHF <- formals(HLfit) ## makes sure about default values 
    names_FHF <- names(FHF)
    if ( ! is.null(mc$resid.formula)) mc$resid.model <- mc$resid.formula
    names_nondefault  <- intersect(names(mc),names_FHF) ## mc including dotlist
    FHF[names_nondefault] <- mc[names_nondefault] ##  full HLfit args
    preprocess.formal.args <- FHF[which(names_FHF %in% names(formals(preprocess)))] 
    preprocess.formal.args$family <- family ## already checked 
    preprocess.formal.args$rand.families <- FHF$rand.family ## because preprocess expects $rand.families 
    preprocess.formal.args$predictor <- FHF$formula ## because preprocess stll expects $predictor 
    preprocess.formal.args$ranFix <- fixed ## because preprocess expects ranFix
    preprocess.formal.args$adjMatrix <- mc$adjMatrix ## because adjMatrix not in formals(HLfit)
    preprocess.formal.args$HLmethod <- HLmethod ## forces evaluation
    #
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
    mc$formula <- NULL
    mc$prior.weights <- NULL
    mc$HLmethod <- NULL ## processed$HL  
    mc$rand.family <- NULL ## processed$rand.families  
    mc$control.glm <- NULL ## processed$control.glm  
    mc$resid.formula <- NULL ## mc$resid.model  
    mc$REMLformula <- NULL ## processed$REMLformula
  }  
  
  mc[[1L]] <- quote(spaMM::fitme_body) 
  return(mc)
}

.fitme_future <- function(formula,data, ## matches minimal call of HLfit
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
                  ... 
) {
  spaMM.options(spaMM_glm_conv_crit=list(max=-Inf),COMP_maxn_warned=FALSE,COMP_geom_approx_warned=FALSE)
  time1 <- Sys.time()
  oricall <- mc <- match.call(expand.dots=TRUE) ## mc including dotlist
  oricall$formula <- .stripFormula(formula)
  if ( ! is.null(form <- resid.model$formula)) {
    # cf comments in corrHLfit wrt to next lines
    resid.model$formula <- .stripFormula(form)
    oricall$resid.model <- resid.model
  }
  mc[[1L]] <- .def_call_fitme_body
  mc <- eval(mc,parent.frame())  
  hlcor <- eval(mc,parent.frame()) 
  .check_conv_glm_reinit()
  attr(hlcor,"fitmecall") <- oricall ## this says the hlcor was returned by fitme
  attr(hlcor,"HLCorcall") <- NULL
  #class(hlcor) <- c(class(hlcor,"fitme"))
  lsv <- c("lsv",ls())
  hlcor$fit_time <- .timerraw(time1)
  rm(list=setdiff(lsv,"hlcor")) ## empties the whole local envir except the return value
  return(hlcor)
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
                  ... 
) {
  spaMM.options(spaMM_glm_conv_crit=list(max=-Inf),COMP_maxn_warned=FALSE,COMP_geom_approx_warned=FALSE)
  time1 <- Sys.time()
  oricall <- mc <- match.call(expand.dots=TRUE) ## mc including dotlist
  oricall$formula <- .stripFormula(formula)
  if ( ! is.null(form <- resid.model$formula)) {
    # cf comments in corrHLfit wrt to next lines
    resid.model$formula <- .stripFormula(form)
    oricall$resid.model <- resid.model
  }
  ## Preventing confusions
  if ( missing(HLmethod)) {
    HLmethod <- method
  } else message("'HLmethod' argument for fitme() may become obsolete: use 'method' instead. ")
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
                names(formals(designL.from.Corr)),names(formals(make_scaled_dist))))  
  dotnames <- setdiff(names(mc)[-1],names(formals(fitme)))
  argcheck <- setdiff(dotnames,HLnames)
  if (length(argcheck)>0) warning(paste("suspect argument(s) ",paste(argcheck, collapse=",")," in fitme call."))
  # 
  if (is.null(processed)) {
    family <- checkRespFam(family)
    FHF <- formals(HLfit) ## makes sure about default values 
    names_FHF <- names(FHF)
    if ( ! is.null(mc$resid.formula)) mc$resid.model <- mc$resid.formula
    names_nondefault  <- intersect(names(mc),names_FHF) ## mc including dotlist
    FHF[names_nondefault] <- mc[names_nondefault] ##  full HLfit args
    preprocess.formal.args <- FHF[which(names_FHF %in% names(formals(preprocess)))] 
    preprocess.formal.args$family <- family ## already checked 
    preprocess.formal.args$rand.families <- FHF$rand.family ## because preprocess expects $rand.families 
    preprocess.formal.args$predictor <- FHF$formula ## because preprocess stll expects $predictor 
    preprocess.formal.args$ranFix <- fixed ## because preprocess expects ranFix
    preprocess.formal.args$adjMatrix <- mc$adjMatrix ## because adjMatrix not in formals(HLfit)
    preprocess.formal.args$HLmethod <- HLmethod ## forces evaluation
    #
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
    mc$formula <- NULL
    mc$prior.weights <- NULL
    mc$HLmethod <- NULL ## processed$HL  
    mc$rand.family <- NULL ## processed$rand.families  
    mc$control.glm <- NULL ## processed$control.glm  
    mc$resid.formula <- NULL ## mc$resid.model  
    mc$REMLformula <- NULL ## processed$REMLformula
  }  
  
  mc[[1L]] <- quote(spaMM::fitme_body) 
  hlcor <- eval(mc,parent.frame()) 
  .check_conv_glm_reinit()
  attr(hlcor,"fitmecall") <- oricall ## this says the hlcor was returned by fitme
  attr(hlcor,"HLCorcall") <- NULL
  #class(hlcor) <- c(class(hlcor,"fitme"))
  lsv <- c("lsv",ls())
  hlcor$fit_time <- .timerraw(time1)
  rm(list=setdiff(lsv,"hlcor")) ## empties the whole local envir except the return value
  return(hlcor)
}

