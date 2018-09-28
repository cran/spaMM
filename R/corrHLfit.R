## wrapper for optimization of HLCor.obj OR (iterateSEMSmooth -> HLCor directly)
.def_call_corrHLfit_body <- function(formula,data, ## matches minimal call of HLfit
                      init.corrHLfit=list(),
                      init.HLfit=list(),
                      ranFix=list(), 
                      lower=list(),upper=list(),
                      objective=NULL, ## return value of HLCor.obj for optim calls... FR->FR meaningless for full SEM
                      resid.model=~1, 
                      control.dist=list(),
                      control.corrHLfit=list(), ## optim.scale, optimizer, <optimizer controls>
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
                names(formals(.spaMM.data$options$mat_sqrt_fn)),names(formals(make_scaled_dist))))  
  dotnames <- setdiff(names(mc)[-1],names(formals(corrHLfit)))
  argcheck <- setdiff(dotnames,HLnames)
  if (length(argcheck)) warning(paste("suspect argument(s) ",paste(argcheck, collapse=",")," in corrHLfit call."))
  # 
  if (is.null(processed)) {
    family <- .checkRespFam(family)
    #family <- .as_call_family(family)
    FHF <- formals(HLfit) ## makes sure about default values 
    names_FHF <- names(FHF)
    names_nondefault  <- intersect(names(mc),names_FHF) ## mc including dotlist
    FHF[names_nondefault] <- mc[names_nondefault] ##  full HLfit args
    preprocess.formal.args <- FHF[which(names_FHF %in% names(formals(.preprocess)))] 
    preprocess.formal.args$For <- "corrHLfit"
    preprocess.formal.args$family <- family ## already checked 
    preprocess.formal.args$rand.families <- FHF$rand.family ## because preprocess expects $rand.families 
    preprocess.formal.args$predictor <- FHF$formula ## because preprocess stll expects $predictor 
    preprocess.formal.args$objective <- objective ## not useful except perhaps to reproduce results from ori paper. 
    preprocess.formal.args$adjMatrix <- mc$adjMatrix ## because adjMatrix not in formals(HLfit)    #
    preprocess.formal.args$corrMatrix <- mc$corrMatrix ## because corrMatrix not in formals(HLfit)    #
    preprocess.formal.args$covStruct <- mc$covStruct ## because covStruct not in formals(HLfit)    #
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
        preprocess.formal.args$data <- do.call(binomialize,c(list(data=data),familyargs)) ## if data not already binomialized
      }     
    }
    mc$processed <- do.call(.preprocess,preprocess.formal.args,envir=parent.frame(1L))
    #mc$ranFix$ranCoefs <- NULL ## but new ranFix can be added by fitme/corrHLfit
    ## removing all elements that duplicate info in processed: 
    pnames <- c("data","family","formula","prior.weights","HLmethod","rand.family","control.glm","REMLformula",
                "resid.model", "verbose")
    for (st in pnames) mc[st] <- NULL 
  }  
  
  mc[[1L]] <- get("corrHLfit_body", asNamespace("spaMM")) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  return(mc)
}

## wrapper for optimization of HLCor.obj OR (iterateSEMSmooth -> HLCor directly)
.corrHLfit_future <- function(formula,data, ## matches minimal call of HLfit
                      init.corrHLfit=list(),
                      init.HLfit=list(),
                      ranFix=list(), 
                      lower=list(),upper=list(),
                      objective=NULL, ## return value of HLCor.obj for optim calls... FR->FR meaningless for full SEM
                      resid.model=~1, 
                      control.dist=list(),
                      control.corrHLfit=list(), ## optim.scale, optimizer, <optimizer controls>
                      processed=NULL, ## added 2014/02 for programming purposes
                      family=gaussian(),
                      nb_cores=NULL,
                      ... ## pb est risque de passer des args mvs genre HL.method et non HLmethod...
) {
  spaMM.options(spaMM_glm_conv_crit=list(max=-Inf))
  time1 <- Sys.time()
  oricall <- match.call(expand.dots=TRUE) ## mc including dotlist
  oricall$formula <- .stripFormula(formula) ## Cf comment in .getValidData
  mc <- oricall
  mc[[1L]] <- get(".def_call_corrHLfit_body", asNamespace("spaMM")) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  mc <- eval(mc,parent.frame())  
  hlcor <- eval(mc,parent.frame()) 
  .check_conv_glm_reinit()
  if (inherits(hlcor,"HLfitlist")) {
    attr(hlcor,"call") <- oricall
  } else hlcor$call <- oricall ## this is a call to corrHLfit()
  attr(hlcor,"HLCorcall") <- NULL
  lsv <- c("lsv",ls())
  if ( ! identical(paste(family[[1L]]),"multi"))  {
    hlcor$how$fit_time <- .timerraw(time1)
    hlcor$fit_time <- structure(hlcor$how$fit_time,
                                message="Please use how(<fit object>)[['fit_time']] to extract this information cleanly.")
  }
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
                       resid.model=~1, 
                       control.dist=list(),
                      control.corrHLfit=list(), ## optim.scale, optimizer, <Optimizer controls>
                      processed=NULL, ## added 2014/02 for programming purposes
                       family=gaussian(),
                       nb_cores=NULL,
                       ... ## pb est risque de passer des args mvs genre HL.method et non HLmethod...
) {
  spaMM.options(spaMM_glm_conv_crit=list(max=-Inf))
  time1 <- Sys.time()
  oricall <- match.call(expand.dots=TRUE) ## mc including dotlist
  oricall$formula <- .stripFormula(formula) ## Cf comment in .getValidData
  mc <- oricall
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
                names(formals(.spaMM.data$options$mat_sqrt_fn)),names(formals(make_scaled_dist))))  
  dotnames <- setdiff(names(mc)[-1],names(formals(corrHLfit)))
  argcheck <- setdiff(dotnames,HLnames)
  if (length(argcheck)) warning(paste("suspect argument(s) ",paste(argcheck, collapse=",")," in corrHLfit call."))
  # 
  if (is.null(processed)) {
    family <- .checkRespFam(family)
    #family <- .as_call_family(family)
    FHF <- formals(HLfit) ## makes sure about default values 
    names_FHF <- names(FHF)
    names_nondefault  <- intersect(names(mc),names_FHF) ## mc including dotlist
    FHF[names_nondefault] <- mc[names_nondefault] ##  full HLfit args
    preprocess.formal.args <- FHF[which(names_FHF %in% names(formals(.preprocess)))] 
    preprocess.formal.args$For <- "corrHLfit"
    preprocess.formal.args$family <- family ## already checked 
    preprocess.formal.args$rand.families <- FHF$rand.family ## because preprocess expects $rand.families 
    preprocess.formal.args$predictor <- FHF$formula ## because preprocess stll expects $predictor 
    preprocess.formal.args$objective <- objective ## now a non-trivial use: cAIC 
    preprocess.formal.args$adjMatrix <- mc$adjMatrix ## because adjMatrix not in formals(HLfit    #
    preprocess.formal.args$corrMatrix <- mc$corrMatrix ## because corrMatrix not in formals(HLfit)    #
    preprocess.formal.args$covStruct <- mc$covStruct ## because covStruct not in formals(HLfit)    #
    preprocess.formal.args$uniqueGeo <- mc$uniqueGeo ## because uniqueGeo not in formals(HLfit)    #
    preprocess.formal.args$distMatrix <- mc$distMatrix ## because distMatrix not in formals(HLfit)    #
    preprocess.formal.args[["control.dist"]] <- control.dist ## because control.dist not in formals(HLfit)    #
    if ( identical(family$family,"multi")) {
      ## then data are reformatted as a list. Both HLCor and HLfit can analyse such lists for given corrPars and return the joint likelihood
      ## By contrast HLCor should not fit different corrPars to each data, so it does not lapply("corrHLfit",...)
      ## Rather, it calls preprocess which will construct a list of processed objects, to be used conjointly with the data list.
      ## But then we should not attempt to modify an element of 'processed' as if it was a single processed object
      ## We must use setProcessed / getProcessed to access list elements.
      if ( ! inherits(data,"list")) {
        familyargs <- family
        familyargs$family <- NULL
        familyargs$binfamily <- NULL
        ## we need the data list in the corrHLfit envir for the call to .makeCheckGeoMatrices
        preprocess.formal.args$data <- do.call(binomialize,c(list(data=data),familyargs)) ## if data not already binomialized
      }     
    }
    mc$processed <- do.call(.preprocess,preprocess.formal.args,envir=parent.frame(1L))
    pnames <- c("data","family","formula","prior.weights","HLmethod","rand.family","control.glm","REMLformula",
                "resid.model", "verbose","distMatrix","uniqueGeo","adjMatrix") 
    for (st in pnames) mc[st] <- NULL 
  }  
  mc[[1L]] <- get("corrHLfit_body", asNamespace("spaMM")) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  hlcor <- eval(mc,parent.frame()) 
  .check_conv_glm_reinit()
  if (inherits(hlcor,"HLfitlist")) {
    attr(hlcor,"call") <- oricall
  } else hlcor$call <- oricall ## this is a call to corrHLfit()
  #attr(hlcor,"HLCorcall") <- NULL # not necess?
  lsv <- c("lsv",ls())
  if ( ! identical(paste(family[[1L]]),"multi") && ! is.call(hlcor) )  {
    hlcor$how$fit_time <- .timerraw(time1)
    hlcor$fit_time <- structure(hlcor$how$fit_time,
                                message="Please use how(<fit object>)[['fit_time']] to extract this information cleanly.")
  }
  rm(list=setdiff(lsv,"hlcor")) 
  return(hlcor)
}

