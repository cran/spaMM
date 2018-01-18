HLCor <- function(formula,
                  data,family=gaussian(),
                  ranPars=NULL, ## all dispersion and correlation params ideally provided through ranPars
                  distMatrix,uniqueGeo=NULL,adjMatrix,corrMatrix, covStruct=NULL,
                  verbose=c(trace=FALSE),control.dist=list(),
                  ...) { 
  spaMM.options(spaMM_glm_conv_crit=list(max=-Inf))
  time1 <- Sys.time()
  oricall <- match.call(expand.dots = TRUE) 
  if ( ! is.null(oricall$ranFix)) { ## avoiding user's confusion
    stop("!From HLCor: ranFix found in '...'. Make sure to use ranPars only")
  }
  if (!is.null(oricall$LamFix)) {
    stop("argument LamFix of HLCor is obsolete")
  }
  if (!is.null(oricall$PhiFix)) {
    stop("argument PhiFix of HLCor is obsolete")
  }  
  # frst steps as in HLFit: (no need to test missing(data) in several functions)
  if (is.null(processed <- oricall$processed)) { ## no 'processed'
    ## FR->FR suggests we should add processed as argument of HLCor...
    oricall$formula <- .stripFormula(formula)
    if ( ! is.null(oricall$resid.formula)) oricall$resid.model <- oricall$resid.formula
    #
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
        locmc$uniqueGeo <- oricall$uniqueGeo[[data_it]]
        eval(locmc)
      }) ## a pure list of HLCor objects
      liks <- sapply(fitlist, function(v) {unlist(v$APHLs)})
      liks <- apply(liks,1,sum)
      attr(fitlist,"APHLs") <- as.list(liks)
      class(fitlist) <- c("HLfitlist",class(fitlist)) 
      return(fitlist) ## list of HLfit object + one attribute
    } else {## there is a single data set, still without processed
      mc <- oricall
      FHF <- formals(HLfit) ## makes sure about default values 
      names_FHF <- names(FHF)
      #if ( ! is.null(mc$resid.formula)) mc$resid.model <- mc$resid.formula
      names_nondefault  <- intersect(names(mc),names_FHF) ## mc including dotlist
      FHF[names_nondefault] <- mc[names_nondefault] ##  full HLfit args
      preprocess.formal.args <- FHF[which(names_FHF %in% names(formals(.preprocess)))] 
      preprocess.formal.args$For <- "HLCor"
      preprocess.formal.args$family <- family ## already checked 
      preprocess.formal.args$rand.families <- FHF$rand.family ## because preprocess expects $rand.families 
      preprocess.formal.args$predictor <- FHF$formula ## because preprocess stll expects $predictor 
      preprocess.formal.args$ranFix <- ranPars ## because preprocess expects ranFix
      if (! missing(adjMatrix)) preprocess.formal.args$adjMatrix <- adjMatrix ## because adjMatrix not in formals(HLfit)
      if (! missing(corrMatrix)) preprocess.formal.args$corrMatrix <- corrMatrix ## because corrMatrix not in formals(HLfit)    #
      preprocess.formal.args$covStruct <- mc$covStruct ## because covStruct not in formals(HLfit)    #
      preprocess.formal.args$uniqueGeo <- mc$uniqueGeo ## because uniqueGeo not in formals(HLfit)    #
      preprocess.formal.args$distMatrix <- mc$distMatrix ## because distMatrix not in formals(HLfit)    #
      preprocess.formal.args[["control.dist"]] <- control.dist ## because control.dist not in formals(HLfit)    #
      mc$processed <- do.call(.preprocess,preprocess.formal.args,envir=parent.frame(1L))
      oricall$resid.model <- mc$processed$residModel
      #mc$ranPars$ranCoefs <- NULL ## but new ranFix can be added by fitme/corrHLfit
      # HLCor_body() called below
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
  pnames <- c("data","family","formula","prior.weights","HLmethod","rand.family","control.glm","resid.formula","REMLformula",
              "resid.model", "verbose","distMatrix","uniqueGeo","adjMatrix") ## try covStruct too...
  for (st in pnames) mc[st] <- NULL 
  mc[[1L]] <- quote(spaMM::HLCor_body)
  hlcor <- eval(mc,parent.frame())
  .check_conv_glm_reinit()
  if ( ! is.null(processed$return_only)) {
    return(hlcor)    ########################   R E T U R N   a list with $APHLs
  }
  # FR->FR 10/2016 remis une ligne supp de 1.9.24 (mais ça doit être supprimé par fitme ou corrHLfit)
  attr(hlcor,"HLCorcall") <- oricall ## potentially used by getCall(object) in update.HL ./., NOT processed
  # ./. and more directly by confint (very convenient)
  if ( ! identical(paste(family[[1L]]),"multi")) hlcor$fit_time <- .timerraw(time1)
  return(hlcor)
}

