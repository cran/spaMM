HLCor <- function(formula,
                  data,family=gaussian(),
                  ranPars=NULL, ## all dispersion and correlation params ideally provided through ranPars
                  distMatrix,uniqueGeo=NULL,adjMatrix,corrMatrix,
                  verbose=c(warn=TRUE,trace=FALSE,summary=FALSE),control.dist=list(),
                  ...) { 
  spaMM.options(spaMM_glm_conv_crit=list(max=-Inf),COMP_maxn_warned=FALSE,COMP_geom_approx_warned=FALSE)
  time1 <- Sys.time()
  oricall <- mc <- match.call(expand.dots = TRUE) 
  if ( ! is.null(mc$ranFix)) { ## avoiding user's confusion
    stop("!From HLCor: ranFix found in '...'. Make sure to use ranPars only")
  }
  if (!is.null(mc$LamFix)) {
    stop("argument LamFix of HLCor is obsolete")
  }
  if (!is.null(mc$PhiFix)) {
    stop("argument PhiFix of HLCor is obsolete")
  }  
  # frst steps as in HLFit: (no need to test missing(data) in several functions)
  if (is.null(processed <- mc$processed)) { ## no 'processed'
    ## FR->FR suggests we should add processed as argument of HLCor...
    family <- checkRespFam(family)
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
      fitlist <- lapply(seq_len(length(data)),function(it){
        locmc <- mc
        if (identical(family$family,"multi")) locmc$family <- family$binfamily
        locmc$data <- data[[it]]
        locmc$distMatrix <- mc$distMatrix[[it]]
        locmc$uniqueGeo <- mc$uniqueGeo[[it]]
        eval(locmc)
      }) ## a pure list of HLCor objects
      liks <- sapply(fitlist,function(v) {unlist(v$APHLs)})
      liks <- apply(liks,1,sum)
      attr(fitlist,"APHLs") <- as.list(liks)
      class(fitlist) <- c("HLfitlist",class(fitlist)) 
      return(fitlist) ## list of HLfit object + one attribute
    } else {## there is a single data set, still without processed
      FHF <- formals(HLfit) ## makes sure about default values 
      names_FHF <- names(FHF)
      if ( ! is.null(mc$resid.formula)) mc$resid.model <- mc$resid.formula
      names_nondefault  <- intersect(names(mc),names_FHF) ## mc including dotlist
      FHF[names_nondefault] <- mc[names_nondefault] ##  full HLfit args
      preprocess.formal.args <- FHF[which(names_FHF %in% names(formals(preprocess)))] 
      preprocess.formal.args$family <- family ## already checked 
      preprocess.formal.args$rand.families <- FHF$rand.family ## because preprocess expects $rand.families 
      preprocess.formal.args$predictor <- FHF$formula ## because preprocess stll expects $predictor 
      preprocess.formal.args$ranFix <- ranPars ## because preprocess expects ranFix
      if (! missing(adjMatrix)) preprocess.formal.args$adjMatrix <- adjMatrix ## because adjMatrix not in formals(HLfit)
      mc$processed <- do.call(preprocess,preprocess.formal.args,envir=parent.frame(1L))
      # HLCor_body() called below
    }
  } else { ## 'processed' is available
    multiple <- attr(processed,"multiple")
    if ( ( ! is.null(multiple)) && multiple)  { ## "multiple" processed list 
      ## RUN THIS LOOP and return
      fitlist <- lapply(seq_len(length(processed)),function(it){
        locmc <- mc
        locmc$processed <- processed[[it]] ## The data are in processed !
        locmc$distMatrix <- distMatrix[[it]] ## but the matrices are not HLfit args hence not in processed ! 
        locmc$uniqueGeo <- uniqueGeo[[it]]
        eval(locmc)
      }) ## a pure list of HLCor objects
      liks <- sapply(fitlist,function(v) {unlist(v$APHLs)})
      liks <- apply(liks,1,sum)
      attr(fitlist,"APHLs") <- as.list(liks)
      class(fitlist) <- c("HLfitlist",class(fitlist)) 
      return(fitlist) ## list of HLfit object + one attribute
    } else { ## there is one processed for a single data set 
      # HLCor_body() called below
    }
  }
  ################# single processed, single data analysis: 
  if (identical(mc$verbose["getCall"][[1L]],TRUE)) return(oricall)
  #
  mc$verbose <- .reformat_verbose(eval(mc$verbose),For="HLCor")
  mc$data <- NULL
  mc$family <- NULL
  mc$formula <- NULL
  mc$prior.weights <- NULL
  mc$HLmethod <- NULL ## processed$HL  
  mc$rand.family <- NULL ## processed$rand.families  
  mc$control.glm <- NULL ## processed$control.glm  
  mc$resid.formula <- NULL ## mc$resid.model  
  mc$REMLformula <- NULL ## processed$REMLformula
  mc[[1L]] <- quote(spaMM::HLCor_body)
  hlcor <- eval(mc,parent.frame())
  if (is.null(processed)) .check_conv_glm_reinit()
  if ( ! is.null(processed$return_only)) {
    return(hlcor)    ########################   R E T U R N   a list with $APHLs
  }
  # FR->FR 10/2016 remis une ligne supp de 1.9.24 (mais ça doit être supprimé par fitme ou corrHLfit)
  attr(hlcor,"HLCorcall") <- oricall ## potentially used by getCall(object) in update.HL ./.
  # ./. and more directly by confint (very convenient)
  if (mc$verbose["HLCorSummary"]) { ## useful in final call from corrHLfit
    summary(hlcor) ## input corr pars have been printed at the beginning...   
  }
  hlcor$fit_time <- .timerraw(time1)
  return(hlcor)
}

