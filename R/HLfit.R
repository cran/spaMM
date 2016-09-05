
## glm convention in binomial models : eta, fitted values describes FREQUENCIES
##                                     linkinv(eta) describes frequencies, but we need mu to scale as y in the code...
## but the input response ar COUNTS

# function to set and modify various controls of HLfit etc. Accepts a single argument
setControl <- function(...) {
  if (nargs() == 0) return(NULL)
  temp <- list(...)
  if (length(temp) == 1 && is.null(names(temp))) {
    arg <- temp[[1]]
    if( ! is.character(arg)) stop("invalid argument: ", sQuote(arg))
    res <- switch(arg,verbose=logical(0), ## something  on which further code in this fn can operate
                  stop("Unhandled argument:",arg) ## 
                  )
  } else {
    arg <- names(temp)
    res <- temp[[1]]
  }
  if (arg=="verbose") { ## default values
    if (is.na(res["trace"])) res["trace"] <- FALSE
    if (is.na(res["warn"])) res["warn"] <- TRUE
    if (is.na(res["summary"])) res["summary"] <- FALSE
    if (is.na(res["SEM"])) res["SEM"] <- FALSE 
  }
  return(res)
}

# local fn defs
# attention au piege vite oublié
# locfn1 <- fn() {... use global bc/def'd in global; e.g. mu}
# locfn2 <- fn() {... modif mu; locfn1()}
# => locfn2->locfn1-> lit mu global pas local a locfn2

HLfit <- function(formula,
                  data,family=gaussian(),rand.family=gaussian(), 
                  resid.model= ~ 1, resid.formula ,REMLformula=NULL,
                  verbose=c(warn=TRUE,trace=FALSE,summary=FALSE),
                  HLmethod="HL(1,1)",
                  control.HLfit=list(),
                  control.glm=list(),
                  init.HLfit = list(), 
                  ranFix=list(), ## phi, lambda, possibly nu, rho if not in init.HLfit
                  # FR->FR ranFix should be able to handle full phi.object and lambda.object ./.
                  # ./. that could be copied in return value.
                  etaFix=list(), ## beta, v_h (or even u_h)
                  prior.weights=NULL, ## I must avoid default argument reference as formals(HLfit) serves as a template for calls to preprocess() 
                  processed=NULL
) {
  oricall <- mc <- match.call()  ## there is no dots in HLfit
  if (is.null(processed)) { 
    if (missing(data)) {
      data <- environment(formula)
      warning("It is _strongly_ recommanded to use the 'data' argument\n for any application beyond a single fit (e.g. for predict(), etc.)")
    }
    family <- checkRespFam(family) ## same, family as HLCor argument ?
    ################### create data list if family is multi #################################
    if (identical(family$family,"multi")) {
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
    if ( inherits(data,"list")) {
      ## RUN THIS LOOP and return
      multiHLfit <- function() {
        fitlist <- lapply(data,function(dt){
          locmc <- mc
          if (identical(family$family,"multi")) locmc$family <- family$binfamily
          locmc$data <- dt
          eval(locmc) ## calls HLfit recursively
        })
        liks <- sapply(fitlist,function(v) {unlist(v$APHLs)})
        liks <- apply(liks,1,sum)
        attr(fitlist,"APHLs") <- as.list(liks)
        attr(fitlist,"sortedTypes") <- attr(data,"sortedTypes")
        attr(fitlist,"responses") <- attr(data,"responses")
        class(fitlist) <- c("HLfitlist",class(fitlist))     
        return(fitlist)
      }
      return(multiHLfit())
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
      mc$processed <- do.call(preprocess,preprocess.formal.args,envir=environment(formula))
      # HLfit_body() called below
    }
  } else { ## 'processed' is available
    multiple <- attr(processed,"multiple")
    if ( ( ! is.null(multiple)) && multiple)  { ## "multiple" processed list 
      ## RUN THIS LOOP and return
      fitlist <- lapply(seq_len(length(processed)),function(it){
        locmc <- mc
        locmc$processed <- processed[[it]] ## The data are in processed !
        eval(locmc)
      }) ## a pure list of HLCor objects
      liks <- sapply(fitlist,function(v) {unlist(v$APHLs)})
      liks <- apply(liks,1,sum)
      attr(fitlist,"APHLs") <- as.list(liks)
      class(fitlist) <- c("HLfitlist",class(fitlist)) 
      return(fitlist) ## list of HLfit object + one attribute
    } else { ## there is one processed for a single data set 
      # mc$processed <- processed
      # HLfit_body() called below
    }
  }
  #
  mc$data <- NULL
  mc$family <- NULL
  mc$formula <- NULL
  mc$prior.weights <- NULL
  mc$HLmethod <- NULL ## processed$HL  
  mc$rand.family <- NULL ## processed$rand.families  
  mc$control.glm <- NULL ## processed$control.glm  
  mc$resid.formula <- NULL ## mc$resid.model  
  mc$REMLformula <- NULL ## processed$REMLformula
  mc$resid.model <- NULL ## info in processed
  mc[[1L]] <- quote(spaMM::HLfit_body)
  hlfit <- eval(mc,parent.frame())
  if ( ! is.null(processed$return_only)) {
    return(hlfit)    ########################   R E T U R N   a list with $APHLs
  }
  hlfit$call <- oricall ## potentially used by getCall(object) in update.HL
  return(hlfit)
}

`HLfit.obj` <- function(ranefParsVec,skeleton,objective=processed$objective,traceFileName=NULL,processed,...) { ## name of first arg MUST differ from names in dotlist...
  mc <- match.call(expand.dots=TRUE) ## (1) expand.dots added 11/04/2014 for the multinomial... eval 

  if (is.null(processed)) { stop("Call to HLfit.obj() without a 'processed' argument is invalid") } 

  multiple <- attr(processed,"multiple")
  if ( ( ! is.null(multiple)) && multiple)  { ## "multiple" processed list 
    ## RUN THIS LOOP and return
    fitlist <- lapply(seq_len(length(processed)),function(it){
      locmc <- mc
      locmc[[1L]] <- as.name("HLfit.obj") ## replaces "f" !
      locmc$ranefParsVec <- ranefParsVec ## replaces "arg" !
      locmc$processed <- processed[[it]] ## The data are in processed !
      eval(locmc)
    }) ## a pure list of HLfit objects
    resu <- sum(unlist(fitlist))
    if(mc$verbose["objective"]) {
      unrelist <- unlist(relist(ranefParsVec,skeleton)) ## handles elements of lemgth>1
      cat(paste(names(forGiven),"=",signif(unrelist,6),sep="",collapse=", "),
          ", ",objective,"=",resu,"\n",sep="")
    }
    if (is.character(traceFileName)) {
      verif <- paste("#global:",ranefParsVec,resu) 
      write(verif,file=traceFileName,append=TRUE) ## the file is unlink'ed in corrHLfit()  
    }
    return(resu)
  } else { ## there is one processed for a single data set 
    family <- processed$family
    data <- processed$data
  }
  
  HLnames <- names(formals(HLfit))
  HLfit.call <- mc[c(1,which(names(mc) %in% HLnames))] ## keep the call structure
  HLfit.call[[1L]] <- quote(spaMM::HLfit)
  forGiven <- relist(ranefParsVec,skeleton) ## given values of the optimized variables
  notlambda <- setdiff(names(forGiven),"lambda")
  HLfit.call$ranFix$lambda[names(forGiven$lambda)] <- forGiven$lambda
  HLfit.call$ranFix[notlambda] <- forGiven[notlambda] ## do not wipe out other fixed, non optimized variables
  types <- attr(skeleton,"type")
  attr(HLfit.call$ranFix,"type")[names(types)] <- types
  if (is.character(traceFileName)) {
    if(.spaMM.data$options$TRACE.UNLINK) unlink("HLfit.call*.RData")
    zut <- paste(ranefParsVec,collapse="")  
    save(HLfit.call,file=paste("HLfit.call",zut,".RData",sep="")) ## for replicating the problem
  }
  hlfit <- eval(HLfit.call)
  aphls <- hlfit$APHLs
  resu <- aphls[[objective]]
  if(mc$verbose["objective"]) {
    unrelist <- unlist(forGiven) ## handles elements of lemgth>1
    cat(paste(names(forGiven),"=",signif(unrelist,6),sep="",collapse=", "),
        ", ",objective,"=",resu,"\n",sep="")
  }
  if (is.character(traceFileName)) {
    readable <- unlist(canonizeRanPars(ranPars=forGiven,corr.model="",checkComplete=FALSE)$ranPars) 
    verif <- c(unlist(aphls),readable,ranefParsVec) ## hlfit$phi may be NULL
    write(verif,file=traceFileName,ncolumns=length(verif),append=TRUE) ## the file is unlink'ed in corrHLfit()  
  }
  return(resu) #
}