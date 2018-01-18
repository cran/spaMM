
## glm convention in binomial models : eta, fitted values describes FREQUENCIES
##                                     linkinv(eta) describes frequencies, but we need mu to scale as y in the code...
## but the input response ar COUNTS


# local fn defs
# attention au piege vite oublié
# locfn1 <- fn() {... use global bc/def'd in global; e.g. mu}
# locfn2 <- fn() {... modif mu; locfn1()}
# => locfn2->locfn1-> lit mu global pas local a locfn2

HLfit <- function(formula,
                  data,family=gaussian(),rand.family=gaussian(), 
                  resid.model= ~ 1, resid.formula ,REMLformula=NULL,
                  verbose=c(trace=FALSE),
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
  spaMM.options(spaMM_glm_conv_crit=list(max=-Inf))
  time1 <- Sys.time()
  oricall <- match.call()  ## there is no dots in HLfit
  if (is.null(processed)) { 
    oricall$formula <- .stripFormula(formula)
    if ( ! missing(resid.formula)) oricall$resid.model <- resid.formula
    if (missing(data)) {
      data <- environment(formula)
      warning("It is _strongly_ recommanded to use the 'data' argument\n for any application beyond a single fit (e.g. for predict(), etc.)")
    }
    ################### create data list if family is multi #################################
    family <- .checkRespFam(family)
    #family <- .as_call_family(family) ## same, family as HLCor argument ?
    ## mix of mc$family and family:
    ## only the evaluated multi() family is certain to have a $binResponse and ad $binfamily
    ## error for missing "npos" typically follows from lack of explicit $binResponse 
    if (identical(paste(family[[1L]]),"multi")) { ## test syntax valid for all formats of 'family'
      if ( ! inherits(data,"list")) {
        if(family$binfamily$family=="binomial") {
          familyargs <- family ## must be the evaluated multi() 
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
        fitlist <- lapply(data, function(dt){
          locmc <- oricall
          if (identical(paste(family[[1L]]),"multi")) 
            locmc$family <- family$binfamily ## typically binomial()
          locmc$data <- dt
          eval(locmc) ## calls HLfit recursively
        })
        liks <- sapply(fitlist, function(v) {unlist(v$APHLs)})
        liks <- apply(liks,1,sum)
        attr(fitlist,"APHLs") <- as.list(liks)
        attr(fitlist,"sortedTypes") <- attr(data,"sortedTypes")
        attr(fitlist,"responses") <- attr(data,"responses")
        class(fitlist) <- c("HLfitlist",class(fitlist))     
        return(fitlist)
      }
      return(multiHLfit())
    } else {## there is a single data set, still without processed
      mc <- oricall
      FHF <- formals(HLfit) ## makes sure about default values 
      names_FHF <- names(FHF)
      #if ( ! is.null(mc$resid.formula)) mc$resid.model <- mc$resid.formula
      names_nondefault  <- intersect(names(mc),names_FHF) ## mc including dotlist
      FHF[names_nondefault] <- mc[names_nondefault] ##  full HLfit args
      preprocess.formal.args <- FHF[which(names_FHF %in% names(formals(.preprocess)))] 
      preprocess.formal.args$For <- "HLfit"
      preprocess.formal.args$family <- family ## already checked 
      preprocess.formal.args$rand.families <- FHF$rand.family ## because preprocess expects $rand.families 
      preprocess.formal.args$predictor <- FHF$formula ## because preprocess stll expects $predictor 
      mc$processed <- do.call(.preprocess,preprocess.formal.args,envir=parent.frame(1L))
      oricall$resid.model <- mc$processed$residModel
      #mc$ranFix$ranCoefs <- NULL ## but new ranFix can be added by fitme/corrHLfit
      # HLfit_body() called below
    }
  } else { ## 'processed' is available
    if (  is.list(processed))  { ## "multiple" processed list 
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
      # mc$processed <- processed
      # HLfit_body() called below
    }
  }
  #
  pnames <- c("data","family","formula","prior.weights","HLmethod","rand.family","control.glm","resid.formula","REMLformula",
              "resid.model","verbose")
  for (st in pnames) mc[st] <- NULL ## info in processed
  mc[[1L]] <- quote(spaMM::HLfit_body)
  hlfit <- eval(mc,parent.frame())
  .check_conv_glm_reinit()
  if ( ! is.null(processed$return_only)) {
    return(hlfit)    ########################   R E T U R N   a list with $APHLs
  }
  hlfit$call <- oricall ## potentially used by getCall(object) in update.HL
  if ( ! identical(paste(family[[1L]]),"multi")) hlfit$fit_time <- .timerraw(time1)
  return(hlfit)
}

`HLfit.obj` <- function(ranefParsVec,skeleton,objective=processed$objective,processed,...) { ## name of first arg MUST differ from names in dotlist...
  mc <- match.call(expand.dots=TRUE) ## (1) expand.dots added 11/04/2014 for the multinomial... eval 

  if (is.null(processed)) { stop("Call to HLfit.obj() without a 'processed' argument is invalid") } 

  if (  is.list(processed))  { ## "multiple" processed list 
    ## RUN THIS LOOP and return
    fitlist <- lapply(seq_len(length(processed)), function(it){
      locmc <- mc
      locmc[[1L]] <- as.name("HLfit.obj") ## replaces "f" !
      locmc$ranefParsVec <- ranefParsVec ## replaces "arg" !
      locmc$processed <- processed[[it]] ## The data are in processed !
      eval(locmc)
    }) ## a pure list of HLfit objects
    resu <- sum(unlist(fitlist))
    return(resu)
  } else { ## there is one processed for a single data set 
    family <- processed$family
    data <- processed$data
  }
  
  HLnames <- names(formals(HLfit))
  HLfit.call <- mc[c(1,which(names(mc) %in% HLnames))] ## keep the call structure
  HLfit.call[[1L]] <- quote(spaMM::HLfit)
  forGiven <- relist(ranefParsVec,skeleton) ## given values of the optimized variables
  if (spaMM.getOption("wDEVEL2")) {
    parlist <- attr(HLfit.call$ranFix,"parlist")
    parlist <- .merge_parlist(parlist,new=forGiven,types="fix")## consistently with previous code
    attr(HLfit.call$ranFix,"parlist") <- parlist
  }
  notlambda <- setdiff(names(forGiven),"lambda")
  HLfit.call$ranFix$lambda[names(forGiven$lambda)] <- forGiven$lambda
  HLfit.call$ranFix[notlambda] <- forGiven[notlambda] ## do not wipe out other fixed, non optimized variables
  #types <- attr(skeleton,"type")
  #attr(HLfit.call$ranFix,"type")[names(types)] <- types
  hlfit <- eval(HLfit.call)
  aphls <- hlfit$APHLs
  resu <- aphls[[objective]]
  lsv <- c("lsv",ls())
  rm(list=setdiff(lsv,"resu"))
  return(resu) #
}