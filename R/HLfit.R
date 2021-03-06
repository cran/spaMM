
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
                  resid.model= ~ 1, REMLformula=NULL,
                  verbose=c(trace=FALSE),
                  HLmethod="HL(1,1)",
                  method="REML",
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
  assign("spaMM_glm_conv_crit",list(max=-Inf) , envir=environment(spaMM_glm.fit))
  time1 <- Sys.time()
  oricall <- match.call()  ## there is no dots in HLfit
  if (is.null(processed)) { 
    if (missing(data)) {
      data <- environment(formula)
      warning("It is _strongly_ recommended to use the 'data' argument\n for any application beyond a single fit (e.g. for predict(), etc.)")
    }
    #oricall$formula <- .preprocess_formula(formula, env=control.HLfit$formula_env)
    ################### create data list if family is multi #################################
    family <- .checkRespFam(family)
    #family <- .as_call_family(family) ## same, family as HLCor argument ?
    ## mix of mc$family and family:
    ## only the evaluated multi() family is certain to have a $binResponse and ad $binfamily
    ## error for missing "npos" typically follows from lack of explicit $binResponse 
    if (identical(paste(family[[1L]]),"multi")) { ## test syntax valid for all formats of 'family' (except function...)
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
      mc[[1L]] <- get(".preprocess_formula", asNamespace("spaMM"), inherits=FALSE)  ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
      oricall$formula <- mc$formula <- eval(mc,parent.frame()) # 
      preprocess_args <- .get_inits_preprocess_args(For="HLfit")
      names_nondefault  <- intersect(names(mc),names(preprocess_args)) ## mc including dotlist
      preprocess_args[names_nondefault] <- mc[names_nondefault] 
      preprocess_args$family <- family ## checked version of 'family'
      if ( ! is.null(mc$rand.family)) preprocess_args$rand.families <- mc$rand.family ## because preprocess expects $rand.families 
      preprocess_args$predictor <- mc$formula ## because preprocess still expects $predictor 
#      preprocess_args$HLmethod <- HLmethod ## forces evaluation
      if ( ! missing(method)) preprocess_args$HLmethod <- method
      mc$processed <- do.call(.preprocess, preprocess_args, envir=parent.frame(1L))
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
  
  pnames <- c("data","family","formula","prior.weights","HLmethod","method","rand.family","control.glm","REMLformula",
              "resid.model","verbose")
  for (st in pnames) mc[st] <- NULL ## info in processed
  mc[[1L]] <- get("HLfit_body", asNamespace("spaMM"), inherits=FALSE) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  if (identical(mc$processed[["verbose"]]["getCall"][[1L]],TRUE)) return(mc) ## returns a call if verbose["getCall"'"] is TRUE
  hlfit <- eval(mc,parent.frame())
  .check_conv_glm_reinit()
  if ( ! is.null(processed$return_only)) {
    return(hlfit)    ########################   R E T U R N   a list with $APHLs
  }
  hlfit$call <- oricall ## potentially used by getCall(object) in update.HL
  if ( inherits(hlfit,"HLfitlist") ) {
    attr(hlfit,"how") <- list(fit_time=.timerraw(time1),fnname="HLfit", spaMM.version=hlfit[[1L]]$how$spaMM.version)
  } else {
    hlfit$how$fit_time <- .timerraw(time1)
    hlfit$how$fnname <- "HLfit"
    hlfit$fit_time <- structure(hlfit$how$fit_time,
                                message="Please use how(<fit object>)[['fit_time']] to extract this information cleanly.")
  }
  return(hlfit)
}

`HLfit.obj` <- function(ranefParsVec, skeleton, objective=processed$objective, processed, ranFix=list(), ...) { ## name of first arg MUST differ from names in dotlist...

  if (  is.list(processed))  { ## "multiple" processed list 
    mc <- match.call(expand.dots=TRUE) 
    ## RUN THIS LOOP and return
    fitlist <- lapply(seq_len(length(processed)), function(it){
      locmc <- mc
      locmc[[1L]] <- get("HLfit.obj", asNamespace("spaMM"), inherits=FALSE) # as.name("HLfit.obj") ## replaces "f" !
      locmc$ranefParsVec <- ranefParsVec ## replaces "arg" !
      locmc$processed <- processed[[it]] ## The data are in processed !
      eval(locmc)
    }) ## a pure list of HLfit objects
    resu <- sum(unlist(fitlist))
    return(resu)
  } #else there is one processed for a single data set 

  ranefParsList <- relist(ranefParsVec, skeleton)
  ranFix <- .modify_list(ranFix, ranefParsList)
  rpType <- .modify_list(attr(ranFix, "type"), attr(skeleton, "type"))
  attr(ranFix, "type") <- rpType
  if (processed$augZXy_cond) { 
    hlfit <- eval(call(.spaMM.data$options$augZXy_fitfn, processed=processed, ranFix=ranFix)) # .HLfit_body_augZXy has only these two arguments
    aphls <- hlfit$APHLs
    resu <- aphls[[objective]]
    if (objective=="cAIC") resu <- - resu ## for minimization of cAIC (private & experimental)
    if (resu>processed$augZXy_env$objective) {
      processed$augZXy_env$objective <- resu
      processed$augZXy_env$phi_est <- aphls[["phi_est"]]
    }
  } else {
    print_phiHGLM_info <- ( ! is.null(processed$residProcessed) && processed$verbose["phifit"])  
    if (print_phiHGLM_info) {
      # set a 'prefix' for the line to be printed for each iteration of the phi fit when outer optimization is used for the main response. 
      # In that case a *distinct line* of the form HLfit for <outer opt pars>: phi fit's iter=<say up to 6>, .phi[1]=... 
      # is written for each call of the outer objfn (=> multi-line output).
      # Currently there is no such 'prefix' for mv (_F I X M E_)
      # That would require checking processed$residProcesseds (with -'s') and some further effort.
      urP <- unlist(.canonizeRanPars(ranefParsList, corr_info=processed$corr_info,checkComplete=FALSE, rC_transf=.spaMM.data$options$rC_transf))
      processed$port_env$prefix <- paste0("HLfit for ", paste(signif(urP,6), collapse=" "), ": ")
    } 
    # since there is a $processed, we can call HLfit_body here (with HLnames <- names(formals(HLfit_body))), rather than HLfit
    # The main difference is a more definite selection of arguments in the HLfit_body() call through HLfit()
    # and the call to .check_conv_glm_reinit()
    mc <- match.call(expand.dots=TRUE) 
    HLnames <- names(formals(HLfit))
    HLfit.call <- mc[c(1L,which(names(mc) %in% HLnames))] ## keep the call structure
    HLfit.call$ranFix <- ranFix
    HLfit.call[[1L]] <- get("HLfit", asNamespace("spaMM"), inherits=FALSE) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
    hlfit <- eval(HLfit.call)
    resu <- hlfit$APHLs[[objective]]
    if (print_phiHGLM_info) cat(paste0(objective,"=",resu)) # verbose["phifit"]
    if (objective=="cAIC") resu <- - resu ## for minimization of cAIC (private & experimental)
  } 
  return(resu) #
}