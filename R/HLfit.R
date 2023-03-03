
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
                  verbose=c(inner=FALSE),
                  HLmethod="HL(1,1)",
                  method="REML",
                  control.HLfit=list(),
                  control.glm=list(),
                  init.HLfit = list(), 
                  fixed=list(), ## phi, lambda, possibly nu, rho if not in init.HLfit
                  # FR->FR ranFix should be able to handle full phi.object and lambda.object ./.
                  # ./. that could be copied in return value.
                  ranFix,
                  etaFix=list(), ## beta, v_h (or even u_h)
                  prior.weights=NULL, weights.form= NULL, # See comments in Z_variants/prior.weights.R
                  X2X=NULL, #  (implemented, but possibly useless, and not tried... remove with care as this would affect preprocessing. See v4.1.23)
                  processed=NULL
) {
  .spaMM.data$options$xLM_conv_crit <- list(max=-Inf)
  time1 <- Sys.time()
  oricall <- match.call()  ## there is no dots in HLfit
  oricall$control.HLfit <- eval(oricall$control.HLfit, parent.frame()) # to evaluate variables in the formula_env, otherwise there are bugs in waiting 
  if ( ! missing(ranFix)) { 
    oricall$fixed <- ranFix
    oricall$ranFix <- NULL
  }
  oricall$fixed <- eval(oricall$fixed, parent.frame()) # allows modif in post-fit code (cf get_HLCorcall) 
  mc <- oricall
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
    if ( ! is.null(weights.form)) mc[["prior.weights"]] <-  weights.form[[2]] # then prior.weights is evaluated to a 'name' and substitute(prior.weights) (elsewhere) does not try to evaluate the variable
    if ( inherits(data,"list")) {
      ## RUN THIS LOOP and return
      multiHLfit <- function() {
        fitlist <- lapply(data, function(dt){
          if (identical(paste(family[[1L]]),"multi")) 
          mc$family <- family$binfamily ## typically binomial()
          mc$data <- dt
          eval(mc) ## calls HLfit recursively
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
      mc[[1L]] <- get(".preprocess_formula", asNamespace("spaMM"), inherits=FALSE)  ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
      oricall$formula <- mc$formula <- eval(mc,parent.frame()) # 
      preprocess_args <- .get_inits_preprocess_args(For="HLfit")
      names_nondefault  <- intersect(names(mc),names(preprocess_args)) ## mc including dotlist
      preprocess_args[names_nondefault] <- mc[names_nondefault] 
      preprocess_args$family <- family ## checked version of 'family'
      preprocess_args$ranFix <- oricall$fixed ## because preprocess expects ranFix
      if ( ! is.null(mc$rand.family)) preprocess_args$rand.families <- mc$rand.family ## because preprocess expects $rand.families 
      preprocess_args$predictor <- mc$formula ## because preprocess still expects $predictor 
#      preprocess_args$HLmethod <- HLmethod ## forces evaluation
      if ( ! missing(method)) preprocess_args$HLmethod <- method
      mc$processed <- processed <- do.call(.preprocess, preprocess_args, envir=parent.frame(1L))
   }
  } else { ## 'processed' is available
    if (  is.list(processed))  { ## "multiple" processed list 
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
      # mc$processed <- processed
      # HLfit_body() called below
    }
  }
  #
  
  # pnames <- c("data","family","formula","prior.weights", "weights.form", "HLmethod","method","rand.family","control.glm","REMLformula",
  #             "resid.model","verbose")
  pnames <- c("data","family","formula","prior.weights", "weights.form","HLmethod","method","rand.family","control.glm","REMLformula",
              "resid.model", "verbose","ranFix") 
  for (st in pnames) mc[st] <- NULL ## info in processed
  mc[[1L]] <- processed$HLfit_body_fn2
  if (.safe_true(processed[["verbose"]]["getCall"][[1L]])) return(mc) ## returns a call if verbose["getCall"'"] is TRUE or 1
  hlfit <- eval(mc,parent.frame())
  .check_conv_dispGammaGLM_reinit()
  if ( ! is.null(processed$return_only)) {
    return(hlfit)    ########################   R E T U R N   a list with $APHLs
  }
  if ( processed$fitenv$prevmsglength) { # there was output for a phi-resid.model. The fit object may then be printed...
    cat("\n")
    processed$fitenv$prevmsglength <- 0L
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

`HLfit.obj` <- function(ranefParsVec, skeleton, objective=processed$objective, processed, fixed=list(), ...) { ## name of first arg MUST differ from names in dotlist...

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
  
  mc <- match.call(expand.dots=TRUE) 
  HLnames <- names(formals(HLfit))
  HLfit.call <- mc[c(1L,which(names(mc) %in% HLnames))] ## keep the call structure

  ranefParsList <- relist(ranefParsVec, skeleton)
  if ( ! is.null(processed$X_off_fn)) { # beta outer-optimisation
    if ( ! is.null(trBeta <- ranefParsList$trBeta)) { # outer beta
      ranefParsList$trBeta <- NULL
      HLfit.call$etaFix$beta <- .betaInv(trBeta)
    } else if ( ! is.null(beta <- ranefParsList$beta)) { # outer beta
      ranefParsList$beta <- NULL
      HLfit.call$etaFix$beta <- beta
    }
  }
  fixed <- .modify_list(fixed, ranefParsList)
  rpType <- .modify_list(attr(fixed, "type"), attr(skeleton, "type"))
  attr(fixed, "type") <- rpType
  if (processed$augZXy_cond) { 
    hlfit <- eval(call(.spaMM.data$options$augZXy_fitfn, processed=processed, fixed=fixed)) # .HLfit_body_augZXy has only these two arguments
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
      # set a 'prefix' for the line to be printed for each iteration of the phi fit when outer optimization is used for the mean response. 
      # In that case a *distinct line* of the form HLfit for <outer opt pars>: phi fit's iter=<say up to 6>, .phi[1]=... 
      # is written for each call of the outer objfn (=> multi-line output).
      # Currently there is no such 'prefix' for mv (_F I X M E_)                           (?)
      # That would require checking processed$residProcesseds (with -'s') and some further effort.
      urP <- unlist(c(HLfit.call$etaFix['beta'], # handles the outer-beta case (spec. when its the only parameter here)
                      .canonizeRanPars(ranefParsList, corr_info=processed$corr_info,checkComplete=FALSE, rC_transf=.spaMM.data$options$rC_transf)))
      processed$port_env$prefix <- paste0("HLfit for ", paste(signif(urP,6), collapse=" "), ": ")
    } 
    # since there is a $processed, we can call HLfit_body here (with HLnames <- names(formals(HLfit_body))), rather than HLfit
    # The main difference is a more definite selection of arguments in the HLfit_body() call through HLfit()
    # and the call to .check_conv_dispGammaGLM_reinit()
    HLfit.call$fixed <- fixed
    HLfit.call[[1L]] <- processed$HLfit
    hlfit <- eval(HLfit.call)
    resu <- hlfit$APHLs[[objective]]
    if (print_phiHGLM_info && processed$verbose["phifit"]>1L) {
      cat(paste0(objective,"=",resu,"\n"))
      processed$fitenv$prevmsglength <- 0L
    }
    if (objective=="cAIC") resu <- - resu ## for minimization of cAIC (private & experimental)
  } 
  return(resu) #
}