.check_args_fitme <- function(...,HLmethod, ranPars=NULL,ranFix=NULL, fixed=list(), init=list(), lower=list(),upper=list(), what_checked="fitme() call") {
  mc <- match.call(expand.dots = TRUE)
  if ( missing(HLmethod)) {
    mc$HLmethod <- mc$method
  } else message(paste0("'HLmethod' argument in ",what_checked," may become obsolete: use 'method' instead. "))  
  if (is.null(fixed)) mc$fixed <- list() ## deep reason is that relist(., fixed) will need a list ## fitme-specific
  ## Preventing confusions
  if (!is.null(mc$ranPars)) {
    stop(paste0("incorrect 'ranPars' argument in ",what_checked,". Use 'fixed' (ranPars is for HLCor() only)"))
  }
  if (!is.null(mc$ranFix)) {
    stop(paste0("incorrect 'ranFix' argument in ",what_checked,". Use 'fixed' (ranFix is for HLfit() and corrHLfit() only)"))
  }
  if ( ! (is.list(lower) && is.list(upper))) {
    wrongclass <- setdiff(unique(c(class(lower),class(upper))),"list")
    stop(paste("'lower' and 'upper' must be of class list, not",paste(wrongclass,collapse=" or ")))
    ## as.list() would flatten rho vectors
  } # else if ((length(lower$phi) || length(upper$phi)) && ! length(init$phi)) {
  #   warning("'lower' or 'upper' specifications without matching 'init' have no effect",immediate. = TRUE)
  # }
  #
  HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                names(formals(mat_sqrt)),names(formals(make_scaled_dist))))  
  dotnames <- setdiff(names(mc)[-1],c(names(formals(fitme)),"what_checked"))
  argcheck <- setdiff(dotnames,HLnames)
  if (length(argcheck)) {
    warning(paste0("suspect argument(s) '",paste(argcheck, sep="'", collapse=","),"' in ",what_checked,"."))
    if ("offset" %in% argcheck) {
      stop("the offset should be a formula term, not a distinct argument.")
    }
  }
  # 
  return(mc)
}

.preprocess_fitme <- function(formula,data,
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
                             nb_cores = NULL, # to be used by SEM...
                             objective=NULL,
                             For="fitme", # alternative is "fitmv"
                             ... 
) {
  # Here if e.g. 'data' is a promise, str(data) is OK but eval(mc$data) fails. We can manipulate mc elements 
  # but cannot evaluate (hence, test) them easily in the current envir (nor in a fn called from here) 
  # Hence args that must be tested in a fn called from here cannot be passed through mc 
  # A fortiori, we must call .preprocess() in the current envir (where we created mc), not in a fn called from here
  mc <- match.call(expand.dots = TRUE)
  if (is.null(processed)) {
    family <- .checkRespFam(family) ## beware negbin not shadowed by mgcv::negbin()
    preprocess_args <- .get_inits_preprocess_args(For=For)
    names_nondefault  <- intersect(names(mc),names(preprocess_args)) ## mc including dotlist
    preprocess_args[names_nondefault] <- mc[names_nondefault] 
    preprocess_args$family <- family ## checked version of 'family'
    if ( ! is.null(mc$rand.family)) preprocess_args$rand.families <- mc$rand.family ## because preprocess expects $rand.families 
    preprocess_args$predictor <- mc$formula ## because preprocess still expects $predictor 
    preprocess_args$ranFix <- fixed ## because preprocess expects ranFix
    preprocess_args$HLmethod <- HLmethod ## forces evaluation
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
        preprocess_args$data <- do.call(binomialize,c(list(data=data),familyargs)) ## if data not already binomialized
      }     
    }
    mc$processed <- do.call(.preprocess, preprocess_args, envir=parent.frame(1L))
    ## removing all elements that are matched in processed:
    # We should remove all processed arguments, in particular those that go into the 'dotlist", otherwise their promises are evaluated again
    ## which is a waste of time (cf corrMatrix=as_precision(...))
    pnames <- c("data","family","formula","prior.weights", "weights.form","HLmethod","method","rand.family","control.glm","REMLformula",
                "resid.model", "verbose","distMatrix","adjMatrix", "control.dist", "corrMatrix","covStruct") 
    # control.HLfit" "init.HLfit"    "etaFix"  remain.
    for (st in pnames) mc[st] <- NULL 
  }  
  return(mc)
}

.preprocess_fixed <- function(fixed) {
  # Whether we can use a chol transfo (or some others) depend on the constraints... so ad hoc code for specific constraints seems justified
  
  if (! is.null(ranCoefs <- fixed$ranCoefs)) {
    for (it in seq_along(ranCoefs)) {
      ranCoef <- ranCoefs[[it]]
      names(ranCoef) <- NULL # In case the users added fancy names to the vector elements, which could break some post-fit name-matching code
      Xi_ncol <- floor(sqrt(length(ranCoef)*2))
      vdiagPos <- cumsum(c(1L,rev(seq(Xi_ncol-1L)+1L))) # diagpos on vector repre of half matrix, not on matrix
      if (is.null(attr(ranCoefs[[it]],"isDiagFamily"))) {
        attr(ranCoefs[[it]],"isDiagFamily") <-  (! anyNA(ranCoef[ - vdiagPos])) & # If all non-diag pos are set by ranCoef, this may be 'is_diag".
          (anyNA(ranCoef[vdiagPos]))
      } # otherwise the user can avoid the ad hoc code for diag family by explitly setting the attribute to FALSE
      attr(ranCoefs[[it]],"Xi_ncol") <- Xi_ncol # used by .constr_ranCoefsInv()
    }
    fixed$ranCoefs <- ranCoefs
  }
  if (! is.null(fixed$lambda) && is.null(names(fixed$lambda))) names(fixed$lambda) <- seq_along(fixed$lambda)
  fixed
}

fitme <- function(formula,data, ## matches minimal call of HLfit
                  family=gaussian(),
                  init=list(),
                  fixed=list(), ## replaces ranFix
                  lower=list(),upper=list(),
                  resid.model=~1,
                  init.HLfit=list(),
                  control=list(), ## optim.scale (private), nloptr, refit... ultimately passed to fitme_body
                  control.dist=list(),
                  method="ML", 
                  HLmethod=method, ## LRT fns assume HLmethod when they are called and when calling
                  processed=NULL, 
                  nb_cores = NULL, # to be used by SEM...
                  objective=NULL,
                  weights.form=NULL,
                  ... # control.HLfit passed through the dots to .preprocess_fitme() -> .preprocess()
) {
  .spaMM.data$options$xLM_conv_crit <- list(max=-Inf)
  time1 <- Sys.time()
  oricall <- match.call(expand.dots=TRUE) ## mc including dotlist
  oricall$"control.HLfit" <- eval(oricall$control.HLfit, parent.frame()) # to evaluate variables in the formula_env, otherwise there are bugs in waiting 
  oricall$fixed <- eval(oricall$fixed, parent.frame()) # allows modif in post-fit code (cf get_HLCorcall) 
  oricall$init <- eval(oricall[["init"]], parent.frame()) # allows modif in post-fit code (cf get_HLCorcall). Better way ? One should be in principle able 
                                                     # to provide the arguments again to the post-fit call => argument o post-fit fn to provide control of eval envir.
                                                     # cf also alternative strategy of trying 3 envirs in lme4:::update.merMod(), incl. sys.frames()[[1]]
  mc <- oricall
  #
  if ( ! is.null(weights.form)) {
    mc[["prior.weights"]] <-  weights.form[[2]]
  } else if ("prior.weights" %in%  evalq(names(substitute(...())))) { # ~ R >= 4.1's ...names() 
                                                                      # substitute(...())  trick found in a Dunlap post to R-devel, 22/05/2020, 20:53
    p_weights <- substitute(alist(...))$prior.weights # necessary when prior weights has been passed to fitme 
    # through the '...' of another function. In that case we reconstruct the call argument as if they had not been passed in this way.
    # is user quoted the pw, the str() of the result of the substitute() calls is language quote(...)  ~  doubly quoted stuff... => eval 
    if ( (inherits(p_weights,"call") && p_weights[[1L]] == "quote") ) p_weights <- eval(p_weights)
    mc[["prior.weights"]] <- p_weights
  }
  #
  mc[[1L]] <- get(".preprocess_formula", asNamespace("spaMM"), inherits=FALSE)  ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  oricall$formula <- mc$formula <- eval(mc,parent.frame()) # 
  ## : among other effects, forces eval of promise for formula, so re-evaluating the call later will work 
  ## [cf probitgem re-evaluating fitme(form,.....) in eval_smoothtest()];
  ## likewise, associates the evaluated formula_env to the formula
  mc[[1L]] <- get(".check_args_fitme", asNamespace("spaMM"), inherits=FALSE) 
  mc <- eval(mc,parent.frame()) # 
  mc[["fixed"]] <- .preprocess_fixed(fixed)
  mc[[1L]] <- get(".preprocess_fitme", asNamespace("spaMM"), inherits=FALSE)
  mc <- eval(mc,parent.frame()) # returns modified call including an element 'processed'
  mc[[1L]] <- get("fitme_body", asNamespace("spaMM"), inherits=FALSE) 
  hlcor <- eval(mc,parent.frame()) 
  .check_conv_dispGammaGLM_reinit()
  if (inherits(hlcor,"HLfitlist")) {
    attr(hlcor,"call") <- oricall
  } else {
    oricall$control.dist <- mc$processed$control_dist ## but never in the fitme_body() call
    hlcor$call <- oricall ## this is a call to fitme()
  }
#  attr(hlcor,"HLCorcall") <- NULL # presumably no more needed
  lsv <- c("lsv",ls())
  if ( ! is.call(hlcor) ) {
    if ( inherits(hlcor,"HLfitlist") ) {
      attr(hlcor,"how") <- list(fit_time=.timerraw(time1),fnname="fitme", spaMM.version=hlcor[[1L]]$how$spaMM.version)
    } else {
      hlcor$how$fit_time <- .timerraw(time1)
      hlcor$how$fnname <- "fitme"
      hlcor$fit_time <- structure(hlcor$how$fit_time,
                                  message="Please use how(<fit object>)[['fit_time']] to extract this information cleanly.")
    }
  }
  rm(list=setdiff(lsv,"hlcor")) ## empties the whole local envir except the return value
  return(hlcor)
}

