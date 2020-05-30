
.Bartlett_robust <- function(LRTobject,robust=TRUE,verbose=F) {
      bootLRTS <- with(LRTobject,2*(bootreps[,1]-bootreps[,2]))## full -null but can be p_v or p_bv
      # plot(qchisq(ppoints(zut),1),sort(zut)) ## QQplot, MASS p. 108
      filter <- bootLRTS[bootLRTS>-1e-08]
      resu <- list()
      if (robust) {
        ## finds the df of the distribution by robust regression (MM: MASS p. 161) of QQplot, 
        robustMean <- MASS::rlm(sort(filter)~qchisq(ppoints(filter),1)-1,method="MM",maxit=200)$coefficients[[1]] 
        ## [[1]] to remove name which otherwise finishes as a rowname in a subsequent dataframe       
        # plot(qchisq(ppoints(filter),1),sort(filter)) ## QQplot, MASS p. 108
        # points(qchisq(ppoints(filter),1),robustMean * qchisq(ppoints(filter),1),pch=".")   
        if (inherits(robustMean,"try-error")) {
          resu$robustMean <- NA
          warning("problem in computation of robustMean.")
        } else resu$robustMean <- robustMean
      }
      resu$meanPosbootLRT <- mean(filter)
      resu$nPosbootLRT <- length(filter) 
      if (verbose) print(unlist(resu))
      return(resu)
}

## fixedLRT is a safe interface for performing tests as described in Ecography paper.
## it does not allow the profiling procedure in .corrMM_LRT
fixedLRT <- function(  ## interface to .LRT or (for devel only) .corrMM_LRT
  null.formula,formula,data,
  method, 
  HLmethod=method,
  REMLformula=NULL,boot.repl=0,
  control="DEPRECATED",
  control.boot="DEPRECATED",
  fittingFunction, 
  #nb_cores=NULL,
  seed=NULL,
  resp_testfn=NULL, 
  # type="marginal", # not necess sing fixedLRT -> .LRT -> spaMM_boot( type="marginal" hard coded)
  ...) {  ## since .corrMM_LRT is not doc'ed, REMLformula=NULL,boot.repl=0 cannot go into '...' 
  if (missing(null.formula)) stop("'null.formula' argument is missing, with no default.")
  if (missing(formula)) stop("'formula' argument is missing, with no default.")
  if (missing(data)) stop("'data' argument is missing, with no default.")
  mc <- match.call(expand.dots = TRUE)
  # method mess
  if (missing(method)) {
    if (missing(HLmethod)) {
      stop("'method' (or 'HLmethod') argument is missing, with no default. Provide 'method'.")
    } else {
      METHOD <- eval(HLmethod,parent.frame())
    }
  } else if ( ! missing(HLmethod)) {
    stop("Don't use both 'method' and 'HLmethod' arguments.")
  } else METHOD <- eval(method,parent.frame())
  #
  ## other possible settings, through iterative fits
  ## We had a potential backward compatiblity problem, since the simulation scripts 
  ## for the Ecography paper assume that the package automatically interpret the model as spatial, even if findSpatial returns NULL
  ## and we no longer want such a behaviour
  ## but fixedLRT is not used in these scripts, so it can make a different assumption
  spatial <- .findSpatial(formula)
  if ( length(spatial)) { ## not list()
    if (missing(fittingFunction)) { ## guessing:
      if (nrow(data)<300L) {
        mc$fittingFunction <- "corrHLfit" ## leverage computation is fast
        mc$method <- NULL
        mc$HLmethod <- METHOD
      } else {
        mc$fittingFunction <- "fitme"
        mc$HLmethod <- NULL
        mc$method <- METHOD
      }
    }
    ## both will use p_v for the optim steps, we need to distinguish whether some REML correction is used in iterative algo :
    if ( METHOD %in% c("ML","PQL/L","SEM") || substr(METHOD,0,2) == "ML") {
      mc[[1L]] <- get(".LRT", asNamespace("spaMM")) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
    } else { ## EQL, REPQL or REML variants: 
      stop(paste0("Likelihood-ratio test not (or no longer) implemented for method ",METHOD)) # .corrMM_LRT long removed
      # FIXME typos (e.g. "PQLL") are not detected...
    }
  } else {
    mc[[1L]] <- get(".LRT", asNamespace("spaMM")) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
    mc$method <- NULL
    mc$HLmethod <- METHOD
    ## No maxit, restarts
    if (missing(fittingFunction)) {
      if (is.null(mc$corrMatrix)) { ## neither explicit spatial nor corrMatrix -> HLfit
        mc$fittingFunction <- "HLfit"
      } else {
        ## corrMatrix -> we need to use HLCor
        mc$fittingFunction <- "HLCor"
      }                 
    }
  }  
  if (mc$fittingFunction=="fitme") {
    if ( ! is.null(mc$init.corrHLfit)) {
      mc[["init"]] <- mc$init.corrHLfit
      mc["init.corrHLfit"] <- NULL
    }
    fixed <- c(mc$ranFix,mc$ranPars) #,mc$etaFix)
    if ( ! is.null(fixed)) {
      mc[["fixed"]] <- fixed
      mc["ranFix"] <- NULL
      mc["ranPars"] <- NULL
      #mc["etaFix"] <- NULL
    }
  }
  eval(mc, parent.frame())
}

corrMM.LRT <- function(...) stop("'corrMM.LRT' is deprecated since version 2.1.66.\n The up-to-date function is 'fixedLRT' (or 'spaMM:::.LRT' for elaborate experiments)")

summary.fixedLRT <- function(object,verbose=TRUE,...) {
  LRT <- object$basicLRT
  print(object$basicLRT)
  #   if (verbose) {
  #     cat(" ========      'full' model:     ========\n")    
  #     summary(object$fullfit,...) 
  #     cat(" ========      'null' model:     ========\n")    
  #     summary(object$nullfit,...) 
  #     cat(" ======== Likelihood ratio test: ========\n")    
  #   }
  #   outst <- paste0(" LR statistic (",object$df," df): ",signif(object$LRTori,3))    
  #   cat(outst)
  bootInfo <- object$bootInfo
  if (!is.null(bootInfo)) {
    cat(" ======== Bootstrap: ========\n")    
    outst <- paste0("Raw simulated p-value: ",signif(object$rawBootLRT$p_value,3))    
    cat(outst)
    cat(paste("\nBartlett-corrected LR test:\n"))
    print(object$BartBootLRT) ## print data frame
    if ( ! is.null(bootInfo$warnlist)) lapply(bootInfo$warnlist, message)
  } 
}

print.fixedLRT <- function(x,...) {
  summary(x,...)
  invisible(x)
}

    
