
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
fixedLRT <- function(  ## interface to spaMMLRT or (for devel only) .corrMM_LRT
  null.formula,formula,data,
  method, 
  HLmethod=method,
  REMLformula=NULL,boot.repl=0,
  control=list(),control.boot=list(),fittingFunction, nb_cores=NULL,
  ...) {  ## since .corrMM_LRT is not doc'ed, REMLformula=NULL,boot.repl=0 cannot go into '...' 
  if (missing(null.formula)) stop("'null.formula' argument is missing, with no default.")
  if (missing(formula)) stop("'formula' argument is missing, with no default.")
  if (missing(data)) stop("'data' argument is missing, with no default.")
  if (missing(fittingFunction)) {
    #if (method!="SEM") 
    # message("Assuming 'fittingFunction' is corrHLfit(), but fitme() could be faster.")
    fittingFunction <- "corrHLfit"
  }
  # construct mc$method from method/HLmethod
  ## see 'lm' code for template
  mc <- match.call(expand.dots = TRUE)
  # method mess
  if (missing(method)) { ## OK for corrHLfit -> HLfit
    if (fittingFunction == "fitme") { ## I could allow HLmethod here but presumably not a good long-term solution
      stop("'method' argument is missing, with no default.")
    } else {
      if (missing(HLmethod)) {
        stop("'method' and 'HLmethod' arguments are missing, with no default. Provide 'method'.")
      } else mc$method <- method <- eval(HLmethod,parent.frame())
    } ##  else fitme expects 'method', not 'HLmethod
  } else if ( ! missing(HLmethod))  stop("Don't use both 'method' and 'HLmethod' arguments.")
  mc$HLmethod <- NULL ## method always overrides HLmethod; HLmethod may be re-created in spaMMLRT()
  #
  if (! is.null(control$profiles)) {
    stop("'fixedLRT' does not allow 'control$profiles'.")
  }
  ## other possible settings, through iterative fits
  ## We had a potential backward compatiblity problem, since the simulation scripts 
  ## for the Ecography paper assume that the package automatically interpret the model as spatial, even if findSpatial returns NULL
  ## and we no longer want such a behaviour
  ## but fixedLRT is not used in these scripts, so it can make a different assumption
  spatial <- .findSpatial(formula)
  if ( ! is.null(spatial)) {
    ## both will use p_v for the optim steps, we need to distinguish whether some REML correction is used in iterative algo :
    if ( method %in% c("ML","PQL/L","SEM") || substr(method,0,2) == "ML") {
      callfn <- ".LRT" # mc[[1L]] <- as.name(".LRT") ## does not (yet) handles well other method's  when eg init.corrHLfit contains lambda
      ## there's no profile etc in spaMMLRT... 
    } else { ## EQL, REPQL or REML variants: profiles then not allowed within .corrMM_LRT!
      # FIXME typos (e.g. "PQLL") are not detected...
      callfn <- ".corrMM_LRT" # mc[[1L]] <- as.name(".corrMM_LRT") ## .corrMM_LRT methods and its options below are best frozen to their v1.0 state
      mc$control<-list(profiles=0,prefits=FALSE) ## default values in call by fixedLRT. corrMM.LRT further has default restarts=TRUE and maxit=1
      mc$control[names(control)] <- control ## overrides with user values
      mc$control.boot <- control.boot ## default values in call by fixedLRT are those of .corrMM_LRT ie prefits=FALSE,profiles=0. We can directly copy user values. 
    }
  } else {
    callfn <- ".LRT" # mc[[1L]] <- as.name(".LRT") 
    ## No profiles, maxit, restarts, prefits
    if (is.null(mc$corrMatrix)) { ## neither explicit spatial nor corrMatrix -> HLfit
        mc$fittingFunction <- "HLfit"
      } else {
      ## corrMatrix -> we need to use HLCor
        mc$fittingFunction <- "HLCor"
      }                 
  }  
  
  # mc[[1L]] must be understood by .eval_boot_replicates and quote() does not work there
  # see further problems with thisFnName <- as.character(sys.call()[[1]])
  # eval with an internal function worked in interacive tests but not in R CMD CHECK...
  do.call(callfn, as.list(mc[-1])) # eval(mc, parent.frame()) ## need to specify an environment where .LRT and fixedLRT args can both be found...
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
  #   outst <- paste(" LR statistic (",object$df," df): ",signif(object$LRTori,3),sep="")    
  #   cat(outst)
  bootInfo <- object$bootInfo
  if (!is.null(bootInfo)) {
    cat(" ======== Bootstrap: ========\n")    
    outst <- paste("Raw simulated p-value: ",signif(object$rawBootLRT$p_value,3),sep="")    
    cat(outst)
    cat(paste("\nBartlett-corrected LR test:\n"))
    print(object$BartBootLRT) ## print data frame
  } 
}

print.fixedLRT <- function(x,...) {
  summary(x,...)
  invisible(x)
}

    
