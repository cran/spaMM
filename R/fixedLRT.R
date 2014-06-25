fixedLRT <-
function(null.formula,formula,data,HLmethod,REMLformula=NULL,boot.repl=0,
                        control=list(),control.boot=list(),...) {  ## since corrMM.LRT is not doc'ed, REMLformula=NULL,boot.repl=0 cannot go into '...' 
  if (missing(null.formula)) stop("'null.formula' argument is missing, with no default.")
  if (missing(formula)) stop("'formula' argument is missing, with no default.")
  if (missing(data)) stop("'data' argument is missing, with no default.")
  if (missing(HLmethod)) stop("'HLmethod' argument is missing, with no default.")
  if (! is.null(control$profiles)) {
    stop("'fixedLRT' does not allow 'control$profiles'.")
  }
  ## see 'lm' code for template
  mc <- match.call(expand.dots = TRUE)
  ## other possible settings, through iterative fits
  #  mc$init.corrHLfit$lambda <- NULL
  #  mc$init.corrHLfit$phi <- NULL
  ## we have a potential backward compatiblity problem, since the simulation scripts 
  ## for the Ecography paper assume that the package automatically interpret the model as spatial, even if findSpatial returns NULL
  ## and we no longer want such a behaviour
  ## but fixedLRT is not used in these scripts, so it can make a different assumption
  spatial <- findSpatial(formula)
  if ( ! is.null(spatial)) {
    mc$method <- "corrHLfit" ## useful for spaMMLRT call but not for corrMM.LRT where it is the default
    ## both will use p_v for the optim steps, we need to distinguish whether some REML correction is used in iterative algo :
    if ( HLmethod %in% c("ML","PQL/L","SEM") || substr(HLmethod,0,2) == "ML") {
      mc[[1L]] <- as.name("spaMMLRT") ## does not (yet) handles well other HLmethod's  when eg init.corrHLfit contains lambda
      ## there's no profile etc in spaMMLRT... 
    } else { ## EQL, REPQL or REML variants: profiles then not allowed within corrMM.LRT!
      mc[[1L]] <- as.name("corrMM.LRT") ## corrMM.LRT methods and its options below are best frozen to their v1.0 state
      mc$control<-list(profiles=0,prefits=FALSE) ## default values in call by fixedLRT. corrMM.LRT further has default restarts=TRUE and maxit=1
      mc$control[names(control)]<-control ## overrides with user values
      mc$control.boot <- control.boot ## default values in call by fixedLRT are those of corrMM.LRT ie prefits=FALSE,profiles=0. We can directly copy user values. 
    }
  } else {
    mc[[1L]] <- as.name("spaMMLRT") 
    ## No profiles, maxit, restarts, prefits
    if (is.null(mc$corrMatrix)) { ## neither explicit spatial nor corrMatrix -> HLfit
        mc$method <- "HLfit"
      } else {
      ## corrMatrix -> we need to use HLCor
        mc$method <- "HLCor"
      }                 
  }  
  eval(mc, parent.frame())
}
