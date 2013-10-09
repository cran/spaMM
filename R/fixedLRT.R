fixedLRT <-
function(null.predictor,predictor,data,HLmethod,REMLformula=NULL,boot.repl=0,
                        control=list(),control.boot=list(),...) {  ## since corrMM.LRT is not doc'ed, REMLformula=NULL,boot.repl=0 cannot go into '...' 
  if (missing(null.predictor)) stop("'null.predictor' argument is missing, with no default.")
  if (missing(predictor)) stop("'predictor' argument is missing, with no default.")
  if (missing(data)) stop("'data' argument is missing, with no default.")
  if (missing(HLmethod)) stop("'HLmethod' argument is missing, with no default.")
  ## see 'lm' code for template
  mc <- match.call(expand.dots = TRUE)
  mc$trace <- F
  mc$control<-list(profiles=0,REMLfits=FALSE,restarts=TRUE,maxit=1) ## default values
  mc$control[names(control)]<-control
  mc$control.boot <- list(REMLfits=FALSE,optimFits=TRUE,profiles=0) ## default values
  mc$control.boot[names(control.boot)]<-control.boot
  ## other possible settings, through iterative fits
  #  mc$init.corrHLfit$lambda <- NULL
  #  mc$init.corrHLfit$phi <- NULL
  mc[[1L]] <- as.name("corrMM.LRT")
  eval(mc, parent.frame())
}
