# define save as a generic with method save.HLfit ?

.stripFormula <- function(formula) {
  form <- as.formula(deparse(formula))
  assign("formula",NULL,environment(form)) 
  ## so that nothing remains in the within-.stripFormula environment attached to form
  form
}

stripHLfit <- function(object,...) {
  if (inherits(object,"HLfit")) {
    ######## clean matrices
    for (st in names(object$envir)) object$envir[[st]] <- NULL
    if (object$spaMM.version <= "1.10.3") { ## for later version, this part of code documents some obscure issues.
      ######### family matters ...
      ## family$variance, $simulate, $validmu, $aic, $dev.resids get an environment that seems to be 
      ##    the same as object formula (unless the formula itself has been cleaned)
      ## sometimes, the following works
      if (FALSE) {
        ## it is sufficient to assign NULL to any one of the following to recover most of the memory   
        assign("stats",NULL,environment(object$family$aic))
        ## let's assign NULL to all others, though.
        assign("simfun",NULL,environment(object$family$aic))
        assign("okLinks",NULL,environment(object$family$aic))
        assign("linktemp",NULL,environment(object$family$aic))
        ## but this does not work on phi.object$glm_phi$family (eg first example(HLfit))
        ## The alterantive code, as side effect, works all the time. 
      } else {
        ## Wonderful, but not understood!
        as.list(environment(object$family$aic))
        ## ergo...
        if ( ! is.null(object$phi.object$glm_phi$family)) 
          as.list(environment(object$phi.object$glm_phi$family$aic))
      }
      parent.env(environment(object$family$aic)) <- environment(stats::Gamma)
      if ( ! is.null(object$phi.object$glm_phi$family)) 
        parent.env(environment(object$phi.object$glm_phi$family$aic)) <- environment(stats::Gamma)
      ######### formula matters; ideally already previously cleaned 
      if ( ! is.null(fitcall <- attr(object,"corrHLfitcall"))) {
        attr(object,"corrHLfitcall")$formula <- .stripFormula(fitcall$formula)  
        if ( ! is.null(form <- fitcall$resid.formula)) 
          attr(object,"corrHLfitcall")$resid.formula <- .stripFormula(form)
        if ( ! is.null(form <- fitcall$resid.model$formula)) 
          attr(object,"corrHLfitcall")$resid.model$formula <- .stripFormula(form)
      }
      if ( ! is.null(fitcall <- attr(object,"fitmecall"))) {
        attr(object,"fitmecall")$formula <- .stripFormula(fitcall$formula)  
        if ( ! is.null(form <- fitcall$resid.model$formula)) 
          attr(object,"fitmecall")$resid.model$formula <- .stripFormula(form)
      }
      ######### optr matters 
      attr(attr(object,"optimInfo")$optim.pars,"optr")$eval_f <- NULL
      attr(attr(object,"optimInfo")$optim.pars,"optr")$nloptr_environment <- NULL
    }
  } else stop("object does not inherit from class 'HLfit'.")
  invisible(object) #" do not print it otherwise some of the matrices are recomputed.
}

## The above reduces the size _on disk_, which can be checked by saveSize():
if (TRUE) {
  .saveSize <- function (object,...) {
    tf <- tempfile(fileext = ".RData")
    on.exit(unlink(tf))
    save(object, file = tf,...)
    file.size(tf)
  }
}
