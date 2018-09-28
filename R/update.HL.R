getCall.HLfit <- function(x,...) { ## FIXME ? getCall()$resid.model does not look like a language object (cf DHGLM)
  # The [stats::: !]getCall.default method cannot be called b/c it is not exported from stats.
  # stats::getCall() with only call getCall.HLfit in an infinite recursion.
  #
  if (x$spaMM.version> "2.4.62") {
    return(x$call)
  } else {
    # O L D E R version
    # only one of these call may be present in the object: HLCorcall is removed by fitme and corrHLfit
    if ( ! is.null(call <- attr(x,"fitmecall"))) return(call) 
    if ( ! is.null(call <- attr(x,"HLCorcall"))) return(call) ## eg confint on an HLCor object
    if ( ! is.null(call <- attr(x,"corrHLfitcall"))) return(call) 
    return(x$call) ## this one is the HLfit call
  }
}

## to get a call with the structure of the final HLCorcall in fitme or corrHLfit
## ranFix is mandatory: Do not set a default value, so that one has to think about the correct value.
## Therefore, the original ranFix of the outer_object is replaced, unless it is explicitly set to getCall(object)$ranFix or $fixed... (in confint.HLfit)
## Parameters not in ranFix are set to the initial value of of the optimization call.
##   
get_HLCorcall <- function(outer_object, ## accepts fit object, or call, or list of call arguments
                          fixed, ## see comments above
                          ... # anything needed to overcome promises in the call
                          ) { 
  
  outer_call <- getCall(outer_object) ## corrHLfit/fitme/HLCor/HLfit call
  outer_call$data <- outer_object$data ## removes dependence on promise
  outer_fn <- paste(outer_call[[1L]])
  if (outer_fn=="fitme") {
    outer_call$fixed <- fixed
  } else if (outer_fn=="HLCor") {
    outer_call$ranPars <- fixed
  } else outer_call$ranFix <- fixed
  verbose <- outer_call$verbose
  verbose["getCall"] <- TRUE
  outer_call$verbose <- verbose
  ## compare to update.default, commented in R language Definition.
  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(outer_call)))
    dotlist <- list(...)
    for (a in names(extras)[existing]) outer_call[[a]] <- dotlist[[a]]
    if (any(!existing)) {
      outer_call <- c(as.list(outer_call), dotlist[!existing])
    }
  }
  #
  HLCorcall <- eval(as.call(outer_call)) ## calls outer fn and bypasses any optimization to get the inner call HLCor/HLfit
  HLCorcall$call <- NULL ## $call kept the outer call! 
  if (outer_fn=="HLfit") {
    HLCorcall[[1L]] <- quote(HLfit)
  } else if (outer_fn=="fitme") {
    if (is.null(HLCorcall$ranPars)) {
      HLCorcall[[1L]] <- quote(HLfit)
    } else HLCorcall[[1L]] <- quote(HLCor)
  } else HLCorcall[[1L]] <- quote(HLCor)
  .assignWrapper(HLCorcall$processed,"verbose['getCall'] <- NA")
  return(HLCorcall)
}

update.HLfit <- function (object, formula., ..., evaluate = TRUE) {
  if (is.null(call <- getCall(object))) 
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.)) {
    oriform <- formula.HLfit(object) ## formula from $call, not $predictor
    ## fixme A long time ago I wrote "does not handle etaFix$beta"... 
    if (is.null(data <- extras$data)) data <- object$data ## fortunately keeping more than the variables required in the original formula
    call$formula <- .update_formula(oriform,formula.) 
  }
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call))) ## which to replace and which to add to the call
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]] ## replace
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing]) ## add
      call <- as.call(call)
    }
  }
  if (evaluate) 
    eval(call, parent.frame())
  else call
}

# Fix intercept issue in local def of stats::terms.formula
.fixFormulaObject <- function (object) {
  Terms <- terms(object)
  tmp <- attr(Terms, "term.labels")
  ind <- grep("|", tmp, fixed = TRUE)
  if (length(ind)) 
    tmp[ind] <- paste("(", tmp[ind], ")")
  if (length(ind <- attr(Terms, "offset"))) {
    tmp2 <- as.character(attr(Terms, "variables"))[-1L]
    tmp <- c(tmp, tmp2[ind])
  }
  rhs <- if (length(tmp)) 
    paste(tmp, collapse = " + ")
  else "1"
  if (attr(Terms, "intercept")) rhs <- paste("1 +", rhs) ## opposite logic to R
  if (length(form <- formula(object)) > 2L) {
    res <- formula(paste("lhs ~", rhs))
    res[[2L]] <- form[[2L]]
    res
  }
  else formula(paste("~", rhs))
}

.update_formula <- function (old, new, ...) { 
  C_updateform <- get("C_updateform",asNamespace("stats")) ## not kocher?
  tmp <- do.call(".Call",list(C_updateform, as.formula(old), as.formula(new))) # circumventing RcppExports' kind bureaucracy...  
  ## at some point I started to write another fn where 'tmp' was actually 'out' and that continued as :
  #HLframes <- .HLframes(formula=out, data=old$data) ## design matrix X, Y... 
  #attributes(out) <- attributes(HLframes$fixef_terms)
  ## Was it useful ?
  out <- formula(terms.formula(tmp, simplify = FALSE))
  out <- .fixFormulaObject(out)
  environment(out) <- environment(tmp)
  if (!inherits(out, "formula")) 
    class(out) <- c(oldClass(out), "formula")
  return(out)
}
