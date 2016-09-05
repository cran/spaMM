getCall.HLfit <- function(x,...) {
  if ( ! is.null(call <- attr(x,"corrHLfitcall"))) return(call) 
  if ( ! is.null(call <- attr(x,"fitmecall"))) return(call) 
  return(x$call)
  # The [stats::: !]getCall.default method cannot be called b/c it is not exported from stats.
  # stats::getCall() with only call getCall.HLfit in an infinite recursion.
}

get_HLCorcall <- function(outercall) { ## to get a callwith the structure of the final HLCorcall in fitme of corrHLfit
  outercall$verbose["getCall"] <- TRUE
  HLCorcall <- eval(as.call(outercall)) ## calls corrHLfit and bypasses optimization to get the call from within the final HLCor
  HLCorcall[[1L]] <- quote(HLCor)
  HLCorcall$verbose["getCall"] <- NA
  #llc$verbose["getCall"] <- NA ## local, useless
  return(HLCorcall)
}


update.HLfit <- function (object, formula., ..., evaluate = TRUE) {
  if (is.null(call <- getCall(object))) 
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.)) {
    predictor <- formula(getCall(object)) ## formula.default gets formula from $call, not from $predictor; and we must use getCall, not $call
    if (inherits(predictor,"predictor")) {
      form <- update.formula(attr(predictor,"oriFormula"),formula.) ## LOSES ALL ATTRIBUTES 
    } else  form <- update.formula(predictor,formula.) 
    ## !!!! FR->FR does not handle etaFix$beta !!!!
    if (! is.null(findOffset(formula.))) {off <- NULL} else { off <- attr(predictor,"offsetObj")$total }
    predArgs <- list(formula=form,
                     LMatrix=attr(predictor,"LMatrix"),
                     AMatrix=attr(predictor,"AMatrix"),
                     ZALMatrix=attr(predictor,"ZALMatrix"), ## see above: *again* from object$call, not from object$predictor
                     offset=off)
    ## attributes BinDenForm and oriFormula will be reconstructed:
    call$formula <- do.call("Predictor",predArgs) ## reconstructs oriFormula... otherwise we have a predictor without it...
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


#`update.HLfit` <- function(object,formula.,...) {update.HL(object=object,formula.=formula.,...)}
#`update.HLCor` <- function(object,formula.,...) {update.HL(object=object,formula.=formula.,...)}

