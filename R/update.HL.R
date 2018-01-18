getCall.HLfit <- function(x,...) {
  # only one of these call may be present in the object: HLCorcall is removed by fitme and corrHLfit
  if ( ! is.null(call <- attr(x,"fitmecall"))) return(call) 
  if ( ! is.null(call <- attr(x,"HLCorcall"))) return(call) ## eg confint on an HLCor object
  if ( ! is.null(call <- attr(x,"corrHLfitcall"))) return(call) 
  return(x$call) ## this one is the HLfit call
  # The [stats::: !]getCall.default method cannot be called b/c it is not exported from stats.
  # stats::getCall() with only call getCall.HLfit in an infinite recursion.
}

## to get a call with the structure of the final HLCorcall in fitme of corrHLfit
## ranFix is mandatory: Do not set a default value, so that one has to think about the correct value.
## Therefore, the original ranFix of the outer_object is replaced, unless it is explicitly set to getCall(object)$ranFix or $fixed... (in confint.HLfit)
## Parameters not in ranFix are set to the initial value of of the optimization call.
##   
get_HLCorcall <- function(outer_object, ## accepts fit object, or call, or list of call arguments
                          fixed, ## see comments above
                          ... # anything needed to overcome promises in the call
                          ) { 
  if (inherits(outer_object,"HLfit")) {
    outer_call <- getCall(outer_object) ## gets a corrHLfit/fitme call, => eval it to get HLCor callS
    outer_call$data <- outer_object$data ## removes dependence on promise
    outer_fn <- paste(outer_call[[1L]])
    if (outer_fn=="fitme") {
      outer_call$fixed <- fixed
    } else if (outer_fn=="HLCor") {
      outer_call$ranPars <- fixed
    } else outer_call$ranFix <- fixed
  }
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
  HLCorcall <- eval(as.call(outer_call)) ## calls corrHLfit and bypasses optimization to get the call from within the final HLCor
  HLCorcall[[1L]] <- quote(HLCor)
  .assignWrapper(HLCorcall$processed,"verbose['getCall'] <- NA")
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
    if (! is.null(.findOffset(formula.))) {off <- NULL} else { off <- attr(predictor,"offsetObj")$total }
    predArgs <- list(formula=form,
                     LMatrix=NULL, ## F I X M E argument for Predictor, to be removed ? (modif function Predictor() ?)
                     AMatrix=attr(predictor,"AMatrix"),
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

