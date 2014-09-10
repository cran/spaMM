update.HLfit <-
function (object, formula., ..., evaluate = TRUE) 
{
  if (is.null(call <- getCallHL(object))) 
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.)) {
    predictor <- formula(object) ## get from $predictor, not from $call
    if ("predictor" %in% class(predictor)) {
      form <- update.formula(attr(predictor,"oriFormula"),formula.) ## LOSES ALL ATTRIBUTES 
    } else  form <- update.formula(predictor,formula.) 
    if (! is.null(findOffset(formula.))) {off <- NULL} else { off <- attr(predictor,"offset")}
    predArgs <- list(formula=form,
                     LMatrix=attr(predictor,"LMatrix"),
                     AMatrix=attr(predictor,"AMatrix"),
                     ZALMatrix=attr(predictor,"ZALMatrix"),
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
