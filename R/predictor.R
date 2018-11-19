Predictor <- function(...) {stop("'Predictor' is deprecated; see help('Predictor').")}

.as_predictor <- function (formula, data) {
  if (inherits(formula,"predictor")) { stop("Do not call '.as_predictor' on a predictor object.")}
  #formula <- terms(formula, data=data) ## that would deserve some explanations... 
  formula <- .stripFormula(formula) ## remove environment
  class(formula) <- c("predictor",class(formula))
  return(formula)
}

## modified from getAnywhere(print.formula)
print.predictor <- function (x, showAttr = FALSE, ...) 
{
  .x <- x ## .x is original, returned invisibly
  if (! showAttr) attributes(x) <- NULL
  print.default(unclass(x), ...) ## from print.formula...
  invisible(.x)
}



