Predictor <- function(...) {stop("'Predictor' is deprecated; see help('Predictor').")}

## modified from getAnywhere(print.formula)
print.predictor <- function (x, showAttr = FALSE, ...) 
{
  .x <- x ## .x is original, returned invisibly
  if (! showAttr) attributes(x) <- NULL
  print.default(unclass(x), ...) ## from print.formula...
  invisible(.x)
}



