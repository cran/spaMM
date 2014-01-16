print.predictor <-
function (x, showAttr = FALSE, ...) 
{
  .x <- x ## .x is original, returned invisibly
  if (! showAttr) attributes(x) <- NULL
  print.default(unclass(x), ...) ## from print.formula...
  invisible(.x)
}
