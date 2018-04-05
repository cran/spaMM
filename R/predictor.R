
## initialiser offset a O ?? pas actuellement pcq bcp de tests qu'il est nul
Predictor <- function (formula, offset=NULL, LMatrix = NULL,  AMatrix = NULL, ZALMatrix = NULL) {
  if (inherits(formula,"predictor")) { stop("Do not call 'Predictor' on a predictor object.")}
  formula <- .stripFormula(formula)
  # we redefine the envir of the formula to be the closure of Predictor, and we will empty before exiting Predictor
  #formula <- .stripFormula(formula)  # env = emptyenv() not possible as model.frame needs a non-empty envir.
  #oriFormula <- formula
  if ( ! ( is.null(AMatrix) || is.list(AMatrix)) ) {
    ## assumes any of a number of matrix classes. Not clear how to test compactly for the soup of Matrix classes
    AMatrix <- list(AMatrix) ## further code expects a list (should ultimately be the same for all matrices...) 
  }
  if ( ! is.null(offset) ) {
    offterm <- .findOffset(formula)
    if ( ! is.null(offterm) ) stop("in 'Predictor', offset should be given EITHER as $formula term OR as $offset element")
  } 
  if ( ! is.null(LMatrix)) {
    if ( ! is.list(LMatrix)) LMatrix <- list(dummyid=LMatrix)
    LMatrix <- lapply(LMatrix, function(lmatrix) {
      ranefs <- attr(lmatrix,"ranefs")
      if (is.null(ranefs)) {
        ranefs <- .parseBars(formula) ## FR->FR or oriformula ??? ## currently, expand does not seem important here
      } else {
        #ranefs <- unlist(lapply(ranefs, function(term) {.parseBars(as.formula(paste("bla~",term)))}))      
        ranefs <- paste("bla~",ranefs)
        ranefs <- sapply(as.formula, ranefs)
        ranefs <- unlist(lapply(ranefs, .parseBars))      
      }
      attr(lmatrix,"ranefs") <- ranefs
      #attr(lmatrix,"set_by_Predictor") <- TRUE ## else remains NULL...
      lmatrix
    })
  }
  res <- formula
  attr(res,"oriFormula") <- formula     ## they may diverge by .noOffset, so they need to be kept distinct
  attr(res,"AMatrix") <- AMatrix
  attr(res,"LMatrix") <- LMatrix ## only usage now is to be copied in AUGI0_ZX$envir$LMatrices, and immediately removed
  attr(res,"offsetObj") <- list(offsetArg=offset,nonZeroInfo= !is.null(offset))
  class(res) <- c("predictor",class(res))
  #lsv <- c("lsv",ls())
  #rm(list=setdiff(lsv,"res")) ## empties the environment pointed by the formula in "res"
  return(res)
}

## modified from getAnywhere(print.formula)
print.predictor <- function (x, showAttr = FALSE, ...) 
{
  .x <- x ## .x is original, returned invisibly
  if (! showAttr) attributes(x) <- NULL
  print.default(unclass(x), ...) ## from print.formula...
  invisible(.x)
}



