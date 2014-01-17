Predictor <-
function (formula, offset=NULL, LMatrix = NULL,  AMatrix = NULL, ZALMatrix = NULL) {
  if ("predictor" %in% class(formula)) {
    pastefrom("Do not call 'Predictor' on a predictor object.",prefix="(!) From ")
  }
  oriFormula <- formula
  if (substr((as.character(formula[2])),1,5)=="cbind") { 
       positives <- strsplit(strsplit(as.character(formula[2]),"\\(")[[1]][2],",")[[1]][1] ## names of variables
       negatives <- strsplit(strsplit(as.character(formula[2]),",")[[1]][2],"\\)")[[1]][1] ## names of variables
       if (length(positives)>1 && length(negatives)>1) {
          stop("For binomial data, please use cbind(<pos>,<neg>) where at least one of <pos> and <neg> is a variable from the data frame")
       } ## because problem bootstrap otherwise...
       ## cbind syntax for binomial model => is converted to alternative syntax; cbind would fail in HLframes as marked there FR->FR (change HLframes some day ?)
       formula <- as.formula(paste(positives,"~",as.character(formula[3])))
       BinDenForm <- paste(positives,"+",negatives)       
  } else {BinDenForm <-NULL}
  if ( ! ( is.null(AMatrix) || is.list(AMatrix)) ) {
    ## assumes any of a number of matrix classes. Not clear how to test compactly for the soup of Matrix classes
    AMatrix <- list(AMatrix) ## further code espects a list (should ultimately be the same for all matrices...) 
  }
  if ( ! is.null(offset) ) {
    offterm <- findOffset(formula)
    if ( ! is.null(offterm) ) stop("in 'Predictor', offset should be given EITHER as $formula term OR as $offset element")
  } 
  #  attr(formula,"offset") <- offset
  res <- formula
  attr(res,"oriFormula") <- oriFormula
  attr(res,"ZALMatrix") <- ZALMatrix
  attr(res,"AMatrix") <- AMatrix
  attr(res,"BinDenForm") <- BinDenForm
  attr(res,"LMatrix") <- LMatrix
  attr(res,"offset") <- offset
  class(res) <- c("predictor",class(res))
  return(res)
}
