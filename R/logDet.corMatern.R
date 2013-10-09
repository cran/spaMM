logDet.corMatern <-
function(object, covariate = getCovariate(object), ...)
{
  if (!is.null(aux <- attr(object, "logDet"))) {
    return(aux)
  }
  if (is.null(aux <- attr(object, "factor"))) {
    ## getting the transpose sqrt factor
    aux <- corMatrix(object, covariate = covariate, corr = FALSE)
  }
  if (is.null(aux1 <- attr(aux, "logDet"))) {
    ## checking for logDet attribute; if not present, get corr matrix
    aux <- corMatrix(object, covariate)
    if (data.class(aux) == "list") {    # by group
      sum(log(abs(unlist(lapply(aux, function(el) svd(el)$d)))))/2
    } else {
      sum(log(abs(svd(aux)$d)))/2
    }
  } else {
    -aux1
  }
}
