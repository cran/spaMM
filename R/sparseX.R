.spaMM_lm.wfit <- function(x, y, offset=NULL,w=NULL) {
  if (!is.null(w)) {
    XtWX <- .ZtWZwrapper(x,w, allow_as_mat=FALSE) # without allow_as_mat, result might be dsy (eg if I forced sparse_X on an inherently dense X)
    rhs <- crossprod(x,w*y)
  } else {
    XtWX <- crossprod(x)
    rhs <- crossprod(x,y)
  }
  chmfactor <- Cholesky(XtWX)
  if (!is.null(offset)) y <- y-offset
  beta <- solve(chmfactor,rhs,system="A")
  fitted <- x %*% beta
  residuals <- y-fitted ## offset removed in each term
  if (!is.null(offset)) fitted <- fitted+offset
  return(list(coefficients=beta[,1], fitted.values=fitted, residuals=residuals, 
              df.residual=nrow(x)-ncol(x) ##assuming rank has been 'preprocessed'
        ))
}
