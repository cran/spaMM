CondNormfn <-
function(LMatrix,lambda) {
  type <- attr(LMatrix,"type")
  if (is.null(type)) {stop("'type' attribute needed for 'LMatrix' argument in 'CondNormfn'.")}
  if (type=="chol") {stop("LMatrix argument in 'CondNormfn' is inadequate because its 'type' attribute is 'chol'.")}
  decomp <- attr(LMatrix,type)
  diago <- decomp$d/(decomp$d+1/lambda)
  sqrtCondCovLv <- sweep(decomp$u,2,sqrt(diago),`*`); ## so that cond Corr = this.t(this)
  condLvReg <- tcrossprodCpp(sqrtCondCovLv) ## conditional regr = cond Corr
  return(list(sqrtCondCovLv=sqrtCondCovLv,condLvReg=condLvReg))
}
