ranef.HLfit <-
function(object,...) {
  object <- getHLfit(object)
  object$ranef #random effects \eqn{u}
}
