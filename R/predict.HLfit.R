predict.HLfit <-
function(object,newX = NULL,...) { ## but not new Y needed; coeffs not meaningful with newX and not useful otherwise
  if (is.null(newX)) { ## trivial but convenient option
    return(object$fv)
  } ## ELSE
  ## reconstructs linear predictor
  newMeanFrames <- HLframes(formula=attr(object$predictor,"oriFormula"),data=newX) ## may need to reconstruct offset using formula term
  X.pv <- newMeanFrames$X  
  if (ncol(X.pv)>0) {eta <- X.pv %*% object$fixef} else {eta <- rep(0,nrow(newX))} 
  eta <- eta + attr(object$predictor,"offset") ## $offset not NULL in a processed predictor
  fv <- object$family$linkinv(eta) ## ! freqs for binomial, counts for poisson
  return(fv)
}
