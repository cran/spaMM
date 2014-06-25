newetaFix <-
function(object,newX) {
  ## [-2] so that HLframes does not try to find the response variables  
  newMeanFrames <- HLframes(formula=attr(object$predictor,"oriFormula")[-2],data=newX) ## may need to reconstruct offset using formula term
  X.pv <- newMeanFrames$X  
  if (ncol(X.pv)>0) {etaFix <- X.pv %*% object$fixef} else {etaFix <- rep(0,nrow(newX))} 
  ## newX -> offset must be recomputed
  off <- model.offset(newMeanFrames$mf) ## Predictor has checked that there was either an offset term XOR an $offset 
  if ( is.null(off) ) { ## no offset oriFormula term
    ## then we check that no non zero $offset was used. This would make prediction generally incorrect
    off <- attr(object$predictor,"offset") ## a PROCESSED predictor always has a non-NULL offset term 
    if (any(range(off)!=c(0,0))) { ## that means there was a non trivial offset in the original call 
      message("Prediction in new design points with an offset from original design points is suspect.")
    }
  } else etaFix <- etaFix + off ## we add a non-trivial offset from the offset formula   
  list(etaFix=etaFix,newMeanFrames=newMeanFrames)
}
