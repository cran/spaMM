newetaFix <-
function(object,newMeanFrames) {
  X.pv <- newMeanFrames$X  
  if (ncol(X.pv)>0) {
    etaFix <- X.pv %*% object$fixef
  } else {
    etaFix <- rep(0,nrow(newMeanFrames$mf)) ## nrow(X.pv)=0
  } 
  ## newX -> offset must be recomputed
  off <- model.offset(newMeanFrames$mf) ## Predictor has checked that there was either an offset term XOR an $offset 
  if ( is.null(off) ) { ## no offset oriFormula term
    ## then we check that no non zero $offset was used. This would make prediction generally incorrect
    off <- attr(object$predictor,"offset") ## a PROCESSED predictor always has a non-NULL offset term 
    if (any(range(off)!=c(0,0))) { ## that means there was a non trivial offset in the original call 
      message("Prediction in new design points with an offset from original design points is suspect.")
    }
  } else etaFix <- etaFix + off ## we add a non-trivial offset from the offset formula   
  return(etaFix)
}
