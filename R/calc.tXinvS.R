calc.tXinvS <-
function(Sig,X.pv,stop.on.error,lambda_est,ranFix) { ## slow... 
  #  qr.Sig <- qr(Sig) ## FR->FR SLOW
  #  XinvS <- safesolve.qr.matrix(qr.Sig,X.pv,stop.on.error=stop.on.error) ## invSig %*% X.pv
  qr.Sig <- Rcpp_QR(Sig) ## FR->FR still SLOW
  XinvS <- solveWrap.matrix(qr.Sig,X.pv,stop.on.error=stop.on.error) ## invSig %*% X.pv
  if (class(XinvS)=="try-error") {
    if (stop.on.error) {
      mess <- pastefrom("the augmented 'Sigma' matrix appears singular. Extreme lambda/phi value and/or extremely correlated random effects?",prefix="(!) From ")
      message(mess)
      cat(paste("max(lambda estimates)=",max(lambda_est)))
      if (length(ranFix)>0) {
        cat("; correlation parameters=")
        cat(paste(names(ranFix),"=",ranFix))
      }
      largeLambdaMessages()
      #      if (is.null(options()$error)) { ## default if not error=recover or error=stop
      #        return(list(error="Singular augmented 'Sigma' matrix")) ## HLCor has code to handle return(list(error=...))
      #      } else stop() ## will call options()$error i.e. (ideally) recover
      stop()
    } else {
      return(XinvS) ## returns a try-error
    }
  } else tXinvS <- t(XinvS) ## we have to transpose either this one or X.pv, which are of the same size 
  return(tXinvS)
}
