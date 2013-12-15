calc.XinvS <-
function(Sig,X.pv,stop.on.error) {
          qr.Sig <- qr(Sig)
          tXinvS <- safesolve.qr.matrix(qr.Sig,X.pv,stop.on.error=stop.on.error) ## invSig %*% X.pv
          if (class(tXinvS)=="try-error") {
            mess <- pastefrom("the augmented 'Sigma' matrix appears singular. Extreme lambda/phi value and/or extremely correlated random effects?",prefix="(!) From ")
            message(mess)
            if (is.null(options()$error)) { ## default
              return(list(error="Singular augmented 'Sigma' matrix")) ## HLCor has code to handle return(list(error=...))
            } else stop() ## will call options()$error i.e. (ideally) recover
          } else XinvS <- t(tXinvS) ## we have to transpose either this one or X.pv, which are of the same size 
          return(XinvS)
}
