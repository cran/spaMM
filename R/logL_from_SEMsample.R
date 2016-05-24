logL_from_SEMsample <- function(betaMat, lambdaVec,  
                                SEMsample, ZAL, nobs, X.pv, off, SEMlogL,
                                pmvnorm_arglist, pMVN_arglist
                               ) {
  beta_eta <- colMeans(betaMat[SEMsample,,drop=FALSE]) 
  lambda <- mean(lambdaVec[SEMsample]) ## lambda_est will be given a different length
  ## estim lik
  fix <- X.pv %*% beta_eta + off
  ZALtZAL <- tcrossprodCpp(ZAL)
  integrationSig <- lambda*ZALtZAL+diag(nobs)  
  if (SEMlogL=="pmvnorm") { 
    pmvnorm_arglist$sigma <- integrationSig
    Lapp <- do.call("pmvnorm",pmvnorm_arglist) 
    # 
    # pmvnorm may fail in various ways. Worst,
    ## pmvnorm bug, difficult to replicate without control of seed:
    ## Lapp my be NaN with "normal completion" message.
    pmvnorm_fatal <- (is.na(Lapp) || Lapp[1]==0 || attr(Lapp,"error")==0) 
    if ( (! pmvnorm_fatal) && attr(Lapp,"msg")=="Completion with error > abseps" ) {
      ## we could increase $maxpts, but this does not appear useful (and must handle the quote object)
      # we could now use $abseps > 0 but given the $releps=0.5 failed, it's not clear what this would improve
      # => do nothing, keep the large Lapp error
    }
  } else pmvnorm_fatal <- FALSE
  if (pmvnorm_fatal || SEMlogL %in% c("halton","pMVN")) { ## input rand_seq or default GHK method
    pMVN_arglist$L <- RcppChol(integrationSig)$L ## t(R::chol)
    #if (SEMlogL=="halton") { ## but avoid the fatal qrng::ghalton() for dim > 360
    #  pMVN_arglist$rand_seq <- randtoolbox::halton(4L*nobs,d=nobs,init=FALSE,usetime=FALSE)
    #  pMVN_arglist$nrep <- NULL
    #} else 
    if (is.null(pMVN_arglist$nrep)) 
      pMVN_arglist$nrep <- quote(200L*nobs) ## cf notes 22/03/2016
    # CALL
    blob <- do.call("pMVN",pMVN_arglist) ## GHK by default
    # 
    logLapp <- blob$logInt
    seInt <- blob$seInt
    attr(logLapp,"method") <- "      logL (GHK)" ## directly usable for screen output
  } else if (SEMlogL=="pmvnorm") {
    logLapp <- log(Lapp[1L])
    seInt <- attr(Lapp,"error")/(2.575829*Lapp) ## = qnorm(0.995) because pmvnorm's error is def'ine'd as bound of 99% interval
    attr(logLapp,"method") <- "  logL (pmvnorm)" ## directly usable for screen output
  } else stop(paste("Unknown integration method",SEMlogL,"requested"))
  ## naive MC removed 04/2016, GHK is now reliable
  return(list(beta_eta=beta_eta,lambda=lambda,logLapp=logLapp,seInt=seInt))
}