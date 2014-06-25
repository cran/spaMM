muetafn <-
function(family,eta,BinomialDen,priorWeights=1) { ## note outer var BinomialDen 
  ## a patch for large eta in poisson case
  if (family$family=="poisson" && family$link =="log") {
    eta <- pmin(30,eta) ## 100 -> mu = 2.688117e+43 ; 30 -> 1.068647e+13
  }
  mu <- family$linkinv(eta) ## linkinv(eta) is FREQS for binomial, COUNTS for poisson...
  if (family$link %in% c("logit","probit","cloglog")) {
        mu[mu > (1-1e-12)] <- (1-1e-12)
        mu[mu < (1e-12)] <- (1e-12)
  }
  dmudeta <- family$mu.eta(eta) ## aberrant at hoc code for cloglog 'elsewhere'...
  Vmu <- family$variance(mu)
  if (family$family=="binomial") {
      Vmu <- Vmu * BinomialDen ## not checked for probit et cloglog
      mu <- mu *BinomialDen
      dmudeta <- dmudeta * BinomialDen
  } 
  GLMweights <- priorWeights * dmudeta^2 /Vmu ## must be O(n) in binomial cases
  if (any(is.nan(GLMweights))) {
#    calls <- sys.calls()
#    ncalls <- length(calls)
#    HLfitcall <- calls[[ncalls-1]]
#    save(HLfitcall,file=generateFileName("HLfitcall")) ## FR->FR debug  code, not for package
    stop("NaN GLMweights generated in 'muetafn'")
  }
return(list(mu=mu,dmudeta=dmudeta,Vmu=Vmu,GLMweights=GLMweights))
}
