simulate.HLfit <-
function(object, nsim = 1, seed = NULL, ...) { ## object must have class HLfit; corr pars are not used, but the ZAL matrix is.
  ## RNG stuff copied from simulate.lm
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else { ## this makes changes to RNG local where 'seed' is used:
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  nr <- ncol(object$ZALMatrix)
  nobs <- nrow(object$ZALMatrix)
  ## rebuild linear predictor in three steps
  if (ncol(object$X)>0) {eta <- object$X %*% object$fixef} else {eta <- rep(0,nobs)}
  ##
  eta <- eta + attr(object$predictor,"offset") ## $offset not NULL in a processed predictor
  ##
  if (object$models[2] != "") { ## i.e. not a GLM 
    newU <- replicate(nsim,{
      switch(tolower(object$rand.family$family),
             gaussian = rnorm(nr,sd=sqrt(object$lambda)),
             gamma = rgamma(nr,shape=1/object$lambda,scale=object$lambda),
             beta = rbeta(nr,1/(2*object$lambda),1/(2*object$lambda)),
             "inverse.gamma" = 1/rgamma(nr,shape=1+1/object$lambda,scale=object$lambda), ## yields inverse gamma (1+1/object$lambda,1/object$lambda)
             stop("(!) random sample from given rand.family not yet implemented")
      )},simplify=TRUE) ## should have nsim columns
    newV <- object$rand.family$linkfun(newU) ## each column a simulation
    if (nsim==1) {
      eta <- eta + object$ZALMatrix %*% newV 
    } else eta <-  matrix(rep(eta,nsim),ncol=nsim) + object$ZALMatrix %*% newV ## nobs rows, nsim col
  }
  mu <- object$family$linkinv(eta) ## ! freqs for binomial, counts for poisson
  ## 
  respv <- function(mu) {switch(tolower(object$family$family),
    gaussian = rnorm(nobs,mean=mu,sd=sqrt(object$phi)),
    poisson = rpois(nobs,mu),
    binomial = rbinom(nobs,size=object$weights,prob=mu),
    gamma = rgamma(nobs,shape= mu^2 / object$phi, scale=object$phi/mu),
    stop("(!) random sample from given family not yet implemented")
  )} ## vector
  if (nsim>1) {
    resu <- apply(mu,2,respv) ## matrix
  } else {
    resu <- respv(mu)
  }
  return(resu)    
}
