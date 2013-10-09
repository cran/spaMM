simulate.HLfit <-
function(object, nsim = 1, seed = NULL, ...) { ## object must have class HLfit; corr pars are not used, but the ZL matrix is.
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
  nr <- ncol(object$ZLMatrix)
  nobs <- nrow(object$ZLMatrix)
  ##FR->FR quickly modified to conform to the generic; but not optimized, and tested only for nsim=1
  newU <- replicate(nsim,{
    switch(tolower(object$rand.family$family),
      gaussian = rnorm(nr,sd=sqrt(object$lambda)),
      gamma = rgamma(nr,shape=1/object$lambda,scale=object$lambda),
      beta = rbeta(nr,1/(2*object$lambda),1/(2*object$lambda)),
      "inverse.gamma" = 1/rgamma(nr,shape=1+1/object$lambda,scale=object$lambda), ## yields inverse gamma (1+1/object$lambda,1/object$lambda)
      stop("(!) random sample from given rand.family not yet implemented")
    )},simplify=TRUE) ## should have nsim columns
  newV <- object$rand.family$linkfun(newU) ## each column a simulation
  ## rebuild linear predictor in three steps
  if (ncol(object$X)>0) {eta <- object$X %*% object$fixef} else {eta <- rep(0,nobs)}
  offset<-object$predictor$offset
  if (!is.null(offset)) eta <- eta + offset 
  if (nsim==1) {
    eta <- eta + object$ZLMatrix %*% newV 
  } else eta <-  matrix(rep(eta,nsim),ncol=nsim) + object$ZLMatrix %*% newV ## nobs rows, nsim col
  mu <- object$family$linkinv(eta) ## ! freqs for binomial, counts for poisson
  ntot <- length(mu) ## total # of elements for a matrix 
  ## the docs of the random functions incorrectly say that the second argument is a 'vector'
  newy <- switch(tolower(object$family$family),
    gaussian = rnorm(ntot,mean=mu,sd=sqrt(object$phi)),
    poisson = rpois(ntot,mu),
    binomial = rbinom(ntot,size=object$weights,prob=mu),
    gamma = rgamma(ntot,shape= mu^2 / object$phi, scale=object$phi/mu),
    stop("(!) random sample from given family not yet implemented")
  ) ## vector
  if (nsim>1) {
    resu <- matrix(newy,ncol=nsim)
  } else {
    resu <- newy
  }
  return(resu)    
}
