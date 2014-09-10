simulate.HLfit <-
function(object, nsim = 1, seed = NULL, newX=NULL, sizes=object$weights,...) { ## object must have class HLfit; corr pars are not used, but the ZAL matrix is.
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
  if (inherits(object,"HLfitlist")) { ## testing for list is not valid since an HLfit object is always a list
    message("simulate does not yet work on list of fits as returned by multinomial fit:")
    message(" run simulate on each of the individual fit in the list")
    stop() ## FR->FR also some basic changes in fixedLRT but more would be needed 
  }  
  if (is.null(newX)) {
    # rebuild linear predictor in three steps eta= X . beta + off + ZAL . RANDOM v
    # hence do not use object$eta which contains PREDICTED V
    nobs <- length(object$y)
    if (ncol(object$X)>0) {eta <- object$X %*% object$fixef} else {eta <- rep(0,nobs)}
    ##
    eta <- eta + attr(object$predictor,"offset") ## $offset not NULL in a processed predictor  
  } else {
    nobs <- nrow(newX)
    ## [-2] so that HLframes does not try to find the response variables  
    allFrames <- HLframes(formula=attr(object$predictor,"oriFormula")[-2],data=newX) ## may need to reconstruct offset using formula term
    eta <- newetaFix(object,newMeanFrames=allFrames)
  }
  ##
  if (any(object$models[["lambda"]] != "")) { ## i.e. not a GLM
    if (is.null(newX)) {
      ZAL <- object$ZALMatrix
      vec_n_u_h <- attr(object$lambda,"n_u_h")
    } else {
      FL <- spMMFactorList(object$predictor, allFrames$mf, 0L, 0L) 
      ZALlist <- compute.ZALlist(LMatrix=NULL,ZAlist=FL$Design,Groupings=FL$Groupings)
      nrand <- length(ZALlist)
      vec_n_u_h <- rep(0, nrand)
      for (i in 1:nrand) vec_n_u_h[i] <- ncol(ZALlist[[i]]) ## nb cols each design matrix = nb realizations each ranef
      ZAL <- do.call(cbind,ZALlist)
    }
    lcrandfamfam <- unlist(lapply(object$rand.families,function(rf) {tolower(rf$family)})) 
    cum_n_u_h <- cumsum(c(0,vec_n_u_h))
    fittedLambda <- object$fittedLambda
    newV <- lapply(seq(length(vec_n_u_h)), function(it) {
      nr <- vec_n_u_h[it]
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      loclambda <- fittedLambda[u.range]
      newU <- replicate(nsim,{
        switch(lcrandfamfam[it], ## remainder of code should be OK for rand.families
               gaussian = rnorm(nr,sd=sqrt(loclambda)),
               gamma = rgamma(nr,shape=1/loclambda,scale=loclambda),
               beta = rbeta(nr,1/(2*loclambda),1/(2*loclambda)),
               "inverse.gamma" = 1/rgamma(nr,shape=1+1/loclambda,scale=loclambda), ## yields inverse gamma (1+1/object$lambda,1/object$lambda)
               stop("(!) random sample from given rand.family not yet implemented")
        )},simplify=TRUE) ## should have nsim columns
      object$rand.families[[it]]$linkfun(newU) 
    }) ## one multi-rand.family simulation
    newV <- do.call(rbind,newV) ## each column a simulation
    if (nsim==1) {
      eta <- eta + ZAL %*% newV 
    } else eta <-  matrix(rep(eta,nsim),ncol=nsim) + ZAL %*% newV ## nobs rows, nsim col
  }
  mu <- object$family$linkinv(eta) ## ! freqs for binomial, counts for poisson
  ## 
  phiW <- object$phi/object$prior.weights ## cf syntax and meaning in Gamma()$simulate / GammaForDispGammaGLM()$simfun
  famfam <- tolower(object$family$family)
#  if (famfam=="binomial" && is.null(size)) size <- object$weights ## original sample size by default
  respv <- function(mu) {switch(famfam,
                                gaussian = rnorm(nobs,mean=mu,sd=sqrt(phiW)),
                                poisson = rpois(nobs,mu),
                                binomial = rbinom(nobs,size=sizes,prob=mu),
                                gamma = rgamma(nobs,shape= mu^2 / phiW, scale=phiW/mu), ## ie shape increase with prior weights, consistent with Gamma()$simulate / GammaForDispGammaGLM()$simfun
                                stop("(!) random sample from given family not yet implemented")
  )} ## vector
  if (nsim>1) {
    resu <- apply(mu,2,respv) ## matrix
  } else {
    resu <- respv(mu)
  }
  return(resu)    
}
