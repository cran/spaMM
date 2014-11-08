## from the devel version of lme4:
##' test for no-random-effect specification: TRUE if NA or ~0, other
##' possibilities are NULL or a non-trivial formula
noReForm <- function(re.form) {
  (!is.null(re.form) && !is(re.form,"formula") && is.na(re.form)) ||
    (is(re.form,"formula") && length(re.form)==2 && identical(re.form[[2]],0))
}

newetaFix <- function(object,newMeanFrames) {
  X.pv <- newMeanFrames$X  
  if (ncol(X.pv)>0) {
    etaFix <- X.pv %*% object$fixef
  } else {
    etaFix <- rep(0,nrow(newMeanFrames$mf)) ## nrow(X.pv)=0
  } 
  ## newX -> offset must be recomputed. 
  off <- model.offset(newMeanFrames$mf) ### look for offset from (ori)Formula 
  if ( is.null(off) ) { ## ## no offset (ori)Formula term. Check attribute (Predictor ensures offset is not both in formula and attr)
    ## then we check that no non zero $offset was used. This would make prediction generally incorrect
    off <- attr(object$predictor,"offsetObj")$vector ## a PROCESSED predictor or resid.predictor always has a non-NULL offset term 
    if (any(range(off)!=c(0,0))) { ## that means there was a non trivial offset in the original call 
      message("Prediction in new design points with an offset from original design points is suspect.")
    }
  } else etaFix <- etaFix + off ## we add a non-trivial offset from the offset formula   
  return(etaFix)
}



# simulate.HLfit(fullm[[2]],newX=fullm[[1]]$data,size=fullm[[1]]$data$total) for multinomial avec binomial nichées de dimension différentes
# FR->FR misses the computation of randoem effects for new spatial positions: cf comments in the code below
simulate.HLfit <- function(object, nsim = 1, seed = NULL, newX=NULL, sizes=object$weights,...) { ## object must have class HLfit; corr pars are not used, but the ZAL matrix is.
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
    if (ncol(object$`X.pv`)>0) {eta <- object$`X.pv` %*% object$fixef} else {eta <- rep(0,nobs)}
    ##
    eta <- eta + attr(object$predictor,"offsetObj")$vector ## a PROCESSED predictor or resid.predictor always has a non-NULL offset term  
  } else {
    nobs <- nrow(newX)
    ## [-2] so that HLframes does not try to find the response variables  
    allFrames <- HLframes(formula=attr(object$predictor,"oriFormula")[-2],data=newX) ## may need to reconstruct offset using formula term
    eta <- newetaFix(object,newMeanFrames=allFrames) ## X . beta + off
  }
  ##
  if (any(object$models[["lambda"]] != "")) { ## i.e. not a GLM
    if (is.null(newX)) {
      ZAL <- object$ZALMatrix
      vec_n_u_h <- attr(object$lambda,"n_u_h")
    } else {
      FL <- spMMFactorList(object$predictor, allFrames$mf, 0L, drop=TRUE) 
      ##### simulation given the observed response ! :
      #       spatial.terms <- findSpatial(locform)
      #       spatial.model <- spatial.terms[[1]] 
      #       if( is.null(spatial.model)) {
      #         uuCnewold <- NULL
      #       } else {
      #         v_h_coeffs <- predictionCoeffs(object) ## changes the coefficients in the right u_range
      #         blob <- calcNewCorrs(object=object,locdata=newX,predVar=FALSE,spatial.model=spatial.model)
      #         uuCnewold <- blob$uuCnewold
      #       }
      #       ZALlist <- compute.ZALlist(CMatrix=uuCnewold,ZAlist=FL$Design,Groupings=FL$Groupings)
      #       (unfinished:) il faut rajouter la conditional variance comme dans le SEM => simuler comme dans le SEM ?
      ##### independent simulation ! : 
      # en fait il faut recycler du code de HLCor...
      ##### the following code with NULL LMatrix ignores spatial effects with newX:
      ZALlist <- compute.ZALlist(LMatrix=NULL,ZAlist=FL$Design,Groupings=FL$Groupings)
      nrand <- length(ZALlist)
      vec_n_u_h <- rep(0, nrand)
      for (i in 1:nrand) vec_n_u_h[i] <- ncol(ZALlist[[i]]) ## nb cols each design matrix = nb realizations each ranef
      ZAL <- do.call(cbind,ZALlist)
    }
    lcrandfamfam <- attr(object$rand.families,"lcrandfamfam") ## unlist(lapply(object$rand.families,function(rf) {tolower(rf$family)})) 
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


simulate.HLfitlist <- function(object,nsim=1,seed=NULL,newX=object[[1]]$data,sizes=object[[1]]$weights,...) {
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
  replicate(nsim, {
    allrownames <- unique(unlist(lapply(object,function(hl){rownames(hl$data)})))
    resu <- matrix(0,nrow=length(allrownames),ncol=length(object)) ## two cols if 3 types
    cumul <- 0
    if (length(sizes) != nrow(newX)) {
      mess <- pastefrom("length(sizes) != nrow(newX).",prefix="(!) From ")
      stop(mess)
    }
    for (it in seq(ncol(resu))) {
      ## it = 1 se ramène à simulate(object[[1]])
      resu[,it] <- simulate(object[[it]],newX=newX,sizes=sizes - cumul)
      cumul <- rowSums(resu)  
    }
    resu <- cbind(resu,sizes - cumul) ## now 3 cols if 3 types
    rownames(resu) <- allrownames
    colnames(resu) <- attr(object,"sortedTypes")
    as.data.frame(resu)
  },simplify=FALSE)
}
  

## there is update.HL for new fits of the same X and same Y...

## there is ?? for new X and new Y... GCV vs Krig....






