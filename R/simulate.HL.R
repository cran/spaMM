## derived from the devel version of lme4:
## tests if argument implies a model without random effects or not
## TRUE if argument is NA or ~0, 
## FALSE if argument is NULL or a non-trivial formula
.noRanef <- function(re.form) {
  (!is.null(re.form) && !is(re.form,"formula") && is.na(re.form)) ||
    (inherits(re.form,"formula") && length(re.form)==2 && identical(re.form[[2]],0))
}

.newetaFix <- function(object, newMeanFrames,validnames=NULL) {
  ## newdata -> offset must be recomputed. 
  off <- model.offset( newMeanFrames$mf) ### look for offset from (ori)Formula 
  if ( is.null(off) ) { ## ## no offset (ori)Formula term. 
    ## then we check that no non zero $offset was used. This would make prediction generally incorrect
    if (! is.null(off <- attr(object$predictor,"offsetObj")$offsetArg)) { ## that means there was a non trivial offset argument in the original Predictor(formula...)  
      message("Prediction in new design points from a fit with formula=Predictor(... non-NULL offset ...) is suspect.")
    } ## but we still proceed with this dubious offset
  }    
  ## dans l'état actuel $fixef et complet, incluant les etaFix$beta: pas besoin de les séparer
  ## mais il peut contenir des NA ! à enlever
  # le newX contient a priori les cols des etaFix$beta, contrairement à object$X.pv => don't use the latter cols 
  if (is.null(validnames)) {
    est_and_fix <- names(which(!is.na(object$fixef))) ## estimated + etaFix$beta
    validnames <- intersect(colnames(newMeanFrames$X) ,est_and_fix) # newX.pv may  contain names of unestimated coeff for the sam reason asthe original MenanFrames$X...
  }
  if (length(validnames)==0L) validnames <- c() ## without this, validnames could be character(0) and [,validnames,drop=FALSE] fails.
  etaFix <-  drop(newMeanFrames$X[,validnames,drop=FALSE] %*% object$fixef[validnames]) ## valide even if ncol(newMeanFrames$X) = 0
  if ( ! is.null(off)) etaFix <- etaFix + off   
  return(etaFix)
}

.r_resid_var <- function(mu,phiW,sizes,COMP_nu,NB_shape,zero_truncated, famfam) { 
  # we cannot use family()$simulate bc it assumes a fit object as input
  switch(famfam,
                    gaussian = rnorm(length(mu),mean=mu,sd=sqrt(phiW)),
                    poisson = .rpois(length(mu),mu,zero_truncated=zero_truncated), 
                    binomial = rbinom(length(mu),size=sizes,prob=mu),
                    Gamma = rgamma(length(mu),shape= mu^2 / phiW, scale=phiW/mu), ## ie shape increase with prior weights, consistent with Gamma()$simulate / spaMM_Gamma()$simulate
                    COMPoisson = sapply(mu, function(muv) {
                      lambda <- family$linkfun(muv,log=FALSE)
                      .COMP_simulate(lambda=lambda,nu=COMP_nu)
                    }),
                    negbin = .rnbinom(length(mu),size=NB_shape, mu_str=mu, zero_truncated=zero_truncated), 
                    stop("(!) random sample from given family not yet implemented")
  )
} ## vector-valued function from vector input

# simulate.HLfit(fullm[[2]],newdata=fullm[[1]]$data,size=fullm[[1]]$data$total) for multinomial avec binomial nichées de dimension différentes
# FR->FR misses the computation of random effects for new spatial positions: cf comments in the code below
simulate.HLfit <- function(object, nsim = 1, seed = NULL, newdata=NULL,
                           type = "marginal", conditional=NULL, verbose=TRUE,
                           sizes=NULL , ...) { ## object must have class HLfit; corr pars are not used, but the ZAL matrix is.
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
  if ( ! is.null(conditional)) {
    warning("argument 'conditional' is obsolete and will be deprecated. Use 'type' instead.")
    if (conditional) {type <- "residual"} else type <- "marginal"
  }
  if (is.null(sizes)) sizes <- .get_BinomialDen(object)
  nrand <- length(object$ZAlist)
  if (nrand==0L) {
    mu <- predict(object,newdata=newdata,binding=NA)
    if (nsim>1) mu <- replicate(nsim,mu,simplify=FALSE) # mu <- matrix(rep(mu,nsim),ncol=nsim)
  } else { ## MIXED MODEL
    if (type=="(ranef|response)") { ## Booth & Hobert approximate approach
      if (verbose) cat("Simulation from conditional random effect distribution | observed response:\n") 
      point_pred_eta <- predict(object,newdata=newdata,variances=list(BH98=TRUE)) ## (with newX.pv=0)
      if (all(attr(object$rand.families,"lcrandfamfam")=="gaussian")){
        rand_eta <- mvrnorm(n=nsim,mu=point_pred_eta[,1L],Sigma=attr(point_pred_eta,"predVar"))
        if (nsim>1L) {
          rand_eta <- t(rand_eta)
        } else rand_eta <- rand_eta[1L,]
        if ( ! is.null(zero_truncated <- object$family$zero_truncated)) {
          mu <- object$family$linkinv(rand_eta,mu_truncated=zero_truncated)
        } else mu <- object$family$linkinv(rand_eta) ## ! freqs for binomial, counts for poisson: suitable for final code
      } else stop("This conditional simulation is not implemented for non-gaussian random-effects")
    } else if ( type=="residual") { ## conditional on predicted ranefs
      if (verbose) cat("Unconditional simulation given predicted random effects:\n") 
      mu <- predict(object,newdata=newdata,binding=NA)
      if (nsim>1) mu <- replicate(nsim,mu,simplify=FALSE) #matrix(rep(mu,nsim),ncol=nsim)
    } else if ( type=="marginal"){ ## unconditional MIXED MODEL
      if (verbose) cat("Unconditional simulation:\n") 
      mu_fixed <- predict(object, newdata=newdata, re.form=NA,binding=NA) ## mu_T
      if ( ! is.null(zero_truncated <- object$family$zero_truncated)) {
        eta_fixed <- object$family$linkfun(mu_fixed, mu_truncated=zero_truncated) ## back to _U to add ranefs to the linear predictor.
      } else eta_fixed <- object$family$linkfun(mu_fixed)
      if (is.null(newdata)) {
        ZAL <- get_ZALMatrix(object)
        cum_n_u_h <- attr(object$lambda,"cum_n_u_h")
        vec_n_u_h <- diff(cum_n_u_h)
      } else { ## unconditional MM with newdata
        #   ## [-2] so that HLframes does not try to find the response variables  
        allFrames <- .HLframes(formula=attr(object$predictor,"oriFormula")[-2],data=newdata) 
        ZALlist <- .spMMFactorList(object$predictor, allFrames$mf, 0L, drop=TRUE,sparse_precision=FALSE) 
        ##### the following code ignores spatial effects with newdata:
        # to overcome this we need to calculate the unconditional covmat including for the (new) positions
        vec_n_u_h <- unlist(lapply(ZALlist,ncol)) ## nb cols each design matrix = nb realizations each ranef
        cum_n_u_h <- cumsum(c(0,vec_n_u_h))
        ZALlist <- lapply(ZALlist,as.matrix)
        ZAL <- do.call(cbind,ZALlist)
      }
      lcrandfamfam <- attr(object$rand.families,"lcrandfamfam") ## unlist(lapply(object$rand.families, function(rf) {tolower(rf$family)})) 
      fittedLambda <- object$lambda.object$lambda_est
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
      if (nsim==1L) {
        eta <- eta_fixed + drop(ZAL %id*% newV) 
        if ( ! is.null(zero_truncated <- object$family$zero_truncated)) {
          mu <- object$family$linkinv(eta,mu_truncated=zero_truncated)
        } else mu <- object$family$linkinv(eta) 
      } else {
        eta <-  matrix(rep(eta_fixed,nsim),ncol=nsim) + as.matrix(ZAL %id*% newV) ## nobs rows, nsim col
        mu <- vector("list",nsim)
        if ( ! is.null(zero_truncated <- object$family$zero_truncated)) {
          for (it in seq_len(nsim)) mu[[it]] <- object$family$linkinv(eta[,it],mu_truncated=zero_truncated)
        } else for (it in seq_len(nsim)) mu[[it]] <- object$family$linkinv(eta[,it]) 
      }
    } else stop("Unknown simulate 'type' value.")
  }
  ## ! mu := freqs for binomial, counts for poisson ## vector or matrix
  phiW <- object$phi/object$prior.weights ## cf syntax and meaning in Gamma()$simulate / spaMM_Gamma$simulate
  family <- object$family
  famfam <- object$family$family
  # resu <- object$family$simulate(object,nsim=nsim) ## could be OK and more efficient for CONDITIONAL simulation
  if (famfam =="COMPoisson") {
    COMP_nu <- environment(object$family$aic)$nu
  } else if (famfam == "negbin") {
    NB_shape <- environment(object$family$aic)$shape
    zero_truncated <- identical(object$family$zero_truncated,TRUE)
  } else if (famfam == "poisson") {
    zero_truncated <- identical(object$family$zero_truncated,TRUE)
  }
  if (nsim>1) { ## matrix output, but list input for mu, as it is easier to keep attributes on mu by using a list.
    resu <- sapply(mu, .r_resid_var, phiW=phiW,sizes=sizes,COMP_nu=COMP_nu,NB_shape=NB_shape, zero_truncated=zero_truncated,
                                         famfam=famfam) 
  } else {
    resu <- .r_resid_var(mu,phiW=phiW,sizes=sizes,COMP_nu=COMP_nu,NB_shape=NB_shape, zero_truncated=zero_truncated,
                         famfam=famfam)
  }
  return(resu)    
}


simulate.HLfitlist <- function(object,nsim=1,seed=NULL,newdata=object[[1]]$data,sizes=NULL,...) {
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
    allrownames <- unique(unlist(lapply(object, function(hl){rownames(hl$data)})))
    resu <- matrix(0,nrow=length(allrownames),ncol=length(object)) ## two cols if 3 types
    cumul <- 0
    if (length(sizes) != nrow(newdata)) stop("length(sizes) != nrow(newdata).")
    for (it in seq(ncol(resu))) {
      ## it = 1 se ramène à simulate(object[[1]])
      if (is.null(sizes)) sizes <- .get_BinomialDen(object[[it]])
      resu[,it] <- simulate(object[[it]],newdata=newdata,sizes=sizes - cumul,verbose=FALSE)
      cumul <- rowSums(resu)  
    }
    resu <- cbind(resu,sizes - cumul) ## now 3 cols if 3 types
    rownames(resu) <- allrownames
    colnames(resu) <- attr(object,"sortedTypes")
    as.data.frame(resu)
  },simplify=FALSE)
}
