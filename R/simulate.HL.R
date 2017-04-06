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
  etaFix <-  drop(newMeanFrames$X[,validnames,drop=FALSE] %*% object$fixef[validnames]) ## valide even if ncol(newMeanFrames$X) = 0
  if ( ! is.null(off)) etaFix <- etaFix + off   
  return(etaFix)
}



# simulate.HLfit(fullm[[2]],newdata=fullm[[1]]$data,size=fullm[[1]]$data$total) for multinomial avec binomial nichées de dimension différentes
# FR->FR misses the computation of random effects for new spatial positions: cf comments in the code below
simulate.HLfit <- function(object, nsim = 1, seed = NULL, newdata=NULL,
                           conditional=FALSE, verbose=TRUE,
                           sizes=object$weights, ...) { ## object must have class HLfit; corr pars are not used, but the ZAL matrix is.
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
  if (conditional) {
    if (verbose) cat("Conditional simulation:\n")
    if (is.null(newdata)) {
      mu <- predict(object)[,1L]
      if (nsim>1) mu <- matrix(rep(mu,nsim),ncol=nsim)
    } else {
      stop("'newdata' not yet implemented in this case.")
      # debut d'implementation Mais il peut y avoir des choses inteessantes dans le code cas non conditionnel
      mu_fixed <- predict(object, newdata=newdata, re.form=NA)
      eta <- object$family$linkfun(mu_fixed[,1L]) ## onlythe fixed part at this point
      ranefs <- attr(object$ZAlist,"ranefs")
      for (rit in seq_len(ranefs)) {
        loc_u_h <- 666 ## for given ranef
        n_u_h <- 666 ## for given ranef
        if (666) { ## test for Gaussian correlated effects
          randcond_bMean <- 666 # Cno Coo^{-1/2} loc_u_h a coder en vecteurs seulement...
          randcond_brand <- 666 # Cnn - Cno Coo^{-1} Con 
          # randb <- rmvnorm(nsim,mean=randcond_bMean,sigma=randcond_brand)
        } else {
          randb <- 666 # hum. Mais il peut y avoir des choses inteessantes dans le code cas non conditionnel
        }
        eta <- eta + object$ZAlist[[ranefs]] %*% randb ## could  be donne once for all ranefs 
      }
      # autre concept: simulation conditionnelle aux données(pas equiv a simul sous le modèe fitté)
      # # il faudrait avoir une decomp de la predVar en ses composantes pour chaque ranef... faisable avec re.form ? 
      # # pas clair, le subrange de beta_w_col ne s'additionnent pas...
      # if (all(attr(object$rand.families,"lcrandfamfam")=="gaussian")) {
      #   mu_with_cov <- predict(object,newdata=newdata,variances=list(linPred=TRUE,cov=TRUE))
      #   Eeta <- object$family$linkfun(mu_with_cov[,1L])
      #   eta <- rmvnorm(nsim,mean=Eeta,sigma=attr(mu_with_cov,"predVar"))
      #   mu <- object$family$linkinv(eta) ## ! freqs for binomial, counts for poisson
    }
  } else {
    if (verbose) cat("Unconditional simulation:\n")
    mu_fixed <- predict(object, newdata=newdata, re.form=NA)
    eta_fixed <- object$family$linkfun(mu_fixed[,1L])
    if (is.null(newdata)) {
    #   # rebuild linear predictor in three steps eta= X . beta + off + ZAL . RANDOM v
    #   # hence do not use object$eta which contains PREDICTED V
      nobs <- length(object$y)
    #   if (ncol(object$`X.pv`)>0) {eta <- drop(object$`X.pv` %*% object$fixef) } else {eta <- rep(0,nobs)}
    #   ##
    #   eta <- eta + attr(object$predictor,"offsetObj")$total ## a PROCESSED predictor or resid.predictor always has a non-NULL offset term  
    } else {
      nobs <- nrow(newdata)
    #   ## [-2] so that HLframes does not try to find the response variables  
    #   allFrames <- HLframes(formula=attr(object$predictor,"oriFormula")[-2],data=newdata) ## may need to reconstruct offset using formula term
    #   eta <- newetaFix(object,allFrames) ## X . beta + off
    }
    ##
    if (object$models[["eta"]] != "etaGLM") {
      if (is.null(newdata)) {
        ZAL <- get_ZALMatrix(object,as_matrix=.eval_as_mat_arg(object))
        cum_n_u_h <- attr(object$lambda,"cum_n_u_h")
        vec_n_u_h <- diff(cum_n_u_h)
      } else {
        #   ## [-2] so that HLframes does not try to find the response variables  
        allFrames <- HLframes(formula=attr(object$predictor,"oriFormula")[-2],data=newdata) 
        FL <- .spMMFactorList(object$predictor, allFrames$mf, 0L, drop=TRUE) 
        ##### the following code with NULL LMatrix ignores spatial effects with newdata:
        ZALlist <- .compute_ZAXlist(XMatrix=NULL,ZAlist=FL$Design)
        nrand <- length(ZALlist)
        # to overcome this we needto calculate the unconditional covmat including for the (nex) positions
        vec_n_u_h <- unlist(lapply(ZALlist,ncol)) ## nb cols each design matrix = nb realizations each ranef
        cum_n_u_h <- cumsum(c(0,vec_n_u_h))
        ZALlist <- lapply(seq_len(length(ZALlist)),as.matrix)
        ZAL <- do.call(cbind,ZALlist)
      }
      lcrandfamfam <- attr(object$rand.families,"lcrandfamfam") ## unlist(lapply(object$rand.families,function(rf) {tolower(rf$family)})) 
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
      } else eta <-  matrix(rep(eta_fixed,nsim),ncol=nsim) + as.matrix(ZAL %id*% newV) ## nobs rows, nsim col
    }
    mu <- object$family$linkinv(eta) ## ! freqs for binomial, counts for poisson
  }
  ## 
  phiW <- object$phi/object$prior.weights ## cf syntax and meaning in Gamma()$simulate / spaMM_Gamma$simulate
  famfam <- tolower(object$family$family)
#  if (famfam=="binomial" && is.null(size)) size <- object$weights ## original sample size by default
  respv <- function(mu) {switch(famfam,
                                gaussian = rnorm(nobs,mean=mu,sd=sqrt(phiW)),
                                poisson = rpois(nobs,mu),
                                binomial = rbinom(nobs,size=sizes,prob=mu),
                                gamma = rgamma(nobs,shape= mu^2 / phiW, scale=phiW/mu), ## ie shape increase with prior weights, consistent with Gamma()$simulate / spaMM_Gamma()$simulate
                                stop("(!) random sample from given family not yet implemented")
  )} ## vector
  if (nsim>1) { ## then mu has several columns
    resu <- apply(mu,2,respv) ## matrix
  } else {
    resu <- respv(mu)
  }
  return(resu)    
}


simulate.HLfitlist <- function(object,nsim=1,seed=NULL,newdata=object[[1]]$data,sizes=object[[1]]$weights,...) {
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
    if (length(sizes) != nrow(newdata)) {
      mess <- pastefrom("length(sizes) != nrow(newdata).",prefix="(!) From ")
      stop(mess)
    }
    for (it in seq(ncol(resu))) {
      ## it = 1 se ramène à simulate(object[[1]])
      resu[,it] <- simulate(object[[it]],newdata=newdata,sizes=sizes - cumul,verbose=FALSE)
      cumul <- rowSums(resu)  
    }
    resu <- cbind(resu,sizes - cumul) ## now 3 cols if 3 types
    rownames(resu) <- allrownames
    colnames(resu) <- attr(object,"sortedTypes")
    as.data.frame(resu)
  },simplify=FALSE)
}
