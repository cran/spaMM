if (FALSE) {
  .dispFn <- function(x) log(0.9+x)
  .dispInv <- function(x) exp(x)-0.9
} else {
  ## see control_nloptr.nb for alternatives
  ## generic .dispInv for devel
  # .dispInv <- function(x,xref=1e-4,xreff=1e-2,xm=1e-3) {
  #   uniroot(function(v) dispFn(v)-x,lower=1e-6,upper=1e6,tol=1e-7)$root
  # }
  ## Still derivable at the threshold xm: 
  .dispFn <- function(x,xref=5e-5,xreff=1e-1,xm=1e-3) { 
    num <- x
    x_lo_xm <- (x<=xm)
    num[which(x_lo_xm)] <- log(xref+x[which(x_lo_xm)])
    num[which( ! x_lo_xm)] <- log(xref+xm)+log((xreff+x[which( ! x_lo_xm)])/(xreff+xm))*(xreff+xm)/(xref+xm)
    #print(num)
    num
  } 
  .dispInv <- function(x,xref=5e-5,xreff=1e-1,xm=1e-3) {
    num <- x
    x_lo_xm <- (x<log(xref+xm))
    num[which(x_lo_xm)] <- exp(x[which(x_lo_xm)])-xref
    fac <- exp((x[which( ! x_lo_xm)]-log(xref+xm))*(xm+xref)/(xm+xreff))
    num[which( ! x_lo_xm)] <-  xm*fac + xreff*(fac-1)
    #print(num)
    num
  }
}
#sapply(sapply(10^(seq(-6,6)),dispFn),dispInv)
## These must be directly applicable to say lambda=c(0.1,0.1) hence need to vectorize the test on x
# but Vectorize is not quite efficient
#.dispFn <- Vectorize(.dispFn,"x")
#.dispInv <- Vectorize(.dispInv,"x")
## for rho/trRho: rho=0 exactly is meaningful in adjacency model. rhoFn and rhoInv should be used only for other models.
.rhoFn <- function(x,RHOMAX) {log(x/(RHOMAX-x))} ## rho should be constrained to <RHOMAX and trRho should diverge as rho approaches RHOMAX
.rhoInv <- function(trRho,RHOMAX) {
  RHOMAX*exp(trRho)/(1+exp(trRho))
} 
.nuFn <- function(nu,rho,NUMAX) {log(nu/(NUMAX-nu))} ## nu should be constrained to <NUMAX and trNu should diverge as nu approaches NUMAX
.nuInv <- function(trNu,trRho,NUMAX) {NUMAX*exp(trNu)/(1+exp(trNu))}
.longdepFn <- function(longdep,LDMAX) {log(longdep/(LDMAX-longdep))} 
.longdepInv <- function(trLongdep,LDMAX) {LDMAX*exp(trLongdep)/(1+exp(trLongdep))}
.kappaFn <- function(kappa,KAPPAMAX) {log(kappa/(KAPPAMAX-kappa))} 
.kappaInv <- function(trKappa,KAPPAMAX) {KAPPAMAX*exp(trKappa)/(1+exp(trKappa))}

## spherical transfo of Pinheiro and Bates 96; see further explanation
# https://math.stackexchange.com/questions/1326462/spherical-parametrization-of-a-cholesky-decomposition/1329660#1329660
# if this is used then 2 lines must be uncommented in .makeLowerUpper() !!! 
.sphFn <- function(covpars,Xi_ncol=NULL, unbounded=TRUE) { ## from *cov* parameter vector to trRancoefs
  if (is.null(Xi_ncol)) Xi_ncol <- floor(sqrt(length(covpars)*2))
  cholmat <- diag(Xi_ncol)
  cholmat[lower.tri(cholmat,diag = TRUE)] <- covpars
  lambdas <- diag(cholmat)
  tocorr <- diag(x=sqrt(1/lambdas))
  cholmat <- tocorr %*% cholmat %*% tocorr
  diag(cholmat) <- diag(cholmat)/2
  cholmat <- t(cholmat) + cholmat
  if (TRUE){ ## svd +qr more precise ? .sphFn is not often called... (only in inits outer optim)
    svdv <- svd(cholmat)
    esys <- with(svdv,list(values=d,vectors=(u+v)/2))
    crossfac <- .Dvec_times_matrix(sqrt(esys$values), t(esys$vectors))
    cholmat <- qr.R(qr(crossfac))
  } else cholmat <- .CHOL(cholmat) ### potential source of problems as crossprod() in .sphInv() may not be the exact inverse
  if (unbounded) {
    corrpars <- numeric(0)
    for (i in 2:Xi_ncol) {
      # operations on col i of the chol 
      aux <- acos(cholmat[1:(i - 1), i]/sqrt(cumsum(cholmat[i:1, i]^2)[i:2]))
      # corr: -1   0   1
      # aux:  pi pi/2  0
      if (.spaMM.data$options$rC_unbounded) {
        corrpars <- c(corrpars, log(aux/(pi - aux))) # log:  Inf  0  -Inf
      } else corrpars <- c(corrpars,2*(aux/pi)-1) # second terms in (-1,1)
    }
    trRancoef <- c(log(lambdas)/2, corrpars) ## trRancoef:= unconstrained parameterization of the cov matrix with log sigma in first positions
    if(any(is.nan(trRancoef))) stop("any(is.nan(trRancoef))")
  } else trRancoef <- c(log(lambdas)/2, cholmat[lower.tri(cholmat)])
  return(trRancoef) 
}

.sphInv <- function(trRancoef,Xi_ncol=NULL) { ## from trRancoefs, with log sigma in first positions, to *cov* parameter vector
  if (is.null(Xi_ncol)) Xi_ncol <- floor(sqrt(length(trRancoef)*2))
  covmat <- .calc_cov_from_trRancoef(trRancoef,Xi_ncol=Xi_ncol)
  return(covmat[lower.tri(covmat,diag = TRUE)]) ## vector representation of the cov matrix using covariances
}
#.sphInv(.sphFn(c(1,1,1,5,5,14)))

.regularized_eigen <- function(compactcovmat, condnum, raw_regul = NULL) {
  if ( ! is.null(raw_regul)) {
    Xi_ncol <- ncol(compactcovmat)
    diagPos <- seq.int(1L,Xi_ncol^2,Xi_ncol+1L)
    compactcovmat[diagPos] <- compactcovmat[diagPos] + raw_regul 
  }
  #if ( ! is.null(condnum)) {
    if (TRUE) {
      svdv <- svd(compactcovmat)
      esys <- with(svdv,list(values=d,vectors=(u+v)/2))  ## a bit heuristic, but this appears more accurate
    } else {
      esys <- eigen(compactcovmat,symmetric = TRUE) ## COV= .ZwZt(esys$vectors,esys$values) ##
    }
    target_min_d <- max(esys$values)/condnum ## so that corrected condition number is at most the denominator: 
    d_corr <- max(c(0,target_min_d-esys$values))
    d_regul <- esys$values+ d_corr # all diag is corrected => added a constant diagonal matrix to compactcovmat
    #
    compactcovmat <- .ZWZt(esys$vectors, d_regul) # at least for return value
    return(structure(compactcovmat, esys=esys, d_regul=d_regul))
  #} else return(compactcovmat)
}

.calc_cov_from_trRancoef <- function(trRancoef, Xi_ncol=NULL) { ## note use in .sphInv
  if (is.null(Xi_ncol)) Xi_ncol <- attr(trRancoef,"Xi_ncol")
  if (is.null(Xi_ncol)) Xi_ncol <- floor(sqrt(length(trRancoef)*2))
  cumnp <- c(0,cumsum(seq(Xi_ncol-1)))
  offdiag <- trRancoef[-seq(Xi_ncol)] #*10
  if (.spaMM.data$options$rC_unbounded) {
    ox <- pi*exp(offdiag)/(1+exp(offdiag)) # PinheiroB96
  } else ox <- pi*((offdiag + 1)/2)
  cholmat <- diag(nrow=Xi_ncol)
  for (col in 2:Xi_ncol) {
    oxrange <- (cumnp[col-1L]+1L):(cumnp[col])
    subox <- cos(ox[oxrange])
    cossign <- sign(subox)
    subox <- - diff(c(1,cumprod(1-subox^2)))
    cholmat[seq(col-1L),col] <- cossign*sqrt(subox)
    cholmat[col,col] <- sqrt(1-sum(subox))
  }
  corrmat <- crossprod(cholmat)  ## .crossprod not interesting for small matrices
  sqrtlam <- diag(x=exp(trRancoef[1:Xi_ncol])) # tried .dispInv
  covmat <- sqrtlam %*% corrmat %*% sqrtlam
  if(any(is.nan(covmat))) stop("any(is.nan(covmat))")
  return(covmat)
}

.calc_invL_from_trRancoef <- function(trRancoef, Xi_ncol=NULL,sing_threshold=.spaMM.data$options$invL_threshold) { ## note use in .sphInv
  if (is.null(Xi_ncol)) Xi_ncol <- attr(trRancoef,"Xi_ncol")
  if (is.null(Xi_ncol)) Xi_ncol <- floor(sqrt(length(trRancoef)*2))
  cumnp <- c(0,cumsum(seq(Xi_ncol-1)))
  offdiag <- trRancoef[-seq(Xi_ncol)] 
  if (.spaMM.data$options$rC_unbounded) {
    ox <- pi*exp(offdiag)/(1+exp(offdiag)) # PinheiroB96
  } else ox <- pi*((offdiag + 1)/2)
  cholmat_corr <- diag(nrow=Xi_ncol)
  for (col in 2:Xi_ncol) {
    oxrange <- (cumnp[col-1L]+1L):(cumnp[col])
    subox <- cos(ox[oxrange])
    cossign <- sign(subox)
    subox <- - diff(c(1,cumprod(1-subox^2)))
    cholmat_corr[seq(col-1L),col] <- cossign*sqrt(subox)
    cholmat_corr[col,col] <- sqrt(1-sum(subox))
  }
  inv_cholmat_corr <- solve(cholmat_corr)
  dvec <- exp(-trRancoef[1:Xi_ncol]) #1/sigma
  crossfac_precmat <- t(.Dvec_times_matrix(dvec, inv_cholmat_corr))
  #print(range(crossfac_precmat))
  #kappa(crossfac_precmat)>sing_threshold
  if (prod(.diagfast(crossfac_precmat,Xi_ncol))>sing_threshold) { ## F I X_invL any(.diagfast(crossfac_precmat,Xi_ncol)>sing_threshold) maybe not the good test
    return(NULL)
  } else return(crossfac_precmat)
}



.ranCoefsFn <- function(vec) { # from canonical vector (var+corr) space
  if ( ! is.null(vec)) {
    transf <- attr(vec,"transf")
    if ( TRUE &&  ! is.null(transf)) { 
      resu <- structure(transf,
                        canon=as.vector(vec), # as.vector() to drop attributes
                        Xi_ncol=attr(vec,"Xi_ncol"))
      return(resu)
    }
    Xi_ncol <- attr(vec,"Xi_ncol")
    if (is.null(Xi_ncol)) Xi_ncol <- floor(sqrt(length(vec)*2))
    #
    compactcovmat <- matrix(0,nrow=Xi_ncol,ncol=Xi_ncol)
    lowerbloc <- lower.tri(compactcovmat,diag=TRUE) ## a matrix of T/F !
    compactcovmat[lowerbloc] <- vec
    diagPos <- seq.int(1L,Xi_ncol^2,Xi_ncol+1L)
    lambdas <- compactcovmat[diagPos]
    # lambdas <- lambdas + 1e-10 # here that woudl affect the fonversion from correlation to covariances
    sigmas <- diag(x=sqrt(lambdas)) 
    #
    compactcovmat[diagPos] <- 1
    compactcovmat <- sigmas %*% compactcovmat %*% sigmas ## using lower tri only
    trRancoef <- structure(.sphFn(compactcovmat[lowerbloc],Xi_ncol=Xi_ncol),
                      canon=vec,
                      Xi_ncol=Xi_ncol)
    return(trRancoef) # to transformed, 'unbounded' parameter space with log sigma in first positions
  } else return(NULL)
}

.calc_cov_from_ranCoef <- function(ranCoef, Xi_ncol=attr(ranCoef,"Xi_ncol")) {
  # assume input is marginal variances + correlation as in user input, but ordered as in lower.tri
  compactcovmat <- matrix(0,nrow=Xi_ncol,ncol=Xi_ncol)
  lowerbloc <- lower.tri(compactcovmat,diag=TRUE) ## a matrix of T/F !
  compactcovmat[lowerbloc] <- ranCoef
  diagPos <- seq.int(1L,Xi_ncol^2,Xi_ncol+1L)
  lambdas <- compactcovmat[diagPos]
  # lambdas <- lambdas + 1e-10 # here that would affect the conversion from correlation to covariances
  sigmas <- diag(x=sqrt(lambdas)) 
  compactcovmat <- (compactcovmat+t(compactcovmat))
  compactcovmat[diagPos] <- 1
  compactcovmat <- sigmas %*% compactcovmat %*% sigmas # cov2cor use distinct, recommended syntax for the last operations
  return(compactcovmat)
}  

.ranCoefsInv <- function(trRancoef) { # from transformed, 'unbounded' parameter space to correlation parameter vector
  if ( ! is.null(trRancoef)) {
    canon <- attr(trRancoef,"canon")
    if (TRUE &&  ! is.null(canon)) {
      resu <- structure(canon,
                        transf=as.vector(trRancoef), # as.vector() to drop attributes
                        Xi_ncol=attr(trRancoef,"Xi_ncol"))
      return(resu)
    }
    Xi_ncol <- attr(trRancoef,"Xi_ncol")
    if (is.null(Xi_ncol)) Xi_ncol <- floor(sqrt(length(trRancoef)*2))
    covpars <- .sphInv(trRancoef, Xi_ncol=Xi_ncol) # to vector of parameters with *cov*ariances
    # from vector with *cov*ariances to canonical vector with *corr*elations:
    rancoefs <- diag(nrow=Xi_ncol)
    rancoefs[lower.tri(rancoefs,diag = TRUE)] <- covpars
    diagPos <- seq.int(1L,Xi_ncol^2,Xi_ncol+1L)
    lambdas <- rancoefs[diagPos]
    # lambdas <- lambdas + 1e-10 # here affecting also the conversion from covariances to correlations
    torancoefs <- diag(x=sqrt(1/lambdas))
    rancoefs <- torancoefs %*% rancoefs %*% torancoefs ## lower tri of corr mat
    rancoefs[diagPos] <- lambdas
    varcorr <- structure(rancoefs[lower.tri(rancoefs,diag = TRUE)],
                      transf=trRancoef,
                      Xi_ncol=Xi_ncol)
    return(varcorr) # to canonical vector (var+corr) space
  } else return(NULL)
}
#(.ranCoefsInv(.ranCoefsFn(c(1.0000000,  0.4472136,  0.2672612,  5.0000000,  0.5976143, 14.0000000))))

.NB_shapeFn <- function(x) {log(log(1+x))} ## drastic handling of flat likelihoods for high shape= LOW variance. .../...
# .../... negbin example in gentle intro is a test (using optimize() -> initial value cannot be controlled)
.NB_shapeInv <- function(x) {exp(exp(x))-1}


## FR->FR rho/sqrt(nu) scaling => should be fixed nu and transformed rho to handle vectorial rho
## thus nu fns should be indep of nu and rho fns should be functions of nu ! :
if (FALSE) {
  futurerhoFn <- function(rho,nu) {
    rhosc <- rho/sqrt(nu)
    log(rhosc/(.spaMM.data$options$RHOMAX-rhosc))
  }
  futurerhoInv <- function(trRhosc,trNu) {
    nu <- .spaMM.data$options$NUMAX*exp(trNu)/(1+exp(trNu))
    rhosc <- .spaMM.data$options$RHOMAX*exp(trRhosc)/(1+exp(trRhosc))
    rhosc*sqrt(nu)
  }
}

.calc_LowUp_trRancoef <- function(trRancoef, ## single vector 
                                  Xi_ncol,
                                  tol_ranCoefs,
                                  adjust=.spaMM.data$options$max_corr_ranCoefs
                                  ) { 
  lower <- c(
    rep(log(tol_ranCoefs["lo_lam"]), Xi_ncol), # log(sigma)
    rep(-1+tol_ranCoefs["corr"], Xi_ncol*(Xi_ncol-1L)/2L) ## better then -1+tol_ranCoefs
  )
  lower <- adjust*lower
  lower <- pmin(lower, trRancoef)
  upper <- c(
    rep(log(tol_ranCoefs["up_lam"]), Xi_ncol), # log(sigma)
    rep(1-tol_ranCoefs["corr"], Xi_ncol*(Xi_ncol-1L)/2L) ## better then -1+tol_ranCoefs
  )
  upper <- adjust*upper
  upper <- pmax(upper, trRancoef)
  return(list(lower=lower,upper=upper)) ## transformed scale with log sigma in first positions
}


.match_coords_in_tUniqueGeo <- function(coords,tUniqueGeo) {any(apply(tUniqueGeo,1L,identical,y=coords))}

.checkDistMatrix <- function(distMatrix,data,coordinates) {
  ## HLCor + non-null distMatrix (cf example ?HLCor MLdistMat) => plantage un peu bizarre dans .checkDistMatrix, 
  #   problème de transposition de matrice uniqueGeo. A priori le code actuel est correct, mais logique pas claire.
  if (inherits(distMatrix,"dist")) {
    usernames <- labels(distMatrix)
  } else if (inherits(distMatrix,"matrix")) {
    usernames <- rownames(distMatrix)
  } else message(paste("(!) 'distMatrix' is neither a 'matrix' or 'dist' object. Check the input. I exit."))
  ## chol() fails on distances matrices with repeated locations (which are pos SD)... but chol() not used by default
  ## the following code assumes that distMatrix deals only with unique locations, and checks this
  ## HENCE ******* distMatrix must refer to unique values of a grouping variable *********
  checknames <- all(sapply(usernames,`%in%`, table=rownames(data))) ## 
  if (!checknames) {
    warning("The rownames of 'distMatrix' are not rownames of the 'data'. Further checking of 'distMatrix' is not possible.")
    nbUnique <- NA
  } else {
    uniqueGeo <- .calcUniqueGeo(data=data[usernames,coordinates,drop=FALSE]) ## check that this corresponds to unique locations
    nbUnique <- nrow(uniqueGeo)
    if (nbUnique != nrow(distMatrix)) {
      stop("The dimension of 'distMatrix' does not match the number of levels of the grouping variable")
    } else { ## check order
      redondGeo <- data[,coordinates,drop=F]
      tUniqueGeo <- t(uniqueGeo)
      designRU <- apply(redondGeo, 1L, .match_coords_in_tUniqueGeo, tUniqueGeo=tUniqueGeo) ## has no names
      # designRU <- apply(redondGeo,1, function(v) {which(apply(v==t(uniqueGeo),2,all))}) ## has no names
      ## eg 1 1 2 2 3 2 3 4 is valid for 8 obs, 4 unique locations
      designRU <- unique(as.vector(designRU)) ## should then be 1 2 3 4
      ## but if distMatrix in reverse order, the first row of redondGeo would match the 4th of uniqueGeo and then the following test is FALSE:
      if ( ! all (designRU==seq_len(length(designRU))) ) {
        stop("The rows of 'distMatrix' are not ordered as rows of the 'data'.")
      }
    } 
  }
  nbUnique ## if stop() did not occur
}

.expand_GeoMatrices <- function(w_uniqueGeo, e_uniqueGeo, coords_nesting, coord_within, dist.method) {
  rownames(e_uniqueGeo) <- seq(nrow(e_uniqueGeo)) ## local & unique rownames
  ## Remarkably the next line works only if the cols are not factors ! Unless we have a fix for this,
  #  uniqueGeo classes should be integer not factor: see instances of as.numeric(levels(fac)) in the sources.
  rows_bynesting <- by(e_uniqueGeo ,e_uniqueGeo[,coords_nesting],rownames) 
  distMatrix <- matrix(Inf,ncol=nrow(e_uniqueGeo),nrow =nrow(e_uniqueGeo))
  rownames(distMatrix) <- colnames(distMatrix) <- rownames(e_uniqueGeo) ## same trivial rownames
  for (lit in seq_len(length(rows_bynesting))) {
    blockrows <- rows_bynesting[[lit]]
    within_values <- e_uniqueGeo[blockrows,coord_within]
    w_dist <- proxy::dist(within_values,method=dist.method)
    distMatrix[blockrows,blockrows] <- as.matrix(w_dist)
  }
  return(list(distMatrix=distMatrix))
}

.calc_fullrho <- function(rho, coordinates,rho_mapping) {
  if ( is.null(rho_mapping) ) {
    if (length(rho)==1L ) rho <- rep(rho,length(coordinates))
    names(rho) <- coordinates
  } else {
    rho <- rho[rho_mapping]
    names(rho) <- names(rho_mapping)
  }  
  rho
}


.reformat_verbose <- function(verbose,For) {
  if (is.null(verbose)) verbose <- logical(0)
  if (is.list(verbose)) stop("is.list(verbose)")
  if (is.na(verbose["TRACE"])) verbose["TRACE"] <- FALSE ## for fitme
  if (is.na(verbose["trace"])) verbose["trace"] <- FALSE
  if (is.na(verbose["SEM"])) verbose["SEM"] <- FALSE
  if (is.na(verbose["iterateSEM"])) verbose["iterateSEM"] <- TRUE ## summary info and plots for each iteration
  if (is.na(verbose["phifit"])) verbose["phifit"] <- TRUE ## DHGLM information
  ## when a fit is by HLCor or HLfit, the serious message at the end of HLfit_body are by default displayed,
  ## but not in the other cases:
  if (is.na(verbose["all_objfn_calls"])) verbose["all_objfn_calls"] <- switch(For, "HLCor" = TRUE, "HLfit" = TRUE, 
                                                        "corrHLfit" = FALSE, "fitme" = FALSE)
  return(verbose)
}

.reformat_LevenbergM <- function(LevenbergM) {
  if (is.list(LevenbergM)) stop("is.list(LevenbergM)")
  ## allows user to control the default globally by spaMM.options(): BUT this should typically be NULL
  default_LM_start <- .spaMM.data$options$LevenbergM
  if ( ! is.logical(default_LM_start)) default_LM_start <- NA ## converts NULL to NA (default resolved further in preprocess)
  default_prefit <- .spaMM.data$options$PQL_prefit
  if ( ! is.logical(default_prefit)) default_prefit <- NA ## converts NULL to NA (default resolved further in preprocess)
  resu <- c(user_LM=LevenbergM["LevenbergM"][[1L]], user_prefit=LevenbergM["PQL_prefit"][[1L]], 
            # [.][[1L]] important to get rid of name yet be valid for implicit NAs
            default_LM_start=default_LM_start, default_prefit=default_prefit, force=NA) 
  if (length(LevenbergM)==1L && is.logical(LevenbergM)) resu["user_LM"] <- LevenbergM ## T, F, or NA : handles non-vector input
  return(resu) 
} 

.determine_LevenbergM <- function(LevenbergM, which="user_LM") {
  if (is.na(info <- LevenbergM["force"])) {
    if (is.na(info <- LevenbergM[which])) {
      if (which=="user_prefit") {
        info <- LevenbergM["default_prefit"]
      } else info <- LevenbergM["default_LM_start"]
    }  
  }
  if (! is.logical(info) || is.na(info)) stop("invalid input: check input (invalid 'which' or invalid user input?)")
  return(info)
} 

.calc_inits_ranCoefs <- function(init,init.optim,init.HLfit,ranFix,user.lower,user.upper) {
  # fixme c/should allow a user-provided init.HLfit$ranCoefs to control what remains in init.optim relative to init.HLfit 
  if ( ! is.null(init.optim$ranCoefs)) { ## should always be a complete ordered list in the case for random-coefficient models
    if (! is.null(init$ranCoefs)) for (st in names(init$ranCoefs)) init.optim$ranCoefs[[st]] <- init$ranCoefs[[st]]
    # fixme ideally the code whould check that st corresponds to a random-coefficient term so that the user does not mess everything....
    init$ranCoefs <- init.optim$ranCoefs
    trRanCoefs <- lapply(init.optim$ranCoefs,.ranCoefsFn)
    init.optim$trRanCoefs <- trRanCoefs
    init.optim$ranCoefs <- NULL
  }
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

.rename_lambda <- function(lambda) {
  if (! is.null(lambda)) {
    if (is.null(names(lambda))) {
      names(lambda) <- paste(seq_len(length(lambda))) 
    } else { ## names are not null, they should contain integer info
      checknames <- strsplit(names(lambda), "[^0-9]")
      for (rd in seq_along(checknames)) {
        if (length(checknames[[rd]])) {
          if (checknames[[rd]][1L]=="") {
            checknames[[rd]] <- list(NULL)
          } else checknames[[rd]] <- checknames[[rd]][1L] # this must be a char_rd
        } else checknames[[rd]] <- list(NULL)
      }
      n_hasint <- length(unames <- unlist(checknames))
      if (n_hasint==0L) {
        names(lambda) <- paste(seq_len(length(lambda))) 
      } else if (n_hasint < length(lambda)) {
        stop(paste0("lambda name(s) '",paste(names(lambda),collapse="', '"),"' have no obvious integer interpretation.\n",
                    "  Check names of initial or fixed lambda values."))
      } else names(lambda) <- unames
    }
  }
  return(lambda)
}

.calc_inits_dispPars <- function(init,init.optim,init.HLfit,ranFix,user.lower,user.upper) { 
  ## does not modify init.HLfit, but keeps its original value. Also useful to keep ranFix for simple and safe coding
  lambda <- init.optim$lambda # <- .rename_lambda(init.optim$lambda) 
  ## : matters ultimately  for optPars names in fitme's refit code; and now for merging with (canon)init from IMRF
  init$lambda[names(lambda)] <- lambda # merging with init$lambda from IMRF; works even if init$lambda was NULL
  ## the following code assumes that any inconsistency between init.optim$lambda and ranFix$trLambda has been handled by .calc_inits_hyper()
#  ranFix$lambda <- .rename_lambda(ranFix$lambda) 
  lambdaNames <- sort(unique(c(names(init.optim$lambda),names(ranFix$lambda))))
  optimVal <- ( ! is.na(init.optim$lambda[lambdaNames])) # not yet up to date: 
  # init.optim$lambda has kept default initial value except those removed by .calc_inits_IMRF()
  fixVal <- ( ! is.na(ranFix$lambda[lambdaNames]))
  conflicts <- ( optimVal & fixVal)
  if (any(conflicts))  stop("'fixed lambda' and 'init lambda' arguments conflict with each other.")
  ##
  if (length(init$lambda)) { ## length rather than !is.null() bc .dispFn(numeric(0)) returns list() bc of the Vectorize code:
    ## zut <- function(x) return(666); zut <- Vectorize(zut); zut(numeric(0)) ## gives list()
    init$lambda <- pmax(init$lambda, 1e-4) ## see remark on init$phi below. Lower user- or bootstrap- inits are not heeded.
    init.optim$trLambda <- .dispFn(init$lambda) 
  }
  init.optim$lambda <- NULL
  if (is.null(.getPar(ranFix,"phi"))) { ## FIXME getPar not useful ?
    init$phi <- init.optim$phi 
    if (!is.null(init$phi)) {
      ## pmax(), really ? FIXME
      init$phi <- pmax(init$phi, 1e-4) 
      init.optim$trPhi <- .dispFn(init$phi)
      init.optim$phi <- NULL
    }
  } else {
    if (!is.null(init.optim$phi)) stop("(!) Arguments 'ranFix$phi' and 'init$phi' conflict with each other.")  
  } 
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

.calc_inits_ARphi <- function(init,init.optim,init.HLfit,ranFix,user.lower=list(),user.upper=list(), char_rd) {
  if (is.null(ranFix$corrPars[[char_rd]][["ARphi"]])) {
    if (! is.numeric(HLfit_ARphi <- init.HLfit$corrPars[[char_rd]][["ARphi"]])) { ## if condition is TRUE then HLfit_ARphi should be NA or NULL
      init_ARphi <- init.optim$corrPars[[char_rd]][["ARphi"]] 
      if (is.null(init_ARphi)) init_ARphi <- (.get_cP_stuff(user.lower,"ARphi",char_rd)+.get_cP_stuff(user.upper,"ARphi",char_rd))/2 
      if ( ! length(init_ARphi)) init_ARphi <- 0. ## may be numeric(0) from previous line
      if (! is.null(HLfit_ARphi)) { ## then typically NA, given that the original init.HLfit$ was not numeric... 
        init.HLfit$corrPars[[char_rd]] <- list(ARphi=init_ARphi) ## numeric
        init.optim$corrPars[[char_rd]] <- NULL ## 'should' not be necessary
      } else { ## then, HLfit_ARphi was NULL and stays NULL
        init.optim$corrPars[[char_rd]] <- list(ARphi=init_ARphi) ## numeric
      } ## now HLfit_ARphi is NULL or numeric, and optim_ARphi is numeric 
      init$corrPars[[char_rd]] <- list(ARphi=init_ARphi)
    } else { ## else init.HLfit$ -> HLfit_ARphi is numeric
      init.optim$corrPars[[char_rd]] <- NULL ## 'should' not be necessary
      init$corrPars[[char_rd]] <- list(ARphi=init_ARphi)
    }
  }
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

.calc_inits_nu <- function(init,init.optim,init.HLfit,ranFix,control_dist_rd,optim.scale,NUMAX,user.lower,user.upper,char_rd) {
  if (is.null(.get_cP_stuff(ranFix,"nu",which=char_rd))) { 
    nu <- .get_cP_stuff(init.optim,"nu",which=char_rd)
    if (is.null(nu)) {
      if ( ! is.null(dm <- control_dist_rd$`dist.method`) && dm %in% c("Geodesic","Earth")) {
        nu <- 0.25  
      } else nu <- 0.5 
    }
    init$corrPars[[char_rd]] <- .modify_list(init$corrPars[[char_rd]], list(nu=nu)) ## canonical scale
    optim_cP <- init.optim$corrPars[[char_rd]]
    if (optim.scale=="transformed") {
      rho <- .get_cP_stuff(ranFix,"rho",which=char_rd)
      if (is.null(rho)) rho <- .get_cP_stuff(init,"rho",which=char_rd)
      optim_cP <- .modify_list(optim_cP, list(trNu=.nuFn(nu,rho,NUMAX=NUMAX))) 
      optim_cP$nu <- NULL
    } else optim_cP <- .modify_list(optim_cP, list(nu=nu)) 
    init.optim$corrPars[[char_rd]] <- optim_cP
  } 
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

.calc_inits_cauchy <- function(init,init.optim,init.HLfit,ranFix,control_dist_rd,optim.scale,LDMAX,user.lower,user.upper,char_rd) {
  if (is.null(.get_cP_stuff(ranFix,"shape",which=char_rd))) { 
    shape <- .get_cP_stuff(init.optim,"shape",which=char_rd)
    if (is.null(shape)) shape <- 1
    init$corrPars[[char_rd]] <- .modify_list(init$corrPars[[char_rd]], list(shape=shape)) ## canonical scale
    optim_cP <- init.optim$corrPars[[char_rd]]
    optim_cP <- .modify_list(optim_cP, list(shape=shape)) 
    init.optim$corrPars[[char_rd]] <- optim_cP
  }
  if (is.null(.get_cP_stuff(ranFix,"longdep",which=char_rd))) { 
    longdep <- .get_cP_stuff(init.optim,"longdep",which=char_rd)
    if (is.null(longdep)) longdep <- 0.5
    init$corrPars[[char_rd]] <- .modify_list(init$corrPars[[char_rd]], list(longdep=longdep)) ## canonical scale
    optim_cP <- init.optim$corrPars[[char_rd]]
    if (optim.scale=="transformed") {
      optim_cP <- .modify_list(optim_cP, list(trLongdep=.longdepFn(longdep,LDMAX=LDMAX))) 
      optim_cP$longdep <- NULL
    } else optim_cP <- .modify_list(optim_cP, list(longdep=longdep)) 
    init.optim$corrPars[[char_rd]] <- optim_cP
  }
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}



.calc_inits_Auto_rho <- function(init,init.optim,init.HLfit,ranFix,rhorange,For,user.lower,user.upper,char_rd) {
  if (is.null(.get_cP_stuff(ranFix,"rho",which=char_rd))) { ## $rho NULL or NA
    init_rho <- .get_cP_stuff(init.HLfit,"rho",which=char_rd)
    if (For=="corrHLfit") { ## with corrHLfit, defaul is outer optim. We provide init.optim
      if (is.null(init_rho)) { ## default init.HLfit
        init_rho <- .get_cP_stuff(init.optim,"rho",which=char_rd)
        if (! is.numeric(init_rho)) init.optim$corrPars[[char_rd]] <- init_rho <- list(rho=mean(rhorange))
      } else if ( ! is.numeric(init_rho) ) { ## non-default: inner optim, but no numeric init  
        init.HLfit$corrPars[[char_rd]] <- init_rho <- list(rho=mean(rhorange))
        init.optim$corrPars[[char_rd]] <- NULL 
      } else { ## non-default: inner optim, already with init value
        init.optim$corrPars[[char_rd]] <- NULL 
      }
    } else if (For=="fitme") { ## with fitme, default is inner optim. We provide init.HLfit
      if (is.null(init_rho)) {
        init_rho <- .get_cP_stuff(init.optim,"rho",which=char_rd)
        if ( ! is.numeric(init_rho) ) { ## default: outer optim, but no numeric init
          init.optim$corrPars[[char_rd]] <- init_rho <- list(rho=mean(rhorange))
          init.HLfit$corrPars[[char_rd]] <- NULL
        } else { ## non-default: inner optim, already with init value
          init.HLfit$corrPars[[char_rd]] <- NULL
        }
      } else if (! is.numeric(init_rho)) init.HLfit$corrPars[[char_rd]] <- init_rho <- list(rho=mean(rhorange))
    }
    init$corrPars[[char_rd]] <- list(rho=init_rho) ## bugs corrected in v2.3.30 and (again) v2.4.51
  }
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

.calc_inits_geostat_rho <- function(init,init.optim,init.HLfit,ranFix,maxrange,optim.scale,RHOMAX,user.lower,user.upper,char_rd) {
  if (is.null(.get_cP_stuff(ranFix,"rho",which=char_rd))) {
    rho <- .get_cP_stuff(init.optim,"rho",which=char_rd)
    if (is.null(rho)) {
      rho <- 30/(2*maxrange)
    } else if (any( narho <- is.na(rho))) rho[narho] <- 30/(2*maxrange[narho]) ## 05/2015 allows init.corrHLfit=list(rho=rep(NA,...
    rho <- min(max(rho, ## checking against user-provided min/max
                   user.lower$corrPars[[char_rd]]$rho)*1.000001, 
               user.upper$corrPars[[char_rd]]$rho/1.000001)
    # Info in canonical scale:
    init$corrPars[[char_rd]] <- .modify_list(init$corrPars[[char_rd]], list(rho=rho)) ## synthesis of user init and default init
    # Template of what will affect init.optim$corrPars in transformed scale:
    optim_cP <- init.optim$corrPars[[char_rd]] 
    if (optim.scale=="transformed") {
      optim_cP <- .modify_list(optim_cP, list(trRho=.rhoFn(rho,RHOMAX=RHOMAX)))
      optim_cP$rho <- NULL
    } else optim_cP <- .modify_list(optim_cP, list(rho=rho))
    init.optim$corrPars[[char_rd]] <- optim_cP
  } 
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

.calc_inits_IMRF <- function(init,init.optim,init.HLfit,ranFix,optim.scale,moreargs_rd,user.lower,user.upper,char_rd) {
  ## This is based on char_rd hence has to ignore hyper, unless contrived use of hyper"s attributes is implemented;
  ## HENCE the return value may contain more values than ultimately retained based on hyper info.
  ## For independent IMRF terms, the user input may be used; then, transformed user input is ignored, as in .calc_inits_dispPars()
  if (is.null(.get_cP_stuff(ranFix,"kappa",which=char_rd))) { 
    kappa <- .get_cP_stuff(init.optim,"kappa",which=char_rd)
    if (is.null(kappa)) { kappa <- 0.5 }
    kappa <- min(max(kappa, ## checking against user-provided min/max
                   user.lower$corrPars[[char_rd]]$kappa)*1.000001, 
               user.upper$corrPars[[char_rd]]$kappa/1.000001)
    # Info in canonical scale:
    init$corrPars[[char_rd]] <- .modify_list(init$corrPars[[char_rd]], list(kappa=kappa)) ## synthesis of user init and default init
    # Template of what will affect init.optim$corrPars in transformed scale:
    optim_cP <- init.optim$corrPars[[char_rd]] 
    if (optim.scale=="transformed") {
      optim_cP <- .modify_list(optim_cP, list(trKappa=.kappaFn(kappa,KAPPAMAX=moreargs_rd$KAPPAMAX)))
      optim_cP$kappa <- NULL
    } else optim_cP <- .modify_list(optim_cP, list(kappa=kappa))
    init.optim$corrPars[[char_rd]] <- optim_cP
  } 
  if (is.null(ranFix$lambda) || is.na(ranFix$lambda[char_rd])) { 
    if (is.null(moreargs_rd$minKappa)) { # slightly cryptic way of excludig spde case
      init$lambda[char_rd] <- 200 # but this will be ignored by optimize (i.e., in the comparison to LatticeKrig)
      init.optim$lambda[char_rd] <- 200 
    } # else control of initial value should be as by default method (~that for Matern)
  }
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

.calc_inits <- function(init.optim,init.HLfit,ranFix,corr_info,
                        moreargs,
                        user.lower,user.upper,
                       optim.scale, For, hyper_info) { 
  # final $ranFix should be suitable for .modify_list()... hence may contain NA's (in trLambda notably: cf mixing of ["1"] optimized, ["2"] fixed)
  # final $init.optim should be suitable for optimizer
  # final $[canon.]init should track $init.optim. It start different from $init.optim, though.
  init.optim$lambda <- .rename_lambda(init.optim$lambda) 
  ranFix$lambda <- .rename_lambda(ranFix$lambda) # Now, to allow name matching in .calc_inits_IMRF()
  inits <- list(init=list(corrPars=list()), ## minimal structure for assignments $corrPars[[char_rd]][[<name>]] to work.
                init.optim=init.optim, init.HLfit=init.HLfit, ranFix=ranFix,
                user.lower=user.lower, user.upper=user.upper, init=list()) 
  corr_types <- corr_info$corr_types
  for (rd in seq_along(corr_types)) {
    corr_type <- corr_types[rd]
    if (! is.na(corr_type)) {
      char_rd <- as.character(rd)
      inits <- corr_info$corr_families[[rd]]$calc_inits(inits=inits, char_rd=char_rd, moreargs_rd=moreargs[[char_rd]], 
                                                        user.lower=user.lower, user.upper=user.upper, optim.scale=optim.scale, 
                                                        ranFix=ranFix, 
                                                        init.optim=init.optim, For=For)
    }
  }
  inits <- .calc_inits_hyper(inits, hyper_info=hyper_info, fixed=inits$ranFix, moreargs=moreargs)
  # phi, lambda
  inits <- .calc_inits_dispPars(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                                user.lower=user.lower,user.upper=user.upper)
  # random coefficients
  inits <- .calc_inits_ranCoefs(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                                user.lower=user.lower,user.upper=user.upper)
  # GLM family parameters
  inits$init$COMP_nu <- inits$init.optim$COMP_nu ## may be NULL. No checks needed
  inits$init$NB_shape <- inits$init.optim$NB_shape ## may be NULL. No checks needed
  if (! is.null(inits$init.optim$NB_shape)) {
    inits$init.optim$trNB_shape <- .NB_shapeFn(inits$init.optim$NB_shape)
    inits$init.optim$NB_shape <- NULL
  }
  inits <- eval(inits) ## not sure why
  if ( ! length(inits[["init.optim"]][["corrPars"]])) { ## remove corrParslist()
    inits[["init.optim"]]["corrPars"] <- NULL
    inits[["init"]]["corrPars"] <- NULL
  }
  return(inits)
}
