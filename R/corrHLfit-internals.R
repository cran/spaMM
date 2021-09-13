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

.rcDispFn <- function(v) log(v) # log(log1p(v)) # .dispFn(v) # 
.rcDispInv <- function(v) exp(v) #  exp(exp(v))-1 # .dispInv(v) # 

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
# Called only by .ranCoefsFn() in devel code, for optim bound, and to convert back the optim results for refit
.sphFn <- function(covpars,Xi_ncol=NULL) { ## from *cov* parameter in lower.tri order vector to trRanCoefs
  if (is.null(Xi_ncol)) Xi_ncol <- floor(sqrt(length(covpars)*2))
  cholmat <- diag(Xi_ncol)
  cholmat[lower.tri(cholmat,diag = TRUE)] <- covpars
  # ~ cov2cor but on half matrix (on which cov2cor() gives an incorrect result)
  lambdas <- diag(cholmat)
  tocorr <- diag(x=sqrt(1/lambdas))
  cholmat <- tocorr %*% cholmat %*% tocorr
  diag(cholmat) <- diag(cholmat)/2
  cholmat <- t(cholmat) + cholmat
  esys <- .eigen_sym(cholmat)
  crossfac <- .Dvec_times_matrix(sqrt(esys$values), t(esys$vectors))
  qrblob <- qr(crossfac)
  cholmat <- qr.R(qrblob) # applying .lmwithQR() systematically (badly) affects numerical precision
  if (! all(unique(diff(qrblob$pivot))==1L)) { # eval an unpermuted triangular R
    cholmat <- .lmwithQR(cholmat[, sort.list(qrblob$pivot)] ,yy=NULL,returntQ=FALSE,returnR=TRUE)$R
  } 
  # sequel assumes UPPER triangular cholmat (cholmat[i:1, i]...)
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
  trRancoef <- c(.rcDispFn(lambdas)/2, corrpars) ## trRancoef:= unconstrained parameterization of the cov matrix with log sigma in first positions
  if(any(is.nan(trRancoef))) stop("any(is.nan(trRancoef))")
  return(trRancoef) 
}

.smooth_regul <- function(covmat, epsi= .spaMM.data$options$tol_ranCoefs_outer["regul"]) {
  es <- .eigen_sym(covmat)
  v <- es$values
  # correction <- v*(v/(4*epsi)-1)+epsi # goes from epsi to 0 as v goes from 0 to 2*epsi. d/dv = -1 in 0 and =0 in 2*epsi.
  # v[v<2*epsi] <- (v+correction)[v<2*epsi] # *increasing* from epsi to 2*epsi as v goes from 0 to 2*epsi, derivable in 2*epsi.
  # i.e.:
  v[v<2*epsi] <- (4*epsi)*v[v<2*epsi]^2+epsi ## inverse is w[w<2*epsi] <- sqrt( (4*epsi)*(w[w<2*epsi]-epsi) )
  es$d_regul <- v
  covmat <- .ZWZt(es$vector,v)
  return(structure(covmat, esys=es))
}

# called only at the begining of inner or outer optim 
.chlFn <- function(covpars,Xi_ncol=NULL) { # from *cov* parameter in lower.tri order vector to trRanCoefs
  if (is.null(Xi_ncol)) Xi_ncol <- floor(sqrt(length(covpars)*2))
  svdv <- diag(Xi_ncol)
  .lower.tri(svdv,diag = TRUE) <- covpars
  svdv[upper.tri(svdv)] <- t(svdv)[upper.tri(svdv)] ## __F I X M E__ ugly...
  crossfac <- .Utri_chol_by_qr(svdv) ## upper.tri crossfac
  return(crossfac[upper.tri(crossfac,diag = TRUE)]) ## returned as vector 
}

.corrFn <- function(covpars,Xi_ncol=NULL) { ## from *cov* parameter in lower.tri order vector to trRanCoefs
  if (is.null(Xi_ncol)) Xi_ncol <- floor(sqrt(length(covpars)*2))
  covcorr <- diag(Xi_ncol)
  .lower.tri(covcorr,diag = TRUE) <- covpars
  lambdas <- diag(covcorr)
  covcorr[upper.tri(covcorr)] <- t(covcorr)[upper.tri(covcorr)]
  covcorr <- cov2cor(covcorr) # cov2cor works only o nthe full matrix
  corrpars <- covcorr[lower.tri(covcorr)]
  trRancoef <- c(.rcDispFn(lambdas)/2, corrpars) ## trRancoef:= constrained parameterization of the cov matrix with log sigma in first positions
  if(any(is.nan(trRancoef))) stop("any(is.nan(trRancoef))")
  return(trRancoef) 
}

.calc_cov_from_trRancoef <- function(trRancoef, Xi_ncol=NULL, rC_transf) { 
  if (is.null(Xi_ncol)) Xi_ncol <- attr(trRancoef,"Xi_ncol")
  if (is.null(Xi_ncol)) Xi_ncol <- floor(sqrt(length(trRancoef)*2L))
  if (rC_transf=="chol") { # trRancoef in upper.tri order
    trRancoef <- trRancoef/.spaMM.data$options$rC_transf_fac
    lampos <- cumsum(seq(Xi_ncol)) 
    trRancoef[lampos] <-  pmax(1e-32, trRancoef[lampos]) # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    cholmat <- matrix(0, nrow=Xi_ncol, ncol=Xi_ncol) 
    .upper.tri(cholmat,diag = TRUE) <- trRancoef
    covmat <- crossprod(cholmat)
    attr(covmat,"chol_crossfac") <- cholmat
    return(covmat)
  } else if (rC_transf=="corr") { # trRancoef in (log(sigma), corr) order
    covmat <- diag(Xi_ncol)/2
    covmat[lower.tri(covmat)] <- trRancoef[-seq(Xi_ncol)]
    covmat <- covmat+t(covmat)
    sigmas <- sqrt(.rcDispInv(2*trRancoef[seq(Xi_ncol)])) 
    ds <- diag(x=sigmas)
    covmat <- ds %*% covmat %*% ds
    return(covmat)
  } else { # "sph"  trRancoef in (log(sigma), corr) order
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
    # 
    # corrmat <- solve(corrmat)
    # 
    sqrtlam <- diag(x=.rcDispInv(trRancoef[1:Xi_ncol])) # (this assumes sqrt() is used in transfo, which is true for sph)
    covmat <- sqrtlam %*% corrmat %*% sqrtlam
    if(any(is.nan(covmat))) stop("any(is.nan(covmat))")
    return(covmat)
  }
}

.calc_invL_from_trRancoef <- function(trRancoef, Xi_ncol=NULL,
                                      sing_threshold=.spaMM.data$options$invL_threshold,
                                      rC_transf) { 
  if (is.null(Xi_ncol)) Xi_ncol <- attr(trRancoef,"Xi_ncol")
  if (is.null(Xi_ncol)) Xi_ncol <- floor(sqrt(length(trRancoef)*2))
  if (rC_transf=="chol") {
    cholmat <- diag(nrow=Xi_ncol)
    .upper.tri(cholmat,diag = TRUE) <- trRancoef## upper tri crossfac (usual chol() convention)
    crossfac_precmat <- t(solve(cholmat)) ## lower tri crossfac
  } else { # "sph"
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
      cholmat_corr[seq(col-1L),col] <- cossign*sqrt(subox) ## upper tri
      cholmat_corr[col,col] <- sqrt(1-sum(subox))
    }
    inv_cholmat_corr <- solve(cholmat_corr) ## upper tri
    dvec <- .rcDispInv(-trRancoef[1:Xi_ncol]) #1/sigma
    crossfac_precmat <- t(.Dvec_times_matrix(dvec, inv_cholmat_corr)) # lower tri
  }
  #print(range(crossfac_precmat))
  #kappa(crossfac_precmat)>sing_threshold
  if (prod(.diagfast(crossfac_precmat,Xi_ncol))>sing_threshold) { ## F I X_invL any(.diagfast(crossfac_precmat,Xi_ncol)>sing_threshold) maybe not the good test
    return(NULL)
  } else return(crossfac_precmat)
}


# __F I X M E__ C version ? inelegant + repetitive call in ..process_ranCoefs() could be avoided sometimes... but profiling shows it's a non-issue
.ranCoefsFn <- function(vec, rC_transf) { # from canonical vector (var+corr) space in lower.tri order
  if ( ! is.null(vec)) {
    transf <- attr(vec,"transf")
    if ( ! is.null(transf)) { 
      attr(transf,"canon") <- as.vector(vec) # as.vector() to drop attributes
      attr(transf,"Xi_ncol") <- attr(vec,"Xi_ncol")
      return(transf)
    }
    Xi_ncol <- attr(vec,"Xi_ncol")
    if (is.null(Xi_ncol)) Xi_ncol <- floor(sqrt(length(vec)*2))
    #
    compactcovmat <- matrix(0,nrow=Xi_ncol,ncol=Xi_ncol)
    lowerbloc <- lower.tri(compactcovmat,diag=TRUE) ## a matrix of T/F !
    compactcovmat[lowerbloc] <- vec
    diagPos <- seq.int(1L,Xi_ncol^2,Xi_ncol+1L)
    lambdas <- compactcovmat[diagPos]
    # lambdas <- lambdas + 1e-10 # here that would affect the conversion from correlation to covariances
    sigmas <- diag(x=sqrt(lambdas)) 
    #
    compactcovmat[diagPos] <- 1
    compactcovmat <- sigmas %*% compactcovmat %*% sigmas ## using lower tri only
    if (rC_transf=="chol") {
      sc_chol <- .chlFn(compactcovmat[lowerbloc],Xi_ncol=Xi_ncol)*.spaMM.data$options$rC_transf_fac
      trRancoef <- structure(sc_chol,
                             canon=vec,
                             Xi_ncol=Xi_ncol)
    } else if (rC_transf=="corr") {
      trRancoef <- structure(.corrFn(compactcovmat[lowerbloc],Xi_ncol=Xi_ncol),
                             canon=vec,
                             Xi_ncol=Xi_ncol)
    } else {
      trRancoef <- structure(.sphFn(compactcovmat[lowerbloc],Xi_ncol=Xi_ncol),
                             canon=vec,
                             Xi_ncol=Xi_ncol)
    }
    
    return(trRancoef) # to transformed, 'unbounded' parameter space with log sigma in first positions
  } else return(NULL)
}

# obsolete: there is a  .C_ version of it
.calc_cov_from_ranCoef <- function(ranCoef, Xi_ncol=attr(ranCoef,"Xi_ncol")) {
  # assume input is marginal variances + correlation as in user input, but ordered as in lower.tri
  compactcovmat <- matrix(0,nrow=Xi_ncol,ncol=Xi_ncol)
  .lower.tri(compactcovmat,diag=TRUE) <- ranCoef
  diagPos <- seq.int(1L,Xi_ncol^2,Xi_ncol+1L)
  lambdas <- pmin(1e12,compactcovmat[diagPos]) # motivated by \ref{ares_bobyqa_ssprec}
  # lambdas <- lambdas + 1e-10 # here that would affect the conversion from correlation to covariances
  sigmas <- sqrt(lambdas) 
  compactcovmat <- (compactcovmat+t(compactcovmat))
  compactcovmat[diagPos] <- 1
  compactcovmat <- sigmas * compactcovmat * rep(sigmas, each = Xi_ncol)
  return(compactcovmat)
}  

# not called in inner optimization; only called is ranCoefsInv()
.tr2cov <- function(trRancoef,Xi_ncol=NULL, rC_transf) { ## from trRanCoefs to *cov* parameter vector
  if (is.null(Xi_ncol)) Xi_ncol <- floor(sqrt(length(trRancoef)*2))
  covmat <- .calc_cov_from_trRancoef(trRancoef,Xi_ncol=Xi_ncol, rC_transf=rC_transf) # must handle all transfos so there is no .chlInv() nor .sphInv()
  covmat <- .smooth_regul(covmat)
  return(covmat[lower.tri(covmat,diag = TRUE)]) ## vector representation of the cov matrix using covariances
}
#.tr2cov(.sphFn(c(1,1,1,5,5,14)))

.ranCoefsInv <- function(trRancoef, rC_transf) { # from transformed parameter space to correlation parameter vector
  if ( ! is.null(trRancoef)) {
    canon <- attr(trRancoef,"canon")
    if ( ! is.null(canon)) {
      attr(canon,"transf") <- as.vector(trRancoef) # as.vector() to drop attributes
      attr(canon,"Xi_ncol") <- attr(trRancoef,"Xi_ncol")
      return(canon)
    }
    if (rC_transf=="chol") return(.rC_inv_chol_cpp(trRancoef))
    Xi_ncol <- attr(trRancoef,"Xi_ncol")
    if (is.null(Xi_ncol)) Xi_ncol <- floor(sqrt(length(trRancoef)*2L))
    diagPos <- seq.int(1L,Xi_ncol^2,Xi_ncol+1L)
    if (FALSE) {
      covpars <- .tr2cov(trRancoef, Xi_ncol=Xi_ncol, rC_transf=rC_transf) # to vector of parameters with *cov*ariances
      # from vector with *cov*ariances to canonical vector with *corr*elations:
      lampos <- rev(length(covpars) -cumsum(seq(Xi_ncol))+1L)
      lambdas <- covpars[lampos]
      rancoefs <- diag(nrow=Xi_ncol)
      rancoefs[lower.tri(rancoefs,diag = TRUE)] <- covpars
      # lambdas <- lambdas + 1e-10 # here affecting also the conversion from covariances to correlations
    } else {
      covmat <- .calc_cov_from_trRancoef(trRancoef, Xi_ncol = Xi_ncol, rC_transf = rC_transf)
      rancoefs <- .smooth_regul(covmat)
      lambdas <- rancoefs[diagPos]
    }
    torancoefs <- sqrt(1/lambdas)
    rancoefs <- torancoefs * rancoefs * rep(torancoefs, each = Xi_ncol) # cf cov2cor()
    rancoefs[diagPos] <- lambdas
    varcorr <- rancoefs[lower.tri(rancoefs,diag = TRUE)]
    if (any( ! is.finite(varcorr))) stop(paste("Random-coefficient fitting failed. Trying spaMM.options(rC_transf='sph')", 
                                               "or spaMM.options(rC_transf_inner='sph') may be useful.", collapse="\n"))
    attr(varcorr,"transf") <- trRancoef
    attr(varcorr,"Xi_ncol") <- Xi_ncol
    return(varcorr) # to canonical vector (var+corr) space
  } else return(NULL)
}
# zut <- c(1.0000000,  0.4472136,  0.2672612,  5.0000000,  0.5976143, 14.0000000)
#(.ranCoefsInv(.ranCoefsFn(zut))[1:6]) # test requires [1:6] otherwise the attribute is used...

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
                                  adjust=.spaMM.data$options$max_bound_ranCoefs,
                                  rC_transf=.spaMM.data$options$rC_transf
                                  ) { 
  # !!!!!!!!! If I change the parametrisation I must also change the .xtol_abs_fn() code !!!!!!!!!
  if ( rC_transf=="chol") {
    if (.spaMM.data$options$augZXy_fitfn==".HLfit_body_augZXy_invL") { # only devel code AFAICS
      bnd1 <- diag(1e-10, Xi_ncol) # bc there is a solve on the implied chol factor
    } else bnd1 <- diag(0, Xi_ncol) # so not really chol factor, but still bounds a space of triangular factor containing chol factors
    if (tol_ranCoefs["inner"]) {
      cholmax <- tol_ranCoefs["cholmax"] ## Inf by default...
    } else cholmax <- Inf
    .upper.tri(bnd1, diag=FALSE) <- rep(-cholmax, Xi_ncol*(Xi_ncol-1L)/2L)  
    lower <- bnd1[upper.tri(bnd1,diag=TRUE)]
    upper <- rep(cholmax, Xi_ncol*(Xi_ncol+1L)/2L) 
    adjust_init <- NULL
  } else { ## "sph" etc.
    lower <- c(
      rep(.rcDispFn(tol_ranCoefs["lo_lam"]), Xi_ncol), # log(sigma)... or log(sigmafac) !!
      rep(-1+tol_ranCoefs["corr"], Xi_ncol*(Xi_ncol-1L)/2L) 
    )
    upper <- c(
      rep(.rcDispFn(tol_ranCoefs["up_lam"]), Xi_ncol), # log(sigma)... or log(sigmafac) !!
      rep(1-tol_ranCoefs["corr"], Xi_ncol*(Xi_ncol-1L)/2L) 
    )
    adjust_init <- list(lower=c(
      rep(.rcDispFn((tol_ranCoefs["lo_lam"])^(2/3)), Xi_ncol), #for log(sigma or sigmafac)
      rep(-1+tol_ranCoefs["corr"], Xi_ncol*(Xi_ncol-1L)/2L) 
    ))
  }
  lower <- adjust*lower
  lower <- pmin(lower, trRancoef)
  upper <- adjust*upper
  upper <- pmax(upper, trRancoef)
  #print(list(lower=lower,upper=upper))
  return(list(lower=lower,upper=upper,adjust_init=adjust_init)) ## transformed scale with log sigma in first positions
}

# obsolete:
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

## replacement for .expand_GeoMatrices() (<= ~3.8.23) producing sparser matrices 
# produces "dsCDIST" matrices where implicit 0's mean infinite distance, while Na on diagonal means 0 distance 
# seems to work for the 3 "algebra"s, but would be useful for spcorr
# Possible problems:
# vector 'rho' scale argument not checked
# production, in e.g. predict() of matrix types other that those handled by MaternCorr methods.
# test-nested-geostat for basic tests:
.sparse_expand_GeoMatrices <- function(w_uniqueGeo, e_uniqueGeo, coords_nesting, coord_within, dist.method,
                                       allow_DIST) {
  rownames(e_uniqueGeo) <- seq(nrow(e_uniqueGeo)) ## local & unique rownames
  ## Remarkably the next line works only if the cols are not factors ! Unless we have a fix for this,
  #  uniqueGeo classes should be integer not factor: see instances of as.numeric(levels(fac)) in the sources.
  rows_bynesting <- by(e_uniqueGeo ,e_uniqueGeo[,coords_nesting],rownames) # disjoints subsets of rownames since rownames are unique
  if (FALSE) { # assignments to blockrows (not in sequential order!) of sparse matrices is quite inefficient
    distMatrix <- Matrix(0,ncol=nrow(e_uniqueGeo),nrow =nrow(e_uniqueGeo)) 
    rownames(distMatrix) <- colnames(distMatrix) <- rownames(e_uniqueGeo) ## same trivial rownames
    for (lit in seq_len(length(rows_bynesting))) {
      blockrows <- rows_bynesting[[lit]]
      within_values <- e_uniqueGeo[blockrows,coord_within]
      w_dist <- proxy::dist(within_values,method=dist.method)
      distMatrix[blockrows,blockrows] <- as.matrix(w_dist)
    }
  } else {
    distMatrix <- vector("list", length(rows_bynesting))
    for (lit in seq_along(rows_bynesting)) {
      blockrows <- rows_bynesting[[lit]]
      within_values <- e_uniqueGeo[blockrows,coord_within]
      w_dist <- proxy::dist(within_values,method=dist.method)
      distMatrix[[lit]] <- as(as.matrix(w_dist),"dsCMatrix")
    }
    distMatrix <- Matrix::bdiag(distMatrix)
    rowperm <- sort.list(as.integer(unlist(rows_bynesting)))
    distMatrix <- distMatrix[rowperm,rowperm]
    rownames(distMatrix) <- colnames(distMatrix) <- rownames(e_uniqueGeo) ## trivial rownames -- not sure we need them now
  }
  if( ! inherits(distMatrix,"dsCMatrix")) stop("code missing in .sparse_expand_GeoMatrices()")
  if (allow_DIST) { # should be TRUE when the @x is to be used directly for conversion from distances to correlations
    diag(distMatrix) <- NA
    attr(distMatrix,"dsCDIST") <- TRUE # should mean 'dsC with NA on the diagonal'
  } else diag(distMatrix) <- 1
  # the diag<- are a bit slow, but there is no point in trying to assign on the blocks; 
  # moreover, proxy::as.matrix(,diag=NA) does not have the desired effect.
  return(list(distMatrix=distMatrix))
}

.expand_grouping <- function(w_uniqueGeo, e_uniqueGeo, coords_nesting, coord_within) {
  rownames(e_uniqueGeo) <- seq(nrow(e_uniqueGeo)) ## local & unique rownames
  ## Remarkably the next line works only if the cols are not factors ! Unless we have a fix for this,
  #  uniqueGeo classes should be integer not factor: see instances of as.numeric(levels(fac)) in the sources.
  rows_bynesting <- by(e_uniqueGeo ,e_uniqueGeo[,coords_nesting],rownames) 
  notSameGrp <- matrix(TRUE,ncol=nrow(e_uniqueGeo),nrow =nrow(e_uniqueGeo))
  rownames(notSameGrp) <- colnames(notSameGrp) <- rownames(e_uniqueGeo) ## same trivial rownames
  for (lit in seq_len(length(rows_bynesting))) {
    blockrows <- rows_bynesting[[lit]]
    notSameGrp[blockrows,blockrows] <- FALSE
  }
  return(notSameGrp)
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
                                                        "corrHLfit" = FALSE, "fitme" = FALSE, "is_separated" = FALSE)
  return(verbose)
}

.calc_inits_ranCoefs <- function(init,init.optim,init.HLfit,ranFix,user.lower,user.upper) {
  if ( ! is.null(init.optim$ranCoefs)) { ## should always be a complete ordered list in the case for random-coefficient models
    if (! is.null(init$ranCoefs)) for (st in names(init$ranCoefs)) init.optim$ranCoefs[[st]] <- init$ranCoefs[[st]]
    init$ranCoefs <- init.optim$ranCoefs
    trRanCoefs <- lapply(init.optim$ranCoefs,.ranCoefsFn, rC_transf=.spaMM.data$options$rC_transf)
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
    } else if (For=="fitme" || For=="fitmv") { ## with fitme, default is inner optim. We provide init.HLfit
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
    # We need to distinguish IMRF(.| model=.) from multIMRF and currently the info seem to be only in $minKappa.
    # I tracked a bug where user's init.lambda was not obeyed bc the wrong moreargs_rd was passed
    # One would have to modify IMRF$calc_moreArgs() to provide extra info allowing to define a more explicit condition?
    if (is.null(moreargs_rd$minKappa)) { # multIMRF: then we provide default init values for lambda.
      # But these default values are overwritten by any init$hyper value (i.e., init=list(hyper=list("1"=list(hy_lam=666))) works)
      # so the user retains control.
      init$lambda[char_rd] <- 200 # but this will be ignored by optimize (i.e., in the comparison to LatticeKrig)
      init.optim$lambda[char_rd] <- 200 
    } # else case IMRF(.| model=.): control of initial lambda value should be done by default method for scalar lambdas (~that for Matern)
  }
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

.calc_inits <- function(init.optim,init.HLfit,ranFix,corr_info,
                        moreargs,
                        user.lower,user.upper,
                       optim.scale, For, hyper_info) { 
  # final $ranFix should be suitable for .modify_list()... hence may contain NA's (in trLambda notably: cf mixing of ["1"] optimized, ["2"] fixed)
  # final $init.optim should be suitable for optimizer
  # final $[canon.]init should track $init.optim. It starts different from $init.optim, though.
  init.optim$lambda <- .rename_lambda(init.optim$lambda) 
  ranFix$lambda <- .rename_lambda(ranFix$lambda) # Now, to allow name matching in IMRF()$calc_inits()
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
  # If inits are not forced to be in user bounds at this step, .safe_opt() will detect that init is 'too close to bounds',
  # and will use bobyqa with locally corrected init that happen to satisfy the bounds. The optimInfo bears hardly any trace of that:
  # It does not show the locally corrected init, and one can only guess from it that nloptr() was not used from the 'optr' format.
  # GLM family parameters
  if (is.null(COMP_nu <- inits[["init"]]$COMP_nu)) COMP_nu <- inits[["init.optim"]]$COMP_nu # 'if user did not provide any, use the default one present in init.optim'
  if ( ! is.null(COMP_nu)) {
    if ( ! is.null(user.upper$COMP_nu)) COMP_nu <- min(user.upper$COMP_nu,COMP_nu) # mv: this is only called on submodels, so pmin, pman NOT needed.
    if ( ! is.null(user.lower$COMP_nu)) COMP_nu <- max(user.lower$COMP_nu,COMP_nu)
    inits[["init.optim"]]$COMP_nu <- inits[["init"]]$COMP_nu <- COMP_nu
  } # else I could add warnings that lower|upper values will be ignored... 
  if (is.null(NB_shape <- inits[["init"]]$NB_shape)) NB_shape <- inits[["init.optim"]]$NB_shape 
  if ( ! is.null(NB_shape)) {
    if ( ! is.null(user.upper$NB_shape)) NB_shape <- min(user.upper$NB_shape,NB_shape)
    if ( ! is.null(user.lower$NB_shape)) NB_shape <- max(user.lower$NB_shape,NB_shape)
    inits[["init"]]$NB_shape <- NB_shape
  }
  if (! is.null(inits[["init"]]$NB_shape)) {
    inits[["init.optim"]]$trNB_shape <- .NB_shapeFn(inits[["init"]]$NB_shape)
    inits[["init.optim"]]$NB_shape <- NULL
  }
  inits <- eval(inits) ## not sure why
  if ( ! length(inits[["init.optim"]][["corrPars"]])) { ## remove corrParslist()
    inits[["init.optim"]]["corrPars"] <- NULL
    inits[["init"]]["corrPars"] <- NULL
  }
  return(inits)
}

.eigen_sym <- function(X) {
  if (FALSE) { # transient effort to improve over eigen()
    svdv <- svd(X)
    chk <- with(svdv,sign(u)*sign(v)) 
    if (any(chk<0L)) { # well, perhaps the matrix is symmetric but not positive definite. Correction code handles this case (only):
      for (colit in seq_len(ncol(chk))) if (all(chk[,colit]<0L)) svdv$d[colit] <- - svdv$d[colit]
    }
    esys <- with(svdv,list(values=d,vectors=(u+sign(u*v)*v)/2)) # a bit heuristic, but this appears more accurate than simpler alternatives.
  } else {
    esys <- eigen(X,symmetric = TRUE) ## COV= .ZwZt(esys$vectors,esys$values) ##
    # => bare-bone version, which R CMD check does not like, but that passes the long tests:
    # z <- .Internal(La_rs(X, FALSE))
    # ord <- rev(seq_along(z$values))
    # esys <- structure(class = "eigen", list(values = z$values[ord], 
    #                                 vectors = z$vectors[, ord, drop = FALSE]))
  }
  return(esys)
}