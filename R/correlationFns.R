.HL_process.args <- function (...) { ## from gsl package
  a <- list(...)
  #attr <- attributes(a[[which.max(unlist(lapply(a, length)))]])
  a <- lapply(a, as.vector)
  out <- do.call(rbind, a)
  out <- list(`1`=out[1,],`2`=out[2,]) ## == out <- split(out, row(out)) but profiling shows the latter is slow
  names(out) <- paste("arg", 1:length(a), sep = "")
  #return(c(out, attr = list(attr)))
  return(out)
}

"MaternCorr" <- function(d, rho=1, smoothness, nu=smoothness, Nugget=0L) UseMethod("MaternCorr") 

#Matern.corr <- MaternCorr ## for back compat as it is in spMMjob.R

MaternCorr.default <- function (d, rho=1, smoothness, nu=smoothness, Nugget=NULL) { ## rho is alpha in fields
  ## ideally (but not necess) on a 'dist' so the diagonal is not  manipulated 
  if (any(d < 0)) 
    stop("distance argument must be nonnegative")
  dscal <- d * rho
  isd0 <- d == 0L
  dscal[isd0] <- 1e-10 ## avoids errors on distance =0; but the resulting corrvals can be visibly < 1 for small nu ## FIXME make value dependent on rho, nu ?
  logcon <- (nu - 1)*log(2)+ lgamma(nu) 
  ## base::besselK() is fast but not very accurate; 
  ## gsl::bessel_lnKnu() [v2.1-6] is sensitively slower; I derived from an older version .bessel_lnKnu() which is in-between;
  ## Bessel::BesselK() is much slower
  # corrvals <- - logcon + nu*log(dscal)+ gsl::bessel_lnKnu(x=dscal, nu=nu)  
  # corrvals <- - logcon + nu*log(dscal)+ log(Bessel::BesselK(z=dscal,nu=nu,expon.scaled = TRUE))-dscal 
  # corrvals <- - logcon + nu*log(dscal)+ .bessel_lnKnu(x=dscal, nu=nu) 
  corrvals <- - logcon + .nuln_plus_bessel_lnKnu(x=dscal, nu=nu) # dist of matrix dscal -> still vector return value
  corrvals <- exp(corrvals) 
  if ( ! is.null(Nugget)) corrvals[!isd0] <- (1-Nugget)* corrvals[!isd0]
  corrvals[isd0] <- 1 ## 
  corrvals[corrvals < 1e-16] <- 0 ## an attempt to deal with problem in chol/ldl/svd which don't like 'nearly-identity' matrices
  attributes(corrvals) <- attributes(dscal) ## subtle :-) handling of both matrix and dist objects 
  attr(corrvals,"corr.model") <- "Matern"
  return(corrvals)
}

.MaternCorr_hasSlotx <- function(d, rho=1, smoothness, nu=smoothness, Nugget=NULL) { ## rho is alpha in fields
  # such that NA's are treated as zero distances and implicit 0's as infinite distances 
  corrvals <- d
  x <- d@x
  x[is.na(x)] <- 0 # for "dsCDIST" case (maybe it should not longer by dsCDIST after this operation?)
  corrvals@x <- MaternCorr.default(x, rho, smoothness, nu, Nugget)
  attr(corrvals,"corr.model") <- "Matern"
  attr(corrvals,"pars") <- list(rho=rho, nu=nu, Nugget=Nugget) # for diagnostic purposes
  return(corrvals)
}

MaternCorr.dsCMatrix <- .MaternCorr_hasSlotx

MaternCorr.dgCMatrix <- .MaternCorr_hasSlotx


## ess <- function(nu,d) {exp(-(d/(2*sqrt(nu)))^2)} ...

"CauchyCorr" <- function(d, rho=1, shape, longdep, Nugget=NULL) UseMethod("CauchyCorr") 

CauchyCorr.default <- function(d, rho=1, shape, longdep, Nugget=NULL) { ## rho is 1/c in Gneiting
  ## ideally (but not necess) on a 'dist' so the diagonal is not  manipulated 
  if (any(d < 0)) 
    stop("distance argument must be nonnegative")
  dscal <- d * rho
  isd0 <- d == 0L
  corrvals <- 1+dscal^shape
  corrvals <- corrvals^(-longdep/shape)
  if ( ! is.null(Nugget)) corrvals[!isd0] <- (1-Nugget)* corrvals[!isd0]
  corrvals[corrvals < 1e-16] <- 0L ## an attempt to deal with problem in chol/ldl/svd which don't like 'nearly-identity' matrices
  attr(corrvals,"corr.model") <- "Cauchy"
  return(corrvals)
}

.CauchyCorr_hasSlotx <- function(d, rho=1, shape, longdep, Nugget=NULL) { ## rho is alpha in fields
  # such that NA's are treated as zero distances and implicit 0's as infinite distances 
  corrvals <- d
  x <- d@x
  x[is.na(x)] <- 0
  corrvals@x <- CauchyCorr.default(d, rho, shape, longdep, Nugget)
  attr(corrvals,"corr.model") <- "Cauchy"
  return(corrvals)
}

CauchyCorr.dsCMatrix <- .CauchyCorr_hasSlotx

CauchyCorr.dgCMatrix <- .CauchyCorr_hasSlotx

#### demo
if (F) {
 mym <- matrix(c(3,2,1,2,3,2,1,2,3),ncol=3)
 cL <- t(chol(mym))
 cL %*% t(cL)
 LDL <- eigen(mym,symmetric=T)
 eL <- LDL$vectors %*% diag(sqrt(LDL$values))    
 eL %*% t(eL)
 seL <- LDL$vectors %*% diag(sqrt(LDL$values)) %*% t(LDL$vectors)   ## symmetric
 seL %*% t(seL)
}

# tcrossprod factor, Ltri or not.
mat_sqrt <- function(m=NULL, # coRRelation matrix
                       symSVD=NULL, try.chol=TRUE, condnum=1e12) { ## tcrossprod(mat_sqrt(X))=X
  ## cf return value: the code must compute 'L', and if the type of L is not chol, also 'corr d' and 'u'
  type <- NULL
  
  if (is.null(symSVD)) {
    dim_m <- dim(m)
    if (dim_m[1L]!=dim_m[2L]) { stop("matrix is not square") }
    if (try.chol) {
      L <- .silent_W_E(.wrap_Ltri_t_chol(m)) # That's t(Matrix::chol(corrMatrix)) for sparse matrix and an Rcpp_chol otherwise
      if (inherits(L,"simpleError")) {
        # assuming the error is due to high condition number:
        EEV <- extreme_eig(m, symmetric=TRUE, required=TRUE) # bc required=TRUE Is what the old mat_sqrt() effectively did
        e1 <- EEV[1]
        en <- EEV[2]
        if (en < -1e-4) {
          mess <- paste0("Matrix has suspiciously large negative eigenvalue(s): is it a valid correlation matrix?")
          warning(mess)
        }
        diagcorr <- (e1 - condnum * en)/(-1 + condnum + e1 - condnum * en) ## so that 
        # corrected condition number is exactly the condnum param, given the following correction.
        # Further if the input matrix has unit diag the corrected one also has unit diagonal.
        # For general covariance matrices here is not such result (but no straightforward rule for a corrected matrix to belong to a declared 'corr'Family )
        diagcorr <- max(0,diagcorr) 
        m <- (1-diagcorr)*m
        # OK for *m*atrix but less clear for sparse Matrix:
        #nc <- ncol(corrMatrix)
        #diagPos <- seq.int(1L,nc^2,nc+1L)
        #m[diagPos] <- m[diagPos] + diagcorr
        diag(m) <- diag(m) + diagcorr ## # all diag is corrected => added a constant diagonal matrix 
        L <- .silent_W_E(try(.wrap_Ltri_t_chol(m)))
      } 
      if ( ! inherits(L,"simpleError") ) type <- "cholL_LLt" ## else type remains NULL      
    }
    if ( is.null(type) ) { ## no chol or failed chol
      decomp <- eigen(m, symmetric=TRUE) ## such that v= t(u) without any sign issue
      svdnames <- names(decomp)
      svdnames[svdnames=="values"] <- "d"
      svdnames[svdnames=="vectors"] <- "u"
      names(decomp) <- svdnames
      type <- "eigen" 
      if (any(decomp$d < -1e-08)) { # wrong test if R's svd had been used as its $d are the eigenvalues up to the sign, 
        # and once found <0 for a positve eigenvalue of a 2x2 matrix.
        message("correlation matrix has suspiciously large negative eigenvalue(s).")
        return(structure("correlation matrix has suspiciously large negative eigenvalue(s).",
                         class="try-error", ## passes control to calling function to print further info
                         call=match.call())) ## not fully conformant to "try-error" class ?"
      } else { 
        decomp$d[decomp$d<1e-16] <- 1e-16
        L <- .ZWZt(decomp$u,sqrt(decomp$d))
      }
    } 
    colnames(L) <- rownames(L) <- colnames(m) ## for checks in HLfit ## currently typically missing from symSVD  
  } else {
    decomp <- symSVD[c("d","u","adjd","eigrange")] # keep any $adjd  useful for SEM CAR; otherwise may be NULL 
    decomp$d[decomp$d<1e-16] <- 1e-16
    L <- .ZWZt(symSVD$u,sqrt(symSVD$d))  
    type <- "symsvd"      
  }
  attr(L,"type") <- type
  attr(L,"corr.model") <- attr(m,"corr.model")
  if (type != "cholL_LLt") attr(L,type) <- decomp
  return(L)
} 

# Speculative: not used. 
.safe_condnum <- function(M, symmetric) {
  condnum <- NULL
  if ((inherits(M,"sparseMatrix") &&  ncol(M)<2000L)
      ||  ncol(M)<1000L) decomp <- .try_RSpectra(M, symmetric=symmetric) # 1000 -> 0.28s
  if (is.null(decomp)) {  # RSpectra was not available or probem was too large
    # kappa() computes by default (an estimate of) the 2-norm condition number of a matrix or of 
    # the R matrix of a QR decomposition, perhaps of a linear fit. The 2-norm condition number can 
    # be shown to be the ratio of the largest to the smallest *non-zero* singular value of the matrix.
    if (ncol(M)<1000L) condnum <- kappa(M) # 1000 -> 0.5s
  } else condnum <- decomp$eigrange[2]/decomp$eigrange[1]
  if (is.null(condnum)) { # too large sparse or dense problem
    # Matrix::condest() use RNG, but we do not want any global effect on RNG! 
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(123) # makes the condest() call deterministic
    condnum <- Matrix::condest(M)$est 
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  condnum
}

extreme_eig <- function(M, symmetric, required=TRUE) {
  EEV <- decomp <- NULL
  if (ncol(M)>2L && (
    required ||
    (inherits(M,"sparseMatrix") &&  ncol(M)<2000L) ||
    ncol(M)<1000L
    ) 
  ) decomp <- .try_RSpectra(M,symmetric=symmetric) # 1000 -> 0.28s   # but increases fast: squirrels...
  if (is.null(decomp)) {  # RSpectra was not available or it failed or problem was too large ... or 2x2 matrix...
    # we cannot use kappa() or Matrix::condest() here since we really want the two extreme eigenvalues, not simply the condnum
    if (required ||
        ncol(M)<1000L) {
      values <- eigen(M,symmetric=symmetric, only.values = TRUE)$values
      EEV <- c(values[1],tail(values,1))
    }
  } else EEV <- c(decomp$eigrange[2],decomp$eigrange[1])
  EEV # in eigen order: largest first. remains NULL only if not 'required' and problem too large
}

regularize <- function(A, EEV=extreme_eig(A,symmetric=TRUE), maxcondnum=1e12#, minmax=0
                       ) {
  if (EEV[1]/EEV[2] > maxcondnum || EEV[2]<0)  {
    diagcorr <- (EEV[1] - maxcondnum * EEV[2])/(-1 + maxcondnum + EEV[1] - maxcondnum * EEV[2])  
    A <- (1-diagcorr)*A
    diag(x=A) <- diag(x=A) + diagcorr ## # all diag is corrected => added a constant diagonal matrix 
#  } else if (EEV[1]<minmax) { ## matrix very close to zero matrix (default is not to correct)
#    diag(x=A) <-  diag(x=A)+1/maxcondnum # effect is different. condnum ~ 1
  }
  A
}

as_precision <- function(corrMatrix, condnum=1e12) {
  # Compared to older as_precision: 
  # * no systematic eigen computation
  # * correction exactly achieves condnum (as corr_sqrt)
  # default condnum 1e12 (as mat_sqrt and corr_sqrt)
  # Compared to mat_sqrt: here input in coRR matrix, output too.
  if (inherits(corrMatrix,"dist")) { corrMatrix <- proxy::as.matrix(corrMatrix, diag=1) }
  corrMatrix <- forceSymmetric(corrMatrix)
  precmat <- tryCatch(chol2inv(chol(corrMatrix)),error=function(e) e)
  if (inherits(precmat,"simpleError")) {
    if (TRUE) {
      EEV <- extreme_eig(corrMatrix, symmetric=TRUE, required=TRUE) # bc required=TRUE Is what the old mat_sqrt() effectively did
      e1 <- EEV[1]
      en <- EEV[2]
      if (en < -1e-4) {
        mess <- paste0("Matrix has suspiciously large negative eigenvalue(s): is it a valid correlation matrix?")
        warning(mess)
      }
      diagcorr <- (e1 - condnum * en)/(-1 + condnum + e1 - condnum * en) ## so that 
      # corrected condition number is exactly the condnum param, given the following correction.
      # Further if the input matrix has unit diag the corrected one also has unit diagonal.
      # For general covariance matrices here is not such result (but no straightforward rule for a corrected matrix to belong to a declared 'corr'Family )
      diagcorr <- max(0,diagcorr) 
    } else { # (corrected) older version
      esys <- eigen(corrMatrix, only.values = TRUE, symmetric=TRUE)
      evalues <- esys$values
      min_d <- evalues[1L]/condnum ## so that corrected condition number is at most the 1e14
      diagcorr <- max(c(0,min_d-evalues)) # SINGLE SCALAR
    }
    corrMatrix <- (1-diagcorr)*corrMatrix
    diag(corrMatrix) <- diag(corrMatrix) + diagcorr ## # all diag is corrected => added a constant diagonal matrix 
    # oldMDCopt <- options(Matrix.warnDeprecatedCoerce = 0) # chol2inv(<dtC>) problem in Matrix v1.4.2
    precmat <- chol2inv(chol(corrMatrix)) # as explained in first version
    # options(oldMDCopt)
  }
  colnames(precmat) <- rownames(precmat) <- colnames(corrMatrix) 
  # drop0 useful to convert to sparseMarix even if no 0 to drop
  return(structure(list(matrix=drop0(precmat)),class=c("list","precision"))) # return value must be sparse, not simply Matrix. 
}


make_scaled_dist <- local({
  Chord_warned <- FALSE
  function(uniqueGeo,uniqueGeo2=NULL,distMatrix,rho,rho.mapping=seq_len(length(rho)),
                               dist.method="Euclidean",return_matrix=FALSE) {
  if (length(rho)>1L && dist.method!="Euclidean") { 
    stop("'rho' length>1 not allowed for non-Euclidean distance.")
  }
  if ( missing(distMatrix) ) { ## ## then provide it (most of this fn's code)
    if ( missing(uniqueGeo) ) {
      stop("missing(distMatrix) && missing(uniqueGeo).")
    } 
    distnrow <- nrow(uniqueGeo)
    distncol <- NROW(uniqueGeo2)
    scaled.dist <- NULL
    if (dist.method=="Euclidean") {
      if (length(rho)==1L) {
        uniqueScal <- uniqueGeo * rho 
      } else if (ncol(uniqueGeo)==length(rho.mapping)) {
        uniqueScal <- t(t(uniqueGeo) * rho[rho.mapping]) ## valid for vectorial rho...
      } else {
        mess  <- paste("Invalid length(rho[rho.mapping]). Length should be either 1 or",ncol(uniqueGeo))
        stop(mess)
      }
      if (! is.null(uniqueGeo2)) {
        if (length(rho)==1L) {
          uniqueScal2 <- uniqueGeo2 * rho  
        } else uniqueScal2 <- t(t(uniqueGeo2) * rho[rho.mapping]) 
      } else {
        uniqueScal2 <- NULL
      }
      scaled.dist <- .dist_fn(x=uniqueScal,y=uniqueScal2, method=dist.method) 
    } else { ## not Euclidean
      if (dist.method=="Chord" && (! Chord_warned) ) {
        warning("NB: using dist's Chord distance on a circle, not EarthChord on a sphere", immediate. = TRUE )
        Chord_warned <<- TRUE
      }
      scaled.dist <- rho * .dist_fn(uniqueGeo,y=uniqueGeo2,method=dist.method)  
    }
  } else { ## distMatrix provided
    scaled.dist <- rho * distMatrix
  }
  ## Here .dist_fn must always returns a dist object, maybe dist(0) 
  ## but rho * dist(0) is numeric(0); we standardise it:
  if ( identical(scaled.dist, numeric(0))) scaled.dist <- dist(0)
  if (return_matrix) {
    if (inherits(scaled.dist,"dist")) {
      scaled.dist <- as.matrix(scaled.dist)
    } else if (inherits(scaled.dist,"crossdist")) scaled.dist <- scaled.dist[] ## []: same effect as what one would expect from non-existent as.matrix.crossdist()
  }
  return(scaled.dist)
  }
}) 

getDistMat <- function(object,scaled=FALSE, which=1L) {
  if (! is.null(msd_arglist <- attr(object$strucList[[which]],"msd.arglist"))) {
    if (is.null(msd_arglist$distMatrix)) { ## we reconstruct it
      olduniqueGeo <- .get_old_info_uniqueGeo(object, char_rd=as.character(which)) 
      coordinates <- colnames(olduniqueGeo)
      msd_arglist$distMatrix <- eval(msd_arglist$distcall,
                                   list(uniqueGeo=olduniqueGeo,dist.method=msd_arglist$dist.method))
      msd_arglist$distcall <- NULL 
    }
    if ( ! scaled)  {
      msd_arglist$rho <- 1 
      msd_arglist$`rho.mapping` <- NULL 
    }
    return(do.call(make_scaled_dist,msd_arglist))
  } else {
    message("no random effect for which a distance matrix can be constructed.") # Matern: seek "msd.arglist"
    return(NULL)
  }
}
