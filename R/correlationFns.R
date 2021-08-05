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

.HL_strictify <- function (val, status) { ## from gsl package
    val[status < 0] <- NaN ##FR inverted the test // gsl package !!
    return(val)
}

.nuln_plus_bessel_lnKnu <- function (nu, x) { 
  jj <- .HL_process.args(nu, x) ## converts (nu, distance matrice to (nu vector, distance vector)
  nu.vec <- jj$arg1
  x.vec <- jj$arg2
  jj <- .nuln_plus_bessel_lnKnu_e(jj$arg1,jj$arg2) 
  # jj <- .nuln_plus_bessel_lnKnu_e_boost(jj$arg1,jj$arg2)  # boost attempt (v3.8.11) was quite correct but slower
  val <- .HL_strictify(jj$val, jj$status)
  attributes(val) <- attributes(x) ## FIXME is this useful here ?
  return(val)
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
  # corrvals <- - logcon + nu*log(dscal)+ log(gsl::bessel_lnKnu(x=dscal, nu=nu)) 
  # corrvals <- - logcon + nu*log(dscal)+ log(Bessel::BesselK(z=dscal,nu=nu,expon.scaled = TRUE))-dscal 
  # corrvals <- - logcon + nu*log(dscal)+ .bessel_lnKnu(x=dscal, nu=nu) 
  corrvals <- - logcon + .nuln_plus_bessel_lnKnu(x=dscal, nu=nu) 
  corrvals <- exp(corrvals) 
  if ( ! is.null(Nugget)) corrvals[!isd0] <- (1-Nugget)* corrvals[!isd0]
  corrvals[isd0] <- 1 ## 
  corrvals[corrvals < 1e-16] <- 0 ## an attempt to deal with problem in chol/ldl/svd which don't like 'nearly-identity' matrices
  attr(corrvals,"corr.model") <- "Matern"
  return(corrvals)
}


## ess <- function(nu,d) {exp(-(d/(2*sqrt(nu)))^2)} ...

CauchyCorr <- function(d, rho=1, shape, longdep, Nugget=NULL) { ## rho is 1/c in Gneiting
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
mat_sqrt <- function(m=NULL, symSVD=NULL, try.chol=TRUE, condnum=1e12) { ## tcrossprod(mat_sqrt(X))=X
  ## cf return value: the code must compute 'L', and if the type of L is not chol, also 'corr d' and 'u'
  type <- NULL
  
  if (is.null(symSVD)) {
    dim_m <- dim(m)
    if (dim_m[1L]!=dim_m[2L]) { stop("matrix is not square") }
    if (try.chol) {
      L <- try(.wrap_Ltri_t_chol(m),silent=TRUE)
      if (inherits(L,"try-error")) {
        blob <- eigen(m, symmetric=TRUE,only.values = TRUE) ## COV= blob$u %*% diag(blob$d) %*% t(blob$u) ##
        if (min(blob$values)< -1e-4) {
          mess <- paste0("Matrix has suspiciously large negative eigenvalue(s): is it a valid correlation matrix?")
          warning(mess)
        }
        min_d <- max(blob$values)/condnum ## so that corrected condition number is at most the denominator
        diagcorr <- max(c(0,min_d-blob$values)) # all diag is corrected => added a constant diagonal matrix to compactcovmat
        nc <- ncol(m)
        diagPos <- seq.int(1L,nc^2,nc+1L)
        m[diagPos] <- m[diagPos] + diagcorr
        L <- try(.wrap_Ltri_t_chol(m),silent=TRUE)
      } 
      if ( ! inherits(L,"try-error") ) type <- "cholL_LLt" ## else type remains NULL      
    }
    if ( is.null(type) ) { ## no chol or failed chol
      decomp <- eigen(m, symmetric=TRUE) ## such that v= t(u) without any sign issue
      svdnames <- names(decomp)
      svdnames[svdnames=="values"] <- "d"
      svdnames[svdnames=="vectors"] <- "u"
      names(decomp) <- svdnames
      type <- "eigen" 
      if (any(decomp$d < -1e-08)) {
        ## could occur for two reasons: wrong input matrix; or problem with R's svd which $d are the eigenvalues up to the sign, 
        ##   and thus $d can be <0 for eigenvalues >0 (very rare, observed in a 2x2 matrix) 
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
      scaled.dist <- proxy::dist(x=uniqueScal,y=uniqueScal2, method=dist.method) 
    } else { ## not Euclidean
      if (dist.method=="Chord" && (! Chord_warned) ) {
        warning("NB: using dist's Chord distance on a circle, not EarthChord on a sphere", immediate. = TRUE )
        Chord_warned <<- TRUE
      }
      scaled.dist <- rho * proxy::dist(uniqueGeo,y=uniqueGeo2,method=dist.method)  
    }
  } else { ## distMatrix provided
    scaled.dist <- rho * distMatrix
  }
  ## Here proxy::dist always returns a dist object, maybe dist(0) 
  ## but rho * dist(0) is numeric(0); we standardise it:
  if ( identical(scaled.dist, numeric(0))) scaled.dist <- dist(0)
  if (return_matrix) {
    if (inherits(scaled.dist,"dist")) {
      scaled.dist <- as.matrix(scaled.dist)
    } else if (inherits(scaled.dist,"crossdist")) scaled.dist <- scaled.dist[] ## []: same effect as what oen would expect from non-existent as.matrix.crossdist()
  }
  return(scaled.dist)
  }
}) 

getDistMat <- function(object,scaled=FALSE, which=1L) {
  if (! is.null(msd_arglist <- attr(object$strucList[[which]],"msd.arglist"))) {
    if (is.null(msd_arglist$distMatrix)) { ## we reconstruct it
      info_olduniqueGeo <- attr(object,"info.uniqueGeo") 
      if ( ! is.array(info_olduniqueGeo)) { ## test TRUE for version > 2.3.18:
        olduniqueGeo <- info_olduniqueGeo[[as.character(which)]] ## spaMM3.0 but could be better documented
      } else olduniqueGeo <- info_olduniqueGeo 
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
