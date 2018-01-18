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

.bessel_lnKnu <- function (nu, x, give = FALSE, strict = TRUE) { ## from bessel_lnKnu in gsl package
    jj <- .HL_process.args(nu, x) ## converts (nu, distance matrice to (nu vector, distance vector)
    nu.vec <- jj$arg1
    x.vec <- jj$arg2
    if (give) attr <- jj$attr
    jj <- .bessel_lnKnu_e(jj$arg1,jj$arg2)
    val <- jj$val
    err <- jj$err
    status <- jj$status
    if (strict) val <- .HL_strictify(val, status)
    attributes(val) <- attributes(x) ## FIXME is this useful here ?
    if (give) {
      attributes(status) <- attr
      attributes(err) <- attr
      return(list(val = val, err = err, status = status))
    } else {
      return(val)
    }
}

"MaternCorr" <- function(d, rho=1, smoothness, nu=smoothness, Nugget=0L) UseMethod("MaternCorr") 

Matern.corr <- MaternCorr ## for back compat as it is in spMMjob.R

MaternCorr.default <- function (d, rho=1, smoothness, nu=smoothness, Nugget=0L) { ## rho is alpha in fields
  ## ideally (but not necess) on a 'dist' so the diagonal is not  manipulated 
  if (any(d < 0)) 
    stop("distance argument must be nonnegative")
  dscal <- d * rho
  isd0 <- d == 0L
  dscal[isd0] <- 1e-10 ## avoids errors on distance =0; but the resulting corrvals can be visibly < 1 for small nu ## FR->FR make value dependent on rho, nu ?
  logcon <- (nu - 1)*log(2)+ lgamma(nu) 
  corrvals <- - logcon + nu*log(dscal)+ .bessel_lnKnu(x=dscal, nu=nu) ## 
  ##    corrvals <- - logcon + nu*log(dscal)+ log(besselK(x=dscal, nu=nu)) ## function from package gsl
  corrvals <- exp(corrvals) 
  corrvals[!isd0] <- (1-Nugget)* corrvals[!isd0]
  corrvals[isd0] <- 1 ## 
  corrvals[corrvals < 1e-16] <- 0L ## an attempt to deal with problem in chol/ldl/svd which don't like 'nearly-identity' matrices
  attr(corrvals,"corr.model") <- "Matern"
  return(corrvals)
}


## ess <- function(nu,d) {exp(-(d/(2*sqrt(nu)))^2)} ...


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

# .designL.from.Qmat <- function(Qmat) {
#   ## Cholesky gives proper LL' (think LDL')  while chol() gives L'L...
#   Q_CHMfactor <- Matrix::Cholesky(Matrix::drop0(Qmat),LDL=FALSE,perm=FALSE)
#   LMatrix  <- solve(Q_CHMfactor,system="Lt") ## solve(t(as(Q_CHMfactor,"sparseMatrix")))
#   # next line adds attributes to an S4 object. str() does not show these attributes... 
#   return(structure(LMatrix,type="invQ_CHMfactor",Q_CHMfactor=Q_CHMfactor,Qmat=Qmat))
# } 

## FR->FR we also use this function once on d2hdv2 in HLfit...
`designL.from.Corr` <- function(m=NULL, symSVD=NULL, try.chol=TRUE, try.eigen=FALSE,threshold=1e-06,SVDfix=1/10) {
  ## cf return value: the code must compute 'L', and if the type of L is not chol, also 'corr d' and 'u'
  type <- NULL
  
  if (is.null(symSVD)) {
    dim_m <- dim(m)
    if (dim_m[1L]!=dim_m[2L]) { stop("matrix is not square") }
    if (try.chol) {
      L <- try(.Cholwrap(m),silent=TRUE)
      if (inherits(L,"try-error")) {
        mreg <- m *(1-1e-08)
        diag(mreg) <- diag(mreg) + 1e-8 * diag(m) ## allows covMatrix; diag(m) has 1's if m is correlation matrix
        L <- try(.Cholwrap(mreg),silent=TRUE)
      } 
      if ( ! inherits(L,"try-error") ) {
        type <- "cholL_LLt"  
      } ## else type remains NULL      
    }
    if ( is.null(type) ) { ## hence if none of the chol algos (nor symSVD input) has been used 
      ## slower by more stable. Not triangular but should not be a pb
      if (try.eigen) { ## R's eigen()
        LDL <- try(eigen(m,symmetric=TRUE),silent=TRUE) ## may _hang_ in R2.15.2 on nearly-I matrices
        if ( ! inherits(LDL,"try-error")) d <- LDL$values 
      }
      if ( (! try.eigen) || inherits(LDL,"try-error") || any(d < -1e-08)) {
        if (.spaMM.data$options$USEEIGEN) { ## see package irlba for SVD of sparse matrices
          symSVD <- sym_eigen(m) ## such that v= t(u) without any sign issue  
          u <- symSVD$u
          d <- symSVD$d  
          type <- "symsvd" ## RcppEigen's SVD: ## "the SVD implementation of Eigen (...) is not a particularly fast SVD method." (RcppEigen vignette)          
        } else { ## buggy code for testing !  $USEEIGEN
          SVD <- try(svd(m)) # R's svd()
          if (inherits(SVD,"try-error")) {
            print("spaMM retries SVD following 'svd' failure.") 
            ## numerically different but otherwise equivalent computation
            m <- diag(rep(1-SVDfix,ncol(m))) + m*SVDfix 
            SVD <- try(svd(m)) 
            if (! inherits(SVD,"try-error") ) SVD$d <- 1+ (SVD$d-1)/SVDfix # valid for covMatrix, but maybe not optimal
          } 
          if (inherits(SVD,"try-error")) {
            ## typically m is I + a few large elements
            ## Most informative post: http://r.789695.n4.nabble.com/Observations-on-SVD-linpack-errors-and-a-workaround-td837282.html
            print("Singular value decomposition failed.") 
            print(" See documentation of the 'SVDfix' argument of 'designL.from.Corr'")
            print("   for ways to handle this.")
            return(try(stop(),silent=TRUE)) ## passes control to calling function
          } else {
            d <- SVD$d
            ## must be valid for sym (semi) PD matrices using U, V being eigenvectors of m %*% t(m)
            ## symm matrix => $u, $v match left and right eigenvectors of original Corr matrix
            ## FR->FR but they can be of opposite sign with negative $d...
            ## => bug in this code not normally used
            u <- SVD$u
            type <- "svd"          
          }         
        } ## svd by R
      } else { ## we have an LDL decomp
        u <- LDL$vectors
        d <- LDL$values
        type <- "eigen"
      }
      if (any(d< -1e-08)) {
        ## could occur for two reasons: wrong input matrix; or problem with R's svd which $d are the eigenvalues up to the sign, 
        ##   and thus $d can be <0 for eigenvalues >0 (very rare, observed in a 2x2 matrix) 
        message("correlation matrix has suspiciously large negative eigenvalue(s).")
        return(try(stop(),silent=TRUE)) ## passes control to calling function to print correlation params
      } else { ## we have a not-too-suspect decomp
        # d[d< threshold]<- threshold ## wrong for corrmats, would be OK for d2hdv2 computation which uses this function 
        if (any(d<threshold)) d <- threshold + (1-threshold) * d ## 17/05/2014 ## maybe not optimal for covMatrix
        L <- .ZWZt(u,sqrt(d))
      }
    } 
    colnames(L) <- rownames(L) <- colnames(m) ## for checks in HLfit ## currently typically missing from symSVD  
    attr(L,"corr.model") <- attr(m,"corr.model")
  } else {
    u <- symSVD$u ## local copy needed for processing attributes at the end of the function
    d <- symSVD$d ## of corr matrix !
    L <- .ZWZt(u,sqrt(d))  
    type <- "symsvd"      
    attr(L,"corr.model") <- symSVD$corr.model
  }
  attr(L,"type") <- type
  ## add the 'sanitized' matrix decomp as attribute of the L matrix
  if (type != "cholL_LLt") {
    decomp <- list(u=u,d=d)
    if ( ! is.null(symSVD)) decomp$adjd <- symSVD$adjd ## useful for SEM CAR; otherwise may be NULL
    attr(L,type) <- decomp
  }
  return(L)
} 

make_scaled_dist <- function(uniqueGeo,uniqueGeo2=NULL,distMatrix,rho,rho.mapping=seq_len(length(rho)),
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
      if (inherits(scaled.dist,"ff_matrix")) {
        scaled.dist[] <- proxy::dist(x=uniqueScal,y=uniqueScal2, method=dist.method) 
      } else scaled.dist <- proxy::dist(x=uniqueScal,y=uniqueScal2, method=dist.method) 
    } else { ## not Euclidean
      if (inherits(scaled.dist,"ff_matrix")) {
        scaled.dist[] <- rho * proxy::dist(uniqueGeo,y=uniqueGeo2,method=dist.method) 
      } else scaled.dist <- rho * proxy::dist(uniqueGeo,y=uniqueGeo2,method=dist.method)  
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

getDistMat <- function(object,scaled=FALSE) {
  if (! is.null(dist_info <- attr(object,"dist_info"))) {
    if (is.null(dist_info$distMatrix)) { ## we reconstruct it
      dist_info$distMatrix <- eval(dist_info$distcall,
                                   list(uniqueGeo=attr(object,"info.uniqueGeo"),dist.method=dist_info$dist.method))
      dist_info$distcall <- NULL 
    }
    if ( ! scaled)  {
      dist_info$rho <- 1 
      dist_info$`rho.mapping` <- NULL 
    }
    return(do.call(make_scaled_dist,dist_info))
  } else {
    message("no Matern-correlated random effects")
    return(NULL)
  }
}
