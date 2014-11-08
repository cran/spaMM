HL_process.args <- function (...) { ## from gsl package
  a <- list(...)
  attr <- attributes(a[[which.max(unlist(lapply(a, length)))]])
  a <- lapply(a, as.vector)
  out <- do.call("rbind", a)
  out <- list(`1`=out[1,],`2`=out[2,]) ## == out <- split(out, row(out)) but profiling shows the latter is slow
  names(out) <- paste("arg", 1:length(a), sep = "")
  return(c(out, attr = list(attr)))
}

HL_strictify <- function (val, status) { ## from gsl package
    val[status < 0] <- NaN ##FR inverted the test // gsl package !!
    return(val)
}

bessel_lnKnu <- function (nu, x, give = FALSE, strict = TRUE) { ## from bessel_lnKnu in gsl package
    jj <- HL_process.args(nu, x)
    nu.vec <- jj$arg1
    x.vec <- jj$arg2
    attr <- jj$attr
    jj <- .C("bessel_lnKnu_e", as.double(nu.vec), as.double(x.vec), 
        as.integer(length(x.vec)), val = as.double(x.vec), err = as.double(x.vec), 
        status = as.integer(0 * x.vec), PACKAGE = "spaMM")
    val <- jj$val
    err <- jj$err
    status <- jj$status
    attributes(val) <- attr
    attributes(err) <- attr
    attributes(status) <- attr
    if (strict) {
        val <- HL_strictify(val, status)
    }
    if (give) {
        return(list(val = val, err = err, status = status))
    }
    else {
        return(val)
    }
}


`Matern.corr`<-function (d, rho=1, smoothness, nu=smoothness, Nugget=0L) { ## rho is alpha in fields
  ## ideally (but not necess) on a 'dist' so the diagonal is not  manipulated 
    if (any(d < 0)) 
        stop("distance argument must be nonnegative")
    dscal <- d * rho
    dscal[d == 0L] <- 1e-10 ## avoids errors on distance =0; but the resulting corrvals can be visibly < 1 for small nu ## FR->FR make value dependent on rho, nu ?
    logcon <- (nu - 1)*log(2)+ lgamma(nu) 
    corrvals <- - logcon + nu*log(dscal)+ bessel_lnKnu(x=dscal, nu=nu) ## 
##    corrvals <- - logcon + nu*log(dscal)+ log(besselK(x=dscal, nu=nu)) ## function from package gsl
    corrvals <- exp(corrvals) 
    corrvals[d != 0L] <- (1-Nugget)* corrvals[d != 0L]
    corrvals[d == 0L] <- 1 ## 
    corrvals[corrvals < 1e-16] <- 0L ## an attempt to deal with problem in chol/ldl/svd which don't like 'nearly-identity' matrices
    return(corrvals)
}

## ess <- function(nu,d) {exp(-(d/(2*sqrt(nu)))^2)} ...


Matern.local <- Matern.corr   ## back.compatibility code

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

`designL.from.Corr` <- function(m,try.chol=TRUE,try.eigen=FALSE,threshold=1e-06,debug=FALSE,SVDfix=1/10) {
  type <- NULL
  if (try.chol) {
    if (.spaMM.data$options$USEEIGEN) {
      if (inherits(m,"blockDiag")) { 
        trychol <- RcppChol.blockDiag(m) ## cf pb RcppEigen / generic
      } else trychol <- RcppChol(m)
      if (trychol$Status==TRUE) { ## if chol computation successful
        L <- trychol$L
        type <- "chol"  
      } else {
        mreg <- m *(1-1e-8)
        diag(mreg) <- diag(mreg) + 1e-8 
        trychol <- RcppChol(mreg)
        if (trychol$Status==TRUE) {
          L <- trychol$L
          type <- "chol"  
        } ## else type remains NULL      
      } 
    } else { ## chol de R
      cholR <- try(chol(m),silent=TRUE) ## dim -> unique geo coordinates
      if (class(cholR)=="try-error") { ## attempt at regularization 050213
        mreg <- m *(1-1e-8)
        diag(mreg) <- diag(mreg) + 1e-8 
        cholR <- try(chol(mreg),silent=TRUE) 
        if (class(cholR) != "try-error") type <- "chol"
      } else type <- "chol"
      if (!is.null(type)) L<-t(cholR) ## such that L %*% t(L) = the correlation matrix in both cases        
    }    
  }
  if ( is.null(type) ) { ## hence if none of the chol algos has been used 
    ## slower by more stable. Not triangular but should not be a pb
    if (try.eigen) {
      LDL <- try(eigen(m,symmetric=TRUE),silent=TRUE) ## may _hang_ in R2.15.2 on nearly-I matrices
      if (class(LDL) != "try-error") d <- LDL$values 
    }
    if ( (! try.eigen) || class(LDL)=="try-error" || any(d < -1e-08)) {
      if (.spaMM.data$options$USEEIGEN) {
        symSVD <- selfAdjointSolverCpp(m) ## such that v= t(u)without any sign issue  
        u <- symSVD$u
        d <- symSVD$d  
        type <- "symsvd"          
      } else {
        SVD <- try(svd(m)) ## "the SVD implementation of Eigen (...) is not a particularly fast SVD method." (RcppEigen vignette)
        if (class(SVD)=="try-error") {
          print("spaMM retries SVD following 'svd' failure.") 
          ## numerically different but otherwise equivalent computation
          m <- diag(rep(1-SVDfix,ncol(m))) + m*SVDfix 
          SVD <- try(svd(m)) 
          if (class(SVD)!="try-error") SVD$d <- 1+ (SVD$d-1)/SVDfix 
        } 
        if (class(SVD)=="try-error") {
          ## typically m is I + a few large elements
          ## Most iformative post: http://r.789695.n4.nabble.com/Observations-on-SVD-linpack-errors-and-a-workaround-td837282.html
          ####################### aide sur le .Rdata :  zut <- as.list(attr(resu,"condition")$call) est un HLCor call
          print("Singular value decomposition failed.") 
          print(" See documentation of the 'SVDfix' argument of 'designL.from.Corr'")
          print("   for ways to handle this.")
          return(try(stop(),silent=T)) ## passes control to calling function
        } else {
          d <- SVD$d
          ## must be valid for sym (semi) PD matrices using U, V being eigenvectors of m %*% t(m)
          ## symm matrix => $u, $v match left and right eigenvectors of original Corr matrix
          ## FR->FR but they can be of opposite sign with negative $d...
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
      mess <- pastefrom("correlation matrix has suspiciously large negative eigenvalue(s).")
      print(mess)
      return(try(stop(),silent=T)) ## passes control to calling function
    } else { ## we have a not-too-suspect decomp
      # d[d< threshold]<- threshold ## no longer a corrMat
      if (any(d<threshold)) d <- threshold + (1-threshold) * d ## 17/05/2014
      ## FR->FR but we also use this function once on d2hdv2 in HLfit...
      L <- ZWZt(u,sqrt(d))
    }
  } 
  colnames(L) <- colnames(m)
  rownames(L) <- colnames(m) ## for checks in HLfit
  attr(L,"type") <- type
  ## add the 'sanitized' matrix decomp as attribute of the L matrix
  if (type != "chol") attr(L,type) <- list(u=u,d=d)
  return(L)
} 

## we want to compute (1) [2e ligne de ChanK97 eq 11] the conditional covariance Cov - Cov inv(Cov+I) Cov  (we return a Cholesky factor of it)
## and (2) the conditional partial regression coeffs for Lv: Cov inv(Cov+I)
##  It turns out that the two are identical : 
## If Corr= LDL', cov= lam LDL', we want  L [lam D - lam^2 D^2/(lam D+I)] L'
## =  L [lam D/(lam D +I)] L' 
CondNormfn <- function(LMatrix,lambda) {
  type <- attr(LMatrix,"type")
  if (is.null(type)) {stop("'type' attribute needed for 'LMatrix' argument in 'CondNormfn'.")}
  if (type=="chol") {stop("LMatrix argument in 'CondNormfn' is inadequate because its 'type' attribute is 'chol'.")}
  decomp <- attr(LMatrix,type)
  diago <- decomp$d/(decomp$d+1/lambda)
  sqrtCondCovLv <- sweep(decomp$u,2,sqrt(diago),`*`); ## so that cond Corr = this.t(this)
  condLvReg <- tcrossprodCpp(sqrtCondCovLv) ## conditional regr = cond Corr
  # FR->FR un probleme est que la repres sqrtCondCovLv est loin d'être "numerically unique". 
  # Donc même si on a des distrib equivalentes pour differents sqrtCondCovLv 
  # (en particulier des condLvReg equivalents)
  # on va avoir sqrtCondCovLv %*% rnorm nettement divergents sous linux vs Windows 
  return(list(sqrtCondCovLv=sqrtCondCovLv,condLvReg=condLvReg))
} 

`make.scaled.dist` <- function(distMatrix,uniqueGeo,rho,rho.mapping=seq_len(length(rho))) {
  if ( missing(distMatrix) ) { ## 
     if ( missing(uniqueGeo) ) {
       mess <- pastefrom("missing(distMatrix) && missing(uniqueGeo).",prefix="(!) From ")
       stop(mess)
     } else {
        if (length(rho)==1) {
           uniqueScal<-uniqueGeo * rho 
        } else if (ncol(uniqueGeo)==length(rho.mapping)) {
           uniqueScal<-t(t(uniqueGeo) * rho[rho.mapping]) ## valid for vectorial rho...
        } else {
          mess <- pastefrom("invalid length(rho[rho.mapping]).",prefix="(!) From ")
          print(mess)
          mess  <- paste("Length should be either 1 or",ncol(uniqueGeo))
          stop(mess)
        }
     }
     scaled.dist <- proxy::dist(uniqueScal) ## one can do a lot with a dist object !
     ## here scaled.dist is always a dist object: if uniqueScal was a single row, it is dist(0) 
  } else if (length(rho)==1) {
     scaled.dist <- rho * distMatrix
     ## here inconsistent behavior: scaled.dist is a dist object except if distMatrix was dist(0), scaled.dist is now numeric(0); we standardise it:
     if ( identical(scaled.dist, numeric(0))) scaled.dist <- dist(0)
  } else {
    mess <- pastefrom("input 'distMatrix' but length(rho)!=1.",prefix="(!) From ")
    stop(mess)
  }
  return(scaled.dist)
}


