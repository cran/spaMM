designL.from.Corr <-
function(m,try.chol=TRUE,try.eigen=FALSE,threshold=1e-06,debug=FALSE,SVDfix=1/10) {
  type <- NULL
  if (try.chol) {
    if (.SpaMM$USEEIGEN) {
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
      if (.SpaMM$USEEIGEN) {
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
