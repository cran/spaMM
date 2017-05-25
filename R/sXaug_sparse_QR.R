# 'constructor' for sXaug_Eigen_sparse_QR (unpivoted QR factorization) object
# from Xaug which already has a *scaled* ZAL 
def_sXaug_EigenSparse_QR_scaled <- function(Xaug,weight_X,w.ranef,H_global_scale) {
  n_u_h <- length(w.ranef)
  Xrows <- n_u_h+seq(length(weight_X))
  #Xaug[Xrows,] <- Diagonal(x=weight_X) %*% Xaug[Xrows,] ## apparemment sweep vraiment pas efficace sur Matrix
  which_i_affected_rows <- Xaug@i>(n_u_h-1L)
  Xaug@x[which_i_affected_rows] <- Xaug@x[which_i_affected_rows]*weight_X[Xaug@i[which_i_affected_rows]-n_u_h+1L] 
  resu <- structure(Xaug,
                    sXaug_EigenSparse_QR_scaled=create_sXaug_EigenSparse_QR_scaled(),
                    get_from="get_from_MME.sXaug_EigenSparse_QR_scaled",
                    w.ranef=w.ranef,
                    n_u_h=n_u_h,
                    pforpv=ncol(Xaug)-n_u_h,
                    H_global_scale=H_global_scale
  )
  ## cannot modify the class attribute... => immediate lumsy code below...
  return( resu ) 
}

# single 'internal method' (with <<- ) for sXaug_Eigen_sparse_QR objects
# if QR not previously available, szAug=NULL: constructs the QR factorization, return value depends on 'which'
# if QR not previously available, szAug ! NULL:   ..........................   returns the solution
# if which="solve_d2hdv2", B=NULL: returns solve(d2hdv2)
# if which="solve_d2hdv2", B ! NULL: returns solve(d2hdv2,B)
create_sXaug_EigenSparse_QR_scaled <- function() {
  BLOB <- list()
  locfn <- function(sXaug,which="",szAug=NULL,B=NULL) {
    if (is.null(BLOB$R_scaled)) {
      blob <- lmwith_sparse_QRp(sXaug,szAug,
                               returntQ=TRUE, ## actually FALSE, cf lmwith_sparse_QRp code
                               returnR=TRUE,pivot=FALSE) 
      ## ... using RcppEigen; szAug may be NULL
      seq_ncol <- seq(ncol(sXaug))
      BLOB$R_scaled <<- blob$R #[seq_ncol,]
      if ( ! is.null(szAug)) return(blob$coef)   
    } 
    if ( is.null(BLOB$t_Q_scaled) && 
         (#! is.null(szAug) || 
           which %in% c("t_Q_scaled","hatval","hatval_Z")) ) {
      BLOB$t_Q_scaled <<- Matrix::drop0(Matrix::solve(t(BLOB$R_scaled),t(sXaug)))
    }
    if ( ! is.null(szAug)) {
      if (is.null(BLOB$t_Q_scaled)) {
        ## many(solve(augmented linear equations) in GLMM, avoid t_Q_computation there)
        return(Matrix::solve(crossprod(BLOB$R_scaled), crossprod(sXaug, szAug)))
      } else {
        return(Matrix::solve(BLOB$R_scaled,BLOB$t_Q_scaled %*% szAug))      }
    }
    # ELSE
    if ( is.null(BLOB$d2hdv2) && 
         (which =="d2hdv2" || (which=="solve_d2hdv2" && is.null(BLOB$inv_d2hdv2) && !is.null(B)) ) 
    ) {
      w.ranef <- attr(sXaug,"w.ranef")
      seq_n_u_h <- seq_len(attr(sXaug,"n_u_h"))
      R_Md2hdv2 <- Matrix_times_Dvec(BLOB$R_scaled[seq_n_u_h, seq_n_u_h], sqrt(w.ranef)) ## : such that R'R gives -d2hdv2
      BLOB$d2hdv2 <<- - Matrix::crossprod(R_Md2hdv2)
    }
    if (which=="t_Q_scaled") { 
      return(BLOB$t_Q_scaled)
    } else if (which=="hatval") {
      if (is.null(BLOB$hatval_ZX) ) {
        tmp <- BLOB$t_Q_scaled
        tmp@x <- tmp@x^2
        BLOB$hatval_ZX <<- colSums(tmp)
      }
      return(BLOB$hatval_ZX)
    } else if (which=="hatval_Z") {
      if (is.null(BLOB$hatval_Z_)) {
        n_u_h <- attr(sXaug,"n_u_h")
        tmp <- BLOB$t_Q_scaled[ seq_len(n_u_h) ,]
        tmp@x <- tmp@x^2
        tmp <- colSums(tmp)
        phipos <- n_u_h+seq_len(nrow(sXaug)-n_u_h)
        BLOB$hatval_Z_ <<-  list(lev_lambda=tmp[-phipos],lev_phi=tmp[phipos])      
      }
      return(BLOB$hatval_Z_)
    } else if (which=="R_scaled") {
      return(BLOB$R_scaled)
    } else if (which=="d2hdv2") {
      return(BLOB$d2hdv2)
    } else if (which=="diag_RtR") {
      if (is.null(BLOB$diag_RtR)) {
        tmp <- BLOB$R_scaled
        tmp@x <- tmp@x^2
        BLOB$diag_RtR <<- colSums(tmp)
      }
      return(BLOB$diag_RtR)
    } else if (which=="sortPerm") {
      return(BLOB$sortPerm) ## NULL !
    } else if (which=="solve_d2hdv2") {
      if (is.null(B)) {
        if (is.null(BLOB$inv_d2hdv2)) {
          w.ranef <- attr(sXaug,"w.ranef")
          seq_n_u_h <- seq_len(attr(sXaug,"n_u_h"))
          R_Md2hdv2 <- Matrix_times_Dvec(BLOB$R_scaled[seq_n_u_h, seq_n_u_h], sqrt(w.ranef)) ## : such that R'R gives -d2hdv2
          BLOB$inv_d2hdv2 <<- - Matrix::chol2inv(R_Md2hdv2) ## valid bc chol has no pivoting
        }
        if (is.null(B)) {
          return(BLOB$inv_d2hdv2)
        } else return(BLOB$inv_d2hdv2 %*% B) ## FR->FR could add speific code for Diagonal B ?
      } else { ## solve (Matrix,vector)
        if (is.null(BLOB$inv_d2hdv2)) {
          return( Matrix::solve(BLOB$d2hdv2, B)) 
        } else return(BLOB$inv_d2hdv2 %*% B)
      }
    } else if (which=="beta_cov") { 
      beta_cov <- Matrix::chol2inv(BLOB$R_scaled) ## valid bc no pivoting
      beta_pos <- attr(sXaug,"n_u_h")+seq_len(attr(sXaug,"pforpv"))
      beta_cov <- beta_cov[beta_pos,beta_pos,drop=FALSE] * attr(sXaug,"H_global_scale")
      colnames(beta_cov) <- colnames(sXaug)[beta_pos]
      return(beta_cov)
    } else if (which=="beta_v_cov_from_wAugX") { ## using a weighted Henderson's augmented design matrix, not a true sXaug  
      beta_v_cov <- as.matrix(Matrix::chol2inv(BLOB$R_scaled))
      # this tends to be dense bc v_h estimates covary (example: wafers)
      # otherwise, dropO(,tol=...), + some fix in summary.HLfit for matrix[] <- Matrix assignment, would be useful.  
      return(beta_v_cov)
    }   
    stop("invalid 'which' value.")
  }
  environment(locfn) <- list2env(list(BLOB=list()),
                                 parent=environment(def_sXaug_EigenSparse_QR_scaled))
  return(locfn)
} 

get_from_MME.sXaug_EigenSparse_QR_scaled <- function(sXaug,which="",szAug=NULL,B=NULL,
                                           damping, LMrhs, ...) {
  R_scaled <- attr(sXaug,"sXaug_EigenSparse_QR_scaled")(sXaug,which="R_scaled")
  resu <- switch(which,
                 "logdet_R_scaled_v" = {
                   sum(log(abs(diag(R_scaled)[ seq_len(attr(sXaug,"n_u_h")) ]))) 
                 },
                 "logdet_R_scaled_b_v" = {
                   sum(log(abs(diag(R_scaled)))) ## 
                 },
                 "LevMar_step" = {
                   diag_RtR <- attr(sXaug,"sXaug_EigenSparse_QR_scaled")(sXaug,which="diag_RtR")
                   dampDpD <- damping*diag_RtR ## NocedalW p. 266
                   if (FALSE) {
                     ## version QR inspree de More (mais sans les Givens bien choisis)
                     Raug <- rbind(R_scaled, Diagonal(x=sqrt(dampDpD)))
                     stop("need asparse version of the following")
                     resu <- LevMar_cpp(Raug,LMrhs=LMrhs)
                     resu$dampDpD <- dampDpD
                   } else {
                     ## version chol assez directe (a traduire en c++?)
                     AtAdDpD <- Matrix::crossprod(R_scaled)
                     nc <- ncol(AtAdDpD)
                     diagPos <- seq.int(1L,nc^2,nc+1L)
                     AtAdDpD[diagPos] <- AtAdDpD[diagPos] + dampDpD 
                     blobR <- Matrix::Cholesky(AtAdDpD,LDL=FALSE)
                     RRsP <- sort.list(as(blobR,"pMatrix")@perm) 
                     RRR <- t(as(blobR,"sparseMatrix"))
                     resu <- list(dVscaled_beta= Matrix::chol2inv(RRR)[RRsP,RRsP] %*% LMrhs,
                                  rhs=LMrhs, dampDpD = dampDpD)
                   }
                   return(resu)
                 },
                 "logdet_sqrt_d2hdv2" = {
                   w.ranef <- attr(sXaug,"w.ranef")
                   sum(log(sqrt(w.ranef) * abs(diag(R_scaled)[ seq_len(attr(sXaug,"n_u_h")) ])))
                 },
                 # "logdet_sqrt_d2hdBV2" = {
                 #   w.ranef <- attr(sXaug,"w.ranef")
                 #   sum(log(abs(diag(R_scaled)))) + sum(log(w.ranef))/2 - (attr(sXaug,"pforpv"))*log(H_global_scale)/2 
                 # },
                 "logdet_sqrt_d2hdbeta2" = {
                   H_global_scale <- attr(sXaug,"H_global_scale")
                   sum(log(abs(diag(R_scaled)[ - seq_len(attr(sXaug,"n_u_h")) ]))) - attr(sXaug,"pforpv")*log(H_global_scale)/2 
                 },
                 ## all other cases:
                 attr(sXaug,"sXaug_EigenSparse_QR_scaled")(sXaug,which=which,szAug=szAug,B=B)
  )
  return(resu)
}
