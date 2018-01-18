# 'constructor' for sXaug_Eigen_QR (unpivoted QR factorization) object
# from Xaug which already has a *scaled* ZAL 
def_sXaug_EigenDense_QRP_Chol_scaled <- function(Xaug,weight_X,w.ranef,H_global_scale) {
  n_u_h <- length(w.ranef)
  Xrows <- n_u_h+seq(length(weight_X))
  Xaug[Xrows,] <- .Dvec_times_matrix(weight_X, Xaug[Xrows,,drop=FALSE]) # sweep(TT,MARGIN=1,sqrt.ww,`*`) # rWW %*%TT
  resu <- structure(Xaug,
                    get_from="get_from_MME.sXaug_EigenDense_QRP_Chol_scaled",
                    BLOB=list2env(list(), parent=environment(.sXaug_EigenDense_QRP_Chol_scaled)),
                    w.ranef=w.ranef,
                    n_u_h=n_u_h,
                    pforpv=ncol(Xaug)-n_u_h,
                    H_global_scale=H_global_scale
  )
  class(resu) <- c(class(resu),"sXaug_EigenDense_QRP_Chol_scaled")
  return( resu ) 
}


.sXaug_EigenDense_QRP_Chol_scaled <- function(sXaug,which="",szAug=NULL,B=NULL) {
  BLOB <- attr(sXaug,"BLOB") ## an environment
  if (is.null(BLOB$R_scaled)) {
    lmwithqrp <- .lmwithQRP(sXaug,yy=szAug,returntQ=FALSE,returnR=TRUE) ## using RcppEigen; szAug may be NULL
    for (st in names(lmwithqrp)) BLOB[[st]] <- lmwithqrp[[st]] ## including BLOB$R_scaled...
    BLOB$perm <- BLOB$perm +1L
    n_u_h <- attr(sXaug,"n_u_h")
    seq_n_u_h <- seq(n_u_h)
    BLOB$sortPerm <- sort.list(BLOB$perm) 
    if ( ! is.null(szAug)) return(BLOB$coef)   
  } 
  if ( ! is.null(szAug)) {
    rhs <- .crossprodCpp(sXaug, szAug)
    rhs <- rhs[BLOB$perm,,drop=FALSE]
    rhs <- backsolve(BLOB$R_scaled, forwardsolve(BLOB$R_scaled, rhs, upper.tri = TRUE, transpose = TRUE))
    return(rhs[BLOB$sortPerm,,drop=FALSE])
  }
  # ELSE
  if ( is.null(BLOB$R_R_v) && which %in% c("absdiag_R_v","d2hdv2","solve_d2hdv2","hatval_Z","R_scaled_v_h_blob")) {
    seq_n_u_h <- seq_len(attr(sXaug,"n_u_h"))
    # remove $u_h_cols_on_left code in version 2.1.61 since there is a bug in it: further code may require 
    #  no nontrivial perm_R_v, so we remove perm_R_v it from further code (as in Matrix_QRP_CHM_scaled version)
    wd2hdv2w <- .crossprodCpp(BLOB$R_scaled[,BLOB$sortPerm[ seq_n_u_h ]], NULL)
    BLOB$R_R_v <- .Rcpp_chol_R(wd2hdv2w)$R
  }
  # 
  if (which=="t_Q_scaled") { 
    if ( is.null(BLOB$t_Q_scaled) ) BLOB$t_Q_scaled <- backsolve(BLOB$R_scaled,t(sXaug[,BLOB$perm]),transpose=TRUE)
    return(BLOB$t_Q_scaled)
  } else if (which=="hatval") {
    if ( is.null(BLOB$t_Q_scaled) ) BLOB$t_Q_scaled <- backsolve(BLOB$R_scaled,t(sXaug[,BLOB$perm]),transpose=TRUE)
    return(colSums(BLOB$t_Q_scaled^2))
  } else if (which=="hatval_Z") { ## Pdiag
    if (is.null(BLOB$hatval_Z_)) {
      ## t(sXaug) is scaled such that the left block is an identity matrix, so we can work 
      ## on two separate blocks if the Cholesky is not permuted. Then 
      lev_lambda <- backsolve(BLOB$R_R_v,diag(ncol(BLOB$R_R_v)),transpose=TRUE)
      lev_lambda <- lev_lambda^2
      lev_lambda <- colSums(lev_lambda)
      n_u_h <- attr(sXaug,"n_u_h")
      phipos <- (n_u_h+1):nrow(sXaug)
      lev_phi <- backsolve(BLOB$R_R_v, t(sXaug[phipos, seq_len(n_u_h) ]),transpose=TRUE) ## fixme the t() may still be costly
      lev_phi <- lev_phi^2
      lev_phi <- colSums(lev_phi)
      BLOB$hatval_Z_ <-  list(lev_lambda=lev_lambda,lev_phi=lev_phi)
    }
    return(BLOB$hatval_Z_)
  } else if (which=="R_scaled") { ## for logdet_r22
    return(BLOB$R_scaled)
  } else if (which=="R_scaled_blob") {
    if (is.null(BLOB$R_scaled_blob)) {
      X <- BLOB$R_scaled[,BLOB$sortPerm]
      BLOB$R_scaled_blob <- list(X=X, diag_pRtRp = colSums(X^2))
    }
    return(BLOB$R_scaled_blob)
  } else if (which=="R_scaled_v_h_blob") {
    if (is.null(BLOB$R_scaled_v_h_blob)) {
      diag_pRtRp_scaled_v_h <- colSums(BLOB$R_R_v^2)
      BLOB$R_scaled_v_h_blob <- list(R_scaled_v_h=BLOB$R_R_v, diag_pRtRp_scaled_v_h=diag_pRtRp_scaled_v_h)
    }
    return(BLOB$R_scaled_v_h_blob)
  } else if (which=="d2hdv2") {
    stop("d2hdv2 requested")
    # don't forget that the factored matrix is not the augmented design matrix ! hence w.ranef needed here
    w.ranef <- attr(sXaug,"w.ranef")
    w_R_R_v <- BLOB$R_R_v %*% diag(sqrt(w.ranef))
    BLOB$d2hdv2 <- - .crossprodCpp(w_R_R_v,yy=NULL)
  } else if (which=="Mg_invH_g") {
    return(sum(solve(BLOB$R_R_v,B)^2))
  } else if (which=="Mg_solve_g") {
    if (is.null(BLOB$invsqrtwranef)) BLOB$invsqrtwranef <- 1/sqrt(attr(sXaug,"w.ranef")) 
    seq_n_u_h <- seq_len(attr(sXaug,"n_u_h"))
    rhs <- B
    rhs[seq_n_u_h] <- BLOB$invsqrtwranef * rhs[seq_n_u_h]
    rhs <- forwardsolve(BLOB$R_scaled, rhs[BLOB$perm], upper.tri = TRUE, transpose = TRUE)
    return(sum(rhs^2))
  } else if (which=="solve_d2hdv2") { ## R'R[seq_n_u_h, seq_n_u_h] gives -d2hd_scaled_v2.
    if ( is.null(B) ) {
      if ( is.null(BLOB$inv_d2hdv2)) {
        if (is.null(BLOB$invsqrtwranef)) BLOB$invsqrtwranef <- 1/sqrt(attr(sXaug,"w.ranef")) 
        inv_d2hdv2 <- chol2inv(BLOB$R_R_v)
        inv_d2hdv2 <- .m_Matrix_times_Dvec(inv_d2hdv2, BLOB$invsqrtwranef) ## _m_atrix sweep 2L
        BLOB$inv_d2hdv2 <- - .Dvec_times_matrix(BLOB$invsqrtwranef,inv_d2hdv2)
      }
      return(BLOB$inv_d2hdv2)
    } else { ## solve (matrix,y)
      if (is.null(BLOB$inv_d2hdv2)) {
        if (is.null(BLOB$invsqrtwranef)) BLOB$invsqrtwranef <- 1/sqrt(attr(sXaug,"w.ranef")) 
        if (is.matrix(B)) {
          rhs <- .Dvec_times_matrix(BLOB$invsqrtwranef,B)
        } else rhs <- BLOB$invsqrtwranef * B
        rhs <- backsolve(BLOB$R_R_v, forwardsolve(BLOB$R_R_v, rhs, upper.tri = TRUE, transpose = TRUE))
        if (is.matrix(rhs)) {
          rhs <- .Dvec_times_matrix(BLOB$invsqrtwranef,rhs)
        } else rhs <- BLOB$invsqrtwranef * rhs
        return( - rhs)
      } else {
        return(BLOB$inv_d2hdv2 %*% B) ## inv_d2hdv2 tends to be dense
      }
    }
  } else if (which %in% c("logdet_R_scaled_b_v")) {
    if (is.null(BLOB$logdet_R_scaled_b_v)) BLOB$logdet_R_scaled_b_v <- sum(log(abs(diag(BLOB$R_scaled))))
    return(BLOB$logdet_R_scaled_b_v)
  } else if (which %in% c("logdet_R_scaled_v")) {
    if (is.null(BLOB$logdet_R_scaled_v)) BLOB$logdet_R_scaled_v <- sum(log(abs(diag(BLOB$R_R_v))))
    return(BLOB$logdet_R_scaled_v)
  } else if (which=="beta_cov") { 
    beta_cov <- chol2inv(BLOB$R_scaled) ## actually augmented beta_v_cov as the following shows
    beta_pos <- attr(sXaug,"n_u_h")+seq_len(attr(sXaug,"pforpv"))
    sP_beta_pos <- BLOB$sortPerm[beta_pos]
    beta_cov <- beta_cov[sP_beta_pos,sP_beta_pos,drop=FALSE] * attr(sXaug,"H_global_scale")
    colnames(beta_cov) <- colnames(sXaug)[beta_pos]
    return(beta_cov)
  } else if (which=="beta_v_cov_from_wAugX") { ## using a weighted Henderson's augmented design matrix, not a true sXaug  
    beta_v_cov <- chol2inv(BLOB$R_scaled)[BLOB$sortPerm,BLOB$sortPerm,drop=FALSE]
    # this tends to be dense bc v_h estimates covary (example: wafers)
    # otherwise, dropO(,tol=...), + some fix in summary.HLfit for matrix[] <- Matrix assignment, would be useful.  
    return(beta_v_cov)
  } else if (which %in% c("absdiag_R_v")) { ## used for logdet_R_scaled_v
    if (is.null(BLOB$absdiag_R_v)) BLOB$absdiag_R_v <- abs(diag(BLOB$R_R_v)) 
    return(BLOB$absdiag_R_v)
  } else if (which=="sortPerm") { ## extracted for LevMar
    return(BLOB$sortPerm)
  }
  stop("invalid 'which' value.")
} 

# trace("get_from_MME.sXaug_EigenDense_QRP_Chol_scaled",print=FALSE, tracer=quote(print(which)),exit=quote(str(resu)))
get_from_MME.sXaug_EigenDense_QRP_Chol_scaled <- function(sXaug,which="",szAug=NULL,B=NULL,
                                    damping, LMrhs, ...) {
  resu <- switch(which,
                 "logdet_sqrt_d2hdv2" = {
                   w.ranef <- attr(sXaug,"w.ranef")
                   absdiag_R_v <- .sXaug_EigenDense_QRP_Chol_scaled(sXaug,which="absdiag_R_v")
                   sum(log(sqrt(w.ranef) * absdiag_R_v))
                 },
                 "logdet_r22" = {
                   H_global_scale <- attr(sXaug,"H_global_scale")
                   R_scaled <- .sXaug_EigenDense_QRP_Chol_scaled(sXaug,which="R_scaled")
                   absdiag_R_v <- .sXaug_EigenDense_QRP_Chol_scaled(sXaug,which="absdiag_R_v")
                   sum(log(abs(diag(R_scaled)))) - sum(log(absdiag_R_v)) - attr(sXaug,"pforpv")*log(H_global_scale)/2 
                 },
                 "LevMar_step" = {
                   R_scaled_blob <- .sXaug_EigenDense_QRP_Chol_scaled(sXaug,which="R_scaled_blob")
                   dampDpD <- damping*R_scaled_blob$diag_pRtRp ## NocedalW p. 266
                   list(dVscaled_beta = .damping_to_solve(X=R_scaled_blob$X,dampDpD=dampDpD,rhs=LMrhs), 
                        dampDpD = dampDpD) 
                 },
                 "LevMar_step_v_h" = {
                   ## FR->FR probably not the most elegant implementation 
                   R_scaled_v_h_blob <- .sXaug_EigenDense_QRP_Chol_scaled(sXaug,which="R_scaled_v_h_blob")
                   dampDpD <- damping*R_scaled_v_h_blob$diag_pRtRp_scaled_v_h ## NocedalW p. 266
                   list(dVscaled = .damping_to_solve(X=R_scaled_v_h_blob$R_scaled_v_h, dampDpD=dampDpD, rhs=LMrhs), 
                        dampDpD = dampDpD) 
                 },
                 ## all other cases:
                 .sXaug_EigenDense_QRP_Chol_scaled(sXaug,which=which,szAug=szAug,B=B)
  )
  return(resu)
}

