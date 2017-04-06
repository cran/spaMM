# 'constructor' for sXaug_Eigen_QR (unpivoted QR factorization) object
# from Xaug which already has a *scaled* ZAL 
def_sXaug_EigenDense_QRP_scaled <- function(Xaug,weight_X,w.ranef,H_global_scale) {
  n_u_h <- length(w.ranef)
  Xrows <- n_u_h+seq(length(weight_X))
  Xaug[Xrows,] <- .sweepZ1Wwrapper(Xaug[Xrows,,drop=FALSE],weight_X) # sweep(TT,MARGIN=1,sqrt.ww,`*`) # rWW %*%TT
  resu <- structure(Xaug,
                    sXaug_EigenDense_QRP_scaled=create_sXaug_EigenDense_QRP_scaled(),
                    get_from="get_from_MME.sXaug_EigenDense_QRP_scaled",
                    w.ranef=w.ranef,
                    n_u_h=n_u_h,
                    pforpv=ncol(Xaug)-n_u_h,
                    H_global_scale=H_global_scale
  )
  class(resu) <- c(class(resu),"sXaug_EigenDense_QRP_scaled")
  return( resu ) 
}

# single 'internal method' (with <<- ) for sXaug_Eigen_QR objects
# if QR not previously available, szAug=NULL: constructs the QR factorization, return value depends on 'which'
# if QR not previously available, szAug ! NULL:   ..........................   returns the solution
# if which="solve_d2hdv2", B=NULL: returns solve(d2hdv2)
# if which="solve_d2hdv2", B ! NULL: returns solve(d2hdv2,B)
create_sXaug_EigenDense_QRP_scaled <- function() {
  BLOB <- list()
  locfn <- function(sXaug,which="",szAug=NULL,B=NULL) {
    if (is.null(BLOB$R_scaled)) {
      BLOB <<- lmwithQRP(sXaug,yy=szAug,returntQ=FALSE,returnR=TRUE) ## using RcppEigen; szAug may be NULL
      BLOB$perm <<- BLOB$perm +1L
      n_u_h <- attr(sXaug,"n_u_h")
      seq_n_u_h <- seq(n_u_h)
      BLOB$sortPerm <<- sort.list(BLOB$perm) 
      BLOB$u_h_cols_on_left <<- (max(BLOB$sortPerm[ seq_n_u_h ])==n_u_h)
      if ( ! is.null(szAug)) return(BLOB$coef)   
    } 
    if ( is.null(BLOB$t_Q_scaled) && 
         ( ! is.null(szAug) || 
           which %in% c("t_Q_scaled","hatval") ||
           (which=="hatval_Z" && BLOB$u_h_cols_on_left)
         ) ) {
      BLOB$t_Q_scaled <<- backsolve(BLOB$R_scaled,t(sXaug[,BLOB$perm]),transpose=TRUE)
    }
    if ( ! is.null(szAug)) {
      return(backsolve(r=BLOB$R_scaled, x=BLOB$t_Q_scaled %*% szAug)[BLOB$sortPerm,,drop=FALSE]) 
    }
    # ELSE
    if (which=="hatval_Z" && BLOB$u_h_cols_on_left) {
      # do nothing, $t_Q_scaled is already available and R_R_v not needed
    } else {
      if ( is.null(BLOB$R_R_v) && which %in% c("absdiag_R_v","d2hdv2","solve_d2hdv2","hatval_Z")) {
        seq_n_u_h <- seq_len(attr(sXaug,"n_u_h"))
        if (BLOB$u_h_cols_on_left) {  # condition is true for small COMPoisson example but not more generally
          BLOB$R_R_v <<- BLOB$R_scaled[seq_n_u_h,seq_n_u_h]
          BLOB$perm_R_v <<- BLOB$perm[ seq_n_u_h ]
        } else {
          Rblob <- lmwithQRP(BLOB$R_scaled[,BLOB$sortPerm[ seq_n_u_h ]],yy=NULL,returntQ=FALSE,returnR=TRUE)
          ## provides a QRP factorization of R_scaled[,sortPerm[seq_n_u_h]] so that 
          ## X[,seq_n_u_h]=Q Q_R_v R_R_v P_R_v
          BLOB$R_R_v <<- Rblob$R_scaled
          BLOB$perm_R_v <<- Rblob$perm +1L
        }
      }
    }
    # test sur B. note that NCOL(NULL)=1. 
    if ( is.null(BLOB$inv_d2hdv2) && which=="solve_d2hdv2" && is.null(B)) {
      w.ranef <- attr(sXaug,"w.ranef")
      w_R_R_v <- BLOB$R_R_v %*% diag(sqrt(w.ranef)[BLOB$perm_R_v])
      if (is.null(BLOB$sortPerm_R_v)) BLOB$sortPerm_R_v <<- sort.list(BLOB$perm_R_v)
      BLOB$inv_d2hdv2 <<- - chol2inv(w_R_R_v)[BLOB$sortPerm_R_v,BLOB$sortPerm_R_v] 
    }
    if ( is.null(BLOB$d2hdv2) && ## important to test this for both 'which'
         (which =="d2hdv2" || (which=="solve_d2hdv2" && is.null(BLOB$inv_d2hdv2) ) ) 
    ) {
      # don't forget that the factored matrix is not the augmented design matrix ! hence w.ranef needed here
      w.ranef <- attr(sXaug,"w.ranef")
      w_R_R_v <- BLOB$R_R_v %*% diag(sqrt(w.ranef)[BLOB$perm_R_v])
      if (is.null(BLOB$sortPerm_R_v)) BLOB$sortPerm_R_v <<- sort.list(BLOB$perm_R_v)
      BLOB$d2hdv2 <<- - crossprodCpp(w_R_R_v,yy=NULL)[BLOB$sortPerm_R_v,BLOB$sortPerm_R_v]
    }
    if (which=="solve_d2hdv2") { ## R'R[seq_n_u_h, seq_n_u_h] gives -d2hd_scaled_v2.
      ## FR->FR there must be better code for CAR models, using a better decomp than R block from QR ?
      if ( is.null(B) ) {
        return(BLOB$inv_d2hdv2)
      } else { ## solve (matrix,y)
        if (is.null(BLOB$inv_d2hdv2)) {
          return( solve(BLOB$d2hdv2, B)) 
        } else {
          return(BLOB$inv_d2hdv2 %*% B) ## inv_d2hdv2 tends to be dense
        }
      }
    } else if (which=="hatval_Z") {
      ## Pdiag
      if (is.null(BLOB$hatval_Z_)) {
        if (TRUE) {
          # X[,cols] = Q R P[,cols] = Q q r p => t(Q q) given by:
          if (BLOB$u_h_cols_on_left) {
            tmp_t_Qq_scaled <- BLOB$t_Q_scaled[seq_len(attr(sXaug,"n_u_h")),]
          } else {
            tmp_t_Qq_scaled <- backsolve(r=BLOB$R_R_v, x=t( sXaug[, seq_len(attr(sXaug,"n_u_h"))[BLOB$perm_R_v] ] ),
                                         transpose=TRUE)
          }
          tmp <- colSums(tmp_t_Qq_scaled^2)
          n_u_h <- attr(sXaug,"n_u_h")
          phipos <- n_u_h+seq_len(nrow(sXaug)-n_u_h)
          BLOB$hatval_Z_ <<- list(lev_lambda=tmp[-phipos],lev_phi=tmp[phipos])
        } else {
          ## FR->FR try two qr.qy instead of one solve as the try in sparse QRP ? 
        }
      }
      return(BLOB$hatval_Z_)
    } else if (which=="hatval") {
      return(colSums(BLOB$t_Q_scaled^2))
    } else if (which=="R_scaled") {
      return(BLOB$R_scaled)
    } else if (which=="t_Q_scaled") { 
      return(BLOB$t_Q_scaled)
    } else if (which %in% c("absdiag_R_v")) { ## used for logdet_R_scaled_v
      if (is.null(BLOB$absdiag_R_v)) BLOB$absdiag_R_v <<- abs(diag(BLOB$R_R_v)) 
      return(BLOB$absdiag_R_v)
    } else if (which=="diag_pRtRp") { ## extracted for LevMar
      if (is.null(BLOB$diag_pRtRp)) BLOB$diag_pRtRp <<- colSums(BLOB$R_scaled^2)[BLOB$sortPerm]
      return(BLOB$diag_pRtRp)
    } else if (which=="sortPerm") { ## extracted for LevMar
      return(BLOB$sortPerm)
    } else if (which=="beta_cov") { 
      beta_cov <- chol2inv(BLOB$R_scaled)
      beta_pos <- attr(sXaug,"n_u_h")+seq_len(attr(sXaug,"pforpv"))
      sP_beta_pos <- BLOB$sortPerm[beta_pos]
      beta_cov <- beta_cov[sP_beta_pos,sP_beta_pos,drop=FALSE] * attr(sXaug,"H_global_scale")
      colnames(beta_cov) <- colnames(sXaug)[beta_pos]
      return(beta_cov)
    } else if (which=="beta_v_cov_from_wAugX") { ## using a weighted Henderson's augmented design matrix, not a true sXaug  
      beta_v_cov <- chol2inv(BLOB$R_scaled)[BLOB$sortPerm,BLOB$sortPerm]
      # this tends to be dense bc v_h estimates covary (example: wafers)
      # otherwise, dropO(,tol=...), + some fix in summary.HLfit for matrix[] <- Matrix assignment, would be useful.  
      return(beta_v_cov)
    } else if (which=="d2hdv2") {
      return(BLOB$d2hdv2)
    }
    stop("invalid 'which' value.")
  }
  environment(locfn) <- list2env(list(BLOB=list()),
                                 parent=environment(def_sXaug_EigenDense_QRP_scaled))
  return(locfn)
} 

# trace("get_from.sXaug_EigenDense_QRP_scaled",tracer=quote(if(which=="hatval_Z") debug(attr(sXaug,"sXaug_EigenDense_QRP_scaled"))))
# trace("get_from.sXaug_EigenDense_QRP_scaled",tracer=quote(print(which)),print=FALSE)
get_from_MME.sXaug_EigenDense_QRP_scaled <- function(sXaug,which="",szAug=NULL,B=NULL,
                                    damping, LMrhs, ...) {
  #R_scaled <- attr(sXaug,"sXaug_EigenDense_QRP_scaled")(sXaug,which="R_scaled")
  resu <- switch(which,
                 "logdet_R_scaled_v" = {
                   absdiag_R_v <- attr(sXaug,"sXaug_EigenDense_QRP_scaled")(sXaug,which="absdiag_R_v")
                   sum(log(absdiag_R_v))
                 },
                 "logdet_R_scaled_b_v" = {
                   R_scaled <- attr(sXaug,"sXaug_EigenDense_QRP_scaled")(sXaug,which="R_scaled")
                   sum(log(abs(diag(R_scaled))))
                 },
                 "logdet_sqrt_d2hdv2" = {
                   w.ranef <- attr(sXaug,"w.ranef")
                   absdiag_R_v <- attr(sXaug,"sXaug_EigenDense_QRP_scaled")(sXaug,which="absdiag_R_v")
                   sum(log(sqrt(w.ranef) * absdiag_R_v))
                 },
                 "logdet_sqrt_d2hdbeta2" = {
                   H_global_scale <- attr(sXaug,"H_global_scale")
                   R_scaled <- attr(sXaug,"sXaug_EigenDense_QRP_scaled")(sXaug,which="R_scaled")
                   absdiag_R_v <- attr(sXaug,"sXaug_EigenDense_QRP_scaled")(sXaug,which="absdiag_R_v")
                   sum(log(abs(diag(R_scaled)))) - sum(log(absdiag_R_v)) - attr(sXaug,"pforpv")*log(H_global_scale)/2 
                 },
                 "LevMar_step" = {
                   ## FR->FR probably not the most elegant implementation 
                   R_scaled <- attr(sXaug,"sXaug_EigenDense_QRP_scaled")(sXaug,which="R_scaled")
                   sortPerm <- attr(sXaug,"sXaug_EigenDense_QRP_scaled")(sXaug,which="sortPerm")
                   diag_pRtRp <- attr(sXaug,"sXaug_EigenDense_QRP_scaled")(sXaug,which="diag_pRtRp")
                   dampDpD <- damping*diag_pRtRp ## NocedalW p. 266
                   Raug <- rbind(R_scaled[,sortPerm], diag(x=sqrt(dampDpD)))
                   RRblob <- lmwithQRP(Raug,yy=NULL,returntQ=FALSE,returnR=TRUE)
                   RRR <- RRblob$R_scaled
                   RRsP <- sort.list(RRblob$perm) 
                   resu <- list(dVscaled_beta= chol2inv(RRR)[RRsP,RRsP] %*% LMrhs,
                                rhs=LMrhs, dampDpD = dampDpD)
                 },
                 ## all other cases:
                 attr(sXaug,"sXaug_EigenDense_QRP_scaled")(sXaug,which=which,szAug=szAug,B=B)
  )
  return(resu)
}
