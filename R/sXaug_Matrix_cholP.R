# 'constructor' for sXaug_Matrix_cholP_scaled object
# from Xaug which already has a *scaled* ZAL 
def_sXaug_Matrix_cholP_scaled <- function(Xaug,weight_X,w.ranef,H_global_scale) {
  n_u_h <- length(w.ranef)
  Xrows <- n_u_h+seq(length(weight_X))
  #Xaug[Xrows,] <- Diagonal(x=weight_X) %*% Xaug[Xrows,] ## apparemment sweep vraiment pas efficace sur Matrix
  which_i_affected_rows <- Xaug@i>(n_u_h-1L)
  Xaug@x[which_i_affected_rows] <- Xaug@x[which_i_affected_rows]*weight_X[Xaug@i[which_i_affected_rows]-n_u_h+1L] 
  resu <- structure(Xaug,
                    sXaug_Matrix_cholP_scaled=create_sXaug_Matrix_cholP_scaled(),
                    get_from="get_from_MME.sXaug_Matrix_cholP_scaled",
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
create_sXaug_Matrix_cholP_scaled <- function() {
  BLOB <- list() ## single list for all objects to be kep in memory
  locfn <- function(sXaug,which="",szAug=NULL,B=NULL) {
    if (is.null(BLOB$R_scaled)) {
      BLOB$CHMfactor <<- Matrix::Cholesky(Matrix::crossprod(sXaug),LDL=FALSE) # P'LL'P = P'R'RP for R=L' OK
      # : stores the whole CHMfactor object which has useful methods including solve()
      BLOB$perm <<- as(BLOB$CHMfactor,"pMatrix")@perm 
      BLOB$sortPerm <<- sort.list(BLOB$perm) 
      BLOB$R_scaled <<- t(as(BLOB$CHMfactor,"sparseMatrix"))
      ## t_Q_scaled needed for the solve(); but we don't request (t) Q from Eigen bc it is terribly slow
      #BLOB$t_Q_scaled <<- blob$tQ #[,seq_ncol]) ## removing cols is essential for all leverage-type comput
      if ( ! is.null(szAug)) return(solve(BLOB$CHMfactor,Matrix::crossprod(sXaug, szAug))) 
    } 
    if ( is.null(BLOB$t_Q_scaled) && 
         ( ## ! is.null(szAug) || 
           which %in% c("t_Q_scaled","hatval")) ) {
      BLOB$t_Q_scaled <<- drop0(Matrix::solve(t(BLOB$R_scaled),t(sXaug[,BLOB$perm])))
    }
    if (is.null(BLOB$sortPerm)) BLOB$sortPerm <<- sort.list(BLOB$perm) ## need for most computations
    if ( ! is.null(szAug)) {
      if (is.null(BLOB$t_Q_scaled)) {
        ## many(solve(augmented linear equations) in GLMM, avoid t_Q_computation there)
        return(Matrix::solve(BLOB$CHMfactor,Matrix::crossprod(sXaug, szAug)))
      } else {
        return(Matrix::solve(BLOB$R_scaled, 
                     (BLOB$t_Q_scaled %*% szAug))[BLOB$sortPerm,,drop=FALSE]) 
      }
      ## valid alternative but requesting the computation of t(Q):
      #return(solve(BLOB$R_scaled[,BLOB$sortPerm], BLOB$t_Q_scaled %*% szAug)) ## sort the vector instaed of the matrix?
    }
    # ELSE
    if ( is.null(BLOB$R_R_v)  && 
         (which %in% c("absdiag_R_v","hatval_Z","d2hdv2","solve_d2hdv2"))
    ) {
      seq_n_u_h <- seq_len(attr(sXaug,"n_u_h"))
      CHMfactor_d2hdv2 <- Matrix::Cholesky(Matrix::crossprod(BLOB$R_scaled[,BLOB$sortPerm[ seq_n_u_h ]] ) )
      ## provides a QRP factorization of R_scaled[,sortPerm[seq_n_u_h]] so that 
      ## X[,seq_n_u_h]=Q Q_R_v R_R_v P_R_v
      BLOB$R_R_v <<- t(as(CHMfactor_d2hdv2,"sparseMatrix"))
      BLOB$perm_R_v <<- as(CHMfactor_d2hdv2,"pMatrix")@perm 
    }
    if ( is.null(BLOB$d2hdv2) && 
         (which =="d2hdv2" || (which=="solve_d2hdv2" && is.null(BLOB$inv_d2hdv2) && !is.null(B)) ) 
    ) {
      # don't forgetthat the factored matrix is not the augmented design matrix ! hence w.ranef needed here
      w.ranef <- attr(sXaug,"w.ranef")
      w_R_R_v <- Matrix_times_Dvec(BLOB$R_R_v,sqrt(w.ranef)[BLOB$perm_R_v])
      if (is.null(BLOB$sortPerm_R_v)) BLOB$sortPerm_R_v <<- sort.list(BLOB$perm_R_v)
      BLOB$d2hdv2 <<- - Matrix::crossprod(w_R_R_v)[BLOB$sortPerm_R_v,BLOB$sortPerm_R_v]
    }
    ## return()'s:
    if (which %in% c("solve_d2hdv2")) {
      seq_n_u_h <- seq_len(attr(sXaug,"n_u_h"))
      #return(solve(BLOB$CHMfactor_d2hdv2, t(BLOB$R_scaled[,BLOB$sortPerm[ seq_n_u_h ]]) %*%  B))
      # don't forget that the factored matrix is not the augmented design matrix ! hence w.ranef needed here
      if (is.null(B)) {
        if (is.null(BLOB$inv_d2hdv2)) {
          w.ranef <- attr(sXaug,"w.ranef")
          if (is.null(BLOB$sortPerm_R_v)) BLOB$sortPerm_R_v <<- sort.list(BLOB$perm_R_v)
          w_R_R_v <- Matrix_times_Dvec(BLOB$R_R_v,sqrt(w.ranef)[sort.list(BLOB$sortPerm_R_v)])
          BLOB$inv_d2hdv2 <<- - Matrix::chol2inv(w_R_R_v)[BLOB$sortPerm_R_v,BLOB$sortPerm_R_v] 
        }
        if (is.null(B)) {
          return(BLOB$inv_d2hdv2)
        } else return(BLOB$inv_d2hdv2 %*% B) ## FR->FR could add speific code for Diagonal B ?
      } else { ## solve (Matrix,vector)
        if (is.null(BLOB$inv_d2hdv2)) {
          return( Matrix::solve(BLOB$d2hdv2, B)) 
        } else return(BLOB$inv_d2hdv2 %*% B)
      } 
    } else if (which=="hatval_Z") {
      if (is.null(BLOB$hatval_Z_)) {
        # X[,cols] = Q R P[,cols] = Q q r p => t(Q q) given by:
        if (is.null(BLOB$t_Qq_scaled)){
          BLOB$t_Qq_scaled <<- solve(t(BLOB$R_R_v),
                                     t(sXaug[, seq_len(attr(sXaug,"n_u_h"))[BLOB$perm_R_v] ]))
        } 
        tmp <- BLOB$t_Qq_scaled
        tmp@x <- tmp@x^2
        tmp <- colSums(tmp)
        n_u_h <- attr(sXaug,"n_u_h")
        phipos <- n_u_h+seq_len(nrow(sXaug)-n_u_h)
        BLOB$hatval_Z_ <<-  list(lev_lambda=tmp[-phipos],lev_phi=tmp[phipos])
      }
      return(BLOB$hatval_Z_)
    } else if (which=="hatval") {
      if (is.null(BLOB$hatval_ZX) ) {
        tmp <- BLOB$t_Q_scaled
        tmp@x <- tmp@x^2
        BLOB$hatval_ZX <<- colSums(tmp)
      }
      return(BLOB$hatval_ZX)
    } else if (which=="t_Q_scaled") {
      return(BLOB$t_Q_scaled)
    } else if (which=="R_scaled") {
      return(BLOB$R_scaled)
    } else if (which=="d2hdv2") {
      return(BLOB$d2hdv2)
    } else if (which=="diag_pRtRp") {
      if (is.null(BLOB$diag_pRtRp)) {
        tmp <- BLOB$R_scaled
        tmp@x <- tmp@x^2
        BLOB$diag_pRtRp <<- colSums(tmp)[BLOB$sortPerm]
      }
      return(BLOB$diag_pRtRp)
    } else if (which=="sortPerm") {
      return(BLOB$sortPerm)
    } else if (which %in% c("absdiag_R_v")) {
      if (is.null(BLOB$absdiag_R_v)) 
        BLOB$absdiag_R_v <<- abs(diag(BLOB$R_R_v)) 
      return(BLOB$absdiag_R_v)
    } else if (which=="beta_cov") { 
      beta_cov <- Matrix::chol2inv(BLOB$R_scaled)
      beta_pos <- attr(sXaug,"n_u_h")+seq_len(attr(sXaug,"pforpv"))
      sP_beta_pos <- BLOB$sortPerm[beta_pos]
      beta_cov <- as.matrix(beta_cov[sP_beta_pos,sP_beta_pos]) * attr(sXaug,"H_global_scale")
      colnames(beta_cov) <- colnames(sXaug)[beta_pos]
      return(beta_cov)
    } else if (which=="beta_v_cov_from_wAugX") { ## using a weighted Henderson's augmented design matrix, not a true sXaug  
      beta_v_cov <- Matrix::chol2inv(BLOB$R_scaled)[BLOB$sortPerm,BLOB$sortPerm]
      # this tends to be dense bc v_h estimates covary (example: wafers)
      # otherwise, dropO(,tol=...), + some fix in summary.HLfit for matrix[] <- Matrix assignment, would be useful.  
      return(beta_v_cov)
    } 
    stop("invalid 'which' value.")
  }
  environment(locfn) <- list2env(list(BLOB=list()),
                                 parent=environment(def_sXaug_Matrix_cholP_scaled))
  return(locfn)
} 

# trace("get_from.sXaug_Matrix_cholP_scaled",tracer=quote(if(which=="hatval_Z") debug(attr(sXaug,"sXaug_Matrix_cholP_scaled"))))
get_from_MME.sXaug_Matrix_cholP_scaled <- function(sXaug,which="",szAug=NULL,B=NULL,
                                        damping, LMrhs, ...) {
  resu <- switch(which,
                 "logdet_R_scaled_v" = {
                   absdiag_R_v <- attr(sXaug,"sXaug_Matrix_cholP_scaled")(sXaug,which="absdiag_R_v")
                   sum(log(absdiag_R_v))
                 },
                 "logdet_R_scaled_b_v" = {
                   R_scaled <- attr(sXaug,"sXaug_Matrix_cholP_scaled")(sXaug,which="R_scaled")
                   sum(log(abs(diag(R_scaled))))
                 },
                 "logdet_sqrt_d2hdv2" = {
                   w.ranef <- attr(sXaug,"w.ranef")
                   absdiag_R_v <- attr(sXaug,"sXaug_Matrix_cholP_scaled")(sXaug,which="absdiag_R_v")
                   sum(log(sqrt(w.ranef) * absdiag_R_v))
                 },
                 "LevMar_step" = {
                   ## FR->FR probably not the most elegant implementation 
                   R_scaled <- attr(sXaug,"sXaug_Matrix_cholP_scaled")(sXaug,which="R_scaled")
                   sortPerm <- attr(sXaug,"sXaug_Matrix_cholP_scaled")(sXaug,which="sortPerm")
                   diag_pRtRp <- attr(sXaug,"sXaug_Matrix_cholP_scaled")(sXaug,which="diag_pRtRp")
                   dampDpD <- damping*diag_pRtRp ## NocedalW p. 266
                   # Extend the X in X'X = P'R'RP: 
                   Raug <- rbind(R_scaled[,sortPerm], Diagonal(x=sqrt(dampDpD)))
                   RRblob <- Matrix::qr(Raug)
                   RRR <- Matrix::qrR(RRblob,backPermute = FALSE)
                   RRsP <- sort.list(RRblob@q) 
                   resu <- list(dVscaled_beta= Matrix::chol2inv(RRR)[RRsP,RRsP] %*% LMrhs,
                                rhs=LMrhs, dampDpD = dampDpD)
                   return(resu)
                 },
                 # "logdet_sqrt_d2hdBV2" = {
                 #   w.ranef <- attr(sXaug,"w.ranef")
                 #   R_scaled <- attr(sXaug,"sXaug_Matrix_cholP_scaled")(sXaug,which="R_scaled")
                 #   sum(log(abs(diag(R_scaled)))) + sum(log(w.ranef))/2 - attr(sXaug,"pforpv")*log(H_global_scale)/2 
                 # },
                 "logdet_sqrt_d2hdbeta2" = {
                   H_global_scale <- attr(sXaug,"H_global_scale")
                   R_scaled <- attr(sXaug,"sXaug_Matrix_cholP_scaled")(sXaug,which="R_scaled")
                   absdiag_R_v <- attr(sXaug,"sXaug_Matrix_cholP_scaled")(sXaug,which="absdiag_R_v")
                   sum(log(abs(diag(R_scaled)))) - sum(log(absdiag_R_v)) - attr(sXaug,"pforpv")*log(H_global_scale)/2 
                 },
                 ## all other cases:
                 attr(sXaug,"sXaug_Matrix_cholP_scaled")(sXaug,which=which,szAug=szAug,B=B)
  )
  return(resu)
}
