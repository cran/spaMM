# 'constructor' for sXaug_EigenSparse_LDL_scaled object
# from Xaug which already has a *scaled* ZAL 
def_sXaug_EigenSparse_LDL_scaled <- function(Xaug,weight_X,w.ranef,H_global_scale) {
  n_u_h <- length(w.ranef)
  Xrows <- n_u_h+seq(length(weight_X))
  #Xaug[Xrows,] <- Diagonal(x=weight_X) %*% Xaug[Xrows,] ## apparemment sweep vraiment pas efficace sur Matrix
  which_i_affected_rows <- Xaug@i>(n_u_h-1L)
  Xaug@x[which_i_affected_rows] <- Xaug@x[which_i_affected_rows]*weight_X[Xaug@i[which_i_affected_rows]-n_u_h+1L] 
  resu <- structure(Xaug,
                    sXaug_EigenSparse_LDL_scaled=create_sXaug_EigenSparse_LDL_scaled(),
                    get_from="get_from_MME.sXaug_EigenSparse_LDL_scaled",
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
create_sXaug_EigenSparse_LDL_scaled <- function() {
  BLOB <- list() ## single list for all objects to be kep in memory
  locfn <- function(sXaug,which="",szAug=NULL,B=NULL) {
    if (is.null(BLOB$D_scaled)) {
      BLOB <<- lmwith_sparse_LDLp(sXaug, szAug, returntQ=FALSE, returnR=TRUE,pivot=FALSE) 
      if ( ! is.null(szAug)) return(BLOB$coef) 
    } 
    if ( is.null(BLOB$t_Q_scaled) && 
         ( ## ! is.null(szAug) || 
           which %in% c("t_Q_scaled","hatval","hatval_Z")) ) {
      # le solve() est bizarrement lent... but then pivoting may be useful...
      BLOB$t_Q_scaled <<- Diagonal(x=1/sqrt(BLOB$D_scaled)) %*% Matrix::drop0(Matrix::solve(t(BLOB$U_scaled), t(sXaug))) 
    }
    if ( ! is.null(szAug)) {
      if (is.null(BLOB$t_Q_scaled)) {
        ## many(solve(augmented linear equations) in GLMM, avoid t_Q_computation there)
        return(Matrix::solve(BLOB$XtX, Matrix::crossprod(sXaug, szAug))) 
      } else {
        return(Matrix::solve(BLOB$U_scaled, 
                              (BLOB$t_Q_scaled %*% szAug)/sqrt(BLOB$D_scaled))) 
      }
    }
    if ( is.null(BLOB$d2hdv2) && 
         (which =="d2hdv2" || (which=="solve_d2hdv2" && is.null(BLOB$inv_d2hdv2) && !is.null(B) ) ) 
    ) {
      w.ranef <- attr(sXaug,"w.ranef")
      seq_n_u_h <- seq_len(attr(sXaug,"n_u_h"))
      BLOB$d2hdv2 <<- - Diagonal(x=sqrt(w.ranef)) %*% Matrix_times_Dvec(BLOB$XtX[seq_n_u_h, seq_n_u_h], sqrt(w.ranef))    
    }
    # ELSE
    if (which=="t_Q_scaled") {
      return(BLOB$t_Q_scaled)
    } else if (which=="hatval") {
      if (is.null(BLOB$hatval_ZX) ) {
        tmp <- BLOB$t_Q_scaled
        tmp@x <- tmp@x^2
        if (is.null(BLOB$hatval_ZX) ) BLOB$hatval_ZX <<-  colSums(tmp)
      }
      return(BLOB$hatval_ZX)
    } else if (which=="hatval_Z") {
      if (is.null(BLOB$hatval_Z_)) {
        tmp <- BLOB$t_Q_scaled[ seq_len(attr(sXaug,"n_u_h")) ,]
        tmp@x <- tmp@x^2
        tmp <- colSums(tmp)
        n_u_h <- attr(sXaug,"n_u_h")
        phipos <- n_u_h+seq_len(nrow(sXaug)-n_u_h)
        BLOB$hatval_Z_ <<-  list(lev_lambda=tmp[-phipos],lev_phi=tmp[phipos])
      }
      return(BLOB$hatval_Z_)
    } else if (which=="D_scaled") { 
      return(BLOB$D_scaled)
    } else if (which=="R_scaled") { 
      if (is.null(BLOB$R_scaled)) BLOB$R_scaled <<- Diagonal(x=sqrt(abs(BLOB$D_scaled))) %*% BLOB$U_scaled 
      return(BLOB$R_scaled)
    } else if (which=="diag_pRtRp") {
      return(diag(BLOB$XtX))
    } else if (which=="sortPerm") {
      return(BLOB$sortPerm)
    } else if (which=="beta_cov") {
      if (is.null(BLOB$R_scaled)) BLOB$R_scaled <<- Diagonal(x=sqrt(abs(BLOB$D_scaled))) %*% BLOB$U_scaled 
      beta_cov <- as.matrix(Matrix::chol2inv(BLOB$R_scaled))
      beta_pos <- attr(sXaug,"n_u_h")+seq_len(attr(sXaug,"pforpv"))
      beta_cov <- beta_cov[beta_pos,beta_pos,drop=FALSE] * attr(sXaug,"H_global_scale")
      colnames(beta_cov) <- colnames(sXaug)[beta_pos]
      return(beta_cov)
    } else if (which=="beta_v_cov_from_wAugX") { ## using a weighted Henderson's augmented design matrix, not a true sXaug  
      if (is.null(BLOB$R_scaled)) BLOB$R_scaled <<- Diagonal(x=sqrt(abs(BLOB$D_scaled))) %*% BLOB$U_scaled 
      beta_v_cov <- as.matrix(Matrix::chol2inv(BLOB$R_scaled))
      # this tends to be dense bc v_h estimates covary (example: wafers)
      # otherwise, dropO(,tol=...), + some fix in summary.HLfit for matrix[] <- Matrix assignment, would be useful.  
      return(beta_v_cov)
    } else if (which %in% c("solve_d2hdv2")) {
      if (is.null(B)) {
        if (is.null(BLOB$inv_d2hdv2)) {
          w.ranef <- attr(sXaug,"w.ranef")
          seq_n_u_h <- seq_len(attr(sXaug,"n_u_h"))
          w_R_R_v <- Matrix_times_Dvec(BLOB$U_scaled[seq_n_u_h,seq_n_u_h], sqrt(w.ranef[seq_n_u_h]))
          w_R_R_v  <- Diagonal(x=sqrt(abs(BLOB$D_scaled[seq_n_u_h]))) %*% w_R_R_v
          BLOB$inv_d2hdv2 <<- - Matrix::chol2inv(w_R_R_v) 
        }
        if (is.null(B)) {
          return(BLOB$inv_d2hdv2)
        } else return(BLOB$inv_d2hdv2 %*% B) ## FR->FR could add speific code for Diagonal B ?
      } else { ## solve (Matrix,vector)
        if (is.null(BLOB$inv_d2hdv2)) {
          return( Matrix::solve(BLOB$d2hdv2, B)) 
        } else return(BLOB$inv_d2hdv2 %*% B)
      }
    } 
    stop("invalid 'which' value.")
  }
  environment(locfn) <- list2env(list(BLOB=list()),
                                 parent=environment(def_sXaug_EigenSparse_LDL_scaled))
  return(locfn)
} 

# trace("get_from.sXaug_EigenSparse_LDL_scaled",tracer=quote(if(which=="hatval_Z") debug(attr(sXaug,"sXaug_EigenSparse_LDL_scaled"))))
get_from_MME.sXaug_EigenSparse_LDL_scaled <- function(sXaug,which="",szAug=NULL,B=NULL,
                                        damping, LMrhs, ...) {
  D_scaled <- attr(sXaug,"sXaug_EigenSparse_LDL_scaled")(sXaug,which="D_scaled")
  resu <- switch(which,
                 "logdet_R_scaled_v" = {
                   sum(log(abs(D_scaled[ seq_len(attr(sXaug,"n_u_h")) ])))/2
                 },
                 "logdet_R_scaled_b_v" = {
                   sum(log(abs(D_scaled)))/2
                 },
                 "logdet_sqrt_d2hdv2" = {
                   w.ranef <- attr(sXaug,"w.ranef")
                   sum(log(w.ranef * D_scaled[ seq_len(attr(sXaug,"n_u_h")) ]))/2
                 },
                 "LevMar_step" = {
                   ## FR->FR probably not the most elegant implementation 
                   R_scaled <- attr(sXaug,"sXaug_EigenSparse_LDL_scaled")(sXaug,which="R_scaled")
                   diag_pRtRp <- attr(sXaug,"sXaug_EigenSparse_LDL_scaled")(sXaug,which="diag_pRtRp")
                   dampDpD <- damping*diag_pRtRp ## NocedalW p. 266
                   # Extend the X in X'X = P'R'RP: 
                   Raug <- rbind(R_scaled, Diagonal(x=sqrt(dampDpD)))
                   RRblob <- Matrix::qr(Raug)
                   RRR <- Matrix::qrR(RRblob,backPermute = FALSE)
                   RRsP <- sort.list(RRblob@q) 
                   resu <- list(dVscaled_beta= Matrix::chol2inv(RRR)[RRsP,RRsP] %*% LMrhs,
                                rhs=LMrhs, dampDpD = dampDpD)
                   return(resu)
                 },
                 "logdet_sqrt_d2hdbeta2" = {
                   H_global_scale <- attr(sXaug,"H_global_scale")
                   sum(log(abs(D_scaled[ - seq_len(attr(sXaug,"n_u_h"))])))/2 - attr(sXaug,"pforpv")*log(H_global_scale)/2 
                 },
                 ## all other cases:
                 attr(sXaug,"sXaug_EigenSparse_LDL_scaled")(sXaug,which=which,szAug=szAug,B=B)
  )
  return(resu)
}
