# 'constructor' for sXaug_Eigen_QR (unpivoted QR factorization) object
# from Xaug which already has a *scaled* ZAL 
def_sXaug_EigenDense_QR_scaled <- function(Xaug,weight_X,w.ranef,H_global_scale) {
  n_u_h <- length(w.ranef)
  Xrows <- n_u_h+seq(length(weight_X))
  Xaug[Xrows,] <- .sweepZ1Wwrapper(Xaug[Xrows,,drop=FALSE],weight_X) # sweep(TT,MARGIN=1,sqrt.ww,`*`) # rWW %*%TT
  resu <- structure(Xaug,
                    get_from="get_from_MME.sXaug_EigenDense_QR_scaled",
                    BLOB=list2env(list(), parent=environment(.sXaug_EigenDense_QR_scaled)),
                    w.ranef=w.ranef,
                    n_u_h=n_u_h,
                    pforpv=ncol(Xaug)-n_u_h,
                    H_global_scale=H_global_scale
  )
  class(resu) <- c(class(resu),"sXaug_EigenDense_QR_scaled")
  return( resu ) 
}

# single 'internal method' for sXaug_Eigen_QR objects
# if QR not previously available, szAug=NULL: constructs the QR factorization, return value depends on 'which'
# if QR not previously available, szAug ! NULL:   ..........................   returns the solution
# if which="solve_d2hdv2", B=NULL: returns solve(d2hdv2)
# if which="solve_d2hdv2", B ! NULL: returns solve(d2hdv2,B)
.sXaug_EigenDense_QR_scaled <- function(sXaug,which="",szAug=NULL,B=NULL) {
  BLOB <- attr(sXaug,"BLOB") ## an environment
  if (is.null(BLOB$R_scaled)) {
    lmwithqr <- lmwithQR(sXaug,szAug,returntQ=TRUE,returnR=TRUE)
    for (st in names(lmwithqr)) BLOB[[st]] <- lmwithqr[[st]]
    if ( ! is.null(szAug)) return(BLOB$coef)   
  } 
  ## else QR was already available
  if ( ! is.null(szAug)) return(backsolve(BLOB$R_scaled,BLOB$t_Q_scaled %*% szAug))
  # ELSE
  # test sur B. note that NCOL(NULL)=1. 
  if ( is.null(BLOB$inv_d2hdv2) && which=="solve_d2hdv2" && is.null(B)) {
    w.ranef <- attr(sXaug,"w.ranef")
    seq_n_u_h <- seq_len(attr(sXaug,"n_u_h"))
    R_Md2hdv2 <- BLOB$R_scaled[seq_n_u_h, seq_n_u_h] %*% diag(sqrt(w.ranef)) ## : such that R'R gives -d2hdv2
    BLOB$inv_d2hdv2 <- - chol2inv(R_Md2hdv2) ## valid bc no pivoting
  }
  if ( is.null(BLOB$d2hdv2) && ## important to test this for both 'which'
       (which =="d2hdv2" || (which=="solve_d2hdv2" && is.null(BLOB$inv_d2hdv2) ) ) 
  ) {
    w.ranef <- attr(sXaug,"w.ranef")
    seq_n_u_h <- seq_len(attr(sXaug,"n_u_h"))
    R_Md2hdv2 <- BLOB$R_scaled[seq_n_u_h, seq_n_u_h] %*% diag(sqrt(w.ranef)) ## : such that R'R gives -d2hdv2
    BLOB$d2hdv2 <- - crossprodCpp(R_Md2hdv2,NULL)  # w r'r w not ZtWZ     
  }
  if (which=="t_Q_scaled") { 
    return(BLOB$t_Q_scaled)
  } else if (which=="hatval") {
    return(colSums(BLOB$t_Q_scaled^2))
  } else if (which=="hatval_Z") {
    if (is.null(BLOB$hatval_Z_) ) {
      tmp <- colSums(BLOB$t_Q_scaled[ seq_len(attr(sXaug,"n_u_h")) ,]^2)
      n_u_h <- attr(sXaug,"n_u_h")
      phipos <- n_u_h+seq_len(nrow(sXaug)-n_u_h)
      BLOB$hatval_Z_ <- list(lev_lambda=tmp[-phipos],lev_phi=tmp[phipos])
    }
    return(BLOB$hatval_Z_)
  } else if (which=="R_scaled") {
    return(BLOB$R_scaled)
  } else if (which=="d2hdv2") {
    return(BLOB$d2hdv2)
  } else if (which=="diag_RtR") {
    if (is.null(BLOB$diag_RtR)) BLOB$diag_RtR <- colSums(BLOB$R_scaled^2)
    return(BLOB$diag_RtR)
  } else if (which=="solve_d2hdv2") { ## R'R[seq_n_u_h, seq_n_u_h] gives -d2hd_scaled_v2.
    ## FR->FR there must be better code for CAR models, using a better decomp than R block from QR ?
    if ( is.null(B) ) {
      return(BLOB$inv_d2hdv2)
    } else { ## solve (matrix,y)
      if (is.null(BLOB$inv_d2hdv2)) {
        resu <- tryCatch(solve(BLOB$d2hdv2, B),error=.minimalErrorHandler)
        if (inherits(resu,"try-error")) {
          w.ranef <- attr(sXaug,"w.ranef")
          seq_n_u_h <- seq_len(attr(sXaug,"n_u_h"))
          R_Md2hdv2 <- BLOB$R_scaled[seq_n_u_h, seq_n_u_h] %*% diag(sqrt(w.ranef)) ## : such that R'R gives -d2hdv2
          BLOB$inv_d2hdv2 <- - chol2inv(R_Md2hdv2) ## valid bc no pivoting
          return(BLOB$inv_d2hdv2 %*% B)
        } else return(resu)
      } else {
        return(BLOB$inv_d2hdv2 %*% B) ## inv_d2hdv2 tends to be dense
      }
    }
  } else if (which=="beta_cov") { 
    beta_cov <- chol2inv(BLOB$R_scaled) ## valid bc no pivoting
    beta_pos <- attr(sXaug,"n_u_h")+seq_len(attr(sXaug,"pforpv"))
    beta_cov <- beta_cov[beta_pos,beta_pos,drop=FALSE] * attr(sXaug,"H_global_scale")
    colnames(beta_cov) <- colnames(sXaug)[beta_pos]
    return(beta_cov)
  } else if (which=="beta_v_cov_from_sXaug") { ## not yet used
    beta_v_cov <- chol2inv(BLOB$R_scaled) ## valid bc no pivoting
    diagscalings <- diag(x=sqrt(c(1/attr(sXaug,"w.ranef"),attr(sXaug,"H_global_scale"))))
    beta_v_cov <- diagscalings %*% beta_v_cov %*% diagscalings
    toXZ0I <- c(attr(sXaug,"n_u_h")+seq_len(attr(sXaug,"pforpv")),seq_len(attr(sXaug,"n_u_h")))
    beta_v_cov <- beta_v_cov[toXZ0I,toXZ0I]
    return(beta_v_cov)
  } else if (which=="beta_v_cov_from_wAugX") { ## using Henderson's augmented design matrix, not a true sXaug  
    beta_v_cov <- chol2inv(BLOB$R_scaled) ## valid bc no pivoting
    return(beta_v_cov)
  }
  stop("invalid 'which' value.")
}
  

# trace("get_from.sXaug_EigenDense_QR_scaled",tracer=quote(if(which=="hatval_Z") debug(attr(sXaug,"sXaug_EigenDense_QR_scaled"))))
# trace("get_from.sXaug_EigenDense_QR_scaled",tracer=quote(print(which)))
get_from_MME.sXaug_EigenDense_QR_scaled <- function(sXaug,which="",szAug=NULL,B=NULL,
                                    damping, LMrhs, ...) {
  R_scaled <- .sXaug_EigenDense_QR_scaled(sXaug,which="R_scaled")
  resu <- switch(which,
                 "logdet_R_scaled_v" = {
                   sum(log(abs(diag(R_scaled)[ seq_len(attr(sXaug,"n_u_h")) ]))) 
                 },
                 "logdet_R_scaled_b_v" = {
                   sum(log(abs(diag(R_scaled))))
                 },
                 "logdet_sqrt_d2hdv2" = {
                   w.ranef <- attr(sXaug,"w.ranef")
                   sum(log(sqrt(w.ranef) * abs(diag(R_scaled)[ seq_len(attr(sXaug,"n_u_h")) ]))) 
                 },
                 "logdet_sqrt_d2hdbeta2" = {
                   H_global_scale <- attr(sXaug,"H_global_scale")
                   sum(log(abs(diag(R_scaled)[ - seq_len(attr(sXaug,"n_u_h")) ]))) - attr(sXaug,"pforpv")*log(H_global_scale)/2 
                 },
                 "LevMar_step" = {
                   diag_RtR <- .sXaug_EigenDense_QR_scaled(sXaug,which="diag_RtR")
                   dampDpD <- damping*diag_RtR ## NocedalW p. 266
                   if (FALSE) {
                     ## version QR inspree de More (mais sans les Givens bien choisis)
                     Raug <- rbind(R_scaled, diag(x=sqrt(dampDpD)))
                     resu <- LevMar_cpp(Raug,LMrhs=LMrhs)
                     resu$dampDpD <- dampDpD
                   } else {
                     ## version chol assez directe (a traduire en c++?)
                     AtAdDpD <- crossprodCpp(R_scaled,NULL)
                     nc <- ncol(AtAdDpD)
                     diagPos <- seq.int(1L,nc^2,nc+1L)
                     AtAdDpD[diagPos] <- AtAdDpD[diagPos] + dampDpD 
                     cholR <- Rcpp_chol_R(AtAdDpD)$R
                     dVscaled_beta <- backsolve(cholR,x=LMrhs,transpose = TRUE)
                     dVscaled_beta <- backsolve(cholR,x=dVscaled_beta)
                     resu <- list(dVscaled_beta=dVscaled_beta, rhs=LMrhs, dampDpD = dampDpD)
                   }
                   return(resu)
                 },
                 ## all other cases:
                 .sXaug_EigenDense_QR_scaled(sXaug,which=which,szAug=szAug,B=B)
  )
  return(resu)
}
