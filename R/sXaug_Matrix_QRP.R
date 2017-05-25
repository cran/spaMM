# 'constructor' for sXaug_Eigen_sparse_QR (unpivoted QR factorization) object
# from Xaug which already has a *scaled* ZAL 
def_sXaug_Matrix_QRP_scaled <- function(Xaug,weight_X,w.ranef,H_global_scale) {
  n_u_h <- length(w.ranef)
  Xrows <- n_u_h+seq(length(weight_X))
  ## Bates https://stat.ethz.ch/pipermail/r-help/2010-December/262365.html
  ## "Assignment of submatrices in a sparse matrix can be slow because there is so much checking that needs to be done."
  # Xaug[Xrows,] <- Diagonal(x=weight_X) %*% Xaug[Xrows,] ## 
  which_i_affected_rows <- Xaug@i>(n_u_h-1L)
  Xaug@x[which_i_affected_rows] <- Xaug@x[which_i_affected_rows]*weight_X[Xaug@i[which_i_affected_rows]-n_u_h+1L] 
  ## FR6>FR mais c'est lent aussi! Revoir.
  resu <- structure(Xaug,
                    get_from="get_from_MME.sXaug_Matrix_QRP_scaled",
                    BLOB=list2env(list(), parent=environment(.sXaug_Matrix_QRP_scaled)),
                    w.ranef=w.ranef,
                    n_u_h=n_u_h,
                    pforpv=ncol(Xaug)-n_u_h,
                    H_global_scale=H_global_scale
  )
  ## cannot modify the class attribute... => immediate clumsy code below...
  return( resu ) 
}

# single 'internal method' (with <<- ) for sXaug_Eigen_sparse_QR objects
# if QR not previously available, szAug=NULL: constructs the QR factorization, return value depends on 'which'
# if QR not previously available, szAug ! NULL:   ..........................   returns the solution
# if which="solve_d2hdv2", B=NULL: returns solve(d2hdv2)
# if which="solve_d2hdv2", B ! NULL: returns solve(d2hdv2,B)
.sXaug_Matrix_QRP_scaled <- function(sXaug,which="",szAug=NULL,B=NULL) {
  BLOB <- attr(sXaug,"BLOB") ## an environment
  if (is.null(BLOB$R_scaled)) {
    BLOB$blob <- Matrix::qr(sXaug)
    # sXaug = t(tQ) %*% R[,sP] but then also sXaug = t(tQ)[,sP'] %*% R[sP',sP] for any sP'
    BLOB$perm <- BLOB$blob@q + 1L
    BLOB$R_scaled <- Matrix::qrR(BLOB$blob,backPermute = FALSE)
    ## t_Q_scaled needed for the solve(); but we don't request (t) Q from Eigen bc it is terribly slow
    #BLOB$t_Q_scaled <<- blob$tQ #[,seq_ncol]) ## removing cols is essential for all leverage-type comput
    n_u_h <- attr(sXaug,"n_u_h")
    seq_n_u_h <- seq(n_u_h)
    BLOB$sortPerm <- sort.list(BLOB$perm) 
    BLOB$u_h_cols_on_left <- (max(BLOB$sortPerm[ seq_n_u_h ])==n_u_h)
    if ( ! is.null(szAug)) return(Matrix::qr.coef(BLOB$blob,szAug))   
  } 
  ##if ( ! is.null(szAug)) return(solve(R_scaled[,sortPerm],t_Q_scaled %*% szAug))
  if ( is.null(BLOB$t_Q_scaled) && 
       ( #! is.null(szAug) || 
         which %in% c("t_Q_scaled","hatval") ||
         (which=="hatval_Z" && BLOB$u_h_cols_on_left)
       ) ) {
    BLOB$t_Q_scaled <- Matrix::solve(t(BLOB$R_scaled),t(sXaug[,BLOB$perm]))
  }
  # if (is.null(BLOB$sortPerm)) BLOB$sortPerm <<- sort.list(BLOB$perm) ## need for most computations
  if ( ! is.null(szAug)) {
    if (is.null(BLOB$t_Q_scaled)) {
      return(Matrix::qr.coef(BLOB$blob,szAug)) # ~ solve 
      ## many(solve(augmented linear equations) in GLMM, avoid t_Q_computation there)
      # return(Matrix::solve(crossprod(BLOB$R_scaled), (t(sXaug) %*% szAug)[BLOB$perm])[BLOB$sortPerm,,drop=FALSE])
    } else {
      return(Matrix::solve(BLOB$R_scaled, BLOB$t_Q_scaled %*% szAug)[BLOB$sortPerm,,drop=FALSE]) 
    }
  }
  # ELSE
  if (which=="hatval_Z" && BLOB$u_h_cols_on_left) {
    # do nothing, $t_Q_scaled is already available and R_R_v not needed
  } else {
    if ( is.null(BLOB$R_R_v) && which %in% c("absdiag_R_v","d2hdv2","solve_d2hdv2","hatval_Z")) {
      seq_n_u_h <- seq_len(attr(sXaug,"n_u_h"))
      if (BLOB$u_h_cols_on_left) {  # condition is true for small COMPoisson example but not more generally
        BLOB$R_R_v <- BLOB$R_scaled[seq_n_u_h,seq_n_u_h]
        BLOB$perm_R_v <- BLOB$perm[ seq_n_u_h ] 
      } else {
        #BLOB$Rblob <<- 
        Rblob <- Matrix::qr(BLOB$R_scaled[,BLOB$sortPerm[ seq_n_u_h ]])
        ## provides a QRP factorization of R_scaled[,sortPerm[seq_n_u_h]] so that 
        ## X[,seq_n_u_h]=Q Q_R_v R_R_v P_R_v
        BLOB$R_R_v <- Matrix::qrR(Rblob,backPermute = FALSE)
        BLOB$perm_R_v <- Rblob@q + 1L 
      }
    }
  }
  if ( is.null(BLOB$inv_d2hdv2) && which=="solve_d2hdv2" && is.null(B)) {
    w.ranef <- attr(sXaug,"w.ranef")
    w_R_R_v <- Matrix_times_Dvec(BLOB$R_R_v,sqrt(w.ranef)[BLOB$perm_R_v])
    if (is.null(BLOB$sortPerm_R_v)) BLOB$sortPerm_R_v <- sort.list(BLOB$perm_R_v)
    BLOB$inv_d2hdv2 <- - Matrix::chol2inv(w_R_R_v)[BLOB$sortPerm_R_v,BLOB$sortPerm_R_v] 
  }
  if ( is.null(BLOB$d2hdv2) && 
       (which =="d2hdv2" || (which=="solve_d2hdv2" && is.null(BLOB$inv_d2hdv2)) ) 
  ) {
    # don't forget that the factored matrix is not the augmented design matrix ! hence w.ranef needed here
    w.ranef <- attr(sXaug,"w.ranef")
    w_R_R_v <- Matrix_times_Dvec(BLOB$R_R_v,sqrt(w.ranef)[BLOB$perm_R_v])
    if (is.null(BLOB$sortPerm_R_v)) BLOB$sortPerm_R_v <- sort.list(BLOB$perm_R_v)
    BLOB$d2hdv2 <- - Matrix::crossprod(w_R_R_v)[BLOB$sortPerm_R_v,BLOB$sortPerm_R_v]
  } 
  ## return()'s
  if (which %in% c("solve_d2hdv2")) {
    ## somewhat involved but important for timings.  
    ## FR->FR there must be better code for CAR models, using a better decomp than R block from QR ?
    ## crossprod(R_R_v P_R_v) must give gives -d2hdv2 and then:
    if ( is.null(B) ) {
      return(BLOB$inv_d2hdv2)
    } else { 
      if (is.null(BLOB$inv_d2hdv2)) {
        return( Matrix::solve(BLOB$d2hdv2, B)) 
      } else return(BLOB$inv_d2hdv2 %*% B)
    }
  } else if (which=="hatval_Z") { ## Pdiag
    if (is.null(BLOB$hatval_Z_)) {
      if (TRUE) {
        # X[,cols] = Q R P[,cols] = Q q r p => t(Q q) given by:
        if (BLOB$u_h_cols_on_left) {
          tmp_t_Qq_scaled <- BLOB$t_Q_scaled[seq_len(attr(sXaug,"n_u_h")),]
        } else {
          tmp_t_Qq_scaled <- Matrix::solve(t(BLOB$R_R_v), t( sXaug[, seq_len(attr(sXaug,"n_u_h"))[BLOB$perm_R_v] ] ))
          ## radically faster than 
          #tmp_t_Qq_scaled <- t(Matrix::qr.Q(Matrix::qr(sXaug[, seq_len(attr(sXaug,"n_u_h")) ])))
          ## but what is needed here is a function for low-rank update of QR factorizations. 
        }
        tmp_t_Qq_scaled@x <- tmp_t_Qq_scaled@x^2
        tmp <- colSums(tmp_t_Qq_scaled)
        n_u_h <- attr(sXaug,"n_u_h")
        phipos <- n_u_h+seq_len(nrow(sXaug)-n_u_h)
        BLOB$hatval_Z_ <-  list(lev_lambda=tmp[-phipos],lev_phi=tmp[phipos])
      } else { ## F A L S E
        if (TRUE) { ## slow !
          ## two qr.qy instead of one solve... but qr.qy tends to produce large dense dgeMatrix... 
          # and... this is slow... the y is always treated as dense...
          q_Rblob <- Matrix::drop0(Matrix::qr.qy(BLOB$Rblob, diag(1, nrow = nrow(BLOB$Rblob), ncol = ncol(BLOB$Rblob))))
          # Important: <dgCMatrix> <- drop0(<dgeMatrix>) allows next line
          q_Rblob@Dim[1L] <- nrow(BLOB$blob) ## adding 0-filled rows in a way sufficient for next computation
          Qq <- Matrix::drop0(Matrix::qr.Q(BLOB$blob,complete=TRUE)) %*% q_Rblob
        } else { ## slow !
          Q_Rblob <- Matrix::drop0(Matrix::qr.Q(BLOB$Rblob,complete=FALSE))
          # Important: <dgCMatrix> <- drop0(<dgeMatrix>) allows next line
          Q_Rblob@Dim[1L] <- nrow(BLOB$blob) ## adding 0-filled rows in a way sufficient for next computation
          Qq <- Matrix::drop0(Matrix::qr.Q(BLOB$blob,complete=TRUE)) %*% Q_Rblob
        }
        BLOB$hatval_Z_ <- rowSums(Qq^2)
      }
    }
    return(BLOB$hatval_Z_)
  } else if (which=="hatval") { ## used in get_hatvalues
    if (is.null(BLOB$hatval_ZX) ) {
      tmp <- BLOB$t_Q_scaled
      tmp@x <- tmp@x^2
      BLOB$hatval_ZX <-  colSums(tmp)
    }
    return(BLOB$hatval_ZX)
  } else if (which=="R_scaled") { ## used for logdet_R_scaled_b_v, logdet_sqrt_d2hdbeta2, LevMar
    return(BLOB$R_scaled)
  } else if (which=="t_Q_scaled") { ## used in get_hatvalues for nonstandard case
    return(BLOB$t_Q_scaled)
  } else if (which %in% c("absdiag_R_v")) { ## used for logdet_R_scaled_v
    if (is.null(BLOB$absdiag_R_v)) 
      BLOB$absdiag_R_v <- abs(diag(BLOB$R_R_v)) 
    return(BLOB$absdiag_R_v)
  } else if (which=="diag_pRtRp") { ## extracted for LevMar
    if (is.null(BLOB$diag_pRtRp)) {
      tmp <- BLOB$R_scaled
      tmp@x <- tmp@x^2
      BLOB$diag_pRtRp <- colSums(tmp)[BLOB$sortPerm]
    }
    return(BLOB$diag_pRtRp)
  } else if (which=="sortPerm") { ## extracted for LevMar
    return(BLOB$sortPerm)
  } else if (which=="beta_cov") { 
    beta_cov <- Matrix::chol2inv(BLOB$R_scaled)
    beta_pos <- attr(sXaug,"n_u_h")+seq_len(attr(sXaug,"pforpv"))
    sP_beta_pos <- BLOB$sortPerm[beta_pos]
    beta_cov <- as.matrix(beta_cov[sP_beta_pos,sP_beta_pos]) * attr(sXaug,"H_global_scale")
    colnames(beta_cov) <- colnames(sXaug)[beta_pos]
    return(beta_cov)
  } else if (which=="beta_v_cov_from_wAugX") { ## using a weighted Henderson's augmented design matrix, not a true sXaug  
    beta_v_cov <- as.matrix(Matrix::chol2inv(BLOB$R_scaled)[BLOB$sortPerm,BLOB$sortPerm])
    # this tends to be dense bc v_h estimates covary (example: wafers)
    # otherwise, dropO(,tol=...), + some fix in summary.HLfit for matrix[] <- Matrix assignment, would be useful.  
    return(beta_v_cov)
  } else if (which=="d2hdv2") { ## does not seem to be used
    return(BLOB$d2hdv2)
  } 
  stop("invalid 'which' value.")
}

# trace("get_from.sXaug_Matrix_QRP_scaled",tracer=quote(if(which=="hatval_Z") debug(attr(sXaug,"sXaug_Matrix_QRP_scaled"))))
get_from_MME.sXaug_Matrix_QRP_scaled <- function(sXaug,which="",szAug=NULL,B=NULL,
                                      damping, LMrhs, ...) {
  resu <- switch(which,
                 "logdet_R_scaled_v" = {
                   absdiag_R_v <- .sXaug_Matrix_QRP_scaled(sXaug,which="absdiag_R_v")
                   sum(log(absdiag_R_v))
                 },
                 "logdet_R_scaled_b_v" = {
                   R_scaled <- .sXaug_Matrix_QRP_scaled(sXaug,which="R_scaled")
                   sum(log(abs(diag(R_scaled))))
                 },
                 "logdet_sqrt_d2hdv2" = {
                   w.ranef <- attr(sXaug,"w.ranef")
                   absdiag_R_v <- .sXaug_Matrix_QRP_scaled(sXaug,which="absdiag_R_v")
                   sum(log(sqrt(w.ranef) * absdiag_R_v))
                 },
                 # "logdet_sqrt_d2hdBV2" = {
                 #   w.ranef <- attr(sXaug,"w.ranef")
                 #   R_scaled <- attr(sXaug,"sXaug_Matrix_QRP_scaled")(sXaug,which="R_scaled")
                 #   sum(log(abs(diag(R_scaled)))) + sum(log(w.ranef))/2 - attr(sXaug,"pforpv")*log(H_global_scale)/2 
                 # },
                 "logdet_sqrt_d2hdbeta2" = {
                   H_global_scale <- attr(sXaug,"H_global_scale")
                   R_scaled <- .sXaug_Matrix_QRP_scaled(sXaug,which="R_scaled")
                   absdiag_R_v <- .sXaug_Matrix_QRP_scaled(sXaug,which="absdiag_R_v")
                   sum(log(abs(diag(R_scaled)))) - sum(log(absdiag_R_v)) - attr(sXaug,"pforpv")*log(H_global_scale)/2 
                 },
                 "LevMar_step" = {
                   ## FR->FR probably not the most elegant implementation 
                   R_scaled <- .sXaug_Matrix_QRP_scaled(sXaug,which="R_scaled")
                   sortPerm <- .sXaug_Matrix_QRP_scaled(sXaug,which="sortPerm")
                   diag_pRtRp <- .sXaug_Matrix_QRP_scaled(sXaug,which="diag_pRtRp")
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
                 ## all other cases:
                 .sXaug_Matrix_QRP_scaled(sXaug,which=which,szAug=szAug,B=B)
  )
  return(resu)
}
