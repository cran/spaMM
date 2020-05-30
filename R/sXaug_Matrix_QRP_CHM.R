# 'constructor' 
# from Xaug which already has a *scaled* ZAL 
def_sXaug_Matrix_QRP_CHM_scaled <- function(Xaug,weight_X,w.ranef,H_global_scale) {
  n_u_h <- length(w.ranef)
  Xrows <- n_u_h+seq(length(weight_X))
  ## Bates https://stat.ethz.ch/pipermail/r-help/2010-December/262365.html
  ## "Assignment of submatrices in a sparse matrix can be slow because there is so much checking that needs to be done."
  Xaug <- .Dvec_times_Matrix_lower_block(weight_X,Xaug,n_u_h)
  resu <- structure(Xaug,
                    get_from="get_from_MME.sXaug_Matrix_QRP_CHM_scaled",
                    BLOB=list2env(list(), parent=environment(.sXaug_Matrix_QRP_CHM_scaled)),
                    w.ranef=w.ranef,
                    n_u_h=n_u_h, # mandatory for all sXaug types
                    pforpv=ncol(Xaug)-n_u_h,  # mandatory for all sXaug types
                    weight_X=weight_X, # new mandatory 08/2018 # not used by the "methods", 
                                       # but by .calc_APHLs_by_augZXy_or_sXaug() -> attributes(sXaug) (wherein it may itself no longer be used ?)
                    H_global_scale=H_global_scale
  )
  ## cannot modify the 'class' attribute... => immediate clumsy code below... and how(.) reports method: dgCMatrix.
  return( resu ) 
}

.sXaug_Matrix_QRP_CHM_scaled <- function(sXaug,which="",szAug=NULL,B=NULL) {
  BLOB <- attr(sXaug,"BLOB") ## an environment
  if (is.null(BLOB$blob)) {
    BLOB$blob <- qr(sXaug) ##  Matrix::qr
    # sXaug = t(tQ) %*% R[,sP] but then also sXaug = t(tQ)[,sP'] %*% R[sP',sP] for any sP'
    BLOB$perm <- BLOB$blob@q + 1L
    BLOB$R_scaled <- qrR(BLOB$blob,backPermute = FALSE)
    n_u_h <- attr(sXaug,"n_u_h")
    seq_n_u_h <- seq(n_u_h)
    BLOB$sortPerm <- sort.list(BLOB$perm) 
    BLOB$u_h_cols_on_left <- (max(BLOB$sortPerm[ seq_n_u_h ])==n_u_h) # e.g., Rasch
    if ( ! is.null(szAug)) return(qr.coef(BLOB$blob,szAug))   
  } 
  if ( is.null(BLOB$t_Q_scaled) && 
       ( which %in% c("t_Q_scaled","hatval") || (which=="hatval_Z" && BLOB$u_h_cols_on_left)
       ) ) {
    if (FALSE) { # after efforts to find an Eigen syntax that works without storing a dense matrix in memory 
      # (solveInPlace may be the only way), we see that it's slow! (e.g test bigranefs)
      # https://stackoverflow.com/questions/53290521/eigen-solver-for-sparse-matrix-product
      BLOB$t_Q_scaled <- .Rcpp_backsolve_M_M(r=as(BLOB$R_scaled,"dgCMatrix"),x=t(sXaug[,BLOB$perm]),transpose=TRUE)
      # base::backsolve() calls as.matrix()...
    } else BLOB$t_Q_scaled <- solve(t(BLOB$R_scaled),t(sXaug[,BLOB$perm])) ## Matrix::solve ## this is faster than qr.Q and returns dgCMatrix rather than qr.Q -> dge !
  }
  if ( ! is.null(szAug)) {
    if (is.null(BLOB$t_Q_scaled)) {
      return(qr.coef(BLOB$blob,szAug)) # avoid t_Q_computation there ## Matrix::qr.coef
    } else {
      return(solve(BLOB$R_scaled, BLOB$t_Q_scaled %*% szAug)[BLOB$sortPerm,,drop=FALSE]) ## Matrix::solve
    }
  }
  # ELSE
  if ( is.null(BLOB$CHMfactor_wd2hdv2w)  && 
       (which %in% c("logdet_R_scaled_v","hatval_Z","d2hdv2","solve_d2hdv2","R_scaled_v_h_blob", "Mg_invH_g"))
  ) {
    seq_n_u_h <- seq_len(attr(sXaug,"n_u_h"))
    wd2hdv2w <- .crossprod(BLOB$R_scaled[,BLOB$sortPerm[ seq_n_u_h ]], allow_as_mat = FALSE ) 
    BLOB$CHMfactor_wd2hdv2w <- Cholesky(wd2hdv2w,LDL=FALSE,
                                                perm=FALSE ) ## perm=FALSE useful for leverage computation as explained below
  }
  ## return()'s
  if (which=="Mg_solve_g") {
    if (is.null(BLOB$invsqrtwranef)) BLOB$invsqrtwranef <- 1/sqrt(attr(sXaug,"w.ranef")) 
    seq_n_u_h <- seq_len(attr(sXaug,"n_u_h"))
    rhs <- B
    rhs[seq_n_u_h] <- BLOB$invsqrtwranef * rhs[seq_n_u_h]
    rhs <- Matrix::solve(BLOB$R_scaled,rhs[BLOB$perm])
    return(sum(rhs^2))
  } 
  if (which=="Mg_invH_g") {
    if (is.null(BLOB$invsqrtwranef)) BLOB$invsqrtwranef <- 1/sqrt(attr(sXaug,"w.ranef")) 
    rhs <- BLOB$invsqrtwranef * B
    if (is.null(BLOB$inv_factor_wd2hdv2w)) { ## can be assigned elsewhere by lev_lambda <- BLOB$inv_factor_wd2hdv2w <- ....
      rhs <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w,rhs,system="L")
    } else rhs <- BLOB$inv_factor_wd2hdv2w %*% rhs
    return(sum(rhs^2))
  } 
  if (which=="Mg_invXtWX_g") { ## 
    if (is.null(BLOB$XtWX)) BLOB$XtWX <- sXaug[-seq_n_u_h,-seq_n_u_h]
    Mg_invXtWX_g <- crossprod(B,solve(BLOB$XtWX,B))[1L] # [1L] drops possible Matrix class...
    return(Mg_invXtWX_g)
  }
  if (which %in% c("solve_d2hdv2")) {
    # don't forget that the factored matrix is not the augmented design matrix ! hence w.ranef needed here
    #if (is.null(BLOB$sortPerm_R_v)) BLOB$sortPerm_R_v <- sort.list(BLOB$perm_R_v)
    if (is.null(BLOB$invsqrtwranef)) BLOB$invsqrtwranef <- 1/sqrt(attr(sXaug,"w.ranef")) 
    if (is.null(B)) {
      if (is.null(BLOB$inv_d2hdv2)) {
        # if (is.null(BLOB$sortPerm_R_v)) BLOB$sortPerm_R_v <- sort.list(BLOB$perm_R_v)
        # w_R_R_v <- .Matrix_times_Dvec(BLOB$R_R_v,sqrt(w.ranef)[BLOB$perm_R_v])
        # BLOB$inv_d2hdv2 <- - Matrix::chol2inv(w_R_R_v)[BLOB$sortPerm_R_v,BLOB$sortPerm_R_v] 
        if (is.null(BLOB$inv_factor_wd2hdv2w)) {
          inv_d2hdv2 <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w, Diagonal(x=BLOB$invsqrtwranef)) 
        } else inv_d2hdv2 <- .Matrix_times_Dvec(BLOB$inv_factor_wd2hdv2w, BLOB$invsqrtwranef)
        BLOB$inv_d2hdv2 <- .Dvec_times_Matrix( - BLOB$invsqrtwranef,inv_d2hdv2)
      }
      return(BLOB$inv_d2hdv2)
    } else { ## solve (Matrix,vector)
      if (is.null(BLOB$inv_d2hdv2)) {
        not_vector <- (( ! is.null(dimB <- dim(B))) && length(dimB)==2L && dimB[2L]>1L) ## more canonical method ?
        if (not_vector) {
          rhs <- .Dvec_times_m_Matrix(BLOB$invsqrtwranef,B)
        } else rhs <- BLOB$invsqrtwranef * B
        rhs <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w,rhs,system="A") ## dge (if rhs is dense, or a vector), or dgC...
        # has dim[2] in all cases
        if (not_vector) { ## is.matrix(rhs) is not the correct test 
          rhs <- .Dvec_times_m_Matrix(BLOB$invsqrtwranef,rhs)
        } else rhs <- BLOB$invsqrtwranef * rhs
        return( - rhs)
      } else return(BLOB$inv_d2hdv2 %*% B)
    } 
  } else if (which=="hatval_Z") { ## Pdiag
    if (is.null(BLOB$hatval_Z_)) {
      # X[,cols] = Q R P[,cols] = Q q r p => t(Q q) given by:
      if (BLOB$u_h_cols_on_left) { # Rasch...
        tmp_t_Qq_scaled <- BLOB$t_Q_scaled[seq_len(attr(sXaug,"n_u_h")),]
        tmp_t_Qq_scaled@x <- tmp_t_Qq_scaled@x^2
        tmp <- colSums(tmp_t_Qq_scaled)
        n_u_h <- attr(sXaug,"n_u_h")
        phipos <- n_u_h+seq_len(nrow(sXaug)-n_u_h)
        BLOB$hatval_Z_ <-  list(lev_lambda=tmp[-phipos],lev_phi=tmp[phipos])
        # previously to v2.1.46 there were comments showing slow code using qr.qy 
      } else {
        # X[,cols] = Q R P[,cols] = Q q r p => t(Q q) given by:
        #t_Qq_scaled <- solve(BLOB$CHMfactor_wd2hdv2w, ## likely bottleneck for large data 
        #                          t(sXaug[, seq_len(attr(sXaug,"n_u_h"))[BLOB$perm_R_v] ]),system="L")
        ## 
        ## but t(sXaug) is scaled such that the left block is an identity matrix, so we can work 
        ## on two separate blocks if the Cholesky is not permuted. Then 
        if (is.null(lev_lambda <- BLOB$inv_factor_wd2hdv2w)) {
          lev_lambda <- BLOB$inv_factor_wd2hdv2w <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w,system="L")
        }  
        lev_lambda@x <- lev_lambda@x^2
        lev_lambda <- colSums(lev_lambda)
        n_u_h <- attr(sXaug,"n_u_h")
        phipos <- (n_u_h+1):nrow(sXaug)
        lev_phi <- .tcrossprod(BLOB$inv_factor_wd2hdv2w, sXaug[phipos, seq_len(n_u_h) ]) # Matrix::solve(BLOB$CHMfactor_wd2hdv2w, t(sXaug[phipos, seq_len(n_u_h) ]),system="L") ## fixme the t() may still be costly
        lev_phi@x <- lev_phi@x^2
        lev_phi <- colSums(lev_phi)
        BLOB$hatval_Z_ <-  list(lev_lambda=lev_lambda,lev_phi=lev_phi)
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
  } else if (which=="R_scaled_blob") { ## used for LevMar
    if (is.null(BLOB$R_scaled_blob)) {
      tmp <- X <- BLOB$R_scaled[ , BLOB$sortPerm]
      tmp@x <- tmp@x^2
      BLOB$R_scaled_blob <- list(X=X, diag_pRtRp=colSums(tmp), XDtemplate=.XDtemplate(X))
    }
    return(BLOB$R_scaled_blob)
  } else if (which=="R_scaled_v_h_blob") {
    if (is.null(BLOB$R_scaled_v_h_blob)) {
      R_scaled_v_h <- t( as(BLOB$CHMfactor_wd2hdv2w,"sparseMatrix") ) ## the t() for .damping_to_solve... (fixme: if we could avoid t()...)
      tmp <- R_scaled_v_h 
      tmp@x <- tmp@x^2
      diag_pRtRp_scaled_v_h <- colSums(tmp)
      BLOB$R_scaled_v_h_blob <- list(R_scaled_v_h=R_scaled_v_h,diag_pRtRp_scaled_v_h=diag_pRtRp_scaled_v_h, XDtemplate=.XDtemplate(R_scaled_v_h))
    }
    return(BLOB$R_scaled_v_h_blob)
  } else if (which=="R_beta_blob") {
    if (is.null(BLOB$R_beta_blob)) {
      n_u_h <- attr(sXaug,"n_u_h")
      seq_n_u_h <- seq(n_u_h)
      X <- as.matrix(sXaug[-seq_n_u_h,-seq_n_u_h]) ## folloing code assuming it is dense...
      R_beta <- .lmwithQR(X,yy=NULL,returntQ=FALSE,returnR=TRUE)$R_scaled
      diag_pRtRp_beta <-  colSums(R_beta^2)
      BLOB$R_beta_blob <- list(R_beta=R_beta,diag_pRtRp_beta=diag_pRtRp_beta, XDtemplate=.XDtemplate(R_beta))
    }
    return(BLOB$R_beta_blob)
  # } else if (which=="sortPerm") { 
  #   return(BLOB$sortPerm)
  } else if (which=="t_Q_scaled") { ## used in get_hatvalues for nonstandard case
    return(BLOB$t_Q_scaled)
  } else if (which %in% c("logdet_R_scaled_b_v")) {
    if (is.null(BLOB$logdet_R_scaled_b_v)) BLOB$logdet_R_scaled_b_v <- sum(log(abs(diag(BLOB$R_scaled))))
    return(BLOB$logdet_R_scaled_b_v)
  } else if (which %in% c("logdet_R_scaled_v")) {
    if (is.null(BLOB$logdet_R_scaled_v)) BLOB$logdet_R_scaled_v <- Matrix::determinant(BLOB$CHMfactor_wd2hdv2w)$modulus[1]
    return(BLOB$logdet_R_scaled_v)
  } else if (which=="beta_cov_info_from_sXaug") {
    return(.calc_beta_cov_info_from_sXaug(BLOB=BLOB, sXaug=sXaug))
  } else if (which=="beta_cov_info_from_wAugX") { ## using a weighted Henderson's augmented design matrix, not a true sXaug  
    if (TRUE) {
      tcrossfac_beta_v_cov <- solve(BLOB$R_scaled)
      tPmat <- sparseMatrix(seq_along(BLOB$sortPerm), BLOB$sortPerm, x=1)
      tcrossfac_beta_v_cov <- as.matrix(tPmat %*% tcrossfac_beta_v_cov)
      rownames(tcrossfac_beta_v_cov) <- colnames(sXaug) ## necessary for summary.HLfit, already lost in BLOB$R_scaled
      #beta_v_cov <- .tcrossprod(tcrossfac_beta_v_cov)
      pforpv <- attr(sXaug,"pforpv")
      seqp <- seq_len(pforpv)
      beta_cov <- .tcrossprod(tcrossfac_beta_v_cov[seqp,,drop=FALSE])
      return( list(beta_cov=beta_cov, 
                   #beta_v_cov=beta_v_cov,
                   tcrossfac_beta_v_cov=tcrossfac_beta_v_cov) )
    }
    ####################
    beta_v_cov <- as.matrix(Matrix::chol2inv(BLOB$R_scaled)[BLOB$sortPerm,BLOB$sortPerm])
    # this tends to be dense bc v_h estimates covary   
    return(beta_v_cov)
  } else if (which=="d2hdv2") { ## does not seem to be used
    ## #return(BLOB$d2hdv2)
    stop("d2hdv2 requested ")
  } 
  stop("invalid 'which' value.")
}

# trace("get_from_MME.sXaug_Matrix_QRP_CHM_scaled",print=FALSE, tracer=quote(print(which)),exit=quote(str(resu)))
get_from_MME.sXaug_Matrix_QRP_CHM_scaled <- function(sXaug,which="",szAug=NULL,B=NULL,
                                      damping, LMrhs, ...) {
  resu <- switch(which,
                 "logdet_sqrt_d2hdv2" = {
                   w.ranef <- attr(sXaug,"w.ranef")
                   logdet_R_scaled_v <- .sXaug_Matrix_QRP_CHM_scaled(sXaug,which="logdet_R_scaled_v")
                   sum(log(w.ranef))/2 + logdet_R_scaled_v
                 },
                 "logdet_r22" = {
                   H_global_scale <- attr(sXaug,"H_global_scale")
                   logdet_R_scaled_b_v <- .sXaug_Matrix_QRP_CHM_scaled(sXaug,which="logdet_R_scaled_b_v") 
                   logdet_R_scaled_v <- .sXaug_Matrix_QRP_CHM_scaled(sXaug,which="logdet_R_scaled_v")
                   logdet_R_scaled_b_v - logdet_R_scaled_v - attr(sXaug,"pforpv")*log(H_global_scale)/2 ## '-', not '<-'
                 },
                 "LevMar_step" = {
                   ## FR->FR probably not the most elegant implementation 
                   R_scaled_blob <- .sXaug_Matrix_QRP_CHM_scaled(sXaug,which="R_scaled_blob")
                   dampDpD <- damping*R_scaled_blob$diag_pRtRp ## NocedalW p. 266
                   # Extend the X in X'X = P'R'RP:
                   list(dVscaled_beta=.damping_to_solve(XDtemplate=R_scaled_blob$XDtemplate, dampDpD=dampDpD, rhs=LMrhs), 
                        dampDpD = dampDpD) 
                 },
                 "LevMar_step_v_h" = {
                   ## FR->FR probably not the most elegant implementation 
                   R_scaled_v_h_blob <- .sXaug_Matrix_QRP_CHM_scaled(sXaug,which="R_scaled_v_h_blob")
                   dampDpD <- damping*R_scaled_v_h_blob$diag_pRtRp_scaled_v_h ## NocedalW p. 266
                   # Extend the X in X'X = P'R'RP: 
                   list(dVscaled = .damping_to_solve(XDtemplate=R_scaled_v_h_blob$XDtemplate, dampDpD=dampDpD, rhs=LMrhs), 
                        dampDpD = dampDpD) 
                 },
                 "LevMar_step_beta" = {
                   if ( ! length(LMrhs)) stop("LevMar_step_beta called with 0-length LMrhs: pforpv=0?")
                   R_beta_blob <- .sXaug_Matrix_QRP_CHM_scaled(sXaug,which="R_beta_blob")
                   dampDpD <- damping*R_beta_blob$diag_pRtRp_beta ## NocedalW p. 266
                   list(dbeta = .damping_to_solve(XDtemplate=R_beta_blob$XDtemplate, dampDpD=dampDpD, rhs=LMrhs), 
                        dampDpD = dampDpD) 
                 } ,
                 ## all other cases:
                 .sXaug_Matrix_QRP_CHM_scaled(sXaug,which=which,szAug=szAug,B=B)
  )
  return(resu)
}
