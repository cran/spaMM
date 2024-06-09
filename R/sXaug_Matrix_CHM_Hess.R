# 'constructor' for sXaug_Matrix_CHM_H_scaled object
# from Xaug which already has a *scaled* ZAL 
def_sXaug_Matrix_CHM_H_scaled <- function(Xaug,weight_X,w.ranef,H_global_scale, force_QRP=.spaMM.data$options$force_LLF_CHM_QRP) {
  H_w.resid <- attr(weight_X,"H_w.resid")
  signs <- attr(H_w.resid,"signs")
  if ( is.null(signs) &&  ! .spaMM.data$options$force_LLM_nosigns_CHM_H) {
    return(def_sXaug_Matrix_QRP_CHM_scaled(Xaug,weight_X,w.ranef,H_global_scale, nonSPD=FALSE)) 
  }
  n_u_h <- length(w.ranef)
  BLOB <- list2env(list(H_w.resid=H_w.resid, # the unscaled version,
                     #weight_X=weight_X, # only for debugging
                     signs=signs), # may be NULL,
                   # not for internal use in MME methods, but for external use of attr(envir$sXaug,"BLOB")$H_w.resid
                   parent=emptyenv())
  sXaug <- .Dvec_times_Matrix_lower_block(weight_X,Xaug,n_u_h) # using signless weights
  if (is.null(signs))  {
    negHess <- .crossprod(sXaug)
  } else {
    BLOB$signed <- .Dvec_times_Matrix_lower_block(signs,sXaug,n_u_h) 
    negHess <- forceSymmetric(.crossprod(BLOB$signed,sXaug))
  }
  AUGI0_ZX_envir <- attr(sXaug,"AUGI0_ZX")$envir # immediately used, but we could imagine extending AUGI0_ZX usage
  # both to block operations here, and to other sXaug methods, which would thus all 
  # have a $AUGI0_ZX (for stable info cross sXaug) and a $BLOB.  
  if (is.null(template <- AUGI0_ZX_envir$template_CHM_negHess)) { 
    BLOB$CHMfactor <- .silent_W_E(Matrix::Cholesky(negHess,LDL=FALSE, perm=TRUE))
    if ((nonSPD <- inherits(BLOB$CHMfactor, "simpleError")) || force_QRP) {
      return(def_sXaug_Matrix_QRP_CHM_scaled(Xaug,weight_X,w.ranef,H_global_scale, nonSPD=nonSPD)) 
      ## lowest <- RSpectra::eigs(as(negHess,"dgCMatrix"), k=1, which="SR", opts=list(retvec=FALSE))$values
      # EEV <- extreme_eig(as(negHess,"dgCMatrix"), symmetric=TRUE, required=TRUE)
      # BLOB$CHMfactor <- Matrix::Cholesky(negHess,LDL=FALSE, perm=TRUE, Imult= (1e-10*EEV[1]-EEV[2])/(1-1e-10) )
    } else {
      BLOB$nonSPD <- FALSE
      if (all(AUGI0_ZX_envir$updateable)) AUGI0_ZX_envir$template_CHM_negHess <- BLOB$CHMfactor
    }
  } else {
    BLOB$CHMfactor <- .silent_W_E(Matrix::.updateCHMfactor(template, negHess, mult=0)) 
    if ((nonSPD <- inherits(BLOB$CHMfactor, "simpleError")) || 
        force_QRP
        # if I force QRP at this stage, in LevM case with signs, the exact H is never used (current QRP_CHM method; modifying its LevM algo seems difficult).
        # Without force_QRP, LevM uses the exact H with signs as long as the H is SPD (CHM_H methods).
        # so it is presumably better not to force_QRP in LeM case. 
        # Hence force_QRP should be set to TRUE only in specific, non-LevM cases.
        # Allowing QRP in those cases was useful to avoid the (then quite costly) hatval computation by .CHM2uhatvals_by_subsetinv: 
        # this is effective when there are signs, as otherwise QRP_CHM is already used by default.
      ) {
      return(def_sXaug_Matrix_QRP_CHM_scaled(Xaug,weight_X,w.ranef,H_global_scale, nonSPD=nonSPD)) 
      ## lowest <- RSpectra::eigs(as(negHess,"dgCMatrix"), k=1, which="SR", opts=list(retvec=FALSE))$values
      # EEV <- extreme_eig(as(negHess,"dgCMatrix"), symmetric=TRUE, required=TRUE)
      # BLOB$CHMfactor <- Matrix::.updateCHMfactor(template, negHess, mult= (1e-10*EEV[1]-EEV[2])/(1-1e-10))
    } else BLOB$nonSPD <- FALSE
  }
  BLOB$signs_in_WLS_mat <- ( ! is.null(BLOB$signs)) # sufficient condition HERE given we are in SPD case
  BLOB$WLS_mat_weights <- H_w.resid # nonSPD always FALSE at this point in code
  resu <- structure(sXaug,
                    # AUGI0_ZX attr already present
                    BLOB=BLOB,
                    get_from="get_from_MME.sXaug_Matrix_CHM_H_scaled",
                    w.ranef=w.ranef,
                    n_u_h=n_u_h,
                    pforpv=ncol(Xaug)-n_u_h,
                    H_global_scale=H_global_scale
  )
  ## cannot modify the class attribute... => immediate clumsy code below...
  return( resu ) 
}

.CHM2uhatval <- function(BLOB, IZ) {
  # Even if we used CHM_H (=> LevM-signs-SPD), this qr-based code is much faster than using subsetinv.
  qrX <- qr(IZ) 
  R_scaled <- qrR(qrX,backPermute = FALSE)
  t_Q <- crossprod(solve(R_scaled),t(IZ[,qrX@q+1]))
  if (is.null(BLOB$signs)) { # if I (stupidly) avoid QRP in that case...
    x <- t_Q@x
    t_Q@x <- x*x
    return(colSums(t_Q))
  } else { # this case makes more sense (LevM-SPD)
    # but so: use QRP as 'main' facto in all cases and a special Cholesky facto within QRP method in (LevM) sign-SPD case
    # => this is not run in test-devel_LLM at least.
    invIm2QtdQ_ZX <- .Rcpp_adhoc_shermanM_sp(t_Q, c(rep(0L,ncol(IZ)),BLOB$signs<0))
    colSums(( invIm2QtdQ_ZX %*% t_Q) * t_Q)
  }  
}

.calc_hatval_Z_by_subsetinv <- function(BLOB, sXaug, n_u_h) {
  
  # tmp <- colSums(t(sXaug[,seq_n_u_h]) * tcrossprod(subsetinv, sXaug[,seq_n_u_h]))
  # faster:
  # tmp <- rowSums(tcrossprod(sXaug[,seq_n_u_h],subsetinv) * sXaug[,seq_n_u_h])
  # tmp <- rowSums(.sparse_cwiseprod(sXaug[,seq_n_u_h], 
  #                                  (sXaug[,seq_n_u_h] %*% subsetinv)))
  seq_n_u_h <- seq(n_u_h)
  lev_Z <- .CHM2uhatval(BLOB=BLOB, IZ=sXaug[, seq_n_u_h]) # \bG corresponds to the true Hessian and sXaug to the absolute weights,
  # and the BLOB$signs are not factored in .CHM2uhatval), so its return value corresponds to the QQ' term in my doc, 
  # and the signs must be added to get the leverages of the true Hessian:
  hatval_Z_ <- list()
  if  ( ! is.null(BLOB$signs)) {
    hatval_Z_$lev_phi <- lev_Z[-seq_n_u_h]  * BLOB$signs
  } else hatval_Z_$lev_phi <- lev_Z[-seq_n_u_h]
  # if ("lambda" %in% B) 
  hatval_Z_$lev_lambda <- lev_Z[seq_n_u_h]
  # if ("phi" %in% B) 
  hatval_Z_
  
}

# if factorization not previously available, szAug=NULL: constructs the factorization, return value depends on 'which'
# if factorization not previously available, szAug ! NULL:   ..........................   returns the solution
# if which="solve_d2hdv2", B=NULL: returns solve(d2hdv2)
# if which="solve_d2hdv2", B ! NULL: returns solve(d2hdv2,B)
.sXaug_Matrix_CHM_H_scaled <- function(sXaug,which="",szAug=NULL,B=NULL) {
  BLOB <- attr(sXaug,"BLOB") ## an environment
  # Peviously there was an optional attr(sXaug,"AUGI0_ZX"), 
  # presumably inherited from the Xaug argument of def_sXaug_Matrix_QRP_CHM_scaled(), itself from 
  # (Xaug=Xscal with attr(Xscal,"AUGI0_ZX") <- processed$AUGI0_ZX in .solve_IRLS_as_ZX(),
  # but this attribute was not necessarily present (in .solve_v_h_IRLS)
  # => we could not generally assume it is present. 
  # In version 3.9.64 I tentatively fixed this by adding
  #     attr(Xscal,"AUGI0_ZX") <- AUGI0_ZX 
  # in .make_Xscal and in a more ad-hoc function.
  # By contrast .sXaug_Matrix_QRP_CHM_scaled() worked bc it did not try to use AUGI0_ZX except as provided by preprocess for aug_ZXy's .get_absdiagR_blocks().
  if (is.null(BLOB$perm)) {
    n_u_h <- attr(sXaug,"n_u_h")
    #
    # : stores the whole CHMfactor object which has useful methods including solve()
    pMat <- as(BLOB$CHMfactor,"pMatrix")
    if (.hasSlot(pMat, "margin") && pMat@margin==2) {
      BLOB$perm <- sort.list(pMat@perm)
    } else BLOB$perm <- pMat@perm 
    delayedAssign("L_scaled", as(BLOB$CHMfactor,"CsparseMatrix"), assign.env = BLOB )
    delayedAssign("invsqrtwranef", 1/sqrt(attr(sXaug,"w.ranef")), assign.env = BLOB )
    delayedAssign("sortPerm", sort.list(BLOB$perm), assign.env = BLOB ) # never NULL
    delayedAssign("u_h_cols_on_left", (max(BLOB$sortPerm[seq_len(n_u_h)])==n_u_h), assign.env = BLOB ) # but currently a single use
    # ! L to be used a tcrossfactor of permuted Hessian, in contract to the crossfact provided by chol() !
    delayedAssign("hatval_Z_by_subsetinv",  .calc_hatval_Z_by_subsetinv(BLOB, sXaug, n_u_h), assign.env = BLOB )
    delayedAssign("t_Q_scaled", {
      warning("Inefficient code in .sXaug_Matrix_CHM_H_scaled", immediate. = TRUE)
      Matrix::drop0(Matrix::solve(BLOB$CHMfactor,t(sXaug[,BLOB$perm]),system="L")) # quite slow...
    } , assign.env = BLOB ) 
    delayedAssign("CHMfactor_wd2hdv2w", {
      seq_n_u_h <- seq_len(n_u_h)
      wd2hdv2w <- Matrix::tcrossprod(BLOB$L_scaled[BLOB$sortPerm[ seq_n_u_h ],] ) # L_scaled is tcrossfac; CHMfactor too
      Cholesky(wd2hdv2w,LDL=FALSE, perm=FALSE ) ##  _F I X M E__?  see comments on updateable etc in QRP_CHM version
    }, assign.env = BLOB )
    delayedAssign("inv_d2hdv2", {        
      if (is.null(BLOB$inv_factor_wd2hdv2w)) {
        inv_d2hdv2 <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w, Diagonal(x=BLOB$invsqrtwranef), system="A") 
      } else inv_d2hdv2 <- .Matrix_times_Dvec(.crossprod(BLOB$inv_factor_wd2hdv2w), BLOB$invsqrtwranef) 
      inv_d2hdv2 <- .Dvec_times_Matrix( - BLOB$invsqrtwranef,inv_d2hdv2)
    }, assign.env = BLOB )
    delayedAssign("logdet_R_scaled_b_v", Matrix::determinant(BLOB$CHMfactor, sqrt=TRUE)$modulus[1], assign.env = BLOB )
    delayedAssign("logdet_R_scaled_v", { Matrix::determinant(BLOB$CHMfactor_wd2hdv2w, sqrt=TRUE)$modulus[1] }, assign.env = BLOB ) 
    delayedAssign("logdet_sqrt_d2hdv2", { sum(log(attr(sXaug,"w.ranef")))/2 + BLOB$logdet_R_scaled_v }, assign.env = BLOB )
    delayedAssign("logdet_r22", { BLOB$logdet_R_scaled_b_v - BLOB$logdet_R_scaled_v - attr(sXaug,"pforpv")*log(attr(sXaug,"H_global_scale"))/2 }, assign.env = BLOB )
    delayedAssign("hatval", { # colSums(t_Q_scaled@x^2)
      tmp <- BLOB$t_Q_scaled
      xx <- tmp@x
      xx <- xx*xx
      tmp@x <- xx
      if  (! is.null(BLOB$signs)) {
        resu <- colSums(tmp)
        n_u_h <- attr(sXaug,"n_u_h")
        phipos <- n_u_h+seq_len(nrow(sXaug)-n_u_h)
        resu[phipos] <- resu[phipos]*BLOB$signs
        resu
      } else colSums(tmp) 
    }, assign.env = BLOB )
    delayedAssign("hatval_Z_u_h_cols_on_left", {
      n_u_h <- attr(sXaug,"n_u_h")
      tmp_t_Qq_scaled <- BLOB$t_Q_scaled[seq_len(n_u_h),]
      xx <- tmp_t_Qq_scaled@x
      xx <- xx*xx
      tmp_t_Qq_scaled@x <- xx
      tmp <- colSums(tmp_t_Qq_scaled)
      phipos <- n_u_h+seq_len(nrow(sXaug)-n_u_h)
      if  (! is.null(BLOB$signs)) {
        hatval_Z_ <-  list(lev_lambda=tmp[-phipos],lev_phi=tmp[phipos]*BLOB$signs) # I checked this by inactivating the BLOB$hatval_Z_by_subsetinv computation
        #   and comparing the results to both those of spprec and those of BLOB$hatval_Z_by_subsetinv.
      } else hatval_Z_ <-  list(lev_lambda=tmp[-phipos],lev_phi=tmp[phipos])
    }, assign.env = BLOB )
    delayedAssign("Z_lev_lambda", {      
      if (is.null(t_Qq_lam_cols <- BLOB$inv_factor_wd2hdv2w)) {
        t_Qq_lam_cols <- BLOB$inv_factor_wd2hdv2w <- solve(BLOB$CHMfactor_wd2hdv2w,system="L", b=attr(sXaug,"AUGI0_ZX")$I)
      }  
      xx <- t_Qq_lam_cols@x
      xx <- xx*xx
      t_Qq_lam_cols@x <- xx
      lev_lambda <- colSums(t_Qq_lam_cols)
    }, assign.env = BLOB )
    delayedAssign("Z_lev_phi", {      
      n_u_h <- attr(sXaug,"n_u_h")
      phipos <- (n_u_h+1L):nrow(sXaug)                 # ZAL block:
      if (is.null(BLOB$inv_factor_wd2hdv2w)) {
        BLOB$inv_factor_wd2hdv2w <- solve(BLOB$CHMfactor_wd2hdv2w,system="L", b=attr(sXaug,"AUGI0_ZX")$I)
      } 
      # there is no $t_Qq but the names here show the parallel with other matrix methods
      t_Qq_phi_cols <- .tcrossprod(BLOB$inv_factor_wd2hdv2w, sXaug[phipos, seq_len(n_u_h) ], chk_sparse2mat = FALSE) # Matrix::solve(BLOB$CHMfactor_wd2hdv2w, t(sXaug[phipos, seq_len(n_u_h) ]),system="L") ## fixme the t() may still be costly
      # :if QRmethod is forced to "sparse" on mathematically dense matrices, we may reach this code yielding a *m* atrix lev_phi unless chk_sparse2mat = FALSE
      xx <- t_Qq_phi_cols@x
      xx <- xx*xx
      t_Qq_phi_cols@x <- xx
      lev_phi <- colSums(t_Qq_phi_cols)
      if  (! is.null(BLOB$signs)) lev_phi <- lev_phi*BLOB$signs
    }, assign.env = BLOB )
  } 
  if ( ! is.null(szAug)) {
    if (.is_evaluated("t_Q_scaled", BLOB)) {
      if  (! is.null(BLOB$signs)) {
        augsigns <- c(rep(1,attr(sXaug,"n_u_h")), BLOB$signs) 
        coefs <- Matrix::solve(BLOB$CHMfactor, 
                               (BLOB$t_Q_scaled %*% (augsigns * szAug)),system="Lt")[BLOB$sortPerm,,drop=FALSE] 
      } else coefs <- Matrix::solve(BLOB$CHMfactor, 
                                            (BLOB$t_Q_scaled %*% szAug),system="Lt")[BLOB$sortPerm,,drop=FALSE] 
    } else if ( ! is.null(BLOB$signs)) {
      ## many(solve(augmented linear equations) in GLMM, avoid t_Q_computation there)
      coefs <- Matrix::solve(BLOB$CHMfactor,crossprod(BLOB$signed, szAug))
    } else coefs <- Matrix::solve(BLOB$CHMfactor,crossprod(sXaug, szAug))
    return(drop(coefs))
    ## valid alternative but requesting the computation of t(Q):
    ###### #return(solve(BLOB$R_scaled[,BLOB$sortPerm], BLOB$t_Q_scaled %*% szAug)) ## sort the vector instead of the matrix?
  }
  # ELSE
  if (which=="Qt_leftcols*B") {
    return(Matrix::solve(BLOB$CHMfactor,.crossprod(sXaug[-seq(attr(sXaug,"n_u_h")),BLOB$perm], B),system="L"))
  } 
  if (which=="Mg_solve_g") {
    seq_n_u_h <- seq_len(attr(sXaug,"n_u_h"))
    rhs <- B
    rhs[seq_n_u_h] <- BLOB$invsqrtwranef * rhs[seq_n_u_h]
    #essai <- solve(BLOB$CHMfactor, rhs, system="A")
    #rhs %*% essai
    rhs <- solve(BLOB$CHMfactor, rhs[BLOB$perm], system="L")
    return(sum(rhs^2))
  } 
  if (which=="Mg_invH_g") { # - sum(B %*% inv_d2hdv2 %*% B) # was oddly inconsistent with the Matrix_QRP_CHM code until change in v3.6.4
    rhs <- BLOB$invsqrtwranef * B
    if (is.null(BLOB$inv_factor_wd2hdv2w)) { ## can be assigned elsewhere by lev_lambda <- BLOB$inv_factor_wd2hdv2w <- ....
      rhs <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w,rhs,system="L")
    } else rhs <- BLOB$inv_factor_wd2hdv2w %*% rhs
    return(sum(rhs^2)) 
  } 
  if (which=="Mg_invXtWX_g") { ## for which_LevMar_step="b", not currently used
    if (is.null(BLOB$XtWX)) BLOB$XtWX <- .crossprod(sXaug[-seq_n_u_h,-seq_n_u_h])
    Mg_invXtWX_g <- crossprod(B,solve(BLOB$XtWX,B))[1L] # [1L] drops possible Matrix class...
    return(Mg_invXtWX_g)
  } 
  if (which=="hatval_Z") { ## hat values ML or calc_sscaled_new -> Pdiag
    
    return(BLOB$hatval_Z_by_subsetinv)
    
    # older agos requiring slow step: computation of t_Q_scaled
    
    # X[,cols] = Q R P[,cols] = Q q r p => t(Q q) given by:
    #t_Qq_scaled <- solve(BLOB$CHMfactor_wd2hdv2w, ## likely bottleneck for large data 
    #                          t(sXaug[, seq_len(attr(sXaug,"n_u_h"))[BLOB$perm_R_v] ]),system="L")
    ## 
    ## but t(sXaug) is scaled such that the left block is an identity matrix, so we can work 
    ## on two separate blocks if the Cholesky is not permuted. Then 
    if (BLOB$u_h_cols_on_left) { # Rasch... bigranefs... 
      return(BLOB$hatval_Z_u_h_cols_on_left) # no clear benefits in separating lev_phi and lev_lambda computations
    } else {
      hatval_Z_ <- list()
      if ("lambda" %in% B) hatval_Z_$lev_lambda <- BLOB$Z_lev_lambda
      if ("phi" %in% B) hatval_Z_$lev_phi <- BLOB$Z_lev_phi
      return(hatval_Z_)
    }
  }
  if (which=="solve_d2hdv2") {
    # don't forget that the factored matrix is not the augmented design matrix ! hence w.ranef needed here
    w.ranef <- attr(sXaug,"w.ranef")
    if (is.null(B)) {
      return(BLOB$inv_d2hdv2)
    } else { ## solve (Matrix,vector)
      if (.is_evaluated("inv_d2hdv2", BLOB)) {
        return(BLOB$inv_d2hdv2 %*% B)
      } else {        
        not_vector <- (( ! is.null(dimB <- dim(B))) && length(dimB)==2L && dimB[2L]>1L) ## more canonical method ?
        if (not_vector) {
          rhs <- .Dvec_times_m_Matrix(BLOB$invsqrtwranef,B)
        } else rhs <- BLOB$invsqrtwranef * B
        if (is.null(BLOB$inv_factor_wd2hdv2w)) {
          rhs <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w,rhs,system="A") ## dge (if rhs is dense, or a vector), or dgC...
        } else rhs <- .crossprod(BLOB$inv_factor_wd2hdv2w, drop(BLOB$inv_factor_wd2hdv2w %*% rhs)) # typical case when solve_d2hdv2 follows hatval_Z in .calc_sscaled_new()
        if (not_vector) {
          rhs <- .Dvec_times_m_Matrix( - BLOB$invsqrtwranef,rhs)
        } else rhs <- - BLOB$invsqrtwranef * rhs
        return( rhs) # note the minus sign on the vector, - BLOB$invsqrtwranef, rather than the final, possibly matrix, rhs
      }
    } 
  } 
  if (which=="hatval") { return(BLOB$hatval) } # REML hatval computation (also named hatval_ZX)
  if (which=="R_scaled_blob") { ## used for LevMar
    if (is.null(BLOB$R_scaled_blob)) {
      tmp <- X <- t(BLOB$L_scaled[BLOB$sortPerm,, drop=FALSE] ) # crossfac needed here
      xx <- tmp@x
      xx <- xx*xx
      tmp@x <- xx
      BLOB$R_scaled_blob <- list(X=X, diag_pRtRp=colSums(tmp), 
                                 XDtemplate=.XDtemplate(X,upperTri=FALSE))
    }
    return(BLOB$R_scaled_blob)
  } 
  if (which=="R_scaled_v_h_blob") {
    if (is.null(BLOB$R_scaled_v_h_blob)) {
      R_scaled_v_h <- t( as(BLOB$CHMfactor_wd2hdv2w,"CsparseMatrix") ) ## the t() for .damping_to_solve... (fixme: if we could avoid t()...)
      tmp <- R_scaled_v_h 
      xx <- tmp@x
      xx <- xx*xx
      tmp@x <- xx
      diag_pRtRp_scaled_v_h <- colSums(tmp)
      BLOB$R_scaled_v_h_blob <- list(R_scaled_v_h=R_scaled_v_h,diag_pRtRp_scaled_v_h=diag_pRtRp_scaled_v_h, 
                                     XDtemplate=.XDtemplate(R_scaled_v_h, upperTri=TRUE))
    }
    return(BLOB$R_scaled_v_h_blob)
  } 
  if (which=="R_beta_blob") {
    if (is.null(BLOB$R_beta_blob)) {
      n_u_h <- attr(sXaug,"n_u_h")
      seq_n_u_h <- seq(n_u_h)
      X <- as.matrix(sXaug[-seq_n_u_h,-seq_n_u_h]) ## The following code assumes it is dense...
      R_beta <- .lmwithQR(X,yy=NULL,returntQ=FALSE,returnR=TRUE)$R_scaled
      diag_pRtRp_beta <-  colSums(R_beta^2)
      BLOB$R_beta_blob <- list(R_beta=R_beta,diag_pRtRp_beta=diag_pRtRp_beta, 
                               XDtemplate=.XDtemplate(R_beta, upperTri=TRUE))
    }
    return(BLOB$R_beta_blob)
  } 
  # if (which=="sortPerm") {
  #   return(BLOB$sortPerm)
  # } 
  if (which=="t_Q_scaled") {
    return(BLOB$t_Q_scaled)
  } 
  if (which=="logdet_R_scaled_b_v") { return(BLOB$logdet_R_scaled_b_v) } 
  if (which=="logdet_sqrt_d2hdv2") { return(BLOB$logdet_sqrt_d2hdv2)} 
  if (which=="logdet_r22") { return(BLOB$logdet_r22) }
  if (which=="beta_cov_info_from_sXaug") { 
    return(.calc_beta_cov_info_from_sXaug(BLOB=BLOB, sXaug=sXaug, 
                                          tcrossfac=solve(BLOB$CHMfactor,system="Lt", b=.sparseDiagonal(n=ncol(BLOB$CHMfactor), shape="g"))))
  } 
  if (which=="beta_cov_info_from_wAugX") { ## using a weighted Henderson's augmented design matrix, not a true sXaug  
    if (TRUE) {
      tcrossfac_beta_v_cov <- solve(BLOB$CHMfactor,system="Lt", b=.sparseDiagonal(n=ncol(BLOB$CHMfactor), shape="g"))
      tPmat <- sparseMatrix(seq_along(BLOB$sortPerm), BLOB$sortPerm, x=1)
      tcrossfac_beta_v_cov <- as.matrix(tPmat %*% tcrossfac_beta_v_cov)
      rownames(tcrossfac_beta_v_cov) <- colnames(sXaug) ## necessary for summary.HLfit, already lost in BLOB$R_scaled
      pforpv <- attr(sXaug,"pforpv")
      seqp <- seq_len(pforpv)
      beta_cov <- .tcrossprod(tcrossfac_beta_v_cov[seqp,,drop=FALSE])
      return( list(beta_cov=beta_cov, 
                   #beta_v_cov=beta_v_cov,
                   tcrossfac_beta_v_cov=tcrossfac_beta_v_cov) )
    }
    ########################
    beta_v_cov <- chol2inv(BLOB$R_scaled)
    beta_v_cov <- beta_v_cov[BLOB$sortPerm,BLOB$sortPerm,drop=FALSE]
    # this tends to be dense bc v_h estimates covary (example: wafers)
    # otherwise, dropO(,tol=...), + some fix in summary.HLfit for matrix[] <- Matrix assignment, would be useful.  
    return(beta_v_cov)
  }  
  if (which=="d2hdv2") {
    # if (is.null(BLOB$d2hdv2)) {
    #     # don't forgetthat the factored matrix is not the augmented design matrix ! hence w.ranef needed here
    #     w.ranef <- attr(sXaug,"w.ranef")
    #     w_R_R_v <- .Matrix_times_Dvec(BLOB$R_R_v,sqrt(w.ranef)[BLOB$perm_R_v])
    #     if (is.null(BLOB$sortPerm_R_v)) BLOB$sortPerm_R_v <- sort.list(BLOB$perm_R_v)
    #     BLOB$d2hdv2 <- - Matrix::crossprod(w_R_R_v)[BLOB$sortPerm_R_v,BLOB$sortPerm_R_v]
    #   }
    #   return(BLOB$d2hdv2)
    stop("d2hdv2 requested")
  } 
  stop("invalid 'which' value.")
} 

# trace("get_from_MME.sXaug_Matrix_CHM_H_scaled",tracer=quote(if(which=="hatval_Z") debug(.sXaug_Matrix_CHM_H_scaled)))
# trace("get_from_MME.sXaug_Matrix_CHM_H_scaled",exit=quote(str(resu)))
# trace("get_from_MME.sXaug_Matrix_CHM_H_scaled",exit=quote({cat(which); str(resu)}))
get_from_MME.sXaug_Matrix_CHM_H_scaled <- function(sXaug,which="",szAug=NULL,B=NULL,
                                        damping, LMrhs, ...) {
  resu <- switch(which,
                 "LevMar_step" = { # test-cloglog...
                   R_scaled_blob <- .sXaug_Matrix_CHM_H_scaled(sXaug,which="R_scaled_blob") 
                   dampDpD <- damping*R_scaled_blob$diag_pRtRp ## NocedalW p. 266
                   # Extend the X in X'X = P'R'RP:
                   list(dVscaled_beta=.damping_to_solve(XDtemplate=R_scaled_blob$XDtemplate, dampDpD=dampDpD, rhs=LMrhs),
                                rhs=LMrhs, dampDpD = dampDpD) 
                 },
                 "LevMar_step_v_h" = {
                   ## FR->FR probably not the most elegant implementation 
                   R_scaled_v_h_blob <- .sXaug_Matrix_CHM_H_scaled(sXaug,which="R_scaled_v_h_blob") 
                   dampDpD <- damping*R_scaled_v_h_blob$diag_pRtRp_scaled_v_h ## NocedalW p. 266
                   # Extend the X in X'X = P'R'RP: 
                   list(dVscaled = .damping_to_solve(XDtemplate=R_scaled_v_h_blob$XDtemplate, dampDpD=dampDpD, rhs=LMrhs), 
                        dampDpD = dampDpD) 
                 },
                 "LevMar_step_beta" = {
                   if ( ! length(LMrhs)) stop("LevMar_step_beta called with 0-length LMrhs: pforpv=0?")
                   R_beta_blob <- .sXaug_Matrix_CHM_H_scaled(sXaug,which="R_beta_blob")
                   dampDpD <- damping*R_beta_blob$diag_pRtRp_beta ## NocedalW p. 266
                   list(dbeta = .damping_to_solve(XDtemplate=R_beta_blob$XDtemplate, dampDpD=dampDpD, rhs=LMrhs), 
                        dampDpD = dampDpD) 
                 } ,
                 ## all other cases:
                 .sXaug_Matrix_CHM_H_scaled(sXaug,which=which,szAug=szAug,B=B)
  )
  return(resu)
}
