def_AUGI0_ZX_sparsePrecision <- function(AUGI0_ZX, corrPars, w.ranef, cum_n_u_h,
                                         # w.resid, 
                                         H_w.resid 
                                         #,weight_X ## currently ignored
                                         #,H_global_scale ## currently ignored
                                         ) {
  BLOB=list2env(list(H_w.resid=H_w.resid), parent=emptyenv())
  if (any(H_w.resid<0)) BLOB$signs <- sign(H_w.resid) # a NULL BLOB$signs value meaning that all signs are >0
  resu <- list(AUGI0_ZX=AUGI0_ZX, BLOB=BLOB)
  attributes(resu) <- list(
    names=names(resu),
    w.ranef=w.ranef,
    cum_n_u_h=cum_n_u_h,
    pforpv=ncol(AUGI0_ZX$X.pv),
    # w.resid= if (is.list(w.resid)) {w.resid$w_resid} else w.resid,
    corrPars=corrPars
    )
  class(resu) <- c("AUGI0_ZX_sparsePrecision","list") # (## )do not define recursively if object is an envir...)
  .init_spprec(resu)
  return( resu ) 
}

.init_spprec <- function(sXaug) {
  AUGI0_ZX <- sXaug$AUGI0_ZX 
  BLOB <- sXaug$BLOB 
  w.ranef <- attr(sXaug,"w.ranef") 
  cum_n_u_h <- attr(sXaug,"cum_n_u_h")
  ## costly, cannot be precomputed
  if (is.null(AUGI0_ZX$envir$precisionFactorList)) stop("is.null(AUGI0_ZX$envir$precisionFactorList)")
  precisionBlocks <- .reformat_Qmat_info(BLOB, AUGI0_ZX$envir, corrPars=attr(sXaug,"corrPars")) ## includes computations of BLOB$chol_Q not function of w.ranef
  for (it in seq_len(length(precisionBlocks))) { ## FIXME respecialize this block ? 
    if (AUGI0_ZX$envir$finertypes[it]=="ranCoefs") {
      ## do not change precisionBlocks[[it]] which in this case already contains the full matrix with lambda 
    } else {
      # precisionBlocks[[it]] * w.ranef[u.range] is a correct way of computing the precision block in two cases:
      # Either a precisionBlocks is that of a non-trivial corr mat of a gaussian ranef => the lambda and w.ranef are constant (even in CAR)
      # or the lambda and w.ranef are not constant but the precisionBlocks is an identity matrix.
      # but not in spprec ranCoefs cas (heterogeneous w.ranef). 
      # .Matrix_times_Dvec() is BOTH faster by directly manipulating slots AND correct in this additional case.
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      precisionBlocks[[it]] <- .Matrix_times_Dvec(precisionBlocks[[it]], 
                                                  w.ranef[u.range],
                                                  check_dsC=FALSE, # risky: use TRUE to detect when the result cannot be dsC
                                                  keep_dsC=TRUE # result MUST be dsC
      ) 
    }
  }
  if (length(precisionBlocks)>1L) {
    #precisionMatrix <- forceSymmetric(Matrix::bdiag(precisionBlocks)) # bdiag degrades dsC into dgC
    precisionMatrix <- .bdiag_dsC(precisionBlocks)
  } else precisionMatrix <- (precisionBlocks[[1L]]) # forceSymmetric removed bc input should in principle be dsC
  if (BLOB$nonSPD <- ! is.null(BLOB$signs)) {
    
    BLOB$sqrtW <- sqrt(abs(BLOB$H_w.resid))
    #
    WZ <- .Dvec_times_Matrix(BLOB$H_w.resid, AUGI0_ZX$ZAfix)
    ZtWZ <- .crossprod(WZ,AUGI0_ZX$ZAfix, allow_as_mat=FALSE, as_sym=TRUE)
    BLOB$Gmat <- .calc_Gmat(ZtWZ=ZtWZ, precisionMatrix) # depends on w.ranef and w.resid
    #
    if (no_template <- is.null(template <- AUGI0_ZX$template_G_CHM)) { ## occurs if $update_CHM is FALSE, OR first comput. of CHM, OR precision factors not yet all $updateable
      BLOB$G_CHMfactor <- suppressWarnings(try(Cholesky(BLOB$Gmat,LDL=FALSE,perm=.spaMM.data$options$perm_G),silent=TRUE))
    } else BLOB$G_CHMfactor <- suppressWarnings(try(Matrix::.updateCHMfactor(template, BLOB$Gmat, mult=0),silent=TRUE)) 
    #
    if (BLOB$nonSPD <- inherits(BLOB$G_CHMfactor, "try-error")) {
      BLOB$WLS_mat_weights <- abs(BLOB$H_w.resid)
      WZ <- .Dvec_times_Matrix(BLOB$WLS_mat_weights, AUGI0_ZX$ZAfix)
      ZtWZ <- .crossprod(WZ,AUGI0_ZX$ZAfix, allow_as_mat=FALSE, as_sym=TRUE)
      BLOB$Gmat <- .calc_Gmat(ZtWZ=ZtWZ, precisionMatrix) # depends on w.ranef and w.resid
      if (no_template) { 
        BLOB$G_CHMfactor <- Cholesky(BLOB$Gmat,LDL=FALSE,perm=.spaMM.data$options$perm_G)
      } else BLOB$G_CHMfactor <- Matrix::.updateCHMfactor(template, BLOB$Gmat, mult=0) 
    }
    #
    if (no_template && .spaMM.data$options$update_CHM && all(AUGI0_ZX$envir$updateable) ) AUGI0_ZX$template_G_CHM <- BLOB$G_CHMfactor
    
  } else {
    BLOB$sqrtW <- sqrt(BLOB$H_w.resid)
    ZtWZ <- .ZtWZwrapper(AUGI0_ZX$ZAfix,sqrtW=BLOB$sqrtW) # should return a symmetric-type matrix (dsC or possibly dsy, not ddi...)
    BLOB$Gmat <- .calc_Gmat(ZtWZ=ZtWZ, precisionMatrix) # depends on w.ranef and w.resid
    if (is.null(template <- AUGI0_ZX$template_G_CHM)) { ## occurs if $update_CHM is FALSE, OR first comput. of CHM, OR precision factors not yet all $updateable
      BLOB$G_CHMfactor <- Cholesky(BLOB$Gmat,LDL=FALSE,perm=.spaMM.data$options$perm_G) ## costly
      if (.spaMM.data$options$update_CHM && all(AUGI0_ZX$envir$updateable)) AUGI0_ZX$template_G_CHM <- BLOB$G_CHMfactor
    } else BLOB$G_CHMfactor <- Matrix::.updateCHMfactor(template, BLOB$Gmat, mult=0) ## If it fails in spde, look for too low kappa.
  }
  if (BLOB$nonSPD) {
    BLOB$signs_in_WLS_mat <- FALSE
    # BLOB$WLS_mat_weights arleady defined
  } else {
    BLOB$signs_in_WLS_mat <- ( ! is.null(BLOB$signs)) # sufficient condition HERE given we are in SPD case
    BLOB$WLS_mat_weights <- BLOB$H_w.resid
  }
  ## with perm=TRUE G=P'LL'P and the P'L (non-triangular) factor is given by solve(<G_CHM>,as(<G_CHM>,"sparseMatrix"),system="Pt")
  # (?__F I X M E__?) It's theo. possible to call .updateCHMfactor on the tcrossfac of Gmat...
  BLOB$pMat_G <- as(BLOB$G_CHMfactor, "pMatrix") # used itself or for its @perm slot, with indices from 1L
  .init_promises_spprec(sXaug)
}


.reformat_Qmat_info <- function(BLOB, # a pointer to an envir !
                                envir, # another pointer to an envir ! 
                                corrPars=NULL) { ## provide the precision factors Q and their CHMfactor from different types of input
  ## used to compute Gmat (all non trivial elements needed) & chol_Q  
  # We want chol_Q to be (dtCMatrix: Csparse triangular) so that efficient solve methods can be used.
  # but eg bdiag(dtCMatrix, Diagonal Matrix) gives a dgCMatrix. 
  # => All of chol_Q_list must be dtCMatrix so that bdiag() gives something that can be converted to dtCMatrix
  precisionBlocks <- envir$precisionFactorList ## local copy so that AUGI0_ZX$envir$precisionFactorList is not modified
  chol_Q_list <- vector("list",length(envir$precisionFactorList)) ## of L in LL' factorisation of precision factor 
  for (rd in seq_len(length(envir$precisionFactorList))) {
    if( ! is.null(envir$precisionFactorList[[rd]])) {
      if (inherits(envir$precisionFactorList[[rd]],c("Matrix","matrix"))) stop("list expected here in .reformat_Qmat_info()") 
      chol_Q_list[[rd]] <- envir$precisionFactorList[[rd]]$chol_Q
      if (envir$finertypes[rd]=="ranCoefs") {
        precisionBlocks[[rd]] <- envir$precisionFactorList[[rd]]$precmat ## full already including lambda
      } else { ## both Q and factorization provided (see .assign_geoinfo_and_LMatrices_but_ranCoefs())
        precisionBlocks[[rd]] <- envir$precisionFactorList[[rd]]$Qmat ## indep of w.ranef (lambda); should be dsC at this point
      }
    } ## ELSE NULL precisionFactorList element <-> chol_Q_list element prefilled
  }
  #BLOB$chol_Q <- as(Matrix::.bdiag(chol_Q_list),"dtCMatrix") # !not dtC by default! (!). Has worked as (triangular) dgCMatrix for a long time.
  BLOB$chol_Q <- .bdiag_dtC(chol_Q_list) # assuming all chol_Q are lower triangular
  return(precisionBlocks)
}

..calc_ZtWX <- function(AUGI0_ZX, WLS_mat_weights, X.pv=AUGI0_ZX$X.pv) {
  if (attr(WLS_mat_weights,"is_unit")) {
    ZtWX <- as.matrix(crossprod(AUGI0_ZX$ZAfix, X.pv)) 
  } else {
    if (methods::.hasSlot(AUGI0_ZX$ZAfix, "x") && # 1st condition is true in all tests. But not clearly enforced. (__F I X M E__ ?) 
        length(AUGI0_ZX$ZAfix@x)>prod(dim(X.pv))) {
      ZtWX <- .Dvec_times_m_Matrix(WLS_mat_weights, X.pv) 
      ZtWX <- as.matrix(crossprod(AUGI0_ZX$ZAfix, ZtWX)) 
    } else {
      ZtWX <- .Dvec_times_m_Matrix(WLS_mat_weights, AUGI0_ZX$ZAfix) 
      ZtWX <- as.matrix(crossprod(ZtWX, X.pv)) 
    }
  }
  ZtWX
}

.calc_ZtWX <- function(sXaug) {
  AUGI0_ZX <- sXaug$AUGI0_ZX
  bloclength <- ceiling(1e7/nrow(AUGI0_ZX$X.pv))
  nc <- ncol(AUGI0_ZX$X.pv)
  if (nc==0L) { ## pforpv=0
    ZtWX <- matrix(0, nrow=ncol(AUGI0_ZX$ZAfix), ncol=0L) # with the 0 the matrix in num (as in general case) rather than logi (logi effect not tested)
  } else if (nc>bloclength) {
    seqncol <- seq(nc)
    blocs <- split(seqncol, ceiling(seq_along(seqncol)/bloclength))
    ZtWX <- vector("list",length(blocs))
    for (it in seq_along(blocs)) {
      ZtWX[[it]] <- ..calc_ZtWX(AUGI0_ZX, WLS_mat_weights=sXaug$BLOB$WLS_mat_weights, X.pv=AUGI0_ZX$X.pv[,blocs[[it]],drop=FALSE])
    }
    ZtWX <- do.call(cbind, ZtWX)
  } else ZtWX <- ..calc_ZtWX(AUGI0_ZX, WLS_mat_weights=sXaug$BLOB$WLS_mat_weights)
  ZtWX ## always *m*atrix
}

# in-fit function
# cf comments in .calc_r22 that also computes crossprod_r22 but returns its chol.
.calc_inv_beta_cov <- function(sXaug) { ## crossprod_r22 without computing r22 (nor r12) when r22 is ot already available
  pforpv <- attr(sXaug, "pforpv") 
  if (pforpv==0L) return(Diagonal(n=0L))
  # ELSE
  ## In the working doc it is shown that the beta,beta block of inv_d2hdbv2 is solve(crossprod(r22)). We compute crossprod(r22) here,
  ## We know this crossprod is X'inv(Sig_e)X - r12'r12 in the working doc (not X'inv(Sig_e)X ! ) since r22 is defined from this. 
  AUGI0_ZX <- sXaug$AUGI0_ZX
  BLOB <- sXaug$BLOB
  bloclength <- ceiling(1e7/nrow(AUGI0_ZX$X.pv)) 
  if (pforpv>bloclength) { ## I felt a need to avoid crossprod_r12 when r12 is tall (even if r12 is slim so that its crossprod is small)
    seqncol <- seq(pforpv)
    row_blocs <- split(seqncol, ceiling(seq_along(seqncol)/bloclength))
    crossprod_r22 <- vector("list",length(row_blocs))
    if ( ! is.null(BLOB$signs)) {
      WX <- as.matrix(.Dvec_times_m_Matrix(BLOB$WLS_mat_weights,AUGI0_ZX$X.pv)) # not sqrt
      for (it in seq_along(row_blocs)) {
        crossprod_r22[[it]] <- crossprod(BLOB$r12, BLOB$r12[,row_blocs[[it]],drop=FALSE]) # don't use .crossprod on BLOB$r12 for numerical precision.
        if(! is.null(BLOB$G_scaling)) crossprod_r22[[it]] <- crossprod_r22[[it]]/BLOB$G_scaling
        crossprod_r22[[it]] <- .crossprod(WX, AUGI0_ZX$X.pv[,row_blocs[[it]],drop=FALSE], as_mat=TRUE) - crossprod_r22[[it]]
      }
    } else {
      sqrtWX <- .Dvec_times_m_Matrix( BLOB$sqrtW,AUGI0_ZX$X.pv) 
      for (it in seq_along(row_blocs)) {
        crossprod_r22[[it]] <- crossprod(BLOB$r12, BLOB$r12[,row_blocs[[it]],drop=FALSE]) # don't use .crossprod on BLOB$r12 for numerical precision.
        if(! is.null(BLOB$G_scaling)) crossprod_r22[[it]] <- crossprod_r22[[it]]/BLOB$G_scaling
        crossprod_r22[[it]] <- .crossprod(sqrtWX, sqrtWX[,row_blocs[[it]],drop=FALSE], as_mat=TRUE) - crossprod_r22[[it]]
      }
    }
    crossprod_r22 <- do.call(cbind, crossprod_r22)
  } else {
    crossprod_r12 <- crossprod(as.matrix(BLOB$r12))  # Xt_X ! # using most precise computation as in .calc_r22()
    crossprod_r22 <- BLOB$XtWX - crossprod_r12 ## XtWX-Xt_X ! Both lines as explained in working doc
  }
  return(crossprod_r22) # dgeMatrix
}



.calc_r22_postfit <- function(X.pv, H_w.resid, r12, XtWX=NULL) { # currently post-fit bc called conditionally on is.null(r22 <- envir$r22)
  ## lots of alternatives here, removed from [v2.6.54. See comments in delayedAssign("r22"...)
  if (ncol(X.pv)) { 
    if (is.null(XtWX))  # currently TRUE because in post-fit calls, XtWX remains NULL (FIXME? Bof)
      XtWX <- .WLS_mat_XtWX(X.pv, WLS_mat_weights=H_w.resid, signs_in_WLS_mat=any(H_w.resid<0))
    crossr22 <- XtWX - as.matrix(crossprod(r12)) 
    return(.wrap_Utri_chol(crossr22)) # rather than (Matrix::)chol()
  } else return(diag(nrow=0L))
}

# returns a *m*atrix bc used in differences of matrices, faster when both matrices are *m*atrices.
# The following as.matrix()'s are a safety for rare case where X.pv is sparse
.WLS_mat_XtWX <- function(X.pv, WLS_mat_weights, signs_in_WLS_mat) {
  if (attr(WLS_mat_weights,"is_unit")) { 
    .crossprod(as.matrix(X.pv), NULL, as_mat = TRUE)
  } else  if (signs_in_WLS_mat) { 
    wX <- .Dvec_times_m_Matrix(WLS_mat_weights, X.pv)
    .crossprod(wX,X.pv, as_mat = TRUE)
  } else as.matrix(.ZtWZwrapper(X.pv,WLS_mat_weights)) # weights>0
}

.solve_crossr22 <- function(BLOB, # must include $r12
                            AUGI0_ZX, sXaug, dbeta_rhs, 
                            use_crossr22=FALSE) { #use_crossr22=TRUE affects mv test 'zut2_testr22' (test-mv-nested, v3.6.46, l.768)
  if (use_crossr22 && ! .is_evaluated("r22",BLOB) ) { # FALSE, hence currently BLOB$r22 is always used 
    dbeta_eta <- solve(BLOB$crossr22, dbeta_rhs)
  } else dbeta_eta <- backsolve(BLOB$r22, backsolve(BLOB$r22, dbeta_rhs, transpose = TRUE)) # ie (chol2inv(BLOB$r22) %*% dbeta_rhs)[,1L]
  dbeta_eta
}

.calc_sum_pwt_Q_y_o_2  <- function(sXaug, pwy_o) { 
  ## Result z, viewed as a col vector,  is such that sum(z^2)=c(0,z)'c(0,z)=c(0,pwy_o)'QQ'c(0,pwy_o)  # ('O' of length n_u_h)
  ##   for the orthog factor Q of the augmented matrix.
  ## The names of the variables come from the analogy with the hatval_ZX code.
  get_from_MME(sXaug, which="initialize") ## sets G_CHMfactor if not already set (possibly useless in current code)
  BLOB <- sXaug$BLOB
  AUGI0_ZX <- sXaug$AUGI0_ZX
  Zt_sqrtw_pwy_o <- crossprod(AUGI0_ZX$ZAfix, sqrt(sXaug$BLOB$H_w.resid)*pwy_o) 
  Xt_sqrtw_pwy_o <- crossprod(AUGI0_ZX$X.pv, sqrt(sXaug$BLOB$H_w.resid)*pwy_o)
  if (.is_evaluated("invL_G.P", BLOB)) { 
    lev_phi_z_pwy_o <- BLOB$invL_G.P %*% Zt_sqrtw_pwy_o
  } else lev_phi_z_pwy_o <- Matrix::solve(BLOB$G_CHMfactor, Zt_sqrtw_pwy_o[BLOB$pMat_G@perm, ], system="L") ## R_11^{-T}.Ztw.pwy_o
  if (attr(sXaug,"pforpv")>0L) {
    if (FALSE) { # correct but maybe slower
      crossprod_r12_z_pwy_o <- crossprod(BLOB$r12,lev_phi_z_pwy_o) ## R_12^T . R_11^{-T}.Ztw.pwy_o
      u_of_quadratic_utAu <- Xt_sqrtw_pwy_o - as.matrix(crossprod_r12_z_pwy_o)
      sum_lev_phi_x_pwy_o_2 <- .crossprod(u_of_quadratic_utAu, 
                                          .solve_crossr22(BLOB, # requires $r12
                                                          AUGI0_ZX, sXaug, u_of_quadratic_utAu, use_crossr22=TRUE))
      return(sum(lev_phi_z_pwy_o^2)+sum_lev_phi_x_pwy_o_2)
    } else { 
      crossprod_r12_z_pwy_o <- crossprod(BLOB$r12,lev_phi_z_pwy_o) ## R_12^T . R_11^{-T}.Ztw.pwy_o
      lev_phi_x_pwy_o <- backsolve(BLOB$r22, Xt_sqrtw_pwy_o - as.matrix(crossprod_r12_z_pwy_o), # as.matrix() has a notable effect on ohio|fitme test, and no visible drawback
                                   transpose=TRUE) ## -  R_22^{-T}.R_12^T . R_11^{-T}.Ztw.pwy_o + R_22^{-T}.Xtw.pwy_o
      return(sum(lev_phi_z_pwy_o^2)+sum(lev_phi_x_pwy_o^2))
    }
  } else return(sum(lev_phi_z_pwy_o^2))
} # "pwt_Q_y_o"

# called only for spprec "hatval", given conditions calling for the use of qrXa 
.calc_spprec_hatval_ZX_by_QR <- function(BLOB, sXaug, AUGI0_ZX, w.ranef) { ## sparse qr is fast
  n_u_h <- ncol(AUGI0_ZX$I)
  phipos <- n_u_h+seq_len(length(tmp)-n_u_h)
  if (BLOB$nonSPD) {
    tmp <- colSums(( BLOB$invIm2QtdQ_ZX %*% BLOB$t_Q) * BLOB$t_Q)
    tmp[phipos] <- tmp[phipos]*BLOB$signs
  } else {
    t_Q <- BLOB$t_Q
    xx <- t_Q@x
    xx <- xx*xx
    t_Q@x <- xx
    tmp <- colSums(t_Q)
  }
  list(lev_lambda=tmp[-phipos],lev_phi=tmp[phipos]) # BLOB$hatval , formely hatval_ZX, for REML
}

.calc_spprec_hatval_ZX_by_r22 <- function(BLOB, sXaug, AUGI0_ZX, w.ranef) {
  if (BLOB$nonSPD) stop("this should not be called in nonSPD case") # nonSPD => QR approach with invIm2Q... correction being the only feasible one.  
  # We use the  Chol factorization of T'T as explained in the working doc
  # We know one of its block is the same as for the ML leverages, we compute two additional blocks: 
  lev_lambda_z <- BLOB$factor_inv_Md2hdv2 ## ul block of R^{-T} as described in working doc
  xx <- lev_lambda_z@x
  xx <- xx*xx
  lev_lambda_z@x <- xx
  lev_lambda_z <- .Matrix_times_Dvec(lev_lambda_z,w.ranef)
  lev_lambda_z <- colSums(lev_lambda_z)
  #
  if (attr(sXaug,"pforpv")>0L) {
    lev_lambda_x <- - backsolve(BLOB$r22,crossprod(BLOB$r12, BLOB$factor_inv_Md2hdv2),transpose=TRUE) ## ll block of R^{-T} as described in working doc
    lev_lambda_x <- lev_lambda_x*lev_lambda_x
    ## one dimension is dropped when (?) X is a single column matrix
    if ( is.matrix(lev_lambda_x)) {  ## m ou M atrix
      lev_lambda_x <- colSums(.m_Matrix_times_Dvec(lev_lambda_x,w.ranef))
    } else {
      lev_lambda_x <- lev_lambda_x * w.ranef
    }
    lev_lambda <- lev_lambda_z+lev_lambda_x
  } 
  if (AUGI0_ZX$is_unitary_ZAfix) {  
    # no complete simplif for AUGI0_ZX$is_unitary_ZAfix bc lev_phi_x is not the diag of a Z...Zt matrix
    # but we can still use p*n matrices instead of r*n matrices to compute the lev_phi_x component
    ## 
    sqrtwZ <- .Dvec_times_Matrix(BLOB$sqrtW, AUGI0_ZX$ZAfix)
    crossprod_r12_z_rows <- crossprod(crossprod(BLOB$invL_G.P, BLOB$r12),t(sqrtwZ)) 
    lev_phi_z <- colSums(BLOB$invL_G.P^2)
    lev_phi_z <- drop(AUGI0_ZX$ZAfix %*% lev_phi_z)*BLOB$WLS_mat_weights  # as in hatval_Z code.
  } else {
    lev_phi_z <- BLOB$inv_L_G_ZtsqrW
    crossprod_r12_z_rows <- crossprod(BLOB$r12,lev_phi_z)
    xx <- lev_phi_z@x
    xx <- xx*xx
    lev_phi_z@x <- xx
    lev_phi_z <- colSums(lev_phi_z) 
    if  ( ! is.null(BLOB$signs)) lev_phi_z <- lev_phi_z*BLOB$signs
  }
  if (attr(sXaug,"pforpv")>0L) {
    sqrtwX <- .Dvec_times_m_Matrix(BLOB$sqrtW,AUGI0_ZX$X.pv)
    lev_phi_x <- backsolve(BLOB$r22, t(sqrtwX) - crossprod_r12_z_rows, transpose=TRUE)
    lev_phi_x <- lev_phi_x^2
    if ( is.matrix(lev_phi_x)) lev_phi_x <- colSums(lev_phi_x) 
    if  ( ! is.null(BLOB$signs)) lev_phi_x <- lev_phi_x*BLOB$signs
    lev_phi <- lev_phi_z+lev_phi_x
    list(lev_lambda=lev_lambda,lev_phi=lev_phi) # BLOB$hatval , formely hatval_ZX, for REML
  } else list(lev_lambda=lev_lambda_z,lev_phi=lev_phi_z) # BLOB$hatval , formely hatval_ZX, for REML
}

.calc_H_dH <- function(BLOB, damping) {
  ## See comments in .calc_G_dG()
  H_dH <- BLOB$Md2hdv2
  nc <- ncol(H_dH)
  diagPos <- seq.int(1L,nc^2,nc+1L)
  dH <- H_dH[diagPos]*damping
  H_dH[diagPos] <- H_dH[diagPos] + dH
  return(list(H_dH=H_dH, dampDpD_2=dH))
}

.calc_G_dG <- function(BLOB, damping) {
  spprec_LevM_D <- .spaMM.data$options$spprec_LevM_D
  ################
  as_sym <- FALSE # If as_sym is TRUE, the final operation BLOB$Gmat + dampdG is dsC + dsC = forceSymmetric(callGeneric(as(e1, "dgCMatrix"), as(e2, "dgCMatrix"))) (with generic for '+')
  # so having previously forced symmetry on dG is a waste of time=> as_sym=FALSE # ___F I X M E___ retry now that we have .dsCsum
  if (spprec_LevM_D=="update") { # experimental. Effect dependent a priori on .spaMM.data$options$perm_G. OK but not faster with default permG (TRUE) 01/2020 
    BLOB$D_Md2hdv2 <- diag(chol2inv(BLOB$chol_Q))
    BLOB$dG <- NaN ## To detect problems in further usages
    G_dG <- Matrix::.updateCHMfactor(BLOB$G_CHMfactor, BLOB$Gmat, mult=damping) # actually method of stats:::update for class CHMfactor (?`CHMfactor-class`) 
    return(list(G_dG=G_dG, dampDpD_2=damping * BLOB$D_Md2hdv2))
  } else {
    if (is.null(BLOB$dG)) {  
      if (spprec_LevM_D=="1") { # default
        BLOB$D_Md2hdv2 <- rep(1,ncol(BLOB$chol_Q))
        BLOB$dG <- drop0(.tcrossprod(BLOB$chol_Q, as_sym=as_sym)) # dsCMatrix if as_sym is TRUE# drop0 makes a diff in LevM.Frailty test # do we have the precisionMatrix somewhere ?
      } else { ## experimental: costly solve() for tcrossfac_Md2hdv2
        ## part of the problem is avoiding tall matrices from tall ZA (and X), but for squarish ZA, using solve(chol_Q,... ) looks complicated.
        ## solve(chol_Q,... ) further assuming that we have not stored the LMatrix
        tmp <- BLOB$tcrossfac_Md2hdv2 ## (Hfac)
        xx <- tmp@x
        xx <- xx*xx
        tmp@x <- xx
        # def of pertubration D_Md2hdv2 affects decimals in test-adjacency-long
        if (spprec_LevM_D=="rowSums") { ## [as originally used for (full) LevMar_step]
          BLOB$D_Md2hdv2 <- rowSums(tmp) # the diagonal elements since *row*Sums = diag Tcrossprod(tcrossfac) (= invQ G Gt invQt)
        } else { 
          BLOB$D_Md2hdv2 <- colSums(tmp) # colSums() seems to give good results
        }
        ## convert diag perturb of Md2hdv2 into a non-diag perturb of G :
        ## since H=invL_Q G t(invL_Q),  dG= L_Q dH t(L_Q) 
        BLOB$dG <- drop0(.ZWZtwrapper(BLOB$chol_Q , BLOB$D_Md2hdv2, as_sym=as_sym)) # dsCMatrix if as_sym is TRUE
      } 
    }
    dampdG <- (damping*BLOB$dG) ## not always diagonal...
    dsC_Gmat <- inherits(BLOB$Gmat,"dsCMatrix")
    dsC_ddG <- inherits(dampdG,"dsCMatrix")
    if (dsC_Gmat && dsC_ddG) {
      G_dG <- .dsCsum(BLOB$Gmat, dampdG) # BLOB$Gmat + dampdG ## probably not so sparse... yet this occurs in the tests
      #G_dG <- .dsC_plus_dsC(BLOB$Gmat,dampdG) ##
    } else {
      if (dsC_Gmat || dsC_ddG) warning("possibly inefficient code in .calc_G_dG()")
      G_dG <- forceSymmetric(BLOB$Gmat + dampdG)
    }
    return(list(G_dG=G_dG, dampDpD_2=damping * BLOB$D_Md2hdv2))
  }
}

.provide_BLOB_hatval_Z_ <- function(sXaug, BLOB=sXaug$BLOB, w.ranef=attr(sXaug,"w.ranef") , AUGI0_ZX=sXaug$AUGI0_ZX , needed=c("lambda","phi")) {
  if (BLOB$nonSPD) { # __F I X M E__ componentwise product presumably not efficient
    hatval_Z_ <- colSums(( BLOB$invIm2QtdQ_Z %*% BLOB$t_Qq) * BLOB$t_Qq)
    seq_n_u_h <- seq_len(ncol(BLOB$Gmat))
    BLOB$hatval_Z_ <- list(lev_lam=hatval_Z_[seq_n_u_h], lev_phi=hatval_Z_[-seq_n_u_h]*BLOB$signs)
  } else { 
    # Chol-based approach H_w.resid=WLS_mat_weights
    if ("lambda" %in% needed) {
      if (is.null(BLOB$hatval_Z_$lev_lambda)) {
        lev_lambda <- BLOB$factor_inv_Md2hdv2
        xx <- lev_lambda@x
        xx <- xx*xx
        lev_lambda@x <- xx
        BLOB$hatval_Z_$lev_lambda <- colSums(lev_lambda) * w.ranef
      }
    }
    if ("phi" %in% needed) {
      if (is.null(BLOB$hatval_Z_$lev_phi)) {
        if (AUGI0_ZX$is_unitary_ZAfix) { ## then only diagonal values of invG matter  ## adjacency-long case...?
          ## invL_G.P presumably needed for lev_lambda from factor_inv_Md2hdv2
          ## (...lev_phi from invL_G.P, lev_lambda from invL_G.P %*% chol_Q ) 
          lev_phi <- colSums(BLOB$invL_G.P^2)
          # __F I X M E__ as.vector() otherwise next line puts the weights' attributes on the lev_phi result. this is ugly 
          # here SPD case, H_w.resid=WLS_mat_weights are signed so they contain BLOB$signs 
          lev_phi <- as.vector(drop(AUGI0_ZX$ZAfix %*% lev_phi) *sXaug$BLOB$WLS_mat_weights) # signs here are checked by compar with SPD spcorr sign() QQ' 
        } else {
          # inv_L_G_ZtsqrW arises as simplif of BLOB$factor_inv_Md2hdv2 %*% t(sqrtwZL) (L and chol_Q cancel each other)
          lev_phi <- BLOB$inv_L_G_ZtsqrW ## solve(as(BLOB$G_CHMfactor,"sparseMatrix"), t(sqrtwZ))
          xx <- lev_phi@x
          xx <- xx*xx
          lev_phi@x <- xx
          lev_phi <- colSums(lev_phi)
          if  (! is.null(BLOB$signs)) lev_phi <- lev_phi*BLOB$signs # signs here consistent with unitary_ZAfix case, already independently checked
        }
        BLOB$hatval_Z_$lev_phi <- lev_phi
      }
    }
  }
  BLOB$hatval_Z_ 
}

.provide_BLOB_factor_inv_Md2hdv2 <- function(BLOB) { ## Code encapsulated in a function for easier profiling
  # .calc_sscaled_new() -> (Pdiag <- get_from_MME(sXaug,"hatval_Z", B=...)) -> appears to be the last big bottleneck for fitme(). 
  # We need a diag(crossprod(factor)) =colSums(factor$x^2) needing the full factor
  LL <- BLOB$invL_G.P %*% BLOB$chol_Q ## dgCMatrix 
  # Even is BLOB$chol_Q is Identity, this this hardly slower than Matrix::solve(BLOB$G_CHMfactor, system="L") so cannot be improved.
  # The pattern Matrix::solve(BLOB$G_CHMfactor, as(BLOB$G_CHMfactor,"pMatrix") %*% ., system="L")  occurs repeatedly and is costly
  #                        late discovery: as(BLOB$G_CHMfactor,"pMatrix") is quite costly
  # so the Q is whether to store the (rel dense) invL_G.P= solve(BLOB$G_CHMfactor, as(BLOB$G_CHMfactor,"pMatrix"))
  # Currently it is computed and stored only for several leverage computation and for factor_inv_Md2hdv2
  #factor_Md2hdv2 <- as.matrix(Matrix::solve(BLOB$chol_Q,as(BLOB$G_CHMfactor,"sparseMatrix")))
  #BLOB$chol_Md2hdv2 <- Rcpp_chol_R(tcrossprod(factor_Md2hdv2)) #A = R'.R as in R's chol()
  factor_inv_Md2hdv2 <- drop0(LL,tol=.Machine$double.eps) ## crossfactor
} # *returns* the factor

.factor_inv_Md2hdv2_times_rhs <- function(BLOB, B) { 
  ## get_from_MME(., "solve_d2hdv2", .) => B generally vector but may be diag, Diagonal 
  ## cf HLfit_body -> .hatvals2std_lev -> -> .calc_dvdloglamMat_new()
  ## or even fuller matrices (neg.d2h0_dv_dlogphi <- .m_Matrix_times_Dvec(t(ZAL), drop(dh0deta)) in .calc_dvdlogphiMat_new())
  if (.is_evaluated("factor_inv_Md2hdv2", BLOB)) {
    BLOB$factor_inv_Md2hdv2 %*% B # __F I X M E__ optimize for diagonal matrices? Not easy (info to be passed efficiently from get_from_MME call)
  } else {
    B <- BLOB$chol_Q %*% B # => which typically yields a dge, so [] is OK:
    B <- B[BLOB$pMat_G@perm, ] # BLOB$pMat_G %*% B
    Matrix::solve(BLOB$G_CHMfactor, B, system="L")  
  }
} # drop0() might be useful (__F I X M E__? have to wait for profiling to point it...)

.crossprod_factor_inv_Md2hdv2_rhs <- function(BLOB, B) { # B is typically grad_v
  if (.is_evaluated("factor_inv_Md2hdv2", BLOB)) {
    .crossprod(BLOB$factor_inv_Md2hdv2, B)
  } else {
    B <- Matrix::solve(BLOB$G_CHMfactor, B, system="Lt")
    B <- B[BLOB$sortPerm, ] # t(pMatrix) %*% vec
    .crossprod(BLOB$chol_Q, B) 
  }
}

.calc_r12 <- function(BLOB, X.pv) {
  
  # Element of the WLS_mat, even in nonSPD case. if there are signs in WLS_mat weights we must recover them.  
  # Use inv_L_G_ZtsqrW when it's available (as needed for certain leverage computations) 
  # Otherwise use ZtWX which is typically much smaller, so that solve(., ZtWX) is faster than solve(., ZtsqrW)
  if (.is_evaluated("inv_L_G_ZtsqrW", BLOB)) { 
    if (BLOB$signs_in_WLS_mat) {
      sqrtwX <- .Dvec_times_m_Matrix(BLOB$signs * BLOB$sqrtW, X.pv)
    } else sqrtwX <- .Dvec_times_m_Matrix(BLOB$sqrtW,X.pv)
    r12 <- BLOB$inv_L_G_ZtsqrW %*% sqrtwX
  } else if (.is_evaluated("invL_G.P", BLOB)) { 
    r12 <- BLOB$invL_G.P %*% BLOB$ZtWX
  } else r12 <-  Matrix::solve(BLOB$G_CHMfactor, BLOB$ZtWX[BLOB$pMat_G@perm, ], system="L")
  # : not computing invL_G.P here is big improvement for test IMRF rawGNIPdataEU.1, and also for pcmatern LMM
  as.matrix(r12) # it tends to be dense (dge) and is used in a few crossprod operations that are faster on matrix (Matrix::crossprod time ~time crossprod(as.matrix)).
  # Using .Rcpp_crossprod() on the result may be even faster, but overall a notable waste of time by loss of numerical precision
  # => DON'T use .crossprod() not .Rcpp_crossprod on it.
  
}

.ad_hoc_dsy_warning <- local({
  dsy_warned <- FALSE
  function() {
    if ( ! dsy_warned) {
      dsy_warned <<- TRUE
      if (is.null(.spaMM.data$options$sparse_precision)) {
        warning("Sparse-precision algorithm (automatically selected) possibly inefficient. Please report to the package maintainer.", immediate.=TRUE)
      } else {
        warning("Sparse-precision algorithm (selected by user) possibly inefficient.", immediate.=TRUE)
      }
    }
  }
})

.init_promises_spprec <- function(sXaug, non_X_ones=TRUE, nullify_X_ones =FALSE, intervalInfo=NULL) {  
  BLOB <- sXaug$BLOB # environment that should contain  'G_CHMfactor', 'chol_Q', 'perm', optionally 'signs' when promises are evaluated;
  AUGI0_ZX <- sXaug$AUGI0_ZX # envir (prime fit) or list (converted to list by confint)
  # 
  if ( non_X_ones ) {
    ### What does NOT depend on X
    delayedAssign("sortPerm", sort.list(BLOB$pMat_G@perm), assign.env = BLOB ) # never NULL
    delayedAssign("invL_G.P", drop0(Matrix::solve(BLOB$G_CHMfactor, BLOB$pMat_G, system="L")), assign.env = BLOB ) # never NULL
    delayedAssign("inv_L_G_ZtsqrW", {
      sqrtwZ <- .Dvec_times_Matrix(BLOB$sqrtW, AUGI0_ZX$ZAfix)
      inv_L_G_ZtsqrW <- tcrossprod(BLOB$invL_G.P, sqrtwZ)
    }, assign.env = BLOB )
    delayedAssign("factor_inv_Md2hdv2",.provide_BLOB_factor_inv_Md2hdv2(BLOB), assign.env = BLOB )
    #delayedAssign("half_logdetQ",  sum(log(diag(x=BLOB$chol_Q))), assign.env = BLOB ) ## currently not used (in variant of LevMar step)
    #delayedAssign("half_logdetG", Matrix::determinant(BLOB$G_CHMfactor)$modulus[1], assign.env = BLOB ) ## currently not used (in variant of LevMar step)
    delayedAssign("logdet_R_scaled_v", { BLOB$logdet_sqrt_d2hdv2 - sum(log(attr(sXaug,"w.ranef")))/2 }, assign.env = BLOB )
    # more or less for speculative code:
    delayedAssign("Md2hdv2", .tcrossprod(BLOB$tcrossfac_Md2hdv2), assign.env = BLOB ) ## currently not used (for H_dH)
    delayedAssign("tcrossfac_Md2hdv2",
                  Matrix::solve(BLOB$chol_Q, crossprod(BLOB$pMat_G, as(BLOB$G_CHMfactor,"sparseMatrix"))), assign.env = BLOB )
    if (BLOB$nonSPD) {
      ## I need the leverages for the gradient, not from the then-unsigned WLS_mat
      # For hatval_Z, construct the t_Qq matrix and invIm2QtdQ_Z, as in sp|decorr nonSPD case
      # For hatval_ZX, a similar approach using qrXa...  below
      delayedAssign("t_Qq", { cbind(tcrossprod(BLOB$invL_G.P, .sparseDiagonal(x= sqrt(attr(sXaug,"w.ranef")), shape="g")), 
                                    BLOB$inv_L_G_ZtsqrW) }, assign.env = BLOB) 
      delayedAssign("invIm2QtdQ_Z", { .Rcpp_adhoc_shermanM_sp(BLOB$t_Qq, c(rep(0L,ncol(BLOB$Gmat)),BLOB$signs<0)) }, assign.env = BLOB)
      delayedAssign("logdet_sqrt_d2hdv2", { # keep in mind that L, u and hlik will differ from correl algos; Lu, mu, clik and p_v will match  
        half_logdetG <- Matrix::determinant(BLOB$G_CHMfactor)$modulus[1]
        half_logdetQ <- sum(log(diag(x=BLOB$chol_Q)))
        half_logdetG - half_logdetQ - determinant(BLOB$invIm2QtdQ_Z)$modulus[1]/2 # logdet of WLS_mat + correction to get logdet of true Hessian 
      }, assign.env = BLOB )
    } else {
      delayedAssign("logdet_sqrt_d2hdv2", { # keep in mind that L, u and hlik will differ from correl algos; Lu, mu, clik and p_v will match  
        half_logdetG <- Matrix::determinant(BLOB$G_CHMfactor)$modulus[1]
        half_logdetQ <- sum(log(diag(x=BLOB$chol_Q)))
        half_logdetG - half_logdetQ 
      }, assign.env = BLOB )
    }
  }
  #
  ### What DOES depend on X:
  #  Most of them assuming unscaled X (for confint code: lost columns; and predVar computations currently all based on unscaled beta's), BUT
  # one exception is qrXa. Ggiven it's needed sometimes in-fit, the post-fit code using it takes the rescaling into account, so we treat it as the non-X ones, 
  if ( non_X_ones ) {
    delayedAssign("qrXa", { ## used for "hatval" in case of failure of computation of BLOB$r22... and for BLOB$nonSPD
      if (.spaMM.data$options$Matrix_old) { # this block appears to evade the long tests
        X.pv <- as(AUGI0_ZX$X.pv,"dgCMatrix")#  .Rcpp_as_dgCMatrix(AUGI0_ZX$X.pv)
      } else X.pv <- as(X.pv <- as(AUGI0_ZX$X.pv,"generalMatrix"),"CsparseMatrix")
      Xscal <- with(sXaug,rbind(cbind(AUGI0_ZX$I, AUGI0_ZX$ZeroBlock), # I0_ZX order
                                cbind(AUGI0_ZX$ZAfix %*% t(.rawsolve(BLOB$chol_Q)), X.pv)))
      Xscal <- .Dvec_times_Matrix(c(sqrt(attr(sXaug,"w.ranef")),BLOB$sqrtW),Xscal)
      qrXa <- qr(Xscal)
    }, assign.env = BLOB )
  }
  if (nullify_X_ones) {
    if (.is_evaluated("r22", BLOB)) {
      # Then we use it for beta_cov but must be careful not to mix it with the other objects bc it is the r22 for a scaled X while the new promises for post fit
      # refer to an unscaled version of X.pv. ./.
      unsc_r22 <-  .m_Matrix_times_Dvec(BLOB$r22, attr(AUGI0_ZX$X.pv,"scaled:scale")) 
      beta_cov <- solve(crossprod(unsc_r22))
      colnames(beta_cov) <- rownames(beta_cov) <- setdiff(colnames(AUGI0_ZX$X.pv), intervalInfo$parm)
      BLOB$beta_cov <- beta_cov # distinct variable from $beta_cov_info$beta_cov
    } # ./. When BLOB$r22 is present and then $r12 is also expected by .old_calc_Md2hdvb2_info_spprec_by_r22() 
    # => if we plan to use BLOB$r22 post fit we must also provide BLOB$r12. Currently we remove them and all X-related promises: 
    BLOB$ZtWX <-BLOB$XtWX <- BLOB$DpD <- BLOB$r12 <- BLOB$r22 <- BLOB$LZtWX <- BLOB$crossr22 <- NULL 
    if (BLOB$nonSPD) t_Q <- invIm2QtdQ_ZX <- logdet_R_scaled_b_v <- NULL
  } else {
    
    delayedAssign("ZtWX", .calc_ZtWX(sXaug), assign.env = BLOB ) # *m*atrix
    delayedAssign("XtWX", .WLS_mat_XtWX(X.pv=AUGI0_ZX$X.pv, WLS_mat_weights=BLOB$WLS_mat_weights, signs_in_WLS_mat=BLOB$signs_in_WLS_mat)
                  , assign.env = BLOB )  ## as(,"matrix") effect included within .WLS_mat_XtWX def.
    delayedAssign("r12", .calc_r12(BLOB, X.pv=AUGI0_ZX$X.pv), assign.env = BLOB )
    delayedAssign("r22",{ # => the former .calc_r22() code
      ## i's vain to try to regularize crossr22 by manipulating diagonal or eigenvalues. 
      #  It would be better to control the accuracy of the computation of crossr22 as a difference.
      if (ncol(AUGI0_ZX$X.pv)) { 
        # Element of the WLS_mat, even in nonSPD case. if there are signs in WLS_mat weights we must recover them => BLOB$XtWX and  BLOB$r12 must include them 
        crossr22 <- BLOB$XtWX - crossprod(BLOB$r12) 
        # .Rcpp_crossprod(r12,NULL,as_mat=TRUE) may be faster but less accurate numerically. Numerical precision important here.
        # + see https://stackoverflow.com/questions/52888650/eigenselfadjointviewrankupdate-slower-than-a-ww-transpose
        # A test is in test-mv-nested: logLik(mrf1T <- fitme(migStatus ~ 1 +multIMRF(1... should give phi=1e-6 and p_v=-11.24422
        # 
        # Different numerics for the same pb: 
        #  .calc_inv_beta_cov(sXaug) computes crossr22 from sXaug, handling some complication, 
        # .get_absdiagR_blocks() face similar pbs and finally calls .Utri_chol_by_qr()
        .wrap_Utri_chol(crossr22) # rather than (Matrix::)chol()
      } else diag(nrow=0L)
      ## lots of alternatives here, removed from [v2.6.54
    }, assign.env = BLOB )
    delayedAssign("DpD", c(BLOB$D_Md2hdv2,diag(x=BLOB$XtWX)), assign.env = BLOB )
    #
    # more or less for speculative code:
    delayedAssign("LZtWX", as(solve(BLOB$chol_Q, BLOB$ZtWX),"matrix"), assign.env = BLOB ) ## currently not used (in variant of LevMar step)
    delayedAssign("crossr22", { # NOT currently used: Used only in-fit -> .calc_sum_pwt_Q_y_o_2() ->  .solve_crossr22(., use_crossr22=TRUE), but use_crossr22 set to FALSE:
      if (ncol(X.pv)) { 
        crossr22 <- BLOB$XtWX - crossprod(as.matrix(BLOB$r12),NULL,as_mat=TRUE) 
      } else crossr22 <- diag(nrow=0L)
    } , assign.env = BLOB )
    
    if (BLOB$nonSPD) {
      delayedAssign("t_Q", { t(qr.Q(BLOB$qrXa)) }, assign.env = BLOB) # (In QRP factos, crossprod X is always eqauld to crossprod Q)
      delayedAssign("invIm2QtdQ_ZX", { .Rcpp_adhoc_shermanM_sp(BLOB$t_Q, c(rep(0L,ncol(BLOB$Gmat)),BLOB$signs<0)) }, assign.env = BLOB)
      delayedAssign("logdet_R_scaled_b_v", {
        R_scaled  <- qrR(BLOB$qrXa,backPermute = FALSE)
        sum(log(abs(diag(x=R_scaled))))
      }, assign.env = BLOB)
      delayedAssign("logdet_r22", { # X_scaled_H_unscaled_logdet_r22
        # some fun here. A QR facto is needed.... I wrote p_bv comput for all methods as a correction of p_v computation
        # by the logdet of an r22 bloc which is not really available in nonSPD case. But we have a similar subcase in the alternative to nonSPD.
        # However, this implies two QR facto ___F I X M E___ is it possible to optimze to a single one?
        logdet_sqrt_d2hdbv2 <- determinant(qrR(BLOB$qrXa,backPermute = FALSE))$modulus[1] - determinant(BLOB$invIm2QtdQ_ZX)$modulus[1]/2
        logdet_sqrt_d2hdbv2 - BLOB$logdet_sqrt_d2hdv2
      }, assign.env = BLOB) 
    } else {
      delayedAssign( "logdet_R_scaled_b_v", { BLOB$logdet_R_scaled_v + BLOB$logdet_r22} , assign.env = BLOB )
      delayedAssign("logdet_r22", { # X_scaled_H_unscaled_logdet_r22
        if (.is_evaluated("r22",BLOB)) {
          logdet_r22 <- determinant(BLOB$r22)$modulus[1]
        } else {
          if (attr(sXaug,"pforpv")==0L) {
            logdet_r22 <- 0 
          } else if (inherits(AUGI0_ZX$X.pv,"sparseMatrix")) { ## sparse qr is fast
            logdet_r22 <- determinant(qrR(BLOB$qrXa,backPermute = FALSE))$modulus[1] - BLOB$logdet_sqrt_d2hdv2
          } else {
            #cat("cross_r22:") # FIXME should use X.Re for non-standard REML? (don't modify this but code which may call for it)
            crossprod_r22 <- .calc_inv_beta_cov(sXaug) # inv_beta_cov is crossprod_r22 as shown in working doc.
            logdet_r22 <- determinant(crossprod_r22)$modulus[1]/2
          }
        } 
        logdet_r22
      } , assign.env = BLOB )
    }

    delayedAssign( "hatval", { # hatval_ZX
      if (identical(.spaMM.data$options$use_spprec_QR,TRUE) || is.null(dim(BLOB$r22)) || BLOB$nonSPD) { # ((currently FALSE) || failure of comput of r22 || nonSPD)
        .calc_spprec_hatval_ZX_by_QR(BLOB=BLOB, AUGI0_ZX=AUGI0_ZX, sXaug=sXaug, w.ranef=attr(sXaug,"w.ranef"))
      } else .calc_spprec_hatval_ZX_by_r22(BLOB=BLOB, AUGI0_ZX=AUGI0_ZX, sXaug=sXaug, w.ranef=attr(sXaug,"w.ranef")) ## fragile bc  BLOB$r22 <- .calc_r22() is fragile...
    } , assign.env = BLOB )
    
  } # end block of promises depending on X
}

.calc_Gmat <- function(ZtWZ, precisionMatrix) {
  if (inherits(ZtWZ,"dsyMatrix")) { 
    # if tmp is small (2*2 and diagonal: ranef with two levels...), then it may have been returned by .ZtWZwrapper -> .crossprod as dsy, not dsC.
    # then, as(tmp,"sparseMatrix") is dgC not dsC
    #  and  as(tmp,"symmetricMatrix") is dsy not dsC
    # so there is no way to make sure that the result is dsC without adding as(.,"dsCMatrix") in .ZtWZwrapper or .crossprod
    .ad_hoc_dsy_warning()
    ZtWZ <- as(ZtWZ,"dsCMatrix")
  }
  if (inherits(ZtWZ,"dsCMatrix") && inherits(precisionMatrix,"dsCMatrix") ) { # test introduced 02/2020. 
    .dsCsum(ZtWZ, precisionMatrix) # faster than tmp + precisionMatrix where '+' calls forceSymmetric(callGeneric(as(e1, "dgCMatrix"), as(e2, "dgCMatrix")))
  } else { 
    stop("Unexpected matrix types in .AUGI0_ZX_sparsePrecision(). Please report to the package maintainer.") 
    # actually there is some other step that is silently wrong if both matrices are not dsC...
    ZtWZ + precisionMatrix
  }
}

# pure debug code 
.check_stepwise_sol <- function(sXaug, zInfo) {
  AUGI0_ZX <- sXaug$AUGI0_ZX # spprec!
  Xaug <- rbind(
    cbind(AUGI0_ZX$I,AUGI0_ZX$ZeroBlock),
    cbind(AUGI0_ZX$ZAfix,AUGI0_ZX$X.pv)
  )
  ww <- c(attr(sXaug,"w.ranef"), sXaug$BLOB$WLS_mat_weights)
  wXaug <- .Dvec_times_m_Matrix(sqrt(ww),Xaug) 
  
  solve(crossprod(wXaug), zInfo$m_grad_obj)
}


.AUGI0_ZX_sparsePrecision <- function(sXaug,which="",z=NULL,B=NULL,
                                      damping,LM_z) {
  BLOB <- sXaug$BLOB
  if (which=="initialize") return(NULL)
  if (which=="hatval_Z") {
    ## calcul correct; 
    # TT <- rbind(diag(sqrt(w.ranef)),diag(sqrt(sXaug$BLOB$H_w.resid)) %*% AUGI0_ZX$ZAfix %*% t(solve(BLOB$chol_Q)) )
    # diag(TT %*% (solve(crossprod(TT))) %*% t(TT)) ## lev_lambda, lev_phi
    ### if (is.null(B)) B <- c("lambda","phi") 
    return(.provide_BLOB_hatval_Z_(sXaug, BLOB=BLOB, needed=B))
  }
  if (which=="Mg_solve_g") { # only written for .pot4improv() in LevM code
    return(sum(B*unlist(get_from_MME(sXaug,szAug=list(m_grad_obj=B))))) # ASSUMING B=zInfo$m_grad_obj 
  } 
  if (which=="Mg_invH_g") { ## - grad_v invD2hdv2 grad_v 
    LLgrad <- .factor_inv_Md2hdv2_times_rhs(BLOB, B) #BLOB$factor_inv_Md2hdv2 %*% B
    return( sum(LLgrad^2))
  }
  if (which=="Mg_invXtWX_g") { ## 
    Mg_invXtWX_g <- crossprod(B,solve(BLOB$XtWX,B))[1L] # [1L] drops possible Matrix class...
    return(Mg_invXtWX_g)
  }
  if (which=="solve_d2hdv2") { 
    if ( is.null(B) ) {
      # stop("B is NULL") # this would actually 'work' bc this case never occurs in current code (v2.6.3)
      if (BLOB$nonSPD) {
        sol <- - .crossprod(BLOB$factor_inv_Md2hdv2, BLOB$invIm2QtdQ_Z %*% BLOB$factor_inv_Md2hdv2) 
      } else sol <-  - .crossprod(BLOB$factor_inv_Md2hdv2)
    } else if (BLOB$nonSPD) {
      sol <- BLOB$factor_inv_Md2hdv2 %*% B
      sol <- BLOB$invIm2QtdQ_Z %*% sol
      sol <- - .crossprod(BLOB$factor_inv_Md2hdv2, sol) 
    } else if (.is_evaluated("factor_inv_Md2hdv2", BLOB)) {
      sol <- - .crossprod(BLOB$factor_inv_Md2hdv2, BLOB$factor_inv_Md2hdv2 %*% B)
    } else {
      B <- BLOB$chol_Q %*% B
      B <- Matrix::solve(BLOB$G_CHMfactor, B, system="A")  
      sol <- - .crossprod(BLOB$chol_Q, B) 
    }
    return(sol)
  }
  #if (which=="half_logdetQ") { return(BLOB$half_logdetQ) }
  #if (which=="half_logdetG") { return(BLOB$half_logdetG) }
  # cases where Sigsolve is needed ## Sig is the Cov(response), tcrossprod(ZAL)/lambda+diag(1/w.resid) 
  #X.pv <- AUGI0_ZX$X.pv
  ## Solving for model coefficients:
  if ( ! is.null(z)) { ## which=""
    grad_v <- z$m_grad_obj[seq(ncol(BLOB$chol_Q))]
    rhs <- .factor_inv_Md2hdv2_times_rhs(BLOB, grad_v)[,1L]
    if (attr(sXaug,"pforpv")) { ## if there are fixed effect coefficients to estimate
      # rhs for beta
      grad_p_v <- z$m_grad_obj[-seq(ncol(BLOB$chol_Q))]
      dbeta_rhs <- grad_p_v - crossprod(BLOB$r12, rhs)[,1L] # numeric vector # using crossprod() rather than .crossprod() -> .Rcpp_crossprod() appears important here ./.
      #                                                      for numerical precision => for overall computation time.     
      dbeta_eta <-  .solve_crossr22(BLOB, AUGI0_ZX=sXaug$AUGI0_ZX, sXaug, dbeta_rhs)
      rhs <- rhs - (BLOB$r12 %*% dbeta_eta)[,1]
    } else {
      dbeta_eta <- numeric(0)
    }
    dv_h <- .crossprod_factor_inv_Md2hdv2_rhs(BLOB, rhs)[,1L]
    return(list(dv_h=dv_h,dbeta_eta=dbeta_eta)) # can be checked by .check_stepwise_sol(sXaug, zInfo)
  }
  if (which=="logdet_sqrt_d2hdv2") { return(BLOB$logdet_sqrt_d2hdv2)} 
  #### REML:
  if (which=="logdet_r22") { return(BLOB$logdet_r22) }
  if (which=="logdet_R_scaled_b_v") { return(BLOB$logdet_R_scaled_b_v)} 
  if (which=="hatval") { return(BLOB$hatval) } ##  REML hatval computation (also named hatval_ZX)
  if (which=="LevMar_step") {  ## LM with diagonal perturbation as in NocedalW p. 266
    #
    if (.spaMM.data$options$use_G_dG) { ## compared 01/2022 on adjacency-long with LevM=TRUE => faster than alternative. Small numerical diff too.
      # uses rather complex code to avoid solve(dense d2hdv2,...) but cannot avoid solve(rather dense G_dG)
      G_dG_blob <- .calc_G_dG(BLOB, damping) # list(G_dG=G_dG, dampDpD_2=damping * BLOB$D_Md2hdv2)
      ## sequel recomputed for each new damping value...
      grad_v <- LM_z$scaled_grad[seq(ncol(BLOB$chol_Q))] 
      G_dG <- G_dG_blob$G_dG
      L_dv_term_from_grad_v <- Matrix::drop(Matrix::solve(G_dG, Matrix::tcrossprod(BLOB$chol_Q, grad_v))) # part from grad_h, to be completed
      if (attr(sXaug,"pforpv")) { ## if there are fixed effect coefficients to estimate
        # rhs for beta
        grad_beta <- LM_z$scaled_grad[-seq(ncol(BLOB$chol_Q))]
        dbeta_rhs <- grad_beta - crossprod(BLOB$ZtWX, L_dv_term_from_grad_v) # - dampDpD_1 * LM_z$v_h_beta$beta_eta one damps the gradient eq, not the one giving v_h_beta
        # matrix for beta
        Xt_X <- crossprod(BLOB$ZtWX, as(Matrix::solve(G_dG, BLOB$ZtWX),"matrix"))
        XX_D <- BLOB$XtWX
        dampDpD_1 <- damping*diag(XX_D)
        diag(XX_D) <- diag(XX_D)+dampDpD_1 ## adds D1=damping*diag(XX_D)
        # see INLA_vs_spaMM_Loaloa.R for one case where the matrix is singular.
        dbeta_eta <- try(solve(XX_D-Xt_X , dbeta_rhs)[,1],silent=TRUE)
        if (inherits(dbeta_eta,"try-error")) dbeta_eta <- (ginv(XX_D-Xt_X) %*% dbeta_rhs)[,1]
        L_dv <- L_dv_term_from_grad_v - drop(solve(G_dG, BLOB$ZtWX %*% dbeta_eta)) # - dampDpD_2 %*% LM_z$v_h_beta$v_h
        dv_h <- Matrix::drop(Matrix::crossprod(BLOB$chol_Q, Matrix::drop(L_dv))) ## inner drop() to avoid signature issue with dgeMatrix dv_rhs...
      } else {
        dbeta_eta <- numeric(0)
        dv_h <- Matrix::drop(Matrix::crossprod(BLOB$chol_Q, Matrix::drop(L_dv_term_from_grad_v))) ## inner drop() to avoid signature issue with dgeMatrix dv_rhs...
      }
    } else { ## may be faster but with useQR set to FALSE 
      ## For each new damping value:
      grad_v <- LM_z$scaled_grad[seq(ncol(BLOB$chol_Q))]
      grad_beta <- LM_z$scaled_grad[-seq(ncol(BLOB$chol_Q))]
      if ((useQR <- FALSE)) {
        if (is.null(BLOB$D_Md2hdv2)) {
          ## t() bc .damping_to_solve implicitly uses (or used) crossprod ( as in solve(X'X)=solve(R'R) ) (or equivalently tcrossprod(solve))
          BLOB$Hfac <- tmp <- t(drop0(BLOB$tcrossfac_Md2hdv2)) ## probably not so sparse...
          xx <- tmp@x
          xx <- xx*xx
          tmp@x <- xx
          BLOB$D_Md2hdv2 <- colSums(tmp) ## colSums bc t() above // rowSums for diag Tcrossprod(invQ G)  = invQ G Gt invQt
        }
        dampDpD_2 <- damping * BLOB$D_Md2hdv2
        DS <- .damping_to_solve(X=BLOB$Hfac, dampDpD=dampDpD_2, rhs=NULL) ## potentially slow chol2inv()
        dv_term_from_grad_v <- drop((DS$inv %*% grad_v[DS$Rperm])[DS$RRsP]) 
      } else {
        H_dH_blob <- .calc_H_dH(BLOB, damping) # list(dVscaled= dv_h, dampDpD = dampDpD_2)
        dv_term_from_grad_v <- drop(solve(H_dH_blob$H_dH,  grad_v))
      }
      if (attr(sXaug,"pforpv")) { ## if there are fixed effect coefficients to estimate
        dbeta_rhs <- grad_beta - crossprod(BLOB$LZtWX, dv_term_from_grad_v) # - dampDpD_1 * LM_z$v_h_beta$beta_eta one damps the gradient eq, not the one giving v_h_beta
        if (useQR) {
          thinsolve <- as((DS$inv %*% BLOB$LZtWX[DS$Rperm,])[DS$RRsP,],"matrix") ## something wrong here ?
        } else thinsolve <- solve(H_dH_blob$H_dH, BLOB$LZtWX) ## thin result, needed for Xt_X
        Xt_X <- crossprod(BLOB$LZtWX, thinsolve)
        XX_D <- BLOB$XtWX
        diag(XX_D) <- diag(XX_D)*(1+damping) ## adds D1=damping*diag(XX_D)
        #
        dbeta_eta <- try(solve(XX_D-Xt_X , dbeta_rhs)[,1],silent=TRUE)
        if (inherits(dbeta_eta,"try-error")) dbeta_eta <- (ginv(XX_D-Xt_X) %*% dbeta_rhs)[,1]
        dv_h <- dv_term_from_grad_v - drop(thinsolve %*% dbeta_eta) 
      } else {
        dbeta_eta <- numeric(0)
        dv_h <- dv_term_from_grad_v
      }
      #invertand <- solve(BLOB$Gmat)+ (.ZtWZwrapper(solve(BLOB$chol_Q),invdH)) ## attendion pas matrix diag en dern pos
      
      # range(BLOB$Gmat+BLOB$chol_Q %*% Diagonal(x=damping * BLOB$D_Md2hdv2) %*% t(BLOB$chol_Q)-G_dG)
      # # 2 representations de H+dH
      # range(solve(BLOB$chol_Q, BLOB$Gmat %*% t(solve(BLOB$chol_Q)))+Diagonal(x=damping * BLOB$D_Md2hdv2)-solve(BLOB$chol_Q, G_dG %*% t(solve(BLOB$chol_Q))))
      # # 2 repres de solve(H+dH)
      # range(solve(solve(BLOB$chol_Q) %*%  BLOB$Gmat %*% t(solve(BLOB$chol_Q))+Diagonal(x=damping * BLOB$D_Md2hdv2))-
      #         t(BLOB$chol_Q) %*% solve(G_dG, BLOB$chol_Q))
      # # Woodbury 
      # invWW <- Diagonal(x=1/(damping * BLOB$D_Md2hdv2))
      # UU <- solve(BLOB$chol_Q)
      # VV <- t(UU)
      # HH <- UU %*%  BLOB$Gmat %*% VV
      # 
      # range(solve(UU %*%  BLOB$Gmat %*% VV+Diagonal(x=damping * BLOB$D_Md2hdv2))-
      #         t(BLOB$chol_Q) %*% solve(G_dG, BLOB$chol_Q))
      # range(invWW-invWW %*% UU %*% solve(solve(BLOB$Gmat)+ (t(solve(BLOB$chol_Q)) %*% invWW %*% solve(BLOB$chol_Q))) %*% VV %*% invWW -
      #         t(BLOB$chol_Q) %*% solve(G_dG, BLOB$chol_Q))
      # range(invWW-invWW %*% UU %*% solve(solve(BLOB$Gmat)+ (.ZtWZwrapper(solve(BLOB$chol_Q),1/(damping * BLOB$D_Md2hdv2)))) %*% VV %*% invWW -
      #         t(BLOB$chol_Q) %*% solve(G_dG, BLOB$chol_Q))
      # range(solve(t(BLOB$chol_Q)) %*% solve(HH+Diagonal(x=damping * BLOB$D_Md2hdv2)) %*% solve(BLOB$chol_Q)-
      #         solve(G_dG))
      
    }
    return(list(dVscaled=dv_h, dbeta_eta=dbeta_eta, dampDpD = damping*BLOB$DpD))
  }
  if (which =="LevMar_step_v_h") {  ## LM with diagonal perturbation as in NocedalW p. 266
    ## Here I deleted from v2.4.100 a variant using .damping_to_solve()
    grad_v <- LM_z$scaled_grad[seq(ncol(BLOB$chol_Q))]
    if (.spaMM.data$options$use_G_dG) {
      G_dG_blob <- .calc_G_dG(BLOB, damping) # list(G_dG=G_dG, dampDpD_2=damping * BLOB$D_Md2hdv2)
      rhs <- BLOB$chol_Q %*% grad_v
      rhs <- Matrix::solve(G_dG_blob$G_dG,rhs)
      dv_h <- drop(crossprod(BLOB$chol_Q,rhs))
      resu <- list(dVscaled= dv_h, dampDpD = G_dG_blob$dampDpD_2)
    } else { ## slower... sigh. Maybe bc the code for (full) LevMar_step uses the other approach.
      H_dH_blob <- .calc_H_dH(BLOB, damping) # list(dVscaled= dv_h, dampDpD = dampDpD_2)
      dv_h <- Matrix::solve(H_dH_blob$H_dH,  grad_v)[,1]
      resu <- list(dVscaled= dv_h, dampDpD = H_dH_blob$dampDpD_2)
    }
    return(resu)
  }
  
  if (which =="LevMar_step_beta") {
    if (attr(sXaug,"pforpv")) { ## if there are fixed effect coefficients to estimate
      dbeta_rhs <- LM_z$scaled_grad[-seq(ncol(BLOB$chol_Q))]
      XX_D <- BLOB$XtWX
      dampDpD <- damping*diag(XX_D)
      diag(XX_D) <- diag(XX_D)+dampDpD
      dbeta_eta <- try(solve(XX_D , dbeta_rhs)[,1],silent=TRUE)
      if (inherits(dbeta_eta,"try-error")) dbeta_eta <- (ginv(XX_D) %*% dbeta_rhs)[,1]
      resu <- list(dbeta_eta=dbeta_eta, dampDpD =dampDpD)
      return(resu)
    } else stop("LevMar_step_beta called with 0-length rhs: pforpv=0?")
  }
  if (which=="beta_cov_info_from_wAugX") { ## for predVar 
    ## not called during the fit, so ideally we store the necessary info in the fit rather than use get_from_MME() 
    stop("Programming error. Use ad hoc .calc_beta_cov_info_spprec() function instead.")
  }
  stop(cat(paste0("'which=\"",which,"\"' is invalid.")))
}

get_from_MME.AUGI0_ZX_sparsePrecision <- function(sXaug, 
                                                  which="",szAug=NULL,B=NULL,
                                                  damping=NULL, LM_z=NULL, ...) {
  if (which=="LevMar_step" && damping==0L) {
    resu <- get_from_MME.AUGI0_ZX_sparsePrecision(sXaug=sXaug, szAug=LM_z)
  } else resu <- .AUGI0_ZX_sparsePrecision(sXaug=sXaug,which=which,z=szAug,B=B,damping=damping,LM_z=LM_z)
  resu
}
