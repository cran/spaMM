def_AUGI0_ZX_sparsePrecision <- function(AUGI0_ZX, corrPars,w.ranef,cum_n_u_h,w.resid 
                                         #,weight_X ## currently ignored
                                         #,H_global_scale ## currently ignored
                                         ) {
  resu <- structure(list(AUGI0_ZX=AUGI0_ZX,
                         BLOB=list2env(list(), parent=environment(.AUGI0_ZX_sparsePrecision))),
               w.ranef=w.ranef,
               cum_n_u_h=cum_n_u_h,
               pforpv=ncol(AUGI0_ZX$X.pv),
               w.resid= if (is.list(w.resid)) {w.resid$w_resid} else w.resid,
               corrPars=corrPars )
  class(resu) <- c("AUGI0_ZX_sparsePrecision","list") # (## )do not define recursively if object is an envir...)
  return( resu ) 
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
  BLOB$chol_Q <- as(Matrix::.bdiag(chol_Q_list),"dtCMatrix") # !not dtC by default! (!). Has worked as (triangular) dgCMatrix for a long time.
  return(precisionBlocks)
}

..calc_ZtWX <- function(AUGI0_ZX, w.resid, X.pv=AUGI0_ZX$X.pv) {
  if (attr(w.resid,"is_unit")) {
    ZtWX <- as.matrix(crossprod(AUGI0_ZX$ZAfix, X.pv)) 
  } else {
    if (methods::.hasSlot(AUGI0_ZX$ZAfix, "x") && # 1st condition is true in all tests. But not clearly enforced. (__F I X M E__ ?) 
        length(AUGI0_ZX$ZAfix@x)>prod(dim(X.pv))) {
      ZtWX <- .Dvec_times_m_Matrix(w.resid, X.pv) 
      ZtWX <- as.matrix(crossprod(AUGI0_ZX$ZAfix, ZtWX)) 
    } else {
      ZtWX <- .Dvec_times_m_Matrix(w.resid, AUGI0_ZX$ZAfix) 
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
      ZtWX[[it]] <- ..calc_ZtWX(AUGI0_ZX, w.resid=attr(sXaug,"w.resid"), X.pv=AUGI0_ZX$X.pv[,blocs[[it]],drop=FALSE])
    }
    ZtWX <- do.call(cbind, ZtWX)
  } else ZtWX <- ..calc_ZtWX(AUGI0_ZX, w.resid=attr(sXaug,"w.resid"))
  ZtWX ## always *m*atrix
}

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
    sqrtWX <- as.matrix(.Dvec_times_m_Matrix( sqrt(attr(sXaug,"w.resid")),AUGI0_ZX$X.pv))
    seqncol <- seq(pforpv)
    row_blocs <- split(seqncol, ceiling(seq_along(seqncol)/bloclength))
    crossprod_r22 <- vector("list",length(row_blocs))
    for (it in seq_along(row_blocs)) {
      crossprod_r22[[it]] <- .crossprod(BLOB$r12, BLOB$r12[,row_blocs[[it]],drop=FALSE])
      if(! is.null(BLOB$G_scaling)) crossprod_r22[[it]] <- crossprod_r22[[it]]/BLOB$G_scaling
      crossprod_r22[[it]] <- .crossprod(sqrtWX, sqrtWX[,row_blocs[[it]],drop=FALSE]) - crossprod_r22[[it]]
    }
    crossprod_r22 <- do.call(cbind, crossprod_r22)
  } else {
    crossprod_r12 <- crossprod(as.matrix(BLOB$r12))  # Xt_X ! # using most precise computation as in .calc_r22()
    crossprod_r22 <- .ZtWZwrapper(AUGI0_ZX$X.pv,attr(sXaug,"w.resid")) - crossprod_r12 ## XtWX-Xt_X ! Both lines as explained in working doc
  }
  return(crossprod_r22) # dgeMatrix
}

.calc_r22 <- function(X.pv, w.resid, r12, XtX=NULL) { 
  ## i's vain to try to regularize crossr22 by manipulating diagonal or eigenvalues. 
  #  It would be better to control the accuracy of the computation of crossr22 as a difference.
  if (ncol(X.pv)) { 
    if (attr(w.resid,"is_unit") && ! is.null(XtX)) { # 2nd test needed bc currently in post-fit calls, XtX remains NULL (__F I X M E__ ?)
      crossr22 <- XtX - crossprod(as.matrix(r12),NULL,as_mat=TRUE) 
    } else crossr22 <- as.matrix(.ZtWZwrapper(X.pv,w.resid))-crossprod(as.matrix(r12),NULL,as_mat=TRUE)# Make sure that '-' operates on *m*atrices 
    # (1)
    # as.matrix(.ZtWZwrapper(X.pv,w.resid)) is then a safety for rare case where X.pv is sparse
    # (2)
    # .Rcpp_crossprod(r12,NULL,as_mat=TRUE) may be faster but less accurate numerically. Numerical precision important here.
    # + see https://stackoverflow.com/questions/52888650/eigenselfadjointviewrankupdate-slower-than-a-ww-transpose
    # A test is in test-mv-nested: logLik(mrf1T <- fitme(migStatus ~ 1 +multIMRF(1... should give phi=1e-6 and p_v=-11.24422
    # 
    # Different numerics for the same pb: 
    #  .calc_inv_beta_cov(sXaug) computes crossr22 from sXaug, handling some complication, 
    # .get_absdiagR_blocks() face similar pbs and finally calls .Utri_chol_by_qr()
    return(.wrap_Utri_chol(crossr22)) # rather than (Matrix::)chol()
  } else return(diag(nrow=0L))
  ## lots of alternatives here, removed from [v2.6.54
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
  Zt_sqrtw_pwy_o <- crossprod(AUGI0_ZX$ZAfix, sqrt(attr(sXaug,"w.resid"))*pwy_o) 
  Xt_sqrtw_pwy_o <- crossprod(AUGI0_ZX$X.pv, sqrt(attr(sXaug,"w.resid"))*pwy_o)
  if (.is_evaluated("invL_G.P", BLOB)) { 
    lev_phi_z_pwy_o <- BLOB$invL_G.P %*% Zt_sqrtw_pwy_o
  } else lev_phi_z_pwy_o <- Matrix::solve(BLOB$G_CHMfactor, 
                                          as(BLOB$G_CHMfactor,"pMatrix") %*% Zt_sqrtw_pwy_o, system="L") ## R_11^{-T}.Ztw.pwy_o
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
  qrQ <- qr.Q(BLOB$qrXa)
  qrQ@x <- qrQ@x^2
  tmp <- rowSums(qrQ)
  n_u_h <- ncol(AUGI0_ZX$I)
  phipos <- n_u_h+seq_len(length(tmp)-n_u_h)
  list(lev_lambda=tmp[-phipos],lev_phi=tmp[phipos]) # BLOB$hatval , formely hatval_ZX, for REML
}

.calc_spprec_hatval_ZX_by_r22 <- function(BLOB, sXaug, AUGI0_ZX, w.ranef) {
  # We use the  Chol factorization of T'T as explained in the working doc
  # We know one of its block is the same as for the ML leverages, we compute two additional blocks: 
  lev_lambda_z <- BLOB$factor_inv_Md2hdv2 ## ul block of R^{-T} as described in working doc
  lev_lambda_z@x <- lev_lambda_z@x^2 
  lev_lambda_z <- .Matrix_times_Dvec(lev_lambda_z,w.ranef)
  lev_lambda_z <- colSums(lev_lambda_z)
  #
  if (attr(sXaug,"pforpv")>0L) {
    lev_lambda_x <- - backsolve(BLOB$r22,crossprod(BLOB$r12, BLOB$factor_inv_Md2hdv2),transpose=TRUE) ## ll block of R^{-T} as described in working doc
    lev_lambda_x <- lev_lambda_x^2
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
    sqrtwZ <- .Dvec_times_Matrix(sqrt(attr(sXaug,"w.resid")), AUGI0_ZX$ZAfix)
    crossprod_r12_z_rows <- crossprod(crossprod(BLOB$invL_G.P, BLOB$r12),t(sqrtwZ)) 
    lev_phi_z <- colSums(BLOB$invL_G.P^2)
    lev_phi_z <- drop(AUGI0_ZX$ZAfix %*% lev_phi_z)*attr(sXaug,"w.resid")
  } else {
    lev_phi_z <- BLOB$inv_L_G_ZtsqrW
    crossprod_r12_z_rows <- crossprod(BLOB$r12,lev_phi_z)
    lev_phi_z@x <- lev_phi_z@x^2
    lev_phi_z <- colSums(lev_phi_z) 
  }
  if (attr(sXaug,"pforpv")>0L) {
    XtW <- t(AUGI0_ZX$X.pv *sqrt(attr(sXaug,"w.resid")))
    lev_phi_x <- backsolve(BLOB$r22, XtW - crossprod_r12_z_rows, transpose=TRUE)
    lev_phi_x <- lev_phi_x^2
    if ( is.matrix(lev_phi_x)) lev_phi_x <- colSums(lev_phi_x) 
    lev_phi <- lev_phi_z+lev_phi_x
    list(lev_lambda=lev_lambda,lev_phi=lev_phi) # BLOB$hatval , formely hatval_ZX, for REML
  } else list(lev_lambda=lev_lambda_z,lev_phi=lev_phi_z) # BLOB$hatval , formely hatval_ZX, for REML
}


.Sigsolve_sparsePrecision <- function(sXaug, rhs) { ## no longer used: compare .calc_inv_beta_cov() when rhs= X.pv. But useful as doc.
  if (is.null(sXaug$BLOB$ZtW)) sXaug$BLOB$ZtW <- t(.Dvec_times_m_Matrix(attr(sXaug,"w.resid"), sXaug$AUGI0_ZX$ZAfix)) 
  v <- Matrix::solve(sXaug$BLOB$G_CHMfactor, sXaug$BLOB$ZtW %*% rhs)
  # implicit .Dvec_times_m_Matrix( on a dgeMatrix:
  invM.z <- attr(sXaug$AUGI0_ZX,"w.resid") * (rhs - sXaug$AUGI0_ZX$ZAfix %*% v)  # W.rhs - W Z (v= L inv(...) L' Z' W rhs) = [W- WZL inv(...) L'Z'W] rhs
  return(invM.z) ## where M=Z.solve(precMat).Zt+diag(~phi)
} # solve(AUGI0_ZX$ZAfix %*% solve(precisionMatrix) %*% t(AUGI0_ZX$ZAfix)+diag(1/attr(sXaug,"w.resid"))) = Sigsolve(diag(x=ncol(ZtW)))

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
  # so having previously forced symmetry on dG is a waste of time=> as_sym=FALSE
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
        tmp@x <- tmp@x^2
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
    if (as_sym) {
      G_dG <- .dsCsum(BLOB$Gmat, dampdG) # BLOB$Gmat + dampdG ## probably not so sparse... 
      #G_dG <- .dsC_plus_dsC(BLOB$Gmat,dampdG) ##
    } else G_dG <- forceSymmetric(BLOB$Gmat + dampdG)
    return(list(G_dG=G_dG, dampDpD_2=damping * BLOB$D_Md2hdv2))
  }
}

.provide_BLOB_hatval_Z_ <- function(BLOB, w.ranef, AUGI0_ZX, sXaug, needed=c("lambda","phi")) {
  if ("lambda" %in% needed) {
    if (is.null(BLOB$hatval_Z_$lev_lambda)) {
      lev_lambda <- BLOB$factor_inv_Md2hdv2
      lev_lambda@x <- lev_lambda@x^2 
      BLOB$hatval_Z_$lev_lambda <- colSums(lev_lambda) * w.ranef
    }
  }
  if ("phi" %in% needed) {
    if (is.null(BLOB$hatval_Z_$lev_phi)) {
      if (AUGI0_ZX$is_unitary_ZAfix) { ## then only diagonal values of invG matter  ## adjacency-long case...?
        ## invL_G.P presumably needed for lev_lambda from factor_inv_Md2hdv2
        ## (...lev_phi from invL_G.P, lev_lambda from invL_G.P %*% chol_Q ) 
        lev_phi <- colSums(BLOB$invL_G.P^2)
        BLOB$hatval_Z_$lev_phi <- drop(AUGI0_ZX$ZAfix %*% lev_phi)*attr(sXaug,"w.resid")
      } else {
        # inv_L_G_ZtsqrW arises as simplif of BLOB$factor_inv_Md2hdv2 %*% t(sqrtwZL) (L and chol_Q cancel each other)
        lev_phi <- BLOB$inv_L_G_ZtsqrW ## solve(as(BLOB$G_CHMfactor,"sparseMatrix"), t(sqrtwZ))
        lev_phi@x <- lev_phi@x^2
        BLOB$hatval_Z_$lev_phi <- colSums(lev_phi)
      }
    }
  }
  BLOB$hatval_Z_ 
}

.provide_BLOB_factor_inv_Md2hdv2 <- function(BLOB) { ## Code encapsulated in a function for easier profiling
  # .calc_sscaled_new() -> (Pdiag <- get_from_MME(sXaug,"hatval_Z")) -> appears to be the last big bottleneck for fitme(). 
  # We need a diag(crossprod(factor)) =colSums(factor$x^2) needing the full factor
  LL <- BLOB$invL_G.P %*% BLOB$chol_Q ## dgCMatrix 
  # Even is BLOB$chol_Q is Identity, this this hardly slower than Matrix::solve(BLOB$G_CHMfactor, system="L") so cannot be improved.
  # The pattern Matrix::solve(BLOB$G_CHMfactor, as(BLOB$G_CHMfactor,"pMatrix") %*% ., system="L")  occurs repeatedly and is costly
  # so the Q is whether to store the (rel dense) invL_G.P= solve(BLOB$G_CHMfactor, as(BLOB$G_CHMfactor,"pMatrix"))
  # Currently it is computed and stored only for several leverage computation and for factor_inv_Md2hdv2
  #factor_Md2hdv2 <- as.matrix(Matrix::solve(BLOB$chol_Q,as(BLOB$G_CHMfactor,"sparseMatrix")))
  #BLOB$chol_Md2hdv2 <- Rcpp_chol_R(tcrossprod(factor_Md2hdv2)) #A = R'.R as in R's chol()
  factor_inv_Md2hdv2 <- drop0(LL,tol=.Machine$double.eps) ## crossfactor
} # *returns* the factor

.factor_inv_Md2hdv2_times_rhs <- function(BLOB, B) { 
  ## get_from_MME(., "solve_d2hdv2", .) => B generally vector but may be diag, Diagonal 
  ## cf HLfit_body -> .calc_lev_from_hat -> -> .calc_dvdloglamMat_new()
  ## or even fuller matrices (neg.d2h0_dv_dlogphi <- .m_Matrix_times_Dvec(t(ZAL), drop(dh0deta)) in .calc_dvdlogphiMat_new())
  if (.is_evaluated("factor_inv_Md2hdv2", BLOB)) {
    BLOB$factor_inv_Md2hdv2 %*% B # __F I X M E__ optimize for diagonal matrices? Not easy (info to be passed efficiently from get_from_MME call)
  } else {
    B <- BLOB$chol_Q %*% B
    B <- B[BLOB$perm, ] # as(BLOB$G_CHMfactor,"pMatrix") %*% vec
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

.init_promises_spprec <- function(sXaug, non_X_ones=TRUE, nullify_X_ones =FALSE) {  
  BLOB <- sXaug$BLOB # environment that should contain  G_CHMfactor, chol_Q, perm when promises are evaluated
  AUGI0_ZX <- sXaug$AUGI0_ZX # envir (prime fit) or list (converted to list by confint)
  # 
  if ( non_X_ones ) {
    ### What does NOT depend on X
    delayedAssign("sortPerm", sort.list(BLOB$perm), assign.env = BLOB ) # never NULL
    delayedAssign("invL_G.P", drop0(Matrix::solve(BLOB$G_CHMfactor, as(BLOB$G_CHMfactor,"pMatrix"), system="L")), assign.env = BLOB ) # never NULL
    delayedAssign("inv_L_G_ZtsqrW", {
      sqrtwZ <- .Dvec_times_Matrix(sqrt(attr(sXaug,"w.resid")), AUGI0_ZX$ZAfix)
      inv_L_G_ZtsqrW <- tcrossprod(BLOB$invL_G.P, sqrtwZ)
    }, assign.env = BLOB )
    delayedAssign("factor_inv_Md2hdv2",.provide_BLOB_factor_inv_Md2hdv2(BLOB), assign.env = BLOB )
    #delayedAssign("half_logdetQ",  sum(log(diag(x=BLOB$chol_Q))), assign.env = BLOB ) ## currently not used (in variant of LevMar step)
    #delayedAssign("half_logdetG", Matrix::determinant(BLOB$G_CHMfactor)$modulus[1], assign.env = BLOB ) ## currently not used (in variant of LevMar step)
    delayedAssign("logdet_sqrt_d2hdv2", { 
      half_logdetG <- Matrix::determinant(BLOB$G_CHMfactor)$modulus[1]
      half_logdetQ <- sum(log(diag(x=BLOB$chol_Q)))
      half_logdetG - half_logdetQ 
    }, assign.env = BLOB )
    delayedAssign("logdet_R_scaled_v", { BLOB$logdet_sqrt_d2hdv2 - sum(log(attr(sXaug,"w.ranef")))/2 }, assign.env = BLOB )
    # more or less for speculative code:
    delayedAssign("Md2hdv2", .tcrossprod(BLOB$tcrossfac_Md2hdv2), assign.env = BLOB ) ## currently not used (for H_dH)
    delayedAssign("tcrossfac_Md2hdv2",
                  Matrix::solve(BLOB$chol_Q, crossprod(as(BLOB$G_CHMfactor,"pMatrix"), as(BLOB$G_CHMfactor,"sparseMatrix"))), assign.env = BLOB )
  }
  #
  ### What DOES depend on X (for confint code: lost columns; and predVar computations currently all based on unscaled beta's)
  # One exception is qrXa, bc the post-fit code using it takes the rescaling into account, so we treat it as the non-X ones 
  if ( non_X_ones ) {
    delayedAssign("qrXa", { ## used for "hatval" in case of failure of computation of BLOB$r22
      X.pv <- as(AUGI0_ZX$X.pv,"dgCMatrix")#  .Rcpp_as_dgCMatrix(AUGI0_ZX$X.pv)
      Xscal <- with(sXaug,rbind(cbind(AUGI0_ZX$I, AUGI0_ZX$ZeroBlock), # I0_ZX order
                                cbind(AUGI0_ZX$ZAfix %*% t(solve(BLOB$chol_Q)), X.pv)))
      Xscal <- .Dvec_times_Matrix(sqrt(c(attr(sXaug,"w.ranef"),attr(sXaug,"w.resid"))),Xscal)
      qrXa <- qr(Xscal)
    }, assign.env = BLOB )
  }
  if (nullify_X_ones) {
    BLOB$ZtWX <- BLOB$XtX <- BLOB$XtWX <- BLOB$r12 <- BLOB$r22 <- BLOB$DpD <- BLOB$LZtWX <- BLOB$crossr22 <- NULL
  } else {
    delayedAssign("ZtWX", .calc_ZtWX(sXaug), assign.env = BLOB ) # *m*atrix
    delayedAssign("XtX", crossprod(as.matrix(AUGI0_ZX$X.pv), NULL, as_mat=TRUE), assign.env=BLOB)
    delayedAssign("XtWX", .ZtWZwrapper(AUGI0_ZX$X.pv,attr(sXaug,"w.resid")), assign.env = BLOB )  ## as(,"matrix") ?
    delayedAssign("r12", {
      # Use inv_L_G_ZtsqrW when it's available (as needed for certain leverage computations) 
      # Otherwise use ZtWX which is typically much smaller, so that solve(., ZtWX) is faster than solve(., ZtsqrW)
      if (.is_evaluated("inv_L_G_ZtsqrW", BLOB)) { 
        sqrtwX <- .Dvec_times_m_Matrix(sqrt(attr(sXaug,"w.resid")),AUGI0_ZX$X.pv)
        r12 <- BLOB$inv_L_G_ZtsqrW %*% sqrtwX
      } else if (.is_evaluated("invL_G.P", BLOB)) { 
        r12 <- BLOB$invL_G.P %*% BLOB$ZtWX
      } else r12 <-  Matrix::solve(BLOB$G_CHMfactor, as(BLOB$G_CHMfactor,"pMatrix") %*% BLOB$ZtWX, system="L")
      # : not computing invL_G.P here is big improvement for test IMRF rawGNIPdataEU.1, and also for pcmatern LMM
    }, assign.env = BLOB )
    delayedAssign("r22", .calc_r22(AUGI0_ZX$X.pv,attr(sXaug,"w.resid"),BLOB$r12, XtX=BLOB$XtX), assign.env = BLOB )
    delayedAssign("DpD", c(BLOB$D_Md2hdv2,diag(x=BLOB$XtWX)), assign.env = BLOB )
    #
    # more or less for speculative code:
    delayedAssign("LZtWX", as(solve(BLOB$chol_Q, BLOB$ZtWX),"matrix"), assign.env = BLOB ) ## currently not used (in variant of LevMar step)
    delayedAssign("crossr22", { # NOT currently used: Used only in-fit -> .calc_sum_pwt_Q_y_o_2() ->  .solve_crossr22(., use_crossr22=TRUE), but use_crossr22 set to FALSE:
      if (ncol(X.pv)) { 
        w.resid <- attr(sXaug,"w.resid")
        if (attr(w.resid,"is_unit")) { 
          crossr22 <- BLOB$XtX - crossprod(as.matrix(BLOB$r12),NULL,as_mat=TRUE) 
        } else crossr22 <- as.matrix(.ZtWZwrapper(AUGI0_ZX$X.pv,w.resid))-crossprod(as.matrix(BLOB$r12),NULL,as_mat=TRUE)# Make sure that '-' operates on *m*atrices 
      } else crossr22 <- diag(nrow=0L)
    } , assign.env = BLOB )
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
    delayedAssign( "logdet_R_scaled_b_v", { BLOB$logdet_R_scaled_v + BLOB$logdet_r22} , assign.env = BLOB )
    delayedAssign( "hatval", { # hatval_ZX
      if (identical(.spaMM.data$options$use_spprec_QR,TRUE) || is.null(dim(BLOB$r22))) { # ((currently FALSE) || failure of comput of r22 )
        .calc_spprec_hatval_ZX_by_QR(BLOB=BLOB, AUGI0_ZX=AUGI0_ZX, sXaug=sXaug, w.ranef=attr(sXaug,"w.ranef"))
      } else .calc_spprec_hatval_ZX_by_r22(BLOB=BLOB, AUGI0_ZX=AUGI0_ZX, sXaug=sXaug, w.ranef=attr(sXaug,"w.ranef")) ## fragile bc  BLOB$r22 <- .calc_r22() is fragile...
    } , assign.env = BLOB )
  } # end block of promises depending on X
}


.AUGI0_ZX_sparsePrecision <- function(sXaug,which="",z=NULL,B=NULL,
                                      damping,LM_z) {
  AUGI0_ZX <- sXaug$AUGI0_ZX 
  BLOB <- sXaug$BLOB ## an environment defined by def_AUGI0_ZX_sparsePrecision
  w.ranef <- attr(sXaug,"w.ranef") 
  cum_n_u_h <- attr(sXaug,"cum_n_u_h")
  if (is.null(BLOB$G_CHMfactor)) { ## costly, cannot be precomputed
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
      precisionMatrix <- forceSymmetric(Matrix::bdiag(precisionBlocks)) # bdiag degrades dsC into dgC
    } else precisionMatrix <- (precisionBlocks[[1L]]) # forceSymmetric removed bc input should in principle be dsC
    tmp <- .ZtWZwrapper(AUGI0_ZX$ZAfix,attr(sXaug,"w.resid")) # should return a symmetric-type matrix (dsC or possibly dsy, not ddi...)
    if (inherits(tmp,"dsyMatrix")) { 
      # if tmp is small (2*2 and diagonal: ranef with two levels...), then it may have been returned by .ZtWZwrapper -> .crossprod as dsy, not dsC.
      # then, as(tmp,"sparseMatrix") is dgC not dsC
      #  and  as(tmp,"symmetricMatrix") is dsy not dsC
      # so there is no way to make sure that the result is dsC without adding as(.,"dsCMatrix") in .ZtWZwrapper or .crossprod
      .ad_hoc_dsy_warning()
      tmp <- as(tmp,"dsCMatrix")
    }
    if (inherits(tmp,"dsCMatrix") && inherits(precisionMatrix,"dsCMatrix") ) { # test introduced 02/2020. 
      tmp <- .dsCsum(tmp, precisionMatrix) # faster than tmp + precisionMatrix where '+' calls forceSymmetric(callGeneric(as(e1, "dgCMatrix"), as(e2, "dgCMatrix")))
    } else { 
      stop("Unexpected matrix types in .AUGI0_ZX_sparsePrecision(). Please report to the package maintainer.") 
      # actually there is some other step that is silently wrong if both matrices are not dsC...
      tmp <- tmp + precisionMatrix
    }
    BLOB$Gmat <- drop0(tmp) ## depends on w.ranef and w.resid
    # if (inherits(BLOB$Gmat,"dgCMatrix")) stop("ICI") 
    if (is.null(template <- AUGI0_ZX$template_G_CHM)) { ## occurs if $update_CHM is FALSE, OR first comput. of CHM, OR precision factors not yet all $updateable
      BLOB$G_CHMfactor <- Cholesky(BLOB$Gmat,LDL=FALSE,perm=.spaMM.data$options$perm_G) ## costly
      if (.spaMM.data$options$update_CHM && all(AUGI0_ZX$envir$updateable)) AUGI0_ZX$template_G_CHM <- BLOB$G_CHMfactor
    } else BLOB$G_CHMfactor <- Matrix::.updateCHMfactor(template, BLOB$Gmat, mult=0) ## If it fails in spde, look for too low kappa.
    ## with perm=TRUE G=P'LL'P and the P'L (non-triangular) factor is given by solve(<G_CHM>,as(<G_CHM>,"sparseMatrix"),system="Pt")
    # (?__F I X M E__?) It's theo. possible to call .updateCHMfactor on the tcrossfac of Gmat...
    BLOB$perm <- as(BLOB$G_CHMfactor, "pMatrix")@perm # from 1L
    .init_promises_spprec(sXaug)
  }
  if (which=="initialize") return(NULL)
  if (which=="hatval_Z") {
    ## calcul correct; 
    # TT <- rbind(diag(sqrt(w.ranef)),diag(sqrt(attr(AUGI0_ZX, "w.resid"))) %*% AUGI0_ZX$ZAfix %*% t(solve(BLOB$chol_Q)) )
    # diag(TT %*% (solve(crossprod(TT))) %*% t(TT)) ## lev_lambda, lev_phi
    if (is.null(B)) B <- c("lambda","phi")
    return(.provide_BLOB_hatval_Z_(BLOB, w.ranef, AUGI0_ZX, sXaug, needed=B))
  }
  if (which=="Mg_solve_g") {
    stop("Mg_solve_g not implemented") # .unlist(dVscaled_beta)*zInfo$m_grad_obj should be OK
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
      sol <-  - .crossprod(BLOB$factor_inv_Md2hdv2)
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
      dbeta_rhs <- grad_p_v - .crossprod(BLOB$r12, rhs)[,1L] # numeric  vector
      dbeta_eta <-  .solve_crossr22(BLOB, AUGI0_ZX, sXaug, dbeta_rhs)
      rhs <- rhs - (BLOB$r12 %*% dbeta_eta)[,1]
    } else {
      dbeta_eta <- numeric(0)
    }
    dv_h <- .crossprod_factor_inv_Md2hdv2_rhs(BLOB, rhs)[,1L]
    return(list(dv_h=dv_h,dbeta_eta=dbeta_eta)) 
  }
  if (which=="logdet_sqrt_d2hdv2") { return(BLOB$logdet_sqrt_d2hdv2)} 
  #### REML:
  if (which=="logdet_r22") { return(BLOB$logdet_r22) }
  if (which=="logdet_R_scaled_b_v") { return(BLOB$logdet_R_scaled_b_v)} 
  if (which=="hatval") { return(BLOB$hatval) } ##  REML hatval computation (also named hatval_ZX)
  if (which %in% c("LevMar_step")) {  ## LM with diagonal perturbation as in NocedalW p. 266
    #
    if (.spaMM.data$options$use_G_dG) { ## not so clear which is faster ## big-ranefs would be a good test F I X M E
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
          BLOB$Hfac <- tmp <- t(drop0(solve(BLOB$chol_Q,  
                                            crossprod(as(BLOB$G_CHMfactor,"pMatrix"), as(BLOB$G_CHMfactor,"sparseMatrix"))))) ## probably not so sparse...
          tmp@x <- tmp@x^2
          BLOB$D_Md2hdv2 <- colSums(tmp) ## colSums bc t() above // rowSums for diag Tcrossprod(invQ G)  = invQ G Gt invQt
        }
        dampDpD_2 <- damping * BLOB$D_Md2hdv2
        DS <- .damping_to_solve(X=BLOB$Hfac, dampDpD=dampDpD_2, rhs=NULL) ## potentially slow chol2inv()
        dv_term_from_grad_v <- drop((DS$inv %*% grad_v[DS$Rperm])[DS$RRsP]) 
      } else {
        H_dH_blob <- .calc_H_dH(BLOB, damping) # list(dVscaled= dv_h, dampDpD = dampDpD_2)
        dv_term_from_grad_v <- solve(H_dH_blob$H_dH,  grad_v)
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
        dv_h <- dv_term_from_grad_v - thinsolve %*% dbeta_eta 
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

if (FALSE) {
  trace(".AUGI0_ZX_sparsePrecision", where=asNamespace("spaMM"), print=FALSE, tracer=quote(print(which)),
        exit=quote(if ( ! is.null(z)) {str(dv_h)} else switch(which,
                                                              "hatval_Z"=str(BLOB$hatval_Z_),
                                                              "solve_d2hdv2"= str(sol), 
                                                              "logdet_r22"=print(logdet_r22,digits = 12),
                                                              "half_logdetQ"=print(half_logdetQ,digits = 12),
                                                              "half_logdetG"=print(half_logdetG,digits = 12),
                                                              "beta_cov"={str(beta_cov);print(eigen(beta_cov,symmetric = TRUE,only.values = TRUE)$values)}
                                                              ) ) )
}
# untrace(".AUGI0_ZX_sparsePrecision", where=asNamespace("spaMM"))

get_from_MME.AUGI0_ZX_sparsePrecision <- function(sXaug, 
                                                  which="",szAug=NULL,B=NULL,
                                                  damping=NULL, LM_z=NULL, ...) {
  if (which=="LevMar_step" && damping==0L) {
    resu <- get_from_MME.AUGI0_ZX_sparsePrecision(sXaug=sXaug, szAug=LM_z)
  } else resu <- .AUGI0_ZX_sparsePrecision(sXaug=sXaug,which=which,z=szAug,B=B,damping=damping,LM_z=LM_z)
  resu
}
