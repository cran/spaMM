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
  # => All of chol_Q_list must be dtCMatrix so that bdiag() gives a dtCMatrix
  precisionBlocks <- envir$precisionFactorList ## local copy so that AUGI0_ZX$envir$precisionFactorList is not modified
  chol_Q_list <- vector("list",length(envir$precisionFactorList)) ## of L in LL' factorisation of precision factor 
  for (rd in seq_len(length(envir$precisionFactorList))) {
    if( ! is.null(envir$precisionFactorList[[rd]])) {
      if (inherits(envir$precisionFactorList[[rd]],c("Matrix","matrix"))) { ## Q_CHMfactor NOT provided
        ## not clear that adjacency or AR1 can occur here
        if (envir$finertypes[rd] %in% c("AR1", "ranCoefs", "adjacency")) { # adjacency-specific code removed [from version 2.5.30
          ## envir$precisionFactorList[[rd]] should be a list hence this block should not be reached
          stop("reached here in .reformat_Qmat_info()") 
        } else {## other models where a single matrix is provided : Qmatrix provided. # is this ever run ?
          if (is.null(template <- envir$precisionFactorList[[rd]]$template)) { 
            Q_CHMfactor <- Cholesky(envir$precisionFactorList[[rd]],LDL=FALSE,perm=FALSE)
            if (identical(.spaMM.data$options$TRY_update,TRUE) 
                && ! isDiagonal(envir$precisionFactorList[[rd]])) { # protection against silly bugs 
              envir$precisionFactorList[[rd]]$template <- Q_CHMfactor
            } 
          } else Q_CHMfactor <- Matrix::update(template, envir$precisionFactorList[[rd]])
          chol_Q_list[[rd]] <- as(Q_CHMfactor, "sparseMatrix")
        }
      } else if (envir$finertypes[rd]=="ranCoefs") {
        chol_Q_list[[rd]] <- envir$precisionFactorList[[rd]]$chol_Q
        precisionBlocks[[rd]] <- envir$precisionFactorList[[rd]]$precmat ## full allready including lambda 
      } else { ## both Q and factorization provided (see .assign_geoinfo_and_LMatrices_but_ranCoefs())
        chol_Q_list[[rd]] <- envir$precisionFactorList[[rd]]$chol_Q
        precisionBlocks[[rd]] <- envir$precisionFactorList[[rd]]$Qmat ## indep of w.ranef (lambda) 
      }
    } ## ELSE NULL precisionFactorList element <-> chol_Q_list element prefilled
  }
  BLOB$chol_Q <- Matrix::bdiag(chol_Q_list) 
  return(precisionBlocks)
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
      #cat(it," ")
      ZtWX[[it]] <- AUGI0_ZX$X.pv[,blocs[[it]],drop=FALSE]
      ZtWX[[it]] <- .Dvec_times_m_Matrix(attr(sXaug,"w.resid"), ZtWX[[it]])
      ZtWX[[it]] <- as.matrix(crossprod(AUGI0_ZX$ZAfix, ZtWX[[it]]))
    }
    ZtWX <- do.call(cbind, ZtWX)
  } else {
    ZtWX <- .Dvec_times_m_Matrix(attr(sXaug,"w.resid"), AUGI0_ZX$X.pv)
    ZtWX <- as.matrix(crossprod(AUGI0_ZX$ZAfix, ZtWX)) 
  }
  return(ZtWX) ## always *m*atrix
}

.get_inv_L_G_ZtsqrW <- function(BLOB, sXaug) {
  if (is.null(BLOB$inv_L_G_ZtsqrW)) {
    sqrtwZ <- .Dvec_times_Matrix(sqrt(attr(sXaug,"w.resid")), sXaug$AUGI0_ZX$ZAfix)
    BLOB$inv_L_G_ZtsqrW <- Matrix::solve(BLOB$G_CHMfactor, tcrossprod(as(BLOB$G_CHMfactor,"pMatrix"), sqrtwZ), system="L")
  }
  return(BLOB$inv_L_G_ZtsqrW)
}

.provide_BLOB_r12 <- function(BLOB, sXaug) {
  if (is.null(BLOB$r12)) {
    # Use inv_L_G_ZtsqrW when it's available (as needed for certain leverage computations) 
    # Otherwise use ZtWX which is typically much smaller so that solve(., ZtWX) is faster than solve(., ZtsqrW)
    if (is.null(BLOB$inv_L_G_ZtsqrW)) { 
      #inv_L_G_ZtsqrW <- .get_inv_L_G_ZtsqrW(BLOB, sXaug)
      if (is.null(BLOB$ZtWX)) BLOB$ZtWX <- .calc_ZtWX(sXaug) # *m*atrix
      BLOB$r12 <- Matrix::solve(BLOB$G_CHMfactor, as(BLOB$G_CHMfactor,"pMatrix") %*% BLOB$ZtWX, system="L")
    } else {
      sqrtwX <- .Dvec_times_m_Matrix(sqrt(attr(sXaug,"w.resid")),sXaug$AUGI0_ZX$X.pv)
      BLOB$r12 <- BLOB$inv_L_G_ZtsqrW %*% sqrtwX
    }
  }
}

.calc_inv_beta_cov <- function(sXaug) { ## crossprod_r22 without computing r22 (nor r12) when r22 is ot already available
  pforpv <- attr(sXaug, "pforpv") 
  if (pforpv==0L) return(Diagonal(n=0L))
  # ELSE
  ## In the working doc it is shown that the beta,beta block of inv_d2hdbv2 is solve(crossprod(r22)). We compute crossprod(r22) here,
  ## We know this crossprod is X'inv(Sig_e)X - r12'r12 in the working doc (not X'inv(Sig_e)X ! ) since r22 is defined from this. 
  AUGI0_ZX <- sXaug$AUGI0_ZX
  BLOB <- sXaug$BLOB
  if (is.null(BLOB$ZtWX)) BLOB$ZtWX <- .calc_ZtWX(sXaug)
  bloclength <- ceiling(1e7/nrow(AUGI0_ZX$X.pv)) 
  if (pforpv>bloclength) { ## I felt a need to avoid crossprod_r12 when r12 is tall (even if r12 is slim so that its crossprod is small)
    .provide_BLOB_r12(BLOB, sXaug) 
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
    .provide_BLOB_r12(BLOB, sXaug) 
    crossprod_r12 <- .crossprod(BLOB$r12)  # Xt_X !
    crossprod_r22 <- .ZtWZwrapper(AUGI0_ZX$X.pv,attr(sXaug,"w.resid")) - as.matrix(crossprod_r12) ## XtWX-Xt_X ! Both lines as explained in working doc
  }
  return(crossprod_r22) # dgeMatrix
}

.calc_r22 <- function(X.pv, w.resid, r12, sXaug) { ## given r12 is known
  ## i's vain to try to regularize crossr22 by manipulating diagonal or eigenvalues. 
  #  It would be better to control the accuracy of the computation of crossr22 as a difference.
  if (ncol(X.pv)==0L) { # sXaug may be NULL bc not accessible in all calls.
    return(diag(nrow=0L))
  } else {
    #crossr22 <- .ZtWZwrapper(X.pv,w.resid)-crossprod(r12) ## typically dsyMatrix (dense symmetric)
    crossr22 <- .ZtWZwrapper(X.pv,w.resid)-crossprod(as.matrix(r12)) ## avoids a lot of bureaucraty in '-' that is useless given what .wrap_Utri_chol() does 
    return(.wrap_Utri_chol(crossr22)) # rather than (Matrix::)chol()
  }
  ## lots of alternatives here, removed from [v2.6.54
}

.provide_assignment_qrXa_or_r22 <- function(BLOB, sXaug, AUGI0_ZX, w.ranef) {
  if (is.null(BLOB$qrXa)) { ## use it whenever available
    if ( ! identical(.spaMM.data$options$use_spprec_QR,TRUE)) {
      if (is.null(BLOB$r22)) {
        # the code assumes that when r22 is available, r12 is, so I cannot simply compute crossprod_r12 here
        .provide_BLOB_r12(BLOB, sXaug) 
        BLOB$r22 <- .calc_r22(AUGI0_ZX$X.pv,attr(sXaug,"w.resid"),BLOB$r12, sXaug)        ## both lines as explained in working doc
      }
    }
    if ( is.null(dim(BLOB$r22))) { ## not attempted to compute or failed computation
      X.pv <- as(AUGI0_ZX$X.pv,"dgCMatrix")#  .Rcpp_as_dgCMatrix(AUGI0_ZX$X.pv)
      Xscal <- with(sXaug,rbind(cbind(AUGI0_ZX$I, AUGI0_ZX$ZeroBlock), # I0_ZX order
                                cbind(AUGI0_ZX$ZAfix %*% t(solve(BLOB$chol_Q)), X.pv)))
      Xscal <- .Dvec_times_Matrix(sqrt(c(w.ranef,attr(sXaug,"w.resid"))),Xscal)
      BLOB$qrXa <- qr(Xscal)
      BLOB$X_scale <- attr(AUGI0_ZX$X.pv,"scaled:scale") ## weakness in conception: X_scaling impacts post-fit computations. 
    }
  } ## then we have either r22 or qrXa
  # no return value
}

.calc_Qt_pwy_o <- function(sXaug, pwy_o) { 
  ## This function copes with the absence of a "t_Q_scaled" extractor in spprec method in the one case we need it.
  ## Result z, viewed as a col vector,  is such that sum(z^2)=c(0,z)'c(0,z)=c(0,pwy_o)'QQ'c(0,pwy_o) 
  ##   for the orthog factor Q of the augmented matrix.
  ## The names of the variables come from the analogy with the hatval_ZX code.
  get_from_MME(sXaug, which="initialize") ## sets G_CHMfactor
  BLOB <- sXaug$BLOB
  AUGI0_ZX <- sXaug$AUGI0_ZX
  Zt_sqrtw_pwy_o <- crossprod(AUGI0_ZX$ZAfix, sqrt(attr(sXaug,"w.resid"))*pwy_o) 
  Xt_sqrtw_pwy_o <- crossprod(AUGI0_ZX$X.pv, sqrt(attr(sXaug,"w.resid"))*pwy_o)
  lev_phi_z_pwy_o <- Matrix::solve(BLOB$G_CHMfactor, as(BLOB$G_CHMfactor,"pMatrix") %*% Zt_sqrtw_pwy_o, system="L") ## R_11^{-T}.Ztw.pwy_o
  if (attr(sXaug,"pforpv")>0L) {
    .provide_assignment_qrXa_or_r22(BLOB, sXaug, AUGI0_ZX, w.ranef=attr(sXaug,"w.ranef")) 
    ## : currently provides r12, r22 ratehr than qrXaand the follwing code assumes so.
    crossprod_r12_z_pwy_o <- crossprod(BLOB$r12,lev_phi_z_pwy_o) ## R_12^T . R_11^{-T}.Ztw.pwy_o
    lev_phi_x_pwy_o <- backsolve(BLOB$r22, Xt_sqrtw_pwy_o - as.matrix(crossprod_r12_z_pwy_o), # as.matrix() has a notable effect on ohio|fitme test, and no visible drawback
                                 transpose=TRUE) ## -  R_22^{-T}.R_12^T . R_11^{-T}.Ztw.pwy_o + R_22^{-T}.Xtw.pwy_o
    return(c(drop(lev_phi_z_pwy_o),lev_phi_x_pwy_o))
  } else return(drop(lev_phi_z_pwy_o))
} # "pwt_Q_y_o"

.calc_spprec_hatval_ZX_by_QR <- function(BLOB, sXaug, AUGI0_ZX, w.ranef) { ## sparse qr is fast
  qrQ <- qr.Q(BLOB$qrXa)
  qrQ@x <- qrQ@x^2
  tmp <- rowSums(qrQ)
  n_u_h <- ncol(AUGI0_ZX$I)
  phipos <- n_u_h+seq_len(length(tmp)-n_u_h)
  BLOB$hatval_ZX <- list(lev_lambda=tmp[-phipos],lev_phi=tmp[phipos])
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
    g <- Matrix::solve(BLOB$G_CHMfactor, as(BLOB$G_CHMfactor,"pMatrix"), system="L") ## fixme rethink this.
    sqrtwZ <- .Dvec_times_Matrix(sqrt(attr(sXaug,"w.resid")), AUGI0_ZX$ZAfix)
    crossprod_r12_z_rows <- crossprod(crossprod(g, BLOB$r12),t(sqrtwZ)) 
    lev_phi_z <- colSums(g^2)
    lev_phi_z <- drop(AUGI0_ZX$ZAfix %*% lev_phi_z)*attr(sXaug,"w.resid")
  } else {
    lev_phi_z <- .get_inv_L_G_ZtsqrW(BLOB, sXaug)
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
    BLOB$hatval_ZX <- list(lev_lambda=lev_lambda,lev_phi=lev_phi)
  } else BLOB$hatval_ZX <- list(lev_lambda=lev_lambda_z,lev_phi=lev_phi_z)
}


.Sigsolve_sparsePrecision <- function(sXaug, rhs) { ## no longer used: compare .calc_inv_beta_cov() when rhs= X.pv. But useful as doc.
  v <- Matrix::solve(sXaug$BLOB$G_CHMfactor, sXaug$BLOB$ZtW %*% rhs) ## BLOB$ZtW removed !
  # implicit .Dvec_times_m_Matrix( on a dgeMatrix:
  invM.z <- attr(sXaug$AUGI0_ZX,"w.resid") * (rhs - sXaug$AUGI0_ZX$ZAfix %*% v)  # W.rhs - W Z (v= L inv(...) L' Z' W rhs) = [W- WZL inv(...) L'Z'W] rhs
  return(invM.z) ## where M=Z.solve(precMat).Zt+diag(~phi)
} # solve(AUGI0_ZX$ZAfix %*% solve(precisionMatrix) %*% t(AUGI0_ZX$ZAfix)+diag(1/attr(sXaug,"w.resid"))) = Sigsolve(diag(x=ncol(ZtW)))

.calc_H_dH <- function(BLOB, damping) {
  ## See comments in .calc_G_dG()
  if (is.null(BLOB$dH)) {
    if (is.null(BLOB$tcrossfac_Md2hdv2)) BLOB$tcrossfac_Md2hdv2 <- Matrix::solve(BLOB$chol_Q,  
                                                                                 crossprod(as(BLOB$G_CHMfactor,"pMatrix"), as(BLOB$G_CHMfactor,"sparseMatrix")))
    BLOB$Md2hdv2 <- .tcrossprod(BLOB$tcrossfac_Md2hdv2)
  }
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
    BLOB$dG <- NaN ## Todetect prolems in further usages
    G_dG <- Matrix::update(BLOB$G_CHMfactor, BLOB$Gmat, damping) # actually method of stats:::update for class CHMfactor (?`CHMfactor-class`) 
    return(list(G_dG=G_dG, dampDpD_2=damping * BLOB$D_Md2hdv2))
  } else {
    if (is.null(BLOB$dG)) {  
      if (spprec_LevM_D=="1") {
        BLOB$D_Md2hdv2 <- rep(1,ncol(BLOB$chol_Q))
        BLOB$dG <- drop0(.tcrossprod(BLOB$chol_Q, as_sym=as_sym)) # dsCMatrix if as_sym is TRUE# drop0 makes a diff in LevM.Frailty test # do we have the precisionMatrix somewhere ?
      } else { ## costly solve()
        ## part of the problem is avoiding tall matrices from tall ZA (and X), but for squarish ZA, using solve(chol_Q,... ) looks complicated.
        ## solve(chol_Q,... ) further assuming that we have not stored the LMatrix
        if (is.null(BLOB$tcrossfac_Md2hdv2)) BLOB$tcrossfac_Md2hdv2 <- Matrix::solve(BLOB$chol_Q,  
                                                                                     crossprod(as(BLOB$G_CHMfactor,"pMatrix"), as(BLOB$G_CHMfactor,"sparseMatrix")))
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

.provide_BLOB_hatval_Z_ <- function(BLOB, w.ranef, AUGI0_ZX, sXaug) {
  if (is.null(BLOB$hatval_Z_)) {
    lev_lambda <- BLOB$factor_inv_Md2hdv2
    lev_lambda@x <- lev_lambda@x^2 
    lev_lambda <- .Matrix_times_Dvec(lev_lambda, w.ranef)
    lev_lambda <- colSums(lev_lambda)
    if (AUGI0_ZX$is_unitary_ZAfix) { ## then only diagonal values of invG matter  ## adjacency-long case...?
      lev_phi <- colSums(Matrix::solve(BLOB$G_CHMfactor, as(BLOB$G_CHMfactor,"pMatrix"), system="L")^2)
      lev_phi <- drop(AUGI0_ZX$ZAfix %*% lev_phi)*attr(sXaug,"w.resid")
    } else {
      # inv_L_G_ZtsqrW arises as simplif of BLOB$factor_inv_Md2hdv2 %*% t(sqrtwZL) (L and chol_Q cancel each other)
      lev_phi <- .get_inv_L_G_ZtsqrW(BLOB, sXaug) ## solve(as(BLOB$G_CHMfactor,"sparseMatrix"), t(sqrtwZ))
      lev_phi@x <- lev_phi@x^2
      lev_phi <- colSums(lev_phi)
    }
    BLOB$hatval_Z_ <- list(lev_lambda=lev_lambda,lev_phi=lev_phi)
  }
}

.provide_BLOB_factor_inv_Md2hdv2 <- function(BLOB) { ## Code encapsulated in a function for easier profiling
  LL <- Matrix::solve(BLOB$G_CHMfactor, as(BLOB$G_CHMfactor,"pMatrix") %*% BLOB$chol_Q,system="L") ## dgCMatrix 
  # Even is BLOB$chol_Q is Identy, this this hardly slower than Matrix::solve(BLOB$G_CHMfactor, system="L") so cannot be improved.
  BLOB$factor_inv_Md2hdv2 <- drop0(LL,tol=.Machine$double.eps) ## crossfactor
  #factor_Md2hdv2 <- as.matrix(Matrix::solve(BLOB$chol_Q,as(BLOB$G_CHMfactor,"sparseMatrix")))
  #BLOB$chol_Md2hdv2 <- Rcpp_chol_R(tcrossprod(factor_Md2hdv2)) #A = R'.R as in R's chol()
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
        # Either a precisionBlocks is that of a non-trivial corr mat of a gaussian ranef => the lambda and w.ranef are constant (even in CAR)
        # or the lambda and w.ranef are not constant but the precisionBlocks is an identity matrix.
        # In both cases, precisionBlocks[[it]] * w.ranef[u.range] is a correct way of computing the precision block
        # .Matrix_times_Dvec() is faster by directly manipulating slots.
        u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
        precisionBlocks[[it]] <- .Matrix_times_Dvec(precisionBlocks[[it]], 
                                                    w.ranef[u.range],
                                                    check_dsC=FALSE, # risky: use TRUE to detect when the result cannot be dsC
                                                    keep_dsC=TRUE # result MUST be dsC
                                                    ) 
      }
    }
    if (length(precisionBlocks)>1L) {
      precisionMatrix <- forceSymmetric(do.call(Matrix::bdiag, precisionBlocks)) # bdiag degrades dsC into dgC
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
      warning("Possibly inefficient code in .AUGI0_ZX_sparsePrecision(). Please report to the package maintainer.", immediate.=TRUE) 
      tmp <- tmp + precisionMatrix
    }
    BLOB$Gmat <- drop0(tmp) ## depends on w.ranef and w.resid
    # if (inherits(BLOB$Gmat,"dgCMatrix")) stop("ICI") 
    if (is.null(template <- AUGI0_ZX$template_G_CHM)) { 
      BLOB$G_CHMfactor <- Cholesky(BLOB$Gmat,LDL=FALSE,perm=.spaMM.data$options$perm_G) ## costly
      if (identical(.spaMM.data$options$TRY_update,TRUE) 
          && ! isDiagonal(precisionMatrix)) { # protection against silly bugs # F I X M E Could we identify cases where prec matrix is always Diagonal ?
        AUGI0_ZX$template_G_CHM <- BLOB$G_CHMfactor
      }
    } else BLOB$G_CHMfactor <- Matrix::update(template, BLOB$Gmat) ## If it fails in spde, look for too low kappa.
    ## with perm=TRUE G=P'LL'P and the P'L (non-triangular) factor is given by solve(<G_CHM>,as(<G_CHM>,"sparseMatrix"),system="Pt")
  }
  if (which=="initialize") return(NULL)
  if ( is.null(BLOB$factor_inv_Md2hdv2)  
       && which %in% c("","solve_d2hdv2","Mg_invH_g","hatval_Z","hatval", "Mg_invH_g")) .provide_BLOB_factor_inv_Md2hdv2(BLOB)
  if (which=="hatval_Z") {
    ## calcul correct; 
    # TT <- rbind(diag(sqrt(w.ranef)),diag(sqrt(attr(AUGI0_ZX, "w.resid"))) %*% AUGI0_ZX$ZAfix %*% t(solve(BLOB$chol_Q)) )
    # diag(TT %*% (solve(crossprod(TT))) %*% t(TT)) ## lev_lambda, lev_phi
    .provide_BLOB_hatval_Z_(BLOB, w.ranef, AUGI0_ZX, sXaug)
    return(BLOB$hatval_Z_)
  }
  if (which=="Mg_solve_g") {
    stop("Mg_solve_g not implemented, using alternative way")
  } 
  if (which=="Mg_invH_g") { ## - grad_v invD2hdv2 grad_v 
    LLgrad <- BLOB$factor_inv_Md2hdv2 %*% B
    return( sum(LLgrad^2))
  }
  if (which=="Mg_invXtWX_g") { ## 
    if (is.null(BLOB$XtWX)) BLOB$XtWX <- .ZtWZwrapper(AUGI0_ZX$X.pv,attr(sXaug,"w.resid")) ## as(,"matrix") ?
    Mg_invXtWX_g <- crossprod(B,solve(BLOB$XtWX,B))[1L] # [1L] drops possible Matrix class...
    return(Mg_invXtWX_g)
  }
  if (which=="solve_d2hdv2") { 
      if ( is.null(B) ) {
        # stop("B is NULL") # this would actually 'work' bc this case never occurs in current code (v2.6.3)
        return( - .crossprod(BLOB$factor_inv_Md2hdv2))
      } else {return( - .crossprod(BLOB$factor_inv_Md2hdv2, BLOB$factor_inv_Md2hdv2 %*% B))}
  }
  if (which=="half_logdetQ") {
    half_logdetQ <-  sum(log(diag(BLOB$chol_Q))) 
    return(half_logdetQ)
  }
  if (which=="half_logdetG") {
    half_logdetG <- Matrix::determinant(BLOB$G_CHMfactor)$modulus[1]
    return(half_logdetG)
  }
  # cases where Sigsolve is needed ## Sig is the Cov(response), tcrossprod(ZAL)/lambda+diag(1/w.resid) 
  #X.pv <- AUGI0_ZX$X.pv
  ## Solving for model coefficients:
  if ( ! is.null(z)) { ## which=""
    grad_v <- z$m_grad_obj[seq(ncol(BLOB$chol_Q))]
    if (attr(sXaug,"pforpv")) { ## if there are fixed effect coefficients to estimate
      # rhs for beta
      grad_p_v <- z$m_grad_obj[-seq(ncol(BLOB$chol_Q))]
      if (is.null(BLOB$r22)) {
        .provide_BLOB_r12(BLOB, sXaug) 
        BLOB$r22 <- .calc_r22(AUGI0_ZX$X.pv,attr(sXaug,"w.resid"), BLOB$r12, sXaug)        ## both lines as explained in working doc
      }
      rhs <- (BLOB$factor_inv_Md2hdv2 %*% grad_v)[,1L]
      dbeta_rhs <- grad_p_v - .crossprod(BLOB$r12, rhs)[,1L] # numeric  vector
      dbeta_eta <- backsolve(BLOB$r22, backsolve(BLOB$r22, dbeta_rhs, transpose = TRUE)) # ie (chol2inv(BLOB$r22) %*% dbeta_rhs)[,1L]
      dv_h <- .crossprod(BLOB$factor_inv_Md2hdv2, rhs - (BLOB$r12 %*% dbeta_eta)[,1])[,1]
    } else {
      dbeta_eta <- numeric(0)
      dv_h <- .crossprod(BLOB$factor_inv_Md2hdv2, (BLOB$factor_inv_Md2hdv2 %*% grad_v)[,1L])[,1L]
    }
    return(list(dv_h=dv_h,dbeta_eta=dbeta_eta)) 
  }
  if (which=="beta_cov") { 
    if (attr(sXaug,"pforpv")==0L) return(list(beta_cov=diag(nrow=0L)))
    # ELSE
    if (is.null(BLOB$r22)) {
      inv_beta_cov <- .calc_inv_beta_cov(sXaug) 
    } else inv_beta_cov <- crossprod(BLOB$r22)
    beta_cov <- try(solve(inv_beta_cov), silent=TRUE) ## chol2inv(Rscaled) in other methods.
    if (inherits(beta_cov,"try-error")) {
      warning("Information matrix is (nearly) singular.")
      beta_cov <- ginv(as.matrix(inv_beta_cov))
    }
    colnames(beta_cov) <- rownames(beta_cov) <- colnames(AUGI0_ZX$X.pv)
    return(list(beta_cov=as.matrix(beta_cov))) # make_beta_table() expects a *m*atrix
  }
  #### REML:
  if (which=="logdet_r22") { 
    if (is.null(BLOB$r22)) {
      if (attr(sXaug,"pforpv")==0L) {
        logdet_r22 <- 0 
      } else if (inherits(AUGI0_ZX$X.pv,"sparseMatrix")) { ## sparse qr is fast
        #cat("qr:")
        Xscal <- with(sXaug,rbind(cbind(AUGI0_ZX$I,AUGI0_ZX$ZeroBlock),cbind(AUGI0_ZX$ZAfix %*% t(solve(BLOB$chol_Q)),AUGI0_ZX$X.pv)))
        Xscal <- .Dvec_times_Matrix(sqrt(c(attr(sXaug,"w.ranef"),attr(sXaug,"w.resid"))),Xscal)
        qrX <- qr(Xscal)
        logdet_r22 <- determinant(qrR(qrX,backPermute = FALSE))$modulus[1] - get_from_MME(sXaug,"logdet_sqrt_d2hdv2")
      } else {
        #cat("cross_r22:") # FIXME should use X.Re for non-standard REML? (don't modify thia but code which may call for it)
        # _FIXME_ is it faster to compute crossprod_r22 or r22 by .calc_r22() ?
        crossprod_r22 <- .calc_inv_beta_cov(sXaug) # inv_beta_cov is crossprod_r22 as shown in working doc.
        logdet_r22 <- determinant(crossprod_r22)$modulus[1]/2
      }
    } else {
      #cat("r22:")
      logdet_r22 <- determinant(BLOB$r22)$modulus[1]
    }
    return(logdet_r22)
  }
  if (which=="hatval") { ##  REML hatval computation 
    if (is.null(BLOB$hatval_ZX)) {
      .provide_assignment_qrXa_or_r22(BLOB=BLOB, AUGI0_ZX=AUGI0_ZX, sXaug=sXaug, w.ranef=w.ranef)
      if ( ! is.null(BLOB$qrXa)) {
        .calc_spprec_hatval_ZX_by_QR(BLOB=BLOB, AUGI0_ZX=AUGI0_ZX, sXaug=sXaug, w.ranef=w.ranef)
      } else .calc_spprec_hatval_ZX_by_r22(BLOB=BLOB, AUGI0_ZX=AUGI0_ZX, sXaug=sXaug, w.ranef=w.ranef) ## fragile bc  BLOB$r22 <- .calc_r22() is fragile...
    }
    return(BLOB$hatval_ZX)
  }
  if (which %in% c("LevMar_step")) {  ## LM with diagonal perturbation as in NocedalW p. 266
    #
    if (.spaMM.data$options$use_G_dG) { ## not so clear which is faster ## big-ranefs would be a good test F I X M E
      # uses rather complex code to avoid solve(dense d2hdv2,...) but cannot avoid solve(rather dense G_dG)
      if (is.null(BLOB$ZtWX)) BLOB$ZtWX <- .calc_ZtWX(sXaug)
      if (is.null(BLOB$XtWX)) BLOB$XtWX <- .ZtWZwrapper(AUGI0_ZX$X.pv,attr(sXaug,"w.resid")) ## as(,"matrix") ?
      G_dG_blob <- .calc_G_dG(BLOB, damping) # list(G_dG=G_dG, dampDpD_2=damping * BLOB$D_Md2hdv2)
      if (is.null(BLOB$DpD)) BLOB$DpD <- c(BLOB$D_Md2hdv2,diag(BLOB$XtWX)) ## store as diagonal perturbation for the return value
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
      if (is.null(BLOB$DpD)) {  
        if (is.null(BLOB$LZtWX)) {
          if (is.null(BLOB$ZtWX)) BLOB$ZtWX <- .calc_ZtWX(sXaug)
          BLOB$LZtWX <- as(solve(BLOB$chol_Q, BLOB$ZtWX),"matrix")
        }
        if (is.null(BLOB$XtWX)) BLOB$XtWX <- .ZtWZwrapper(AUGI0_ZX$X.pv,attr(sXaug,"w.resid")) 
        BLOB$DpD <- c(BLOB$D_Md2hdv2,diag(BLOB$XtWX)) ## store as diagonal perturbation for the retrun value
      } 
      ## sequel recomputed for each new damping value...
      grad_v <- LM_z$scaled_grad[seq(ncol(BLOB$chol_Q))]
      grad_beta <- LM_z$scaled_grad[-seq(ncol(BLOB$chol_Q))]
      if ((useQR <- FALSE)) {
        if (is.null(BLOB$D_Md2hdv2)) {
          ## t() bc .damping_to_solve implicitly uses crossprod ( as in solve(X'X)=solve(R'R) ) (or equivalently tcrossprod(solve))
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
  if (which=="Qty") { ## 
    if (is.null(BLOB$hatval_ZX)) {
      .provide_assignment_qrXa_or_r22(BLOB=BLOB, AUGI0_ZX=AUGI0_ZX, sXaug=sXaug, w.ranef=w.ranef)
      if ( ! is.null(BLOB$qrXa)) {
        .calc_spprec_hatval_ZX_by_QR(BLOB=BLOB, AUGI0_ZX=AUGI0_ZX, sXaug=sXaug, w.ranef=w.ranef)
      } else .calc_spprec_hatval_ZX_by_r22(BLOB=BLOB, AUGI0_ZX=AUGI0_ZX, sXaug=sXaug, w.ranef=w.ranef) ## fragile bc  BLOB$r22 <- .calc_r22() is fragile...
    }
    return(BLOB$hatval_ZX)
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
  } else resu <- switch(which,
                        "logdet_sqrt_d2hdv2" = { ## sqrt( logdet info(v) = -d2hdv2 ) 
                          half_logdetQ <- .AUGI0_ZX_sparsePrecision(sXaug=sXaug,which="half_logdetQ")
                          (- half_logdetQ + .AUGI0_ZX_sparsePrecision(sXaug=sXaug,which="half_logdetG"))
                        }, ## all other cases:
                        .AUGI0_ZX_sparsePrecision(sXaug=sXaug,which=which,z=szAug,B=B,damping=damping,LM_z=LM_z)
                       )
  return(resu)
  
}
