def_AUGI0_ZX_sparsePrecision <- function(AUGI0_ZX, corrPars,w.ranef,cum_n_u_h,w.resid, 
                                         #weight_X, ## currently ignored
                                         H_global_scale ## currently ignored
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
  for (it in seq_len(length(envir$precisionFactorList))) {
    if( ! is.null(envir$precisionFactorList[[it]])) {
      if (inherits(envir$precisionFactorList[[it]],c("Matrix","matrix"))) { ## Q_CHMfactor NOT provided
        ## not clear that adjacency or AR1 can occur here
        if (envir$finertypes[it] %in% c("AR1", "ranCoefs", "adjacency")) { # adjacency-specific code removed [from version 2.5.30
          ## envir$precisionFactorList[[it]] should be a list hence this block should not be reached
          stop("reached here in .reformat_Qmat_info()") 
        } else {## other models where a single matrix is provided : Qmatrix provided.
          chol_Q_list[[it]] <- as(Cholesky(envir$precisionFactorList[[it]],LDL=FALSE,perm=FALSE), "sparseMatrix")
        }
      } else if (envir$finertypes[it]=="ranCoefs") {
        chol_Q_list[[it]] <- envir$precisionFactorList[[it]]$chol_Q
        precisionBlocks[[it]] <- envir$precisionFactorList[[it]]$precmat ## full allready including lambda 
      } else { ## both Q and factorization provided (see .assign_geoinfo_and_LMatrices_but_ranCoefs())
        chol_Q_list[[it]] <- envir$precisionFactorList[[it]]$chol_Q
        precisionBlocks[[it]] <- envir$precisionFactorList[[it]]$Qmat ## indep of w.ranef (lambda) 
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

.calc_inv_beta_cov <- function(sXaug) { ## crossprod_r22
  ## In the working doc it is shown that the beta,beta block of inv_d2hdbv2 is solve(crossprod(r22)). We compute crossprod(r22) here,
  ## We know this crossprod is X'inv(Sig_e)X - r12'r12 in the working doc (not X'inv(Sig_e)X ! ) since r22 is defined from this. 
  AUGI0_ZX <- sXaug$AUGI0_ZX
  BLOB <- sXaug$BLOB
  if (is.null(BLOB$ZtWX)) BLOB$ZtWX <- .calc_ZtWX(sXaug)
  bloclength <- ceiling(1e7/nrow(AUGI0_ZX$X.pv)) 
  pforpv <- attr(sXaug, "pforpv") 
  if (pforpv==0L) {
    return(Diagonal(n=0L))
  } else if (pforpv>bloclength) {
    #cat("_<_")
    seqncol <- seq(pforpv)
    blocs <- split(seqncol, ceiling(seq_along(seqncol)/bloclength))
    cross_r22 <- vector("list",length(blocs))
    for (it in seq_along(blocs)) {
      #cat(it," ")
      cross_r22[[it]] <- BLOB$ZtWX[,blocs[[it]],drop=FALSE] ##*m*atrix
      cross_r22[[it]] <- solve(BLOB$G_CHMfactor, cross_r22[[it]],system="A") ## ~v ## dgeMatrix <- *m*atrix
      if(! is.null(BLOB$G_scaling)) cross_r22[[it]] <- cross_r22[[it]]/BLOB$G_scaling
      cross_r22[[it]] <- AUGI0_ZX$X.pv[,blocs[[it]],drop=FALSE] - AUGI0_ZX$ZAfix %*% cross_r22[[it]]
      cross_r22[[it]] <- .Dvec_times_Matrix( attr(sXaug,"w.resid"),cross_r22[[it]]) ## ~invM.z
      cross_r22[[it]] <- crossprod(AUGI0_ZX$X.pv, cross_r22[[it]])
    }
    cross_r22 <- do.call(cbind, cross_r22)
    #cat("_>_")
  } else {
    cross_r22 <- solve(BLOB$G_CHMfactor, BLOB$ZtWX,system="A") ## dgeMatrix <- *m*atrix
    if(! is.null(BLOB$G_scaling)) cross_r22 <- cross_r22/BLOB$G_scaling
    cross_r22 <- AUGI0_ZX$X.pv - AUGI0_ZX$ZAfix %*% cross_r22 ## still dgeMatrix
    cross_r22 <- .Dvec_times_Matrix( attr(sXaug,"w.resid"),cross_r22) ## invM.X
    cross_r22 <- crossprod(AUGI0_ZX$X.pv, cross_r22) ## Xt.invM.X
    # # That is equivalent to
    # r12 <- solve(BLOB$G_CHMfactor, BLOB$ZtWX,system="L") 
    # h22 <- .ZtWZwrapper(AUGI0_ZX$X.pv,attr(sXaug,"w.resid")) - crossprod(r12) ## both lines as explained in working doc
    # FIXME keep the most efficient algo ! But not obvious difference in efficiency
  }
  return(cross_r22) # dgeMatrix
}

.calc_r22 <- function(X.pv, w.resid, r12, sXaug) {
  ## i's vain to try to regularize crossr22 by manipulating diagonal or eigenvalues. 
  #  It would be better to control the accuracy of the computation of crossr22 as a difference.
  if (ncol(X.pv)==0L) { # sXaug may be NULL bc not accessible in all calls.
    return(diag(nrow=0L))
  } else {
    crossr22 <- .ZtWZwrapper(X.pv,w.resid)-crossprod(r12) ## typically dsyMatrix (dense symmetric)
    return(.wrap_Utri_chol(crossr22)) # rather than (Matrix::)chol()
  }
  
  
  ## lots of alternatives here:
  
  
  if (FALSE) {
    scaling <- 1/sqrt(colSums(r12^2)) ## i.e. control diag(crossprod(r12))
    locX <- .m_Matrix_times_Dvec(X.pv, scaling)
    locr12 <- .m_Matrix_times_Dvec(r12, scaling)
    crossscaledr22 <- .ZtWZwrapper(locX,w.resid)-crossprod(locr12)
    cholscaled <- chol(crossscaledr22)
    return(.m_Matrix_times_Dvec(cholscaled, 1/scaling))
  } else if (TRUE) {
    crossr22 <- .ZtWZwrapper(X.pv,w.resid)-crossprod(r12)
    resu <- try(chol(crossr22),silent=TRUE)
    if (inherits(resu,"try-error")) {
      BLOB <- sXaug$BLOB
      # AUGI0_ZX <- sXaug$AUGI0_ZX
      scaling <- 1/sqrt(colSums(r12^2)) ## i.e. control diag(crossprod(r12))
      locX <- .m_Matrix_times_Dvec(X.pv, scaling)
      # locZtWX <-  .Dvec_times_m_Matrix(attr(sXaug,"w.resid"),locX)
      # locZtWX <- as.matrix(crossprod(AUGI0_ZX$ZAfix, locZtWX)) 
      # cross_r22 <- solve(BLOB$G_CHMfactor, locZtWX) ## dgeMatrix <- *m*atrix
      # cross_r22 <- locX - AUGI0_ZX$ZAfix %*% cross_r22 ## still dgeMatrix
      # cross_r22 <- .Dvec_times_Matrix( attr(sXaug,"w.resid"),cross_r22) ## invM.X
      # cross_r22 <- crossprod(locX, cross_r22) ## Xt.invM.X
      cross_r22 <- t(locX) %*% BLOB$precisionMatrix %*% locX
      cholscaled <- chol(cross_r22)
      return(.m_Matrix_times_Dvec(cholscaled, 1/scaling))
    } else return(resu)
  } else {
    crossr22 <- .ZtWZwrapper(X.pv,w.resid)-crossprod(r12)
    resu <- try(chol(crossr22),silent=TRUE)
    if (inherits(resu,"try-error")) {
      crossr22 <- .calc_inv_beta_cov(sXaug)
    }
    return(chol(crossr22))
  }
}

.provide_assignment_qrXa_or_r22 <- function(BLOB, sXaug, AUGI0_ZX, w.ranef) {
  if (is.null(BLOB$qrXa)) { ## use it whenever available
    if ( ! identical(.spaMM.data$options$use_spprec_QR,TRUE)) {
      if (is.null(BLOB$r22)) {
        if (is.null(BLOB$ZtWX)) BLOB$ZtWX <- .calc_ZtWX(sXaug)
        BLOB$r12 <- solve(BLOB$G_CHMfactor, BLOB$ZtWX,system="L") 
        BLOB$r22 <- try(.calc_r22(AUGI0_ZX$X.pv,attr(sXaug,"w.resid"),BLOB$r12, sXaug),silent=TRUE)         ## both lines as explained in working doc
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
  lev_phi_z_pwy_o <- Matrix::solve(BLOB$G_CHMfactor, Zt_sqrtw_pwy_o, system="L") ## R_11^{-T}.Ztw.pwy_o
  if (attr(sXaug,"pforpv")>0L) {
    .provide_assignment_qrXa_or_r22(BLOB, sXaug, AUGI0_ZX, w.ranef=attr(sXaug,"w.ranef")) 
    ## : currently provides r12, r22 ratehr than qrXaand the follwing code assumes so.
    crossprod_r12_z_pwy_o <- crossprod(BLOB$r12,lev_phi_z_pwy_o) ## R_12^T . R_11^{-T}.Ztw.pwy_o
    lev_phi_x_pwy_o <- backsolve(BLOB$r22, Xt_sqrtw_pwy_o - crossprod_r12_z_pwy_o, 
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
  #sqrtwZ <- AUGI0_ZX$ZAfix
  #sqrtwZ@x <- sqrtwZ@x * sqrt(attr(sXaug,"w.resid"))[sqrtwZ@i+1L]
  sqrtwZ <- .Dvec_times_Matrix(sqrt(attr(sXaug,"w.resid")), AUGI0_ZX$ZAfix)
  XtW <- t(AUGI0_ZX$X.pv *sqrt(attr(sXaug,"w.resid")))
  if (AUGI0_ZX$is_unitary_ZAfix) {
    # no complete simplif for AUGI0_ZX$is_unitary_ZAfix bc lev_phi_x is not the diag of a Z...Zt matrix
    # but we can still use p*n matrices instead of r*n matrices to compute the lev_phi_x component
    g <- Matrix::solve(BLOB$G_CHMfactor, system="L") ## fixme rethink this.
    crossprod_r12_z_rows <- suppressMessages(crossprod(crossprod(g, BLOB$r12),t(sqrtwZ))) 
    lev_phi_z <- colSums(g^2)
    lev_phi_z <- drop(AUGI0_ZX$ZAfix %*% lev_phi_z)*attr(sXaug,"w.resid")
  } else {
    lev_phi_z <- Matrix::solve(BLOB$G_CHMfactor, t(sqrtwZ), system="L") ## solve(as(BLOB$G_CHMfactor,"sparseMatrix"), t(sqrtwZ))
    crossprod_r12_z_rows <- crossprod(BLOB$r12,lev_phi_z)
    lev_phi_z@x <- lev_phi_z@x^2
    lev_phi_z <- colSums(lev_phi_z) 
  }
  if (attr(sXaug,"pforpv")>0L) {
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

# (FIXME) Matrix::solve not used economically in this code => long time of test-Rasch when sparse precision is used.
# this is why spprec is useful only when prec mat is quite sparser than corr mat. 
# It involves the following calls ( * => mMatrix rhs ) :
# * solve(BLOB$G_CHMfactor, chol_Q, "L")
# * solve(BLOB$G_CHMfactor, t(sqrtwZ), "L")
# solve(BLOB$G_CHMfactor, rhs, "A")
# solve(BLOB$G_CHMfactor, (chol_Q %*% grad_v), "A")
# * solve(BLOB$G_CHMfactor, ZtWX, "A")
# solve(BLOB$G_CHMfactor, ZtWX %*% dbeta_eta, "A")
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
        u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
        precisionBlocks[[it]] <- .Matrix_times_Dvec(precisionBlocks[[it]], w.ranef[u.range]) ## more vanilla code actually degrades a dsC into dgC 
      }
    }
    if (length(precisionBlocks)>1L) {
      precisionMatrix <- forceSymmetric(do.call(Matrix::bdiag, precisionBlocks)) # bdiag degrades dsC into dgC
    } else precisionMatrix <- (precisionBlocks[[1L]]) # forceSymmetric removed bc input should in principle be dsC
    tmp <- .ZtWZwrapper(AUGI0_ZX$ZAfix,attr(sXaug,"w.resid")) # should return a symmetric-type matrix (dsC or ddi...)
    tmp <- tmp + precisionMatrix 
    BLOB$Gmat <- drop0(tmp) ## depends on w.ranef and w.resid
    # if (inherits(BLOB$Gmat,"dgCMatrix")) stop("ICI") 
    BLOB$G_CHMfactor <- Cholesky(BLOB$Gmat,LDL=FALSE,perm=FALSE) ## costly
    ## with perm=TRUE G=P'LL'P and the P'L (non-triangular) factor is given by solve(<G_CHM>,as(<G_CHM>,"sparseMatrix"),system="Pt")
  }
  if (which=="initialize") return(NULL)
  if ( is.null(BLOB$factor_inv_Md2hdv2)  && 
       #(which %in% c("logdet_R_scaled_v","hatval_Z","d2hdv2","solve_d2hdv2","R_scaled_v_h_blob"))
       (which %in% c("solve_d2hdv2","Mg_invH_g","hatval_Z","hatval"))
  ) {
    LL <- Matrix::solve(BLOB$G_CHMfactor, BLOB$chol_Q,system="L") ## dgCMatrix 
    BLOB$factor_inv_Md2hdv2 <- drop0(LL,tol=.Machine$double.eps)
    #factor_Md2hdv2 <- as.matrix(Matrix::solve(BLOB$chol_Q,as(BLOB$G_CHMfactor,"sparseMatrix")))
    #BLOB$chol_Md2hdv2 <- Rcpp_chol_R(tcrossprod(factor_Md2hdv2)) #A = R'.R as in R's chol()
    ####BLOB$inv_Md2hdv2 <- Matrix::crossprod(BLOB$factor_inv_Md2hdv2)
    ####BLOB$CHMfactor_inv_Md2hdv2 <- Matrix::Cholesky(BLOB$inv_Md2hdv2,LDL=FALSE,perm=FALSE)
  }
  
  if (which=="hatval_Z") {
    ## calcul correct; 
    # TT <- rbind(diag(sqrt(w.ranef)),diag(sqrt(attr(AUGI0_ZX, "w.resid"))) %*% AUGI0_ZX$ZAfix %*% t(solve(BLOB$chol_Q)) )
    # diag(TT %*% (solve(crossprod(TT))) %*% t(TT)) ## lev_lambda, lev_phi
    if (is.null(BLOB$hatval_Z_)) {
      lev_lambda <- BLOB$factor_inv_Md2hdv2
      lev_lambda@x <- lev_lambda@x^2 
      lev_lambda <- .Matrix_times_Dvec(lev_lambda, w.ranef)
      lev_lambda <- colSums(lev_lambda)
      if (AUGI0_ZX$is_unitary_ZAfix) { ## then only diagonal values of invG matter
        lev_phi <- colSums(Matrix::solve(BLOB$G_CHMfactor, system="L")^2)
        lev_phi <- drop(AUGI0_ZX$ZAfix %*% lev_phi)*attr(sXaug,"w.resid")
      } else {
        #sqrtwZ <- AUGI0_ZX$ZAfix
        #sqrtwZ@x <- sqrtwZ@x * sqrt(attr(sXaug,"w.resid"))[sqrtwZ@i+1L]
        sqrtwZ <- .Dvec_times_Matrix(sqrt(attr(sXaug,"w.resid")), AUGI0_ZX$ZAfix)
        # next line is simplif of BLOB$factor_inv_Md2hdv2 %*% t(sqrtwZL) (L and chol_Q cancel each other)
        lev_phi <- Matrix::solve(BLOB$G_CHMfactor, t(sqrtwZ), system="L") ## solve(as(BLOB$G_CHMfactor,"sparseMatrix"), t(sqrtwZ))
        lev_phi@x <- lev_phi@x^2
        lev_phi <- colSums(lev_phi)
      }
      BLOB$hatval_Z_ <- list(lev_lambda=lev_lambda,lev_phi=lev_phi)
    }
    return(BLOB$hatval_Z_)
  }
  if (which=="Mg_solve_g") {
    stop("Mg_solve_g not implemented, using alternative way")
  } else if (which=="Mg_invH_g") { ## - grad_v invD2hdv2 grad_v 
    LLgrad <- BLOB$factor_inv_Md2hdv2 %*% B
    return( sum(LLgrad^2))
  }
  if (which=="solve_d2hdv2") { 
      if ( is.null(B) ) {
        stop("B is NULL")
      } else {return( - Matrix::crossprod(BLOB$factor_inv_Md2hdv2, BLOB$factor_inv_Md2hdv2 %*% B))}
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
    #L_dv_term_from_grad_v <- Matrix::solve(BLOB$G_CHMfactor, BLOB$factor_inv_Md2hdv2 %*% grad_v, system="Lt")
    L_dv_term_from_grad_v <- Matrix::solve(BLOB$G_CHMfactor, BLOB$chol_Q %*% grad_v,system="A") # part from grad_h, to be completed
    if (attr(sXaug,"pforpv")) { ## if there are fixed effect coefficients to estimate
      if (is.null(BLOB$ZtWX)) BLOB$ZtWX <- .calc_ZtWX(sXaug)
      XtWX <- .ZtWZwrapper(AUGI0_ZX$X.pv,attr(sXaug,"w.resid")) 
      # rhs for beta
      grad_p_v <- z$m_grad_obj[-seq(ncol(BLOB$chol_Q))]
      dbeta_rhs <- Matrix::drop(grad_p_v - crossprod(BLOB$ZtWX, L_dv_term_from_grad_v)) 
      # matrix for beta
      Xt_X <- crossprod(BLOB$ZtWX, solve(BLOB$G_CHMfactor, BLOB$ZtWX, system="A"))
      dbeta_eta <- solve(XtWX-Xt_X , dbeta_rhs)[,1] ## wher Xt(W-_)X = X'SigSolve X
      # modif rhs for v_h :
      L_dv <- L_dv_term_from_grad_v - solve(BLOB$G_CHMfactor, BLOB$ZtWX %*% dbeta_eta, system="A")
      dv_h <- Matrix::drop(Matrix::crossprod(BLOB$chol_Q, Matrix::drop(L_dv))) ## inner drop() to avoid signature issue with dgeMatrix dv_rhs...
    } else {
      dbeta_eta <- numeric(0)
      dv_h <- Matrix::drop(Matrix::crossprod(BLOB$chol_Q, Matrix::drop(L_dv_term_from_grad_v))) ## inner drop() to avoid signature issue with dgeMatrix dv_rhs...
    }
    return(list(dv_h=dv_h,dbeta_eta=dbeta_eta)) 
  }
  if (which=="beta_cov") { 
    if (attr(sXaug,"pforpv")==0L) return(list(beta_cov=diag(nrow=0L)))
    # ELSE
    inv_beta_cov <- .calc_inv_beta_cov(sXaug) 
    beta_cov <- try(solve(inv_beta_cov), silent=TRUE) ## chol2inv(Rscaled) in other methods.
    if (inherits(beta_cov,"try-error")) {
      warning("Information matrix is (nearly) singular.")
      beta_cov <- ginv(as.matrix(inv_beta_cov))
    }
    colnames(beta_cov) <- rownames(beta_cov) <- colnames(AUGI0_ZX$X.pv)
    return(list(beta_cov=as.matrix(beta_cov)))
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
        #cat("cross_r22:") # FIXME should use X.Re for non-standard REML?
        crossprod_r22 <- .calc_inv_beta_cov(sXaug) # inv_beta_cov is crossprod(diagonal block r22 of cholesky factor) as shown in working doc.
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
    if (TRUE) { ## not so clear which is faster
      # uses rather complex code to avoid solve(dense d2hdv2,...) but cannot avoid solve(rather dense G_dG)
      if (is.null(BLOB$dG)) {  
        tmp <- Matrix::solve(BLOB$chol_Q, as(BLOB$G_CHMfactor,"sparseMatrix"))
        tmp@x <- tmp@x^2
        BLOB$D_Md2hdv2 <- rowSums(tmp) ## rowSums = diag Tcrossprod = invQ G Gt invQt
        ## convert diag perturb of Md2hdv2 into a non-diag perturb of G :
        ## since H=invL_Q G t(invL_Q),  dG= L_Q dH t(L_Q) 
        BLOB$dG <- .ZWZtwrapper(BLOB$chol_Q , BLOB$D_Md2hdv2)
        if (is.null(BLOB$ZtWX)) BLOB$ZtWX <- .calc_ZtWX(sXaug)
        BLOB$XtWX <- .ZtWZwrapper(AUGI0_ZX$X.pv,attr(sXaug,"w.resid")) ## as(,"matrix") ?
        BLOB$DpD <- c(BLOB$D_Md2hdv2,diag(BLOB$XtWX)) ## store as diagonal perturbation for the return value
      } 
      ## sequel recomputed for each new damping value...
      # if (is.null(BLOB$Gmat_dense)) BLOB$Gmat_dense <- as.matrix(BLOB$Gmat)
      # G_dG <- BLOB$Gmat_dense + damping*as.matrix(BLOB$dG) ## probably not so sparse...
      # nevertheless the sparse code is faster
      grad_v <- LM_z$scaled_grad[seq(ncol(BLOB$chol_Q))] 
      dampDpD_2 <- (damping*BLOB$dG) ## not diagonal...
      G_dG <- drop0(BLOB$Gmat + dampDpD_2) ## probably not so sparse...
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
        #
        dbeta_eta <- solve(XX_D-Xt_X , dbeta_rhs)[,1]
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
        BLOB$XtWX <- .ZtWZwrapper(AUGI0_ZX$X.pv,attr(sXaug,"w.resid")) 
        BLOB$DpD <- c(BLOB$D_Md2hdv2,diag(BLOB$XtWX)) ## store as diagonal perturbation for the retrun value
      } 
      ## sequel recomputed for each new damping value...
      grad_v <- LM_z$scaled_grad[seq(ncol(BLOB$chol_Q))]
      grad_beta <- LM_z$scaled_grad[-seq(ncol(BLOB$chol_Q))]
      if ((useQR <- FALSE)) {
        if (is.null(BLOB$D_Md2hdv2)) {
          ## t() bc .damping_to_solve implictly uses crossprod ( as in solve(X'X)=solve(R'R) )
          BLOB$Hfac <- tmp <- t(drop0(solve(BLOB$chol_Q, as(BLOB$G_CHMfactor,"sparseMatrix")))) ## probably not so sparse...
          tmp@x <- tmp@x^2
          BLOB$D_Md2hdv2 <- colSums(tmp) ## colSums bc t() above // rowSums for diag Tcrossprod(invQ G)  = invQ G Gt invQt
        }
        dampDpD_2 <- damping * BLOB$D_Md2hdv2
        DS <- .damping_to_solve(X=BLOB$Hfac, dampDpD=dampDpD_2, rhs=NULL) 
        dv_term_from_grad_v <- drop((DS$inv %*% grad_v[DS$Rperm])[DS$RRsP]) 
      } else {
        if (is.null(BLOB$Md2hdv2)) BLOB$Md2hdv2 <- as(tcrossprod(Matrix::solve(BLOB$chol_Q, as(BLOB$G_CHMfactor,"sparseMatrix"))),"matrix")
        H_dH <- BLOB$Md2hdv2 ## fixme inefficace en soi d'autant plus qu'on a un rhs
        diag(H_dH) <- diag(H_dH)*(1+damping)
        dv_term_from_grad_v <- solve(H_dH, grad_v)
      }
      if (attr(sXaug,"pforpv")) { ## if there are fixed effect coefficients to estimate
        dbeta_rhs <- grad_beta - crossprod(BLOB$LZtWX, dv_term_from_grad_v) # - dampDpD_1 * LM_z$v_h_beta$beta_eta one damps the gradient eq, not the one giving v_h_beta
        if (useQR) {
          thinsolve <- as((DS$inv %*% BLOB$LZtWX[DS$Rperm,])[DS$RRsP,],"matrix") ## something wrong here ?
        } else thinsolve <- solve(H_dH, BLOB$LZtWX) ## thin result, needed for Xt_X
        Xt_X <- crossprod(BLOB$LZtWX, thinsolve)
        XX_D <- BLOB$XtWX
        diag(XX_D) <- diag(XX_D)*(1+damping) ## adds D1=damping*diag(XX_D)
        #
        dbeta_eta <- solve(XX_D-Xt_X , dbeta_rhs)[,1]
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
    if (TRUE) {
      if (is.null(BLOB$dG)) {  
        tmp <- Matrix::solve(BLOB$chol_Q, as(BLOB$G_CHMfactor,"sparseMatrix")) ## (Hfac)
        tmp@x <- tmp@x^2
        BLOB$D_Md2hdv2 <- colSums(tmp)
        ## convert diag perturb of Md2hdv2 into a non-diag perturb of G :
        ## since H=invL_Q G t(invL_Q),  dG= L_Q dH t(L_Q) 
        BLOB$dG <- .ZWZtwrapper(BLOB$chol_Q , BLOB$D_Md2hdv2)
      }
      dampdG <- (damping*BLOB$dG) ## not diagonal...
      G_dG <- drop0(BLOB$Gmat + dampdG) ## probably not so sparse...
      rhs <- BLOB$chol_Q %*% grad_v
      rhs <- Matrix::solve(G_dG,rhs)
      dv_h <- drop(crossprod(BLOB$chol_Q,rhs))
      dampDpD_2 <- damping * BLOB$D_Md2hdv2 ## the diagonal perturbation of the full H matrix for the return value
    } else { ## pas vu bcp de difference en terme de vitesse
      ## part of the problem is avoiding tall matrices from tall ZA (and X), but for squarish ZA, using solve(chol_Q,... ) looks complicated
      if (is.null(BLOB$Md2hdv2)) BLOB$Md2hdv2 <- as(tcrossprod(Matrix::solve(BLOB$chol_Q, as(BLOB$G_CHMfactor,"sparseMatrix"))),"matrix")
      H_dH <- BLOB$Md2hdv2 
      dampDpD_2 <- diag(H_dH)*damping
      diag(H_dH) <- diag(H_dH)+dampDpD_2
      dv_h <- solve(H_dH,  grad_v)
    }
    resu <- list(dVscaled= dv_h, dampDpD = dampDpD_2)
    return(resu)
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
  if (which=="beta_v_cov_from_wAugX") { ## for predVar 
    ## not called during the fit, so ideally we store the necessary info in the fit rather than use get_from_MME() 
    stop("Programming error. Use ad hoc .calc_beta_cov_info_spprec() function instead.")
  }
  stop(cat("'which=\"",which,"\"' is invalid."))
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
