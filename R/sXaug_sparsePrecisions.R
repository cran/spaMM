def_AUGI0_ZX_sparsePrecision <- function(AUGI0_ZX,corrPars,unique_w.ranef=NULL,w.ranef,cum_n_u_h,w.resid) {
  resu <- structure(AUGI0_ZX,
                    BLOB=list2env(list(), parent=environment(.AUGI0_ZX_sparsePrecision)),
                    unique_w.ranef=unique_w.ranef,
                    w.ranef=w.ranef,
                    cum_n_u_h=cum_n_u_h,
                    pforpv=ncol(AUGI0_ZX$X.pv),
                    w.resid=w.resid,
                    corrPars=corrPars )
  class(resu) <- c("AUGI0_ZX_sparsePrecision",class(resu))
  return( resu ) 
}

.reformat_Qmat_info <- function(BLOB, # a pointer to an envir !
                                envir, # another pointer to an envir ! 
                                corrPars=NULL) { ## provide the precision factors Q and their CHMfactor from different types of input
  ## used to compute Gmat (all non trivial elements needed) & chol_Q  
  # we want chol_Q to be (dtCMatrix: Csparse triangular) so that efficient solve methods can be used.
  # but eg bdiag(dtCMatrix, Diagonal Matrix) gives a dgCMatrix. 
  # All of chol_Q_list must be dtCMatrix so that bdiag() gives a dtCMatrix
  BLOB$precisionBlocks <- envir$precisionFactorList ## local copy so that AUGI0_ZX$envir$precisionFactorList is not modified
  BLOB$chol_Q_list <- vector("list",length(envir$precisionFactorList)) ## of L in LL' factorisation of precision factor 
  for (it in seq_len(length(envir$precisionFactorList))) {
    if( ! is.null(envir$precisionFactorList[[it]])) {
      if (inherits(envir$precisionFactorList[[it]],c("Matrix","matrix"))) { ## Q_CHMfactor NOT provided
        ## not clear that adjacency or AR1 can occur here
        if (envir$types[it]=="adjacency") { ## adjacency matrix provided
          Qmat <- -  corrPars$rho * envir$precisionFactorList[[it]] ## un seul rho pour l'instant 
          diag(Qmat) <- diag(Qmat)+1
          BLOB$precisionBlocks[[it]] <- Qmat
          BLOB$chol_Q_list[[it]] <- as(Matrix::Cholesky(envir$precisionFactorList[[it]],LDL=FALSE,perm=FALSE), "sparseMatrix")
        } else if (envir$types[it]=="AR1") {
          ## envir$precisionFactorList[[it]] should be a list hence this block should not be reached
          stop("reached here") 
        } else {## other models where a single matrix is provided : Qmatrix provided.
          BLOB$chol_Q_list[[it]] <- as(Matrix::Cholesky(envir$precisionFactorList[[it]],LDL=FALSE,perm=FALSE), "sparseMatrix")
        }
      } else { ## both Q and factorization provided 
        BLOB$chol_Q_list[[it]] <- envir$precisionFactorList[[it]]$chol_Q
        BLOB$precisionBlocks[[it]] <- envir$precisionFactorList[[it]]$Qmat ## indep of w.ranef (lambda) 
      }
    } ## ELSE NULL precisionFactorList element <-> chol_Q_list element prefilled
  }
  BLOB$chol_Q <- Matrix::bdiag(BLOB$chol_Q_list) 
  ## no return value :-)
}

.Sigsolve_sparsePrecision <- function(AUGI0_ZX, rhs) {
  BLOB <- attr(AUGI0_ZX,"BLOB")
  v <- Matrix::solve(BLOB$G_CHMfactor, BLOB$ZtW %*% rhs) 
  invM.z <- attr(AUGI0_ZX,"w.resid")*(rhs - AUGI0_ZX$ZAfix %*% v)  # W.rhs - W Z (v= L inv(...) L' Z' W rhs) = [W- WZL inv(...) L'Z'W] rhs
  return(invM.z) ## where M=Z.solve(precMat).Zt+diag(~phi)
} # solve(AUGI0_ZX$ZAfix %*% solve(precisionMatrix) %*% t(AUGI0_ZX$ZAfix)+diag(1/attr(AUGI0_ZX,"w.resid"))) = Sigsolve(diag(x=ncol(ZtW)))

# FIXME Matrix::solve not used economically in this code => long time of test-Rasch when sparse precision is used.
# It incolves the following calls ( * => mMatrix rhs ) :
# * solve(BLOB$BLOB$G_CHMfactor, chol_Q, "L")
# * solve(BLOB$BLOB$G_CHMfactor, t(sqrtwZ), "L")
# solve(BLOB$BLOB$G_CHMfactor, rhs, "A")
# solve(BLOB$BLOB$G_CHMfactor, (chol_Q %*% grad_v), "A")
# * solve(BLOB$BLOB$G_CHMfactor, ZtWX, "A")
# solve(BLOB$BLOB$G_CHMfactor, ZtWX %*% dbeta_eta, "A")
.AUGI0_ZX_sparsePrecision <- function(AUGI0_ZX,which="",z=NULL,B=NULL,
                                      damping,LM_z) {
  time666 <- Sys.time()
  verbose <- FALSE
  BLOB <- attr(AUGI0_ZX,"BLOB") ## an environment defined by def_AUGI0_ZX_sparsePrecision
  unique_w.ranef <- attr(AUGI0_ZX,"unique_w.ranef") 
  if ( is.null(unique_w.ranef)) w.ranef <- attr(AUGI0_ZX,"w.ranef") 
  cum_n_u_h <- attr(AUGI0_ZX,"cum_n_u_h")
  if (is.null(BLOB$G_CHMfactor)) { ## costly, cannot be precomputed
    if (is.null(AUGI0_ZX$envir$precisionFactorList)) stop("is.null(AUGI0_ZX$envir$precisionFactorList)")
    .reformat_Qmat_info(BLOB, AUGI0_ZX$envir, corrPars=attr(AUGI0_ZX,"corrPars")) ## includes computations of BLOB$chol_Q not function of w.ranef
    if (verbose) print(paste("aa",.timerraw(time666)))
    WZ <- AUGI0_ZX$ZAfix
    WZ <- .Dvec_times_Matrix(attr(AUGI0_ZX,"w.resid"),WZ) 
    BLOB$ZtW <- t(WZ)
    precisionBlocks <-  BLOB$precisionBlocks 
    # Either a precisionBlocks is that of a non-trivial corr mat of a gaussian ranef => the lambda and w.ranef are constant
    # or the lambda and w.ranef are not constant but the precisionBlocks is an identity matrix.
    # In both cases, precisionBlocks[[it]] * w.ranef[u.range] is a correct way of computing the precision block
    if ( ! is.null(unique_w.ranef)) {
      for (it in seq_len(length(precisionBlocks))) {
        precisionBlocks[[it]] <- precisionBlocks[[it]] * unique_w.ranef
      }
    } else for (it in seq_len(length(precisionBlocks))) {
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      precisionBlocks[[it]] <- precisionBlocks[[it]] * w.ranef[u.range]
    }
    if (length(precisionBlocks)>1L) {
      precisionMatrix <- forceSymmetric(do.call(Matrix::bdiag, precisionBlocks))
    } else precisionMatrix <- forceSymmetric(precisionBlocks[[1L]])
    ## FIXME more efficient code for case where isDiagonal(precisionMatrix) ?
    if (verbose) print(paste("ab",.timerraw(time666)))
    BLOB$Gmat <- Matrix::drop0(forceSymmetric(BLOB$ZtW %*% AUGI0_ZX$ZAfix) + precisionMatrix) ## depends on w.ranef and w.resid
    if (verbose) print(paste("ac",.timerraw(time666)))
    BLOB$G_CHMfactor <- Matrix::Cholesky(BLOB$Gmat,LDL=FALSE,perm=FALSE) ## costly
    ## with perm=TRUE G=P'LL'P and the P'L (non-triangular) factor is given by solve(<G_CHM>,as(<G_CHM>,"sparseMatrix"),system="Pt")
    if (verbose) print(paste("ad",.timerraw(time666)))
  }
  if ( is.null(BLOB$factor_inv_Md2hdv2)  && 
       #(which %in% c("logdet_R_scaled_v","hatval_Z","d2hdv2","solve_d2hdv2","R_scaled_v_h_blob"))
       (which %in% c("solve_d2hdv2","Mg_invH_g","hatval_Z","hatval"))
  ) {
    LL <- Matrix::solve(BLOB$G_CHMfactor, BLOB$chol_Q,system="L") ## dgCMatrix 
    BLOB$factor_inv_Md2hdv2 <- Matrix::drop0(LL,tol=.Machine$double.eps)
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
      if ( is.null(unique_w.ranef)) {
        lev_lambda <- .Matrix_times_Dvec(lev_lambda, w.ranef)
      } else lev_lambda <- lev_lambda * unique_w.ranef
      lev_lambda <- colSums(lev_lambda)
      if (AUGI0_ZX$is_unitary_ZAfix) { ## then only diagonal values of invG matter
        lev_phi <- colSums(Matrix::solve(BLOB$G_CHMfactor, system="L")^2)
        lev_phi <- drop(AUGI0_ZX$ZAfix %*% lev_phi)*attr(AUGI0_ZX,"w.resid")
      } else {
        #sqrtwZ <- AUGI0_ZX$ZAfix
        #sqrtwZ@x <- sqrtwZ@x * sqrt(attr(AUGI0_ZX,"w.resid"))[sqrtwZ@i+1L]
        sqrtwZ <- .Dvec_times_Matrix(sqrt(attr(AUGI0_ZX,"w.resid")), AUGI0_ZX$ZAfix)
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
    return(Matrix::determinant(BLOB$G_CHMfactor)$modulus[1])
  }
  # cases where Sigsolve is needed ## Sig is the Cov(response), tcrossprod(ZAL)/lambda+diag(1/w.resid) 
  X.pv <- AUGI0_ZX$X.pv
  ## Solving for model coefficients:
  if ( ! is.null(z)) { ## which=""
    if (verbose) print(paste("e",.timerraw(time666)))
    grad_v <- z$LMrhs[seq(ncol(BLOB$chol_Q))]
    #L_dv_term_from_grad_v <- Matrix::solve(BLOB$G_CHMfactor, BLOB$factor_inv_Md2hdv2 %*% grad_v, system="Lt")
    L_dv_term_from_grad_v <- Matrix::solve(BLOB$G_CHMfactor, BLOB$chol_Q %*% grad_v,system="A") # part from grad_h, to be completed
    if (attr(AUGI0_ZX,"pforpv")>0L) { ## if there are fixed effect coefficients to estimate
      ZtWX <- BLOB$ZtW %*% X.pv ## copy of ZtWXin BLOB not useful as it is used only once per instance of AUGI0_ZX (except when LevM is used)
      XtWX <- .ZtWZ(X.pv,attr(AUGI0_ZX,"w.resid")) 
      # rhs for beta
      grad_p_v <- z$LMrhs[-seq(ncol(BLOB$chol_Q))]
      dbeta_rhs <- Matrix::drop(grad_p_v - crossprod(ZtWX, L_dv_term_from_grad_v)) 
      # matrix for beta
      Xt_X <- crossprod(ZtWX, Matrix::solve(BLOB$G_CHMfactor, ZtWX, system="A"))
      dbeta_eta <- solve(XtWX-Xt_X , dbeta_rhs)[,1] ## wher Xt(W-_)X = X'SigSolve X
      # modif rhs for v_h :
      L_dv <- L_dv_term_from_grad_v - Matrix::solve(BLOB$G_CHMfactor, ZtWX %*% dbeta_eta, system="A")
      #v_h <- LM_z$v_h+dv_h
    } else dbeta_eta <- numeric(0)
    dv_h <- Matrix::drop(Matrix::crossprod(BLOB$chol_Q, Matrix::drop(L_dv))) ## inner drop() to avoid signature issue with dgeMatrix dv_rhs...
    return(list(dv_h=dv_h,dbeta_eta=dbeta_eta)) 
  }
  if (which=="beta_cov") { 
    beta_cov <- solve(crossprod(X.pv, .Sigsolve_sparsePrecision(AUGI0_ZX=AUGI0_ZX, rhs=X.pv))) 
    colnames(beta_cov) <- colnames(X.pv)
    return(as.matrix(beta_cov))
  }
  #### REML:
  if (which=="logdet_sqrt_d2hdbeta2") { 
    M_d2hdbeta2 <- as.matrix(crossprod(X.pv, .Sigsolve_sparsePrecision(AUGI0_ZX=AUGI0_ZX, rhs=X.pv))) ## FIXME should use X.Re for non-standard REML...
    logdet_sqrt_d2hdbv2 <- determinant(M_d2hdbeta2)$modulus[1]/2
    return(logdet_sqrt_d2hdbv2)
  }
  #### REML with inner estimation:
  if (which=="hatval") {
    if (is.null(BLOB$hatval_ZX)) {
      # We use the  Chol factorization of T'T as explained in the working doc
      # We know one of its block is the same as for the ML leverages, we compute two additional blocks: 
      r12 <- Matrix::solve(BLOB$G_CHMfactor, BLOB$ZtW %*% AUGI0_ZX$X.pv,system="L") ## solve(as(BLOB$G_CHMfactor,"sparseMatrix"), BLOB$ZtW %*% AUGI0_ZX$X.pv)
      r22 <- chol(.ZtWZ(AUGI0_ZX$X.pv,attr(AUGI0_ZX,"w.resid"))-crossprod(r12)) ## both lines as explained in working doc
      lev_lambda_z <- BLOB$factor_inv_Md2hdv2 ## ul block of R^{-T} as described in working doc
      lev_lambda_z@x <- lev_lambda_z@x^2 
      if (is.null(unique_w.ranef)) {
        lev_lambda_z <- .Matrix_times_Dvec(lev_lambda_z,w.ranef)
      } else {
        lev_lambda_z <- lev_lambda_z * unique_w.ranef
      }
      lev_lambda_z <- colSums(lev_lambda_z)
      #
      lev_lambda_x <- - backsolve(r22,crossprod(r12, BLOB$factor_inv_Md2hdv2),transpose=TRUE) ## ll block of R^{-T} as described in working doc
      lev_lambda_x <- lev_lambda_x^2
      ## one dimension is dropped when (?) X is a single column matrix
      if ( is.matrix(lev_lambda_x)) {  ## m ou M atrix
        if (is.null(unique_w.ranef)) {
          lev_lambda_x <- colSums(.m_Matrix_times_Dvec(lev_lambda_x,w.ranef))
        } else {
          lev_lambda_x <- colSums(lev_lambda_x) * unique_w.ranef
        }
      } else {
        if (is.null(unique_w.ranef)) {
          lev_lambda_x <- lev_lambda_x * w.ranef
        } else {
          lev_lambda_x <- lev_lambda_x * unique_w.ranef
        }
      }
      lev_lambda <- lev_lambda_z+lev_lambda_x
      #sqrtwZ <- AUGI0_ZX$ZAfix
      #sqrtwZ@x <- sqrtwZ@x * sqrt(attr(AUGI0_ZX,"w.resid"))[sqrtwZ@i+1L]
      sqrtwZ <- .Dvec_times_Matrix(sqrt(attr(AUGI0_ZX,"w.resid")), AUGI0_ZX$ZAfix)
      XtW <- t(AUGI0_ZX$X.pv *sqrt(attr(AUGI0_ZX,"w.resid")))
      if (AUGI0_ZX$is_unitary_ZAfix) {
        # no complete simplif for AUGI0_ZX$is_unitary_ZAfix bc lev_phi_x is not the diag of a Z...Zt matrix
        # but we can still use p*n matrices instead of r*n matrices to compute the lev_phi_x component
        g <- Matrix::solve(BLOB$G_CHMfactor, system="L") ## fixme rethink this.
        crossprod_r12_z_rows <- suppressMessages(crossprod(crossprod(g, r12),t(sqrtwZ))) 
        lev_phi_z <- colSums(g^2)
        lev_phi_z <- drop(AUGI0_ZX$ZAfix %*% lev_phi_z)*attr(AUGI0_ZX,"w.resid")
      } else {
        lev_phi_z <- Matrix::solve(BLOB$G_CHMfactor, t(sqrtwZ), system="L") ## solve(as(BLOB$G_CHMfactor,"sparseMatrix"), t(sqrtwZ))
        crossprod_r12_z_rows <- crossprod(r12,lev_phi_z)
        lev_phi_z@x <- lev_phi_z@x^2
        lev_phi_z <- colSums(lev_phi_z) 
      }
      lev_phi_x <- backsolve(r22, XtW - crossprod_r12_z_rows, transpose=TRUE)
      lev_phi_x <- lev_phi_x^2
      if ( is.matrix(lev_phi_x)) lev_phi_x <- colSums(lev_phi_x) 
      lev_phi <- lev_phi_z+lev_phi_x
      BLOB$hatval_ZX <- list(lev_lambda=lev_lambda,lev_phi=lev_phi)
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
        BLOB$ZtWX <- as(BLOB$ZtW %*% X.pv,"matrix") 
        BLOB$XtWX <- .ZtWZ(X.pv,attr(AUGI0_ZX,"w.resid")) 
        BLOB$DpD <- c(BLOB$D_Md2hdv2,diag(BLOB$XtWX)) ## store as diagonal perturbation for the retrun value
      } 
      ## sequel recomputed for each new damping value...
      # if (is.null(BLOB$Gmat_dense)) BLOB$Gmat_dense <- as.matrix(BLOB$Gmat)
      # G_dG <- BLOB$Gmat_dense + damping*as.matrix(BLOB$dG) ## probably not so sparse...
      # nevertheless the sparse code is faster
      grad_v <- LM_z$LMrhs[seq(ncol(BLOB$chol_Q))] 
      dampDpD_2 <- (damping*BLOB$dG) ## not diagonal...
      G_dG <- Matrix::drop0(BLOB$Gmat + dampDpD_2) ## probably not so sparse...
      L_dv_term_from_grad_v <- Matrix::drop(Matrix::solve(G_dG, Matrix::tcrossprod(BLOB$chol_Q, grad_v))) # part from grad_h, to be completed
      if (attr(AUGI0_ZX,"pforpv")>0L) { ## if there are fixed effect coefficients to estimate
        # rhs for beta
        grad_beta <- LM_z$LMrhs[-seq(ncol(BLOB$chol_Q))]
        dbeta_rhs <- grad_beta - crossprod(BLOB$ZtWX, L_dv_term_from_grad_v) # - dampDpD_1 * LM_z$v_h_beta$beta_eta one damps the gradient eq, not the one giving v_h_beta
        # matrix for beta
        Xt_X <- crossprod(BLOB$ZtWX, as(Matrix::solve(G_dG, BLOB$ZtWX),"matrix"))
        XX_D <- BLOB$XtWX
        dampDpD_1 <- damping*diag(XX_D)
        diag(XX_D) <- diag(XX_D)+dampDpD_1 ## adds D1=damping*diag(XX_D)
        #
        dbeta_eta <- solve(XX_D-Xt_X , dbeta_rhs)[,1]
        L_dv <- L_dv_term_from_grad_v - Matrix::drop(Matrix::solve(G_dG, BLOB$ZtWX %*% dbeta_eta)) # - dampDpD_2 %*% LM_z$v_h_beta$v_h
        dv_h <- Matrix::drop(Matrix::crossprod(BLOB$chol_Q, Matrix::drop(L_dv))) ## inner drop() to avoid signature issue with dgeMatrix dv_rhs...
      } else {
        dbeta_eta <- numeric(0)
        dv_h <- Matrix::drop(Matrix::crossprod(BLOB$chol_Q, Matrix::drop(L_dv))) ## inner drop() to avoid signature issue with dgeMatrix dv_rhs...
      }
    } else { ## may be faster but with useQR set to FALSE
      if (is.null(BLOB$DpD)) {  
        if (is.null(BLOB$LZtWX)) {
          BLOB$ZtWX <- BLOB$ZtW %*% X.pv 
          BLOB$LZtWX <- as(Matrix::solve(BLOB$chol_Q, BLOB$ZtWX),"matrix")
        }
        BLOB$XtWX <- .ZtWZ(X.pv,attr(AUGI0_ZX,"w.resid")) 
        BLOB$DpD <- c(BLOB$D_Md2hdv2,diag(BLOB$XtWX)) ## store as diagonal perturbation for the retrun value
      } 
      ## sequel recomputed for each new damping value...
      grad_v <- LM_z$LMrhs[seq(ncol(BLOB$chol_Q))]
      grad_beta <- LM_z$LMrhs[-seq(ncol(BLOB$chol_Q))]
      if ((useQR <- FALSE)) {
        if (is.null(BLOB$D_Md2hdv2)) {
          ## t() bc .damping_to_solve implictly uses crossprod ( as in solve(X'X)=solve(R'R) )
          BLOB$Hfac <- tmp <- t(Matrix::drop0(solve(BLOB$chol_Q, as(BLOB$G_CHMfactor,"sparseMatrix")))) ## probably not so sparse...
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
      if (attr(AUGI0_ZX,"pforpv")>0L) { ## if there are fixed effect coefficients to estimate
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
    return(list(dv_h=dv_h, dbeta_eta=dbeta_eta, dampDpD = damping*BLOB$DpD))
  }
  if (which =="LevMar_step_v_h") {  ## LM with diagonal perturbation as in NocedalW p. 266
    if (FALSE) { ## 
      if (is.null(BLOB$D_Md2hdv2)) {
        BLOB$Hfac <- tmp <- Matrix::drop0(Matrix::solve(BLOB$chol_Q, as(BLOB$G_CHMfactor,"sparseMatrix"))) ## probably not so sparse...
        tmp@x <- tmp@x^2
        BLOB$D_Md2hdv2 <- colSums(tmp)
      }
      dampDpD <- damping * BLOB$D_Md2hdv2
      dv_h <- .damping_to_solve(X=BLOB$Hfac, dampDpD=dampDpD, LM_z$LMrhs[seq(length(dampDpD))]) 
    } else {
      grad_v <- LM_z$LMrhs[seq(ncol(BLOB$chol_Q))]
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
        G_dG <- Matrix::drop0(BLOB$Gmat + dampdG) ## probably not so sparse...
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
    }
    resu <- list(dv_h= dv_h, dampDpD = dampDpD_2)
    return(resu)
  }
  stop(cat("'which=\"",which,"\"' is invalid.",sep=""))
}

get_from_MME.AUGI0_ZX_sparsePrecision <- function(sXaug, ## really AUGI0_ZX: the generic's argnames are a bit too restrictive
                                                  which="",szAug=NULL,B=NULL,
                                                  damping=NULL, LM_z=NULL, ...) {
  if (which=="LevMar_step" && damping==0L) {
    resu <- get_from_MME.AUGI0_ZX_sparsePrecision(sXaug=sXaug, szAug=LM_z)
  } else resu <- switch(which,
                 "logdet_sqrt_d2hdv2" = { ## sqrt( logdet info(v) = -d2hdv2 ) 
                   if (is.null(d <- attr(sXaug$adjMatrix,"symSVD")$d)) {
                     half_logdetQ <- .AUGI0_ZX_sparsePrecision(AUGI0_ZX=sXaug,which="half_logdetQ")
                   } else half_logdetQ <- sum(log(abs(1-attr(sXaug,"corrPars")$rho*d)))/2
                   - half_logdetQ + .AUGI0_ZX_sparsePrecision(AUGI0_ZX=sXaug,which="half_logdetG")
                 },
                 ## all other cases:
                 .AUGI0_ZX_sparsePrecision(AUGI0_ZX=sXaug,which=which,z=szAug,B=B,damping=damping,LM_z=LM_z)
  )
  return(resu)
  
}
