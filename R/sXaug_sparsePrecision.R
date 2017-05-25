
def_AUGI0_ZX_sparsePrecision <- function(AUGI0_ZX,rho,lambda,w.resid) {
  resu <- structure(AUGI0_ZX,
                    BLOB=list2env(list(), parent=environment(.AUGI0_ZX_sparsePrecision)),
                    lambda=lambda,
                    pforpv=ncol(AUGI0_ZX$X.pv),
                    w.resid=w.resid,
                    rho=rho )
  class(resu) <- c("AUGI0_ZX_sparsePrecision",class(resu))
  return( resu ) 
}

.AUGI0_ZX_sparsePrecision <- function(AUGI0_ZX,which="",z=NULL,B=NULL,LAPACK=TRUE,rho,
                                      damping,LM_z) {
  BLOB <- attr(AUGI0_ZX,"BLOB") ## an environment defined by def_AUGI0_ZX_sparsePrecision
  BLOB$Q_CHMfactor <- AUGI0_ZX$envir$Q_CHMfactor ## may have been computed before call of def_AUGI0_ZX_sparsePrecision
  if (is.null(BLOB$Q_CHMfactor)) { 
    stop("is.null(BLOB$Q_CHMfactor)") ## should not occur in current code
    Qmat <- - attr(AUGI0_ZX,"rho")* AUGI0_ZX$adjMatrix
    diag(Qmat) <- diag(Qmat)+1
    BLOB$Qmat <- Qmat
    BLOB$Q_CHMfactor <- Matrix::Cholesky(BLOB$Qmat,LDL=FALSE,perm=FALSE)
    ## $envir$LMatrix should also be missing
  } else BLOB$Qmat <- AUGI0_ZX$envir$Qmat ## A L S O $envir$LMatrix USED BELOW
  if (is.null(BLOB$G_CHMfactor)) { ## cannot be precomputed
    WZ <- AUGI0_ZX$ZAfix
    WZ@x <- WZ@x * attr(AUGI0_ZX,"w.resid")[WZ@i+1L]
    BLOB$ZtW <- t(WZ)
    precisionMatrix <-  BLOB$Qmat/attr(AUGI0_ZX,"lambda") 
    BLOB$Gmat <- Matrix::drop0(BLOB$ZtW %*% AUGI0_ZX$ZAfix + precisionMatrix)
    BLOB$G_CHMfactor <- Matrix::Cholesky(BLOB$Gmat,LDL=FALSE,perm=FALSE) 
  }
  # ZAL <- AUGI0_ZX$ZAfix %*% LMatrix
  # Xaug <- Matrix::drop0(Matrix::rbind2(
  #   Matrix::cbind2(AUGI0_ZX$I,AUGI0_ZX$ZeroBlock),
  #   Matrix::cbind2(ZAL,AUGI0_ZX$X)
  # ))
  # wXaug <- diag(c(rep(1/sqrt(attr(AUGI0_ZX,"lambda")),10),sqrt(attr(AUGI0_ZX,"w.resid")))) %*% Xaug
  #
  if (which=="hatval_Z") {
    if (is.null(BLOB$hatval_Z_)) {
      if (is.null(BLOB$LL)) BLOB$LL <- solve(BLOB$G_CHMfactor, as(BLOB$Q_CHMfactor,"sparseMatrix"),system="L")
      lev_lambda <- BLOB$LL
      lev_lambda@x <- lev_lambda@x^2
      lev_lambda <- colSums(lev_lambda)/attr(AUGI0_ZX,"lambda")
      if (AUGI0_ZX$is_unitary_ZAfix) { ## then only diagonal values of invG matter
        lev_phi <- colSums(solve(BLOB$G_CHMfactor, system="L")^2)
        lev_phi <- drop(AUGI0_ZX$ZAfix %*% lev_phi)*attr(AUGI0_ZX,"w.resid")
      } else {
        sqrtwZ <- AUGI0_ZX$ZAfix
        sqrtwZ@x <- sqrtwZ@x * sqrt(attr(AUGI0_ZX,"w.resid"))[sqrtwZ@i+1L]
        lev_phi <- solve(BLOB$G_CHMfactor, t(sqrtwZ), system="L") ## solve(as(BLOB$G_CHMfactor,"sparseMatrix"), t(sqrtwZ))
        lev_phi@x <- lev_phi@x^2
        lev_phi <- colSums(lev_phi)
      }
      BLOB$hatval_Z_ <- list(lev_lambda=lev_lambda,lev_phi=lev_phi)
      # zut <- colSums((solve(as(BLOB$G_CHMfactor,"sparseMatrix") , cbind2(solve(t(LMatrix)),t(diag(sqrt(attr(AUGI0_ZX,"w.resid"))) %*% (AUGI0_ZX$ZAfix)))) )^2)
      # if (max(abs(BLOB$hatval_Z_-zut))>1e-6) stop("check2")
    }
    return(BLOB$hatval_Z_)
  }
  if (which=="solve_d2hdv2") {
    if (is.null(BLOB$inv_d2hdv2)) { ## FIXME possibly better implementation assuming B will alway be there ?
      if (is.null(BLOB$LL)) BLOB$LL <- solve(BLOB$G_CHMfactor, as(BLOB$Q_CHMfactor,"sparseMatrix"),system="L")
      BLOB$inv_d2hdv2 <- - crossprod(BLOB$LL) ## dense result but all sparse operations
      # seq_n_u_h <- seq(nrow(LMatrix))
      # zut <- -solve(Matrix::crossprod(wXaug)[seq_n_u_h,seq_n_u_h])
      # if (max(abs(BLOB$inv_d2hdv2-zut))>1e-6) stop("check2")
    }
    if ( is.null(B) ) {
      return(BLOB$inv_d2hdv2)
    } else {
      return(BLOB$inv_d2hdv2 %*% B) ## FIXME (?) method without comput inv_d2hdv2 ? But BLOB$Q_CHMfactor needed elsewhere
    }
  }
  if (which=="half_logdetQ") {
    half_logdetQ <-  Matrix::determinant(BLOB$Q_CHMfactor)$modulus[1] 
    return(half_logdetQ)
  }
  if (which=="half_logdetG") {
    return(Matrix::determinant(BLOB$G_CHMfactor)$modulus[1])
  }
  # cases where Sigsolve is needed ## Sig is the Cov(response), tcrossprod(ZAL)/lambda+diag(1/w.resid) 
  Sigsolve <- function(rhs) {
    v <- Matrix::solve(BLOB$G_CHMfactor, BLOB$ZtW %*% rhs)
    invM.z <- attr(AUGI0_ZX,"w.resid")*(rhs - AUGI0_ZX$ZAfix %*% v)
    return(invM.z) ## where M=Z.solve(precMat).Zt+diag(~phi)
  } # solve(AUGI0_ZX$ZAfix %*% solve(precisionMatrix) %*% t(AUGI0_ZX$ZAfix)+diag(1/attr(AUGI0_ZX,"w.resid"))) = Sigsolve(diag(x=ncol(ZtW)))
  X.pv <- AUGI0_ZX$X.pv
  ## Solving for model coefficients:
  if ( ! is.null(z)) {
    if (is.null(BLOB$inv_d2hdv2)) { ## not yet computed in PQL/L FIXME: optimize which=="solve_d2hdv2"?
      if (is.null(BLOB$LL)) BLOB$LL <- solve(BLOB$G_CHMfactor, as(BLOB$Q_CHMfactor,"sparseMatrix"),system="L")
      BLOB$inv_d2hdv2 <- - crossprod(BLOB$LL) ## dense result of sparse operations
    }
    # computation a(1) here for a(0)=0 (LeeL12 p. 963 col 1 l 2)
    ## Raw inefficient version
    # Sig <- tcrossprod(ZAL)*attr(AUGI0_ZX,"lambda")+diag(1/attr(AUGI0_ZX,"w.resid")) ## Sigsolve(Sig)=I 
    # a <- drop(Sig %*% (attr(AUGI0_ZX,"w.resid") * z$sscaled))
    # invM.z <- Sigsolve(rhs=z$z1-a)
    ## there is SigSolve(Sig...) here! This simplify as
    invM.z <- Sigsolve(rhs=z$z1)-(attr(AUGI0_ZX,"w.resid") * z$sscaled)
    rhs <- Matrix::crossprod(X.pv,invM.z)
    # inv(X'.invSig.X).
    beta_eta <- solve(crossprod(X.pv, Sigsolve(rhs=X.pv)), rhs)[,1]
    ## Raw inefficient version
    # v_h <- - (BLOB$inv_d2hdv2 %*% (
    #   t(ZAL) %*% ((z$z1- X.pv %*% beta_eta)*attr(AUGI0_ZX,"w.resid"))+ z$z2/attr(AUGI0_ZX,"lambda") 
    # ))[,1]
    zz <- (z$z1- X.pv %*% beta_eta)*attr(AUGI0_ZX,"w.resid")
    ## [,1] bc crossprod(AUGI0_ZX$ZAfix,zz) is dense (dgeMatrix) and then dgC#dge siganture message 
    zz <- crossprod(AUGI0_ZX$envir$LMatrix, crossprod(AUGI0_ZX$ZAfix,zz)[,1] ) ## matrix times vector only
    v_h <- - (BLOB$inv_d2hdv2 %*% ( zz + z$z2/attr(AUGI0_ZX,"lambda") ))[,1] ## with z2=0 for GLMM...
    #
    # ZAL <- AUGI0_ZX$ZAfix %*% LMatrix
    # rhs_v <- z$z2+ drop((attr(AUGI0_ZX,"lambda") * t(ZAL)) %*% (z$sscaled * attr(AUGI0_ZX,"w.resid" )))
    # zut <- stats::lm.fit(as.matrix(wXaug),c(rhs_v,(z$z1-z$sscaled)*sqrt(attr(AUGI0_ZX,"w.resid"))))$coef
    # seq_n_u_h <- seq(nrow(LMatrix))
    # if (max(abs(v_h-zut[seq_n_u_h]))>1e-6) stop("check3")
    #
    return(list(v_h=v_h,beta_eta=beta_eta)) 
  }
  if (which=="beta_cov") { 
    beta_cov <- solve(crossprod(X.pv,Sigsolve(X.pv))) 
    colnames(beta_cov) <- colnames(X.pv)
    return(as.matrix(beta_cov))
  }
  #### REML:
  if (which=="logdet_sqrt_d2hdbeta2") { 
    M_d2hdbeta2 <- as.matrix(crossprod(X.pv,Sigsolve(X.pv))) ## FIXME should use X.Re for non-standard REML...
    logdet_sqrt_d2hdbv2 <- determinant(M_d2hdbeta2)$modulus[1]/2
    return(logdet_sqrt_d2hdbv2)
  }
  #### REML with inner estimation:
  if (which=="hatval") {
    if (is.null(BLOB$hatval_ZX)) {
      r12 <- solve(BLOB$G_CHMfactor, BLOB$ZtW %*% AUGI0_ZX$X.pv,system="L") ## solve(as(BLOB$G_CHMfactor,"sparseMatrix"), BLOB$ZtW %*% AUGI0_ZX$X.pv)
      r22 <- chol(ZtWZ(AUGI0_ZX$X.pv,attr(AUGI0_ZX,"w.resid"))-crossprod(r12))
      if (is.null(BLOB$LL)) BLOB$LL <- solve(BLOB$G_CHMfactor, as(BLOB$Q_CHMfactor,"sparseMatrix"),system="L")
      lev_lambda_z <- BLOB$LL
      lev_lambda_x <- - backsolve(r22,crossprod(r12,lev_lambda_z),transpose=TRUE)
      lev_lambda_z@x <- lev_lambda_z@x^2
      lev_lambda_x <- lev_lambda_x^2
      if ( is.matrix(lev_lambda_x)) {
        lev_lambda <- (colSums(lev_lambda_z)+colSums(lev_lambda_x))/attr(AUGI0_ZX,"lambda")
      } else {
        lev_lambda <- (colSums(lev_lambda_z)+lev_lambda_x)/attr(AUGI0_ZX,"lambda")
      }
      sqrtwZ <- AUGI0_ZX$ZAfix
      sqrtwZ@x <- sqrtwZ@x * sqrt(attr(AUGI0_ZX,"w.resid"))[sqrtwZ@i+1L]
      XtW <- t(AUGI0_ZX$X.pv *sqrt(attr(AUGI0_ZX,"w.resid")))
      if (AUGI0_ZX$is_unitary_ZAfix) {
        # no complete simplif for AUGI0_ZX$is_unitary_ZAfix bc lev_phi_x is not the diag of a Z...Zt matrix
        # but we can still use p*n matrices instead of r*n matrices to compute the lev_phi_x component
        g <- solve(BLOB$G_CHMfactor, system="L") 
        crossprod_r12_z_rows <- crossprod(crossprod(g, r12),t(sqrtwZ)) 
        lev_phi_z <- colSums(g^2)
        lev_phi_z <- drop(AUGI0_ZX$ZAfix %*% lev_phi_z)*attr(AUGI0_ZX,"w.resid")
      } else {
        lev_phi_z <- solve(BLOB$G_CHMfactor, t(sqrtwZ), system="L") ## solve(as(BLOB$G_CHMfactor,"sparseMatrix"), t(sqrtwZ))
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
  if (which=="LevMar_step") {  ## LM with diagonal perturbation as in NocedalW p. 266
    ## uses rather complex code to avoid solve(dense d2hdv2,...) but cannot avoid solve(rather dense G_dG)
    if (is.null(BLOB$dG)) {  
      tmp <- solve(BLOB$Q_CHMfactor, as(BLOB$G_CHMfactor,"sparseMatrix"),system="L")
      tmp@x <- tmp@x^2
      D_Md2hdv2 <- colSums(tmp)
      ## convert diag perturb of Md2hdv2 into a non-diag perturb of G :
      ## since H=invL_Q G t(invL_Q),  dG= L_Q dH t(L_Q) 
      BLOB$dG <- .ZWZtwrapper(as(BLOB$Q_CHMfactor,"sparseMatrix") , D_Md2hdv2)
      BLOB$ZtWX <- BLOB$ZtW %*% X.pv 
      BLOB$XtWX <- ZtWZ(X.pv,attr(AUGI0_ZX,"w.resid")) 
      BLOB$DpD <- c(D_Md2hdv2,diag(BLOB$XtWX))
    } 
    ## sequel recomputed for each new damping value...
    # if (is.null(BLOB$Gmat_dense)) BLOB$Gmat_dense <- as.matrix(BLOB$Gmat)
    # G_dG <- BLOB$Gmat_dense + damping*as.matrix(BLOB$dG) ## probably not so sparse...
    # nevertheless the sparse code is faster
    G_dG <- Matrix::drop0(BLOB$Gmat + damping*BLOB$dG) ## probably not so sparse...
    if (attr(AUGI0_ZX,"pforpv")>0L) {
      # RHS of eq for beta:
      XtWze <- crossprod(X.pv, attr(AUGI0_ZX,"w.resid") * LM_z$z1_eta) ## component Xt W (z1-eta) in Xt (invSigma=W-WZ()Z'W) z1
      Xt_ze <- crossprod(BLOB$ZtWX, Matrix::solve(G_dG, BLOB$ZtW %*% LM_z$z1_eta)) ## the other component 
      Xt_v <- crossprod(BLOB$ZtWX, Matrix::solve(G_dG, LM_z$phi_v)) 
      XtWs <- crossprod(X.pv, attr(AUGI0_ZX,"w.resid") * LM_z$sscaled)
      rhs <- XtWze-Xt_ze+Xt_v-XtWs
      # matrix for beta
      Xt_X <- crossprod(BLOB$ZtWX, Matrix::solve(G_dG, BLOB$ZtWX))
      XX_D <- BLOB$XtWX
      diag(XX_D) <- (1+damping)*diag(XX_D) ## adds D1=damping*diag(XX_D)
      dbeta_eta <- solve(XX_D-Xt_X , rhs)[,1]
      # 
      dv_rhs <- Matrix::solve(G_dG, 
                              - LM_z$phi_v + BLOB$ZtW %*% (LM_z$z1_eta - X.pv %*% dbeta_eta))
      dv_h <- crossprod(as(BLOB$Q_CHMfactor,"sparseMatrix"),dv_rhs)[,1]
      #v_h <- LM_z$v_h+dv_h
      resu <- list(dv_h=dv_h, dbeta_eta=dbeta_eta, dampDpD = damping*BLOB$DpD)
    } else {
      resu <- list(v_h=crossprod(as(BLOB$Q_CHMfactor,"sparseMatrix"), Matrix::solve(G_dG, BLOB$ZtW %*% LM_z$z1))[,1], 
                   beta_eta=numeric(0), dampDpD = damping*BLOB$DpD)
    }
    return(resu)
  }
  stop(cat("'which=\"",which,"\"' is invalid.",sep=""))
}

get_from_MME.AUGI0_ZX_sparsePrecision <- function(sXaug, ## really AUGI0_ZX: the generic's argnames are a bit too restrictive
                                                  which="",szAug=NULL,B=NULL,
                                                  damping=NULL, LM_z=NULL, ...) {
  resu <- switch(which,
                 "logdet_sqrt_d2hdv2" = { ## sqrt( logdet info(v) = -d2hdv2 ) 
                   if (is.null(d <- attr(sXaug$adjMatrix,"symSVD")$d)) {
                     half_logdetQ <- .AUGI0_ZX_sparsePrecision(AUGI0_ZX=sXaug,which="half_logdetQ")
                   } else half_logdetQ <- sum(log(abs(1-attr(sXaug,"rho")*d)))/2
                   - half_logdetQ + .AUGI0_ZX_sparsePrecision(AUGI0_ZX=sXaug,which="half_logdetG")
                 },
                 ## all other cases:
                 .AUGI0_ZX_sparsePrecision(AUGI0_ZX=sXaug,which=which,z=szAug,B=B,damping=damping,LM_z=LM_z)
  )
  return(resu)
  
}