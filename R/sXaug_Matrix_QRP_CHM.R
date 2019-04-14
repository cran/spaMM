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
                    weight_X=weight_X, # new mandatory 08/2018
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
    BLOB$u_h_cols_on_left <- (max(BLOB$sortPerm[ seq_n_u_h ])==n_u_h)
    if ( ! is.null(szAug)) return(qr.coef(BLOB$blob,szAug))   
  } 
  if ( is.null(BLOB$t_Q_scaled) && 
       ( which %in% c("t_Q_scaled","hatval") || (which=="hatval_Z" && BLOB$u_h_cols_on_left)
       ) ) {
    BLOB$t_Q_scaled <- solve(t(BLOB$R_scaled),t(sXaug[,BLOB$perm])) ## Matrix::solve ## this is faster than qr.Q and returns dgCMatrix rather than qr.Q -> dge !
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
       (which %in% c("logdet_R_scaled_v","hatval_Z","d2hdv2","solve_d2hdv2","R_scaled_v_h_blob"))
  ) {
    seq_n_u_h <- seq_len(attr(sXaug,"n_u_h"))
    wd2hdv2w <- crossprod(BLOB$R_scaled[,BLOB$sortPerm[ seq_n_u_h ]] ) # Matrix::crossprod
    BLOB$CHMfactor_wd2hdv2w <- Cholesky(wd2hdv2w,LDL=FALSE,
                                                perm=FALSE ) ## perm=FALSE useful for leverage computation as explained below
  }
  ## return()'s
  if (which=="Mg_solve_g") {
    if (is.null(BLOB$invsqrtwranef)) BLOB$invsqrtwranef <- 1/sqrt(attr(sXaug,"w.ranef")) 
    seq_n_u_h <- seq_len(attr(sXaug,"n_u_h"))
    rhs <- B
    rhs[seq_n_u_h] <- BLOB$invsqrtwranef * rhs[seq_n_u_h]
    rhs <- Matrix::solve(BLOB$R_scaled,rhs[BLOB$perm],system="L")
    return(sum(rhs^2))
  } else if (which=="Mg_invH_g") {
    if (is.null(BLOB$invsqrtwranef)) BLOB$invsqrtwranef <- 1/sqrt(attr(sXaug,"w.ranef")) 
    rhs <- BLOB$invsqrtwranef * B
    rhs <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w,rhs,system="L")
    return(sum(rhs^2))
  } else if (which %in% c("solve_d2hdv2")) {
    # don't forget that the factored matrix is not the augmented design matrix ! hence w.ranef needed here
    #if (is.null(BLOB$sortPerm_R_v)) BLOB$sortPerm_R_v <- sort.list(BLOB$perm_R_v)
    if (is.null(BLOB$invsqrtwranef)) BLOB$invsqrtwranef <- 1/sqrt(attr(sXaug,"w.ranef")) 
    if (is.null(B)) {
      if (is.null(BLOB$inv_d2hdv2)) {
        # if (is.null(BLOB$sortPerm_R_v)) BLOB$sortPerm_R_v <- sort.list(BLOB$perm_R_v)
        # w_R_R_v <- .Matrix_times_Dvec(BLOB$R_R_v,sqrt(w.ranef)[BLOB$perm_R_v])
        # BLOB$inv_d2hdv2 <- - Matrix::chol2inv(w_R_R_v)[BLOB$sortPerm_R_v,BLOB$sortPerm_R_v] 
        inv_d2hdv2 <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w, Diagonal(x=BLOB$invsqrtwranef)) 
        BLOB$inv_d2hdv2 <- - .Dvec_times_Matrix(BLOB$invsqrtwranef,inv_d2hdv2)
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
      if (BLOB$u_h_cols_on_left) {
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
        lev_lambda <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w,system="L")
        lev_lambda@x <- lev_lambda@x^2
        lev_lambda <- colSums(lev_lambda)
        n_u_h <- attr(sXaug,"n_u_h")
        phipos <- (n_u_h+1):nrow(sXaug)
        lev_phi <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w, t(sXaug[phipos, seq_len(n_u_h) ]),system="L") ## fixme the t() may still be costly
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
      BLOB$R_scaled_blob <- list(X=X, diag_pRtRp=colSums(tmp))
    }
    return(BLOB$R_scaled_blob)
  } else if (which=="R_scaled_v_h_blob") {
    if (is.null(BLOB$R_scaled_v_h_blob)) {
      R_scaled_v_h <- t( as(BLOB$CHMfactor_wd2hdv2w,"sparseMatrix") ) ## the t() for .damping_to_solve... (fixme: if we could avoid t()...)
      tmp <- R_scaled_v_h 
      tmp@x <- tmp@x^2
      diag_pRtRp_scaled_v_h <- colSums(tmp)
      BLOB$R_scaled_v_h_blob <- list(R_scaled_v_h=R_scaled_v_h,diag_pRtRp_scaled_v_h=diag_pRtRp_scaled_v_h)
    }
    return(BLOB$R_scaled_v_h_blob)
  # } else if (which=="R_beta_blob") {
  #   if (is.null(BLOB$R_beta_blob)) {
  #     n_u_h <- attr(sXaug,"n_u_h")
  #     seq_n_u_h <- seq(n_u_h)
  #     X <- as.matrix(sXaug[-seq_n_u_h,-seq_n_u_h]) ## follwoing code assuming it is dense...
  #     R_beta <- .lmwithQR(X,yy=NULL,returntQ=FALSE,returnR=TRUE)$R_scaled
  #     diag_pRtRp_beta <-  colSums(R_beta^2)
  #     BLOB$R_beta_blob <- list(R_beta=R_beta,diag_pRtRp_beta=diag_pRtRp_beta)
  #   }
  #   return(BLOB$R_beta_blob)
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
    if (TRUE) {
      tcrossfac_v_beta_cov <- solve(BLOB$R_scaled) # *M*atrix
      tPmat <- sparseMatrix(seq_along(BLOB$sortPerm), BLOB$sortPerm, x=1)
      pforpv <- attr(sXaug,"pforpv")
      X_scaling <- sqrt(rep(attr(sXaug,"H_global_scale"),pforpv))
      if ( ! is.null(B)) X_scaling <- X_scaling/B
      diagscalings <- c(1/sqrt(attr(sXaug,"w.ranef")), X_scaling)
      tPmat <- .Dvec_times_Matrix(diagscalings, tPmat) # Pmat <- .Matrix_times_Dvec(Pmat,diagscalings)
      tcrossfac_v_beta_cov <- as.matrix(tPmat %*% tcrossfac_v_beta_cov)
      rownames(tcrossfac_v_beta_cov) <- colnames(sXaug) ## necessary for summary.HLfit, already lost in BLOB$R_scaled
      seqp <- seq_len(pforpv)
      beta_pos <- attr(sXaug,"n_u_h")+seqp
      beta_v_order <- c(beta_pos,seq(attr(sXaug,"n_u_h")))
      tcrossfac_beta_v_cov <- tcrossfac_v_beta_cov[beta_v_order,,drop=FALSE]
      beta_cov <- .tcrossprod(tcrossfac_beta_v_cov[seqp,,drop=FALSE])
      return( list(beta_cov=beta_cov, 
                   #beta_v_cov=beta_v_cov,
                   tcrossfac_beta_v_cov=tcrossfac_beta_v_cov) )
    }
    ########################
    if (FALSE) {
      ## First version quite slow (test-negbin1): chol2inv, [BLOB$sortPerm,BLOB$sortPerm] on spare matrix, .Dvec_times_matrix 
      beta_cov <- Matrix::chol2inv(BLOB$R_scaled) ## actually augmented v_beta_cov as the following shows
      pforpv <- attr(sXaug,"pforpv")
      v_beta_cov <- as.matrix(beta_cov[BLOB$sortPerm,BLOB$sortPerm])
      diagscalings <- sqrt(c(1/attr(sXaug,"w.ranef"),rep(attr(sXaug,"H_global_scale"),pforpv)))
      v_beta_cov <- .Dvec_times_matrix(diagscalings, .m_Matrix_times_Dvec(v_beta_cov, diagscalings))
    } else {
      ## Second version faster, but could we avoid as.matrix ? (and t() ?). If not as.matrix, avoid [.,.] ?
      Pmat <- sparseMatrix(seq_along(BLOB$perm), BLOB$perm, x=1)
      pforpv <- attr(sXaug,"pforpv")
      diagscalings <- sqrt(c(1/attr(sXaug,"w.ranef"),rep(attr(sXaug,"H_global_scale"),pforpv)))
      Pmat <- .Matrix_times_Dvec(Pmat,diagscalings)
      v_beta_cov <- as.matrix(crossprod(solve(t(BLOB$R_scaled),Pmat)))
    }
    colnames(v_beta_cov) <- colnames(sXaug) ## necessary for summary.HLfit
    beta_pos <- attr(sXaug,"n_u_h")+seq_len(pforpv)
    beta_v_order <- c(beta_pos,seq(attr(sXaug,"n_u_h"))) ## 
    return( list(beta_cov=v_beta_cov[beta_pos,beta_pos,drop=FALSE], beta_v_cov=v_beta_cov[beta_v_order,beta_v_order]) )
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
                   list(dVscaled_beta=.damping_to_solve(X=R_scaled_blob$X, dampDpD=dampDpD, rhs=LMrhs), 
                        dampDpD = dampDpD) 
                 },
                 "LevMar_step_v_h" = {
                   ## FR->FR probably not the most elegant implementation 
                   R_scaled_v_h_blob <- .sXaug_Matrix_QRP_CHM_scaled(sXaug,which="R_scaled_v_h_blob")
                   dampDpD <- damping*R_scaled_v_h_blob$diag_pRtRp_scaled_v_h ## NocedalW p. 266
                   # Extend the X in X'X = P'R'RP: 
                   list(dVscaled = .damping_to_solve(X=R_scaled_v_h_blob$R_scaled_v_h, dampDpD=dampDpD, rhs=LMrhs), 
                        dampDpD = dampDpD) 
                 },
                 # "LevMar_step_beta" = {
                 #   R_beta_blob <- .sXaug_Matrix_QRP_CHM_scaled(sXaug,which="R_beta_blob")
                 #   dampDpD <- damping*R_beta_blob$diag_pRtRp_beta ## NocedalW p. 266
                 #   # Extend the X in X'X = R'R: 
                 #   Raug <- rbind(R_beta_blob$R_beta, diag(x=sqrt(dampDpD),nrow = length(dampDpD)))
                 #   RRblob <- .lmwithQRP(Raug,yy=NULL,returntQ=FALSE,returnR=TRUE)
                 #   RRR <- RRblob$R_scaled
                 #   RRsP <- sort.list(RRblob@perm) 
                 #   resu <- list(dVscaled=chol2inv(RRR)[RRsP,RRsP] %*% LMrhs, 
                 #                dampDpD = dampDpD) 
                 #   return(resu)
                 # } ,
                 ## all other cases:
                 .sXaug_Matrix_QRP_CHM_scaled(sXaug,which=which,szAug=szAug,B=B)
  )
  return(resu)
}
