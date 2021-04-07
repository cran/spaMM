# 'constructor' 
# from Xaug which already has a *scaled* ZAL 
def_sXaug_Matrix_QRP_CHM_scaled <- function(Xaug,weight_X,w.ranef,H_global_scale) {
  n_u_h <- length(w.ranef)
  Xrows <- n_u_h+seq(length(weight_X))
  ## Bates https://stat.ethz.ch/pipermail/r-help/2010-December/262365.html
  ## "Assignment of submatrices in a sparse matrix can be slow because there is so much checking that needs to be done."
  Xaug <- .Dvec_times_Matrix_lower_block(weight_X,Xaug,n_u_h)
  attr(Xaug, "get_from") <- "get_from_MME.sXaug_Matrix_QRP_CHM_scaled"
  attr(Xaug, "BLOB") <- list2env(list(), parent=environment(.sXaug_Matrix_QRP_CHM_scaled))
  attr(Xaug, "w.ranef") <- w.ranef
  attr(Xaug, "n_u_h") <- n_u_h # mandatory for all sXaug types
  attr(Xaug, "pforpv") <- ncol(Xaug)-n_u_h # mandatory for all sXaug types
  attr(Xaug, "weight_X") <- weight_X # new mandatory 08/2018
  attr(Xaug, "H_global_scale") <- H_global_scale
  ## cannot modify the 'class' attribute... => immediate clumsy code below... and how(.) reports method: dgCMatrix.
  return( Xaug ) 
}

# trace(get_from_MME, print=FALSE, tracer=quote(cat("'",crayon::yellow(which),"'")))
# trace(spaMM:::.sXaug_Matrix_QRP_CHM_scaled, print=FALSE, tracer=quote(cat("'",which,"'")))
#
.sXaug_Matrix_QRP_CHM_scaled <- function(sXaug,which="",szAug=NULL,B=NULL) {
  BLOB <- attr(sXaug,"BLOB") ## an environment
  if (is.null(BLOB$blob)) {
    BLOB$blob <- qr(sXaug) ##  Matrix::qr
    # sXaug = t(tQ) %*% R[,sP] but then also sXaug = t(tQ)[,sP'] %*% R[sP',sP] for any sP'
    BLOB$perm <- BLOB$blob@q + 1L
    BLOB$R_scaled <- qrR(BLOB$blob,backPermute = FALSE) # crossfac
    ##########################
    delayedAssign("sortPerm", sort.list(BLOB$perm), assign.env = BLOB ) # never NULL
    delayedAssign("seq_n_u_h", seq(attr(sXaug,"n_u_h")), assign.env = BLOB ) # repeatedly used for levM
    delayedAssign("sortPerm_u_h", BLOB$sortPerm[ BLOB$seq_n_u_h], assign.env = BLOB ) # repeatedly used for levM
    delayedAssign("u_h_cols_on_left", (max(BLOB$sortPerm_u_h)==attr(sXaug,"n_u_h")), assign.env = BLOB ) 
    delayedAssign("use_R_block", (BLOB$u_h_cols_on_left && FALSE) , assign.env = BLOB ) 
    # : TRUE for Rasch, bigranefs... 
    # : FALSE for cloglog (also test-ranCoefs but .sXaug_Matrix_QRP_CHM_scaled is used only in refit). 
    delayedAssign("invsqrtwranef", 1/sqrt(attr(sXaug,"w.ranef")), assign.env = BLOB )
    delayedAssign("t_Q_scaled", solve(t(BLOB$R_scaled),t(sXaug[,BLOB$perm])), assign.env = BLOB )  ## Matrix::solve ## this is faster than qr.Q and returns dgCMatrix rather than qr.Q -> dge ! 
    delayedAssign("CHMfactor_wd2hdv2w", {
      wd2hdv2w <- .crossprod(BLOB$R_scaled[,BLOB$sortPerm_u_h, drop=FALSE], allow_as_mat = FALSE ) # R_scaled is crossfac; CHMfactor ~ tcrossfac
      Cholesky(wd2hdv2w,LDL=FALSE, perm=FALSE ) ## perm=TRUE seems simple to implement except 
      #    for which="R_scaled_v_h_blob". Ascertaining when updating might be correct is not obvious. __F I X M E__
    }, assign.env = BLOB )
    delayedAssign("inv_factor_wd2hdv2w", { if (BLOB$use_R_block) { 
      sortPerm_u_h <- BLOB$sortPerm_u_h
      solve(t(BLOB$R_scaled))[sortPerm_u_h,sortPerm_u_h, drop=FALSE] # triangular solve remains sparse... 
      #  !! this is totally wrong if ! u_h_cols_on_left
      ## equivalences:
      # str(a <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w, system="A") ) # inv_d2dhdv2 # includes permuations
      # str(b <- crossprod(Matrix::solve(BLOB$CHMfactor_wd2hdv2w,system="L"))) # Matrix::solve(BLOB$CHMfactor_wd2hdv2w,system="L")
      # str(b <- drop0(Matrix::tcrossprod(solve(BLOB$R_scaled[BLOB$sortPerm_u_h,BLOB$sortPerm_u_h, drop=FALSE]))))
      # range(a-b)
    } else Matrix::solve(BLOB$CHMfactor_wd2hdv2w,system="L") } , assign.env = BLOB )  # crossfac
    delayedAssign("inv_d2hdv2", {        
      if (.is_evaluated("inv_factor_wd2hdv2w",BLOB)) {
        inv_d2hdv2 <- .Matrix_times_Dvec(.crossprod(BLOB$inv_factor_wd2hdv2w), BLOB$invsqrtwranef) 
      } else inv_d2hdv2 <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w, Diagonal(x=BLOB$invsqrtwranef), system="A") 
      inv_d2hdv2 <- .Dvec_times_Matrix( - BLOB$invsqrtwranef,inv_d2hdv2)
    }, assign.env = BLOB )
    delayedAssign("logdet_R_scaled_b_v", sum(log(abs(diag(x=BLOB$R_scaled)))), assign.env = BLOB )  # not .diagfast() on sparse matrix
    delayedAssign("logdet_R_scaled_v", {
      if (BLOB$use_R_block) { # but BLOB$logdet_R_scaled_v may not be requested in that case
        sum(log(abs(diag(x=BLOB$R_scaled)[BLOB$seq_n_u_h])))
      } else Matrix::determinant(BLOB$CHMfactor_wd2hdv2w)$modulus[1]
    }, assign.env = BLOB ) 
    delayedAssign("logdet_r22", {
      # the R's are H-scaled but r22 is H-unscaled... tricky!
      if (BLOB$u_h_cols_on_left) { # default
        sum(log(abs(diag(x=BLOB$R_scaled)[-BLOB$seq_n_u_h]))) - attr(sXaug,"pforpv")*log(attr(sXaug,"H_global_scale"))/2 # '-' not '<-'
      } else {
        BLOB$logdet_R_scaled_b_v - BLOB$logdet_R_scaled_v - attr(sXaug,"pforpv")*log(attr(sXaug,"H_global_scale"))/2 ## '-', not '<-'
      }
    } , assign.env = BLOB )
    delayedAssign("logdet_sqrt_d2hdv2", { sum(log(attr(sXaug,"w.ranef")))/2 + BLOB$logdet_R_scaled_v }, assign.env = BLOB )
    delayedAssign("hatval", { # colSums(t_Q_scaled@x^2)
      tmp <- BLOB$t_Q_scaled
      tmp@x <- tmp@x^2
      colSums(tmp)
    } , assign.env = BLOB )
    ##############################
    if ( ! is.null(szAug)) return(qr.coef(BLOB$blob,szAug))   
  } 
  # In .calc_sscaled_new, I compute the hatval_Z then solve_d2hdv2. 
  # CHMfactor_wd2hdv2w is needed in all cases for solve_d2hdv2.
  # It *may* also be useful to compute the hatval_Z. To compute hatval_Z
  # if (u_h_cols_on_left) {
  #   then t_Q_scaled has been computed to compute leverages
  #   neither inv_factor_wd2hdv2w nor CHMfactor_wd2hdv2w are yet computed
  # } else {
  #   either I precompute CHMfactor_wd2hdv2w (and optionally inv_factor_wd2hdv2w) and it seems logical to use them
  #   or I could use an .lmwithQR call as shown in the code below.
  #   The first option seems faster, but this is not commonly profileable. u_h_cols are generally on the left. An exception occurs
  #   in test-cloglog.R#57: fitme(cbind(Dead, Alive) ~ ... confirming the expectation
  # }
  if (FALSE) { # after efforts to find an Eigen syntax that works without storing a dense matrix in memory 
    # (solveInPlace may be the only way), we see that it's slow! (e.g test bigranefs)
    # https://stackoverflow.com/questions/53290521/eigen-solver-for-sparse-matrix-product
    BLOB$t_Q_scaled <- .Rcpp_backsolve_M_M(r=as(BLOB$R_scaled,"dgCMatrix"),x=t(sXaug[,BLOB$perm]),transpose=TRUE)
    # base::backsolve() calls as.matrix()... 
    # .tcrossprod(t(solve(BLOB$R_scaled))[,BLOB$sortPerm], sXaug)
    # 
    #   We will use:
    # BLOB$t_Q_scaled <- solve(t(BLOB$R_scaled),t(sXaug[,BLOB$perm])) ## Matrix::solve ## this is faster than qr.Q and returns dgCMatrix rather than qr.Q -> dge 
  } 
  if ( ! is.null(szAug)) {
    if (.is_evaluated("t_Q_scaled", BLOB)) {
      return(solve(BLOB$R_scaled, BLOB$t_Q_scaled %*% szAug)[BLOB$sortPerm,,drop=FALSE]) ## Matrix::solve
    } else {
      return(qr.coef(BLOB$blob,szAug)) # avoid t_Q_computation there ## Matrix::qr.coef
    }
  }
  # ELSE 
  if (which=="Qt_leftcols*B") {
    return(solve(t(BLOB$R_scaled), .crossprod(sXaug[-BLOB$seq_n_u_h,BLOB$perm], B)))
  }
  if (which=="Mg_solve_g") {
    seq_n_u_h <- BLOB$seq_n_u_h
    rhs <- B
    rhs[seq_n_u_h] <- BLOB$invsqrtwranef * rhs[seq_n_u_h]
    rhs <- Matrix::solve(BLOB$R_scaled,rhs[BLOB$perm])
    return(sum(rhs^2))
  } 
  if (which=="Mg_invH_g") {
    rhs <- BLOB$invsqrtwranef * B
    if (.is_evaluated("inv_factor_wd2hdv2w",BLOB)) { ## can be assigned elsewhere by lev_lambda <- BLOB$inv_factor_wd2hdv2w <- ....
      rhs <- BLOB$inv_factor_wd2hdv2w %*% rhs
    } else rhs <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w,rhs,system="L")
    return(sum(rhs^2))
  } 
  if (which=="Mg_invXtWX_g") { ## for which_LevMar_step="b", not currently used
    seq_n_u_h <- BLOB$seq_n_u_h
    if (is.null(BLOB$XtWX)) BLOB$XtWX <- .crossprod(sXaug[-seq_n_u_h,-seq_n_u_h])
    Mg_invXtWX_g <- crossprod(B,solve(BLOB$XtWX,B))[1L] # [1L] drops possible Matrix class...
    return(Mg_invXtWX_g)
  }
  if (which=="hatval_Z") { ## sscaled Pdiag
    if (is.null(BLOB$hatval_Z_)) {
      # X[,cols] = Q R P[,cols] = Q q r p => t(Q q) given by:
      n_u_h <- attr(sXaug,"n_u_h")
      if (BLOB$u_h_cols_on_left) { # Rasch... bigranefs... and that is remarkably fast compared to alternative software
        # as well as faster than the alternative code when u_h cols are not on the left. 
        # So, forming inv_factor_wd2hdv2w seems slow compared to this + only using $CHMfactor_wd2hdv2w
        tmp_t_Qq_scaled <- BLOB$t_Q_scaled[BLOB$seq_n_u_h,]
        tmp_t_Qq_scaled@x <- tmp_t_Qq_scaled@x^2
        tmp <- colSums(tmp_t_Qq_scaled)
        phipos <- n_u_h+seq_len(nrow(sXaug)-n_u_h)
        BLOB$hatval_Z_ <-  list(lev_lambda=tmp[-phipos],lev_phi=tmp[phipos])
        # previously to v2.1.46 there were comments showing slow code using qr.qy 
      } else if (FALSE) { #  correct but only for pedagogy. 
        #  Cf comments above. So forming inv_factor_wd2hdv2w, albeit slow,, is faster than a new QR factorization. 
        tq <- drop0(.lmwithQR(as.matrix(BLOB$R_scaled[,BLOB$sortPerm][,BLOB$seq_n_u_h]) ,yy=NULL,returntQ=TRUE,returnR=FALSE)$t_Q_scaled,
                    tol=1e-16)
        tmp_t_Qq_scaled <- tq %*% BLOB$t_Q_scaled
        tmp_t_Qq_scaled@x <- tmp_t_Qq_scaled@x^2
        tmp <- colSums(tmp_t_Qq_scaled)
        phipos <- n_u_h+seq_len(nrow(sXaug)-n_u_h)
        BLOB$hatval_Z_ <-  list(lev_lambda=tmp[-phipos],lev_phi=tmp[phipos])
      } else {
        # X[,cols] = Q R P[,cols] = Q q r p => t(Q q) given by:
        #t_Qq_scaled <- solve(BLOB$CHMfactor_wd2hdv2w, ## likely bottleneck for large data 
        #                          t(sXaug[, seq_len(attr(sXaug,"n_u_h"))[BLOB$perm_R_v] ]),system="L")
        ## 
        ## t(sXaug[,u_h cols]= (I, scaled t(ZAL)) i.e. is scaled such that the left block is an identity matrix, so we can work 
        ## on two separate blocks if the Cholesky is not permuted. Then 
        lev_lambda <- BLOB$inv_factor_wd2hdv2w
        lev_lambda@x <- lev_lambda@x^2
        lev_lambda <- colSums(lev_lambda)
        phipos <- (n_u_h+1L):nrow(sXaug)                 # ZAL block:
        lev_phi <- .tcrossprod(BLOB$inv_factor_wd2hdv2w, sXaug[phipos, BLOB$seq_n_u_h ], allow_as_mat = FALSE) # Matrix::solve(BLOB$CHMfactor_wd2hdv2w, t(sXaug[phipos, seq_len(n_u_h) ]),system="L") ## fixme the t() may still be costly
        # :if QRmethod is forced to "sparse" on mathematically dense matrices, we may reach this code yielding a *m* atrix lev_phi unless allow_as_mat = FALSE
        lev_phi@x <- lev_phi@x^2
        lev_phi <- colSums(lev_phi)
        BLOB$hatval_Z_ <-  list(lev_lambda=lev_lambda,lev_phi=lev_phi)
      }
    }
    return(BLOB$hatval_Z_)
  } 
  if (which=="solve_d2hdv2") {
    # don't forget that the factored matrix is not the augmented design matrix ! hence w.ranef needed here
    #if (is.null(BLOB$sortPerm_R_v)) BLOB$sortPerm_R_v <- sort.list(BLOB$perm_R_v)
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
        # In .calc_sscaled_new, I compute the hatval_Z then solve_d2hdv2. To compute hatval_Z
        # if (u_h_cols_on_left) {
        #   then t_Q_scaled has been computed to compute leverages
        #   neither inv_factor_wd2hdv2w nor CHMfactor_wd2hdv2w are yet computed
        # } else {
        #   CHMfactor_wd2hdv2w and inv_factor_wd2hdv2w has been evaluated to compute leverages.
        # }
        # Thus, when I reach here, I may or may not have inv_factor_wd2hdv2w available
        if ( BLOB$use_R_block || .is_evaluated("inv_factor_wd2hdv2w",BLOB)) {
          rhs <- .crossprod(BLOB$inv_factor_wd2hdv2w, drop(BLOB$inv_factor_wd2hdv2w %*% rhs)) # typical case when solve_d2hdv2 follows hatval_Z in .calc_sscaled_new()
        } else rhs <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w,rhs,system="A") ## dge (if rhs is dense, or a vector), or dgC...
        if (not_vector) { ## is.matrix(rhs) is not the correct test 
          rhs <- .Dvec_times_m_Matrix(BLOB$invsqrtwranef,rhs)
        } else rhs <- BLOB$invsqrtwranef * rhs
        return( - rhs)
      }
    } 
  } 
  if (which=="hatval") {return(BLOB$hatval) }  ## used in get_hatvalues  ##  REML hatval computation (also named hatval_ZX)
  if (which=="R_scaled_blob") { ## used for LevMar
    if (is.null(BLOB$R_scaled_blob)) {
      tmp <- X <- BLOB$R_scaled[ , BLOB$sortPerm] ## crossfac
      tmp@x <- tmp@x^2
      BLOB$R_scaled_blob <- list(X=X, diag_pRtRp=colSums(tmp), 
                                 XDtemplate=structure(.XDtemplate(X),upperTri=FALSE))
    }
    return(BLOB$R_scaled_blob)
  } 
  if (which=="R_scaled_v_h_blob") {
    if (is.null(BLOB$R_scaled_v_h_blob)) {
      if (BLOB$use_R_block) {
        R_scaled_v_h <- t(BLOB$R_scaled[BLOB$sortPerm_u_h,BLOB$sortPerm_u_h, drop=FALSE]) ## the t() for .damping_to_solve... (fixme: if we could avoid t()...)
      } else {
        R_scaled_v_h <- t( as(BLOB$CHMfactor_wd2hdv2w,"sparseMatrix") ) ## the t() for .damping_to_solve... (fixme: if we could avoid t()...)
      }
      tmp <- R_scaled_v_h 
      tmp@x <- tmp@x^2
      diag_pRtRp_scaled_v_h <- colSums(tmp)
      BLOB$R_scaled_v_h_blob <- list(R_scaled_v_h=R_scaled_v_h,diag_pRtRp_scaled_v_h=diag_pRtRp_scaled_v_h,
                                     XDtemplate=structure(.XDtemplate(R_scaled_v_h), upperTri=NA))
    }
    return(BLOB$R_scaled_v_h_blob)
  } 
  if (which=="R_beta_blob") {
    if (is.null(BLOB$R_beta_blob)) {
      n_u_h <- attr(sXaug,"n_u_h")
      seq_n_u_h <- BLOB$seq_n_u_h
      X <- as.matrix(sXaug[-seq_n_u_h,-seq_n_u_h]) ## The following code assumes it is dense...
      R_beta <- .lmwithQR(X,yy=NULL,returntQ=FALSE,returnR=TRUE)$R_scaled
      diag_pRtRp_beta <-  colSums(R_beta^2)
      BLOB$R_beta_blob <- list(R_beta=R_beta,diag_pRtRp_beta=diag_pRtRp_beta,  
                               XDtemplate=structure(.XDtemplate(R_beta), upperTri=TRUE))
    }
    return(BLOB$R_beta_blob)
  # } else if (which=="sortPerm") { 
  #   return(BLOB$sortPerm)
  } 
  if (which=="logdet_sqrt_d2hdv2") { return(BLOB$logdet_sqrt_d2hdv2)} 
  if (which=="logdet_r22") { return(BLOB$logdet_r22) } 
  if (which=="t_Q_scaled") { return(BLOB$t_Q_scaled) } ## used in get_hatvalues for nonstandard REML case
  if (which=="logdet_R_scaled_b_v") { return(BLOB$logdet_R_scaled_b_v)} 
  if (which=="beta_cov_info_from_sXaug") {
    return(.calc_beta_cov_info_from_sXaug(BLOB=BLOB, sXaug=sXaug, tcrossfac=solve(BLOB$R_scaled)))
  } 
  if (which=="beta_cov_info_from_wAugX") { ## using a weighted Henderson's augmented design matrix, not a true sXaug  
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
  } 
  if (which=="d2hdv2") { ## does not seem to be used
    #   # don't forgetthat the factored matrix is not the augmented design matrix ! hence w.ranef needed here
    #   w.ranef <- attr(sXaug,"w.ranef")
    #   w_R_R_v <- .Matrix_times_Dvec(BLOB$R_R_v,sqrt(w.ranef)[BLOB$perm_R_v])
    #   if (is.null(BLOB$sortPerm_R_v)) BLOB$sortPerm_R_v <- sort.list(BLOB$perm_R_v)
    #   BLOB$d2hdv2 <- - Matrix::crossprod(w_R_R_v)[BLOB$sortPerm_R_v,BLOB$sortPerm_R_v]
    #   return(BLOB$d2hdv2)
    stop("d2hdv2 requested ")
  } 
  stop("invalid 'which' value.")
}

# trace("get_from_MME.sXaug_Matrix_QRP_CHM_scaled",print=FALSE, tracer=quote(print(which)),exit=quote(str(resu)))
get_from_MME.sXaug_Matrix_QRP_CHM_scaled <- function(sXaug,which="",szAug=NULL,B=NULL,
                                      damping, LMrhs, ...) {
  resu <- switch(which,
                 "LevMar_step" = {
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
