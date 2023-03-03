# 'constructor' 
# from Xaug which already has a *scaled* ZAL 
def_sXaug_Matrix_QRP_CHM_scaled <- function(Xaug,weight_X,w.ranef,H_global_scale, 
                                            nonSPD=NULL,
                                            force_QRP=NULL # ignored
                                            ) {
  n_u_h <- length(w.ranef)
  Xrows <- n_u_h+seq(length(weight_X))
  ## Bates https://stat.ethz.ch/pipermail/r-help/2010-December/262365.html
  ## "Assignment of submatrices in a sparse matrix can be slow because there is so much checking that needs to be done."
  Xaug <- .Dvec_times_Matrix_lower_block(weight_X,Xaug,n_u_h)
  attr(Xaug, "get_from") <- "get_from_MME.sXaug_Matrix_QRP_CHM_scaled"
  H_w.resid <- attr(weight_X,"H_w.resid")
  BLOB <- list2env(list(H_w.resid=H_w.resid,
                                      signs=attr(H_w.resid,"signs"),
                                      nonSPD=identical(nonSPD,TRUE)), parent=emptyenv())
  if (BLOB$nonSPD) {
    BLOB$WLS_mat_weights <- abs(H_w.resid)
  } else BLOB$WLS_mat_weights <- H_w.resid # same as in CHM_H methods when the two can be compared (all SPD cases => do not condition on BLOB$signs)
  BLOB$signs_in_WLS_mat <- FALSE
  
  attr(Xaug, "BLOB") <- BLOB
  attr(Xaug, "w.ranef") <- w.ranef
  attr(Xaug, "n_u_h") <- n_u_h # mandatory for all sXaug types
  attr(Xaug, "pforpv") <- ncol(Xaug)-n_u_h # mandatory for all sXaug types
  attr(Xaug, "weight_X") <- weight_X # new mandatory 08/2018
   # =>  something like sqrt(H_global_scale*{ w.resid[$w_resid], potentially corrected for Hobs})
   # so the { w.resid[$w_resid], potentially corrected for Hobs} = weight_X^2 / H_global_scale
  attr(Xaug, "H_global_scale") <- H_global_scale
  ## cannot modify the 'class' attribute... => immediate clumsy code below... and how(.) reports method: dgCMatrix.
  return( Xaug ) 
}

.calc_t_Q_scaled <- function(BLOB, sXaug) {
  # after efforts to find an Eigen syntax that works without storing a dense matrix in memory # (solveInPlace may be the only way), 
  #   we see that it's slow! (e.g test bigranefs)
  # https://stackoverflow.com/questions/53290521/eigen-solver-for-sparse-matrix-product
  #### BLOB$t_Q_scaled <- .Rcpp_backsolve_M_M(r=as(BLOB$R_scaled,"dgCMatrix"),x=t(sXaug[,BLOB$perm]),transpose=TRUE)
  # base::backsolve() calls as.matrix()... 
  # Do not try Matrix::qr.qty() on large RHS! (not memory efficient).
  # 
  if (eval(.spaMM.data$options$presolve_cond) && .calc_denseness(sXaug,relative = TRUE)>0.004 ) { # nested_Matern
    # Below, solve(t(BLOB$R_scaled),t(sXaug[,BLOB$perm])) is slow when the result is rather dense. 
    #   It calls .sortCsparse(.Call(dtCMatrix_sparse_solve, a, b)), and its difficult to see whether the .sortCSparse itself is the culprit
    # nested Matern has denseness = 0.00509... and presolving is better
    # bigranefs has denseness 0.0001865483 and presolving is worse
    #
    #Matrix::crossprod(BLOB$solve_R_scaled, t(sXaug[,BLOB$perm])) # Rather avoid transpose on the larger matrix:
    Matrix::tcrossprod(t(BLOB$solve_R_scaled), sXaug[,BLOB$perm])
  } else { # bigranefs
    ## Matrix::solve ## this is faster than qr.Q and returns dgCMatrix rather than qr.Q -> dge ! 
    solve(t(BLOB$R_scaled),t(sXaug[,BLOB$perm])) # a bit ugly, but (1) see comment wrt .Rcpp_backsolve_M_M() call;
    # (2) for dense vector LHS, there is a simple function cs_ltsolve() in suiteSparse, but otherwise see the complex code
    # of (suiteSparse) Cholesky/cholmod_spsolve as integrated in Matrix   (and not even handling the t()) 
  }
}

.calc_hatval_Z_u_h_cols_on_left <- function(BLOB, sXaug) {
  # previously to v2.1.46 there were comments showing slow code using qr.qy 
  n_u_h <- attr(sXaug,"n_u_h") 
  phipos <- seq(n_u_h+1L,nrow(sXaug)) # corresponding to rows of *Q*, which are never permuted
  # getting t_Qq_scaled is cheap when u_h_cols_on_left, so there in no point avoiding its computation as in .calc_Z_lev_[lambda|phi]
  tmp_t_Qq_scaled <- BLOB$t_Qq_scaled # on the left of Q => on the first rows of t_Q => and order among these rows does not matter after colSums(x^2) 
  if  (is.null(BLOB$signs)) {
    xx <- tmp_t_Qq_scaled@x 
    xx <- xx*xx # appears faster than accessing @x twice !
    tmp_t_Qq_scaled@x <- xx
    tmp <- colSums(tmp_t_Qq_scaled)
    hatval_Z_ <-  list(lev_lambda=tmp[-phipos],lev_phi=tmp[phipos])
  } else {
    tmp <- colSums(( BLOB$invIm2QtdQ_Z %*% tmp_t_Qq_scaled) * (tmp_t_Qq_scaled))
    hatval_Z_ <-  list(lev_lambda=tmp[-phipos],lev_phi=tmp[phipos]*BLOB$signs) 
  }
  hatval_Z_
}

.calc_Z_lev_lambda <- function(BLOB) {      
  t_Qq_lam_cols <- BLOB$inv_factor_wd2hdv2w # using sXaug unpermuted col order # think of R^-T %*% I-block-of-sXaug...
  if (is.null(BLOB$signs)) {
    xx <- t_Qq_lam_cols@x
    xx <- xx*xx
    t_Qq_lam_cols@x <- xx
    lev_lambda <- colSums(t_Qq_lam_cols)
  } else {
    lev_lambda <- colSums(( BLOB$invIm2QtdQ_Z %*% t_Qq_lam_cols) * t_Qq_lam_cols)
  }
  lev_lambda
}

.calc_Z_lev_phi <- function(BLOB, sXaug) {      
  n_u_h <- attr(sXaug,"n_u_h")
  phipos <- (n_u_h+1L):nrow(sXaug)                 # ZAL block:
  if (is.null(BLOB$signs)) {
    #lev_phi <- .tcrossprod(BLOB$inv_factor_wd2hdv2w, .leftcols_Csp(sXaug, n_u_h, keep_names=FALSE)[phipos, ], chk_sparse2mat = FALSE) # Matrix::solve(BLOB$CHMfactor_wd2hdv2w, t(sXaug[phipos, seq_len(n_u_h) ]),system="L") 
    # :if QRmethod is forced to "sparse" on mathematically dense matrices, we may reach this code yielding a *m* atrix lev_phi unless chk_sparse2mat = FALSE
    t_Qq_phi_cols <- .tcrossprod(BLOB$inv_factor_wd2hdv2w, sXaug[phipos, BLOB$seq_n_u_h ], chk_sparse2mat = FALSE) # Matrix::solve(BLOB$CHMfactor_wd2hdv2w, t(sXaug[phipos, seq_len(n_u_h) ]),system="L")
    xx <- t_Qq_phi_cols@x
    xx <- xx*xx
    t_Qq_phi_cols@x <- xx
    lev_phi <- colSums(t_Qq_phi_cols)
  } else { # $t_Qq_scaled used for $invIm2QtdQ_Z anyway
    t_Qq_phi_cols <- BLOB$t_Qq_scaled[,phipos] # .tcrossprod(BLOB$inv_factor_wd2hdv2w, sXaug[phipos, BLOB$seq_n_u_h ], chk_sparse2mat = FALSE) # Matrix::solve(BLOB$CHMfactor_wd2hdv2w, t(sXaug[phipos, seq_len(n_u_h) ]),system="L")
    lev_phi <- colSums(( BLOB$invIm2QtdQ_Z %*% t_Qq_phi_cols) * t_Qq_phi_cols)*BLOB$signs
  }
  lev_phi
}


.calc_inv_d2hdv2_QRP_CHM <- function(BLOB) {
   if (.is_evaluated("inv_factor_wd2hdv2w",BLOB)) {
    rhs <- .Matrix_times_Dvec(.crossprod(BLOB$inv_factor_wd2hdv2w), BLOB$invsqrtwranef) 
  } else rhs <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w, Diagonal(x=BLOB$invsqrtwranef), system="A") 
  .Dvec_times_Matrix( - BLOB$invsqrtwranef,rhs)
}

.calc_inv_d2hdv2_QRP_CHM_signs <- function(BLOB) {
  if (.is_evaluated("inv_factor_wd2hdv2w",BLOB)) {
    rhs <-  .Matrix_times_Dvec(BLOB$inv_factor_wd2hdv2w, BLOB$invsqrtwranef)
    rhs <- BLOB$invIm2QtdQ_Z %*% rhs
    rhs <- .crossprod(BLOB$inv_factor_wd2hdv2w, rhs) 
  } else {
    rhs <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w, Diagonal(x=BLOB$invsqrtwranef),system="Lt")
    rhs <- BLOB$invIm2QtdQ_Z %*% rhs
    rhs <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w, rhs, system="L")
  }
  .Dvec_times_Matrix( - BLOB$invsqrtwranef,rhs)
}

# trace(get_from_MME, print=FALSE, tracer=quote(cat("'",crayon::yellow(which),"'")))
# trace(spaMM:::.sXaug_Matrix_QRP_CHM_scaled, print=FALSE, tracer=quote(cat("'",which,"'")))
#
.sXaug_Matrix_QRP_CHM_scaled <- function(sXaug,which="",szAug=NULL,B=NULL) {
  BLOB <- attr(sXaug,"BLOB") ## an environment
  # attr(sXaug,"AUGI0_ZX") was always present in in .solve_IRLS_as_ZX(), being inherited from the Xaug argument of def_sXaug_Matrix_QRP_CHM_scaled
  # but was not present in .solve_v_h_IRLS(). Now it should be always present.
  if (is.null(BLOB$QRsXaug)) {
    BLOB$QRsXaug <- qr(sXaug) ##  Matrix::qr
    # sXaug = t(tQ) %*% R[,sP] but then also sXaug = t(tQ)[,sP'] %*% R[sP',sP] for any sP'
    BLOB$perm <- BLOB$QRsXaug@q + 1L
    BLOB$R_scaled <- qrR(BLOB$QRsXaug,backPermute = FALSE) # crossfac
    ##########################
    # the promises below have an evaluation environment. E.g.
    # ls(rlang:::promise_env("Z_lev_phi",BLOB)) lists  "B"        "BLOB"     "sXaug"    "szAug"    "wd2hdv2w" "which"   
    #              [wd2hdv2 bc CHMfactor_wd2hdv2w promise has been evaluated in BLOB!?]
    # and 
    # > rlang::env_parents(rlang:::promise_env("Z_lev_phi",BLOB))
    # [[1]] $ <env: namespace:spaMM>  # what's defined by spaMM
    # [[2]] $ <env: imports:spaMM>  # what's imported by spaMM
    # [[3]] $ <env: namespace:base>
    # [[4]] $ <env: global> 
    # but the contents of spaMM namespace... and the global envir are not saved in the object: 
    # when loaded in another session (with another spaMM version etc), this points to the local versions of all objects  
    delayedAssign("sortPerm", sort.list(BLOB$perm), assign.env = BLOB ) # never NULL
    delayedAssign("seq_n_u_h", seq(attr(sXaug,"n_u_h")), assign.env = BLOB ) # repeatedly used for levM
    delayedAssign("sortPerm_u_h", BLOB$sortPerm[ BLOB$seq_n_u_h], assign.env = BLOB ) # repeatedly used for levM
    delayedAssign("u_h_cols_on_left", (max(BLOB$sortPerm_u_h)==attr(sXaug,"n_u_h")), assign.env = BLOB ) 
    # : TRUE for Rasch, bigranefs... 
    # : FALSE for cloglog (also test-ranCoefs but .sXaug_Matrix_QRP_CHM_scaled is used only in refit). 
    delayedAssign("use_R_block", (BLOB$u_h_cols_on_left && FALSE) , assign.env = BLOB ) 
    delayedAssign("invsqrtwranef", 1/sqrt(attr(sXaug,"w.ranef")), assign.env = BLOB )
    delayedAssign("solve_R_scaled", 
      # solve(BLOB$R_scaled) is as(.Call(dtCMatrix_sparse_solve, a, .trDiagonal(n, unitri = FALSE)), "dtCMatrix")
      # as(.rawsolve(BLOB$R_scaled),"triangularMatrix"), 
      as(Matrix::solve(BLOB$R_scaled, attr(sXaug,"AUGI0_ZX")$trDiag),"triangularMatrix"), 
      assign.env = BLOB ) # Matrix::solve
    delayedAssign("t_Q_scaled", .calc_t_Q_scaled(BLOB, sXaug), assign.env = BLOB )
    # 
    ## to be used either when BLOB$u_h_cols_on_left or when there are $signs:
    delayedAssign("t_Qq_scaled", {
      if (BLOB$u_h_cols_on_left) {
        tmp_t_Qq_scaled <- BLOB$t_Q_scaled[BLOB$seq_n_u_h,] # on the left of Q => on the first rows of t_Q => and order among these rows does not matter after colSums(x^2) 
      } else .tcrossprod(BLOB$inv_factor_wd2hdv2w, sXaug[, BLOB$seq_n_u_h ], chk_sparse2mat = FALSE) # Matrix::solve(BLOB$CHMfactor_wd2hdv2w, t(sXaug[phipos, seq_len(n_u_h) ]),system="L")
    }, assign.env = BLOB )
    #
    delayedAssign("CHMfactor_wd2hdv2w", {
      wd2hdv2w <- .crossprod(BLOB$R_scaled[,BLOB$sortPerm_u_h, drop=FALSE], allow_as_mat = FALSE ) # R_scaled is crossfac; CHMfactor ~ tcrossfac
      Cholesky(wd2hdv2w,LDL=FALSE, perm=FALSE ) ## perm=TRUE seems 'simple'(*) to implement except 
      #    for which="R_scaled_v_h_blob". whether perm=TRUE (for updating) might be correct is not obvious:
      # There, R_scaled_v_h <- t( as(BLOB$CHMfactor_wd2hdv2w,"sparseMatrix") ) must be triangular and without permutation of the v's
      # (*) the SPPREC code seems better structured to reach the 'updateable' info.
      #  __F I X M E__? progress is not obvious: would need to make updateable info accessible, and even so it might not be useful.
    }, assign.env = BLOB )
    #
    delayedAssign("inv_factor_wd2hdv2w", {   # alwas unpermuted user's column order, notably bc CHMfactor_wd2hdv2w is Cholesky(permuted back to user order,,perm=FALSE)
      if (BLOB$use_R_block) { # always FALSE... 
        sortPerm_u_h <- BLOB$sortPerm_u_h
        # solve(t(BLOB$R_scaled))[sortPerm_u_h,sortPerm_u_h, drop=FALSE] # triangular solve remains sparse...
        t(BLOB$solve_R_scaled)[sortPerm_u_h,sortPerm_u_h, drop=FALSE] # triangular solve remains sparse... 
        #  !! this is totally wrong if ! u_h_cols_on_left
        ## equivalences:
        # str(a <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w, system="A") ) # inv_d2dhdv2 # includes permuations
        # str(b <- crossprod(Matrix::solve(BLOB$CHMfactor_wd2hdv2w,system="L"))) # Matrix::solve(BLOB$CHMfactor_wd2hdv2w,system="L")
        # str(b <- drop0(Matrix::tcrossprod(solve(BLOB$R_scaled[BLOB$sortPerm_u_h,BLOB$sortPerm_u_h, drop=FALSE]))))
        # range(a-b)
      } else Matrix::solve(BLOB$CHMfactor_wd2hdv2w,system="L", b=attr(sXaug,"AUGI0_ZX")$I) 
    } , assign.env = BLOB )  # crossfac # gives correct results for $signs...
    #
    delayedAssign("logdet_R_scaled_b_v", sum(log(abs(diag(x=BLOB$R_scaled)))), assign.env = BLOB )  # not .diagfast() on sparse matrix
    #
    delayedAssign("logdet_R_scaled_v", {
      # the tests led to use_R_block being FALSE in all cases. Gain is obviously not here, 
      # so it must be through making CHMfactor_wd2hdv2w available and efficient for other operations... not transparent.
      if (BLOB$use_R_block) { # but BLOB$logdet_R_scaled_v may not be requested in that case
        sum(log(abs(diag(x=BLOB$R_scaled)[BLOB$seq_n_u_h])))
      } else Matrix::determinant(BLOB$CHMfactor_wd2hdv2w)$modulus[1]
    }, assign.env = BLOB )
    #
    delayedAssign("logdet_r22", {
      # the R's are H-scaled but r22 is H-unscaled... tricky!
      if (BLOB$u_h_cols_on_left) { # default
        sum(log(abs(diag(x=BLOB$R_scaled)[-BLOB$seq_n_u_h]))) - attr(sXaug,"pforpv")*log(attr(sXaug,"H_global_scale"))/2 # '-' not '<-'
      } else {
        BLOB$logdet_R_scaled_b_v - BLOB$logdet_R_scaled_v - attr(sXaug,"pforpv")*log(attr(sXaug,"H_global_scale"))/2 ## '-', not '<-'
      }
    } , assign.env = BLOB )
    #
    if  (is.null(BLOB$signs)) {
      delayedAssign("inv_d2hdv2", .calc_inv_d2hdv2_QRP_CHM(BLOB), assign.env = BLOB )
      
      delayedAssign("logdet_sqrt_d2hdv2", { sum(log(attr(sXaug,"w.ranef")))/2 + BLOB$logdet_R_scaled_v }, assign.env = BLOB )
      #
      delayedAssign("hatval", { # colSums(t_Q_scaled@x^2)
        tmp <- BLOB$t_Q_scaled
        xx <- tmp@x
        xx <- xx*xx
        tmp@x <- xx
        colSums(tmp)
      } , assign.env = BLOB )
      
    } else {  # QR with signs -> typically nonSPD    
      #
      ## The WLS matrix used to fit in that case is the one for the "regularized" QR facto one. 
      # So its inverse should not involve invIm2QtdQ_Z. 
      # But solve_d2hdv2 is used for gradient -> vecdi2,vecdi3 (and also further for non-gaussian ranefs), 
      # and then we need the exact Hessian. To test this code, compare (SPD with sign) chol-based code to QRP-based one.
      # i.e., switch force_LLF_CHM_QRP.
      delayedAssign("inv_d2hdv2", .calc_inv_d2hdv2_QRP_CHM_signs(BLOB), assign.env = BLOB )
      
      delayedAssign("logdet_sqrt_d2hdv2", {
        logdet_sqrt_RtR <- sum(log(attr(sXaug,"w.ranef")))/2 + BLOB$logdet_R_scaled_v 
        logdet_sqrt_RtR - determinant(BLOB$invIm2QtdQ_Z)$modulus[1]/2
      }, assign.env = BLOB )
      #
      delayedAssign("invIm2QtdQ_Z", { .Rcpp_adhoc_shermanM_sp(BLOB$t_Qq_scaled, c(rep(0L,attr(sXaug,"n_u_h")),BLOB$signs<0)) }, assign.env = BLOB)
      #
      delayedAssign("invIm2QtdQ_ZX", .Rcpp_adhoc_shermanM_sp(BLOB$t_Q_scaled, c(rep(0L,attr(sXaug,"n_u_h")),BLOB$signs<0)), assign.env = BLOB)
      #
      delayedAssign("hatval", {
        colSums(( BLOB$invIm2QtdQ_ZX %*% BLOB$t_Q_scaled) * BLOB$t_Q_scaled)
      } , assign.env = BLOB )
      
    }
    #
    delayedAssign("Z_lev_lambda", .calc_Z_lev_lambda(BLOB), assign.env = BLOB )
    #
    delayedAssign("Z_lev_phi", .calc_Z_lev_phi(BLOB, sXaug), assign.env = BLOB )
    #
    delayedAssign("hatval_Z_u_h_cols_on_left", .calc_hatval_Z_u_h_cols_on_left(BLOB, sXaug), assign.env = BLOB )
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
  if ( ! is.null(szAug)) {
    # if (.is_evaluated("solve_R_scaled", BLOB)) {
    #return((BLOB$solve_R_scaled %*% (BLOB$t_Q_scaled %*% szAug))[BLOB$sortPerm,,drop=FALSE]) 
    #return((BLOB$solve_R_scaled %*% (Matrix::qr.qty(BLOB$QRsXaug,szAug)[seq(ncol(BLOB$solve_R_scaled))]))[BLOB$sortPerm,,drop=FALSE]) 
    # } else if (.is_evaluated("t_Q_scaled", BLOB)) { 
    # even in that case qr.coef seems faster(cf nested_Matern example: poisson -> Pdiag for gradient calculation -> hatval_Z -> t_Q_scaled is computed)
    #return(solve(BLOB$R_scaled, BLOB$t_Q_scaled %*% szAug)[BLOB$sortPerm,,drop=FALSE]) # => one solve R_scaled for t_Q, one in this line  
    # } else {
    if (is.null(BLOB$signs)) { # virtual augsigns=1
      return(qr.coef(BLOB$QRsXaug,szAug)) # avoid t_Q_computation there ## Matrix::qr.coef
    } else{ 
      augsigns <- c(rep(1,attr(sXaug,"n_u_h")), BLOB$signs)
      if (BLOB$nonSPD) { # solves the regularized system because the non-regularized one would as well minimize rather than maximize logL
        return(qr.coef(BLOB$QRsXaug, augsigns*szAug)) # avoid t_Q_computation there ## Matrix::qr.coef
      } else { # Case where use of QRP is forced on (weights with $signs but negHess is SPD)
        # allows to find the correct solution to maximize lik when one forces the use of this sXaug when the negHess is SPD
        return(drop(solve( qrR(BLOB$QRsXaug), BLOB$invIm2QtdQ_ZX %*% (BLOB$t_Q_scaled %*% (augsigns*szAug))))) 
      } 
      #
    }
    # }
  }
  # ELSE 
  if (which=="Qt_leftcols*B") { # only for an y-augmented algo
    if (.is_evaluated("solve_R_scaled", BLOB)) {
      # not so many tests? HLfit(distance ~ age + (age | Subject), data = Orthodont, HLmethod = "REML"): then solve_R_scaled is evaluated
      # corresponding fitme() oes not reach here
      return(.crossprod(BLOB$solve_R_scaled, .crossprod(sXaug[-BLOB$seq_n_u_h,BLOB$perm], B)))
    } else return(solve(t(BLOB$R_scaled), .crossprod(sXaug[-BLOB$seq_n_u_h,BLOB$perm], B))) 
  }
  if (which=="Mg_solve_g") {
    seq_n_u_h <- BLOB$seq_n_u_h
    rhs <- B
    rhs[seq_n_u_h] <- BLOB$invsqrtwranef * rhs[seq_n_u_h]
    # Originally forgot to transpose; fixed 2022/08/05 v3.12.38" and retested between v[3.12.54--3.13.0[.
    # General impact of bug was a much larger value than the true Mg_solve_g, with confusing effects on nloptr tests (eg #362) wrongly suggesting old code was better. 
    # Had to adjust pot4improv by Mg_solve_g_fac to compensate this effect.
    if (.is_evaluated("solve_R_scaled", BLOB)) {
      rhs <- .crossprod(BLOB$solve_R_scaled, rhs[BLOB$perm]) # sum(crossprod(BLOB$solve_R_scaled,B[BLOB$perm])^2) = sum(crossprod(BLOB$solve_R_scaled[BLOB$sortPerm,],B)^2)
    } else rhs <- Matrix::solve(t(BLOB$R_scaled),rhs[BLOB$perm]) 
    if  (is.null(BLOB$signs) || BLOB$nonSPD) { # two cases where the WLS matrix is crossprod(R) 
      return(sum(rhs^2))
    } else return( sum((rhs)* drop(BLOB$invIm2QtdQ_ZX %*% (rhs)))) 
  } 
  if (which=="Mg_invH_g") {
    rhs <- BLOB$invsqrtwranef * B
    if (.is_evaluated("inv_factor_wd2hdv2w",BLOB)) { ## can be assigned elsewhere by lev_lambda <- BLOB$inv_factor_wd2hdv2w <- ....
      rhs <- drop(BLOB$inv_factor_wd2hdv2w %*% rhs)
    } else rhs <- drop(Matrix::solve(BLOB$CHMfactor_wd2hdv2w,rhs,system="L"))
    # No sign => no  invIm2QtdQ_Z factor
    # nonSPD => we use the regularized WLS_mat so do not consider invIm2QtdQ_Z factor here even though it is defined
    # sign SPD => still using QR and exact Hessian here so we use the invIm2QtdQ_Z factor. (vs decorr where we use chol)
    if (is.null(BLOB$signs) || BLOB$nonSPD) { # both cases where the WLS_matrix actually used has no invIm2Q... factor
      return(sum(rhs^2))
    } else return(sum(rhs * drop(BLOB$invIm2QtdQ_Z %*% rhs)))
  } 
  if (which=="Mg_invXtWX_g") { ## for which_LevMar_step="b", not currently used
    seq_n_u_h <- BLOB$seq_n_u_h
    if (is.null(BLOB$XtWX)) BLOB$XtWX <- .crossprod(sXaug[-seq_n_u_h,-seq_n_u_h])
    Mg_invXtWX_g <- crossprod(B,solve(BLOB$XtWX,B))[1L] # [1L] drops possible Matrix class...
    return(Mg_invXtWX_g)
  }
  if (which=="hatval_Z") { ## sscaled Pdiag (or HLfit ML )
    # X[,cols] = Q R P[,cols] = Q q r p => t(Q q) given by:
    if (BLOB$u_h_cols_on_left) { # Rasch... bigranefs... and that is remarkably fast compared to alternative software
      # as well as faster than the alternative code when u_h cols are not on the left. 
      # So, forming inv_factor_wd2hdv2w seems slow compared to this + only using $CHMfactor_wd2hdv2w
      return(BLOB$hatval_Z_u_h_cols_on_left) # no clear benefits in separating lev_phi and lev_lambda computations
    } else if (FALSE) { #  correct but only for pedagogy. 
      if (is.null(BLOB$hatval_Z_)) {
        #  Cf comments above. So forming inv_factor_wd2hdv2w, albeit slow,, is faster than a new QR factorization. 
        tq <- drop0(.lmwithQR(as.matrix(BLOB$R_scaled[,BLOB$sortPerm][,BLOB$seq_n_u_h]) ,yy=NULL,returntQ=TRUE,returnR=FALSE)$t_Q_scaled,
                    tol=1e-16)
        tmp_t_Qq_scaled <- tq %*% BLOB$t_Q_scaled
        xx <- tmp_t_Qq_scaled@x
        xx <- xx*xx
        tmp_t_Qq_scaled@x <- xx
        tmp <- colSums(tmp_t_Qq_scaled)
        n_u_h <- attr(sXaug,"n_u_h")
        phipos <- n_u_h+seq_len(nrow(sXaug)-n_u_h)
        BLOB$hatval_Z_ <-  list(lev_lambda=tmp[-phipos],lev_phi=tmp[phipos]) # BLOB$signs not taken into account in this 'documentation' code
      }
      return(BLOB$hatval_Z_)
    } else {
        # X[,cols] = Q R P[,cols] = Q q r p => t(Q q) given by:
        #t_Qq_scaled <- solve(BLOB$CHMfactor_wd2hdv2w, ## likely bottleneck for large data 
        #                          t(sXaug[, seq_len(attr(sXaug,"n_u_h"))[BLOB$perm_R_v] ]),system="L")
        ## 
        ## t(sXaug[,u_h cols]= (I, scaled t(ZAL)) i.e. is scaled such that the left block is an identity matrix, so we can work 
        ## on two separate blocks if the Cholesky is not permuted. Then 
        hatval_Z_ <- list()
        if ("lambda" %in% B) hatval_Z_$lev_lambda <- BLOB$Z_lev_lambda 
        if ("phi" %in% B) hatval_Z_$lev_phi <- BLOB$Z_lev_phi 
      return(hatval_Z_)
    }
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
          if (.spaMM.data$options$Matrix_old) { # ugly... but such versions do not handle as(, "dMatrix"))
            rhs <- .Dvec_times_m_Matrix(BLOB$invsqrtwranef,B)
          } else rhs <- as(.Dvec_times_m_Matrix(BLOB$invsqrtwranef,B),"generalMatrix")
        } else rhs <- BLOB$invsqrtwranef * B
        # In .calc_sscaled_new, I compute the hatval_Z then solve_d2hdv2. To compute hatval_Z
        # if (u_h_cols_on_left) {
        #   then t_Q_scaled has been computed to compute leverages
        #   neither inv_factor_wd2hdv2w nor CHMfactor_wd2hdv2w are yet computed
        # } else {
        #   CHMfactor_wd2hdv2w and inv_factor_wd2hdv2w has been evaluated to compute leverages.
        # }
        # Thus, when I reach here, I may or may not have inv_factor_wd2hdv2w available
        if ( ! is.null(BLOB$signs)) {
          if (.is_evaluated("inv_factor_wd2hdv2w",BLOB)) {
            rhs <- BLOB$inv_factor_wd2hdv2w %*% rhs
            rhs <- BLOB$invIm2QtdQ_Z %*% rhs
            rhs <- .crossprod(BLOB$inv_factor_wd2hdv2w, rhs) 
          } else {
            rhs <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w,rhs,system="Lt")
            rhs <- BLOB$invIm2QtdQ_Z %*% rhs
            rhs <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w, rhs, system="L")
          }
        } else if ( BLOB$use_R_block || .is_evaluated("inv_factor_wd2hdv2w",BLOB)) {
          rhs <- .crossprod(BLOB$inv_factor_wd2hdv2w, drop(BLOB$inv_factor_wd2hdv2w %*% rhs)) # typical case when solve_d2hdv2 follows hatval_Z in .calc_sscaled_new()
        } else rhs <- Matrix::solve(BLOB$CHMfactor_wd2hdv2w,rhs,system="A") ## dge (if rhs is dense, or a vector), or dgC...
        if (not_vector) { ## is.matrix(rhs) is not the correct test 
          rhs <- .Dvec_times_m_Matrix( - BLOB$invsqrtwranef,rhs)
        } else rhs <- - BLOB$invsqrtwranef * rhs
        return(rhs) # note the minus sign on the vector, - BLOB$invsqrtwranef, rather than the final, possibly matrix, rhs.
      }
    } 
  } 
  if (which=="hatval") {return(BLOB$hatval) }  ## used in get_hatvalues  ##  REML hatval computation (also named hatval_ZX)
  if (which=="R_scaled_blob") { ## used for LevMar
    if (is.null(BLOB$R_scaled_blob)) {
      tmp <- X <- BLOB$R_scaled[ , BLOB$sortPerm] ## crossfac
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
      if (BLOB$use_R_block) { # currently FALSE => CHMfactor_wd2hdv2w is used
        R_scaled_v_h <- t(BLOB$R_scaled[BLOB$sortPerm_u_h,BLOB$sortPerm_u_h, drop=FALSE]) ## the t() for .damping_to_solve... (fixme: if we could avoid t()...)
      } else {
        R_scaled_v_h <- t( as(BLOB$CHMfactor_wd2hdv2w,"sparseMatrix") ) ## the t() for .damping_to_solve... (fixme: if we could avoid t()...)
      }
      tmp <- R_scaled_v_h 
      xx <- tmp@x
      xx <- xx*xx
      tmp@x <- xx
      diag_pRtRp_scaled_v_h <- colSums(tmp)
      BLOB$R_scaled_v_h_blob <- list(R_scaled_v_h=R_scaled_v_h,diag_pRtRp_scaled_v_h=diag_pRtRp_scaled_v_h,
                                     XDtemplate=.XDtemplate(R_scaled_v_h, upperTri=NA))
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
                               XDtemplate=.XDtemplate(R_beta, upperTri=TRUE))
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
    return(.calc_beta_cov_info_from_sXaug(BLOB=BLOB, sXaug=sXaug, tcrossfac=BLOB$solve_R_scaled)) 
  } 
  if (which=="beta_cov_info_from_wAugX") { ## using a weighted Henderson's augmented design matrix, not a true sXaug  
    tPmat <- sparseMatrix(seq_along(BLOB$sortPerm), BLOB$sortPerm, x=1)
    tcrossfac_beta_v_cov <- as.matrix(tPmat %*% BLOB$solve_R_scaled)
    rownames(tcrossfac_beta_v_cov) <- colnames(sXaug) ## necessary for summary.HLfit, already lost in BLOB$R_scaled
    pforpv <- attr(sXaug,"pforpv")
    seqp <- seq_len(pforpv)
    beta_cov <- .tcrossprod(tcrossfac_beta_v_cov[seqp,,drop=FALSE])
    # beta_v_cov <- as.matrix(Matrix::chol2inv(BLOB$R_scaled)[BLOB$sortPerm,BLOB$sortPerm])
    # this tends to be dense bc v_h estimates covary   
    return( list(beta_cov=beta_cov, 
                 #beta_v_cov=beta_v_cov,
                 tcrossfac_beta_v_cov=tcrossfac_beta_v_cov) )
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
# trace("get_from_MME.sXaug_Matrix_QRP_CHM_scaled",exit=quote(str(resu)))
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
