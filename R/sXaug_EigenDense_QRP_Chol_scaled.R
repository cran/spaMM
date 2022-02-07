# 'constructor' for sXaug_Eigen_QR (unpivoted QR factorization) object
# from Xaug which already has a *scaled* ZAL (it has no name in the doc but appears in the eq defining mathcal{W}_w)
def_sXaug_EigenDense_QRP_Chol_scaled <- function(Xaug, # already ZAL_scaled
                                                 weight_X,w.ranef,H_global_scale) {
  n_u_h <- length(w.ranef)
  Xrows <- n_u_h+seq(length(weight_X)) 
  Xaug[Xrows,] <- .Dvec_times_matrix(weight_X, Xaug[Xrows,,drop=FALSE]) ## applying def of mathcal{W}_w in the doc
  attr(Xaug, "get_from") <- "get_from_MME.sXaug_EigenDense_QRP_Chol_scaled"
  attr(Xaug, "BLOB") <- list2env(list(), parent=environment(.sXaug_EigenDense_QRP_Chol_scaled))
  attr(Xaug, "w.ranef") <- w.ranef
  attr(Xaug, "n_u_h") <- n_u_h # mandatory for all sXaug types
  attr(Xaug, "pforpv") <- ncol(Xaug)-n_u_h # mandatory for all sXaug types
  attr(Xaug, "weight_X") <- weight_X # new mandatory 08/2018
  attr(Xaug, "H_global_scale") <- H_global_scale
  class(Xaug) <- c("sXaug_EigenDense_QRP_Chol_scaled", class(Xaug))
  return( Xaug ) 
}

# trace(get_from_MME, print=FALSE, tracer=quote(cat("'",crayon::yellow(which),"'")))
# trace(spaMM:::.sXaug_EigenDense_QRP_Chol_scaled, print=FALSE, tracer=quote(cat("'",which,"'")))
#
.sXaug_EigenDense_QRP_Chol_scaled <- function(sXaug,which="",szAug=NULL,B=NULL) {
  BLOB <- attr(sXaug,"BLOB") ## an environment
  if (is.null(BLOB$R_scaled)) {
    delayedAssign("seq_n_u_h", seq_len(attr(sXaug,"n_u_h")), assign.env = BLOB )
    delayedAssign("invsqrtwranef", 1/sqrt(attr(sXaug,"w.ranef")), assign.env = BLOB )
    delayedAssign("t_Q_scaled", {
      if ( ! is.null(BLOB$perm)) {
        backsolve(BLOB$R_scaled,t(sXaug[,BLOB$perm]),transpose=TRUE)
      } else backsolve(BLOB$R_scaled,t(sXaug),transpose=TRUE)}, assign.env = BLOB )  ## but repetitive usage is minimal?
    delayedAssign("inv_d2hdv2", {
      if (.is_evaluated("inv_factor_wd2hdv2w", BLOB)) { 
        inv_d2hdv2 <- .crossprod(BLOB$inv_factor_wd2hdv2w) # inv_factor_wd2hdv2w is crossfactor
      } else inv_d2hdv2 <- chol2inv(BLOB$R_R_v)
      inv_d2hdv2 <- .m_Matrix_times_Dvec(inv_d2hdv2, BLOB$invsqrtwranef) ## _m_atrix sweep 2L
      inv_d2hdv2 <- - .Dvec_times_matrix(BLOB$invsqrtwranef,inv_d2hdv2)
    }, assign.env = BLOB )
    delayedAssign("R_R_v", { # in "absdiag_R_v", "solve_d2hdv2","hatval_Z","R_scaled_v_h_blob", "logdet_R_scaled_v","Mg_invH_g"
      seq_n_u_h <- BLOB$seq_n_u_h
      # remove $u_h_cols_on_left code in version 2.1.61 since there is a bug in it: further code may require 
      #  no nontrivial perm_R_v, so we remove perm_R_v from further code (as in Matrix_QRP_CHM_scaled version)
      if ( is.null(BLOB$sortPerm)) { # default
        R_R_v <- BLOB$R_scaled[seq_n_u_h,seq_n_u_h, drop=FALSE]
      } else {
        wd2hdv2w <- .crossprodCpp_d(BLOB$R_scaled[,BLOB$sortPerm[ seq_n_u_h ], drop=FALSE], NULL)
        R_R_v <- .Rcpp_chol_R(wd2hdv2w)$R
      }
    }, assign.env = BLOB )
    delayedAssign("absdiag_R_v",if ( is.null(BLOB$sortPerm)) { # default
      abs(.diagfast(x=BLOB$R_scaled)[BLOB$seq_n_u_h])
    } else { # non-default, but to optimize this, perhaps try to get the 
      #  determinant(BLOB$R_scaled[,BLOB$sortPerm[ seq_n_u_h ], drop=FALSE])
      # rather than go through the chol crossprod stuff ? 
      abs(.diagfast(x=BLOB$R_R_v))
    } , assign.env = BLOB )
    delayedAssign("inv_factor_wd2hdv2w", {.Rcpp_backsolve(BLOB$R_R_v, NULL, # NULL i.e., identity, but without creating it (but this does not explain most of the time gain) in Rcpp
                                     transpose=TRUE)}, assign.env = BLOB )
    delayedAssign("logdet_R_scaled_b_v", sum(log(abs(.diagfast(x=BLOB$R_scaled)))), assign.env = BLOB )
    delayedAssign("logdet_R_scaled_v", sum(log(BLOB$absdiag_R_v)), assign.env = BLOB )
    delayedAssign("logdet_r22", {
      # the R's are H-scaled but r22 is H-unscaled... tricky!
      if ( is.null(BLOB$sortPerm)) { # default
        sum(log(abs(.diagfast(x=BLOB$R_scaled)[-BLOB$seq_n_u_h]))) - attr(sXaug,"pforpv")*log(attr(sXaug,"H_global_scale"))/2 # '-' not '<-'
      } else {
        BLOB$logdet_R_scaled_b_v - BLOB$logdet_R_scaled_v - attr(sXaug,"pforpv")*log(attr(sXaug,"H_global_scale"))/2 # '-' not '<-'
      }
    } , assign.env = BLOB )
    delayedAssign("logdet_sqrt_d2hdv2", { sum(log( sqrt(attr(sXaug,"w.ranef"))*BLOB$absdiag_R_v )) } , assign.env = BLOB )
    delayedAssign("hatval", {colSums(BLOB$t_Q_scaled^2)} , assign.env = BLOB ) # also named hatval_ZX
    delayedAssign("Z_lev_lambda", {      
      lev_lambda <- BLOB$inv_factor_wd2hdv2w 
      lev_lambda <- lev_lambda^2
      lev_lambda <- colSums(lev_lambda)
    }, assign.env = BLOB )
    delayedAssign("Z_lev_phi", {      
      n_u_h <- attr(sXaug,"n_u_h")
      phipos <- (n_u_h+1L):nrow(sXaug) #                get the scaled-ZAL block:
      lev_phi <- .tcrossprod(BLOB$inv_factor_wd2hdv2w, sXaug[phipos, BLOB$seq_n_u_h ]) # backsolve(BLOB$R_R_v, t(sXaug[phipos, seq_len(n_u_h) ]),transpose=TRUE) 
      #lev_phi <- .tcrossprod(BLOB$inv_factor_wd2hdv2w, .leftcols.matrix(sXaug,n_u_h)[phipos, ]) # backsolve(BLOB$R_R_v, t(sXaug[phipos, seq_len(n_u_h) ]),transpose=TRUE)
      lev_phi <- lev_phi^2
      lev_phi <- colSums(lev_phi)
    }, assign.env = BLOB )
    ##############################
    EigenDense_QRP_method <- .spaMM.data$options$EigenDense_QRP_method ## pre-09/2018 used .lmwithQRP
    if (EigenDense_QRP_method==".lmwithQR") {     # .lmwithQR() is fast - much faster than qr() and even than chol(crossprod())... 
      lmwithqr <- .lmwithQR(sXaug,yy=szAug,returntQ=FALSE,returnR=TRUE) ## using RcppEigen; szAug may be NULL
      ## we don't request (t) Q from Eigen bc it is comparatively slow 
      # for (st in names(lmwithqr)) BLOB[[st]] <- lmwithqr[[st]] ## "R_scaled" and optionally "coef"
      BLOB$R_scaled <- lmwithqr$R_scaled
      # perm and sortPerm remain NULL
      if ( ! is.null(szAug)) return(lmwithqr$coef)   
      # BLOB$sortPerm will be NULL
    } else if (EigenDense_QRP_method=="qr") {
      ## ( this slightly affects predVar in onelambda vs twolambda)
      QRsXaug <- qr(as(sXaug[],"dgCMatrix")) ## QRsXaug <- qr(.Rcpp_as_dgCMatrix(sXaug))  ## Matrix::qr
      BLOB$perm <- QRsXaug@q + 1L
      BLOB$R_scaled <- as.matrix(qrR(QRsXaug,backPermute = FALSE))
      delayedAssign("sortPerm", sort.list(BLOB$perm), assign.env = BLOB ) # never NULL
      if ( ! is.null(szAug)) return(qr.coef(QRsXaug,szAug))   
    } else {
      lmwithqrp <- .lmwithQRP(sXaug,yy=szAug,returntQ=FALSE,returnR=TRUE) ## using RcppEigen; szAug may be NULL
      ## we don't request (t) Q from Eigen 
      # for (st in names(lmwithqrp)) BLOB[[st]] <- lmwithqrp[[st]] ## "perm", "R_scaled" and optionally "coef"
      BLOB$R_scaled <- lmwithqrp$R_scaled
      BLOB$perm <- BLOB$perm +1L
      delayedAssign("sortPerm", sort.list(BLOB$perm), assign.env = BLOB ) # never NULL
      if ( ! is.null(szAug)) return(lmwithqrp$coef)   
    }
    # any delayedAssign() here would be shadowed by the above return()'s
  } 
  if ( ! is.null(szAug)) {
    rhs <- .crossprodCpp_d(sXaug, szAug)
    if ( ! is.null(BLOB$perm)) rhs <- rhs[BLOB$perm,,drop=FALSE]
    # test-adjacency-corrMatrix (R_scaled of dim 114): replicate(100,{source(...)},simplify=FALSE) finds 5e-3s advantage per replicate for backsolve over .Rcpp_backsolve .
    # I.e., fairly NS; and other tests are not more discriminating
    rhs <- backsolve(BLOB$R_scaled, backsolve(BLOB$R_scaled, rhs, transpose = TRUE))
    if ( is.null(BLOB$sortPerm)) { # default
      return(rhs)
    } else return(rhs[BLOB$sortPerm,,drop=FALSE])
  }
  # ELSE
  if (which=="Qt_leftcols*B") { # use tQ=R^{-T}sX^T in a context where we need only the leftcols of sX^T ie ZX rows of sX
    # if (.is_evaluated("t_Q_scaled",BLOB)) browser() # then we should surely use it, but this does not happen
    if ( ! is.null(BLOB$perm)) {
      return(backsolve(BLOB$R_scaled,.crossprod(sXaug[-BLOB$seq_n_u_h, BLOB$perm], B),transpose=TRUE))
    } else return(backsolve(BLOB$R_scaled,.crossprod(sXaug[-BLOB$seq_n_u_h,], B),transpose=TRUE))
  }
  if (which=="logdet_sqrt_d2hdv2") { return(BLOB$logdet_sqrt_d2hdv2)} 
  if (which=="logdet_r22") { return(BLOB$logdet_r22) }
  if (which=="hatval_Z") { ## Pdiag; note that these leverages are constant wrt to any diagonal rescaling of the Z block
    ## As for the sparse version:
    ## t(sXaug[,u_h cols]= (I, scaled t(ZAL)) i.e. is scaled such that the left block is an identity matrix, so we can work 
    ## on two separate blocks if the Cholesky is not permuted. Then 
    #lev_lambda <- BLOB$inv_factor_wd2hdv2w <- backsolve(BLOB$R_R_v,diag(ncol(BLOB$R_R_v)),transpose=TRUE)
    # test-spaMM # measurable time gain by .Rcpp_backsolve() possibly bc it's the full inverse which is computed. 
    if (FALSE) { # that is correct when .lmwithQR was used but not .lmwithQRP. 
      # Further, I would need an u_h_cols_on_left instead of BLOB$perm for maximum use of this.  
      tmp_t_Qq_scaled <- BLOB$t_Q_scaled[BLOB$seq_n_u_h,]
      tmp_t_Qq_scaled <- tmp_t_Qq_scaled^2
      tmp <- colSums(tmp_t_Qq_scaled)
      phipos <- n_u_h+seq_len(nrow(sXaug)-n_u_h)
      hatval_Z_ <-  list(lev_lambda=tmp[-phipos],lev_phi=tmp[phipos])
    } else {
      hatval_Z_ <- list()
      if ("lambda" %in% B) hatval_Z_$lev_lambda <- BLOB$Z_lev_lambda
      if ("phi" %in% B) hatval_Z_$lev_phi <- BLOB$Z_lev_phi
    }
    return(hatval_Z_)
  } 
  if (which=="solve_d2hdv2") { ## R'R[seq_n_u_h, seq_n_u_h] gives -d2hd_scaled_v2.
    if ( is.null(B) ) {
      return(BLOB$inv_d2hdv2)
    } else { ## solve (matrix,y)
      if (.is_evaluated("inv_d2hdv2", BLOB)) {
        return(BLOB$inv_d2hdv2 %*% B) ## inv_d2hdv2 tends to be dense
      } else {
        if (is.matrix(B)) {
          rhs <- .Dvec_times_matrix(BLOB$invsqrtwranef,B)
        } else rhs <- BLOB$invsqrtwranef * B
        if (.is_evaluated("inv_factor_wd2hdv2w", BLOB)) {
          # test-spaMM Nugget multinomial inverse-Gamma .... no clear benefit in using .Rcpp_chol2solve()
          rhs <- .crossprod(BLOB$inv_factor_wd2hdv2w, drop(BLOB$inv_factor_wd2hdv2w %*% rhs)) # typical sscaled computation
        } else rhs <- backsolve(BLOB$R_R_v, backsolve(BLOB$R_R_v, rhs, transpose = TRUE))
        if (is.matrix(rhs)) {
          rhs <- .Dvec_times_matrix(BLOB$invsqrtwranef,rhs)
        } else rhs <- BLOB$invsqrtwranef * rhs
        return( - rhs)
      }
    }
  } 
  if (which=="R_scaled_blob") {
    if (is.null(BLOB$R_scaled_blob)) {
      if ( ! is.null(BLOB$sortPerm)) {
        X <- BLOB$R_scaled[,BLOB$sortPerm,drop=FALSE] ## has colnames
      } else X <- BLOB$R_scaled
      BLOB$R_scaled_blob <- list(X=X, diag_pRtRp = colSums(X^2), 
                                 XDtemplate=.XDtemplate(X,upperTri=is.null(BLOB$sortPerm)))
    }
    return(BLOB$R_scaled_blob)
  } 
  if (which=="R_scaled_v_h_blob") {
    if (is.null(BLOB$R_scaled_v_h_blob)) {
      diag_pRtRp_scaled_v_h <- colSums(BLOB$R_R_v^2)
      BLOB$R_scaled_v_h_blob <- list(R_scaled_v_h=BLOB$R_R_v, diag_pRtRp_scaled_v_h=diag_pRtRp_scaled_v_h,
                                     XDtemplate=.XDtemplate(BLOB$R_R_v,upperTri=TRUE))
    }
    return(BLOB$R_scaled_v_h_blob)
  } 
  if (which=="R_beta_blob") {
    if (is.null(BLOB$R_beta_blob)) {
      X <- as.matrix(sXaug[-BLOB$seq_n_u_h,-BLOB$seq_n_u_h]) ## following code assuming it is dense...
      R_beta <- .lmwithQR(X,yy=NULL,returntQ=FALSE,returnR=TRUE)$R_scaled
      diag_pRtRp_beta <-  colSums(R_beta^2)
      BLOB$R_beta_blob <- list(R_beta=R_beta,diag_pRtRp_beta=diag_pRtRp_beta, 
                               XDtemplate=.XDtemplate(R_beta, upperTri=TRUE)) 
      ## seems to be tested only by testmv (one of the independent-fit tests on adjacency)
    }
    return(BLOB$R_beta_blob)
  } 
  if (which=="hatval") {return(BLOB$hatval) }  ##  REML hatval computation (also named hatval_ZX)
  if (which=="Mg_invH_g") { # - sum(B %*% inv_d2hdv2 %*% B) # was oddly inconsistent with the Matrix_QRP_CHM code until change in v3.6.4
    rhs <- BLOB$invsqrtwranef * B
    if (.is_evaluated("inv_factor_wd2hdv2w", BLOB)) {
      rhs <- BLOB$inv_factor_wd2hdv2w %*% rhs
    } else rhs <- .Rcpp_backsolve(BLOB$R_R_v,rhs,transpose=TRUE)
    return(sum(rhs^2))
  } 
  if (which=="Mg_solve_g") {
    rhs <- B
    rhs[BLOB$seq_n_u_h] <- BLOB$invsqrtwranef * rhs[BLOB$seq_n_u_h]
    if (is.null(BLOB$perm)) {
      rhs <- backsolve(BLOB$R_scaled, rhs, transpose = TRUE)
    } else rhs <- backsolve(BLOB$R_scaled, rhs[BLOB$perm], transpose = TRUE)
    return(sum(rhs^2))
  } 
  if (which=="Mg_invXtWX_g") { ## 
    if (is.null(BLOB$XtWX)) BLOB$XtWX <- .crossprod(sXaug[-BLOB$seq_n_u_h,-BLOB$seq_n_u_h])
    Mg_invXtWX_g <- crossprod(B,solve(BLOB$XtWX,B))[1L] # [1L] drops possible Matrix class...
    return(Mg_invXtWX_g)
  } 
  if (which=="logdet_R_scaled_b_v") {return(BLOB$logdet_R_scaled_b_v)} 
  if (which=="beta_cov_info_from_sXaug") {  
    return(.calc_beta_cov_info_from_sXaug(BLOB=BLOB, sXaug=sXaug, tcrossfac=solve(BLOB$R_scaled)))
  } 
  if (which=="beta_cov_info_from_wAugX") { ## using a weighted Henderson's augmented design matrix, not a true sXaug  
    if (TRUE) {
      if (ncol(BLOB$R_scaled)) { # excludes GLMs with fully fixed fixef. See comments in .calc_beta_cov_info_others()
        tcrossfac_beta_v_cov <- solve(BLOB$R_scaled) # solve(as(BLOB$R_scaled,"dtCMatrix"))
        if ( ! is.null(BLOB$sortPerm)) { # depending on method used for QR facto
          tPmat <- sparseMatrix(seq_along(BLOB$sortPerm), BLOB$sortPerm, x=1)
          tcrossfac_beta_v_cov <- as.matrix(tPmat %*% tcrossfac_beta_v_cov)
        } else tcrossfac_beta_v_cov <- as.matrix(tcrossfac_beta_v_cov)
        rownames(tcrossfac_beta_v_cov) <- colnames(sXaug) ## necessary for summary.HLfit, already lost in BLOB$R_scaled
        # extract beta cols:
        seqp <- seq_len( attr(sXaug,"pforpv"))
        beta_cov <- .tcrossprod(tcrossfac_beta_v_cov[seqp,,drop=FALSE])
        return(list(beta_cov=beta_cov, 
                    #beta_v_cov=beta_v_cov,
                    tcrossfac_beta_v_cov=tcrossfac_beta_v_cov)) 
      } else return(list(beta_cov=matrix(ncol=0,nrow=0), 
                         tcrossfac_beta_v_cov=matrix(ncol=0,nrow=0))) 
    } else {
      beta_v_cov <- chol2inv(BLOB$R_scaled)
      if ( ! is.null(BLOB$sortPerm)) beta_v_cov <- beta_v_cov[BLOB$sortPerm,BLOB$sortPerm,drop=FALSE]
      # this tends to be dense bc v_h estimates covary (example: wafers)
      # otherwise, dropO(,tol=...), + some fix in summary.HLfit for matrix[] <- Matrix assignment, would be useful.  
      return(beta_v_cov)
    }
    ########################
  } else if (which %in% c("absdiag_R_v")) { 
    return(BLOB$absdiag_R_v)
  # } else if (which=="sortPerm") { 
  #   return(BLOB$sortPerm)
  }
  if (which=="t_Q_scaled") { # residual call for non-standard REML 
    return(BLOB$t_Q_scaled)
  }
  if (which=="d2hdv2") {
    stop("d2hdv2 requested")
    # don't forget that the factored matrix is not the augmented design matrix ! hence w.ranef needed here
    w.ranef <- attr(sXaug,"w.ranef")
    w_R_R_v <- .Matrix_times_Dvec(BLOB$R_R_v, sqrt(w.ranef))  #BLOB$R_R_v %*% diag(sqrt(w.ranef))
    BLOB$d2hdv2 <- - .crossprodCpp_d(w_R_R_v,yy=NULL)
    return(BLOB$d2hdv2)
  } 
  stop("invalid 'which' value.")
} 

# trace("get_from_MME.sXaug_EigenDense_QRP_Chol_scaled",print=FALSE, tracer=quote(print(which)),exit=quote(str(resu)))
get_from_MME.sXaug_EigenDense_QRP_Chol_scaled <- function(sXaug,which="",szAug=NULL,B=NULL,
                                    damping, LMrhs, ...) {
  resu <- switch(which,
                 "LevMar_step" = {
                   R_scaled_blob <- .sXaug_EigenDense_QRP_Chol_scaled(sXaug,which="R_scaled_blob")
                   dampDpD <- damping*R_scaled_blob$diag_pRtRp ## NocedalW p. 266
                   list(dVscaled_beta = .damping_to_solve(XDtemplate=R_scaled_blob$XDtemplate, dampDpD=dampDpD,rhs=LMrhs), 
                        dampDpD = dampDpD) 
                 },
                 "LevMar_step_v_h" = {
                   ## FR->FR probably not the most elegant implementation 
                   R_scaled_v_h_blob <- .sXaug_EigenDense_QRP_Chol_scaled(sXaug,which="R_scaled_v_h_blob")
                   dampDpD <- damping*R_scaled_v_h_blob$diag_pRtRp_scaled_v_h ## NocedalW p. 266
                   list(dVscaled = .damping_to_solve(XDtemplate=R_scaled_v_h_blob$XDtemplate, dampDpD=dampDpD, rhs=LMrhs), 
                        dampDpD = dampDpD) 
                 },
                 "LevMar_step_beta" = {
                   if ( ! length(LMrhs)) stop("LevMar_step_beta called with 0-length LMrhs: pforpv=0?")
                   R_beta_blob <- .sXaug_EigenDense_QRP_Chol_scaled(sXaug,which="R_beta_blob")
                   dampDpD <- damping*R_beta_blob$diag_pRtRp_beta
                   list(dbeta = .damping_to_solve(XDtemplate=R_beta_blob$XDtemplate, dampDpD=dampDpD, rhs=LMrhs), 
                        dampDpD = dampDpD) 
                 },
                 ## all other cases:
                 .sXaug_EigenDense_QRP_Chol_scaled(sXaug,which=which,szAug=szAug,B=B)
  )
  return(resu)
}

