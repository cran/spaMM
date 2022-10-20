# 'constructor' for sXaug_Eigen_QR (unpivoted QR factorization) object
# from Xaug which already has a *scaled* ZAL (it has no name in the doc but appears in the eq defining mathcal{W}_w)
#
# This properly handles Hobs negative weights: track all code dependent on BLOB$nonSPD and BLOB$signs.
#
def_sXaug_EigenDense_QRP_Chol_scaled <- function(Xaug, # already ZAL_scaled
                                                 weight_X,w.ranef,H_global_scale,
                                                 force_QRP=NULL # ignored
                                                 ) {
  n_u_h <- length(w.ranef)
  Xrows <- n_u_h+seq(length(weight_X)) 
  Xaug[Xrows,] <- .Dvec_times_matrix(weight_X, Xaug[Xrows,,drop=FALSE]) ## applying def of mathcal{W}_w in the doc
  attr(Xaug, "get_from") <- "get_from_MME.sXaug_EigenDense_QRP_Chol_scaled"
  H_w.resid <- attr(weight_X,"H_w.resid")
  BLOB <- list2env(list(H_w.resid=H_w.resid,
               signs=attr(H_w.resid,"signs")), 
          parent=emptyenv())
  if ( BLOB$nonSPD <- ! is.null(BLOB$signs)) {
    augsigns <- c(rep(1,n_u_h), BLOB$signs)
    BLOB$signed <- .Dvec_times_matrix(augsigns,Xaug) 
    negHess <- crossprod(BLOB$signed,Xaug)
    BLOB$R_scaled <- try(chol(negHess), silent=TRUE)
    if (BLOB$nonSPD <- inherits(BLOB$R_scaled, "try-error")) {
      BLOB$R_scaled <- BLOB$signed <- NULL
    }
  } 
  if (BLOB$nonSPD) {
    BLOB$WLS_mat_weights <- abs(H_w.resid)
  } else BLOB$WLS_mat_weights <- H_w.resid
  attr(Xaug, "BLOB") <- BLOB
  attr(Xaug, "w.ranef") <- w.ranef
  attr(Xaug, "n_u_h") <- n_u_h # mandatory for all sXaug types
  attr(Xaug, "pforpv") <- ncol(Xaug)-n_u_h # mandatory for all sXaug types
  attr(Xaug, "weight_X") <- weight_X # new mandatory 08/2018
  attr(Xaug, "H_global_scale") <- H_global_scale
  .init_promises_decorr(sXaug=Xaug)
  class(Xaug) <- c("sXaug_EigenDense_QRP_Chol_scaled", class(Xaug))
  return( Xaug ) 
}

.init_promises_decorr <- function(sXaug) {
  BLOB <- attr(sXaug,"BLOB") 
  delayedAssign("seq_n_u_h", seq_len(attr(sXaug,"n_u_h")), assign.env = BLOB )
  delayedAssign("invsqrtwranef", 1/sqrt(attr(sXaug,"w.ranef")), assign.env = BLOB )
  delayedAssign("t_Q_scaled", {
    if ( ! is.null(BLOB$perm)) {
      backsolve(BLOB$R_scaled,t(sXaug[,BLOB$perm]),transpose=TRUE)
    } else backsolve(BLOB$R_scaled,t(sXaug),transpose=TRUE)}, assign.env = BLOB )  ## but repetitive usage is minimal?
  #
  delayedAssign("t_Qq_scaled", {
    n_u_h <- attr(sXaug,"n_u_h")
    phipos <- (n_u_h+1L):nrow(sXaug)
    .tcrossprod(BLOB$inv_factor_wd2hdv2w, sXaug[, BLOB$seq_n_u_h ])
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
  delayedAssign("absdiag_R_v", {
    if ( is.null(BLOB$sortPerm)) { # default
      abs(.diagfast(x=BLOB$R_scaled)[BLOB$seq_n_u_h])
    } else { # non-default, but to optimize this, perhaps try to get the 
      #  determinant(BLOB$R_scaled[,BLOB$sortPerm[ seq_n_u_h ], drop=FALSE])
      # rather than go through the chol crossprod stuff ? 
      abs(.diagfast(x=BLOB$R_R_v))
    }
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
  
  if (BLOB$nonSPD) { # then QR and $signs should be used
    delayedAssign("inv_d2hdv2", .calc_inv_d2hdv2_decorr_nonSPD(BLOB), assign.env = BLOB )
    delayedAssign("invIm2QtdQ_ZX", .shermanM_de(Qt=BLOB$t_Q_scaled, indic = c(rep(0,attr(sXaug,"n_u_h")),BLOB$signs<0) ), assign.env = BLOB)
    delayedAssign("invIm2QtdQ_Z", .shermanM_de(Qt=BLOB$t_Qq_scaled, indic = c(rep(0,attr(sXaug,"n_u_h")),BLOB$signs<0)), assign.env = BLOB)
    # redefines:
    delayedAssign("hatval", { 
      tmp <- colSums(( BLOB$invIm2QtdQ_ZX %*% BLOB$t_Q_scaled) * BLOB$t_Q_scaled) 
      # the lev_lam are signed too so we must be careful not changing their sign
      tmp[-BLOB$seq_n_u_h] <- tmp[-BLOB$seq_n_u_h]*BLOB$signs
      tmp
    } , assign.env = BLOB ) 
    delayedAssign("Z_lev_lambda", { 
      t_Qq_lam_cols <- BLOB$inv_factor_wd2hdv2w #  think of R^-T %*% I-block-of-sXaug...
      colSums(( BLOB$invIm2QtdQ_Z %*% t_Qq_lam_cols) * t_Qq_lam_cols)
    } , assign.env = BLOB )
    #
    delayedAssign("Z_lev_phi", { 
      n_u_h <- attr(sXaug,"n_u_h")
      phipos <- (n_u_h+1L):nrow(sXaug)                 
      t_Qq_phi_cols <- BLOB$t_Qq_scaled[,phipos] 
      lev_phi <- colSums(( BLOB$invIm2QtdQ_Z %*% t_Qq_phi_cols) * t_Qq_phi_cols)*BLOB$signs
    } , assign.env = BLOB )
    #
    ## Speculative distinct code for nonSPD case; but the H used to fit in that case is the one for the "regularized" QR facto one. 
    # So its inverse should not involve invIm2QtdQ_Z. Only if inv_d2hdv2 were used for gradient -> hat values, or likelihood calculations,
    # we would need to correct from inv_factor_wd2hdv2w to recover correct values of the latter variables.
    # delayedAssign("inv_d2hdv2", {
    #   rhs <- .Matrix_times_Dvec(BLOB$inv_factor_wd2hdv2w, BLOB$invsqrtwranef)
    #   inv_d2hdv2 <- - .crossprod(rhs, BLOB$invIm2QtdQ_Z %*% rhs)
    # }, assign.env = BLOB )
    #
    delayedAssign("logdet_sqrt_d2hdv2", { 
      logdet_sqrt_RtR <- sum(log( sqrt(attr(sXaug,"w.ranef"))*BLOB$absdiag_R_v )) 
      logdet_sqrt_RtR - determinant(BLOB$invIm2QtdQ_Z)$modulus[1]/2
    } , assign.env = BLOB )
  } else {
    delayedAssign("inv_d2hdv2", .calc_inv_d2hdv2_decorr(BLOB), assign.env = BLOB )
    delayedAssign("logdet_sqrt_d2hdv2", { 
      sum(log( sqrt(attr(sXaug,"w.ranef"))*BLOB$absdiag_R_v )) 
    } , assign.env = BLOB )
    delayedAssign("hatval", {
      hatval <- colSums(BLOB$t_Q_scaled^2)
      if  (! is.null(BLOB$signs)) {
        n_u_h <- attr(sXaug,"n_u_h")
        phipos <- n_u_h+seq_len(nrow(sXaug)-n_u_h)
        hatval[phipos] <- hatval[phipos]*BLOB$signs  # correct only for chol (SPD) case
      } 
      hatval
    } , assign.env = BLOB ) # also named hatval_ZX
    delayedAssign("Z_lev_lambda", {      
      lev_lambda <- BLOB$inv_factor_wd2hdv2w 
      lev_lambda <- lev_lambda^2
      lev_lambda <- colSums(lev_lambda)
    }, assign.env = BLOB )
    delayedAssign("Z_lev_phi", {      
      # n_u_h <- attr(sXaug,"n_u_h")
      # phipos <- (n_u_h+1L):nrow(sXaug) #                get the scaled-ZAL block:
      # lev_phi <- .tcrossprod(BLOB$inv_factor_wd2hdv2w, sXaug[phipos, BLOB$seq_n_u_h ]) # backsolve(BLOB$R_R_v, t(sXaug[phipos, seq_len(n_u_h) ]),transpose=TRUE) 
      # lev_phi <- .tcrossprod(BLOB$inv_factor_wd2hdv2w, .leftcols.matrix(sXaug,n_u_h)[phipos, ]) # backsolve(BLOB$R_R_v, t(sXaug[phipos, seq_len(n_u_h) ]),transpose=TRUE)
      n_u_h <- attr(sXaug,"n_u_h")
      phipos <- (n_u_h+1L):nrow(sXaug)                 
      t_Qq_phi_cols <- BLOB$t_Qq_scaled[,phipos] 
      lev_phi <- t_Qq_phi_cols^2
      lev_phi <- colSums(lev_phi)
      if  (! is.null(BLOB$signs)) lev_phi <- lev_phi*BLOB$signs # correct only for chol (SPD) case
      lev_phi
    }, assign.env = BLOB )
    
  }
}

# trace(get_from_MME, print=FALSE, tracer=quote(cat("'",crayon::yellow(which),"'")))
# trace(spaMM:::.sXaug_EigenDense_QRP_Chol_scaled, print=FALSE, tracer=quote(cat("'",which,"'")))
#
.sXaug_EigenDense_QRP_Chol_scaled <- function(sXaug,which="",szAug=NULL,B=NULL) { 
  BLOB <- attr(sXaug,"BLOB") 
  
  ###
  
  if (is.null(BLOB$R_scaled)) { # may occur if no signs, or signs and (failed chol <=> non SPD) 
    if ( ! is.null(BLOB$signs)) { # => failed Cholesky => use QR
      # cf explanations of the use of signs in the post-factorization code, below.
      if (BLOB$nonSPD) { # => this must be TRUE
        if ( ! is.null(szAug)) {
          augsigns <- c(rep(1,attr(sXaug,"n_u_h")), BLOB$signs)
          lmwithqr <- .lmwithQR(sXaug,yy=augsigns*szAug,returntQ=FALSE,returnR=TRUE) ## using RcppEigen; szAug may be NULL (but avoid numeic(0) !)
        } else lmwithqr <- .lmwithQR(sXaug,yy=NULL,returntQ=FALSE,returnR=TRUE) ## szAug may be NULL (but avoid augsigns*NULL =numeric(0) !)
        BLOB$R_scaled <- lmwithqr$R_scaled
        # perm and sortPerm remain NULL
        if ( ! is.null(szAug)) return(lmwithqrp$coef)   # return solution of regularized system in non-SPD case
      } else { # there are $signs but negHess was SPD => $R_scaded must already be present => this alternative never occurs.
        # we could imagine forcing QR, for devel purposes, here, as for spcorr case. But this is already complicated enough
        # augsigns <- c(rep(1,attr(sXaug,"n_u_h")), BLOB$signs)
        # coefs <- backsolve( BLOB$R_scaled, BLOB$invIm2QtdQ_ZX %*% (BLOB$t_Q_scaled %*% (augsigns*szAug)))
        if ( ! is.null(szAug)) {
          rhs <- crossprod(BLOB$signed,  szAug) # H beta + (grad = X ; dlogcLdeta)
          coef <- backsolve(BLOB$R_scaled, backsolve(BLOB$R_scaled, rhs, transpose = TRUE)) # solve(mH, mH beta + grad) = beta + mH grad = beta - D2logLDbeta2 . DlogLDbeta
          return(drop(coef))   # return for $signs, but Hessian still SPD (=> chol facto)
        }
      }
    } else { # NULL $signs
      EigenDense_QRP_method <- .spaMM.data$options$EigenDense_QRP_method ## pre-09/2018 used .lmwithQRP
      if (EigenDense_QRP_method==".lmwithQR") {     # DEFAULT .lmwithQR() is fast - much faster than qr() and even than chol(crossprod())... 
        lmwithqr <- .lmwithQR(sXaug,yy=szAug,returntQ=FALSE,returnR=TRUE) ## using RcppEigen; szAug may be NULL
        ## we don't request (t) Q from Eigen bc it is comparatively slow 
        # for (st in names(lmwithqr)) BLOB[[st]] <- lmwithqr[[st]] ## "R_scaled" and optionally "coef"
        BLOB$R_scaled <- lmwithqr$R_scaled
        # perm and sortPerm remain NULL
        if ( ! is.null(szAug)) return(lmwithqr$coef)   # return for NO $signs
        # BLOB$sortPerm will be NULL
      } else if (EigenDense_QRP_method=="qr") {
        ## ( this slightly affects predVar in onelambda vs twolambda)
        if (.spaMM.data$options$Matrix_old) { # this block appears to evade the long tests
          QRsXaug <- qr(as(sXaug[],"dgCMatrix")) ## QRsXaug <- qr(.Rcpp_as_dgCMatrix(sXaug))  ## Matrix::qr
        } else  QRsXaug <- qr(as(as(sXaug[],"generalMatrix"),"CsparseMatrix")) ## QRsXaug <- qr(.Rcpp_as_dgCMatrix(sXaug))  ## Matrix::qr
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
    } 
  } 
  
  ###
  
  if ( ! is.null(szAug)) { # then reach within this block when $R_scaled was already available, so the coefs were not provided by the previous block
    # Not using inv(R) %*% tQ %*% szAug here... presumably bc tQ not computed by default, spec. with default .lmwithQR method.
    if (is.null(BLOB$signs)) {
      rhs <- .crossprodCpp_d(sXaug, szAug)
      if ( ! is.null(BLOB$perm)) rhs <- rhs[BLOB$perm,,drop=FALSE]
      # test-adjacency-corrMatrix (R_scaled of dim 114): replicate(100,{source(...)},simplify=FALSE) finds 5e-3s advantage per replicate for backsolve over .Rcpp_backsolve .
      # I.e., fairly NS; and other tests are not more discriminating
      rhs <- backsolve(BLOB$R_scaled, backsolve(BLOB$R_scaled, rhs, transpose = TRUE))
      if ( is.null(BLOB$sortPerm)) { # default
        return(rhs)
      } else return(rhs[BLOB$sortPerm,,drop=FALSE])
    } else {
      # no BLOB$perm... are needed here bc unpermuted factos are used, never the experimental permuted ones
      if (BLOB$nonSPD) { # We then solve the regularized system 
        augsigns <- c(rep(1,attr(sXaug,"n_u_h")), BLOB$signs) # there is no BLOB$signed
        rhs <- .crossprodCpp_d(sXaug, augsigns*szAug) # as in spcorr case 
        return(backsolve(BLOB$R_scaled, backsolve(BLOB$R_scaled, rhs, transpose = TRUE)))
      } else { # Typically does not happen: in SPD signs case I used chol rather than QR, so there is no BLOB$invIm2QtdQ_ZX here.
        # I could imagine forcing QR in signed SPD case as in spcorr methods, but not now. Then I would need ~
        #   return(backsolve(BLOB$R_scaled, BLOB$invIm2QtdQ_ZX %*% backsolve(BLOB$R_scaled, rhs, transpose = TRUE)))
        rhs <- .crossprodCpp_d(BLOB$signed,  szAug) # H beta + (grad = X ; dlogcLdeta), correct for chol case.
        # which should be equivalent to .crossprodCpp_d(sXaug, augsigns*szAug)
        return(backsolve(BLOB$R_scaled, backsolve(BLOB$R_scaled, rhs, transpose = TRUE)))
        # Again if I used QR here I would need an BLOB$invIm2QtdQ_ZX matrix, but again chol() is usd  
      }
    }
  }
  
  ### ELSE
  
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
    if (FALSE) { # that is correct when .lmwithQR was used but not .lmwithQRP. (nor for signed weights)
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
        if (BLOB$nonSPD) {
          rhs <- .crossprod(BLOB$inv_factor_wd2hdv2w, drop(BLOB$invIm2QtdQ_Z %*% (BLOB$inv_factor_wd2hdv2w %*% rhs)))
        } else if (.is_evaluated("inv_factor_wd2hdv2w", BLOB)) {
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
  if (which=="Mg_invH_g") { # - sum(B %*% inv_d2hdv2 %*% B)
    rhs <- BLOB$invsqrtwranef * B
    if (.is_evaluated("inv_factor_wd2hdv2w", BLOB)) {
      rhs <- BLOB$inv_factor_wd2hdv2w %*% rhs
    } else rhs <- .Rcpp_backsolve(BLOB$R_R_v,rhs,transpose=TRUE)
    # No sign => no  invIm2QtdQ_Z factor
    # nonSPD => we use the regularized WLS_mat so do not consider invIm2QtdQ_Z factor here even though it is defined
    # sign SPD => we currently use a Chol facto here so no invIm2QtdQ_Z factor.
    #if (is.null(BLOB$signs) || BLOB$nonSPD) {
      return(sum(rhs^2))
    #} else return(sum(rhs * drop(BLOB$invIm2QtdQ_Z %*% rhs))) 
  } 
  if (which=="Mg_solve_g") {
    rhs <- B
    rhs[BLOB$seq_n_u_h] <- BLOB$invsqrtwranef * rhs[BLOB$seq_n_u_h]
    if (is.null(BLOB$perm)) {
      rhs <- backsolve(BLOB$R_scaled, rhs, transpose = TRUE)
    } else rhs <- backsolve(BLOB$R_scaled, rhs[BLOB$perm], transpose = TRUE)
    ## Same comment as on "inv_d2hdv2"...
    return(sum(rhs^2)) # correct even in the signed SPD case (no invIm2Q... correction is implemented bc "EigenDense_QRP" actually uses chol in that case)
  } 
  if (which=="Mg_invXtWX_g") { ## 
    if (is.null(BLOB$XtWX)) BLOB$XtWX <- .crossprod(sXaug[-BLOB$seq_n_u_h,-BLOB$seq_n_u_h])
    Mg_invXtWX_g <- crossprod(B,solve(BLOB$XtWX,B))[1L] # [1L] drops possible Matrix class...
    return(Mg_invXtWX_g)
  } 
  if (which=="logdet_R_scaled_b_v") {return(BLOB$logdet_R_scaled_b_v)} 
  if (which=="beta_cov_info_from_sXaug") {  
    # backsolve() more robust numerically than solve(). Avoids stop on an otherwise quite poor fit.
    return(.calc_beta_cov_info_from_sXaug(BLOB=BLOB, sXaug=sXaug, tcrossfac=backsolve(BLOB$R_scaled, diag(ncol(sXaug))))) 
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

.calc_inv_d2hdv2_decorr <- function(BLOB) {
  if (.is_evaluated("inv_factor_wd2hdv2w", BLOB)) { 
    inv <- .crossprod(BLOB$inv_factor_wd2hdv2w) # inv_factor_wd2hdv2w is crossfactor
  } else inv <- chol2inv(BLOB$R_R_v)
  # from inv_wd2hdv2w to inv_d2hdv2:
  inv <- .m_Matrix_times_Dvec(inv, BLOB$invsqrtwranef) ## _m_atrix sweep 2L
  inv <- - .Dvec_times_matrix(BLOB$invsqrtwranef,inv)
  inv
}

.calc_inv_d2hdv2_decorr_nonSPD <- function(BLOB) {
  inv <- .crossprod(BLOB$inv_factor_wd2hdv2w, BLOB$invIm2QtdQ_Z %*% BLOB$inv_factor_wd2hdv2w) 
  # from inv_wd2hdv2w to inv_d2hdv2:
  inv <- .m_Matrix_times_Dvec(inv, BLOB$invsqrtwranef) ## _m_atrix sweep 2L
  inv <- - .Dvec_times_matrix(BLOB$invsqrtwranef,inv)
  inv
}
