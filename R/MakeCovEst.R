.upper_tri_tcrossfactorize <- function( esys, d_regul=esys$d_regul, compactcovmat) {
  # SPPREC or use triangular facto for NON-spprec   # currently with d=1
  ## Competing for clumsy code prize... 
  ## need triangular factor for spprec case in particular hence the default value of use_tri_Nspprec
  crossfac_prec <- .Dvec_times_matrix(sqrt(1/d_regul), t(esys$vectors))
  qrblob <- qr(crossfac_prec)
  crossfac_prec <- qr.R(qrblob) # applying .lmwithQR() systematically (badly) affects numerical precision
  if (! all(unique(diff(qrblob$pivot))==1L)) { # eval an unpermuted triangular R
    crossfac_prec <- .lmwithQR(crossfac_prec[, sort.list(qrblob$pivot)] ,yy=NULL,returntQ=FALSE,returnR=TRUE)$R
  } 
  tcrossfac_prec <- t(crossfac_prec) ## tcrossprod(tcrossfac_prec)=compactprecmat
  # Correct signs:
  diagsigns <- sign(.diagfast(tcrossfac_prec))
  tcrossfac_prec <- .m_Matrix_times_Dvec(tcrossfac_prec,diagsigns)
  ## now with triangular factor: .gmp_solve() spares alternative code, but regularization remains necessary anyway.
  # invdcp <- 1/.diagfast(x=tcrossfac_prec)
  # safe_solvand <- .m_Matrix_times_Dvec(tcrossfac_prec, invdcp)
  # list(design_u=t(solve(safe_solvand)) %*% diag(x=invdcp),
  list(design_u=t(.gmp_solve(tcrossfac_prec)), # tcrossprod factor, not orthonormed, UPPER tri (while eigen provides an orthonormed "full" U matrix)
       # radical removal of d : v3.9.19. See notably comments in .post_process_v_h_LMatrices() for good reasons not to reintroduce them.
       compactcovmat=compactcovmat, ## not used for the fit
       compactchol_Q=tcrossfac_prec ## for sparse precision algo. ## compactprecmat=.ZWZtwrapper(chol_Q,1/invdcp^2)
  ) 
}

.lower_tri_tcrossfactorize <- function( esys, d_regul=esys$d_regul, compactcovmat) { ## # LOWER tri tcrossprod factor. 
  # Simple code and this may well be a correct alternative to .upper_tri_tcrossfactorize() as rescue code:
  # (1) checks for v3.8.39 suggests that this works as a general replacement for chol() (although slower).
  # (2) it passes the rC_transf check.
  #
  # qr always produces a crossprod factor so if we want a tcrossprod factor, it must be t(qr.R()) hence lower tri...); unless we qr() the preci mat...
  qrand <- t(.m_Matrix_times_Dvec(esys$vectors,sqrt(d_regul)))
  qrblob <- qr(qrand)
  crossfac <- qr.R(qrblob) # applying .lmwithQR() systematically (badly) affects numerical precision
  if (! all(unique(diff(qrblob$pivot))==1L)) { # eval an unpermuted triangular R
    crossfac <- .lmwithQR(crossfac[, sort.list(qrblob$pivot)] ,yy=NULL,returntQ=FALSE,returnR=TRUE)$R
  } 
  tcrossfac <- t(crossfac) # lower.tri
  list(design_u=tcrossfac, # LOWER tri tcrossprod factor
       # radical removal of d : v3.9.19. See notably comments in .post_process_v_h_LMatrices() for good reasons not to reintroduce them.
       compactcovmat=compactcovmat ## not used for the fit
  ) ## and one might use (solve(design_u)) as a $crossfac_Q precision factor
}

.calc_latentL <- function(compactcovmat, use_tri_CORREL=TRUE, spprecBool, trDiag) {
  # returns a *t*crossfactor 
  ## L and L_Q are distinct matrices used jointly in a fit (L in ZAL...) and must be deduced from each other.
  if (spprecBool) {
    if (regularize <- TRUE) {
      # always regularize first
      compactcovmat <- .smooth_regul(compactcovmat, epsi=.spaMM.data$options$tol_ranCoefs_inner["regul"]) 
      #   value of epsi ultimately also affects whether bobyqa has to be run 
      esys <- attr(compactcovmat,"esys") 
      blob <- .upper_tri_tcrossfactorize( esys, d_regul=esys$d_regul, compactcovmat=compactcovmat)
    } else {
      esys <- .eigen_sym(compactcovmat)
      blob <- .upper_tri_tcrossfactorize( esys, d_regul=esys$values, compactcovmat=compactcovmat)
    }
    ## for sparse precision we want chol_Q to be (dtCMatrix: Csparse triangular) so that efficient solve methods can be used.
    ## (1) direct conversion from matrix to dtC is remarkably slower than going through the as(.,"sparseMatrix") step !
    ## (2) Even with two steps as(.,"dtCMatrix") is costly hence the spprecBool condition (dtCMatrix useful only for chol_Q)
    blob$compactchol_Q <- as(as(blob$compactchol_Q,"sparseMatrix"),"triangularMatrix") # precision factor for sparse_precision algo; # compactcovmat=.ZtWZwrapper(solve(compactchol_Q),d  )
    blob$trDiag <- trDiag
    #  We also want chol_Q it to be lower triangular (dtCMatrix can be Up or Lo) 
    #    to obtain a dtCMatrix by bdiag(list of lower tri dtCMatrices). It already is lower tri.
    # Beyond this, design_u may be lower or upper tri
    if (.spaMM.data$options$Matrix_old) { # ugly... but such versions do not handle as(, "dMatrix"))
      evec <- as(esys$vectors, "dgCMatrix") 
    } else evec <- as(as(esys$vectors,"generalMatrix"),"CsparseMatrix") # the aim of converison to dgC is that ZWZt is automatically dsC through call to Matrix::tcrossprod
    if (regularize) {
      blob$compactprecmat <- .ZWZtwrapper(evec, 1/esys$d_regul) # dsCMatrix
    } else blob$compactprecmat <- .ZWZtwrapper(evec, 1/esys$values) # dsCMatrix
    # if (regularize) {
    #   compactprecmat <- .ZWZtwrapper(evec, 1/esys$d_regul) 
    # } else compactprecmat <- .ZWZtwrapper(evec, 1/esys$values) 
    # blob$compactprecmat <- as(forceSymmetric(compactprecmat),"CsparseMatrix")
  } else { # CORREL algos
    if (is.null(chol_crossfac <- attr(compactcovmat, "chol_crossfac")))  chol_crossfac <- tryCatch(chol(compactcovmat),error=function(e) e)
    if ( ! inherits(chol_crossfac, "simpleError")) {
      blob <- list(design_u=t(chol_crossfac), 
                   # radical removal of d : v3.9.19
                   compactcovmat=compactcovmat)
    } else { # RESCUE CODE
      if (regularize <- TRUE) { # necess bc eigensystem is not accurate enough (it may have slightly negative eigenvalues)
        # regularize the matrix, then choose between different representations of it.
        compactcovmat <- .smooth_regul(compactcovmat, epsi=.spaMM.data$options$tol_ranCoefs_inner["regul"]) 
        #   value of epsi utlmately also affects whether bobyqa has to be run 
        esys <- attr(compactcovmat,"esys")
        if (use_tri_CORREL) {  
          blob <- .upper_tri_tcrossfactorize( esys, d_regul=esys$d_regul, compactcovmat=compactcovmat)
        } else { ## # LOWER tri tcrossprod factor. 
          blob <- .lower_tri_tcrossfactorize( esys, d_regul=esys$d_regul, compactcovmat=compactcovmat)
        } 
      } else {
        esys <- .eigen_sym(compactcovmat)
        if (use_tri_CORREL) {  
          blob <- .upper_tri_tcrossfactorize( esys, d_regul=esys$values, compactcovmat=compactcovmat)
        } else { ## # LOWER tri tcrossprod factor. 
          blob <- .lower_tri_tcrossfactorize( esys, d_regul=esys$values, compactcovmat=compactcovmat)
        } 
      }
      #
    }
  }
  
  blob
}

.gmp.kronecker2D <- function(X,Y) { # X, Y are 2D arrays of gmp type or not
  # cf base::.kronecker, but we perform operations on seq(length(p)) 
  # because they cannot be performed directly on the gmp:: object
  p <- outer(X,Y)
  i <- seq(length(p))
  dim(i) <- c(dim(X),dim(Y))
  dp <- as.vector(t(matrix(1L:4L, ncol = 2L)[, 2:1])) 
  pi <- aperm(i,dp)
  dim(pi) <- dim(X)*dim(Y)
  pp <- p[[pi]]
  dim(pp) <- dim(pi)
  pp
}

# kronecker(matrix(seq(6),nrow=2,ncol=3),matrix(seq(20),nrow=4,ncol=5))
# .gmp.kronecker2D(as.bigq(matrix(seq(6),nrow=2,ncol=3)),as.bigq(matrix(seq(20),nrow=4,ncol=5)))
# .gmp.kronecker2D(as.bigq(matrix(seq(6),nrow=2,ncol=3)),matrix(seq(20),nrow=4,ncol=5))
# .gmp.kronecker2D(matrix(seq(6),nrow=2,ncol=3), as.bigq(matrix(seq(20),nrow=4,ncol=5)))

.makelong_bigq <- function(Lcompact, longsize, kron_Y) { ## ad hoc version with restricted usage
  if ( ! is.null(kron_Y)) return(.gmp.kronecker2D(Lcompact,kron_Y)) 
  # ELSE
  n_levels <- longsize/ncol(Lcompact)
  longLv <- gmp::as.bigq(diag(nrow=longsize))
  pos_template <- matrix(FALSE,nrow=longsize,ncol=longsize)
  for (it in seq_len(ncol(Lcompact))) {
    urange1 <- (it-1)*n_levels + seq(n_levels)
    positions <- pos_template
    positions[cbind(urange1,urange1)] <- TRUE
    longLv[[]][which(as.vector(positions))] <- Lcompact[it,it]
    for (jt in seq_len(it-1)) {
      urange2 <- (jt-1)*n_levels + seq(n_levels)
      positions <- pos_template
      positions[cbind(urange1,urange2)] <- TRUE
      longLv[[]][which(as.vector(positions))] <- Lcompact[it,jt]
      positions <- pos_template
      positions[cbind(urange2,urange1)] <- TRUE
      longLv[[]][which(as.vector(positions))] <- Lcompact[jt,it]
    }
  }
  return(longLv) 
} ## end def .makelong_bigq

.makelong_kronprod <- function(Lcompact, kron_Y) {
  if (methods::.hasSlot(Lcompact,"x") && any(Lcompact@x==0)) Lcompact <- drop0(Lcompact) # Seems cheap
  kronprod <- kronecker(Lcompact, kron_Y) # sparse as inferred from argument structure
  # so if the Lcompact has explicit zeros in @x, a drop0() would be useful. But do it on the LHS, not on the product!
  if ( inherits(Lcompact,c("dtCMatrix","ltCMatrix")) && 
       inherits(kron_Y,c("dtCMatrix","ltCMatrix"))) return(as(as(kronprod,"triangularMatrix"),"CsparseMatrix"))
  if ( inherits(Lcompact,c("dsCMatrix","lsCMatrix", "lsyMatrix", "dsyMatrix")) && 
       inherits(kron_Y,c("dsCMatrix","lsCMatrix"))) return(as(forceSymmetric(kronprod),"CsparseMatrix")) # forceSymmetric(dgT)  previously returned dsC, now dsT argh) 
  if ( inherits(kronprod,c("dtTMatrix","dgTMatrix"))) return(as(kronprod,"CsparseMatrix")) 
  return(kronprod) # *this* return occurs when LHS is full (square, not triangular) map matrix
}

## Build Lmatrix between all pairs of u_h (nr*nr) from parameter estimates (2*2) for the design matrix, longsize is final dim of matrix
.makelong <- function(Lcompact, ## may also be a compact precmat, symmetric rather than triangular
                      longsize,as_matrix=FALSE,template=NULL, kron_Y, kron_long=TRUE,
                      drop0template =TRUE # __F I X M E___ use =FALSE more
                      ) { ## always returns a Matrix unless explicitly as_matrix
  if ( ! is.null(kron_Y)) {
    if (kron_long) {
      return(.makelong_kronprod(Lcompact, kron_Y))
    } else return(.def_Kronfacto(lhs=Lcompact,rhs=kron_Y))
  }
  if (inherits(Lcompact, "bigq")) return(.makelong_bigq(Lcompact, longsize, kron_Y=kron_Y))
  if ( ! is.null(template)) { ## allows efficient filling of long template
    template@x <- Lcompact[template@x] ## template must be *M*atrix 
    if (as_matrix) {
      template <- as.matrix(template)
    } else {
      if (drop0template && .nonzeros(Lcompact)<ncol(Lcompact)^2) template <- drop0(template) 
      # hmf: Lcompact is dense so why the next code ? Lcompact is never a CsparseMatrix in the long routine tests. 
      if ( inherits(Lcompact,"dtCMatrix")) {
        template <- as(as(template,"triangularMatrix"), "CsparseMatrix")
      } else if ( inherits(Lcompact,"dsCMatrix")) template <- forceSymmetric(template)
    } 
    return(template) 
  } else { # no template... (both cases occur)
    longLv <- .C_makelong(Lcompact, longsize) ## efficient, "sparse matrix code" (input Lcompact not S4!) replaces in v3.6.39 previous R code.
    #  Matrix::kronecker(Lcompact,<Diagonal matrix>) can be used to obtain the 'same' (dgT...) Matrix but it's slower (even with precomputed Diag)
    if (as_matrix) {
      longLv <- as.matrix(longLv)
    } else {
      if ( inherits(Lcompact,"dtCMatrix")) { 
        longLv <- as(as(longLv,"triangularMatrix"), "CsparseMatrix")
      } else if ( inherits(Lcompact,"dsCMatrix")) { 
        longLv <- forceSymmetric(longLv)
      }
      #if (! is.null(template) && diff(range(template-longLv))>0) stop("zut") 
    }
    return(longLv) 
  }
} ## end def makelong

.makeCovEst1 <- function(u_h,ZAlist,cum_n_u_h,prev_LMatrices,
                         var_ranCoefs,
                         H_w.resid, 
                         processed,phi_est,
                         #family, ## ignored
                         as_matrix,v_h, MakeCovEst_pars_not_ZAL_or_lambda,
                         init_ranCoefs,
                         prev_lambda_est,
                         w.resid
) {
  rC_transf_inner <- .spaMM.data$options$rC_transf_inner
  nrand <- length(ZAlist)
  locX.Re <- processed$X.Re ## may be NULL 
  if (is.null(locX.Re)) locX.Re <- processed$AUGI0_ZX$X.pv
  if ( ncol(locX.Re) ) { lik_obj <- "p_bv"} else lik_obj <- "p_v" 
  updated_LMatrices <- working_LMatrices <- prev_LMatrices
  Xi_cols <- attr(ZAlist,"Xi_cols")
  loc_lambda_est <- prev_lambda_est # estimates not provided by the loop are necessary when (! augZXy_cond) -> .solve_IRLS_as_ZX()
  # currently prev_lambda_est holds non-unit elements at least of first iteration
  spprecBool <- processed$is_spprec
  use_tri_CORREL <- .spaMM.data$options$use_tri_for_makeCovEst
  verbose <- processed$verbose["TRACE"]
  augZXy_cond <- attr(processed$augZXy_cond,"inner")
  test <- FALSE
  #test <- TRUE
  pars_for_conv_corr <- vector("list",nrand)
  H_global_scale <- .calc_H_global_scale(H_w.resid)
  for (rt in seq_len(nrand)) {
    if ( var_ranCoefs[rt]) { ## inner estimation of cov mat of u_h 
      Xi_ncol <- Xi_cols[rt]
      n_levels <- ncol(ZAlist[[rt]])/Xi_ncol 
      COVpredUnderHuncorr <- matrix(0,ncol=Xi_ncol,nrow=Xi_ncol) ## var on diag, corr outside diag
      ##prevL <- attr(prev_LMatrices[[rt]],"latentL_blob")$design_u
      u.range <- (cum_n_u_h[rt]+1L):(cum_n_u_h[rt+1L])
      loc_lambda_est[u.range] <- 1 # not so on first call of .makeCovEst1() at least. Could it be set before calling .makeCovEst1()? 
      ranCoefs_blob <- processed$ranCoefs_blob
      ########## brute force optimization
      objfn <- function(trRancoef) { ## to be MINimized
        compactcovmat <- .calc_cov_from_trRancoef(trRancoef, Xi_ncol, rC_transf=rC_transf_inner)
        latentL_blob <- .calc_latentL(compactcovmat, use_tri_CORREL=use_tri_CORREL, spprecBool=spprecBool, trDiag=ranCoefs_blob$trDiags[[rt]])
        # Build Xscal
        working_LMatrices[[rt]] <- .makelong(latentL_blob$design_u,longsize=ncol(ZAlist[[rt]]), 
                                             template=processed$ranCoefs_blob$longLv_templates[[rt]], 
                                             kron_Y=attr(processed$corr_info$cov_info_mats[[rt]],"blob")$Lunique 
        ) ## the variances are taken out in $d
        ###### attr(working_LMatrices[[rt]],"ranefs") <- attr(ZAlist,"exp_ranef_strings")[[rt]] 
        locZAL <- .compute_ZAL(XMatrix=working_LMatrices, ZAlist=ZAlist, as_matrix=as_matrix, force_bindable=FALSE) 
        if (augZXy_cond || test) {
          ####################################################################################################
          # we don't want anything specific on u_h values:
          w.ranef <- 1/loc_lambda_est # call to .updateW_ranefS() reduced to this for v3.6.39
          ZAL_scaling <- sqrt(loc_lambda_est/H_global_scale) # sqrt(w.ranef*H_global_scale) ## Q^{-1/2}/s
          Xscal <- .make_Xscal(ZAL=locZAL, ZAL_scaling = ZAL_scaling, processed=processed)
          weight_X <- .calc_weight_X(Hobs_w.resid=H_w.resid,
                                     H_global_scale=H_global_scale, obsInfo=processed$how$obsInfo) ## sqrt(s^2 W.resid) ## should not affect the result up to precision
          sXaug <- do.call(processed$corr_method,
                           list(Xaug=Xscal, weight_X=weight_X, w.ranef=w.ranef, H_global_scale=H_global_scale))
          ####################################################################################################
        } else sXaug <- NULL
        if ( ! augZXy_cond || test) {
          ####################################################################################################
          locw.ranefSblob <- processed$updateW_ranefS(u_h=u_h,v_h=v_h,lambda=loc_lambda_est)
          locarglist <- c(MakeCovEst_pars_not_ZAL_or_lambda, 
                          list(ZAL=locZAL, lambda_est=loc_lambda_est, wranefblob=locw.ranefSblob, H_global_scale=H_global_scale, 
                               H_w.resid=H_w.resid, w.resid=w.resid))
          # it would be really nonsense to compute the objective for constant u_h;
          ## the u_h should be evaluated at the BLUPs (or equivalent) for given parameters
          auglinmodblob <- do.call(".solve_IRLS_as_ZX",locarglist)
          ## FR->FR but its not clear that non-standard REML is handled by calc_APHLs_from_ZX !!
          ####################################################################################################
          obj_augZX <- .calc_APHLs_from_ZX(sXaug=NULL,
                                           auglinmodblob=auglinmodblob, # including $lambda_est, needed here
                                           phi_est=phi_est,
                                           processed=processed,
                                           which=lik_obj)[[lik_obj]]
        } 
        if ( augZXy_cond || test) {
          obj_augZXy <- .calc_APHLs_by_augZXy_or_sXaug(sXaug=sXaug, # with this we don't need a lambda_est argument.
                                           auglinmodblob=NULL, 
                                           phi_est=phi_est,
                                           processed=processed,
                                           which=lik_obj)[[lik_obj]]
          
        }
        if (test) .test_augZXy(obj_augZXy, obj_augZX, phi_est=phi_est)
        if (augZXy_cond) {objective <- obj_augZXy} else objective <- obj_augZX
        return( - objective)
      } 
      ####
      init_trRancoef <- attr(prev_LMatrices[[rt]],"trRancoef")
      if (is.null(init_trRancoef)) {
        init_ranCoef <- init_ranCoefs[[as.character(rt)]]
        if (is.null(init_ranCoef)) {
          init_ranCoef <- diag(x=rep(1, Xi_ncol))
          #init_ranCoef[lower.tri(init_ranCoef,diag=FALSE)] <- 0.001
        } 
        init_trRancoef <- .ranCoefsFn(init_ranCoef[lower.tri(init_ranCoef,diag=TRUE)], rC_transf=rC_transf_inner) # rep(0,Xi_ncol*(Xi_ncol+1L)/2L)
      }
      trRancoef_LowUp <- .calc_LowUp_trRancoef(init_trRancoef, Xi_ncol=Xi_ncol,
                                               tol_ranCoefs=.spaMM.data$options$tol_ranCoefs_inner,
                                               rC_transf=rC_transf_inner)
      lowerb <- trRancoef_LowUp$lower
      upperb <- trRancoef_LowUp$upper
      Optimizer <- .spaMM.data$options$optim_inner
      if (Optimizer==".safe_opt") { # note that optimizer is not preprocessed (and that would require some thoughts)
        optr <- .safe_opt(init=init_trRancoef, objfn=objfn, lower=lowerb, upper=upperb, verbose=verbose,
                          adjust_init=trRancoef_LowUp$adjust_init, LowUp=list(lower=list(trRanCoefs=list("1"=trRancoef_LowUp$lower)),
                                                                              upper=list(trRanCoefs=list("1"=trRancoef_LowUp$upper)))) 
      } else if (Optimizer=="nloptr") {
        optr <- .optim_by_nloptr(initvec=init_trRancoef,lowerb=lowerb,upperb=upperb,objfn_locoptim=objfn,
                                 local_control=NULL,LowUp=trRancoef_LowUp)
      } else if (Optimizer=="bobyqa") {
        optr <- .optim_by_bobyqa(initvec=init_trRancoef,lowerb=lowerb,upperb=upperb,objfn_locoptim=objfn,
                                 local_control=NULL)
        optr$objective <- optr$fval
        optr$solution <- optr$par
      }
      ################# 
      ## reproduces representation in objfn
      compactcovmat <- .calc_cov_from_trRancoef(optr$solution, Xi_ncol, rC_transf=rC_transf_inner)
      latentL_blob <- .calc_latentL(compactcovmat, use_tri_CORREL=use_tri_CORREL, spprecBool=spprecBool, trDiag=ranCoefs_blob$trDiags[[rt]])
      next_LMatrix <- .makelong(latentL_blob$design_u,longsize=ncol(ZAlist[[rt]]), 
                                template=processed$ranCoefs_blob$longLv_templates[[rt]], 
                                kron_Y=attr(processed$corr_info$cov_info_mats[[rt]],"blob")$Lunique 
                                ) ## il faut updater pour estimer les ranef correctement...
      attr(next_LMatrix,"latentL_blob") <- latentL_blob ## kept for updating in next iteration and for output
      attr(next_LMatrix,"trRancoef") <- optr$solution ## kept for updating in next iteration and for output
      ranef_info <- attr(ZAlist,"exp_ranef_strings")[rt]
      attr(ranef_info, "type") <- attr(ZAlist,"exp_ranef_types")[rt]## type is "(.|.)" if LMatrix is for random slope ## ajout 2015/06
      #####  attr(next_LMatrix,"ranefs") <-  ranef_info
      attr(next_LMatrix, "corr.model") <- "random-coef"
      updated_LMatrices[[rt]] <- next_LMatrix # (fixme ??) working_LMatrices[[rt]] <-
      ### in the HLfit6 example, one element of compactcovmat appears to move without affecting logL, preventing compactcovmat convergence
      ## the optpars for "chol" appear to converge, but not those for "sph"
      #pars_for_conv_corr[[rt]] <- compactcovmat[lower.tri(compactcovmat,diag=TRUE)] 
      pars_for_conv_corr[[rt]] <- optr$solution 
    } else updated_LMatrices[rt] <- list(NULL) ## this erases Matern, AR1, and other fixed LMatrices, so updated_LMatrices is not to be used in objfn()   
  } ## loop on rt = ranefs
  if (verbose["TRACE"]>1L) print(optr)
  return(list(updated_LMatrices=updated_LMatrices, 
              next_lambda_est=loc_lambda_est,
              info_for_conv_rC=list(obj=optr$objective, ranCoefs=unlist(pars_for_conv_corr))))
} ## end def makeCovEst1
