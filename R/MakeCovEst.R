# HLfit(Reaction ~ 0 + (Days|Subject), data = sleepstudy,HLmethod="ML")
# HLfit(Reaction ~ 1 + (Days|Subject), data = sleepstudy,HLmethod="ML")
# HLfit(Reaction ~ Days + (Days|Subject), data = sleepstudy,HLmethod="ML")

.calc_latentL <- function(compactcovmat, triangularL=TRUE, spprecBool) {
  ## the tcrossproduct factorization is not unique, so there's no reason for the factorization of compactcovmat and compactprecmat to match
  #  unless one is deduced from the other. L and L_Q are distinct matrices used jointly in a fit (L in ZAL...).
  #  Hence it's pathetic to recompute a Chol from the cov somewhere and another from the prec elsewhere.
  #  In sparse precision code We want chol_Q to be lower triangular (dtCMatrix can be Up or Lo) 
  #    to obtain a dtCMatrix by bdiag(list of lower tri dtCMatrices).
  # Beyond this, design_u may be lower or upper tri
  d_regul <- attr(compactcovmat,"d_regul")
  esys <- attr(compactcovmat,"esys")
  if (spprecBool || triangularL) { ## Competing for clumsy code prize... ## sparse version elsewhere
    ## need triangular factor for spprec case in particular hence the default value of triangularL
    crossfac_prec <- .Dvec_times_matrix(sqrt(1/d_regul), t(esys$vectors))
    crossfac_prec <- qr.R(qr(crossfac_prec))
    tcrossfac_prec <- t(crossfac_prec)  ## equivalent to chol up to signs... # tcrossprod(tcrossfac_prec)=compactprecmat
    # Correct signs:
    diagsigns <- sign(.diagfast(tcrossfac_prec))
    tcrossfac_prec <- .m_Matrix_times_Dvec(tcrossfac_prec,diagsigns)
    # now with triangular factor:
    invdcp <- 1/.diagfast(x=tcrossfac_prec)
    tcrossfac_Q <- .m_Matrix_times_Dvec(tcrossfac_prec, invdcp) ## for design_u, and independently in sparse precision algo. ## compactprecmat=.ZWZtwrapper(chol_Q,1/invdcp^2)
    tcrossfac_corr <- t(solve(tcrossfac_Q)) ## upper tri (while eigen provides a "full" U matrix)
    blob <- list(design_u=tcrossfac_corr, # UPPER tri tcrossprod factor
                 d=invdcp^2, # given variances of latent independent ranefs outer optim algo and iterative algo
                 compactcovmat=compactcovmat ## not used for the fit
    ) 
    if (spprecBool) {
      ## (1) direct conversion from matrix to dtC is remarkably slower than going through the as(.,"sparseMatrix") step !
      ## (2) Even with two steps as(.,"dtCMatrix") is costly hence the spprecBool condition (dtCMatrix useful only for chol_Q)
      tcrossfac_Q <- as(tcrossfac_Q,"sparseMatrix") # precision factor for sparse_precision algo
      blob$compactchol_Q <- as(tcrossfac_Q,"dtCMatrix") # precision factor for sparse_precision algo
      blob$compactprecmat <- .ZWZt(esys$vectors, 1/d_regul) 
    }
  } else if (triangularL) { ## never run as triangularL handled by first case, but would be much nicer if it was OK for HLfit6
    # qr always produces a crossprod factor so if we want a tcrossprod factor, it must be t(qr.R()) hence lower tri...); unless we qr() the preci mat...
    design_u <- t(qr.R(qr(t(.m_Matrix_times_Dvec(esys$vectors,sqrt(d_regul)))))) ## lower tri contrary to clumsy method
    sqrt_d <- .diagfast(x=design_u) 
    tcrossfac_corr <- .m_Matrix_times_Dvec(design_u, 1/sqrt_d)
    blob <- list(design_u=tcrossfac_corr, # LOWER tri tcrossprod factor
                 d= sqrt_d^2, # rep(1,length(d_regul)), # given variances of latent independent ranefs outer optim algo and iterative algo
                 compactcovmat=compactcovmat ## not used for the fit
    ) ## and one might use (solve(design_u)) as a $crossfac_Q precision factor
  } else {blob <- list(design_u=esys$vectors, 
                       d=d_regul, compactcovmat=compactcovmat)}
  return(blob)
  ## for sparse precision we want chol_Q to be (dtCMatrix: Csparse triangular) so that efficient solve methods can be used.
  ## This is not provided by this function
}

## Build Lmatrix between all pairs of u_h (nr*nr) from parameter estimates (2*2) for the design matrix, longsize is final dim of matrix
.makelong <- function(Lcompact,longsize,as_matrix=FALSE,template=NULL) { ## always returns a Matrix unless explicitly as_matrix
  if ( ! is.null(template)) {
    template@x <- Lcompact[template@x]
    if (as_matrix) {
      template <- as.matrix(template)
    } else {
      template <- drop0(template) ## LHS is *M*atrix in all cases
      if ( inherits(Lcompact,"dtCMatrix")) template <- as(template,"dtCMatrix")
    } 
    return(template) 
  } else {
    n_levels <- longsize/ncol(Lcompact)
    if (.spaMM.data$options$Zcolsbyrows) {
      longLv <- do.call(Matrix::bdiag,rep(list(Lcompact),n_levels)) # more diagonal... no benefits (yet)
    } else {
      if (longsize>180L) {
        longLv <- Diagonal(n=longsize) ## declaration ## FIXME with some effort, we could produce better sparse code.
      } else longLv <- diag(nrow=longsize) ## declaration
      for (it in seq_len(ncol(Lcompact))) {
        urange1 <- (it-1)*n_levels + seq(n_levels)
        longLv[cbind(urange1,urange1)] <- Lcompact[it,it]
        for (jt in seq_len(it-1)) {
          urange2 <- (jt-1)*n_levels + seq(n_levels)
          longLv[cbind(urange1,urange2)] <- Lcompact[it,jt]
          longLv[cbind(urange2,urange1)] <- Lcompact[jt,it]
        }
      }
    }
    if (as_matrix) {
      longLv <- as.matrix(longLv)
    } else {
      longLv <- drop0(longLv) ## drop0() returns a *M*atrix
      if ( inherits(Lcompact,"dtCMatrix")) { longLv <- as(longLv,"dtCMatrix")} 
      #if (! is.null(template) && diff(range(template-longLv))>0) stop("zut") 
    }
    return(longLv) 
  }
  return(longLv) 
} ## end def makelong

.makeCovEst1 <- function(u_h,ZAlist,cum_n_u_h,prev_LMatrices,
                         var_ranCoefs,w.resid,processed,phi_est,
                         #family, ## ignored
                        as_matrix,v_h, MakeCovEst_pars_not_ZAL_or_lambda
) {
  nrand <- length(ZAlist)
  locX.Re <- processed$X.Re ## may be NULL 
  if (is.null(locX.Re)) locX.Re <- processed$AUGI0_ZX$X.pv
  if ( ncol(locX.Re) ) { lik_obj="p_bv"} else lik_obj="p_v" 
  updated_LMatrices <- working_LMatrices <- prev_LMatrices
  Xi_cols <- attr(ZAlist,"Xi_cols")
  loc_lambda_est <- numeric(length(u_h))
  for (rt in seq_len(nrand)) {
    if ( var_ranCoefs[rt]) { ## inner estimation of cov mat of u_h 
      Xi_ncol <- Xi_cols[rt]
      n_levels <- ncol(ZAlist[[rt]])/Xi_ncol 
      COVpredUnderHuncorr <- matrix(0,ncol=Xi_ncol,nrow=Xi_ncol) ## var on diag, corr outside diag
      ##prevL <- attr(prev_LMatrices[[rt]],"latentL_blob")$design_u
      u.range <- (cum_n_u_h[rt]+1L):(cum_n_u_h[rt+1L])
      ########## brute force optimization
      #
      augZXy_cond <- attr(processed$augZXy_cond,"inner")
      test <- FALSE
      #test <- TRUE
      objfn <- function(trRancoef) {
        compactcovmat <- .calc_cov_from_trRancoef(trRancoef, Xi_ncol)
        compactcovmat <- .regularized_eigen(compactcovmat, condnum=.spaMM.data$options$condnum_for_latentL_inner)
        latentL_blob <- .calc_latentL(compactcovmat, triangularL=.spaMM.data$options$use_tri_for_makeCovEst,
                                      spprecBool=processed$sparsePrecisionBOOL)
        # Build Xscal
        working_LMatrices[[rt]] <- .makelong(latentL_blob$design_u,longsize=ncol(ZAlist[[rt]]), 
                                             template=processed$ranCoefs_blob$longLv_templates[[rt]]) ## the variances are taken out in $d
        attr(working_LMatrices[[rt]],"ranefs") <- attr(ZAlist,"exp_ranef_strings")[[rt]] 
        locZAL <- .compute_ZAL(XMatrix=working_LMatrices, ZAlist=ZAlist, as_matrix=as_matrix) 
        # w.ranef argument
        loc_lambda_est[u.range] <- .make_long_lambda(latentL_blob$d, n_levels, Xi_ncol) 
        loc_lambda_est[loc_lambda_est<1e-08] <- 1e-08 ## arbitrarily small eigenvalue is possible for corr=+/-1 even for 'large' parvec
        if (augZXy_cond || test) {
          ####################################################################################################
          n_u_h <- length(u_h)
          # we don't want anything specific on u_h values:
          locwranefblob <- .updateW_ranefS(processed$cum_n_u_h, processed$rand.families, lambda=loc_lambda_est, 
                                           u_h=rep(NA,n_u_h),v_h=rep(NA,n_u_h)) ## indeed
          H_global_scale <- .calc_H_global_scale(w.resid)
          ZAL_scaling <- 1/sqrt(locwranefblob$w.ranef*H_global_scale) ## Q^{-1/2}/s
          Xscal <- .make_Xscal(ZAL=locZAL, ZAL_scaling = ZAL_scaling, AUGI0_ZX=processed$AUGI0_ZX)
          if (inherits(Xscal,"Matrix")) { # same type as ZAL
            mMatrix_method <- .spaMM.data$options$Matrix_method
          } else {
            mMatrix_method <- .spaMM.data$options$matrix_method
          }
          weight_X <- .calc_weight_X(w.resid, H_global_scale) ## sqrt(s^2 W.resid) ## should not affect the result up to precision
          sXaug <- do.call(mMatrix_method,
                           list(Xaug=Xscal, weight_X=weight_X, w.ranef=locwranefblob$w.ranef, H_global_scale=H_global_scale))
          ####################################################################################################
        } else sXaug <- NULL
        if ( ! augZXy_cond || test) {
          ####################################################################################################
          locw.ranefSblob <- .updateW_ranefS(cum_n_u_h,processed$rand.families,lambda=loc_lambda_est,u_h,v_h) 
          locarglist <- c(MakeCovEst_pars_not_ZAL_or_lambda, list(ZAL=locZAL, lambda_est=loc_lambda_est, wranefblob=locw.ranefSblob))
          # it would be really nonsens to compute the objective for constant u_h;
          ## the u_h should be evaluated at the BLUPs (or equivalent) for given parameters
          auglinmodblob <- do.call(".solve_IRLS_as_ZX",locarglist)
          ## FR->FR but its not clear that non-standard REML is handled by calc_APHLs_from_ZX !!
          ####################################################################################################
          obj_augZX <- .calc_APHLs_from_ZX(sXaug=NULL, # may be NULL
                                           auglinmodblob=auglinmodblob,  # may be NULL
                                           phi_est=phi_est,
                                           processed=processed,
                                           which=lik_obj)[[lik_obj]]
        } 
        if ( augZXy_cond || test) {
          obj_augZXy <- .calc_APHLs_by_augZXy_or_sXaug(sXaug=sXaug, 
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
      if (is.null(init_trRancoef)) init_trRancoef <- rep(0,Xi_ncol*(Xi_ncol+1L)/2L)
      trRancoef_LowUp <- .calc_LowUp_trRancoef(init_trRancoef, Xi_ncol=Xi_ncol,
                                               tol_ranCoefs=.spaMM.data$options$tol_ranCoefs)
      lowerb <- trRancoef_LowUp$lower
      upperb <- trRancoef_LowUp$upper
      if (TRUE) {
        nloptr_controls <- spaMM.getOption("nloptr") 
        if (is.null(nloptr_controls$maxeval)) nloptr_controls$maxeval <-  eval(.spaMM.data$options$maxeval,list(initvec=init_trRancoef))
        if (is.null(nloptr_controls$xtol_abs)) {
          nloptr_controls$xtol_abs <- eval(.spaMM.data$options$xtol_abs, 
                                           list(LowUp=trRancoef_LowUp)) #, factors=c(rcLam=5e-8,rcCor=5e-7,others=5e-11))) 
        }
        optr <- nloptr::nloptr(x0=init_trRancoef,eval_f=objfn,lb=lowerb,ub=upperb,
                       opts=nloptr_controls)
        while (optr$status==5L) { ## 5 => termination bc maxeval has been reached
          prevlik <- optr$objective
          reinit <- pmax(lowerb,pmin(upperb,optr$solution))
          optr <- nloptr::nloptr(x0=reinit,eval_f=objfn,lb=lowerb,ub=upperb,
                                 opts=nloptr_controls)
          loc_ftol <- max(1e-8, optr$options$ftol_abs)
          if (- optr$objective < - prevlik+loc_ftol) break ## no progress in <= maxeval iterations
        }
      } else { ## deprecated:
        ################# OPTIM
        parscale <- (upperb-lowerb)
        optr <- optim(init_trRancoef,objfn,lower=lowerb,upper=upperb,method="L-BFGS-B", ## deprecated
                      control=list(parscale=parscale))
        optr$solution <- optr$par
      }
      ################# 
      ## reproduces representation in objfn
      compactcovmat <- .calc_cov_from_trRancoef(optr$solution, Xi_ncol)
      compactcovmat <- .regularized_eigen(compactcovmat, condnum=.spaMM.data$options$condnum_for_latentL_inner)
      latentL_blob <- .calc_latentL(compactcovmat, triangularL=.spaMM.data$options$use_tri_for_makeCovEst,
                            spprecBool=processed$sparsePrecisionBOOL) 
      loc_lambda_est[u.range] <- .make_long_lambda(latentL_blob$d, n_levels, Xi_ncol) 
      loc_lambda_est[loc_lambda_est<1e-08] <- 1e-08 
      next_LMatrix <- .makelong(latentL_blob$design_u,longsize=ncol(ZAlist[[rt]]), 
                                template=processed$ranCoefs_blob$longLv_templates[[rt]]) ## il faut updater pour estimer les ranef correctement...
      attr(next_LMatrix,"latentL_blob") <- latentL_blob ## kept for updating in next iteration and for output
      attr(next_LMatrix,"trRancoef") <- optr$solution ## kept for updating in next iteration and for output
      attr(next_LMatrix,"ranefs") <-  structure(attr(ZAlist,"exp_ranef_strings")[rt], 
                                                ## type is "(.|.)" if LMatrix is for random slope ## ajout 2015/06
                                                type= attr(ZAlist,"exp_ranef_types")[rt] )
      attr(next_LMatrix, "corr.model") <- "random-coef"
      updated_LMatrices[[rt]] <- next_LMatrix # (fixme ??) working_LMatrices[[rt]] <- 
    } else updated_LMatrices[rt] <- list(NULL) ## this erases Matern, AR1, and other fixed LMatrices, so updated_LMatrices is not to be used in objfn()   
  } ## loop on rt = ranefs
  return(list(updated_LMatrices=updated_LMatrices, 
              next_lambda_est=loc_lambda_est,
              optr_par=optr$solution))
} ## end def makeCovEst1
