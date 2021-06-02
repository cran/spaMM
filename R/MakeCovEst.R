# HLfit(Reaction ~ 0 + (Days|Subject), data = sleepstudy,HLmethod="ML")
# HLfit(Reaction ~ 1 + (Days|Subject), data = sleepstudy,HLmethod="ML")
# HLfit(Reaction ~ Days + (Days|Subject), data = sleepstudy,HLmethod="ML")

.calc_latentL <- function(compactcovmat, use_tri_Nspprec=TRUE, spprecBool, try_chol=FALSE) {
  # returns a *t*crossfactor 
  ## L and L_Q are distinct matrices used jointly in a fit (L in ZAL...) andmust be deduced from each other.
  #  In sparse precision code We want chol_Q to be lower triangular (dtCMatrix can be Up or Lo) 
  #    to obtain a dtCMatrix by bdiag(list of lower tri dtCMatrices).
  # Beyond this, design_u may be lower or upper tri
  if (try_chol) {
    if (is.null(chol_crossfac <- attr(compactcovmat, "chol_crossfac"))) {
      chol_crossfac <- try(chol(compactcovmat), silent=TRUE)
    }
    if ( ! inherits(chol_crossfac, "try-error")) return(list(design_u=t(chol_crossfac), 
                                                             d=rep(1,ncol(compactcovmat)), 
                                                             compactcovmat=compactcovmat))
  }
  # ELSE:
  compactcovmat <- .smooth_regul(compactcovmat, epsi=.spaMM.data$options$tol_ranCoefs_inner["regul"]) 
  #   value of epsi utlmately also affects whether bobyqa has to be run 
  esys <- attr(compactcovmat,"esys") 
  d_regul <- esys$d_regul
  if (spprecBool || use_tri_Nspprec) { ## Competing for clumsy code prize... ## sparse version elsewhere
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
    # now with triangular factor:
    invdcp <- 1/.diagfast(x=tcrossfac_prec)
    compactchol_Q <- .m_Matrix_times_Dvec(tcrossfac_prec, invdcp) ## for design_u, and independently in sparse precision algo. ## compactprecmat=.ZWZtwrapper(chol_Q,1/invdcp^2)
    blob <- list(design_u=t(solve(compactchol_Q)), # tcrossprod factor, not orthonormed, UPPER tri (while eigen provides an orthonormed "full" U matrix)
                 d=invdcp^2, # given variances of latent independent ranefs outer optim algo and iterative algo
                 compactcovmat=compactcovmat ## not used for the fit
    ) 
    if (spprecBool) {
      ## (1) direct conversion from matrix to dtC is remarkably slower than going through the as(.,"sparseMatrix") step !
      ## (2) Even with two steps as(.,"dtCMatrix") is costly hence the spprecBool condition (dtCMatrix useful only for chol_Q)
      compactchol_Q <- as(compactchol_Q,"sparseMatrix") # precision factor for sparse_precision algo
      blob$compactchol_Q <- as(compactchol_Q,"dtCMatrix") # precision factor for sparse_precision algo; # compactcovmat=.ZtWZwrapper(solve(compactchol_Q),d  )
      evec <- as(esys$vectors,"dgCMatrix")
      blob$compactprecmat <- .ZWZtwrapper(evec, 1/d_regul) # dsCMatrix
    } else blob$compactchol_Q <- compactchol_Q ## not used for the fit
    #
    # if (identical(spaMM.data$options$safe_spprec,TRUE)) { # debug algo: d=1 and design_u is t(solve(compactchol_prec)) not t(solve(compactchol_Q))
    #   cat(crayon::red("TRUE"))
    #   blob$design_u <- blob$design_u %*% diag(x=invdcp) # OK for corr not for prec
    #   blob$d <- rep(1,length(invdcp))
    #   blob$compactchol_Q <- as(as(tcrossfac_prec,"sparseMatrix"),"dtCMatrix") # precision factor for sparse_precision algo; # compactcovmat=.ZtWZwrapper(solve(compactchol_Q),d  )
    #   blob$compactprecmat <- solve(blob$compactcovmat)
    # } else {
    #   cat(crayon::green("FALSE"))
    #   # design_u is t(solve(tcrossfac_Q))=t(solve(compactchol_Q)) and $d is not necess 1
    # }
  } else if (FALSE) { ## Would be much nicer if it was OK for HLfit
    # qr always produces a crossprod factor so if we want a tcrossprod factor, it must be t(qr.R()) hence lower tri...); unless we qr() the preci mat...
    qrand <- t(.m_Matrix_times_Dvec(esys$vectors,sqrt(d_regul)))
    qrblob <- qr(qrand)
    crossfac <- qr.R(qrblob) # applying .lmwithQR() systematically (badly) affects numerical precision
    if (! all(unique(diff(qrblob$pivot))==1L)) { # eval an unpermuted triangular R
      crossfac <- .lmwithQR(crossfac[, sort.list(qrblob$pivot)] ,yy=NULL,returntQ=FALSE,returnR=TRUE)$R
    } 
    tcrossfac <- t(crossfac) # lower.tri
    sqrt_d <- .diagfast(x=tcrossfac) 
    design_u <- .m_Matrix_times_Dvec(tcrossfac, 1/sqrt_d)
    blob <- list(design_u=design_u, # LOWER tri tcrossprod factor
                 d= sqrt_d^2, # rep(1,length(d_regul)), # The given variances of latent independent ranefs outer optim algo and iterative algo
                 compactcovmat=compactcovmat ## not used for the fit
    ) ## and one might use (solve(design_u)) as a $crossfac_Q precision factor
  } else {blob <- list(design_u=esys$vectors, 
                       d=d_regul, # *!* spprec case where d!=1; occurs when Cholesky failed.
                       compactcovmat=compactcovmat)}
  return(blob)
  ## for sparse precision we want chol_Q to be (dtCMatrix: Csparse triangular) so that efficient solve methods can be used.
  ## This is not provided by this function
}

.makelong_bigq <- function(Lcompact, longsize) { ## ad hoc version with restricted usage
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

## Build Lmatrix between all pairs of u_h (nr*nr) from parameter estimates (2*2) for the design matrix, longsize is final dim of matrix
.makelong <- function(Lcompact, ## may also be a compact precmat, symmetric rather than triangular
                      longsize,as_matrix=FALSE,template=NULL) { ## always returns a Matrix unless explicitly as_matrix
  
  if (inherits(Lcompact, "bigq")) return(.makelong_bigq(Lcompact, longsize))
  if ( ! is.null(template)) { ## allows efficient filling of long template
    template@x <- Lcompact[template@x]
    if (as_matrix) {
      template <- as.matrix(template)
    } else {
      #if ( ! inherits(template,"sparseMatrix")) stop()
      if (any(Lcompact==0L)) template <- drop0(template) ## LHS is *M*atrix in all cases    
      if ( inherits(Lcompact,"dtCMatrix")) {
        template <- as(template,"dtCMatrix")
      } else if ( inherits(Lcompact,"dsCMatrix")) template <- forceSymmetric(template)
    } 
    return(template) 
  } else { # no template... (both cases occur)
    longLv <- .C_makelong(Lcompact, longsize) ## efficient, "sparse matrix code" (input Lcompact not S4!) replaces in v3.6.39 previous R code.
    if (as_matrix) {
      longLv <- as.matrix(longLv)
    } else {
      if ( inherits(Lcompact,"dtCMatrix")) { 
        longLv <- as(longLv,"dtCMatrix")
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
                         w.resid, # not variable within the function
                         processed,phi_est,
                         #family, ## ignored
                        as_matrix,v_h, MakeCovEst_pars_not_ZAL_or_lambda,
                        init_ranCoefs,
                        prev_lambda_est
) {
  H_global_scale <- .calc_H_global_scale(w.resid)
  weight_X <- .calc_weight_X(w.resid, H_global_scale) ## sqrt(s^2 W.resid) ## should not affect the result up to precision
  rC_transf_inner <- .spaMM.data$options$rC_transf_inner
  nrand <- length(ZAlist)
  locX.Re <- processed$X.Re ## may be NULL 
  if (is.null(locX.Re)) locX.Re <- processed$AUGI0_ZX$X.pv
  if ( ncol(locX.Re) ) { lik_obj <- "p_bv"} else lik_obj <- "p_v" 
  updated_LMatrices <- working_LMatrices <- prev_LMatrices
  Xi_cols <- attr(ZAlist,"Xi_cols")
  loc_lambda_est <- prev_lambda_est # estimates not provided by the loop are necessary when (! augZXy_cond) -> .solve_IRLS_as_ZX()
  spprecBool <- processed$is_spprec
  use_tri_Nspprec <- .spaMM.data$options$use_tri_for_makeCovEst
  verbose <- processed$verbose["TRACE"]
  try_chol <- ( ! spprecBool ) # && (.spaMM.data$options[["rC_transf_inner"]]=="chol") ## "sph" works with chol latentL now
  augZXy_cond <- attr(processed$augZXy_cond,"inner")
  test <- FALSE
  #test <- TRUE
  pars_for_conv_corr <- vector("list",nrand)
  for (rt in seq_len(nrand)) {
    if ( var_ranCoefs[rt]) { ## inner estimation of cov mat of u_h 
      Xi_ncol <- Xi_cols[rt]
      n_levels <- ncol(ZAlist[[rt]])/Xi_ncol 
      COVpredUnderHuncorr <- matrix(0,ncol=Xi_ncol,nrow=Xi_ncol) ## var on diag, corr outside diag
      ##prevL <- attr(prev_LMatrices[[rt]],"latentL_blob")$design_u
      u.range <- (cum_n_u_h[rt]+1L):(cum_n_u_h[rt+1L])
      ########## brute force optimization
      objfn <- function(trRancoef) { ## to be MINimized
        compactcovmat <- .calc_cov_from_trRancoef(trRancoef, Xi_ncol, rC_transf=rC_transf_inner)
        latentL_blob <- .calc_latentL(compactcovmat, use_tri_Nspprec=use_tri_Nspprec, spprecBool=spprecBool, try_chol= try_chol)
        # Build Xscal
        working_LMatrices[[rt]] <- .makelong(latentL_blob$design_u,longsize=ncol(ZAlist[[rt]]), 
                                             template=processed$ranCoefs_blob$longLv_templates[[rt]]) ## the variances are taken out in $d
        attr(working_LMatrices[[rt]],"ranefs") <- attr(ZAlist,"exp_ranef_strings")[[rt]] 
        locZAL <- .compute_ZAL(XMatrix=working_LMatrices, ZAlist=ZAlist, as_matrix=as_matrix) 
        # w.ranef argument
        loc_lambda_est[u.range] <- .make_long_lambda(latentL_blob$d, n_levels, Xi_ncol) 
        ## latentL_blob$d may be a rep(1,...) in all calls to this fn. In which cas, below,
        ## locwranefblob and ZAL_scaling are constant; non-constant code begisn with Xscal(locZAL...)
        #if (.spaMM.data$options$rC_transf_inner!="chol") loc_lambda_est[loc_lambda_est<1e-08] <- 1e-08 # (two occurrences)
        if (augZXy_cond || test) {
          ####################################################################################################
          n_u_h <- length(u_h)
          # we don't want anything specific on u_h values:
          w.ranef <- 1/loc_lambda_est # call to .updateW_ranefS() reduced to this for v3.6.39
          ZAL_scaling <- sqrt(loc_lambda_est/H_global_scale) # sqrt(w.ranef*H_global_scale) ## Q^{-1/2}/s
          Xscal <- .make_Xscal(ZAL=locZAL, ZAL_scaling = ZAL_scaling, AUGI0_ZX=processed$AUGI0_ZX)
          if (inherits(Xscal,"Matrix")) { # same type as ZAL
            mMatrix_method <- .spaMM.data$options$Matrix_method
          } else {
            mMatrix_method <- .spaMM.data$options$matrix_method
          }
          sXaug <- do.call(mMatrix_method,
                           list(Xaug=Xscal, weight_X=weight_X, w.ranef=w.ranef, H_global_scale=H_global_scale))
          ####################################################################################################
        } else sXaug <- NULL
        if ( ! augZXy_cond || test) {
          ####################################################################################################
          W_ranefS_constant_args <- processed$reserve$W_ranefS_constant_args
          locw.ranefSblob <- do.call(".updateW_ranefS",c(W_ranefS_constant_args, list(u_h=u_h,v_h=v_h,lambda=loc_lambda_est)))
          locarglist <- c(MakeCovEst_pars_not_ZAL_or_lambda, list(ZAL=locZAL, lambda_est=loc_lambda_est, wranefblob=locw.ranefSblob))
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
      latentL_blob <- .calc_latentL(compactcovmat, use_tri_Nspprec=use_tri_Nspprec, spprecBool=spprecBool, try_chol= try_chol)
      loc_lambda_est[u.range] <- .make_long_lambda(latentL_blob$d, n_levels, Xi_ncol) 
      #if (.spaMM.data$options$rC_transf_inner!="chol") loc_lambda_est[loc_lambda_est<1e-08] <- 1e-08 
      next_LMatrix <- .makelong(latentL_blob$design_u,longsize=ncol(ZAlist[[rt]]), 
                                template=processed$ranCoefs_blob$longLv_templates[[rt]]) ## il faut updater pour estimer les ranef correctement...
      attr(next_LMatrix,"latentL_blob") <- latentL_blob ## kept for updating in next iteration and for output
      attr(next_LMatrix,"trRancoef") <- optr$solution ## kept for updating in next iteration and for output
      ranef_info <- attr(ZAlist,"exp_ranef_strings")[rt]
      attr(ranef_info, "type") <- attr(ZAlist,"exp_ranef_types")[rt]## type is "(.|.)" if LMatrix is for random slope ## ajout 2015/06
      attr(next_LMatrix,"ranefs") <-  ranef_info
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
