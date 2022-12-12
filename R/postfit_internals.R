.calc_pd_product <- function(tcrossfac_v_beta_cov, Md2clikdvb2, blockSize=5000L) {
  if ((nc <- ncol(tcrossfac_v_beta_cov))>(blockSize)) {
    ## We reached this point by sparse matrix computations. We need to save memory in the following dense computations,
    message(paste0("Conditional-AIC computation requires operations on a large matrix (square with dimension ",nc,"),\n",
                   "which may take a bit of time. Use 'also_cAIC=FALSE' to avoid it.")) 
    slices <- unique(c(seq(0L,nc,blockSize),nc))
    nslices <- length(slices)-1L
    it <- 0L ## 'global definition' for Rcmd check
    foreach_args <- list(it = seq_len(nslices), .combine = "sum")
    foreach_blob <- do.call(foreach::foreach,foreach_args)
    abyss <- foreach::`%do%`(foreach_blob, Sys.setenv(LANG = "en"))
    barstyle <- eval(spaMM.getOption("barstyle"))
    progrbar_setup <- .set_progrbar(max = nslices, style = barstyle, char="s")
    pd <- foreach::`%do%`(foreach_blob, {
      slice <- (slices[it]+1L):slices[it+1L]
      tmp <- t(.crossprod(tcrossfac_v_beta_cov, Md2clikdvb2[,slice]))
      if (barstyle) progrbar_setup$progress(it)
      return(sum(tcrossfac_v_beta_cov[slice,] * tmp)) 
    })
    if (barstyle) close(progrbar_setup$pb)
    # We could parallelize using %dopar% (twice) but is that worth the overhead? 
  } else {
    # logic of following code is
    # pd = sum(diag(solve(Md2hdbv2,Md2clikdbv2[c(113:114,1:112),]))) 
    #    = sum(diag((tcrossprod(R_invMd2hdvb2)[c(113:114,1:112),] %*% Md2clikdbv2[sort.list(c(113:114,1:112)),]))) 
    #    = sum(diag(R_invMd2hdvb2[c(113:114,1:112),] %*% (crossprod(R_invMd2hdvb2, Md2clikdbv2[sort.list(c(113:114,1:112)),])))) 
    #    = sum((R_invMd2hdvb2[c(113:114,1:112),] * t(crossprod(R_invMd2hdvb2, Md2clikdbv2[sort.list(c(113:114,1:112)),])))) 
    # but we directly use v,b order rather than b,v
    pd <- t(.crossprod(tcrossfac_v_beta_cov,Md2clikdvb2))
    pd <- sum(tcrossfac_v_beta_cov * pd) 
    # that is actually the logic of .traceAB:
    # pd = .traceAB(tcrossfac_v_beta_cov, t(tcrossfac_v_beta_cov),Md2clikdvb2,diag(nrow =ncol(Md2clikdvb2)))
    # but the matrices are square so no use here
  }
  return(pd)
}


.calc_cAIC_pd_spprec <- function(object) {
  H_w.resid <- .get_H_w.resid(object)
  #ZAL <- .compute_ZAL(XMatrix=object$strucList, ZAlist=object$ZAlist,as_matrix=.eval_as_mat_arg.HLfit(object)) 
  ZAL <- get_ZALMatrix(object, force_bind=FALSE) # allows ZAXlist; but if a non-ZAXlist is already in the $envir, no effect; + we will solve(chol_Q, Diagonal()) so the gain is not obvious
  if ( ncol(object$X.pv) ) {
    M12 <- .crossprod(ZAL, .Dvec_times_m_Matrix(H_w.resid, object$X.pv), as_mat=TRUE)
    Md2clikdvb2 <- rbind2(cbind2(as.matrix(.ZtWZwrapper(ZAL,H_w.resid)), M12), ## this .ZtWZwrapper() takes time
                          cbind2(t(M12), as.matrix(.ZtWZwrapper(object$X.pv,H_w.resid)))) 
    # _FIXME_ any way to avoid formation of this matrix ? Or otherwise message() ?           
  } else {
    Md2clikdvb2 <-  as.matrix(.ZtWZwrapper(ZAL,H_w.resid))
  }
  tcrossfac_v_beta_cov <- .calc_Md2hdvb2_info_spprec(X.pv=object$X.pv, envir=object$envir, 
                                                     which="tcrossfac_v_beta_cov") 
  # not triang if we used sparse QR. Following code should not assume triangularity
  pd <- .calc_pd_product(tcrossfac_v_beta_cov, Md2clikdvb2)
  return(pd)
}

.calc_cAIC_pd_from_sXaug <- function(object) {
  if (is.matrix(beta_cov_info <- object$envir$beta_cov_info) || ## matrix is old format, should be a list now
      is.null(tcrossfac_beta_v_cov <- beta_cov_info$tcrossfac_beta_v_cov)) {
    tcrossfac_beta_v_cov <- .get_beta_cov_info(object)$tcrossfac_beta_v_cov
  }    
  pforpv <- ncol(object$X.pv)
  nc <- ncol(tcrossfac_beta_v_cov)
  n_u_h <- nc - pforpv
  seqp <- seq_len(pforpv)
  perm <- c(pforpv+seq_len(n_u_h), seqp)
  tcrossfac_v_beta_cov <- tcrossfac_beta_v_cov[perm,,drop=FALSE] # useful to keep it for predVar computations?
  
  H_w.resid <- .get_H_w.resid(object)
  #ZAL <- .compute_ZAL(XMatrix=obje17ct$strucList, ZAlist=object$ZAlist,as_matrix=.eval_as_mat_arg.HLfit(object)) 
  ZAL <- get_ZALMatrix(object, force_bind=FALSE) # allows ZAXlist; but if a non-ZAXlist is already in the $envir, no effect; + we will solve(chol_Q, Diagonal()) so the gain is not obvious
  if ( ncol(object$X.pv) ) {
    M12 <- .crossprod(ZAL, .Dvec_times_m_Matrix(H_w.resid, object$X.pv), as_mat=TRUE)
    Md2clikdvb2 <- rbind2(cbind2(as.matrix(.ZtWZwrapper(ZAL,H_w.resid)), M12), ## this .ZtWZwrapper() takes time
                          cbind2(t(M12), as.matrix(.ZtWZwrapper(object$X.pv,H_w.resid)))) 
    # _FIXME_ any way to avoid formation of this matrix ? Or otherwise message() ?           
  } else {
    Md2clikdvb2 <-  as.matrix(.ZtWZwrapper(ZAL,H_w.resid))
  }
  # not triang if we used sparse QR. Following code should not assume triangularity
  pd <- .calc_pd_product(tcrossfac_v_beta_cov, Md2clikdvb2)
  return(pd)
}

.calc_cAIC_pd_others <- function(X.pv, ZAL, w.resid, d2hdv2, blockSize=1000L) {
  if ( ncol(X.pv) ) { ## the projection matrix for the response always includes X even for REML!
    hessnondiag <- .crossprod(ZAL, .Dvec_times_m_Matrix(w.resid, X.pv))
    Md2hdbv2 <- as.matrix(rbind2(cbind2(.ZtWZwrapper(X.pv,w.resid), t(hessnondiag)),
                                 cbind2(hessnondiag, - d2hdv2))) 
    Md2clikdbv2 <- as.matrix(rbind2(cbind2(.ZtWZwrapper(X.pv,w.resid), t(hessnondiag)),
                                    cbind2(hessnondiag, .ZtWZwrapper(ZAL,w.resid))))            
  } else {
    Md2hdbv2 <- - d2hdv2 
    Md2clikdbv2 <-  as.matrix(.ZtWZwrapper(ZAL,w.resid))
  }
  if (inherits(Md2hdbv2,"diagonalMatrix")) {
    pd <- sum(diag(Md2clikdbv2)/diag(Md2hdbv2)) ## is sum(diag(solve(Md2hdbv2,Md2clikdbv2)))
  } else {
    ## dans un LMM avec estimation ML, pd = sum(lev_phi), mais pas de simplif plus generale 
    if ((nc <- ncol(Md2hdbv2))>(blockSize)) {
      message(paste0("Conditional-AIC computation requires operations on a large dense matrix (square with dimension ",nc,"),\n",
                     "which may take a lot of time. Use 'also_cAIC=FALSE' to avoid it."))
    }
    ## if we reached this point with a huge dense matrix, then we have huge memory, hence saving memory may not be the issue,
    #  But we might save computation by computing only the requested diagonal (one qr, one backsolve, one sum(. * .))
    # using Md2hdbv2=QRP, tr=sum_i (invP invR Qt(Md2hdbv2))_ii
    qr.Md2hdbv2 <- try(qr(Md2hdbv2))
    if (inherits(qr.Md2hdbv2,"try-error")) {
      warning("Computation of cAIC/GoF df's failed because the information matrix appears singular.")
      pd <- NA
    } else {
      solveR <- try(backsolve(qr.R(qr.Md2hdbv2),diag(nrow =nc)))
      if (inherits(solveR,"try-error")) {
        warning("Computation of cAIC/GoF df's failed because the information matrix appears singular.")
        # determinant(qrR,logarithm=FALSE)$modulus > 1e-14
        pd <- NA
      } else {
        # using inv(RP)= inv(R[,perm <- sort.list(.$pivot)]) = inv(R)[perm,] (as in h9[,perm] %*% solve(h9)[perm,] for any perm)
        # though use of pivot may occur only for nearly singular matrices (poorly doc)
        pd <- sum(solveR[sort.list(qr.Md2hdbv2$pivot),] * t(qr.qty(qr.Md2hdbv2,Md2clikdbv2)))      
      }
    }
  }
  return(pd)
}

# Consistent with Saefken et al, only in terms of theta and mu. The linear predictor eta and link do not appear.
# Simulations for that paper use predictions eta only in the poisson(log) case;
# In the exponential (->Gamma(log)) they do use the theta deduced as -1/mu, not the eta. => cAIC4:::conditionalBootstrap is odd...
.calc_boot_AIC_dfs <- function (object, nsim, type="residual", seed=NULL, # (___F I X M E___) The whole boot procedure does not handle mv fits ?
                           nb_cores=NULL, fit_env=NULL) {
  if ( ! object$family$flags$exp) 
    stop("Bootstrap bias correction not implemented for families not from GLM (exponential family) class.") # assuming $exp methods are always available for GLMs (but see negbin2_dvl)
  bootsims <- simulate(object, nsim = nsim, type = type, verbose=FALSE, seed=seed) # the corresponding lmer code returns a data frame
  muFREQS <- dopar(bootsims, function(x) {
    predict(update_resp(object, newresp=x), type="response") # predict(refit(object, newresp = x))
  }, nb_cores=nb_cores, fit_env=fit_env)
  thetas <- .theta.mu.canonical(muFREQS,family(object)) # would return mu for non-GLMs
  # if (is.factor(bootsims[1])) dataMatrix <- as.numeric(dataMatrix) - 1 # in cAIC4:::conditionalBootstrap
  bootsims <- bootsims - rowMeans(bootsims)
  phis <- residVar(object, which="phi") # For Gamma(), the phis are still those of he canonical form of the exponential family, 
  # not the full residual variance (name 'get_residVar' is ambiguous)
  bootBC <- sum(colSums(thetas * bootsims/phis))/(nsim-1) # Handles the case Where the phis are heteroscedastic, 
  return(bootBC)
}

.calc_p_phi <- function(object, dfs=object$dfs) {
  p_phi <- dfs[["p_fixef_phi"]]
  if  ( ! is.null(resid_fits <- object$resid_fits)) { 
    p_phi <- sum(na.omit(unlist(p_phi)),
                 sum(unlist(lapply(resid_fits, `[[`, x="dfs"), recursive = TRUE, use.names = FALSE)) )
  } else if  ( ! is.null(resid_fit <- object$resid_fit)) { 
    # input p_phi (above) is typically set to NA, and will be ignored
    p_phi <- sum(.unlist(resid_fit$dfs)) ## phi_pd is relevant only for measuring quality of prediction by the resid_fit! 
  } else p_phi <- sum(.unlist(dfs[["p_fixef_phi"]])) # .unlist() for mv
  names_est_ranefPars <- unlist(.get_methods_disp(object))  
  fam_disp_parsnames <- intersect(names_est_ranefPars,c("NB_shape","COMP_nu","beta_prec"))
  if (length(fam_disp_parsnames)) {
    p_GLM_family <- length( # compatible with mv:
      .unlist(.get_outer_inits_from_fit(object, keep_canon_user_inits = FALSE)[fam_disp_parsnames]))
    p_phi <- p_phi+p_GLM_family ## effectively part of the model for residual error structure
  }
  p_phi
}

.get_info_crits <- function(object, also_cAIC=TRUE, nsim=0L, ...) { 
  if (is.null(info_crits <- object$envir$info_crits) || (also_cAIC && is.null(info_crits[["cAIC"]]))) { 
    dfs <- object$dfs
    if ( ! inherits(dfs,"list")) dfs <- as.list(dfs) ## back compatibility
    pforpv <- dfs[["pforpv"]]
    p_phi <- .calc_p_phi(object, dfs)
    p_lambda <- dfs[["p_lambda"]]
    APHLs <- object$APHLs
    H_w.resid <- .get_H_w.resid(object)
    info_crits <- list()
    # poisson-Gamma and negbin should have similar similar mAIC => NB_shape as one df or lambda as one df   
    forAIC <- APHLs
    if (object$models[[1]]=="etaHGLM") {
      if (object$HL[1]=="SEM") {
        forAIC <- list(p_v=APHLs$logLapp,p_bv=APHLs$logLapp,clik=APHLs$clik)
      } 
      # if standard ML: there is an REMLformula ~ 0 (or with ranefs ?); processed$X.Re is 0-col matrix
      # if standard REML: REMLformula is NULL: $X.Re is X.pv, processed$X.Re is NULL
      # non standard REML: other REMLformula: $X.Re and processed$X.Re identical, and may take essentially any value
      # if (identical(attr(object$REMLformula,"isML"),TRUE)) {
      #   Md2hdbv2 <- - d2hdv2 
      #   Md2clikdbv2 <-  as.matrix(.ZtWZwrapper(ZAL,w.resid))
      # } else {
      #   ## REML standard || REML non standard
      #   X.Re <- object$distinctX.Re ## null if not distinct from X.pv
      #   if (is.null(X.Re)) X.Re <- object$X.pv ## standard REML
      #   ## diff de d2hdbv2 slmt dans dernier bloc (-> computation pd)
      #   hessnondiag <- .crossprod(ZAL, .Dvec_times_m_Matrix(w.resid, X.pv))
      #   Md2hdbv2 <- as.matrix(rbind2(cbind2(.ZtWZwrapper(X.Re,w.resid), t(hessnondiag)),
      #                                cbind2(hessnondiag, - d2hdv2))) 
      #   Md2clikdbv2 <- as.matrix(rbind2(cbind2(.ZtWZwrapper(X.Re,w.resid), t(hessnondiag)),
      #                                   cbind2(hessnondiag, .ZtWZwrapper(ZAL,w.resid))))            
      # }
      X.pv <- object$X.pv
      if (is.null(p_corrPars <- dfs[["p_corrPars"]])) { ## back compatibility code, 
        # old code, presumably missing inner-estimated adjacency params 
        #corrPars <- get_ranPars(object,which="corrPars")
        #p_corrPars <- length(intersect(names_est_ranefPars,names(corrPars)))
        # Here inner-estimated adjacency params are in the "vars" 
        p_corrPars <- length(which(unlist(attr(object$CorrEst_and_RanFix,"type")$corrPars, use.names = FALSE) %in% c("outer","var")))
        # This is not the full count for up-to date spaMM (hyper param are missing), 
        # but should be OK for old objects to which this back compat code applies.
      }
      info_crits$mAIC <- -2*forAIC$p_v + 2 *(pforpv+p_lambda+p_corrPars+p_phi)
      info_crits$dAIC <- -2*forAIC$p_bv + 2 * (p_lambda+p_phi+p_corrPars) ## HaLM07 (eq 10) focussed for dispersion params
      #                                                                             including the rho param of an AR model
      if (also_cAIC) {
        if (.is_spprec_fit(object)) {
          pd <- .calc_cAIC_pd_spprec(object)
        } else if ( ! is.null(object$envir$sXaug)) { 
          pd <- .calc_cAIC_pd_from_sXaug(object)
        } else { # SEM
          ZAL <- .compute_ZAL(XMatrix=object$strucList, ZAlist=object$ZAlist,as_matrix=.eval_as_mat_arg.HLfit(object)) 
          d2hdv2 <- .calcD2hDv2(ZAL,H_w.resid,object$w.ranef) ## update d2hdv2= - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
          pd <- .calc_cAIC_pd_others(X.pv, ZAL, H_w.resid, d2hdv2)
        }
        info_crits$GoFdf <- length(object$y) - pd ## <- nobs minus # df absorbed in inference of ranefs
        ## eqs 4,7 in HaLM07
        info_crits$cAIC <- -2*forAIC$clik + 2*(pd+p_phi) ## no p_lambda !
      }
      # print(c(pd,p_phi))
      # but Yu and Yau then suggest caicc <- -2*clik + ... where ... involves d2h/db d[disp params] and d2h/d[disp params]2
    } else { ## fixed effect model
      info_crits$cAIC <- info_crits$mAIC <- -2*forAIC$p_v+2*(pforpv+p_phi) 
      # => sets cAIC so that (also_cAIC && is.null(info_crits[["cAIC"]])) becomes FALSE although also_cAIC is true by default...
    }
    object$envir$info_crits <- info_crits
  } else p_phi <- NULL
  if ( nsim>0L ) {
    if (is.null(p_phi)) p_phi <- .calc_p_phi(object)
    b_pd <- .calc_boot_AIC_dfs(object, nsim=nsim, ...)
    if (object$models[[1]]=="etaHGLM") {
      object$envir$info_crits$b_cAIC <- -2*object$APHLs$clik + 2*(b_pd+p_phi)
    } else object$envir$info_crits$b_cAIC <- -2*object$APHLs$p_v + 2*(b_pd+p_phi) # as documented
  } 
  return(object$envir$info_crits)
}

.calc_logdispObject <- function(object, envir=object$envir, force_fixed) {
  
  dvdloglamMat <- envir$dvdloglamMat
  dvdloglamMat_needed <- ( is.null(dvdloglamMat) && 
                             # (comment this => allows random slope)  all(unlist(attr(object$ZAlist,"namesTerms"))=="(Intercept)") && ## (1|.) or CAR or Matern
                             (force_fixed || 
                                any( ! object$lambda.object$type %in% c("fixed","fix_ranCoefs","fix_hyper"))) ) ## some lambda params were estimated
  dvdlogphiMat <- envir$dvdlogphiMat
  dvdlogphiMat_needed <- (is.null(dvdlogphiMat) && 
                           ( any((phimodel <- object$models[["phi"]])=="phiScal") ||
                               force_fixed)) 
  dvdlogphiMat_needed <- dvdlogphiMat_needed || identical(envir$forcePhiComponent,TRUE) ## hack for code testing !
  if (dvdloglamMat_needed || dvdlogphiMat_needed) {
    ZAL <- get_ZALMatrix(object, force_bind = ! (.is_spprec_fit(object)) )     
    d2hdv2_info <- .calc_d2hdv2_info(object, ZAL) # may be a qr object, or not (SPPREC). F I X M E a gentle message for long computations ? 
  } 
  if (dvdloglamMat_needed) { 
    cum_n_u_h <- attr(.get_u_h(object),"cum_n_u_h")
    psi_M <- rep(attr(object$rand.families,"unique.psi_M"),diff(cum_n_u_h))
    dlogfthdth <- (psi_M - .get_u_h(object))/object$lambda.object$lambda_est ## the d log density of th(u)
    neg.d2f_dv_dloglam <- .calc_neg_d2f_dv_dloglam(dlogfthdth, cum_n_u_h, 
                                                   lcrandfamfam=attr(object$rand.families,"lcrandfamfam"), 
                                                   rand.families=object$rand.families, u_h=.get_u_h(object))
    dvdloglamMat <- .calc_dvdloglamMat_new(neg.d2f_dv_dloglam,
                                           d2hdv2_info=d2hdv2_info) ## d2hdv2_info is either a qr factor or the inverse as a matrix or an environment
  }
  if (dvdlogphiMat_needed) {
    muetablob <- object$muetablob
    # .get_H_w.resid() rather than .get_w.resid here. See section mentioning "dvdloglamMat" in the long doc. 
    if ( ! is.null(envir$G_CHMfactor)) { # spprec; possibly generalisable code not using math-dense ZAL
      # rhs <- .Matrix_times_Dvec(t(envir$sXaug$AUGI0_ZX$ZAfix), - dh0deta) # efficient
      # rhs <- solve(envir$G_CHMfactor,rhs,system="A") # efficient
      unW_dh0deta <- (object$y-muetablob$mu)/muetablob$dmudeta ## (soit Bin -> phi fixe=1, soit BinomialDen=1)
      if (is.null(envir$invG_ZtW)) {
        if (is.null(envir$ZtW)) {
          # cf comments in .old_calc_Md2hdvb2_info_spprec_by_r22()
          envir$ZtW <- t(.Dvec_times_m_Matrix(.get_H_w.resid(envir=envir), envir$sXaug$AUGI0_ZX$ZAfix))
        }
        envir$invG_ZtW <- solve(envir$G_CHMfactor, 
                                envir$ZtW, system="A") # hardly avoidable has there is no comparable operation elsewhere (check "A")
      }
      rhs <- .Matrix_times_Dvec(envir$invG_ZtW, -unW_dh0deta) # efficient
      dvdlogphiMat <- .crossprod(envir$chol_Q, rhs) # _FIXME_ bottleneck in large spprec but .crossprodCpp_d not useful here 
    } else {
      dh0deta <- ( .get_H_w.resid(object) *(object$y-muetablob$mu)/muetablob$dmudeta ) ## (soit Bin -> phi fixe=1, soit BinomialDen=1) 
      dvdlogphiMat  <- .calc_dvdlogphiMat_new(dh0deta=dh0deta, ZAL=ZAL,
                                              d2hdv2_info=d2hdv2_info ## either a qr factor or a matrix inverse or envir
      )
    }
  }
  invV_factors <- .calc_invV_factors(object) ## n_x_r and r_x_n in repres of invV as diag(w.resid)- [n_x_r %*% r_x_n]
  if (length(phimodel)>1L) {
    .calc_logdisp_cov_mv(object, dvdloglamMat=dvdloglamMat, ## square matrix, by  the formulation of the algo 
                                                dvdlogphiMat=dvdlogphiMat, invV_factors=invV_factors,
                         force_fixed=force_fixed)
  } else
    .calc_logdisp_cov(object, dvdloglamMat=dvdloglamMat, ## square matrix, by  the formulation of the algo 
                                             dvdlogphiMat=dvdlogphiMat, invV_factors=invV_factors,
                      force_fixed=force_fixed)
  
}

.get_logdispObject <- function(object) { ## 
  envir <- object$envir
  if (is.null(envir$logdispObject) && object$models[["eta"]]=="etaHGLM" ) {
    envir$logdispObject <- .calc_logdispObject(object, envir=envir, force_fixed=FALSE)
  } 
  return(envir$logdispObject)
} # if dvdloglamMat or dvdlogphiMat were computed ex-tempo, they are NOT saved.

## This provides factor n_x_r and r_x_n of the representation of invV as diag(w.resid)- [n_x_r %*% r_x_n = t(Ztw) %*% invG.ZtW]  (nXr  %*% rxn)
## slow computation the one time .get_logdispObject() is called, for variances$disp (no need to store the result in an $envir)
.calc_invV_factors <- function(object) { ## used by .get_logdispObject
  ## Store inv( G=[{precmat=inv(L invWranef Lt)} +ZtWZ] ) as two matrix nXr and rXn rather than their nXn product
  # F I X M E yet there will be cases where n<r and then it's better to store the n x n product !  
  if (.is_spprec_fit(object)) {
    envir <- object$envir
    ## code clearly related to .Sigsolve_sparsePrecision() algorithm:
    if (is.null(envir$invG_ZtW)) {
      if (is.null(envir$ZtW)) {
        # cf comments in .old_calc_Md2hdvb2_info_spprec_by_r22()
        H_w.resid <- .get_H_w.resid(envir=envir)
        envir$ZtW <- t(.Dvec_times_m_Matrix(H_w.resid, envir$sXaug$AUGI0_ZX$ZAfix))
      }
      object$envir$invG_ZtW <- solve(envir$G_CHMfactor, envir$ZtW, system="A") # hardly avoidable has there is no comparable operation elsewhere (check "A")
    }    
    RES <- list(n_x_r=t(envir$ZtW), r_x_n=as.matrix(envir$invG_ZtW)) # requires both being kept in the envir 
    ZAfix <- .get_ZAfix(object, as_matrix=FALSE)
    RES$r_x_r <- object$envir$invG_ZtW %*% ZAfix
    return(RES)
  } else {
    ZAfix <- .get_ZAfix(object, as_matrix=FALSE)
    H_w.resid <- .get_H_w.resid(object)
    if (.is_identity(ZAfix)) {
      # FIXME inelegant ZAL computation only to test if it is diagonal
      ZAL <- .compute_ZAL(XMatrix=object$strucList, ZAlist=object$ZAlist,as_matrix=.eval_as_mat_arg.HLfit(object)) 
      if (inherits(ZAL,"diagonalMatrix")) { 
        RES <- list(n_x_r=Diagonal(x = H_w.resid),
                    r_x_n=Diagonal(x = H_w.resid/(H_w.resid + object$w.ranef/diag(ZAL)^2))) 
        return(RES)
      } ## ELSE
    } ## ELSE
    #
    invL <- .get_invL_HLfit(object) # t(tcrossfac(precision_matrix))
    wrZ <- .Dvec_times_m_Matrix(H_w.resid, ZAfix) # suppressMessages(sweep(t(ZA), 2L, w.resid,`*`)) 
    if (TRUE) { 
      ZtwrZ <- .crossprod(ZAfix, wrZ, as_sym=TRUE) ## seems more precise and we must compute wrZ anyway
    } else ZtwrZ <- .ZtWZwrapper(ZAfix,.get_H_w.resid(object)) ## -> calls .crossprod(., y=NULL) affects numerical precision of twolambda test in test-predVar.R
    ## is general ZtwrZ should be sparse, while invL may not. Large dense invL will away be a problem
    ## for small Z, ZtwrZ may bedsy, but this an un-intersting subcase that does not call for a special treatment
    if (inherits(ZtwrZ,"dsCMatrix") || inherits(ZtwrZ,"dsyMatrix")) {
      # (NB: dsy+dsy faster than dsC+dsy but maybe just because of conversion from dsC to dsy ?) => no obvious improvement
      if (.is_identity(invL)) {
        precmat <- .symDiagonal(x=object$w.ranef) # dsC+dsC
      } else precmat <- .ZtWZwrapper(invL,object$w.ranef) # dsC+whatever 
    } else if (inherits(ZtwrZ,"ddiMatrix")) {
      if (.is_identity(invL)) {
        precmat <- Diagonal(x=object$w.ranef) # ddi+ddi...
      } else precmat <- .ZtWZwrapper(invL,object$w.ranef) # ddi+whatever
    } else { # ("matrix") may never occur 
      if (.is_identity(invL)) {
        precmat <- diag(x=object$w.ranef)
      } else precmat <- .ZtWZwrapper(invL,object$w.ranef)
      message("Possibly inefficient code in .calc_invV_factors().")
    } 
    ## try to sum dsC or to sum dense matrix but not to mix types..., and to avoid formation of a large nxn matrix:
    if (inherits(precmat,"dsCMatrix") && inherits(ZtwrZ,"dsCMatrix")) {
      Gmat <- .dsCsum(precmat,ZtwrZ)
    } else Gmat <- precmat + ZtwrZ  
    invG_ZtW <- try(solve(Gmat, t(wrZ)),silent=TRUE)
    if (inherits(invG_ZtW,"try-error")) { ## but that should be well behaved when precmat is.
      invG <- ginv(as.matrix(Gmat)) ## FIXME quick patch at least
      invG_ZtW <- .tcrossprod(invG, wrZ)
    }  
    RES <- list(n_x_r=wrZ, r_x_n=invG_ZtW)
    #RES$r_x_r <- invG_ZtW %*% ZAfix  ## only for spprec (above) or  TRY_dense_iVZA (FALSE)
    return(RES)
  }
}

.calc_beta_cov_info_others <- function(wAugX=NULL, AUGI0_ZX, ZAL, ww) { ## post-fit fn
  if (is.null(wAugX)) {
    if (is.null(ZAL)) { # GLM... or LLM
      if (any(ww<0)) { # ... LLM
        negHess <- crossprod(AUGI0_ZX$X.pv, .Dvec_times_m_Matrix( ww, AUGI0_ZX$X.pv))
        cholH <- try(chol(negHess), silent=TRUE)
        if (inherits(cholH, "try-error")) { # => brute regularization
          warning("logLik presumably not maximized (information matrix is not positive definite at attained estimates).", immediate. = TRUE)
          beta_cov <- tcrossfac_beta_v_cov <- matrix(NA, ncol=ncol(negHess), nrow=ncol(negHess))
        } else {
          beta_cov <- chol2inv(cholH)
          tcrossfac_beta_v_cov <- t(cholH)
        }
        colnames(beta_cov) <- rownames(beta_cov) <- colnames(negHess) # sigh
        return(list(beta_cov=beta_cov, tcrossfac_beta_v_cov=tcrossfac_beta_v_cov))
      } else wAugX <- .calc_wAugX(XZ_0I=AUGI0_ZX$X.pv,sqrt.ww=sqrt(ww)) 
      # => ... may have zero cols... X.pv being a zero-col *m*atrix, in which case .spaMM.data$options$matrix_method is selected below 
      # and then get_from_MME(<matrix_method>, "beta_cov_info_from_wAugX") must handle this case.
    } else {
      XZ_0I <- .calc_XZ_0I(AUGI0_ZX=AUGI0_ZX,ZAL) # ZAL is not ZAXlist since .calc_beta_cov_info_others not called in spprec
      wAugX <- .calc_wAugX(XZ_0I=XZ_0I,sqrt.ww=sqrt(ww))
    }
  } ## wAugX is in XZ_OI order 
  if (inherits(wAugX,"Matrix")) {
    corr_method <- .spaMM.data$options$Matrix_method 
  } else corr_method <- .spaMM.data$options$matrix_method
  suppressMessages(try(untrace(corr_method, where=asNamespace("spaMM")),silent=TRUE)) # try() bc this fails when called by Infusion 
  # test: Infusion tests -> ... -> .predict_body -> ;.. -> .calc_beta_cov_info_others -> untrace -> def_sXaug_EigenDense_QRP_Chol_scaled not found
  # Note that the function is in exportPattern, explicitly exporting/impporting it does not help; attaching spaMM seems required to avoid untrace's error..
  # hack to recycle sXaug code; all weights are 1 or unit vectors as the order is not that assumed by corr_method. 
  wAugX <- do.call(corr_method,list(Xaug=wAugX, weight_X=rep(1,nrow(AUGI0_ZX$X.pv)), 
                                       w.ranef=rep(1,ncol(AUGI0_ZX$I)), ## we need at least its length for get_from Matrix methods
                                       H_global_scale=1))
  beta_cov_info <- get_from_MME(wAugX,"beta_cov_info_from_wAugX") 
  return(beta_cov_info)
}



.calc_beta_cov_info_from_sXaug <- function(BLOB, sXaug, tcrossfac) { 
  pforpv <- attr(sXaug,"pforpv")
  X_scaling <- sqrt(rep(attr(sXaug,"H_global_scale"),pforpv))
  X_scale <- attr(sXaug,"scaled:scale") # this is (! spprec) code and the X.pv attribute has been copied here
  if ( ! is.null(X_scale)) X_scaling <- X_scaling/X_scale
  diagscalings <- c(1/sqrt(attr(sXaug,"w.ranef")), X_scaling)
  if ( ! is.null(BLOB$sortPerm)) { # depending on method used for QR facto
    sctPmat <- sparseMatrix(seq_along(BLOB$sortPerm), BLOB$sortPerm, x=1)
    sctPmat <- .Dvec_times_Matrix(diagscalings, sctPmat) # Pmat <- .Matrix_times_Dvec(Pmat,diagscalings)
    tcrossfac_v_beta_cov <- sctPmat %*% tcrossfac
  } else tcrossfac_v_beta_cov <- .Dvec_times_m_Matrix(diagscalings, tcrossfac) ## loses colnames...
  dgC_good <- (inherits(tcrossfac_v_beta_cov, "dgCMatrix") && 
                 .calc_denseness(tcrossfac_v_beta_cov, relative=TRUE)<0.35 )
  if ( ! dgC_good) tcrossfac_v_beta_cov <- as.matrix(tcrossfac_v_beta_cov) # bigranefs.R shows that conversion is not always good.
  rownames(tcrossfac_v_beta_cov) <- colnames(sXaug) ## necessary for summary.HLfit, already lost in BLOB$R_scaled
  seqp <- seq_len(pforpv)
  beta_pos <- attr(sXaug,"n_u_h")+seqp
  beta_v_order <- c(beta_pos,seq(attr(sXaug,"n_u_h")))
  tcrossfac_beta_v_cov <- tcrossfac_v_beta_cov[beta_v_order,,drop=FALSE]
  if (inherits(sXaug,"dtCMatrix")) tcrossfac_beta_v_cov <- as(tcrossfac_beta_v_cov, "sparseMatrix")
  #beta_v_cov <- .tcrossprod(tcrossfac_beta_v_cov)
  beta_cov <- as.matrix(.tcrossprod(tcrossfac_beta_v_cov[seqp,,drop=FALSE])) ## assignment in .make_beta_table() assumes a dense matrix
  return( list(beta_cov=beta_cov, 
               #beta_v_cov=beta_v_cov,
               tcrossfac_beta_v_cov=tcrossfac_beta_v_cov) )
}



.get_beta_cov_info <- function(res) { 
  # Provide list(beta_cov=., tcrossfac_beta_v_cov=.)
  if (.is_spprec_fit(res)) {
    return(.calc_beta_cov_info_spprec(X.pv=res$X.pv, envir=res$envir)) 
  } else if (! is.null(res$envir$sXaug)) { # excludes SEM *and* fixed-effect models 
    if (prod(dim(res$envir$sXaug))>1e7) message("[one-time computation of covariance matrix, which may be slow]")
    if (res$spaMM.version<"2.7.34") attr(res$envir$sXaug,"scaled:scale") <- attr(res$X.pv,"scaled:scale")
    res$envir$beta_cov_info <- get_from_MME(res$envir$sXaug,which="beta_cov_info_from_sXaug") 
  } else { ## older code, generic for non-spprec.
    # What this [block -> .calc_beta_cov_info_others()] provides is a beta_cov + full beta_v_cov
    if (is.matrix(beta_cov_info <- res$envir$beta_cov_info) || ## old format (old fit object), should be a list now
        is.null(beta_cov_info$tcrossfac_beta_v_cov)) { ## Use .calc_beta_cov_info_others -> get_from_MME()
      # F I X M E test useless if call to .get_beta_cov_info() already dependent on the test
      # Further, can $tcrossfac_beta_v_cov be NULL if $beta_cov_info is list ?
      ZAL <- get_ZALMatrix(res, force_bind=FALSE) # note that .calc_beta_cov_info_others -> ... -> rbind2() but ZAL should not be a ZAXlist here
      nrd <- length(res$w.ranef)
      pforpv <- ncol(res$X.pv)
      if (inherits(ZAL,"Matrix")) {
        AUGI0_ZX <- list(I=.trDiagonal(n=nrd), # Ilarge=.trDiagonal(n=nrd+pforpv),
                         ZeroBlock=Matrix(0,nrow=nrd,ncol=pforpv),X.pv=res$X.pv)
      } else {
        AUGI0_ZX <- list(I=diag(nrow=nrd),ZeroBlock=matrix(0,nrow=nrd,ncol=pforpv),X.pv=res$X.pv)
      }
      res$envir$beta_cov_info <- .calc_beta_cov_info_others(AUGI0_ZX=AUGI0_ZX,ZAL=ZAL,ww=c(.get_H_w.resid(res),res$w.ranef)) ## with beta_v_cov attribute
    }
  }
  return(res$envir$beta_cov_info)
}

.calc_newZACvar <- function(newZAlist,cov_newLv_oldv_list) {
  newZACvarlist <- vector("list",length(newZAlist))
  for (new_rd in seq_along(newZAlist)) newZACvarlist[[new_rd]] <- newZAlist[[new_rd]] %ZA*gI% cov_newLv_oldv_list[[new_rd]]
  newZAC <- .ad_hoc_cbind(newZACvarlist, as_matrix=FALSE)
  if (inherits(newZAC,"sparseMatrix") && .calc_denseness(newZAC, relative=TRUE)>0.35) newZAC <- as.matrix(newZAC)
  newZAC
}

## Aggregate info on corrpars, inner-estimated and inner-fixed.
## $corrPars is only for info in messages() and return value, (?!)
.get_CorrEst_and_RanFix <- function(ranFix, ## has "fix", "outer", and also "var" values ! code corrently assumes "var" <=> corr_est
                                    corr_est 
) {
  # When fitting function was HLCor, the type attribute has not been added to the corrPars
  if (is.null(attr(ranFix,"type"))) attr(ranFix,"type") <- .relist_rep("fix", ranFix)
  if ( ! is.null(corr_est)) {
    ranFix <- structure(.modify_list(ranFix,corr_est), 
                        type=.modify_list(attr(ranFix,"type"),.relist_rep("var",corr_est)))
  } 
  return(ranFix) ## correlation params + whatever else was in ranFix
}

# Has become a postfit fn; to recycle processed$u_h_v_h_from_v_h one would have to remove processed from its defining envir.
.u_h_v_h_from_v_h <- function(v_h, rand.families, cum_n_u_h, lower.v_h, upper.v_h, 
                              u_list=vector("list", length(rand.families))) {
  if(!is.null(lower.v_h)) {v_h[v_h<lower.v_h] <- lower.v_h}
  if(!is.null(upper.v_h)) {v_h[v_h>upper.v_h] <- upper.v_h}
  nrand <- length(rand.families)
  for (it in seq_len(nrand)) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    u_list[[it]] <- rand.families[[it]]$linkinv(v_h[u.range])
    if (any(is.infinite(u_list[[it]]))) {
      warning("infinite random values ('u_h') were constrained to finite range.") 
      u_list[[it]] <- pmin(.Machine$double.xmax, pmax(-.Machine$double.xmax,u_list[[it]]) )
    }
  }
  u_h <- .unlist(u_list)
  ## if there were box constr, v_h may have been modified, we put it in return value
  if ( ! (is.null(lower.v_h) && is.null(upper.v_h))) attr(u_h,"v_h") <- v_h
  return(u_h)
}

