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
    pb <- txtProgressBar(max = nslices, style = 3, char="s")
    pd <- foreach::`%do%`(foreach_blob, {
      slice <- (slices[it]+1L):slices[it+1L]
      tmp <- t(.crossprod(tcrossfac_v_beta_cov, Md2clikdvb2[,slice]))
      setTxtProgressBar(pb, it)
      return(sum(tcrossfac_v_beta_cov[slice,] * tmp)) 
    })
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
  w.resid <- object$w.resid
  #ZAL <- .compute_ZAL(XMatrix=object$strucList, ZAlist=object$ZAlist,as_matrix=.eval_as_mat_arg(object)) 
  ZAL <- get_ZALMatrix(object) # allows ZAXlist; but if a non-ZAXlist is already in the $envir, no effect; + we will solve(chol_Q, Diagonal()) so the gain is not obvious
  if ( ncol(object$X.pv) ) {
    M12 <- as.matrix(.crossprod(ZAL, .Dvec_times_m_Matrix(w.resid, object$X.pv)))
    Md2clikdvb2 <- rbind2(cbind2(as.matrix(.ZtWZwrapper(ZAL,w.resid)), M12), ## this .ZtWZwrapper() takes time
                          cbind2(t(M12), as.matrix(.ZtWZwrapper(object$X.pv,w.resid)))) 
    # _FIXME_ any way to avoid formation of this matrix ? Or otherwise message() ?           
  } else {
    Md2clikdvb2 <-  as.matrix(.ZtWZwrapper(ZAL,w.resid))
  }
  tcrossfac_v_beta_cov <- .calc_Md2hdvb2_info_spprec(X.pv=object$X.pv, envir=object$envir, w.resid=w.resid, 
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
  
  w.resid <- object$w.resid
  #ZAL <- .compute_ZAL(XMatrix=obje17ct$strucList, ZAlist=object$ZAlist,as_matrix=.eval_as_mat_arg(object)) 
  ZAL <- get_ZALMatrix(object) # allows ZAXlist; but if a non-ZAXlist is already in the $envir, no effect; + we will solve(chol_Q, Diagonal()) so the gain is not obvious
  if ( ncol(object$X.pv) ) {
    M12 <- as.matrix(.crossprod(ZAL, .Dvec_times_m_Matrix(w.resid, object$X.pv)))
    Md2clikdvb2 <- rbind2(cbind2(as.matrix(.ZtWZwrapper(ZAL,w.resid)), M12), ## this .ZtWZwrapper() takes time
                          cbind2(t(M12), as.matrix(.ZtWZwrapper(object$X.pv,w.resid)))) 
    # _FIXME_ any way to avoid formation of this matrix ? Or otherwise message() ?           
  } else {
    Md2clikdvb2 <-  as.matrix(.ZtWZwrapper(ZAL,w.resid))
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



.get_info_crits <- function(object, also_cAIC=TRUE) {
  if (is.null(info_crits <- object$envir$info_crits) || (also_cAIC && is.null(info_crits$cAIC))) { 
    pforpv <- object$dfs[["pforpv"]]
    p_phi <- object$dfs[["p_fixef_phi"]]
    p_lambda <- object$dfs[["p_lambda"]]
    APHLs <- object$APHLs
    w.resid <- object$w.resid
    info_crits <- list()
    if  ( ! is.null(resid_fit <- object$resid_fit)) { ## fit includes a resid_fit
      # input p_phi (above) is typically set to NA, and will be ignored
      p_phi <- sum(resid_fit$dfs) ## phi_pd is relevant only for measuring quality of prediction by the resid_fit!
    } else p_phi <- object$dfs[["p_fixef_phi"]]
    names_est_ranefPars <- unlist(.get_methods_disp(object))  
    p_GLM_family <- length(intersect(names_est_ranefPars,c("NB_shape","NU_COMP")))
    p_phi <- p_phi+p_GLM_family ## effectively part of the model for residual error structure
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
      #corrPars <- get_ranPars(object,which="corrPars")
      p_corrPars <- object$dfs[["p_corrPars"]] # length(intersect(names_est_ranefPars,names(corrPars)))
      info_crits$mAIC <- -2*forAIC$p_v + 2 *(pforpv+p_lambda+p_corrPars+p_phi)
      info_crits$dAIC <- -2*forAIC$p_bv + 2 * (p_lambda+p_phi+p_corrPars) ## HaLM07 (eq 10) focussed for dispersion params
      #                                                                             including the rho param of an AR model
      if (also_cAIC) {
        if ( "AUGI0_ZX_sparsePrecision" %in% object$MME_method) {
          pd <- .calc_cAIC_pd_spprec(object)
        } else if ( ! is.null(object$envir$sXaug)) {
          pd <- .calc_cAIC_pd_from_sXaug(object)
        } else {
          ZAL <- .compute_ZAL(XMatrix=object$strucList, ZAlist=object$ZAlist,as_matrix=.eval_as_mat_arg(object)) 
          d2hdv2 <- .calcD2hDv2(ZAL,w.resid,object$w.ranef) ## update d2hdv2= - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
          pd <- .calc_cAIC_pd_others(X.pv, ZAL, w.resid, d2hdv2)
        }
        info_crits$GoFdf <- length(object$y) - pd ## <- nobs minus # df absorbed in inference of ranefs
        ## eqs 4,7 in HaLM07
        info_crits$cAIC <- -2*forAIC$clik + 2*(pd+p_phi) ## no p_lambda !
      }
      # print(c(pd,p_phi))
      # but Yu and Yau then suggest caicc <- -2*clik + ... where ... involves d2h/db d[disp params] and d2h/d[disp params]2
    } else { ## fixed effect model
      info_crits$mAIC <- -2*forAIC$p_v+2*(pforpv+p_phi) 
    }
    object$envir$info_crits <- info_crits
  } 
  return(info_crits)
}

.get_logdispObject <- function(object) { ## 
  if (is.null(object$envir$logdispObject) && object$models[["eta"]]=="etaHGLM" ) { 
    dvdloglamMat <- object$envir$dvdloglamMat
    dvdloglamMat_needed <- ( is.null(dvdloglamMat) && 
                               # (comment this => allows random slope)  all(unlist(attr(object$ZAlist,"namesTerms"))=="(Intercept)") && ## (1|.) or CAR or Matern
                               any( ! object$lambda.object$type %in% c("fixed","fix_ranCoefs","fix_hyper")) ) ## some lambda params were estimated
    dvdlogphiMat <- object$envir$dvdlogphiMat
    dvdlogphiMat_needed <- (is.null(dvdlogphiMat) && 
                              object$models[["phi"]]=="phiScal") ## cf comment in calc_logdisp_cov
    dvdlogphiMat_needed <- dvdlogphiMat_needed || identical(object$envir$forcePhiComponent,TRUE) ## hack for code testing !
    if (dvdloglamMat_needed || dvdlogphiMat_needed) {
      ZAL <- get_ZALMatrix(object)     
      d2hdv2_info <- .calc_d2hdv2_info(object, ZAL) # F I X M E a gentle message for long computations ? 
    } 
    if (dvdloglamMat_needed) { 
      cum_n_u_h <- attr(object$ranef,"cum_n_u_h")
      psi_M <- rep(attr(object$rand.families,"unique.psi_M"),diff(cum_n_u_h))
      dlogfthdth <- (psi_M - object$ranef)/object$lambda.object$lambda_est ## the d log density of th(u)
      neg.d2f_dv_dloglam <- .calc_neg_d2f_dv_dloglam(dlogfthdth, cum_n_u_h, 
                                                     lcrandfamfam=attr(object$rand.families,"lcrandfamfam"), 
                                                     rand.families=object$rand.families, u_h=object$ranef)
      dvdloglamMat <- .calc_dvdloglamMat_new(neg.d2f_dv_dloglam,
                                             d2hdv2_info=d2hdv2_info) ## d2hdv2_info is either a qr factor or the inverse as a matrix or an environment
    }
    if (dvdlogphiMat_needed) {
      muetablob <- object$muetablob
      if ( ! is.null(object$envir$G_CHMfactor)) { # possibly generalisable code not using math-dense ZAL
        # rhs <- .Matrix_times_Dvec(t(object$envir$ZAfix), - dh0deta) # efficient
        # rhs <- solve(object$envir$G_CHMfactor,rhs,system="A") # efficient
        unW_dh0deta <- (object$y-muetablob$mu)/muetablob$dmudeta ## (soit Bin -> phi fixe=1, soit BinomialDen=1)
        if (is.null(object$envir$invG_ZtW)) object$envir$invG_ZtW <- solve(object$envir$G_CHMfactor, 
                                                                           object$envir$ZtW, system="A") # hardly avoidable has there is no comparable operation elsewhere (check "A")
        rhs <- .Matrix_times_Dvec(object$envir$invG_ZtW, -unW_dh0deta) # efficient
        dvdlogphiMat <- .crossprod(object$envir$chol_Q, rhs) # _FIXME_ bottleneck in large spprec but .crossprodCpp not useful here 
      } else {
        dh0deta <- ( object$w.resid *(object$y-muetablob$mu)/muetablob$dmudeta ) ## (soit Bin -> phi fixe=1, soit BinomialDen=1)
        dvdlogphiMat  <- .calc_dvdlogphiMat_new(dh0deta=dh0deta, ZAL=ZAL,
                                                d2hdv2_info=d2hdv2_info, ## either a qr factor or a matrix inverse or envir
                                                stop.on.error=TRUE)
      }
    }
    invV_factors <- .calc_invV_factors(object) ## of invV as w.resid- [n_x_r %*% r_x_n]
    if ( (dvdloglamMat_needed || dvdlogphiMat_needed)  && inherits(ZAL,"ZAXlist")) { # I cannot test ZAL directly bc it has not necessarily been computed
      object$envir$logdispObject <- .calc_logdisp_cov_ZAX(object, dvdloglamMat=dvdloglamMat, ## square matrix, by  the formulation of the algo
                                                          dvdlogphiMat=dvdlogphiMat, invV_factors=invV_factors)
    } else
      object$envir$logdispObject <- .calc_logdisp_cov(object, dvdloglamMat=dvdloglamMat, ## square matrix, by  the formulation of the algo 
                                                      dvdlogphiMat=dvdlogphiMat, invV_factors=invV_factors)
  } 
  return(object$envir$logdispObject)
} # if dvdloglamMat or dvdlogphiMat were computed ex-tempo, they are NOT saved.

## This use the representation of invV as w.resid- [n_x_r %*% r_x_n = t(Ztw) %*% invG.ZtW]  (nXr  %*% rxn)
## slow computation the one time .get_logdispObject() is called, for variances$disp (no need to store the result in an $envir)
.calc_invV_factors <- function(object) { ## used by .get_logdispObject
  ## Store Compute inv( G=[{precmat=inv(L invWranef Lt)} +ZtWZ] ) as two matrix nXr and rXn rather than their nXn product
  # F I X M E yet there will be cases where n<r and then it's better to store the n x n product !  
  if ("AUGI0_ZX_sparsePrecision" %in% object$MME_method) {
    ## code clearly related to .Sigsolve_sparsePrecision() algorithm:
    if (is.null(object$envir$invG_ZtW)) object$envir$invG_ZtW <- solve(object$envir$G_CHMfactor, object$envir$ZtW, system="A") # hardly avoidable has there is no comparable operation elsewhere (check "A")
    RES <- list(n_x_r=t(object$envir$ZtW),
                r_x_n=as.matrix(object$envir$invG_ZtW)) 
    # if as(,"dgeMatrix"), .calc_denseness() -> drop0() wastes time.
    ZAfix <- .get_ZAfix(object, as_matrix=FALSE)
    RES$r_x_r <- object$envir$invG_ZtW %*% ZAfix
    return(RES)
  } else {
    # FIXME inelegant ZAL computation only to test if it is diagonal
    ZAL <- .compute_ZAL(XMatrix=object$strucList, ZAlist=object$ZAlist,as_matrix=.eval_as_mat_arg(object)) 
    if (inherits(ZAL,"diagonalMatrix")) { # direct but ad hoc
      w.resid <- object$w.resid
      RES <- list(n_x_r=Diagonal(x = w.resid),
                  r_x_n=Diagonal(x = w.resid/(w.resid + object$w.ranef/diag(ZAL)^2))) 
    } else {
      ZAfix <- .get_ZAfix(object, as_matrix=FALSE)
      wrZ <- .Dvec_times_m_Matrix(object$w.resid, ZAfix) # suppressMessages(sweep(t(ZA), 2L, w.resid,`*`)) 
      invL <- .get_invL_HLfit(object) # equivalent to t(chol_Q)
      if (.is_identity(invL)) {
        precmat <- diag(x=object$w.ranef)
      } else precmat <- .ZtWZwrapper(invL,object$w.ranef)
      ZtwrZ <- .crossprod(ZAfix, wrZ) 
      ## avoid formation of a large nxn matrix:
      Gmat <- precmat+as.matrix(ZtwrZ)
      invG_ZtW <- try(solve(Gmat, t(wrZ)),silent=TRUE)
      if (inherits(invG_ZtW,"try-error")) { ## but that should be well behaved when precmat is.
        invG <- ginv(as.matrix(Gmat)) ## FIXME quick patch at least
        invG_ZtW <- .tcrossprod(invG, wrZ)
      }  
      RES <- list(n_x_r=wrZ, r_x_n=invG_ZtW)
    } 
    ZAfix <- .get_ZAfix(object, as_matrix=FALSE)
    RES$r_x_r <- invG_ZtW %*% ZAfix
    return(RES)
  }
}

.calc_beta_cov_info_others <- function(wAugX=NULL, AUGI0_ZX, ZAL, ww) { ## post-fit fn
  if (is.null(wAugX)) {
    if (is.null(ZAL)) {
      wAugX <- .calc_wAugX(XZ_0I=AUGI0_ZX$X.pv,sqrt.ww=sqrt(ww))
    } else {
      XZ_0I <- .calc_XZ_0I(AUGI0_ZX=AUGI0_ZX,ZAL) # ZAL is not ZAXlist since .calc_beta_cov_info_others not called in spprec
      wAugX <- .calc_wAugX(XZ_0I=XZ_0I,sqrt.ww=sqrt(ww))
    }
  } ## wAugX is in XZ_OI order 
  if (inherits(wAugX,"Matrix")) {
    mMatrix_method <- .spaMM.data$options$Matrix_method 
  } else mMatrix_method <- .spaMM.data$options$matrix_method
  suppressMessages(try(untrace(mMatrix_method, where=asNamespace("spaMM")),silent=TRUE)) # try() bc this fails when called by Infusion 
  # test: Infusion tests -> ... -> .predict_body -> ;.. -> .calc_beta_cov_info_others -> untrace -> def_sXaug_EigenDense_QRP_Chol_scaled not found
  # Note that the function is in exportPattern, explicitly exporting/impporting it does not help; attaching spaMM seems required to avoid untrace's error..
  # hack to recycle sXaug code; all weights are 1 or unit vectors as the order is not that assumed by mMatrix_method. 
  wAugX <- do.call(mMatrix_method,list(Xaug=wAugX, weight_X=rep(1,nrow(AUGI0_ZX$X.pv)), 
                                       w.ranef=rep(1,ncol(AUGI0_ZX$I)), ## we need at least its length for get_from Matrix methods
                                       H_global_scale=1))
  beta_cov_info <- get_from_MME(wAugX,"beta_cov_info_from_wAugX") 
  return(beta_cov_info)
}



.calc_beta_cov_info_from_sXaug <- function(BLOB, sXaug, B) {
  tcrossfac_v_beta_cov <-  solve(BLOB$R_scaled) # solve(as(BLOB$R_scaled,"dtCMatrix"))
  pforpv <- attr(sXaug,"pforpv")
  X_scaling <- sqrt(rep(attr(sXaug,"H_global_scale"),pforpv))
  if ( ! is.null(B)) X_scaling <- X_scaling/B
  diagscalings <- c(1/sqrt(attr(sXaug,"w.ranef")), X_scaling)
  if ( ! is.null(BLOB$sortPerm)) { # depending on method used for QR facto
    tPmat <- sparseMatrix(seq_along(BLOB$sortPerm), BLOB$sortPerm, x=1)
    tPmat <- .Dvec_times_Matrix(diagscalings, tPmat) # Pmat <- .Matrix_times_Dvec(Pmat,diagscalings)
    tcrossfac_v_beta_cov <- tPmat %*% tcrossfac_v_beta_cov
  } else tcrossfac_v_beta_cov <- .Dvec_times_m_Matrix(diagscalings, tcrossfac_v_beta_cov) ## loses colnames...
  dgC_good <- (inherits(tcrossfac_v_beta_cov, "dgCMatrix") && 
                 .calc_denseness(tcrossfac_v_beta_cov)/prod(dim(tcrossfac_v_beta_cov))<0.35 )
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
  if ( "AUGI0_ZX_sparsePrecision" %in% res$MME_method) {
    return(.calc_beta_cov_info_spprec(X.pv=res$X.pv, envir=res$envir, w.resid=res$w.resid)) 
  } else if (! is.null(res$envir$sXaug)) { # excludes SEM *and* fixed-effect models 
    if (prod(dim(res$envir$sXaug))>1e7) message("[one-time computation of covariance matrix, which may be slow]")
    res$envir$beta_cov_info <- get_from_MME(res$envir$sXaug,which="beta_cov_info_from_sXaug", 
                                  B=attr(res$envir$sXaug,"scaled:scale")) 
  } else { ## older code, generic for non-spprec.
    # What this [block -> .calc_beta_cov_info_others()] provides is a beta_cov + full beta_v_cov
    if (is.matrix(beta_cov_info <- res$envir$beta_cov_info) || ## old format (old fit object), should be a list now
        is.null(beta_cov_info$tcrossfac_beta_v_cov)) { ## Use .calc_beta_cov_info_others -> get_from_MME()
      # F I X M E test useless if call to .get_beta_cov_info() already dependent on the test
      # Further, can $tcrossfac_beta_v_cov be NULL if $beta_cov_info is list ?
      ZAL <- get_ZALMatrix(res) # note that .calc_beta_cov_info_others -> ... -> rbind2() but ZAL should not be a ZAXlist here
      nrd <- length(res$w.ranef)
      pforpv <- ncol(res$X.pv)
      if (inherits(ZAL,"Matrix")) {
        AUGI0_ZX <- list(I=.trDiagonal(n=nrd),
                         ZeroBlock=Matrix(0,nrow=nrd,ncol=pforpv),X.pv=res$X.pv)
      } else {
        AUGI0_ZX <- list(I=diag(nrow=nrd),ZeroBlock=matrix(0,nrow=nrd,ncol=pforpv),X.pv=res$X.pv)
      }
      res$envir$beta_cov_info <- .calc_beta_cov_info_others(AUGI0_ZX=AUGI0_ZX,ZAL=ZAL,ww=c(res$w.resid,res$w.ranef)) ## with beta_v_cov attribute
    }
  }
  return(res$envir$beta_cov_info)
}

.calc_newZACvar <- function(newZAlist,cov_newLv_oldv_list) {
  newZACvarlist <- vector("list",length(newZAlist))
  for (new_rd in seq_along(newZAlist)) {
    terme <- newZAlist[[new_rd]] %ZA*gI% cov_newLv_oldv_list[[new_rd]]
    terme <- as.matrix(terme)
    newZACvarlist[[new_rd]] <- as.matrix(terme)
  }
  return(do.call(cbind,newZACvarlist))
}

## Aggregate info on corrpars, inner-estimated and inner-fixed.
## $corrPars is only for info in messages() and return value, (?!)
.get_CorrEst_and_RanFix <- function(ranFix, ## has "fix", "outer", and also "var" values ! code corrently assumes "var" <=> corr_est
                                    corr_est 
) {
  if ( ! is.null(corr_est)) {
    ranFix <- structure(.modify_list(ranFix,corr_est), 
                        type=.modify_list(attr(ranFix,"type"),relist(rep("var",length(unlist(corr_est))),corr_est)))
  } 
  return(ranFix) ## correlation params + whatever else was in ranFix
}

