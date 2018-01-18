.do_damped_WLS <- function(sXaug, zInfo, # fixed descriptors of system to be solved
                         old_Vscaled_beta,
                         oldAPHLs,
                         APHLs_args,
                         damping,
                         dampingfactor=2, ## no need to change the input value
                         ypos,off,GLMMbool,etaFix,constant_u_h_v_h_args,
                         updateW_ranefS_constant_arglist,
                         wranefblob,seq_n_u_h,ZAL_scaling,
                         processed, trace=FALSE, Xscal,mMatrix_method,
                         phi_est, H_global_scale, n_u_h, ZAL, which_i_affected_rows,
                         which_LevMar_step
) {
  newXscal <- Xscal ## template
  initdamping <- damping
  gainratio_grad <- zInfo$gainratio_grad
  # grad wrt scaled v = d f / d (v/ZAL_scaling) = ZAL_scaling * d f / d v
  use_heuristic <- TRUE
  while ( TRUE ) {
    if (processed$HL[1L]==1L) { ## ML fit 
      Vscaled_beta <- old_Vscaled_beta
      ## maximize p_v wrt beta only
      if ( which_LevMar_step=="v_b") { ## note tests on summand too !!!!!!!
        LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step", LMrhs=zInfo$scaled_grad, damping=damping)
        Vscaled_beta <- Vscaled_beta + LevMarblob$dVscaled_beta 
      } else if ( which_LevMar_step=="b") {          
        LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step_beta", LMrhs=zInfo$scaled_grad[-seq_n_u_h], damping=damping)
        Vscaled_beta[-seq_n_u_h] <- Vscaled_beta[-seq_n_u_h] + LevMarblob$dbeta 
      } else if ( which_LevMar_step=="v") { ## v_h estimation given beta 
        LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step_v_h", LMrhs=zInfo$scaled_grad[seq_n_u_h], damping=damping)
        Vscaled_beta[seq_n_u_h] <- Vscaled_beta[seq_n_u_h] + LevMarblob$dVscaled 
      }
    } else { ## joint hlik maximization
      LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step", LMrhs=zInfo$scaled_grad, damping=damping)
      Vscaled_beta <- old_Vscaled_beta + LevMarblob$dVscaled_beta 
    }
    fitted <- drop(Xscal %*% Vscaled_beta) ## length nobs+nr ! 
    fitted[ypos] <- eta <- fitted[ypos] + off
    
    newmuetablob <- .muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed) 
    neww.resid <- .calc_w_resid(newmuetablob$GLMweights,phi_est)
    newweight_X <- .calc_weight_X(neww.resid, H_global_scale) ## sqrt(s^2 W.resid)
    
    if (is.null(etaFix$v_h)) { 
      v_h <- Vscaled_beta[seq_n_u_h] * ZAL_scaling ## use original scaling!
      if (GLMMbool) {
        u_h <- v_h 
        newwranefblob <- wranefblob ## keep input wranefblob since GLMM and lambda_est not changed
      } else {
        u_h_v_h_from_v_h_args <- c(constant_u_h_v_h_args,list(v_h=v_h))
        u_h <- do.call(".u_h_v_h_from_v_h",u_h_v_h_from_v_h_args)
        if ( ! is.null(attr(u_h,"v_h"))) { ## second test = if u_h_info$upper.v_h or $lower.v_h non NULL
          v_h <- attr(u_h,"v_h")
        }
        ## update functions u_h,v_h
        newwranefblob <- do.call(".updateW_ranefS",c(updateW_ranefS_constant_arglist,list(u_h=u_h,v_h=v_h)))
      } 
    } else newwranefblob <- wranefblob
    
    mMatrix_arglist <- list(weight_X=newweight_X, w.ranef=newwranefblob$w.ranef, H_global_scale=H_global_scale)
    if ( ! GLMMbool ) { 
      # newZAL_scaling necessary to get the correct logdet_sqrt_d2hdv2 for newsXaug
      newZAL_scaling <- 1/sqrt(newwranefblob$w.ranef*H_global_scale) ## Q^{-1/2}/s
      ## used only to compute a likelihood, not to update a system to be solved.
      if (TRUE) {
        ## alternative clause shows the meaning. This is distinctly faster. 
        scaledZAL <- .m_Matrix_times_Dvec(ZAL, newZAL_scaling)
        if (inherits(ZAL,"Matrix")) {
          # Next line assumes Xscal (IO_ZX) has no @diag component
          newXscal@x[which_i_affected_rows] <- scaledZAL@x ## should create an error if some elements are stored in @diag
        } else {
          newXscal[n_u_h+seq(nrow(scaledZAL)),seq_n_u_h] <- scaledZAL
        }
      } else newXscal <- .make_Xscal(ZAL, ZAL_scaling = newZAL_scaling, AUGI0_ZX=processed$AUGI0_ZX, n_u_h=n_u_h)
      mMatrix_arglist$Xaug <- newXscal 
    } else mMatrix_arglist$Xaug <- Xscal ##not distinct from the 'resident' Xscal
    newsXaug <- do.call(mMatrix_method, mMatrix_arglist)
    
    APHLs_args$sXaug <- newsXaug
    APHLs_args$dvdu <- newwranefblob$dvdu
    APHLs_args$u_h <- u_h 
    APHLs_args$mu <- newmuetablob$mu
    newAPHLs <- do.call(".calc_APHLs_from_ZX", APHLs_args)
    if (processed$HL[1L]==1L) { 
      if (which_LevMar_step=="v") { 
        newlik <- newAPHLs[["hlik"]]
        oldlik <- oldAPHLs[["hlik"]]
        if (damping==0L) {
          if (trace) print(paste("IRLS step for v_h, hlik=",newlik))
          break
        }
      } else {
        newlik <- newAPHLs[["p_v"]]
        oldlik <- oldAPHLs[["p_v"]]
        if (damping==0L) {
          if (trace) print(paste("IRLS step for (beta,v_h); p_v=",newlik))
          break
        }
      } 
    } else {
      newlik <- newAPHLs[["hlik"]]
      oldlik <- oldAPHLs[["hlik"]]
      if (damping==0L) {
        if (trace) print(paste("IRLS step, hlik=",newlik))
        break
      }
    }
    ## ELSE
    gainratio <- (newlik!=-Inf) ## -Inf occurred in binary probit with extreme eta... 
    if (gainratio) { 
      if (processed$HL[1L]==1L) { ## ML fit 
        if (which_LevMar_step=="v_b") { 
          tempodvhbeta <- LevMarblob$dVscaled_beta
          tempodvhbeta[seq_n_u_h] <- tempodvhbeta[seq_n_u_h]*ZAL_scaling
          summand <- tempodvhbeta*(gainratio_grad+ LevMarblob$dampDpD * tempodvhbeta)
        } else if (which_LevMar_step=="b") {
          summand <- LevMarblob$dbeta*(gainratio_grad[-seq_n_u_h]+ LevMarblob$dampDpD * LevMarblob$dbeta) 
        } else if (which_LevMar_step=="v") { ## v_h estimation given beta (FIXME can surely be made more exact)
          tempodvh <- LevMarblob$dVscaled*ZAL_scaling
          summand <- tempodvh*(gainratio_grad[seq_n_u_h]+ LevMarblob$dampDpD * tempodvh) 
        }
      } else { ## joint hlik maximization
        summand <- LevMarblob$dVscaled_beta*(gainratio_grad+ LevMarblob$dampDpD * LevMarblob$dVscaled_beta) 
      }
      ## The two terms of the summand should be positive. In part. conv_dbetaV*rhs should be positive. 
      ## However, numerical error may lead to <0 or even -Inf
      ## Further, if there are both -Inf and +Inf elements the sum is NaN.
      summand[summand<0] <- 0
      denomGainratio <- sum(summand)
      dlogL <- newlik-oldlik
      conv_logL <- abs(dlogL)/(1+abs(newlik))
      gainratio <- 2*dlogL/denomGainratio ## cf computation in MadsenNT04 below 3.14, but generalized for D' != I ; numerator is reversed for maximization
    } else { ## 2017/10/16 a patch to prevent a stop in this case, but covers up dubious computations (FIXME)
      newlik <- -.Machine$double.xmax
      dlogL <- newlik-oldlik
      conv_logL <- abs(dlogL)/(1+abs(newlik))
    }
    if (trace) print(paste("dampingfactor=",dampingfactor,#"innerj=",innerj,
                           "damping=",damping,"gainratio=",gainratio,"oldlik=",oldlik,"newlik=",newlik))
    if (is.nan(gainratio)) {
      # if the initial logL is the solution logL, then damping diverges 
      # it is then possible that some element of dVscaled_beta =0 and some of dampDpD =Inf
      # then the summand has some NaN
      # At the same time not all elements of dVscaled_beta need be 0 (eg small changes in eta for mu=0 or 1 in binomial models)
      # so testing dVscaled_beta is not sufficient to stop the algo
      # (LevenbergM is quite slow in such cases)
      break
    }
    if (gainratio > 0) { 
      ## cf Madsen-Nielsen-Tingleff again, and as in levmar library by Lourakis
      damping <- damping * max(1/3,1-(2*gainratio-1)^3)  
      # damping <- min(1000,damping )
      break 
    } else if (dampingfactor>4 ## implies iter>2
               && conv_logL <1e-8 && abs(prev_conv_logL) <1e-8) {
      damping <- initdamping ## bc presumably damping has diverged unusefully
      break ## apparently flat likelihood; this has occurred when fit_as_ZX used a wrong initial (beta,v_h), but still occurs in /tests
    } else { ## UNsuccessful step
      prev_conv_logL <- conv_logL
      damping <- damping*dampingfactor
      dampingfactor <- dampingfactor*2
    }
    if (damping>1e100) stop("reached damping=1e100")
  } 
  RESU <- list(lik=newlik,APHLs=newAPHLs,damping=damping,sXaug=newsXaug,
               fitted=fitted, eta=eta, muetablob=newmuetablob, wranefblob=newwranefblob, 
               v_h=v_h,u_h=u_h, w.resid=neww.resid)
  if ( ! GLMMbool ) {
    RESU$ZAL_scaling <- newZAL_scaling
    RESU$Xscal <- newXscal
    Vscaled_beta[seq_n_u_h] <- v_h/newZAL_scaling ## represent solution in new scaling...
  } 
  RESU$Vscaled_beta <- Vscaled_beta 
  return(RESU)
  
}

.solve_IRLS_as_ZX <- function(X.pv, 
                      ZAL, 
                      y, ## could be taken fom processed ? 
                      n_u_h, H_global_scale, lambda_est, muetablob=NULL ,off, maxit.mean, etaFix,
                      wranefblob, processed,
                      ## supplement for ! LMM
                      phi_est, 
                      ## supplement for for LevenbergM or ! GLMM
                      eta=NULL, w.resid=NULL, 
                      ## supplement for ! GLMM
                      u_h, v_h, for_init_z_args, 
                      ## supplement for LevenbergM
                      beta_eta,
                      ## supplement for intervals
                      for_intervals,
                      trace=FALSE
) {
  pforpv <- ncol(X.pv)
  nobs <- length(y)
  seq_n_u_h <- seq_len(n_u_h)
  ypos <- n_u_h+seq_len(nobs)
  lcrandfamfam <- attr(processed$rand.families,"lcrandfamfam")
  LMMbool <- processed$LMMbool
  GLMMbool <- processed$GLMMbool
  LevenbergM <- (.determine_LevenbergM(processed$LevenbergM) && is.null(for_intervals))
  is_HL1_1 <- (processed$HL[1L]==1L)
  which_LevMar_step <- "v_b"
  old_relV_beta <- NULL
  not_moving <- FALSE
  damped_WLS_blob <- NULL
  Ftol_LM <- processed$spaMM_tol$Ftol_LM
  if ( LevenbergM) { 
    damping <- 1e-7
    loc_Xtol_rel <- 1e-03 ## maybe good compromise between optim accuracy and time. 
  } else {
    damping <- 0L ## indicator for early exit from .do_damped_WLS without a check for increase
    loc_Xtol_rel <- processed$spaMM_tol$Xtol_rel/10
  }  
  if ( ! LMMbool) {
    constant_zAug_args <- list(n_u_h=n_u_h, nobs=nobs, pforpv=pforpv, y=y, off=off, ZAL=ZAL, processed=processed)
    if ( ! GLMMbool) {
      constant_init_z_args <- c(list(lcrandfamfam=lcrandfamfam, nobs=nobs, lambda_est=lambda_est, ZAL=ZAL),  
                                # fit_as_ZX args specific for ! GLMM:
                                for_init_z_args,
                                #
                                mget(c("cum_n_u_h","rand.families","stop.on.error"),envir=processed))
      constant_u_h_v_h_args <- c(mget(c("cum_n_u_h","rand.families"),envir=processed),
                                 processed$u_h_info, ## elements of u_h_info as elements of constant_u_h_v_h_args  
                                 list(lcrandfamfam=lcrandfamfam))
      updateW_ranefS_constant_arglist <- c(mget(c("cum_n_u_h","rand.families"),envir=processed),list(lambda=lambda_est))
    } 
  } 
  
  ##### initial sXaug
  ZAL_scaling <- 1/sqrt(wranefblob$w.ranef*H_global_scale) ## Q^{-1/2}/s
  Xscal <- .make_Xscal(ZAL, ZAL_scaling = ZAL_scaling,
                       AUGI0_ZX=processed$AUGI0_ZX,n_u_h=n_u_h)
  if (inherits(Xscal,"Matrix")) { # same type as ZAL
    ## def_sXaug_Eigen_sparse_QR calls lmwith_sparse_QRp(,pivot=FALSE)
    ## def_sXaug_Eigen_sparse_QRP calls lmwith_sparse_QRp(,pivot=TRUE)
    mMatrix_method <- .spaMM.data$options$Matrix_method
    
    #@p[c] must contain the index _in @x_ of the first nonzero element of column c, x[p[c]] in col c and row i[p[c]])  
    elmts_affected_cols <- seq_len(Xscal@p[n_u_h+1L]) ## corresponds to cols seq_n_u_h
    which_i_affected_rows <- which(Xscal@i[elmts_affected_cols]>(n_u_h-1L))    
  } else {
    mMatrix_method <- .spaMM.data$options$matrix_method
    which_i_affected_rows <- NULL
  }
  if ( ! is.null(for_intervals)) {
    Vscaled_beta <- c(v_h/ZAL_scaling ,for_intervals$beta_eta)
  } else if (LevenbergM) {
    Vscaled_beta <- c(v_h/ZAL_scaling ,beta_eta)
  } # else is NULL 
  if (is.null(eta)) { ## NULL input eta allows NULL input muetablob
    eta  <- off + (Xscal %*% c(v_h/ZAL_scaling ,beta_eta))[ypos]
    muetablob <- .muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed) 
  }
  ## weight_X and Xscal varies within loop if ! LMM since at least the GLMweights in w.resid change
  if ( is.null(w.resid) ) w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est)
  weight_X <- .calc_weight_X(w.resid, H_global_scale) ## sqrt(s^2 W.resid)
  sXaug <- do.call(mMatrix_method,
                   list(Xaug=Xscal, weight_X=weight_X, w.ranef=wranefblob$w.ranef, H_global_scale=H_global_scale))
  allow_LM_restart <- ( ! LMMbool && ! LevenbergM && is.null(for_intervals) && is.na(processed$LevenbergM["user_LM"]) )
  if (allow_LM_restart) {
    keep_init <- new.env()
    #names_keep <- ls()  
    names_keep <- c("sXaug","wranefblob","muetablob","u_h","w.resid","eta","v_h","ZAL_scaling","weight_X","Xscal","beta_eta")
    for (st in names_keep) keep_init[[st]] <- environment()[[st]]
  }
  LMcond <- - 10. 
  ################ L O O P ##############
  for (innerj in 1:maxit.mean) {
    if( ! LevenbergM && allow_LM_restart) { ## FIXME the next step improvement would be 
      #  ./. to keep track of lowest lambda that created problem and use LM by default then
      if (innerj>3) {
        LMcond <- LMcond + mean(abs_d_relV_beta/(old_abs_d_relV_beta+1e-8))
        ## cat(mean(abs_d_relV_beta/old_abs_d_relV_beta)," ")
        # cat(LMcond/innerj," ")
        if (LMcond/innerj>0.5) {
          if (trace) cat("!LM")
          for (st in names_keep) assign(st,keep_init[[st]])
          LevenbergM <- TRUE
          Vscaled_beta <- c(v_h/ZAL_scaling ,beta_eta) ## bc initialized only | LevenbergM 
          damping <- 1e-7
          loc_Xtol_rel <- 1e-03 ## maybe good compromise between optim accuracy and time. 
          damped_WLS_blob <- NULL
          allow_LM_restart <- FALSE 
        }
      }
      if (innerj>2) old_abs_d_relV_beta <- abs_d_relV_beta
    }
    ##### get the lik of the current state
    if ( ! is.null(for_intervals)) {
      loc_logLik_args <- list(sXaug=sXaug, processed=processed, phi_est=phi_est,
                              lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, mu=muetablob$mu)
      oldlik <- unlist(do.call(".calc_APHLs_from_ZX",loc_logLik_args)[for_intervals$likfn]) # unlist keeps name
    } else {
      if (LevenbergM) {
        if (is.null(damped_WLS_blob)) {
          oldAPHLs <- .calc_APHLs_from_ZX(sXaug=sXaug, processed=processed, phi_est=phi_est, which=processed$p_v_obj, 
                                          lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, mu=muetablob$mu)
        } else { ## Levenberg and innerj>1
          oldAPHLs <- damped_WLS_blob$APHLs
        }
      } 
    }
    #####
    
    ##### RHS
    if (LMMbool) {
      wzAug <- c(rep(0,n_u_h),(y-off)*weight_X)
    } else {
      if ( ! GLMMbool) {
        # arguments for init_resp_z_corrections_new called in calc_zAug_not_LMM
        init_z_args <- c(constant_init_z_args,
                         list(w.ranef=wranefblob$w.ranef, u_h=u_h, v_h=v_h, dvdu=wranefblob$dvdu, 
                              sXaug=sXaug, w.resid=w.resid))
      } else init_z_args <- NULL
      calc_zAug_args <- c(constant_zAug_args,
                          list(eta=eta, muetablob=muetablob, dlogWran_dv_h=wranefblob$dlogWran_dv_h, 
                               sXaug=sXaug, ## .get_qr may be called for Pdiag calculation
                               w.ranef=wranefblob$w.ranef, 
                               w.resid=w.resid,
                               init_z_args=init_z_args) )
      zInfo <- do.call(".calc_zAug_not_LMM",calc_zAug_args) ## dlogfvdv is 'contained' in $z2
      wzAug <- c(zInfo$y2_sscaled/ZAL_scaling, (zInfo$z1_sscaled)*weight_X) 
    }
    ## keep name 'w'zAug to emphasize the distinct weightings  of zaug and Xaug (should have been so everywhere)
    #####
    ##### improved  Vscaled_beta   
    if ( ! is.null(for_intervals)) {
      currentDy <- (for_intervals$fitlik-oldlik)
      if (currentDy < -1e-4) .warn_intervalStep(oldlik,for_intervals)
      intervalBlob <- .intervalStep_ZX(old_Vscaled_beta=Vscaled_beta,
                                       sXaug=sXaug,szAug=wzAug,
                                       for_intervals=for_intervals,
                                       currentlik=oldlik,currentDy=currentDy)
      damped_WLS_blob <- NULL
      Vscaled_beta <- intervalBlob$Vscaled_beta
    } else if (LevenbergM) { ## excludes IRLS
      ## (w)zAug is all what is needed for the direct solution of the extended system. in GLMM case
      # Hence wZaug contains Phi z_2 including (Phi v^0 +dlogfvdv)/ZAL_scaling (from components of hlik)
      ## now we want the LHS of a d_beta_v solution
      etamo <- eta - off
      zInfo$z1_eta <- zInfo$z1-etamo
      z1_sscaled_eta <- zInfo$z1_sscaled - etamo # zAug[-seq_n_u_h]-etamo # z_1-sscaled-etamo
      if (GLMMbool) {
        zInfo$dlogfvdv <-  - v_h * wranefblob$w.ranef
      } else zInfo$dlogfvdv <- (zInfo$z2 - v_h) * wranefblob$w.ranef
      ## the gradient for -p_v (or -h), independent of the scaling
      if (is.list(w.resid)) {
        m_grad_obj <- c( ## drop() avoids c(Matrix..) 
          m_grad_v <- drop(.crossprod(ZAL, w.resid$w_resid * zInfo$z1_eta) + zInfo$dlogfvdv), # Z'W(z_1-eta)+ dlogfvdv
          drop(.crossprod(X.pv, w.resid$w_resid * z1_sscaled_eta)) # X'W(z_1-sscaled-eta)
        )
      } else {
        m_grad_obj <- c( ## drop() avoids c(Matrix..) 
          m_grad_v <- drop(.crossprod(ZAL, w.resid * zInfo$z1_eta) + zInfo$dlogfvdv), # Z'W(z_1-eta)+ dlogfvdv
          drop(.crossprod(X.pv, w.resid * z1_sscaled_eta)) # X'W(z_1-sscaled-eta)
        )
      }
      if (not_moving && is_HL1_1) { ## not_moving TRUE may occur when we are out of solution space. Hence test Mg_solve_g
        Mg_solve_g <- get_from_MME(sXaug=sXaug, which="Mg_solve_g", B=m_grad_obj)
        if (Mg_solve_g < Ftol_LM/2) {
          if (trace>1L) {"break bc Mg_solve_g<1e-6"}
          break
        }
      } ## else not_moving was a break condition elsewhere in code
      zInfo$gainratio_grad <- m_grad_obj ## before recaling
      # gradient for scaled system from gradient of objective
      scaled_grad <- H_global_scale * m_grad_obj
      scaled_grad[seq_n_u_h] <- scaled_grad[seq_n_u_h] * ZAL_scaling 
      zInfo$scaled_grad <- scaled_grad
      if (trace>1L )  { ## only tracing
        maxs_grad <- c(max(abs(m_grad_obj[seq_n_u_h])),max(abs(m_grad_obj[-seq_n_u_h])))
        cat("iter=",innerj,", max(|grad|): v=",maxs_grad[1L],"beta=",maxs_grad[2L],";")
      }
      constant_APHLs_args <- list(processed=processed, which=processed$p_v_obj, sXaug=sXaug, phi_est=phi_est, lambda_est=lambda_est)
      # the following block needs m_grad_v the new m_grad_v hence its position
      if (is_HL1_1 && ! is.null(damped_WLS_blob)) {
        if (which_LevMar_step=="v") { ## 
          hlik_stuck <- (damped_WLS_blob$APHLs$hlik < oldAPHLs$hlik + Ftol_LM/10)
          if ( ! hlik_stuck) need_v_step <- (get_from_MME(sXaug=damped_WLS_blob$sXaug, which="Mg_invH_g", B=m_grad_v) > Ftol_LM/2) 
          if ( hlik_stuck || ! need_v_step) { ## LevMar apparently maximized h wrt v after several iterations
            ## if hlik has not recently moved or has moved but reached a point where the h gradient is negligible 
            if (trace>2L) print("switch from v to v_b")
            old_relV_beta <- relV_beta ## serves to assess convergence !!!!!!!!!!!!!!!!!!!
            which_LevMar_step <- "v_b" 
          } else {
            if (trace>2L) print("still v")
            ## v_h estimation not yet converged, continue with it
          }
        } else { ## performed one improvement of p_v by new v_b, 
          # indirect way of checking Mg_solve_g:
          # FIXME I made not_moving a sufficent condition fro break below !
          if (not_moving) { ## if we reach this point, Mg_solve_g (tested above) was too large, we must be out of solution space
            # need_v_step <- TRUE ## implicit meaning
          } else {
            p_v_stuck <- (damped_WLS_blob$APHLs$p_v < oldAPHLs$p_v + Ftol_LM/10) ## test whether LevMar apparently solved (v,beta) equations after several iterations
            if ( ! p_v_stuck) need_v_step <- (get_from_MME(sXaug=damped_WLS_blob$sXaug, which="Mg_invH_g", B=m_grad_v) > Ftol_LM/2) 
              ## we have identified two gradient cases where we must check v: Mg_solve_g>0 or (if estimates have just moved) Mg_invH_g>0 
          }
          if ( not_moving || p_v_stuck || need_v_step) { ## logically we may not need p_v_stuck, but this condition is faster to evaluate
            # p_v_stuck is analogous to (not_moving BUT large Mg_solve_g), checking that lik and estimates do not change 
            if (trace>2L) print("switch from v_b to v")
            which_LevMar_step <- "v"
          } else {
            if (trace>2L) print("still v_b")
          }
        }
      } 
      
      damped_WLS_blob <- .do_damped_WLS(sXaug=sXaug, zInfo=zInfo, 
                                        old_Vscaled_beta=Vscaled_beta,
                                        oldAPHLs=oldAPHLs,
                                        APHLs_args = constant_APHLs_args,
                                        damping=damping,
                                        ypos=ypos,off=off,
                                        GLMMbool=GLMMbool,etaFix=etaFix,
                                        constant_u_h_v_h_args=constant_u_h_v_h_args,
                                        updateW_ranefS_constant_arglist=updateW_ranefS_constant_arglist,
                                        wranefblob=wranefblob,seq_n_u_h=seq_n_u_h,ZAL_scaling=ZAL_scaling,
                                        processed=processed, Xscal=Xscal,mMatrix_method = mMatrix_method,
                                        phi_est=phi_est, H_global_scale=H_global_scale, n_u_h=n_u_h, ZAL=ZAL,
                                        which_i_affected_rows=which_i_affected_rows,
                                        which_LevMar_step = which_LevMar_step
      ) 
      ## LevM PQL
      if (! is_HL1_1) {
        if (damped_WLS_blob$lik < oldAPHLs$hlik) { ## if LevM step failed to find a damping that increases the lik
          ## This occurs inconspiscuously in the PQL_prefit providing a bad starting point for ML fit
          damped_WLS_blob <- NULL
          wzAug <- c(zInfo$y2_sscaled/ZAL_scaling, (zInfo$z1_sscaled)*weight_X)
          Vscaled_beta <- get_from_MME(sXaug,szAug=wzAug) # vscaled= v scaling so that v has 'scale' H_global_scale
          LevenbergM <- FALSE ## D O N O T set it to TRUE again !
        } 
      }
    } else { ## IRLS: always accept new v_h_beta
        damped_WLS_blob <- NULL
        Vscaled_beta <- get_from_MME(sXaug,szAug=wzAug) # vscaled= v scaling so that v has 'scale' H_global_scale
    }
    
    ######
    
    ##### Everything that is needed for 
    #  (1) assessment of convergence: c(v_h*sqrt(wranefblob$w.ranef),beta_eta)
    #  (2) all return elements are updated as function of the latest Vscaled_beta.
    #      In particular We need muetablob and (if ! LMM) sXaug, hence a lot of stuff.
    #  Hence, the following code is useful whether a break occurs or not. 
    if ( ! is.null(damped_WLS_blob) ) {
      Vscaled_beta <- damped_WLS_blob$Vscaled_beta
      eta <- damped_WLS_blob$eta
      wranefblob <- damped_WLS_blob$wranefblob
      v_h <- damped_WLS_blob$v_h
      u_h <- damped_WLS_blob$u_h
      muetablob <- damped_WLS_blob$muetablob
      w.resid <- damped_WLS_blob$w.resid ## !important! cf test-adjacency-corrMatrix.R
      # there's no new weight_X in damped_WLS_blob as the weights are purposefully kept constant there
      # same for new 'fitted' 
      sXaug <- damped_WLS_blob$sXaug
      if ( ! GLMMbool ) {
        Xscal <- damped_WLS_blob$Xscal
        ZAL_scaling <- damped_WLS_blob$ZAL_scaling
      }
    } else {
      # Vscaled_beta must have been provided by somethin else than damped_WLS_blob
      # drop, not as.vector(): names are then those of (final) eta and mu -> used by predict() when no new data
      fitted <- drop(Xscal %*% Vscaled_beta) ## length nobs+nr ! 
      fitted[ypos] <- eta <- fitted[ypos] + off
      if (is.null(etaFix$v_h)) { 
        v_h <- Vscaled_beta[seq_n_u_h] * ZAL_scaling
        if (GLMMbool) {
          u_h <- v_h ## keep input wranefblob since lambda_est not changed
        } else {
          u_h_v_h_from_v_h_args <- c(constant_u_h_v_h_args,list(v_h=v_h))
          u_h <- do.call(".u_h_v_h_from_v_h",u_h_v_h_from_v_h_args)
          if ( ! is.null(attr(u_h,"v_h"))) { ## second test = if u_h_info$upper.v_h or $lower.v_h non NULL
            v_h <- attr(u_h,"v_h")
          }
          ## update functions u_h,v_h
          wranefblob <- do.call(".updateW_ranefS",c(updateW_ranefS_constant_arglist,list(u_h=u_h,v_h=v_h)))
          if ( ! GLMMbool) { 
            ZAL_scaling <- 1/sqrt(wranefblob$w.ranef*H_global_scale) ## Q^{-1/2}/s
            if (TRUE) {
              ## alternative clause shows the meaning. This is distinctly faster. 
              scaledZAL <- .m_Matrix_times_Dvec(ZAL, ZAL_scaling)
              if (inherits(ZAL,"Matrix")) {
                # This block of code assumes Xscal (IO_ZX) has no @diag component
                Xscal@x[which_i_affected_rows] <- scaledZAL@x ## should create an error if some elements are stored in @diag
              } else {
                Xscal[n_u_h+seq(nobs),seq_n_u_h] <- scaledZAL
              }
            } else Xscal <- .make_Xscal(ZAL, ZAL_scaling = ZAL_scaling, AUGI0_ZX=processed$AUGI0_ZX,n_u_h=n_u_h)
          }
        }
      }
      muetablob <- .muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed) 
      if ( ! LMMbool ) {
        ## weight_X and Xscal vary within loop if ! LMM since at least the GLMweights in w.resid change
        w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est)
        weight_X <- .calc_weight_X(w.resid, H_global_scale) ## sqrt(s^2 W.resid)
        sXaug <- do.call(mMatrix_method,
                         list(Xaug=Xscal, weight_X=weight_X, w.ranef=wranefblob$w.ranef, H_global_scale=H_global_scale))
      } ## ergo sXaug is not updated for LMM (no need to)
    }
    beta_eta <- Vscaled_beta[n_u_h+seq_len(pforpv)]
    
    ##### assessment of convergence
    if (innerj<maxit.mean) {
      relV_beta <- c(v_h*sqrt(wranefblob$w.ranef),beta_eta)  ## convergence on v_h relative to sqrt(w.ranef)
      abs_d_relV_beta <- abs(relV_beta - old_relV_beta)
      not_moving <- ( ( ! is.null(old_relV_beta)) && mean(abs_d_relV_beta) < loc_Xtol_rel )
      if (is.na(not_moving)) {
        if (anyNA(relV_beta)) {
          if ( ! is.null(damped_WLS_blob)) {
            message(paste("innerj=",innerj,"damping=",damping,"lik=",damped_WLS_blob$lik))
            stop("Numerical problem despite Levenberg algorithm being used: complain.")
          } else stop("Numerical problem: try control.HLfit=list(LevenbergM=TRUE)")
        } else stop("Error in evaluating break condition")
      } 
      if (not_moving) break ## sufficient condition here
      if ( ! (is_HL1_1 && LevenbergM)) { ## possible reversal of condition from F to T in  LevM PQL !!!!
        old_relV_beta <- relV_beta
      } ## ELSE old_relV_beta controlled in block for which_LevMar_step !!
    } else break
  } ################ E N D LOOP ##############
  if (trace>1L && (LevenbergM))  { ## only tracing
    maxs_grad <- c(max(abs(m_grad_obj[seq_n_u_h])),max(abs(m_grad_obj[-seq_n_u_h])))
    cat("iter=",innerj,", max(|grad|): v=",maxs_grad[1L],"beta=",maxs_grad[2L],";")
  }
  
  
  names(beta_eta) <- colnames(X.pv)
  if (! is.null(damped_WLS_blob)) {
    fitted <- damped_WLS_blob$fitted
    weight_X <- damped_WLS_blob$weight_X
  } 
  RESU <- list(sXaug=sXaug, 
               ## used by calc_APHLs_from_ZX: (in particular can use Vscaled values contained in fitted)
               fitted=fitted, #zAug=zAug, 
               weight_X=weight_X, nobs=nobs, pforpv=pforpv, seq_n_u_h=seq_n_u_h, u_h=u_h, 
               muetablob=muetablob, 
               lambda_est=lambda_est,
               phi_est=phi_est,
               ## used by other code
               beta_eta=beta_eta, w.resid=w.resid, wranefblob=wranefblob, 
               v_h=v_h, eta=eta, innerj=innerj)
  return(RESU)
}

.intervalStep_ZX <- function(old_Vscaled_beta,sXaug,szAug,currentlik,for_intervals,currentDy) {
  #print((control.HLfit$intervalInfo$fitlik-currentlik)/(control.HLfit$intervalInfo$MLparm-old_betaV[parmcol]))
  ## voir code avant 18/10/2014 pour une implem rustique de VenzonM pour debugage  
  ## somewhat more robust algo (FR->FR: still improvable ?), updates according to a quadratic form of lik near max
  ## then target.dX = (current.dX)*sqrt(target.dY/current.dY) where dX,dY are relative to the ML x and y 
  ## A nice thing of this conception is that if the target lik cannot be approached, 
  ##   the inferred x converges to the ML x => this x won't be recognized as a CI bound (important for locoptim) 
  parmcol_ZX <- for_intervals$parmcol_ZX
  Vscaled_beta <- rep(NA,length(old_Vscaled_beta))
  if (currentDy <0) { 
    Vscaled_beta[parmcol_ZX] <- old_Vscaled_beta[parmcol_ZX]
  } else {
    currentDx <- (old_Vscaled_beta[parmcol_ZX]-for_intervals$MLparm)
    targetDy <- (for_intervals$fitlik-for_intervals$targetlik)
    Dx <- currentDx*sqrt(targetDy/currentDy)
    ## pb is if Dx=0 , Dx'=0... and Dx=0 can occur while p_v is still far from the target, because other params have not converged.
    ## FR->FR patch:
    if (currentDy<targetDy) { ## we are close to the ML: we extrapolate a bit more confidently
      min_abs_Dx <- for_intervals$asympto_abs_Dparm/1000
    } else min_abs_Dx <- 1e-6 ## far from ML: more cautious move our of Dx=0
    Dx <- sign(currentDx)*max(abs(Dx),min_abs_Dx)
    Vscaled_beta[parmcol_ZX] <- for_intervals$MLparm + Dx 
  }
  locsXaug <- sXaug[,-(parmcol_ZX),drop=FALSE]
  locszAug <- as.matrix(szAug-sXaug[,parmcol_ZX]*Vscaled_beta[parmcol_ZX])
  Vscaled_beta[-(parmcol_ZX)] <- get_from_MME(locsXaug,szAug=locszAug) 
  return(list(Vscaled_beta=Vscaled_beta)) # levQ ispresumably always dense
}

.intervalStep_glm <- function(old_beta,sXaug,szAug,currentlik,for_intervals,currentDy) {
  #print((control.HLfit$intervalInfo$fitlik-currentlik)/(control.HLfit$intervalInfo$MLparm-old_betaV[parmcol]))
  ## voir code avant 18/10/2014 pour une implem rustique de VenzonM pour debugage  
  ## somewhat more robust algo (FR->FR: still improvable ?), updates according to a quadratic form of lik near max
  ## then target.dX = (current.dX)*sqrt(target.dY/current.dY) where dX,dY are relative to the ML x and y 
  ## A nice thing of this conception is that if the target lik cannot be approached, 
  ##   the inferred x converges to the ML x => this x won't be recognized as a CI bound (important for locoptim) 
  parmcol_X <- for_intervals$parmcol_X
  beta <- rep(NA,length(old_beta))
  if (currentDy <0) { 
    beta[parmcol_X] <- old_beta[parmcol_X]
  } else {
    currentDx <- (old_beta[parmcol_X]-for_intervals$MLparm)
    targetDy <- (for_intervals$fitlik-for_intervals$targetlik)
    Dx <- currentDx*sqrt(targetDy/currentDy)
    ## pb is if Dx=0 , Dx'=0... and Dx=0 can occur while p_v is still far from the target, because other params have not converged.
    ## FR->FR patch:
    if (currentDy<targetDy) { ## we are close to the ML: we extrapolate a bit more confidently
      min_abs_Dx <- for_intervals$asympto_abs_Dparm/1000
    } else min_abs_Dx <- 1e-6 ## far from ML: more cautious move our of Dx=0
    Dx <- sign(currentDx)*max(abs(Dx),min_abs_Dx)
    beta[parmcol_X] <- for_intervals$MLparm + Dx 
  }
  locsXaug <- sXaug[,-(parmcol_X),drop=FALSE]
  locszAug <- as.matrix(szAug-sXaug[,parmcol_X]*beta[parmcol_X])
  beta[-(parmcol_X)] <- get_from_MME(locsXaug,szAug=locszAug) 
  return(list(beta=beta)) # levQ ispresumably always dense
}
