# evaluate tr(A %*% B)= sum(A*B) where A and B are large matrices but each of the form l %*% r for 'narrow' l
# this function avoids the formation of the large matrices, using a form of commutation of trace arguments.
.traceAB <- function(lA,rA,lB,rB) {
  ll <- .crossprod(lA, lB)
  rr <- .tcrossprod(rA, rB)
  return(sum(ll*rr))
}
# and same concept for trace( D %*% B)
.traceDB <- function(dD,lB,rB) { sum( sweep(lB * t(rB),1L,dD,`*`) )}
# lA <- matrix(runif(6),ncol=2)
# rA <- matrix(runif(6),ncol=3)
# lB <- matrix(runif(6),ncol=2)
# rB <- matrix(runif(6),ncol=3)
# .traceAB(lA,rA,lB,rB)
# dD <- runif(3)
# sum(diag(diag(dD) %*% lB %*% rB))
# .traceDB(dD,lB,rB)



.calc_invV.dV_info <- function(strucList,lambda.object,corrPars,asDmLR_invV,ZAL) {
  nrand <- length(strucList)
  cum_n_u_h <- attr(lambda.object$lambda,"cum_n_u_h")
  RES <- list(lambda=numeric(nrand),cum_n_u_h=cum_n_u_h)
  ZAL_to_ZALd_vec <- rep(1,cum_n_u_h[nrand+1L])
  for (randit in seq_along(strucList)) {
    lmatrix <- strucList[[randit]]
    corr.model <- attr(lmatrix,"corr.model")
    if (identical(corr.model,"adjacency") || identical(corr.model,"ar1")) { ## allows only one AR model
      glmit <- lambda.object$rand_to_glm_map[randit]
      RES$lambda[[randit]] <- with(lambda.object,linkinvS[[glmit]](coefficients_lambdaS[[randit]][1]))
      if ("adjd" %in% lambda.object$namesTerms[[randit]]) {
        RES$rho <- - with(lambda.object,coefficients_lambdaS[[randit]][2]/coefficients_lambdaS[[randit]][1])
      } else { RES$rho <- corrPars$rho }
      adjd <- attr(lmatrix,"symSVD")$adjd
      denom <- 1-RES$rho*adjd
      RES$adjd_denom2 <- adjd/(denom^2) 
      u.range <- (cum_n_u_h[randit]+1L):(cum_n_u_h[randit+1L])
      ZAL_to_ZALd_vec[u.range] <- ZAL_to_ZALd_vec[u.range]/sqrt(denom)
    } else if (identical(corr.model,"random-coef")) {
      stop("code missing for random-coefficient models.")
    } else { ## standard lamScal model
      coeff <- lambda.object$coefficients_lambdaS[[randit]] ## corrHLfit or fitme with control$refit=TRUE
      if (is.null(coeff)) {
        RES$lambda[[randit]] <- lambda.object$lambda[[randit]] ## basic fitme (may numerically differ)
      }  else RES$lambda[[randit]] <- exp(coeff)
    }
  }
  ZALd <- ZAL %id*id% Diagonal(x=ZAL_to_ZALd_vec)
  ## next 2 lines uses invV= invD- n_x_r %*% r_x_n
  # and the lhs of dvdlam = ZALd %*% t(ZALd)
  lhs_invV.dVdlam <- asDmLR_invV$n_x_r %*% (asDmLR_invV$r_x_n %*% ZALd)  
  RES$lhs_invV.dVdlam <- sweep( ZALd,1L,asDmLR_invV$invD,`*`) - lhs_invV.dVdlam # (invD- n_x_r %*% r_x_n) %*% ZALd 
  RES$rhs_invV.dVdlam <- t(ZALd)
  # lambda^2 *sum(invV.dVdlam* t(invV.dVdlam)) :
  return(RES)
}

.calc_loglamInfo <- function(invV.dV_info,which) {
  lambda <- invV.dV_info$lambda[which]
  nrand <- length(lambda)
  cum_n_u_h <- invV.dV_info$cum_n_u_h
  loglamInfo <- matrix(ncol=length(which),nrow=length(which))
  for (randit in which) {
    iloc <- rep(0,cum_n_u_h[nrand+1L])
    u.range <- (cum_n_u_h[randit]+1L):(cum_n_u_h[randit+1L])
    iloc[u.range] <- 1L
    i_rhs_invV.dVdlam <- sweep( invV.dV_info$rhs_invV.dVdlam,1L,iloc,`*`)
    loglamInfo[randit,randit] <- .traceAB(invV.dV_info$lhs_invV.dVdlam,i_rhs_invV.dVdlam,
                                          t(i_rhs_invV.dVdlam),t(invV.dV_info$lhs_invV.dVdlam))
    for (randjt in intersect(which,seq_len(randit-1L))) {
      jloc <- rep(0,cum_n_u_h[nrand+1L])
      u.range <- (cum_n_u_h[randjt]+1L):(cum_n_u_h[randjt+1L])
      jloc[u.range] <- 1
      j_rhs_invV.dVdlam <- sweep( invV.dV_info$rhs_invV.dVdlam,1L,jloc,`*`)
      loglamInfo[randjt,randit] <- .traceAB(invV.dV_info$lhs_invV.dVdlam,i_rhs_invV.dVdlam, 
                                            t(j_rhs_invV.dVdlam),t(invV.dV_info$lhs_invV.dVdlam))
      loglamInfo[randit,randjt] <- loglamInfo[randjt,randit]
    }
  }
  loglamInfo <- loglamInfo * (lambda[which] %*% t(lambda[which]))
  return(loglamInfo)
} 


## version with strucList (and other modifs), called for recent fits
.calc_logdisp_cov_new <- function(object, dvdloglamMat=NULL, dvdlogphiMat=NULL, asDmLR_invV=NULL) { 
  lambda.object <- object$lambda.object
  strucList <- object$strucList
  dwdlogphi <- dwdloglam <- NULL ## always a cbind at the end of calc_logdisp_cov
  dispcolinfo <- list()
  problems <- list() ## Its elements' names are tested in calcPredVar, and the strings are 'development info'
  # detect if fitme() without refit 
  #  if (any(lambda.object$type=="outer")) {
  #    ## then we need to refit to get a more complete lambda.object
  #  }
  #coefficients_lambdaS <- lambda.object$coefficients_lambdaS
  if (any(lambda.object$type!="fixed")) {
    corr.models <- lapply(strucList,attr,which="corr.model") ## not unlist bc it may contain NULLs
    checkadj <- unlist(lapply(corr.models,identical,y="adjacency"))
    if(any(checkadj)) problems$warnadj <- warning("lambda dispVar component not implemented for adjacency model.")
    checkar1 <- unlist(lapply(corr.models,identical,y="ar1"))
    if(any(checkar1)) problems$warnauto <- warning("lambda dispVar component not implemented for autoregressive model.")
    checkrancoef <- unlist(lapply(corr.models,identical,y="random-coef"))
    if(any(checkrancoef)) problems$warnauto <- warning("lambda dispVar component not implemented for random-coefficients model.")
    others <- ! (checkadj | checkar1 | checkrancoef)
    if (any(others)) {
      if (is.null(dvdloglamMat)) {
        ## note that .get_logdispObject is computed on request by .get_logdispObject()
        # but only if all lambda coeffs are Intercept coeffs.
        problems$stopmiss <- warning("is.null(dvdloglamMat) in a case where it should be available.") 
      } else {
        checklambda <- (others & lambda.object$type!="fixed")
        if (any(checklambda)) dispcolinfo$loglambda <- "loglambda"
      }
      
    } 
    dvdloglam <- matrix(0,ncol=length(checklambda),nrow=nrow(dvdloglamMat))
    cum_n_u_h <- attr(lambda.object$lambda,"cum_n_u_h")
    nrand <- length(strucList)
    for (randit in which(checklambda)) {
      u.range <- (cum_n_u_h[randit]+1L):(cum_n_u_h[randit+1L])
      dvdloglam[u.range,randit] <- rowSums(dvdloglamMat[u.range,]) ## assuming each lambda_i = lambda in each block
    }
    dwdloglam <- .calc_invL_coeffs(object,dvdloglam) ## one matrix for all ranefs 
  } else lambda <- NULL # no ranefs
  ###
  phimodel <- object$models[["phi"]]
  if (phimodel=="phiScal") { ## semble impliquer pas outer phi.Fix... => no need to test object$phi.object$phi_outer,"type")
    phi_est <- object$phi ## no need to get object$phi.object$phi_outer
    if (length(phi_est)!=1L) problems$stopphi <- warning("phimodel=\"phiScal\" but length(phi_est)!=1L.")
    if ( ! is.null(dvdlogphiMat)) {
      dvdlogphi <- rowSums(dvdlogphiMat) ## using each phi_i = phi
      dwdlogphi <- .calc_invL_coeffs(object,dvdlogphi)
      dispcolinfo$logphi <- "logphi"
    } else if (object$models[["eta"]]=="etaHGLM") stop("phimodel=='phiScal' but is.null(dvdlogphiMat)")
  } else {  ## else phimodel="", e.g. binomial
    # if binomial or poisson, phimodel=""; warning for other phimodels
    if (phimodel!="") {
      problems$structphi <- "phi dispVar component not yet available for phi model != ~1."
      if ( ! identical(spaMM.getOption("phi_dispVar_comp_warned"),TRUE)) {
        warning(problems$structphi)
        spaMM.options(phi_dispVar_comp_warned=TRUE)
      }
    }
  }
  dwdlogdisp <- cbind(dwdloglam,dwdlogphi) ## typically nobs * 2
  col_info <- list(nrand=nrand, ranef_ids=c(), phi_cols=0L)
  if (!is.null(dwdloglam)) {
    col_info$ranef_ids <- which(checklambda) ## indices of ranefs, not cols of ranefs
    col_info$cum_n_u_h <- cum_n_u_h
  }
  if (!is.null(dwdlogphi)) col_info$phi_cols=length(col_info$ranef_ids)+seq_len(NCOL(dwdlogphi)) ## cols indices for phi 
  attr(dwdlogdisp,"col_info") <- col_info
  ## compute info matrix:
  if ((length(dispcolinfo))==0L) {
    return(list(problems=problems))
  } else {
    # cf my documentation, based on McCullochSN08 6.62 and 6.74
    # lambda and phi factors enter in dV/dlog(.), computed instead of dV/d(.) to match dwdlog(.) vectors.
    #
    # use repres of two matrices large A and B, each as (thin) lhs %*% (flat) rhs   
    ZAL <- get_ZALMatrix(object, as_matrix=.eval_as_mat_arg(object))
    if ("loglambda" %in% names(dispcolinfo) || "rho" %in% names(dispcolinfo)) {
      invV.dV_info <- .calc_invV.dV_info(strucList=strucList[checklambda],lambda.object=lambda.object,
                                         corrPars=object$corrPars, asDmLR_invV=asDmLR_invV, ZAL=ZAL)
      dispcolinfo$loglambda <- rep("loglambda",length(invV.dV_info$lambda[checklambda]))
    }
    cum_n_disp_pars <- cumsum(c(0,lapply(dispcolinfo,length))) # #ncols for phi, lambda[checklambda], etc.
    dispcols <- lapply(seq_along(dispcolinfo), function(varit) {
      cum_n_disp_pars[varit]+ seq_along(dispcolinfo[[varit]])
    }) ## col ranges for phi, lambda[checklambda], etc
    names(dispcols) <- dispnames <- names(dispcolinfo) ## list names
    nrc <- cum_n_disp_pars[length(cum_n_disp_pars)]
    #
    logdispInfo <- matrix(NA,nrow=nrc,ncol=nrc)
    colnames(logdispInfo) <- rownames(logdispInfo) <- unlist(dispcolinfo)
    if ("loglambda" %in% dispnames) { 
      loglamInfo <- .calc_loglamInfo(invV.dV_info,which=which(checklambda))
      logdispInfo[dispcols$loglambda,dispcols$loglambda] <- loglamInfo 
    }
    if ("rho" %in% dispnames) {
      # no use of sqrt because adjd can be negative
      #invV.dVdrho <- (invV %id*id% ZAL) %*% ( Diagonal(x=lambda*adjd/(denom^2)) %id*id% t(ZAL))
      lhs_invV.dVdrho <- asDmLR_invV$n_x_r %*% (asDmLR_invV$r_x_n %*% ZAL)  
      lhs_invV.dVdrho <- sweep( ZAL,1L,asDmLR_invV$invD,`*`) - lhs_invV.dVdrho
      rhs_invV.dVdrho <- ( Diagonal(x=lambda*invV.dV_info$adjd_denom2) %id*id% t(ZAL)) ## FIXME curently meaningful for only one lambda element
      #logdispInfo["rho","rho"] <- sum(invV.dVdrho*t(invV.dVdrho))
      logdispInfo[dispcols$rho,dispcols$rho] <- .traceAB(lhs_invV.dVdrho,rhs_invV.dVdrho,t(rhs_invV.dVdrho),t(lhs_invV.dVdrho))
      if ("loglambda" %in% dispnames) {
        for (randit in which(checklambda)) {
          iloc <- rep(0,cum_n_u_h[nrand+1L])
          u.range <- (cum_n_u_h[randit]+1L):(cum_n_u_h[randit+1L])
          iloc[u.range] <- 1L
          i_rhs_invV.dVdlam <- sweep( invV.dV_info$rhs_invV.dVdlam,1L,iloc,`*`)
          logdispInfoBlock[randit] <- .traceDB(asDmLR_invV$invD, t(rhs_invV.dVdrho),t(lhs_invV.dVdrho)) -
            .traceAB(invV.dV_info$lhs_invV.dVdlam, i_rhs_invV.dVdlam,t(rhs_invV.dVdrho),t(lhs_invV.dVdrho))
        }
        logdispInfoBlock <- logdispInfoBlock *invV.dV_info$lambda  #lambda * sum(invV.dVdlam*t(invV.dVdrho))
        logdispInfo[dispcols$loglambda,dispcols$rho] <- 
          logdispInfo[dispcols$rho,dispcols$loglambda] <- logdispInfoBlock
      }
    } 
    ## if (! is.null(dwdlogphi)) { ## currently length(phi)==1L && ! is.null(dvdlogphiMat)
    if ("logphi" %in% dispnames) { ## more transparent, but error if mismatch of conditions
      ## next lines assume that  the design matrix for the residual error is I
      logdispInfo[dispcols$logphi,dispcols$logphi] <- phi_est^2 * (
        sum(asDmLR_invV$invD^2) -2 * .traceDB(asDmLR_invV$invD,asDmLR_invV$n_x_r, asDmLR_invV$r_x_n) + 
          .traceAB(asDmLR_invV$n_x_r, asDmLR_invV$r_x_n, asDmLR_invV$n_x_r, asDmLR_invV$r_x_n)
      ) # phi_est^2 * sum(invV^2)
      if ("loglambda" %in% dispnames) {
        nlambda <- length(invV.dV_info$lambda) ## nrand ?
        logdispInfoBlock <- numeric(nlambda)
        cum_n_u_h <- invV.dV_info$cum_n_u_h
        for (randit in which(checklambda)) {
          iloc <- rep(0,cum_n_u_h[nrand+1L])
          u.range <- (cum_n_u_h[randit]+1L):(cum_n_u_h[randit+1L])
          iloc[u.range] <- 1L
          i_rhs_invV.dVdlam <- sweep( invV.dV_info$rhs_invV.dVdlam,1L,iloc,`*`)
          logdispInfoBlock[randit] <- .traceDB(asDmLR_invV$invD, invV.dV_info$lhs_invV.dVdlam, i_rhs_invV.dVdlam) -
            .traceAB(invV.dV_info$lhs_invV.dVdlam, i_rhs_invV.dVdlam, asDmLR_invV$n_x_r, asDmLR_invV$r_x_n)
        }
        logdispInfoBlock <- logdispInfoBlock *invV.dV_info$lambda * phi_est # lambda * phi_est * sum(invV.dVdlam * invV)
        logdispInfo[dispcols$loglambda,dispcols$logphi] <- 
          logdispInfo[dispcols$logphi,dispcols$loglambda] <- logdispInfoBlock
      }
      if ("rho" %in% dispnames) {
        logdispInfo[dispcols$rho,dispcols$logphi] <- 
          logdispInfo[dispcols$logphi,dispcols$rho] <- phi_est * .traceAB(lhs_invV.dVdrho,rhs_invV.dVdrho, 
                                                                          asDmLR_invV$n_x_r, asDmLR_invV$r_x_n)  
        # phi_est * sum(invV.dVdrho * invV)  
      }
    } 
    logdispInfo <- logdispInfo/2
    logdisp_cov <- try(solve(logdispInfo),silent=TRUE)
    if (inherits(logdisp_cov,"try-error")) logdisp_cov <- MASS::ginv(logdispInfo) ## quick patch for uninteresting case
    return(list(dwdlogdisp=dwdlogdisp,logdisp_cov=logdisp_cov,problems=problems)) ## more compact than storing ww %*% logdisp_cov %*% t(ww) which is nobs*nobs 
  }
}

calc_logdisp_cov <- function(object, dvdloglamMat=NULL, dvdlogphiMat=NULL, asDmLR_invV=NULL) {
  if (object$spaMM.version>"1.11.60") {
    .calc_logdisp_cov_new(object=object, dvdloglamMat=dvdloglamMat, dvdlogphiMat=dvdlogphiMat, asDmLR_invV=asDmLR_invV)
  } else {
    .calc_logdisp_cov_old(object=object, dvdloglamMat=dvdloglamMat, dvdlogphiMat=dvdlogphiMat, asDmLR_invV=asDmLR_invV, stop.on.error=TRUE)
  }
} 
