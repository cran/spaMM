.calc_logdisp_cov_old <- function(object, dvdloglamMat=NULL, dvdlogphiMat=NULL, asDmLR_invV=NULL, stop.on.error) { 
  
  # evaluate tr(A %*% B)= sum(A*B) where A and B are large matrices but each of the form l %*% r for 'narrow' l
  # this function avoids the formation of the large matrices, using a form of commutation of trace arguments.
  traceAB <- function(lA,rA,lB,rB) {
    ll <- .crossprod(lA, lB)
    rr <- .tcrossprod(rA, rB)
    return(sum(ll*rr))
  }
  # and same concept for trace( D %*% B)
  traceDB <- function(dD,lB,rB) { sum( sweep(lB * t(rB),1L,dD,`*`) )}
  # lA <- matrix(runif(6),ncol=2)
  # rA <- matrix(runif(6),ncol=3)
  # lB <- matrix(runif(6),ncol=2)
  # rB <- matrix(runif(6),ncol=3)
  # traceAB(lA,rA,lB,rB)
  # dD <- runif(3)
  # sum(diag(diag(dD) %*% lB %*% rB))
  # traceDB(dD,lB,rB)
  
  lambda.object <- object$lambda.object
  LMatrix <- attr(object$predictor,"LMatrix") ## back ompat given we are in .calc_logdisp_cov_old()
  adjacency <- identical(attr(LMatrix,"corr.model"),"adjacency")
  ar1 <- identical(attr(LMatrix,"corr.model"),"ar1")
  ## debut de code de detection parmi +s ranefs
  #whichadj <- attr(attr(ZAlist,"ranefs"),"type")=="adjacency"  
  #notscalar <- unlist(lapply(coefficients_lambdaS,length)==2L)
  dwdlogphi <- dwdloglam <- NULL ## always a cbind at the end of calc_logdisp_cov
  dispnames <- c()
  problems <- list() ## Its elements' names are tested in calcPredVar, and the strings are 'development info'
  if (length(lambda.object$coefficients_lambda)>0L) { ## not lambda.Fix
    ## in the first cases, a dvdloglamMat may exist (eg salamander, male and female effects)
    if (adjacency || ar1) {
      ## cf summary.HLfit for extracting info from object
      locit <- 1L
      it <- 1L
      if ("adjd" %in% object$lambda.object$namesTerms[[it]]) {
        rho <- - with(lambda.object,coefficients_lambda[locit+1L]/coefficients_lambda[locit])
        lambda <- with(lambda.object,linkinvS[[rand_to_glm_map[it]]](coefficients_lambda[locit]))
      } else {
        lambda <- with(lambda.object,linkinvS[[rand_to_glm_map[it]]](coefficients_lambda[locit]))        
        rho <- object$corrPars$rho
      }
      adjd <- attr(LMatrix,"symSVD")$adjd
      denom <- 1-rho*adjd
      problems$warnadj <- warning("lambda dispVar component not implemented for adjacency model.") 
      # il me manque dwdrho (et meme dwdloglam pour ce modele ?) donc on inactive les lignes suivantes:
      #       if (is.null(lambda.object$lambda.fix)) dispnames <- c(dispnames,"loglambda")
      #       corrFixNames <- names(unlist(object$corrPars[which(attr(corrPars,"type")=="fix")]))
      #       if (! ("rho" %in% corrFixNames) ) dispnames <- c(dispnames,"rho")
    } else if (length(lambda.object$coefficients_lambda)>1L) { 
      problems$warnmult <- warning("lambda dispVar component not implemented for model with several lambda parameters.")
    } else if (is.null(dvdloglamMat)) {
      problems$stopmiss <- warning("is.null(dvdloglamMat) in a case where it should be available.") 
    }else {
      dvdloglam <- rowSums(dvdloglamMat) ## assuming each lambda_i = lambda
      dwdloglam <- .calc_invL_coeffs(object,dvdloglam)
      if ( all(lambda.object$type=="inner")) dispnames <- c(dispnames,"loglambda")
    }
  } else lambda <- NULL # no ranefs
  
  phimodel <- object$models[["phi"]]
  if (phimodel=="phiScal") { ## semble impliquer pas outer phi.Fix... => no need to test object$phi.object$phi_outer,"type")
    phi_est <- object$phi ## no need to get object$phi.object$phi_outer
    if (length(phi_est)!=1L) problems$stopphi <- warning("phimodel=\"phiScal\" but length(phi_est)!=1L.")
    if ( ! is.null(dvdlogphiMat)) {
      dvdlogphi <- rowSums(dvdlogphiMat) ## using each phi_i = phi
      dwdlogphi <- .calc_invL_coeffs(object,dvdlogphi)
      dispnames <- c(dispnames,"logphi")
    } else if (object$models[["eta"]]=="etaHGLM") stop("phimodel=='phiScal' but is.null(dvdlogphiMat)")
  } else {  ## else phimodel="", e.g. binomial
    # if binomial or poisson, phimodel=""; warning for other phimodels
    if (phimodel!="") problems$structphi <- warning("phi dispVar component not yet available for phi model != ~1.")
  }
  ## compute info matrix:
  if ((nrc <- length(dispnames))==0L) {
    return(list(problems=problems))
  } else {
    # cf my documentation, based on McCullochSN08 6.62 and 6.74
    # lambda and phi factors enter in dV/dlog(.), computed instead of dV/d(.) to match dwdlog(.) vectors.
    #
    # use repres of two matrices large A and B, each as (thin) lhs %*% (flat) rhs   
    ZAL <- get_ZALMatrix(object,as_matrix=.eval_as_mat_arg(object))
    logdispInfo <- matrix(NA,nrow=nrc,ncol=nrc)
    colnames(logdispInfo) <- rownames(logdispInfo) <- dispnames
    if ("loglambda" %in% dispnames) {
      if (adjacency || ar1) { 
        #  lambda already evaluated 
        ZALd <- ZAL %id*id% Diagonal(x=sqrt(1/denom))
        lhs_invV.dVdlam <- asDmLR_invV$n_x_r %*% (asDmLR_invV$r_x_n %*% ZALd)  
        lhs_invV.dVdlam <- sweep( ZALd,1L,asDmLR_invV$invD,`*`) - lhs_invV.dVdlam
        rhs_invV.dVdlam <- t(ZALd)
      } else { ## standard lamScal model
        lambda <- exp(lambda.object$coefficients_lambda)
        lhs_invV.dVdlam <- asDmLR_invV$n_x_r %*% (asDmLR_invV$r_x_n %*% ZAL)  
        lhs_invV.dVdlam <- sweep( ZAL,1L,asDmLR_invV$invD,`*`) - lhs_invV.dVdlam
        rhs_invV.dVdlam <- t(ZAL)
      }
      # lambda^2 *sum(invV.dVdlam* t(invV.dVdlam)) :
      logdispInfo["loglambda","loglambda"] <- lambda^2 *traceAB(lhs_invV.dVdlam,rhs_invV.dVdlam,
                                                                t(rhs_invV.dVdlam),t(lhs_invV.dVdlam)) 
    }
    if ("rho" %in% dispnames) {
      # no use of sqrt because adjd can be negative
      #invV.dVdrho <- (invV %id*id% ZAL) %*% ( Diagonal(x=lambda*adjd/(denom^2)) %id*id% t(ZAL))
      lhs_invV.dVdrho <- asDmLR_invV$n_x_r %*% (asDmLR_invV$r_x_n %*% ZAL)  
      lhs_invV.dVdrho <- sweep( ZAL,1L,asDmLR_invV$invD,`*`) - lhs_invV.dVdrho
      rhs_invV.dVdrho <- ( Diagonal(x=lambda*adjd/(denom^2)) %id*id% t(ZAL))
      #logdispInfo["rho","rho"] <- sum(invV.dVdrho*t(invV.dVdrho))
      logdispInfo["rho","rho"] <- traceAB(lhs_invV.dVdrho,rhs_invV.dVdrho,t(rhs_invV.dVdrho),t(lhs_invV.dVdrho))
      if ("loglambda" %in% dispnames) {
        logdispInfo["loglambda","rho"] <- 
          logdispInfo["rho","loglambda"] <- lambda* ( 
            traceDB(asDmLR_invV,t(rhs_invV.dVdrho),t(lhs_invV.dVdrho)) -
              traceAB(lhs_invV.dVdlam,rhs_invV.dVdlam,t(rhs_invV.dVdrho),t(lhs_invV.dVdrho)) 
          ) #lambda * sum(invV.dVdlam*t(invV.dVdrho))
      }
    } 
    ## if (! is.null(dwdlogphi)) { ## currently length(phi)==1L && ! is.null(dvdlogphiMat)
    if ("logphi" %in% dispnames) { ## more transparent, but error if mismatch of conditions
      ## next lines assume that  the design matrix for the residual error is I
      logdispInfo["logphi","logphi"] <- phi_est^2 * (
        sum(asDmLR_invV$invD^2) -2*traceDB(asDmLR_invV$invD,asDmLR_invV$n_x_r, asDmLR_invV$r_x_n) + 
          traceAB(asDmLR_invV$n_x_r, asDmLR_invV$r_x_n, asDmLR_invV$n_x_r, asDmLR_invV$r_x_n)
      ) # phi_est^2 * sum(invV^2)
      if ("loglambda" %in% dispnames) {
        logdispInfo["loglambda","logphi"] <- 
          logdispInfo["logphi","loglambda"] <- lambda * phi_est * (
            traceDB(asDmLR_invV$invD, lhs_invV.dVdlam, rhs_invV.dVdlam) -
              traceAB(lhs_invV.dVdlam, rhs_invV.dVdlam, asDmLR_invV$n_x_r, asDmLR_invV$r_x_n)
          ) # lambda * phi_est * sum(invV.dVdlam * invV)
      }
      if ("rho" %in% dispnames) {
        logdispInfo["rho","logphi"] <- 
          logdispInfo["logphi","rho"] <- phi_est * traceAB(lhs_invV.dVdrho,rhs_invV.dVdrho, asDmLR_invV$n_x_r, asDmLR_invV$r_x_n)  
        # phi_est * sum(invV.dVdrho * invV)  
      }
    } 
    logdispInfo <- logdispInfo/2
    logdisp_cov <- try(solve(logdispInfo),silent=TRUE)
    if (inherits(logdisp_cov,"try-error")) logdisp_cov <- MASS::ginv(logdispInfo) ## quick patch for uninteresting case
    dwdlogdisp <- cbind(dwdloglam,dwdlogphi) ## typically nobs * 2
    return(list(dwdlogdisp=dwdlogdisp,logdisp_cov=logdisp_cov,problems=problems)) ## more compact than storing ww %*% logdisp_cov %*% t(ww) which is nobs*nobs 
  }
}