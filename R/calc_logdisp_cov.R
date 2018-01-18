# evaluate tr(A %*% B)= sum(A*B) where A and B are large matrices but each of the form l %*% r for 'narrow' l
# this function avoids the formation of the large 'n x n' matrices, using a form of commutation of trace arguments.
# However the computation of 'r x r' matrix by the crossproducts may still be quite long
# e.g. 792*62388 matrices => (792^2)*62388 = 39 133 746 432 flop'product' and quite a few 'sum'
# FIXME that's something that 'should' be easy to paralellise
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

.calc_lhs_invV.dVdlam <- function(object, ZALd, asDmLR_invV) {
  if ("AUGI0_ZX_sparsePrecision" %in% object$MME_method) {
    ## ~ .Sigsolve_sparsePrecision !
    invG_ZtW_ <- Matrix::solve(object$envir$G_CHMfactor, object$envir$ZtW %*% ZALd)
    ZinvG_ZtW_ <- asDmLR_invV$ZAfix %*% invG_ZtW_
    lhs_invV.dVdlam <- object$w.resid*(ZALd - ZinvG_ZtW_) ## implicit .Dvec_times_mMatrix
  } else {
    ## next 2 lines uses invV= w.resid- n_x_r %*% r_x_n
    W_ZinvG_ZtW_ <- asDmLR_invV$n_x_r %*% (asDmLR_invV$r_x_n %*% ZALd)  
    lhs_invV.dVdlam <- .Dvec_times_m_Matrix(object$w.resid, ZALd) - W_ZinvG_ZtW_ # (w.resid- n_x_r %*% r_x_n) %*% ZALd 
  }
  return(lhs_invV.dVdlam)
}  

.calc_invV.dV_info <- function(object,checklambda,ZAL,asDmLR_invV) {
  strucList <- object$strucList
  lambda.object <- object$lambda.object
  cum_n_u_h <- attr(lambda.object$lambda,"cum_n_u_h")
  lambda_list <- vector("list",length(strucList))
  RES <- list(cum_n_u_h=cum_n_u_h)
  ZAL_to_ZALd_vec <- rep(0,cum_n_u_h[length(cum_n_u_h)]) ## rep(666,.) shoudl be equivalent (as this is *0) (but NA*0 =NA)
  for (randit in which(checklambda)) {
    lmatrix <- strucList[[randit]]
    corr.model <- attr(lmatrix,"corr.model")
    u.range <- (cum_n_u_h[randit]+1L):(cum_n_u_h[randit+1L])
    if (identical(corr.model,"adjacency") &&
        ( ! is.null(glmit <- lambda.object$rand_to_glm_map[randit]))
       ) { 
      lambda_list[[randit]] <- with(lambda.object,linkinvS[[glmit]](coefficients_lambdaS[[randit]][1]))
      if ("adjd" %in% lambda.object$print_namesTerms[[randit]]) {
        RES$rho <- - with(lambda.object,coefficients_lambdaS[[randit]][2]/coefficients_lambdaS[[randit]][1])
      } else { RES$rho <- object$corrPars$rho }
      adjd <- attr(lmatrix,"symSVD")$adjd
      denom <- 1-RES$rho*adjd
      RES$adjd_denom2 <- adjd/(denom^2) 
      ZAL_to_ZALd_vec[u.range] <- ZAL_to_ZALd_vec[u.range]/sqrt(denom)
    } else { ## standard lamScal model or ranCoefs model # F I X M E check that this is OK for ranCoefs
      ZAL_to_ZALd_vec[u.range] <- 1
      coeff <- lambda.object$coefficients_lambdaS[[randit]] ## corrHLfit or fitme with control$refit=TRUE
      if (is.null(coeff)) {
        lambda_list[[randit]] <- lambda.object$lambda[[randit]] ## basic fitme (may numerically differ)
      } else lambda_list[[randit]] <- exp(coeff)
    }
  }
  RES$lambda_list <- lambda_list
  ZALd <- .m_Matrix_times_Dvec(ZAL, ZAL_to_ZALd_vec)
  # The rhs of dvdlam = ZALd %*% t(ZALd) is split between lhs_invV.dVdlam and rhs_invV.dVdlam
  RES$lhs_invV.dVdlam <- .calc_lhs_invV.dVdlam(object, ZALd, asDmLR_invV)
  RES$rhs_invV.dVdlam <- t(ZALd)
  return(RES)
}

.fill_rhs_invV.dVdlam <- function(template, urange, rhs_invV.dVdlam) { ## to localise template and urange
  template[urange] <- 1L
  return(.Dvec_times_m_Matrix(template,rhs_invV.dVdlam)) # sweep( invV.dV_info$rhs_invV.dVdlam,1L,iloc,`*`))
}

.calc_loglamInfo <- function(invV.dV_info,which) { ## called by .calc_logdisp_cov_new(), using the result of .calc_invV.dV_info()
  lambda_list <- invV.dV_info$lambda_list ## for ranCoefs, must be a list with Xi_cols-matching elements 
  Xi_cols <- sapply(lambda_list,length)
  n_sublambda <- sum(Xi_cols[which]) 
  loglamInfo <- matrix(ncol=n_sublambda,nrow=n_sublambda)
  cum_n_u_h <- invV.dV_info$cum_n_u_h
  zerotemplate <- rep(0,cum_n_u_h[length(cum_n_u_h)]) 
  cum_Xi_cols <- cumsum(c(0,Xi_cols))
  for (randit in which) { # say '2' for second ranef
    u.range <- (cum_n_u_h[randit]+1L):(cum_n_u_h[randit+1L])
    Xi_ncol <- Xi_cols[randit] # say '1 2' for ranCoefs
    uirange <- matrix(u.range,ncol=Xi_ncol)
    for (ilam in seq_len(Xi_ncol)) { 
      i_rhs_invV.dVdlam <- .fill_rhs_invV.dVdlam(template=zerotemplate, urange=uirange[,ilam], invV.dV_info$rhs_invV.dVdlam)
      colit <- cum_Xi_cols[randit]+ilam
      loglamInfo[colit,colit] <- .traceAB(invV.dV_info$lhs_invV.dVdlam,i_rhs_invV.dVdlam,
                                          t(i_rhs_invV.dVdlam),t(invV.dV_info$lhs_invV.dVdlam))
      for (jlam in seq_len(ilam-1L)) { ## WITHIN a randit
        j_rhs_invV.dVdlam <- .fill_rhs_invV.dVdlam(template=zerotemplate, urange=uirange[,jlam], invV.dV_info$rhs_invV.dVdlam)
        coljt <- cum_Xi_cols[randit]+jlam
        loglamInfo[coljt,colit] <- .traceAB(invV.dV_info$lhs_invV.dVdlam,i_rhs_invV.dVdlam, 
                                            t(j_rhs_invV.dVdlam),t(invV.dV_info$lhs_invV.dVdlam))
        loglamInfo[colit,coljt] <- loglamInfo[coljt,colit]
      }
      for (randjt in intersect(which,seq_len(randit-1L))) {
        u.range <- (cum_n_u_h[randjt]+1L):(cum_n_u_h[randjt+1L])
        Xj_ncol <- Xi_cols[randjt] # say '1 2' for ranCoefs
        ujrange <- matrix(u.range,ncol=Xj_ncol)
        for (jlam in seq_len(Xj_ncol)) {
          j_rhs_invV.dVdlam <- .fill_rhs_invV.dVdlam(template=zerotemplate, urange=ujrange[,jlam], invV.dV_info$rhs_invV.dVdlam)
          coljt <- cum_Xi_cols[randjt]+jlam
          loglamInfo[coljt,colit] <- .traceAB(invV.dV_info$lhs_invV.dVdlam,i_rhs_invV.dVdlam, 
                                              t(j_rhs_invV.dVdlam),t(invV.dV_info$lhs_invV.dVdlam))
          loglamInfo[colit,coljt] <- loglamInfo[coljt,colit]
        }
      }
    }
  }
  sub_lambda_vec <- unlist(lambda_list[which]) 
  loglamInfo <- loglamInfo * (sub_lambda_vec %*% t(sub_lambda_vec))
  return(loglamInfo)
} 

.calc_logdisp_cov <- function(object, dvdloglamMat=NULL, dvdlogphiMat=NULL, asDmLR_invV=NULL) { 
  if (object$spaMM.version<="1.11.60") stop("objects created with spaMM versions <= 1.11.60 are no longer supported.")
  lambda.object <- object$lambda.object
  strucList <- object$strucList
  dwdlogphi <- dwdloglam <- NULL ## always a cbind at the end of calc_logdisp_cov
  dispcolinfo <- list()
  problems <- list() ## Its elements' names are tested in calcPredVar, and the strings are 'development info'
  nrand <- length(strucList)
  col_info <- list(nrand=nrand, ranef_ids=c(), phi_cols=0L)
  if (any(lambda.object$type!="fixed")) {
    corr.models <- lapply(strucList,attr,which="corr.model") ## not unlist bc it may contain NULLs
    checkadj <- unlist(lapply(corr.models,identical,y="adjacency"))
    if(any(checkadj)) {
      ## several blocks of code are "maintained" below for a future dispVar computation for rho
      # il me manque dwdrho (et meme dwdloglam pour ce modele ?) donc on inactive les lignes suivantes:
      #       if (is.null(lambda.object$lambda.fix)) dispnames <- c(dispnames,"loglambda")
      #       corrFixNames <- names(unlist(object$corrPars[which(attr(corrPars,"type")=="fix")]))
      #       if (! ("rho" %in% corrFixNames) ) dispnames <- c(dispnames,"rho")
    }
    checklambda <- ( ! (lambda.object$type %in% c("fixed","fix_ranCoefs"))) 
    if (any(checklambda)) {
      if (is.null(dvdloglamMat)) {
        ## note that .get_logdispObject is computed on request by .get_logdispObject()
        problems$stopmiss <- warning("is.null(dvdloglamMat) in a case where it should be available.") 
      }
      dispcolinfo$loglambda <- "loglambda"
    }
    Xi_cols <- attr(object$ZAlist, "Xi_cols")
    dvdloglam <- matrix(0,ncol=sum(Xi_cols),nrow=NROW(dvdloglamMat))
    cum_Xi_cols <- cumsum(c(0,Xi_cols))
    cum_n_u_h <- attr(lambda.object$lambda,"cum_n_u_h")
    for (randit in which(checklambda)) {
      u.range <- (cum_n_u_h[randit]+1L):(cum_n_u_h[randit+1L])
      Xi_ncol <- Xi_cols[randit]
      u.range <- matrix(u.range,ncol=Xi_ncol)
      for (colit in seq_len(Xi_ncol)) {
        urange <- u.range[,colit]
        dvdloglam[urange,cum_Xi_cols[randit]+colit] <- rowSums(dvdloglamMat[urange,]) ## assuming each lambda_i = lambda in each block
      }
    }
    dwdloglam <- .calc_invL_coeffs(object,dvdloglam) ## one matrix for all ranefs 
    if (!is.null(dwdloglam)) {
      col_info$ranef_ids <- which(rep(checklambda,Xi_cols)) ## (repeated) indices of ranefs, not cols of ranefs
      col_info$cum_n_u_h <- cum_n_u_h
    }
  } 
  ###
  phimodel <- object$models[["phi"]]
  if (phimodel=="phiScal") { ## semble impliquer pas outer phi.Fix... => no need to test object$phi.object$phi_outer,"type")
    phi_est <- object$phi ## no need to get object$phi.object$phi_outer
    if (length(phi_est)!=1L) problems$stopphi <- warning("phimodel=\"phiScal\" but length(phi_est)!=1L.")
    if ( ! is.null(dvdlogphiMat)) {
      dvdlogphi <- rowSums(dvdlogphiMat) ## using each phi_i = phi
      dwdlogphi <- .calc_invL_coeffs(object,dvdlogphi)
      if (!is.null(dwdlogphi)) col_info$phi_cols=length(col_info$ranef_ids)+seq_len(NCOL(dwdlogphi)) ## cols indices for phi 
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
  ## compute info matrix:
  if ((length(dispcolinfo))==0L) {
    return(list(problems=problems))
  } else {
    dwdlogdisp <- cbind(dwdloglam,dwdlogphi) ## typically nobs * 2
    attr(dwdlogdisp,"col_info") <- col_info
    # cf my documentation, based on McCullochSN08 6.62 and 6.74
    # lambda and phi factors enter in dV/dlog(.), computed instead of dV/d(.) to match dwdlog(.) vectors.
    #
    # use repres of two matrices large A and B, each as (thin) lhs %*% (flat) rhs   
    ZAL <- get_ZALMatrix(object)
    if ("loglambda" %in% names(dispcolinfo) || "rho" %in% names(dispcolinfo)) {
      invV.dV_info <- .calc_invV.dV_info(object, checklambda, asDmLR_invV=asDmLR_invV, ZAL=ZAL)
      sublambda <- unlist(invV.dV_info$lambda_list[checklambda])
      dispcolinfo$loglambda <- rep("loglambda",length(sublambda))
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
    if ("rho" %in% dispnames) { ## will occur only when if (any(checkadj)) {...} block above is fixed and active. 
      # no use of sqrt because adjd can be negative
      #invV.dVdrho <- (invV %id*id% ZAL) %*% ( Diagonal(x=lambda*adjd/(denom^2)) %id*id% t(ZAL))
      lhs_invV.dVdrho <- .calc_lhs_invV.dVdlam(object, ZAL, asDmLR_invV) # sweep( ZAL,1L,object$w.resid,`*`) - lhs_invV.dVdrho
      lambda_adjd <- invV.dV_info$lambda_list[[which(checkadj)]] ## asumes single adjd
      rhs_invV.dVdrho <- ( Diagonal(x=lambda_adjd*invV.dV_info$adjd_denom2) %id*id% t(ZAL)) ## FIXME curently meaningful for only one lambda element
      #logdispInfo["rho","rho"] <- sum(invV.dVdrho*t(invV.dVdrho))
      logdispInfo[dispcols$rho,dispcols$rho] <- .traceAB(lhs_invV.dVdrho,rhs_invV.dVdrho,t(rhs_invV.dVdrho),t(lhs_invV.dVdrho))
      if ("loglambda" %in% dispnames) {
        sublambda <- unlist(invV.dV_info$lambda)
        logdispInfoBlock <- numeric(length(sublambda))
        cum_n_u_h <- invV.dV_info$cum_n_u_h
        zerotemplate <- rep(0,cum_n_u_h[nrand+1L])
        for (randit in which(checklambda)) {
          u.range <- (cum_n_u_h[randit]+1L):(cum_n_u_h[randit+1L])
          Xi_ncol <- Xi_cols[randit] # say '1 2' for ranCoefs
          uirange <- matrix(u.range,ncol=Xi_ncol)
          for (ilam in seq_len(Xi_ncol)) { 
            i_rhs_invV.dVdlam <- .fill_rhs_invV.dVdlam(template=zerotemplate, urange=uirange[,ilam], invV.dV_info$rhs_invV.dVdlam)
            colit <- cum_Xi_cols[randit]+ilam
            logdispInfoBlock[colit] <- .traceDB(object$w.resid, t(rhs_invV.dVdrho),t(lhs_invV.dVdrho)) -
              .traceAB(invV.dV_info$lhs_invV.dVdlam, i_rhs_invV.dVdlam,t(rhs_invV.dVdrho),t(lhs_invV.dVdrho))
          }
        }
        logdispInfoBlock <- logdispInfoBlock * sublambda  #lambda * sum(invV.dVdlam*t(invV.dVdrho))
        logdispInfo[dispcols$loglambda,dispcols$rho] <- 
          logdispInfo[dispcols$rho,dispcols$loglambda] <- logdispInfoBlock
      }
    } 
    ## if (! is.null(dwdlogphi)) { ## currently length(phi)==1L && ! is.null(dvdlogphiMat)
    if ("logphi" %in% dispnames) { ## more transparent, but error if mismatch of conditions
      ## next lines assume that  the design matrix for the residual error is I
      logdispInfo[dispcols$logphi,dispcols$logphi] <- phi_est^2 * (
        sum(object$w.resid^2) -2 * .traceDB(object$w.resid,asDmLR_invV$n_x_r, asDmLR_invV$r_x_n) + 
          .traceAB(asDmLR_invV$n_x_r, asDmLR_invV$r_x_n, asDmLR_invV$n_x_r, asDmLR_invV$r_x_n)
      ) # phi_est^2 * sum(invV^2)
      if ("loglambda" %in% dispnames) {
        sublambda <- unlist(invV.dV_info$lambda)
        logdispInfoBlock <- numeric(length(sublambda))
        cum_n_u_h <- invV.dV_info$cum_n_u_h
        zerotemplate <- rep(0,cum_n_u_h[nrand+1L])
        for (randit in which(checklambda)) {
          u.range <- (cum_n_u_h[randit]+1L):(cum_n_u_h[randit+1L])
          Xi_ncol <- Xi_cols[randit] # say '1 2' for ranCoefs
          uirange <- matrix(u.range,ncol=Xi_ncol)
          for (ilam in seq_len(Xi_ncol)) { 
            i_rhs_invV.dVdlam <- .fill_rhs_invV.dVdlam(template=zerotemplate, urange=uirange[,ilam], invV.dV_info$rhs_invV.dVdlam)
            colit <- cum_Xi_cols[randit]+ilam
            logdispInfoBlock[colit] <- .traceDB(object$w.resid, invV.dV_info$lhs_invV.dVdlam, i_rhs_invV.dVdlam) -
              .traceAB(invV.dV_info$lhs_invV.dVdlam, i_rhs_invV.dVdlam, asDmLR_invV$n_x_r, asDmLR_invV$r_x_n)
          }
        }
        logdispInfoBlock <- logdispInfoBlock * sublambda * phi_est # lambda * phi_est * sum(invV.dVdlam * invV)
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
    if (inherits(logdisp_cov,"try-error")) logdisp_cov <- ginv(logdispInfo) ## quick patch for uninteresting case
    if (any( ! checklambda )) { ## if cols missing from logdisp_cov compared to dwdlogdisp
      ncd <- ncol(dwdlogdisp)
      full_logdisp_cov <- matrix(0,ncd,ncd)
      cols_in_logdisp_cov <- c(rep(checklambda,Xi_cols), ## which cols in dwdloglam match loglambda col in logdisp_cov
                               rep(TRUE,max(ncol(dwdlogphi),0)))  ## cols of dwdlogphi
      full_logdisp_cov[cols_in_logdisp_cov,cols_in_logdisp_cov] <- logdisp_cov
      return(list(dwdlogdisp=dwdlogdisp,logdisp_cov=full_logdisp_cov,problems=problems)) 
    } else return(list(dwdlogdisp=dwdlogdisp,logdisp_cov=logdisp_cov,problems=problems)) 
    ## more compact than storing ww %*% logdisp_cov %*% t(ww) which is nobs*nobs 
  }
}

