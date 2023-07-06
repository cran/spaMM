# evaluate tr(A %*% B)= sum(A*B) where A and B are large matrices but each of the form l %*% r for 'narrow' l
# this function avoids the formation of the large 'n x n' matrices, using a form of commutation of trace arguments.
# However the computation of 'r x r' matrix by the crossproducts may still be quite long
# e.g. 792*62388 matrices => (792^2)*62388 = 39 133 746 432 flop'product' and quite a few 'sum'
# That's something that 'should' be easy to paralellise
.traceAB <- function(lA,rA,lB,rB, B_is_tA=FALSE) {
  if (nrow(lA)>ncol(lA)) { # lA is typically nXr so this occurs when more obs than ranef levels
    if (B_is_tA) {
      rr <- rA %*% lA # not A= lA %*% rA
      return(sum(t(rr)*rr)) ## not sum(rr^2) which is the result when B=A (as used below)
    } else { # more subtly handling the case of NULL B than in  'more ranefs' case (the crossprds are trivial)
      if (.spaMM.data$options$Matrix_old) { # ugly... but such versions do not handle as(, "generalMatrix"))
        ll <- .crossprod(lA, lB)
        rr <- .tcrossprod(rA, rB) ### dsC or dpo or...?    # slower if *both* matrices have sparse storage though being dense
      } else {
        ll <- as(.crossprod(lA, lB),"generalMatrix")
        rr <- as(.tcrossprod(rA, rB),"generalMatrix") ## slower if *both* matrices have sparse storage though being dense
      }
      return(sum(ll*rr)) # elementwise product of dsC if no as(.,"generalMatrix"). Matrix v1.4-2 might complain.
    }
  } else {
    A <- lA %*% rA
    if (B_is_tA) {
      return(sum(A*t(A))) 
    } else if (is.null(lB) && is.null(rB)) {  # meaning that B=A
      return(sum(A^2)) # sum(diag(tcrossprod(A)))
    } else {
      B <- lB %*% rB 
      return(sum(A*B)) # sum(diag(tcrossprod(A,B)))
    }
  }
}

.gmp_traceAB <- function(lA,rA,lB,rB) {
  lA <- gmp::as.bigq(as.matrix(lA))
  rA <- gmp::as.bigq(as.matrix(rA))
  lB <- gmp::as.bigq(as.matrix(lB))
  rB <- gmp::as.bigq(as.matrix(rB))
  ll <- gmp::crossprod(lA, lB)
  rr <- gmp::tcrossprod(rA, rB)
  return(as.numeric(gmp::sum.bigq(ll*rr)))
}



# and same concept for trace( D %*% B)
.traceDB <- function(dD,lB,rB) sum(.Dvec_times_m_Matrix(dD, lB * t(rB))) #  { sum( sweep(lB * t(rB),1L,dD,`*`) )}
# lA <- matrix(runif(6),ncol=2)
# rA <- matrix(runif(6),ncol=3)
# lB <- matrix(runif(6),ncol=2)
# rB <- matrix(runif(6),ncol=3)
# sum((lA %*% rA) * (lB %*% rB))
# .traceAB(lA,rA,lB,rB)
# dD <- runif(3)
# sum(diag(diag(dD) %*% lB %*% rB))
# .traceDB(dD,lB,rB)

.calc_lhs_invV.dVdlam <- function(object, ZALd, invV_factors) { 
  if (.is_spprec_fit(object)) { # alternative code is always valid! BUT....
    # dgeMatrix is more efficient in products using the result 'lhs_invV.dVdlam':
    # Further, we have dgeMatrices where the alternative code produces dgCMatrices, 
    # W_ZinvG_ZtW_ZA <- invV_factors$n_x_r =W Z       
    #                                  %*% invV_factors$r_x_n=Matrix::solve(object$envir$G_CHMfactor, object$envir$ZtW) 
    #                                                          %*% ZA
    if (FALSE) {
      ZAfix <- .get_ZAfix(object)
      W_ZinvG_ZtW_ZA <- invV_factors$n_x_r %*% as(invV_factors$r_x_n %*% ZAfix, "dgeMatrix")
    } else W_ZinvG_ZtW_ZA <- invV_factors$n_x_r %*% invV_factors$r_x_r  # precomput r_x_r controls its type relative to $r_x_n
    lhs_invV.dVdlam <- invV_factors$n_x_r - W_ZinvG_ZtW_ZA # (w.resid- n_x_r %*% r_x_n) %*% ZA 
  } else if (missing(ZALd)) { # To produce "iVZA|L" factorization 
    # ZALd missing eiher for spprec (above or TRY_dense_iVZA (FALSE)); r_x_r used only in these cases
    W_ZinvG_ZtW_ZA <- invV_factors$n_x_r %*% invV_factors$r_x_r  # precomput r_x_r controls its type relative to $r_x_n
    lhs_invV.dVdlam <- invV_factors$n_x_r - W_ZinvG_ZtW_ZA # (w.resid- n_x_r %*% r_x_n) %*% ZA 
  } else {
    ## next lines use invV= w.resid- n_x_r %*% r_x_n
    # W_ZinvG_ZtW_ <- invV_factors$n_x_r =W Z       
    #                                  %*% invV_factors$r_x_n=Matrix::solve(object$envir$G_CHMfactor, object$envir$ZtW) 
    #                                                          %*% ZALd
    W_ZinvG_ZtW_ <- invV_factors$n_x_r %*% (invV_factors$r_x_n %*% ZALd)  
    lhs_invV.dVdlam <- .Dvec_times_m_Matrix(.get_H_w.resid(object), ZALd) - W_ZinvG_ZtW_ # (w.resid- n_x_r %*% r_x_n) %*% ZALd 
  }
  
  return(lhs_invV.dVdlam)
}  

.calc_invV.dV_info <- function(object,checklambda,ZAL,
                               invV_factors
                               ) {
  strucList <- object$strucList
  lambda.object <- object$lambda.object
  cum_n_u_h <- attr(lambda.object$lambda_list,"cum_n_u_h")
  lambda_list <- vector("list",length(strucList))
  RES <- list(cum_n_u_h=cum_n_u_h)
  ZAL_to_ZALd_vec <- rep(0,cum_n_u_h[length(cum_n_u_h)]) 
  for (randit in which(checklambda)) {
    lmatrix <- strucList[[randit]]
    corr.model <- attr(lmatrix,"corr.model")
    u.range <- (cum_n_u_h[randit]+1L):(cum_n_u_h[randit+1L])
    if (identical(corr.model,"adjacency") &&
        ( ! is.null(glmit <- lambda.object$rand_to_glm_map[randit])) && ## means inner estimation of lambda
        "adjd" %in% lambda.object$print_namesTerms[[randit]] ## means inner estimation of adj_rho
       ) { 
      ## we end here if not sparse (sparse => outer optim of rho) and inner optim of lambda was allowed
      lambda_list[[randit]] <- with(lambda.object,linkinvS[[glmit]](coefficients_lambdaS[[randit]][1]))
      RES$rho <- - with(lambda.object,coefficients_lambdaS[[randit]][2]/coefficients_lambdaS[[randit]][1])
      adjd <- attr(lmatrix,"symsvd")$adjd ## svd not SVD because went through mat_sqrt which changed attributes
      denom <- 1-RES$rho*adjd
      RES$adjd_denom2 <- adjd/(denom^2) 
      ZAL_to_ZALd_vec[u.range] <- 1/sqrt(denom)
    } else { ## (1) standard lamScal model or ranCoefs model, as results from outer optimization; lmatrix has no attributes; or
             ## (2) For adjacency case, sparse_precision  forces outer optim of rho and thus we end here too
      ZAL_to_ZALd_vec[u.range] <- 1
      coeff <- lambda.object$coefficients_lambdaS[[randit]] ## corrHLfit or fitme with control$refit=TRUE
      # there are several coeffs, computed in .bloc_lambda(), for ranCoefs terms.
      if (is.null(coeff) || anyNA(coeff)) { # NA for multIMRF and NULL for other outer optim 
        lambda_list[[randit]] <- lambda.object$lambda_list[[randit]] ## basic fitme (may numerically differ)
      } else lambda_list[[randit]] <- exp(coeff) # inner estimation.
    }
  }
  RES$lambda_list <- lambda_list
  if (.is_spprec_fit(object)) {
    RES$lhs_invV.dVdlam <- .calc_lhs_invV.dVdlam(object, invV_factors=invV_factors) ## invV %*% ZA
    RES$envir <- object$envir # to use $chol_Q and $ZAfix without any copy here 
    RES$ZAL_to_ZALd_vec <- ZAL_to_ZALd_vec
    RES$type <- "|L" # "iVZA | (Ld!dL)AZ" 
    ## traceAB will use iVZA as lhs and Ld!dL'A'Z' as rhs: we store iVZA and A'Z and .fill_rhs_invV.dVdlam() factors by (Ld!dL)
  } else {
    if (identical(.spaMM.data$options$TRY_dense_iVZA,TRUE)) { # (ZALd argument missing in next call) # only for deve purpose
      # I may not yet have good tests of the efficiency of this code, but I keep it ready for use ## F I X M E Not convincing. Retry
      RES$lhs_invV.dVdlam <- .calc_lhs_invV.dVdlam(object, invV_factors=invV_factors) ## invV %*% ZA 
      #
      ZAphant <- object$ZAlist
      for (rd in seq_along(ZAphant)) ZAphant[[rd]] <- Diagonal(n=ncol(ZAphant[[rd]]))
      ZAXlist <- .compute_ZAXlist(XMatrix=object$strucList, ZAlist=ZAphant, force_bindable=TRUE)
      RES$Lmatrix <- do.call(Matrix::bdiag,ZAXlist) # it would be nice to avoid this
      #
      .get_ZAfix(object) ## makes sure that ZAfix is in the object's environment
      RES$envir <- object$envir # to use $ZAfix without any copy here
      RES$ZAL_to_ZALd_vec <- ZAL_to_ZALd_vec
      RES$type <- "iVZA|L" # "iVZA | (Ld!dL)AZ" 
    } else {
      ZALd <- .m_Matrix_times_Dvec(ZAL, ZAL_to_ZALd_vec)
      # in invV.dVdlam the rhs of dvdlam = ZALd %*% t(ZALd) is split between lhs_invV.dVdlam and rhs_invV.dVdlam
      #                                                               [               rhs             ]       lhs
      # invV.dVdlam = (w.resid- n_x_r %*% r_x_n).(ZALd %*% t(ZALd)) = ((w.resid- n_x_r %*% r_x_n).(ZALd) %*% t(ZALd)
      RES$lhs_invV.dVdlam <- .calc_lhs_invV.dVdlam(object, ZALd, invV_factors) ## invV %*% ZALd
      RES$rhs_rhs_invV.dVdlam <- t(ZALd)
      RES$type <- "L|"        # "iVZALd | (!)dLAZ" 
      ## traceAB will use iVZALd as lhs and !dL'A'Z' as rhs: we store iVZALd and dL'A'Z' and .fill_rhs_invV.dVdlam() factors by (!).
    }
    return(RES)
  }
  return(RES) 
}

.fill_rhs_invV.dVdlam <- function(template, urange, invV.dV_info) { ## to localise template and urange
  template[urange] <- 1L
  if (invV.dV_info$type=="|L") { # only spprec
    if ( ! is.null(latent_d_list <- invV.dV_info$envir$sXaug$AUGI0_ZX$envir$latent_d_list)) {
      chol_Q_w <- .Matrix_times_Dvec(t(invV.dV_info$envir$chol_Q),1/sqrt(.unlist(latent_d_list)))
    } else chol_Q_w <- t(invV.dV_info$envir$chol_Q)
    # ZAL_to_ZALd_vec is a correction for adjacency not for ranCoefs
    tcrossfac <- solve(chol_Q_w, Diagonal(x=invV.dV_info$ZAL_to_ZALd_vec * template)) # Ld
    lhs <- .tcrossprod(tcrossfac) # LddL'
    return(as.matrix(.tcrossprod(lhs, invV.dV_info$envir$sXaug$AUGI0_ZX$ZAfix))) # LddL'A'Z'
    #lhs <- .ZWZtwrapper(invV.dV_info$LMatrix, (invV.dV_info$ZAL_to_ZALd_vec^2 * template)) 
    #return(.tcrossprod(lhs, invV.dV_info$ZAfix)) ## effectively dense lhs => slow
  } else if (invV.dV_info$type=="iVZA|L") { # this is only devel code
    warning("devel code not tested after redef of structList for ranCoefs")
    # => now, this seems correct
    # Further, we won't see any pb as long as the ranCoefs $d are 1, which is so if *!*spprec and chol() worked in .calc_latentL()
    # Thus this devel code must be OK 
    tcrossfac <- .m_Matrix_times_Dvec(invV.dV_info$Lmatrix, invV.dV_info$ZAL_to_ZALd_vec * template) # Ld
    lhs <- .tcrossprod(tcrossfac) # LddL'
    return(as.matrix(.tcrossprod(lhs, invV.dV_info$envir$sXaug$AUGI0_ZX$ZAfix))) # LddL'A'Z'
  } else {
    # (FIXME) I could add a drop0 when there are several lambda's and matrices are sparse
    return(.Dvec_times_m_Matrix(template, invV.dV_info$rhs_rhs_invV.dVdlam)) # sweep( invV.dV_info$rhs_invV.dVdlam,1L,iloc,`*`))
  }
}

.calc_loglamInfo <- function(invV.dV_info,which) { ## called by .calc_logdisp_cov(), using the result of .calc_invV.dV_info()
  lambda_list <- invV.dV_info$lambda_list ## for ranCoefs, must be a list with Xi_cols-matching elements 
  Xi_cols <- sapply(lambda_list,length)
  n_sublambda <- sum(Xi_cols[which]) 
  loglamInfo <- matrix(ncol=n_sublambda,nrow=n_sublambda)
  cum_n_u_h <- invV.dV_info$cum_n_u_h
  zerotemplate <- rep(0,cum_n_u_h[length(cum_n_u_h)]) 
  cum_Xi_cols <- cumsum(c(0,Xi_cols))
  rhs_invV.dVdlam_list <- list()
  for (randit in which) { # say '2' for second ranef
    u.range <- (cum_n_u_h[randit]+1L):(cum_n_u_h[randit+1L])
    Xi_ncol <- Xi_cols[randit] # say '1 2' for ranCoefs
    uirange <- matrix(u.range,ncol=Xi_ncol)
    for (ilam in seq_len(Xi_ncol)) { 
      i_rhs_invV.dVdlam <- .fill_rhs_invV.dVdlam(template=zerotemplate, urange=uirange[,ilam], invV.dV_info)
      rhs_invV.dVdlam_list[[paste0(randit,"_",ilam)]] <- i_rhs_invV.dVdlam
      colit <- cum_Xi_cols[randit]+ilam
      loglamInfo[colit,colit] <- .traceAB(lA=invV.dV_info$lhs_invV.dVdlam, rA=i_rhs_invV.dVdlam,
                                          #lB=t(i_rhs_invV.dVdlam), rB=t(invV.dV_info$lhs_invV.dVdlam)
                                          B_is_tA=TRUE # (lA is n x r, rA is r x n, so A is square as this cases requires)
                                          )
      for (jlam in seq_len(ilam-1L)) { ## WITHIN a randit
        j_rhs_invV.dVdlam <-  rhs_invV.dVdlam_list[[paste0(randit,"_",jlam)]] #.fill_rhs_invV.dVdlam(template=zerotemplate, urange=uirange[,jlam], invV.dV_info)
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
          j_rhs_invV.dVdlam <- rhs_invV.dVdlam_list[[paste0(randjt,"_",jlam)]] #.fill_rhs_invV.dVdlam(template=zerotemplate, urange=ujrange[,jlam], invV.dV_info)
          coljt <- cum_Xi_cols[randjt]+jlam
          loglamInfo[coljt,colit] <- .traceAB(invV.dV_info$lhs_invV.dVdlam,i_rhs_invV.dVdlam, 
                                              t(j_rhs_invV.dVdlam),t(invV.dV_info$lhs_invV.dVdlam))
          loglamInfo[colit,coljt] <- loglamInfo[coljt,colit]
        }
      }
    }
  }
  sub_lambda_vec <- .unlist(lambda_list[which]) 
  loglamInfo <- loglamInfo * (sub_lambda_vec %*% t(sub_lambda_vec))
  return(list(loglamInfo=loglamInfo,rhs_invV.dVdlam_list=rhs_invV.dVdlam_list))
}

.wrap_solve_logdispinfo <- function(logdispInfo, object) {
  logdisp_cov <- try(solve(logdispInfo), silent = TRUE)
  problem <- inherits(logdisp_cov, "try-error")
  if (!problem) {
    problem <- any(diag(logdisp_cov) < 0) # strngly suggest major inaccuracy or bug in computing the matrix
    if (problem) {
      warning(paste("Numerical precision issue or something else?\n", 
                    "  Information matrix for dispersion parameters does not seem positive-definite.\n", 
                    "  The prediction variance may be inaccurate."))
    }
  } else { # solve() failed
    lambdas <- VarCorr(object, add_residVars=FALSE)[,"Variance"]
    if (any(lambdas<1e-6)) { ## may not be appropriate for all models (cf non-gaussian ranefs?)
      # we could message() that a ranef has low variance, but why do that here?
      # .force_solve will use regularization
      # regularization should lead to underestimate the logdisp_cov where it is tiny, so it should be OK
      # *If* there is still a potential for inaccuracies, this should be dealt at the fix_predVar level ?
    } else if (all(lambdas>1e-4) && any(eigen(logdispInfo, only.values=TRUE)$values==0)) { # suggests exactly singular matrix i.e. redundant ranefs (twolambda example)
      message("Suspiciously-looking information matrix for dispersion parameters: maybe redundant random effects?")
    } else warning(paste("Numerical precision issue in computation of the information matrix for dispersion parameters:\n", 
                         "  the prediction variance may be inaccurate."))
  }
  if (problem) {
    logdisp_cov <- .force_solve(logdispInfo)
  }
  return(list(logdisp_cov=logdisp_cov, problem=problem))
}

.calc_logdisp_cov <- function(object, dvdloglamMat=NULL, dvdlogphiMat=NULL, invV_factors=NULL,
                              force_fixed=FALSE # if TRUE -> force computation for fixed pars too 
                              ) { 
  if (object$spaMM.version<="1.11.60") stop("objects created with spaMM versions <= 1.11.60 are no longer supported.")
  lambda.object <- object$lambda.object
  strucList <- object$strucList
  dispcolinfo <- list()
  problems <- list() ## Its elements' names are tested in calcPredVar, and the strings are 'development info'
  nrand <- length(strucList)
  col_info <- list(nrand=nrand, phi_cols=NULL) 
  Xi_cols <- attr(object$ZAlist, "Xi_cols")
  dwdloglam <- matrix(0,ncol=sum(Xi_cols),nrow=length(object$v_h)) # cols will remain 0 for fixed lambda params
  if (force_fixed) {
    checklambda <- rep(TRUE, length(lambda.object$type))
  } else checklambda <- ( ! (lambda.object$type %in% c("fixed","fix_ranCoefs","fix_hyper"))) 
  if (any(checklambda)) {
    exp_ranef_types <- attr(object$ZAlist,"exp_ranef_types")
    checkadj <- (exp_ranef_types=="adjacency")
    if(any(checkadj)) {
      ## several blocks of code are "maintained" below for a future dispVar computation for rho
      # il me manque dwdrho (et meme dwdloglam pour ce modele ?) donc on inactive les lignes suivantes:
      #       if (is.null(lambda.object$lambda.fix)) dispnames <- c(dispnames,"loglambda")
      #       corrFixNames <- names(unlist(object$corrPars[which(attr(corrPars,"type")=="fix")]))
      #       if (! ("rho" %in% corrFixNames) ) dispnames <- c(dispnames,"rho")
    }
    
    if (is.null(dvdloglamMat)) {
      ## note that .get_logdispObject is computed on request by .get_logdispObject()
      problems$stopmiss <- warning("is.null(dvdloglamMat) in a case where it should be available.") 
    }
    dispcolinfo$loglambda <- "loglambda"
    #dvdloglam <- matrix(0,nrow=NROW(dvdloglamMat), ncol=sum(Xi_cols))
    strucList <- object$strucList
    cum_n_u_h <- attr(lambda.object$lambda_list,"cum_n_u_h")
    n_u_h <- diff(cum_n_u_h)
    cum_Xi_cols <- cumsum(c(0,Xi_cols))
    ## dwdloglam will include cols of zeros for fixed lambda; matching with reduced logdisp_cov is performed at the end of the function.
    for (randit in seq_len(nrand)) { ## ALL ranefs!
      range_in_dw <- (cum_n_u_h[randit]+1L):(cum_n_u_h[randit+1L])
      if ( inherits(strucList[[randit]],"dCHMsimpl")) { # e.g. get_predVar(adjfitsp)
        for_dw_i <- as(strucList[[randit]], "CsparseMatrix") %*% dvdloglamMat[range_in_dw,  ] # i.e L_Q %*% lignes de (t(L_Q) %*% invG %*% L_Q %*% some rhs) 
      } else if ( ! is.null(lmatrix <- strucList[[randit]])) {
        for_dw_i <- solve(t(lmatrix),dvdloglamMat[range_in_dw,]) ## f i x m e for efficiency ? store info about solve(t(lmatrix)) in object ? 
      } else { ## implicit identity lmatrix
        for_dw_i <- dvdloglamMat[range_in_dw,] ## assuming each lambda_i = lambda in each block
      }
      nblocks_randit <- Xi_cols[randit]
      rowranges_in_dw_i <- matrix(seq(n_u_h[randit]),ncol=nblocks_randit) ## this _splits_ seq(n_u_h[randit]) over two columns for a random-slope model
      for (row_block in seq_len(nblocks_randit)) { ## half-ranges for random-slope model
        rowrange_in_dw_i <- rowranges_in_dw_i[,row_block]
        cum_rowrange_in_dw <- rowrange_in_dw_i + cum_n_u_h[randit]
        for (randjt in which(checklambda)) { ## NOT all ranefs! # iteraction over ranefs, not lambdas ranCoef: one ranef, several lambdas)
          nblocks_randjt <- Xi_cols[randjt]
          cum_colrange_in_dw_i <- (cum_n_u_h[randjt]+1L):(cum_n_u_h[randjt+1L])
          cum_colranges_in_dw_i <- matrix(cum_colrange_in_dw_i,ncol=nblocks_randjt) ## this _splits_ seq(n_u_h[randit]) over two columns for a random-slope model
          for (col_in_colranges_dw_i in seq(nblocks_randjt)) { ## half-ranges for random-slope model; trivial bug here detected 08/2021 thanks to test-composite.R
            cum_col_in_dw <- cum_Xi_cols[randjt]+col_in_colranges_dw_i
            cum_cols_in_dw_i <- cum_colranges_in_dw_i[,col_in_colranges_dw_i] 
            dwdloglam[cum_rowrange_in_dw, cum_col_in_dw] <- rowSums(for_dw_i[rowrange_in_dw_i, cum_cols_in_dw_i,drop=FALSE])  
          }
        }
      }
    }
    ## dwdloglam includes cols of zeros for fixed lambda; matching with reduced logdisp_cov is performed at the end of the function.
    ranef_ids <- rep(seq_len(nrand),Xi_cols) ## (repeated for ranCoefs) indices of ranefs, not cols of ranefs
  } else ranef_ids <- NULL
  ###
  phimodel <- object$models[["phi"]]
  if (force_fixed ||
      phimodel=="phiScal" ||  
      identical(object$envir$forcePhiComponent,TRUE) ## old hack for code testing: force dispVar computation as if phi was not fixed.
      ) { ## semble impliquer pas outer phi.Fix... => no need to test object$phi.object$phi_outer,"type")
    phi_est <- object$phi ## no need to get object$phi.object$phi_outer
    if (length(phi_est)!=1L) problems$stopphi <- warning("phimodel=\"phiScal\" but length(phi_est)!=1L.")
    if ( ! is.null(dvdlogphiMat)) {
      dvdlogphi <- rowSums(dvdlogphiMat) ## using each phi_i = phi # always a vector, even from 0-col matrix
      # => r-vector over r v's with element k = sum_over_responses (dv_k/d phi_i)
      # for distinct phi_i we would have a correspondingly expanded information matrix 
      # and then dvdlogphi and dwdlogphi would retain distinct columns
      dwdlogphi <- .calc_invL_coeffs(object,dvdlogphi) # input always a vector, output always a vector
      col_info$phi_cols <- length(ranef_ids)+1L ## cols indices for phi 
      dispcolinfo$logphi <- "logphi"
    } else if (object$models[["eta"]]=="etaHGLM") stop("phimodel=='phiScal' but is.null(dvdlogphiMat)")
  } else {  ## else phimodel="", e.g. binomial
    # if binomial or poisson, phimodel=""; warning for other phimodels
    if (any(phimodel!="")) {
      problems$structphi <- "phi dispVar component not yet available for phi model != ~1."
      if ( ! identical(spaMM.getOption("phi_dispVar_comp_warned"),TRUE)) {
        warning(problems$structphi)
        .spaMM.data$options$phi_dispVar_comp_warned <- TRUE
      }
    }
    dwdlogphi <- NULL
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
    ZAL <- get_ZALMatrix(object, force_bind = ! (.is_spprec_fit(object)) )
    if ("loglambda" %in% names(dispcolinfo) || "rho" %in% names(dispcolinfo)) {
      invV.dV_info <- .calc_invV.dV_info(object, checklambda, invV_factors=invV_factors, ZAL=ZAL) ## $lhs= invV %*% ZALd and $lhs= t(ZALd)
      sublambda <- .unlist(invV.dV_info$lambda_list[checklambda])
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
    colnames(logdispInfo) <- rownames(logdispInfo) <- .unlist(dispcolinfo)
    if ("loglambda" %in% dispnames) { 
      loglamInfo_blob <- .calc_loglamInfo(invV.dV_info,which=which(checklambda))
      logdispInfo[dispcols$loglambda,dispcols$loglambda] <- loglamInfo_blob$loglamInfo 
    }
    if ("rho" %in% dispnames) { ## will occur only when if (any(checkadj)) {...} block above is fixed and active. 
      # no use of sqrt because adjd can be negative
      #invV.dVdrho <- (invV %id*id% ZAL) %*% ( Diagonal(x=lambda*adjd/(denom^2)) %id*id% t(ZAL))
      lhs_invV.dVdrho <- .calc_lhs_invV.dVdlam(object, ZAL, invV_factors) # sweep( ZAL,1L,object$w.resid,`*`) - lhs_invV.dVdrho
      lambda_adjd <- invV.dV_info$lambda_list[[which(checkadj)]] ## asumes single adjd
      rhs_invV.dVdrho <- ( Diagonal(x=lambda_adjd*invV.dV_info$adjd_denom2) %id*id% t(ZAL)) ## FIXME curently meaningful for only one lambda element
      #logdispInfo["rho","rho"] <- sum(invV.dVdrho*t(invV.dVdrho))
      logdispInfo[dispcols$rho,dispcols$rho] <- .traceAB(lhs_invV.dVdrho,rhs_invV.dVdrho,t(rhs_invV.dVdrho),t(lhs_invV.dVdrho))
      if ("loglambda" %in% dispnames) {
        sublambda <- .unlist(invV.dV_info$lambda)
        logdispInfoBlock <- numeric(nrand)
        cum_n_u_h <- invV.dV_info$cum_n_u_h
        zerotemplate <- rep(0,cum_n_u_h[nrand+1L])
        for (randit in which(checklambda)) {
          u.range <- (cum_n_u_h[randit]+1L):(cum_n_u_h[randit+1L])
          Xi_ncol <- Xi_cols[randit] # say '1 2' for ranCoefs
          uirange <- matrix(u.range,ncol=Xi_ncol)
          for (ilam in seq_len(Xi_ncol)) { 
            i_rhs_invV.dVdlam <- loglamInfo_blob$rhs_invV.dVdlam_list[[paste0(randit,"_",ilam)]]  #.fill_rhs_invV.dVdlam(template=zerotemplate, urange=uirange[,ilam], invV.dV_info)
            colit <- cum_Xi_cols[randit]+ilam
            logdispInfoBlock[colit] <- .traceDB(.get_H_w.resid(object), t(rhs_invV.dVdrho),t(lhs_invV.dVdrho)) -
              .traceAB(invV.dV_info$lhs_invV.dVdlam, i_rhs_invV.dVdlam,t(rhs_invV.dVdrho),t(lhs_invV.dVdrho))
          }
        }
        logdispInfoBlock <- logdispInfoBlock[which(checklambda)] * sublambda  #lambda * sum(invV.dVdlam*t(invV.dVdrho))
        logdispInfo[dispcols$loglambda,dispcols$rho] <- 
          logdispInfo[dispcols$rho,dispcols$loglambda] <- logdispInfoBlock
      }
    } 
    ## if (! is.null(dwdlogphi)) { ## currently length(phi)==1L && ! is.null(dvdlogphiMat)
    if ("logphi" %in% dispnames) { ## more transparent, but error if mismatch of conditions
      ## next lines assume that  the design matrix for the residual error is I
      # using the pattern (D-nXr.rXn)^2 = D^2 - 2 D nXr.rXn + (nXr.rXn)^2
      if (.is_spprec_fit(object)) {
        #A <- solve(object$envir$G_CHMfactor, .tcrossprod(object$envir$ZtW), system="A")
        A <- invV_factors$r_x_n %*% invV_factors$n_x_r
        trAB <- sum(A^2)
        H_w.resid <- .get_H_w.resid(object)
        logdispInfo[dispcols$logphi,dispcols$logphi] <- phi_est^2 * (
          sum(H_w.resid^2) -2 * .traceDB(H_w.resid,invV_factors$n_x_r, invV_factors$r_x_n) + 
            trAB
        ) # phi_est^2 * sum(invV^2)
      } else {
        H_w.resid <- .get_H_w.resid(object)
        logdispInfo[dispcols$logphi,dispcols$logphi] <- phi_est^2 * (
          sum(H_w.resid^2) -2 * .traceDB(H_w.resid,invV_factors$n_x_r, invV_factors$r_x_n) + 
            .traceAB(lA=invV_factors$n_x_r, rA=invV_factors$r_x_n, lB=NULL, rB=NULL) # ie B=A
        ) # phi_est^2 * sum(invV^2)
      }
      if ("loglambda" %in% dispnames) {
        sublambda <- .unlist(invV.dV_info$lambda)
        logdispInfoBlock <- numeric(nrand)
        cum_n_u_h <- invV.dV_info$cum_n_u_h
        zerotemplate <- rep(0,cum_n_u_h[nrand+1L])
        for (randit in which(checklambda)) {
          u.range <- (cum_n_u_h[randit]+1L):(cum_n_u_h[randit+1L])
          Xi_ncol <- Xi_cols[randit] # say '1 2' for ranCoefs
          uirange <- matrix(u.range,ncol=Xi_ncol)
          for (ilam in seq_len(Xi_ncol)) { 
            i_rhs_invV.dVdlam <- loglamInfo_blob$rhs_invV.dVdlam_list[[paste0(randit,"_",ilam)]] #.fill_rhs_invV.dVdlam(template=zerotemplate, urange=uirange[,ilam], invV.dV_info)
            colit <- cum_Xi_cols[randit]+ilam
            # sum(diag(diag(w.resid) %*% lhs_invV.dVdlam)) - sum(diag(lhs_invV.dVdlam %*% invV_factors))
            logdispInfoBlock[colit] <- .traceDB(.get_H_w.resid(object), invV.dV_info$lhs_invV.dVdlam, i_rhs_invV.dVdlam) -
              .traceAB(invV.dV_info$lhs_invV.dVdlam, i_rhs_invV.dVdlam, invV_factors$n_x_r, invV_factors$r_x_n) 
            # The tcrossprod (                       i_rhs_invV.dVdlam,                     invV_factors$r_x_n) is a bottleneck for large (n>r)
            # while crossprod(           .                              lB=object$envir$ZtW                   ) has quite sparse lB
                      }
        }
        logdispInfoBlock <- logdispInfoBlock[which(checklambda)] * sublambda * phi_est # lambda * phi_est * sum(invV.dVdlam * invV)
        logdispInfo[dispcols$loglambda,dispcols$logphi] <- 
          logdispInfo[dispcols$logphi,dispcols$loglambda] <- logdispInfoBlock
      }
      if ("rho" %in% dispnames) {
        logdispInfo[dispcols$rho,dispcols$logphi] <- 
          logdispInfo[dispcols$logphi,dispcols$rho] <- phi_est * .traceAB(lhs_invV.dVdrho,rhs_invV.dVdrho, 
                                                                          invV_factors$n_x_r, invV_factors$r_x_n)  
        # phi_est * sum(invV.dVdrho * invV)  
      }
    } 
    logdispInfo <- logdispInfo/2
    resu <- .wrap_solve_logdispinfo(logdispInfo, object)
    if (any( ! checklambda )) { ## if (lambda) cols missing from logdisp_cov compared to dwdlogdisp
      #ncd <- ncol(dwdlogdisp)
      #full_logdisp_cov <- matrix(0,ncd,ncd)
      # : alternative ncd below should be equivalent but more self contained
      cols_in_logdisp_cov <- rep(checklambda,Xi_cols) ## which cols in dwdloglam match loglambda col in logdisp_cov
      if ( ! is.null(dwdlogphi)) cols_in_logdisp_cov <- c(cols_in_logdisp_cov,TRUE)  ## col for dwdlogphi
      ncd <- length(cols_in_logdisp_cov)
      full_logdisp_cov <- matrix(0,ncd,ncd)
      full_logdisp_cov[cols_in_logdisp_cov,cols_in_logdisp_cov] <- resu$logdisp_cov
      resu$logdisp_cov <- full_logdisp_cov
    }  
    # Example(s):
    #
    # mrf <- fitme(migStatus ~ 1 + (1|pos) + multIMRF(1|longitude+latitude,margin=2,levels=2, coarse=4) ...
    # with fixed lambda for (1|pos). logdispInfo will be 3*3, 3 being the number of parameters whether fixed or not 
    #  [one lambda for (.|.), one lambda hyper-par for IMRF, one phi]
    # but full_logdisp_cov is 4*4, (from IMRF(...levels=2)), with a summingMat later operating the match.
    # The first row and column of full_logdisp_cov remain zero bc the first lambda is fixed.
    #
    resu$dwdlogdisp <- dwdlogdisp
    return(resu)
    ## more compact than storing ww %*% logdisp_cov %*% t(ww) which is nobs*nobs 
  }
}
