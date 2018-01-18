`%id*%` <- function(A,b) {
  if (.is_identity(A)) return(b)
  return(A %*% b)
}

`%*id%` <- function(a,B) {
  if (.is_identity(B)) return(a)
  return(a %*% B)
}


`%id*id%` <- function(A,B) {
  if (.is_identity(A)) return(B)
  if (.is_identity(B)) return(A)
  return(A %*% B)
}

.calcZWZt_mat_or_diag <- function(Z,W,returnMat) { ## fixefVar or fixefVar + a bit of ranefVar
  ## bottleneck for large predVar with large n_u_h (W=beta_w_cov being a matrix)
  if (returnMat) {
    return(Z[] %id*id% W[] %id*id% t(Z)[]) ## not ZWZT functions bc W not diag a priori if (returnMat)
  } else { ## returns only the diagonal
    if (is.vector(W)) {
      return(rowSums(.Matrix_times_Dvec(Z[],W[]) * Z[]))
    } else {
      rownames(W) <- colnames(W) <- NULL ## inhibits a useless warning from Matrix:: dimNamesCheck 
      premul <- Z[] %id*id% W[]
      return(rowSums(suppressMessages(premul * Z[]))) ## suppress message("method with signature...") [found by debug(message)]
    }
  }
}

.calc_Evar <- function(newZAlist, newinold, cov_newLv_oldv_list, lambda, invCov_oldLv_oldLv_list, covMatrix,## arguments for all calls
                       # arguments for standard diag or full Evar matrix:
                       cov_newLv_newLv_list, 
                       ## arguments for non-diagonal block of Evar matrix:
                       cov_newLv_fixLv_list, cov_fixLv_oldv_list, fixZAlist=NULL 
) {
  newnrand <- length(newinold)
  Evarlist <- vector("list",newnrand)
  for (new_rd in seq_along(Evarlist)) {
    old_rd <- newinold[new_rd]
    isEachNewLevelInOld <- attr(cov_newLv_oldv_list[[new_rd]],"isEachNewLevelInOld") ## non autocorr effect: a vector of booleans indicating whether new is in old (qualifies cols of Cnewold)
    if ( ! is.null(isEachNewLevelInOld)) { ## non spatial effect: a vector of booleans indicating whether new is in old (qualifies cols of Cnewold)
      Evar <- Diagonal(x=lambda[old_rd]*as.numeric(! isEachNewLevelInOld))
      terme <- .calcZWZt_mat_or_diag(newZAlist[[new_rd]], Evar, covMatrix)
    } else { ## autocorr effect (spatial, ranCoefs...)
      if ( ! is.null(fixZAlist)) {
        if ( ! is.null(isEachNewLevelInOld)) { ## non spatial effect: a vector of booleans indicating whether new is in old (qualifies cols of Cnewold)
          Evar <- Diagonal(x=lambda[old_rd]*as.numeric(! isEachNewLevelInOld))
        } else { ## spatial effect
          Cno_InvCoo_Cof <- (cov_newLv_oldv_list[[new_rd]])[] %id*id% (invCov_oldLv_oldLv_list[[old_rd]])[] %id*id% t(cov_fixLv_oldv_list[[new_rd]])[]
          Evar <- lambda[old_rd] * (cov_newLv_fixLv_list[[new_rd]] - Cno_InvCoo_Cof)
        } 
        terme <- newZAlist[[new_rd]][] %id*id% Evar[] %id*id% t(fixZAlist[[new_rd]])[]
      } else {
        if (is.null(cov_newLv_newLv_list[[new_rd]])) { ## we compute only the var vector
          Cno_InvCoo_Con <- rowSums( ( (cov_newLv_oldv_list[[new_rd]])[] %id*id% (invCov_oldLv_oldLv_list[[old_rd]])[]) * cov_newLv_oldv_list[[new_rd]] ) 
          Evar <- lambda[old_rd] * (1 - Cno_InvCoo_Con)
        } else { ## full cov matrix
          Cno_InvCoo_Con <- (cov_newLv_oldv_list[[new_rd]])[] %id*id% (invCov_oldLv_oldLv_list[[old_rd]])[] %id*id% t(cov_newLv_oldv_list[[new_rd]])[]
          Evar <- lambda[old_rd] * (cov_newLv_newLv_list[[new_rd]] - Cno_InvCoo_Con)
        }
        terme <- .calcZWZt_mat_or_diag(newZAlist[[new_rd]],Evar,covMatrix)
      }
    }
    if (covMatrix) {
      Evarlist[[new_rd]] <- as.matrix(terme)
    } else Evarlist[[new_rd]] <- terme
  }
  Evar <- Reduce("+",Evarlist)
  return(Evar)
}

.submatrices_for_disp_effect <- function(logdispObject, col_info, newinold) { ## logdisp_cov is complete, newinold informs about retained ranefs
  logdisp_cov <- logdispObject$logdisp_cov
  dwdlogdisp <- logdispObject$dwdlogdisp
  cum_n_u_h <- col_info$cum_n_u_h
  #
  which_ranefs <- intersect(col_info$ranef_ids,newinold)
  u_range <- vector("list",length = length(which_ranefs))
  for (it in which_ranefs) u_range[[it]] <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
  u_range <- unlist(u_range)
  #
  which_ranef_cols <- which(col_info$ranef_ids %in% newinold) ## for ranCoefs, several elements of col_info$ranef_ids can match one in newinold
  whichcols <- c(which_ranef_cols,col_info$phi_cols)
  dwdlogdisp <- dwdlogdisp[u_range,whichcols]
  logdisp_cov <- logdisp_cov[whichcols,whichcols]
  return(list(dwdlogdisp=dwdlogdisp,logdisp_cov=logdisp_cov))
}

.calc_sliceVar <- function(slice, cov_newLv_oldv_list,
                           X.pv, newZAlist_slice, 
                           beta_w_cov,
                           cov_newLv_newLv_list,
                           invCov_oldLv_oldLv_list, ## may be null is no newdata.
                           lambda, ## may be null if no newdata.
                           logdispObject, ## should remain NULL is disp not requested
                           newinold, 
                           covMatrix,blockSize) {
  ## here  the problem is that newZA should map the new response levels 
  ## to the 'new' levels of random effects  
  newnrand <- length(newZAlist_slice) ## or of any other of the lists of matrices
  templateList <- vector("list", length = newnrand)
  ## subset new levels in all relevant matrices. The oldv levels must be untouched !
  locnewZA <- templateList
  locCov_newLv_oldv_list <- templateList
  if ( ! is.null(cov_newLv_newLv_list)) {locCov_newLv_newLv_list <- templateList} else locCov_newLv_newLv_list <- NULL
  for (new_rd in seq_len(newnrand)) {
    req_levels <- which(colSums(newZAlist_slice[[new_rd]])) 
    locnewZA[[new_rd]] <- newZAlist_slice[[new_rd]][ , req_levels, drop=FALSE] 
    locCov_newLv_oldv_list[[new_rd]] <- structure(
      cov_newLv_oldv_list[[new_rd]][req_levels, , drop=FALSE],
      isEachNewLevelInOld = attr(cov_newLv_oldv_list[[new_rd]],"isEachNewLevelInOld")[req_levels] ## for non-spatial effects; (qualifies sub cols of sub Cnewold)
    )
    if ( ! is.null(locCov_newLv_newLv_list)) locCov_newLv_newLv_list[[new_rd]] <- cov_newLv_newLv_list[[new_rd]][req_levels, req_levels, drop=FALSE] 
  }
  .calcPredVar(cov_newLv_oldv_list=locCov_newLv_oldv_list, ## needed only if newdata in predict() call
               X.pv=X.pv,## problem is that this creates the appearence of new data end more calculations
               newZAlist=locnewZA, ## either newZA or newZAC needed even if no newdata in predict() call
               beta_w_cov=beta_w_cov,cov_newLv_newLv_list=locCov_newLv_newLv_list,
               invCov_oldLv_oldLv_list=invCov_oldLv_oldLv_list, ## needed only if newdata in predict() call
               lambda=lambda,
               logdispObject=logdispObject,
               newinold=newinold,
               covMatrix=covMatrix,blockSize=blockSize)
}

.calcPredVar <- function(cov_newLv_oldv_list,
                         X.pv,newZAlist,beta_w_cov,
                         #newZACpplist=NULL,
                         cov_newLv_newLv_list=NULL,
                         invCov_oldLv_oldLv_list=NULL, ## may be null is no newdata.
                        lambda=NULL, ## may be null if no newdata.
                        logdispObject=NULL, ## should remain NULL is disp not requested
                        newinold, 
                        covMatrix=FALSE,blockSize) {
  # predVar (in observed points or elsewhere) uses C rather than L hence we need to compute ZA.C in all cases ('newZAC")
  # and then we need newZA and cov_newLv_oldv_list.
  nrX <-  nrow(X.pv)
  if ( ( ! covMatrix ) && nrX > blockSize) {
    ### this part of code is tested by the test-predVar code on Loaloa data
    # et par test geostat dans probitgem (iterateSEMSmooth -> .sampleNextPars -> .spaMM_rhullByEI)
    slices <- unique(c(seq(0L,nrX,blockSize),nrX))
    nslices <- length(slices)-1L
    predVar <- vector("list",nslices)
    newZAlist_slice <- vector("list",length(newZAlist))
    if (showpbar <- interactive()) pb <- txtProgressBar(style = 3)
    for (it in seq_len(nslices)) {
      slice <- (slices[it]+1L):slices[it+1L]
      X.pv_slice <- X.pv[slice,,drop=FALSE]
      for (rt in seq_along(newZAlist)) newZAlist_slice[[rt]] <- newZAlist[[rt]][slice,,drop=FALSE]
      predVar[[it]] <- .calc_sliceVar( ## -> recursive call to .calcPredVar 
                                      slice=slice,
                                      cov_newLv_oldv_list=cov_newLv_oldv_list,
                                      X.pv=X.pv_slice,
                                      newZAlist_slice=newZAlist_slice,
                                      beta_w_cov=beta_w_cov,
                                      cov_newLv_newLv_list=cov_newLv_newLv_list,
                                      invCov_oldLv_oldLv_list=invCov_oldLv_oldLv_list, ## may be null is no newdata.
                                      lambda=lambda, ## may be null if no newdata.
                                      logdispObject=logdispObject, ## should remain NULL is disp not requested
                                      newinold=newinold, 
                                      covMatrix=covMatrix,blockSize=blockSize)
      if (showpbar) setTxtProgressBar(pb, slices[it+1L]/nrX) ## update progress bar
    }
    if (showpbar) close(pb)
    return(unlist(predVar))
  }
  ############################################################
  newZACvar <- .calc_newZACvar(newZAlist,cov_newLv_oldv_list)
  newAugX <- cbind2(X.pv,newZACvar) ## mais en fait pas un AugX since it uses C (in C.w) rather than L (in L.v)
  ## First component of predVar
  # variance of expectation of Xbeta+Zb due to var of (hat(beta),hat(v)) using E[b] as function of hat(v)
  predVar <- .calcZWZt_mat_or_diag(newAugX,beta_w_cov,covMatrix) ## component for linPred=TRUE,disp=FALSE when newdata=ori data
  #### next line of code is an imperfect but useful patch
  # possible problems:
  # * (unhandled) Evar can also be (numerically) negative
  # * (quick-patched) same for min_eigen 
  if ( is.null(dim(predVar))) {
    predVar <- pmax(apply(newAugX^2,1L,sum)*max(0,attr(beta_w_cov,"min_eigen")), predVar)
  } else predVar <- pmax(.tcrossprod(newAugX)*max(0,attr(beta_w_cov,"min_eigen")), predVar)
  ## Second component of predVar:
  # Evar: expect over distrib of (hat(beta),hat(v)) of [variance of Xbeta+Zb given (hat(beta),hat(v))]
  if ( ! is.null(invCov_oldLv_oldLv_list)) { ## must be equivalent to the presence of newdata
    Evar <- .calc_Evar(newZAlist=newZAlist,newinold=newinold, cov_newLv_oldv_list=cov_newLv_oldv_list, 
               lambda=lambda, invCov_oldLv_oldLv_list=invCov_oldLv_oldLv_list, 
               cov_newLv_newLv_list=cov_newLv_newLv_list, covMatrix=covMatrix)
    # Evar must be 0 when newdata=ori data
    predVar <- predVar + Evar ## component for linPred=TRUE,disp=FALSE whether newdata=ori data or not
  }
  # If components for uncertainty in dispersion params were requested,
  #   logdispObject is not NULL
  # If some components ere computable, $$dwdlogdisp should not be NULL
  # Former approach (changed 08/2016) was to test logdispObject and then 
  #   for any 'problem'. But there may some 'problem' and still a valid logdispObject
  # => In new version, dwdlogdisp should be either NULL or a conforming matrix;
  #  'problems" should not be tested.
  if ( ! is.null(logdispObject$dwdlogdisp) ) {
    col_info <- attr(logdispObject$dwdlogdisp,"col_info") ## ranefs for which there is a col in dwdlogdisp
    if ( length(newinold) != col_info$nrand) { ## selection of blocks for re.form ranefs # fixme: this code is not (testthat-)checked
      submatrices <- .submatrices_for_disp_effect(logdispObject, col_info, newinold) 
    } else submatrices <- logdispObject
    newZACw <- newZACvar %*% submatrices$dwdlogdisp ## typically (nnew * n_u_h) %*% (n_u_h * 2) = nnew * 2 hence small 
    if (covMatrix) {
      disp_effect_on_newZACw <- newZACw %*% submatrices$logdisp_cov %*% t(newZACw)  
    } else {
      premul <- newZACw %*% submatrices$logdisp_cov
      disp_effect_on_newZACw <- rowSums(premul * newZACw)
    }
    predVar <- predVar + disp_effect_on_newZACw
  }
  return(predVar) ## may be a Matrix
}


.calcResidVar <- function(object,newdata=NULL) {
  phi.object <- object$phi.object
  if (is.null(phi_outer <- phi.object$phi_outer)) { ## valid whether newdata are NULL or not:
    glm_phi <- phi.object[["glm_phi"]]
    if (is.null(glm_phi)) glm_phi <- .get_glm_phi(object)
    residVar <- predict(glm_phi, newdata=newdata, type="response")
  } else { ## phi, but not glm_phi
    if (length(phi_outer)==1L) {
      if (is.null(newdata)) {
        residVar <- rep(phi_outer,nrow(object$X.pv)) ## assumes (length(phi_outer)==1L)           
      } else residVar <- rep(phi_outer,nrow(newdata))
    } else stop("Unable to compute 'residVar' given length(phi_outer)!=1L.") ## and no glm_phi
    ## FR->FR we could improve this if we get a glm_phi when phi was estimed by outer iterations
  }
  residVar
} 

.make_new_corr_lists <- function(object,locdata, which_mats, newZAlist, newinold, fix_info=NULL) {
  if (which_mats$no) {
    cov_newLv_oldv_list <- vector("list",length(newZAlist)) ## declar-initalization, will be filled in the loop
  } else cov_newLv_oldv_list <- .make_corr_list(object$strucList,newZAlist=NULL) # we always need a non trivial value
  cov_newLv_newLv_list <- vector("list",length(cov_newLv_oldv_list))
  ranefs <- attr(newZAlist,"exp_ranef_strings") ## only the matrices that get a "ranefs" attribute will be used in .compute_ZAXlist()
  if (any(unlist(which_mats))) {
    strucList <- object$strucList
    spatial_terms <- attr(object$ZAlist,"exp_spatial_terms")
    for (new_rd in seq_along(cov_newLv_oldv_list)) {
      old_rd <- newinold[new_rd]
      if ( which_mats$no || which_mats$nn[new_rd]) {
        corr.model <- attr(strucList[[old_rd]],"corr.model")
        if (is.null(corr.model)) {
          ## the list elements remain NULL
        } else if ( corr.model=="") {
          stop("corr.model==''") ## FIXME remove this test ?
        } else if ( corr.model=="random-coef") {
          newZA <- newZAlist[[new_rd]]  
          namesTerms <- attr(newZAlist,"namesTerms")[[new_rd]]
          if ( ! is.null(fix_info)) {
            oldlevels <- colnames(fix_info$newZAlist[[old_rd]]) ## old_rd despite the "newZAlist" name: specific to fix_info
          } else oldlevels <- colnames(object$ZAlist[[old_rd]])
          newlevels <- colnames(newZA)
          oldornew <- unique(c(oldlevels,newlevels))
          numnamesTerms <- length(namesTerms) ## 2 for random-coef
          design_u <- attr(strucList[[old_rd]],"latentL_blob")$u
          newoldC <- .makelong(.tcrossprod(design_u), longsize=length(oldornew)*numnamesTerms)
          repnames <- rep(oldornew, numnamesTerms)
          newcols <- repnames %in% newlevels ## handle replicates (don't `[` newoldC using names !)
          oldcols <- repnames %in% oldlevels
          colnames(newoldC) <- rownames(newoldC) <- repnames ## these will be needed by .match_old_new_levels()
          if (which_mats$no) cov_newLv_oldv_list[[new_rd]] <- structure(newoldC[newcols,oldcols,drop=FALSE],ranefs=ranefs[[new_rd]])
          if (which_mats$nn[new_rd]) cov_newLv_newLv_list[[new_rd]] <- newoldC[newcols,newcols,drop=FALSE]
        } else {
          if ( ! is.null(fix_info)) {
            olduniqueGeo <- fix_info$newuniqueGeo
          } else olduniqueGeo <- attr(object,"info.uniqueGeo")
          geonames <- colnames(olduniqueGeo)
          newuniqueGeo <- locdata[,geonames,drop=FALSE] ## includes nesting factor
          ## distance matrix and then call to correl fn:
          if (corr.model=="AR1") {
            ### older, non nested code:
            # if (which$no) resu$uuCnewold <- proxy::dist(newuniqueGeo,olduniqueGeo) 
            # if (which$nn) resu$uuCnewnew <- proxy::dist(newuniqueGeo)
            ### new code recycling .get_dist_nested_or_not, but not fastest nor most transparent (fixme?)
            onGeo <- rbind(newuniqueGeo,olduniqueGeo) # includes nesting factor
            if (object$spaMM.version < "2.2.118") {
              blob <- .get_dist_nested_or_not(object$spatial_term, data=onGeo, distMatrix=NULL, uniqueGeo=NULL, 
                                              dist.method=object$control.dist$dist.method, as_matrix=TRUE,
                                              needed=c(distMatrix=TRUE))
            } else {
              blob <- .get_dist_nested_or_not(spatial_terms[[old_rd]], data=onGeo, distMatrix=NULL, uniqueGeo=NULL, 
                                            dist.method=object$control.dist$dist.method, as_matrix=TRUE,
                                            needed=c(distMatrix=TRUE)) ## FIXME provide uniqueGeo to save time ?
            }
            ## we merged old and new so need to get the respective cols (which may overlap) 
            uli_onGeo <- .ULI(onGeo) # this should give row and colums in the blob ## FIXME how to make sure of that? .get_dist_nested_or_not must use .ULI()
            uli_new <- uli_onGeo[seq(nrow(newuniqueGeo))]
            uli_old <- uli_onGeo[-seq(nrow(newuniqueGeo))]
            if (which_mats$no) uuCnewold <- blob$distMatrix[uli_new,uli_old,drop=FALSE] ## rows match the newZAlist, cols match th u_h 
            if (which_mats$nn[new_rd]) uuCnewnew <- blob$distMatrix[uli_new,uli_new,drop=FALSE]
          } else {
            ### rho only used to compute scaled distances
            rho <- .getPar(object$ranFix,"rho")
            if ( ! is.null(rho_mapping <- attr(object,"dist_info")$rho.mapping)
                 && length(rho)>1L ) rho <- .calc_fullrho(rho=rho,coordinates=geonames,rho_mapping=rho_mapping)
            ## rows from newuniqueGeo, cols from olduniqueGeo:
            msd.arglist <- list(uniqueGeo=newuniqueGeo,uniqueGeo2=olduniqueGeo,
                                rho=rho,return_matrix=TRUE)
            if ( ! is.null(dist.method <- object$control.dist$dist.method)) msd.arglist$dist.method <- dist.method
            if (which_mats$no) uuCnewold <- do.call(make_scaled_dist,msd.arglist) ## ultimately allows products with Matrix ## '*cross*dist' has few methods, not even as.matrix
            if (which_mats$nn[new_rd])  {
              msd.arglist$uniqueGeo2 <- NULL
              if (nrow(msd.arglist$uniqueGeo)==1L) {
                uuCnewnew <- matrix(0)
              } else uuCnewnew <- do.call(make_scaled_dist,msd.arglist) 
            }
          }
          if (which_mats$no) cov_newLv_oldv_list[[new_rd]] <- structure(.calc_corr_from_dist(uuCnewold, object, corr.model),
                                                                    ranefs=ranefs[[new_rd]])
          if (which_mats$nn[new_rd]) cov_newLv_newLv_list[[new_rd]] <- .calc_corr_from_dist(uuCnewnew, object, corr.model)
        }
      }
    }
  }
  return(list(cov_newLv_oldv_list=cov_newLv_oldv_list, cov_newLv_newLv_list=cov_newLv_newLv_list))
}

.process_variances <- function(variances) {
  if (identical(variances$BH98,TRUE)) { ## my interpretation of BH98
    variances$predVar <- TRUE ## implies $linPred <- TRUE by default
    #variances$respVar <- FALSE ## default
    variances$disp <- FALSE ## non-default
    variances$cov <- TRUE ## non-default
  } else variances$BH98 <- FALSE
  # $respVar implies components
  if (is.null(variances$predVar)) variances$predVar <- variances$respVar ## may still be NULL
  if (is.null(variances$residVar)) variances$residVar <- variances$respVar ## may still be NULL
  if (is.null(variances$respVar)) variances$respVar <- FALSE 
  # $predVar implies components
  if (is.null(variances$linPred)) variances$linPred <- variances$predVar ## may still be NULL
  if (is.null(variances$disp)) variances$disp <- variances$predVar ## may still be NULL
  if (is.null(variances$predVar)) variances$predVar <- FALSE 
  # Do not let any component empty
  if (is.null(variances$fixefVar)) variances$fixefVar <- FALSE 
  if (is.null(variances$linPred)) variances$linPred <- FALSE 
  if (is.null(variances$disp)) variances$disp <- FALSE ## uncertaintly on dispersion parameters
  if (is.null(variances$residVar)) variances$residVar <- FALSE ## uncertaintly on dispersion parameters
  if (is.null(variances$cov)) variances$cov <- FALSE
  ##
  return(variances)
}

.match_old_new_levels <- function(new_rd, old_cum_n_u_h, newinold, spatial_old_rd, w_h_coeffs, subZAlist, newZACpplist, lcrandfamfam, object) {
  old_rd <- newinold[new_rd]
  oldu.range <- (old_cum_n_u_h[old_rd]+1L):(old_cum_n_u_h[old_rd+1L])
  if (old_rd %in% spatial_old_rd) { 
    return(w_h_coeffs[oldu.range])          
  } else {
    oldlevels <- colnames(subZAlist[[new_rd]])
    newlevels <- colnames(newZACpplist[[new_rd]])
    interlevels <- intersect(oldlevels,newlevels)
    oldpos <- which(oldlevels %in% interlevels) ## positions: handle replicates for random-coef
    newpos <- which(newlevels %in% interlevels)
    oldv <- w_h_coeffs[oldu.range]
    names(oldv) <- oldlevels
    psi_M <- switch(lcrandfamfam[new_rd], 
                    gaussian = 0,
                    gamma = 1, 
                    beta = 1/2, 
                    "inverse.gamma" = 1
    )
    vpsi_M <- object$rand.families[[old_rd]]$linkfun(psi_M) 
    ## since vpsi_M can be non-zero, the expectation of the response can be modified in a re.form model compared to the original
    newv <- rep(vpsi_M,length(newlevels)) ## fills new levels with psi_M
    names(newv) <- newlevels
    newv[newpos] <- oldv[oldpos] 
    return(newv)
  }
}

.point_predict <- function(object, newdata, variances, eta_fix, re.form,
                           newinold, spatial_old_rd, subZAlist, newZACpplist ) {
  if (.noRanef(re.form)) {
    if ( ! is.null(zero_truncated <- object$family$zero_truncated)) {
      fv <- object$family$linkinv(eta_fix, mu_truncated=zero_truncated)
    } else fv <- object$family$linkinv(eta_fix) 
    return(list(fv=fv,eta=eta_fix))
  } else if ( is.null(newdata) && ! inherits(re.form,"formula")) {
    return(list(fv=object$fv)) ## eta will be reconstructed from fv on request
  } else { ## 
    newnrand <- length(subZAlist)
    if ( newnrand==0L ) {
      ZACw <- 0 
    } else {
      ## for random-slope model the original eta's can be recomputed as 
      # object$X.pv %*% fixef(object) + 
      #    object$ZAlist[[1]] %*% object$strucList[[1]] %*% object$v_h
      #
      #### precomputation of coeffs
      ## on the gaussian scale, L.v_ori ~ lam C (lam C + phi I)^{-1}y 
      ## new random autocorr term ~ lam c (lam C + phi I)^{-1}y = c C^{-1} L_ori.v_ori = c [t(L_ori)]^{-1} v_ori
      ## [t(L_ori)]^{-1} v_ori can be computed once for all predictions => 'w_h_coeffs'
      if (is.null(w_h_coeffs <- object$envir$w_h_coeffs)) { 
        object$envir$w_h_coeffs <- w_h_coeffs <- .calc_invL_coeffs(object,object$v_h)
      }
      lcrandfamfam <- attr(object$`rand.families`,"lcrandfamfam")[newinold]
      
      augm_w_h_coeffs <- lapply(seq_len(newnrand),.match_old_new_levels,
                                old_cum_n_u_h=attr(object$lambda,"cum_n_u_h"), newinold=newinold, 
                                spatial_old_rd=spatial_old_rd, w_h_coeffs=w_h_coeffs, 
                                subZAlist=subZAlist, newZACpplist=newZACpplist, lcrandfamfam=lcrandfamfam, object=object)
      if (newnrand>1L) {
        ZACw <- vector("list", length=newnrand )
        for (it in seq_len(newnrand)) ZACw[[it]] <- drop(newZACpplist[[it]][] %*% augm_w_h_coeffs[[it]])
        ZACw <- do.call(cbind,ZACw)
        ZACw <- rowSums(ZACw)
      } else ZACw <- drop(newZACpplist[[1]][] %*% augm_w_h_coeffs[[1]])
    }
    eta <- eta_fix + ZACw ## (length(eta)) col vector from coeffs = length(eta) row vector...
    # done with eta
    if (variances$BH98) {
      fv <- eta
    } else if ( ! is.null(zero_truncated <- object$family$zero_truncated)) {
      fv <- object$family$linkinv(eta,mu_truncated=zero_truncated)
    } else fv <- object$family$linkinv(eta) ## ! freqs for binomial, counts for poisson
    return(list(fv=fv,eta=eta))
  }
} 

.predict_body <- function(object, newdata, re.form,
                          variances, binding, intervals, level, blockSize) {
  new_X_ZACblob <- .calc_new_X_ZAC(object=object, newdata=newdata, re.form = re.form,
                                   variances=variances)
  newinold <- new_X_ZACblob$newinold ## says which ranef is kept by re.form
  #
  ## (1) computes fv (2) compute predVar
  ##### fv
  ppblob <- .point_predict(object, newdata, variances, 
                       eta_fix=new_X_ZACblob$etaFix, ## may be NULL. addition of random-effect terms in the function
                       re.form, newinold=newinold, spatial_old_rd= new_X_ZACblob$spatial_old_rd,
                       subZAlist=new_X_ZACblob$subZAlist, ## (subset of) the original object$ZAlist 
                       newZACpplist=new_X_ZACblob$newZACpplist )
  resu <- ppblob$fv
  if ( ! is.na(binding)) resu <- as.matrix(resu) ## suitable for objective function of optim() etc ## matrix ! maybe more suitable than data frame as objective function
  # if (identical(object$family$zero_truncated,TRUE)) {
  #   attr(resu,"p0") <- attr(ppblob$fv,"p0")
  #   attr(resu,"mu_U") <- attr(ppblob$fv,"mu_U")
  # }
  #\item 'make.names(.,unique= TRUE)' is applied to the input data, so the row names may be modified, if the input data did not contain unique, syntactically valid row names as defined by 'help(make.names)'.
  #rownames(resu) <- make.names(rownames(resu),unique = TRUE)
  locdata <- new_X_ZACblob$locdata
  if ( ! is.logical(binding) ) { ## expecting a string
    binding <- .makenewname(base=binding,varnames=colnames(locdata)) ## 09/11/2014 = posterior to CRAN 1.4.1 
    resu <- data.frame(resu)
    colnames(resu) <- binding
    resu <- cbind(locdata,resu) ## becomes a data frame !
    attr(resu,"fittedName") <- binding
  } else { ## expecting binding= FALSE
    if (! is.na(binding) && ncol(locdata)>0)  attr(resu,"frame") <- locdata 
  }
  ##### (2) predVar
  newnrand <- length(new_X_ZACblob$newZAlist) ## may be reduced if non trivial re.form
  newX.pv <- new_X_ZACblob$newX.pv
  if(variances$linPred) {
    beta_cov <- .get_beta_cov_any_version(object)
    beta_w_cov <- .get_beta_w_cov(object)
    if (inherits(re.form,"formula")) {
      # identifies an selects columns for the [retained ranefs, which are given by newinold 
      subrange <- vector("list",length = length(newinold))
      old_cum_n_u_h <- attr(object$lambda,"cum_n_u_h")
      for (it in newinold) subrange[[it]] <- (old_cum_n_u_h[it]+1L):(old_cum_n_u_h[it+1L])
      subrange <- unlist(subrange)
      Xncol <- ncol(beta_cov)
      subrange <- c(seq_len(Xncol),subrange + Xncol)
      beta_w_cov <- beta_w_cov[subrange,subrange]
    }
    if (newnrand) {
      if ( is.null(newdata) && ! inherits(re.form,"formula")) {
        newZAlist <- new_X_ZACblob$subZAlist ## (subset of) the original object$ZAlist 
      } else { ## 
        newZAlist <- new_X_ZACblob$newZAlist ## a new ZAlist built using .spMMFactorList(new formula...)
      }
    }
    if ( ! is.null(newdata)) {
      invCov_oldLv_oldLv_list <- .get_invColdoldList(object)
    } else {
      invCov_oldLv_oldLv_list <- NULL
    }
    if (newnrand) {
      if (variances$BH98) newX.pv[] <- 0
      loclist <- list(X.pv=newX.pv,beta_w_cov=beta_w_cov,covMatrix=variances$cov,lambda=object$lambda,
                      newinold=new_X_ZACblob$newinold,blockSize=blockSize, 
                      newZAlist=newZAlist,
                      cov_newLv_oldv_list=new_X_ZACblob$cov_newLv_oldv_list,
                      ## list for Cnewnew, which enters in  newZA %*% Cnewnew %*% tnewZA, hence should not represent newZA itself 
                      cov_newLv_newLv_list=new_X_ZACblob$cov_newLv_newLv_list)
      if ( ! is.null(newdata)) loclist$invCov_oldLv_oldLv_list <- invCov_oldLv_oldLv_list
      if (variances$disp) loclist$logdispObject <- .get_logdispObject(object)
      respVar <- do.call(.calcPredVar,loclist) 
      if (variances$cov) {
        respVar <- as.matrix(respVar) ## matrix, not Matrix (assumed below)
        rownames(respVar) <- colnames(respVar) <- rownames(locdata)
      } else {
        names(respVar) <- rownames(locdata)
      }
    } else {
      if (variances$cov) {
        respVar <- matrix(0,nrow=nrow(locdata),ncol=nrow(locdata))
      } else respVar <- rep(0,nrow(locdata))
    }
  } else if (any(unlist(variances))) {
    respVar <- rep(0,nrow(locdata))
  } else respVar <- NULL 
  beta_cov <- .get_beta_cov_any_version(object)
  if (! is.null(beta_cov)) {
    if ( variances$fixefVar || (newnrand==0L && variances$linPred) ) {
      fixefcov <- newX.pv %*% beta_cov %*% t(newX.pv)
      if (variances$cov) {
        attr(resu,"fixefVar") <- fixefcov 
      } else attr(resu,"fixefVar") <- diag(fixefcov)
      if (newnrand==0L) { ## otherwise there is already such a term in predVar
        respVar <- respVar + attr(resu,"fixefVar") 
      }
    }
  }
  attr(resu,"predVar") <- respVar ## vector or matrix
  if ( ! is.null(respVar) && object$family$link!="identity") {
    if (is.null(eta <- ppblob$eta)) eta <- object$family$linkfun(ppblob$fv) 
    dmudeta <- object$family$mu.eta(eta)
    if (!is.null(dim(respVar))) {
      respVar <- sweep(respVar, MARGIN = 1, dmudeta, `*`) ## premul
      respVar <- sweep(respVar, MARGIN = 2, dmudeta, `*`) ## postmul
    } else respVar <- respVar*(dmudeta^2)
  }
  if (variances$residVar) {
    pw <- object$prior.weights
    if ( ! (attr(pw,"unique") && pw[1L]==1L)) {
      if (! identical(spaMM.getOption("prior_weights_residvar_warned"),TRUE)) {
        warning("Prior weights are not taken in account in residVar computation.")
        spaMM.options(prior_weights_residvar_warned=TRUE)
      }
    }
    if (object$family$family %in% c("poisson","binomial","COMPoisson","negbin")) {
      attr(resu,"residVar") <- object$family$variance(ppblob$fv)
    } else attr(resu,"residVar") <- .calcResidVar(object,newdata=locdata) 
    if (inherits(respVar,"matrix")) {
      nc <- ncol(respVar)
      diagPos <- seq.int(1L,nc^2,nc+1L)
      respVar[diagPos] <- respVar[diagPos] + attr(resu,"residVar")
    } else respVar <- respVar + attr(resu,"residVar")
  }
  if (variances$respVar) attr(resu,"respVar") <- respVar
  if ( is.matrix(resu) && NCOL(resu)==1L) {
    class(resu) <- c("predictions",class(resu))
  } ## for print.predictions method which expects a 1-col matrix
  # intervals
  checkVar <- setdiff(intervals,names(attributes(resu)))
  if (length(checkVar)) {
    warning(paste("Variances",paste(checkVar,collapse=", "),
                  "not available for interval computation.\n Check arguments."))
    intervals <- intersect(intervals,names(attributes(resu)))
  } 
  if(length(intervals)) {
    intervalresu <- NULL
    for (st in intervals) {
      varcomp <- attr(resu,st)
      if (is.null(varcomp)) warning(paste("Prediction variance component",st,"requested but not available: check input."))
      if (is.matrix(varcomp)) varcomp <- diag(varcomp)
      eta <- object$family$linkfun(resu[,1L])
      pv <- 1-(1-level)/2
      ## special case for simple LM
      if (length(object$rand.families)==0L && # not mixed
          object$family$family=="gaussian" &&
          deparse(object$resid.predictor)=="~1" # not heteroscedastic
      ) { 
        nobs <- length(object$y)
        resdf <- nobs - ncol(object$X.pv) ## don't use fixef here, that contains bot NAs and argument etaFix$beta! 
        is_REML <- ( .REMLmess(object,return_message=FALSE))
        if ( ! is_REML) {
          vart <- varcomp*nobs/resdf
        } else vart <- varcomp
        ## FR->FR (use a type attribute for fixef ?)
        sd <- stats::qt(pv,df=resdf)*sqrt(vart)
      } else {
        sd <- qnorm(pv)*sqrt(varcomp)
      }
      if (! is.null(zero_truncated <- object$family$zero_truncated)) {
        interval <- cbind(object$family$predict(eta-sd, mu_truncated=zero_truncated),
                          object$family$predict(eta+sd, mu_truncated=zero_truncated))
      } else interval <- cbind(object$family$linkinv(eta-sd),object$family$linkinv(eta+sd))
      colnames(interval) <- paste(st,c(signif(1-pv,4),signif(pv,4)),sep="_")
      intervalresu <- cbind(intervalresu,interval)
    }
    attr(resu,"intervals") <- intervalresu
  }
  return(resu)
}

.unlist_with_attributes <- function(somelist) {
  if (is.matrix(somelist[[1L]])) {
    res <- do.call(rbind,somelist)
  } else res <- unlist(somelist)
  mostAttrs <- names(attributes(somelist[[1L]]))
  mostAttrs <- setdiff(mostAttrs,c("class","dim","dimnames","names","row.names"))
  tmpattr <- vector("list",length(somelist))
  for (attrname in mostAttrs) {
    for (it in seq_along(somelist)) tmpattr[[it]] <- attr(somelist[[it]],attrname)
    if (is.matrix(tmpattr[[1L]])) {
      attr(res,attrname) <- do.call(rbind,tmpattr)
    } else attr(res,attrname) <- unlist(tmpattr)
  }
  return(res)
}

## (1) for surface prediction: (developed in InferentialSimulation/InferentialSimulation.R)
## (2) But also for generation of fixed effects in simulation of nested-effects models
predict.HLfit <- function(object, newdata = newX, newX = NULL, re.form = NULL,
                          variances=list(),
                          binding = FALSE, 
                          intervals = NULL,
                          level = 0.95,
                          blockSize = 2000L,
                          ...) { ## but not new Y
  if ( ! is.null(variances$ranef)) {
    message("'variances$ranef' is obsolete: I interpret this as 'variances$linPred'.")
    variances$linPred <- variances$ranef
    variances$ranef <- NULL
  }
  if ( ! is.null(variances$sum)) {
    message("'variances$sum' is obsolete: I interpret this as 'variances$respVar'.")
    variances$respVar <- variances$sum
    variances$sum <- NULL
  }
  if (is.null(object$envir)) object$envir <- list2env(list(), ## back-compatibility fix for old objects
                                                     parent=environment(HLfit_body))
  ## the final components returned as attributes have names ...Var, other terms should be named differently
  checkIntervals <- (substr(x=intervals, nchar(intervals)-2, nchar(intervals))=="Var")
  if (any(!checkIntervals)) warning("Element(s)",intervals[!checkIntervals],"are suspect, not ending in 'Var'.")
  # possible elements in return value: fixefVar, predVar, residVar, respVar
  variances[intervals] <- TRUE 
  variances <- .process_variances(variances)
  nrX <-  NROW(newdata)
  ############################## if (nrX>0L) newdata <- droplevels(newdata) FIXME perhaps here ? 
  if ( (! variances$cov) && nrX > blockSize) {
    ### this part of code is tested by the test-predVar code on Loaloa data
    # et par test geostat dans probitgem (iterateSEMSmooth -> .sampleNextPars -> .spaMM_rhullByEI)
    ## newdata <- droplevels(newdata) ## potential gain of time for droplevels(newdata_slice)
    slices <- unique(c(seq(0L,nrX,blockSize),nrX))
    nslices <- length(slices)-1L
    res <- vector("list",nslices)
    if (showpbar <- interactive()) pb <- txtProgressBar(style = 3)
    for (it in seq_len(nslices)) {
      slice <- (slices[it]+1L):slices[it+1L]
      newdata_slice <- newdata[slice,,drop=FALSE]
      ## newdata_slice <- droplevels(newdata_slice) 
      res[[it]] <- .predict_body(object=object, newdata=newdata_slice, re.form = re.form, variances=variances, 
                                 binding=binding, intervals=intervals, level=level, blockSize=blockSize) ## blocksize shuld not be useful *here*
      if (showpbar) setTxtProgressBar(pb, slices[it+1L]/nrX) ## update progress bar
    }
    if (showpbar) close(pb)
    res <- .unlist_with_attributes(res)
    return(res)
  } else .predict_body(object=object, newdata=newdata, re.form = re.form,
                variances=variances, binding=binding,
                intervals=intervals, level=level, blockSize=blockSize)  ## but blockSize could be useful *here* if newdata was NULL
}

print.vcov.HLfit <- function(x, expanded=FALSE, ...) {
  a <- attributes(x)
  attr(x,"beta_v_cov") <- NULL  
  print.default(x)
  cat("with additional attribute(s):")
  std.attr <- c("names","dim","dimnames","class") ## attributes not to be shown
  nam <- names(a)
  if (expanded) { # shows structure of attributes as in utils:::str.default
    cat("\n")
    nest.lev <- 0
    indent.str <- paste(rep.int(" ", max(0, nest.lev + 1)), collapse = "..")
    strO <- getOption("str")
    strO <- modifyList(strOptions(), strO) ## seems to provide a format.fun function
    `%w/o%` <- function(x, y) x[is.na(match(x, y))]
    nfS <- names(fStr <- formals())
    ## this scans the substructure of each attribute
    strSub <- function(obj, ...) {
      nf <- nfS %w/o% c("object", "give.length", "comp.str", 
                        "no.list", names(match.call())[-(1:2)], "...")
      aList <- as.list(fStr)[nf]
      aList[] <- lapply(nf, function(n) eval(as.name(n)))
      strObj <- function(...) str(obj, ...)
      do.call(strObj, c(aList, list(...)), quote = TRUE)
    }
    for (i in seq_along(a)) if (all(nam[i] != std.attr)) {
      cat(indent.str, paste0("- attr(*, \"", nam[i], "\")="), 
          sep = "")
      strSub(a[[i]], give.length = TRUE, indent.str = paste(indent.str, 
                                                            ".."), nest.lev = nest.lev + 1)
    }
  } else {
    cat(" ")
    nam <- setdiff(nam,std.attr)
    cat(paste(nam,collapse=", "))  
    cat("\n")
  }
  invisible() ## do not return x since it has lost a useful attribute
}


`[.predictions` <- function (x, i, j, 
                             drop = TRUE ## by default, this function will return scalar/vector  
                             ) {
  class(x) <- "matrix" ## removes "predictions" => set back later
  #   if (is.data.frame(x)) {
  #     resu <- x[i,j]
  #   } else 
  resu <- x[i,j,drop=drop]
  if ( ! drop) {
    fixefVar <- attr(x, "fixefVar")
    if ( ! is.null(fixefVar)) {
      if (is.null(dim(fixefVar))) {
        fixefVar <- fixefVar[x]
      } else fixefVar <- fixefVar[x,x,drop=FALSE]
    }
    predVar <- attr(x, "predVar")
    if ( ! is.null(predVar)) {
      if (is.null(dim(predVar))) {
        predVar <- predVar[x]
      } else predVar <- predVar[x,x,drop=FALSE]
    }
    frame <- attr(x, "frame")
    if ( ! is.null(frame)) frame <- frame[x,] ## dataframe => nodrop
    residVar <- attr(x, "residVar")
    if ( ! is.null(frame)) residVar <- residVar[x,drop=FALSE]
    respVar <- attr(x, "respVar")
    if ( ! is.null(respVar)) {
      if (is.null(dim(respVar))) {
        respVar <- respVar[x]
      } else respVar <- respVar[x,x,drop=FALSE]
    }
    class(resu) <- c("predictions","matrix")
    structure(resu,fixefVar=fixefVar,predVar=predVar,residVar=residVar,frame=frame,fittedName=attr(x, "fittedName"))
  } else return(resu)
} # Use unlist() to remove attributes from the return value

print.predictions <- function (x, expanded=FALSE, ...) {
  asvec <- as.vector(x) ## important to remove names and keep them separately
  rnames <- rownames(x)
  toolong <- nchar(rnames)>9
  rnames[toolong] <- paste(substr(rnames[toolong],0,8),".",sep="")
  names(asvec) <- rnames
  cat("Point predictions:\n")
  print(asvec)
  cat("*stored as* 1-col matrix with attributes:")
  std.attr <- c("names","dim","dimnames","class") ## attributes not to be shown
  a <- attributes(x)
  nam <- names(a)
  if (expanded) { # shows structure of attributes as in utils:::str.default
    cat("\n")
    nest.lev <- 0
    indent.str <- paste(rep.int(" ", max(0, nest.lev + 1)), collapse = "..")
    strO <- getOption("str")
    strO <- modifyList(strOptions(), strO) ## seems to provide a format.fun function
    `%w/o%` <- function(x, y) x[is.na(match(x, y))]
    nfS <- names(fStr <- formals())
    ## this scans the substructure of each attribute
    strSub <- function(obj, ...) {
      nf <- nfS %w/o% c("object", "give.length", "comp.str", 
                        "no.list", names(match.call())[-(1:2)], "...")
      aList <- as.list(fStr)[nf]
      aList[] <- lapply(nf, function(n) eval(as.name(n)))
      strObj <- function(...) str(obj, ...)
      do.call(strObj, c(aList, list(...)), quote = TRUE)
    }
    for (i in seq_along(a)) if (all(nam[i] != std.attr)) {
      cat(indent.str, paste0("- attr(*, \"", nam[i], "\")="), 
          sep = "")
      strSub(a[[i]], give.length = TRUE, indent.str = paste(indent.str, 
                                                            ".."), nest.lev = nest.lev + 1)
    }
  } else {
    cat(" ")
    nam <- setdiff(nam,std.attr)
    cat(paste(nam,collapse=", "))  
    cat("\n")
  }
  invisible()
}