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


dimnames.bigq <- function(x) { # colnames() and rownames() will use this for bigq objects
  attr(x, "DIMNAMES")
}

.sparseMatrix_bigq <- function(X) {
  dimNames <- dimnames(X) # using my "DIMNAMES" attribute
  nc <- ncol(X)
  nr <- nrow(X)
  Xv <- X[[]]
  positions <- (Xv!="0")
  if (any(positions)) {
    i <- rep_len(rep.int(seq_len(nc), rep.int(1, nc)), nc*nr)[positions] # cf gl() code
    j <- rep_len(rep.int(seq_len(nc), rep.int(nr, nc)), nc*nr)[positions]
    sparseMatrix(i=i,j=j,x=gmp::asNumeric(Xv[positions]), dims=c(nr,nc), dimnames=dimNames)
  } else sparseMatrix(i=integer(),j=integer(),x=numeric(), dims=c(nr,nc), dimnames=dimNames) # zero matrix in sparse format
}

.mMatrix_bigq <- function(X) {
  dimNames <- dimnames(X) # using my "DIMNAMES" attribute
  nc <- ncol(X)
  nr <- nrow(X)
  Xv <- X[[]]
  positions <- (Xv!="0")
  denseness <- length(which(positions))/length(positions)
  if (denseness>0.3) {
    return(structure(gmp::asNumeric(X), dimnames=dimNames))
  } else if (denseness>0) {
    i <- rep_len(rep.int(seq_len(nc), rep.int(1, nc)), nc*nr)[positions] # cf gl() code
    j <- rep_len(rep.int(seq_len(nc), rep.int(nr, nc)), nc*nr)[positions]
    return(sparseMatrix(i=i,j=j,x=gmp::asNumeric(Xv[positions]), dims=c(nr,nc), dimnames=dimNames))
  } else return(sparseMatrix(i=integer(),j=integer(),x=numeric(), dims=c(nr,nc), dimnames=dimNames)) # zero matrix in sparse format
}

`%ZA*gI%` <- function(A,B) {
  if (inherits(B,"bigq")) B <- .mMatrix_bigq(B)
  return(A %id*id% B)
}

# 'gI' for Id and gmp
`%gI*gI%` <- local({
  gmp_warned <- FALSE
  function(A,B) { 
    A_bigq <- inherits(A,"bigq")
    B_bigq <- inherits(B,"bigq")
    if (A_bigq != B_bigq) stop("code missing here")
    if (A_bigq) return(gmp::`%*%`(A,B))
    if (.is_identity(A)) return(B)
    if (.is_identity(B)) return(A)
    return(A %*% B)
  }
})


.calcZWZt_mat_or_diag <- function(Z,
                                  W, # never identity in the actual calls
                                  cholW,
                                  tcrossfac_W=NULL,
                                  returnMat) { ## fixefVar or fixefVar + a bit of ranefVar
  ## bottleneck for large predVar with large n_u_h (W being a matrix in most cases)
  if (returnMat) {
    if (is.null(tcrossfac_W)) {
      if (inherits(cholW,"try-error")) {
        return(Z %id*% W %*id% t(Z)) ## not ZWZT functions bc W not diag a priori if (returnMat)
      } else { ## better ensures SPDefiniteness
        rhs <- .tcrossprod(cholW,Z) # cW Zt 
        return(.crossprod(rhs)) #Z cW' cW Zt # here cholW is a *c*rossprod factor
      } 
    } else { ## .tcrossprod() is memory bottleneck (for lhs nobs* ~n_u_h) => use as_tcrossfac_list to bypass .calcZWZt_mat_or_diag()
      lhs <- Z %*% tcrossfac_W # Z tcW 
      if (inherits(lhs,"dgeMatrix")) { # Z sparse, but decorr was used => tcrossfac_W is matrix (totally OK)
        return(.tcrossprodCpp(as.matrix(lhs),yy=NULL))
      } else return(.tcrossprod(lhs)) # Z tcW tcW' Z' 
    }
  } else { ## returns only the diagonal
    if (is.vector(W)) {
      return(rowSums(.Matrix_times_Dvec(Z,W) * Z))
    } else { # even if returnMat is false, we need the full beta_w_cov if it is not Diagonal.
      if (is.null(tcrossfac_W)) {
        rownames(W) <- colnames(W) <- NULL ## inhibits a useless warning from Matrix:: dimNamesCheck 
        premul <- Z %id*% W # it would be nice if W was provided as a *t*crossprod factor
        return(rowSums(suppressMessages(premul * Z))) ## suppress message("method with signature...") [found by debug(message)]
      } else {
        lhs <- Z %*% tcrossfac_W # Z tcW 
        return(rowSums(lhs^2)) 
      }
    }
  }
}

.calc_Evar <- function(newZAlist, newinold, cov_newLv_oldv_list, invCov_oldLv_oldLv_list, covMatrix,## arguments for all calls
                       # arguments for standard diag or full Evar matrix:
                       cov_newLv_newLv_list, 
                       ## arguments for non-diagonal block of Evar matrix:
                       cov_newLv_fixLv_list, cov_fixLv_oldv_list, fixZAlist=NULL ,
                       diag_cov_newLv_newLv_list,
                       object,
                       as_tcrossfac_list=FALSE
) {
  if (as_tcrossfac_list && ! is.null(fixZAlist)) {stop("Incompatible arguments: as_tcrossfac_list=TRUE and use of fix_X_ZAC object")}
  newnrand <- length(newinold)
  Evarlist <- vector("list",newnrand)
  lambda_list <- object$lambda.object$lambda_list
  for (new_rd in seq_along(Evarlist)) {
    old_rd <- newinold[new_rd]
    exp_ranef_types <- attr(newZAlist,"exp_ranef_types")
    isEachNewLevelInOld <- attr(cov_newLv_oldv_list[[new_rd]],"isEachNewLevelInOld") ## non autocorr effect: a vector of booleans indicating whether new is in old (qualifies cols of Cnewold)
    if (exp_ranef_types[new_rd] %in% c("IMRF","adjacency","corrMatrix")) { ## in that case the newLv = the oldLv Cno=Coo... and actually the Evar term is zero 
      if ( ! is.null(fixZAlist)) {
        terme <- Matrix(0,nrow=nrow(newZAlist[[new_rd]]),ncol=nrow(fixZAlist[[new_rd]]))
      } else {
        terme <- rep(0, nrow(newZAlist[[new_rd]]))
        if (as_tcrossfac_list) {
          terme <- NULL
        } else if (covMatrix) terme <- diag(x=terme, nrow=length(terme))
      }
    } else if ( ! is.null(isEachNewLevelInOld)) { ## non correlated effect: a vector of booleans indicating whether new is in old (qualifies cols of Cnewold)
      if ( ! is.null(fixZAlist)) { 
        Cno_InvCoo_Cof <- cov_newLv_oldv_list[[new_rd]] %id*id% t(cov_fixLv_oldv_list[[new_rd]])
        Evar <- lambda_list[[old_rd]] * (cov_newLv_fixLv_list[[new_rd]] - Cno_InvCoo_Cof) # !=0 only for (new=fiw) not in old
        terme <- newZAlist[[new_rd]] %id*id% Evar %id*id% t(fixZAlist[[new_rd]])
      } else {
        Evardiag <- lambda_list[[old_rd]] * as.numeric(! isEachNewLevelInOld)
        if (is.null(cov_newLv_newLv_list[[new_rd]])) { ## we compute only the var vector
          terme <- .calcZWZt_mat_or_diag(Z=newZAlist[[new_rd]], W=Evardiag, cholW=sqrt(Evardiag), returnMat=covMatrix)
        } else {
          if (as_tcrossfac_list) {
            terme <- newZAlist[[new_rd]] %*% Diagonal(x=sqrt(Evardiag))
          } else terme <- .calcZWZt_mat_or_diag(Z=newZAlist[[new_rd]], W=Diagonal(x=Evardiag), cholW=Diagonal(x=sqrt(Evardiag)), returnMat=covMatrix)
        }
      }
    } else { ## autocorr effect (spatial, ranCoefs...) EXCEPT IMRF (first case above)
      if ((lam_len <- length(lambda_list[[old_rd]]))>1L && any(lambda_list[[old_rd]]!=1)) { 
        # this must be a ranCoef fitted by spprec. Then we have a specific problem if there are new levels 
        # (so that diag_cov_newLv_newLv_list[[new_rd]] - diag_Cno_InvCoo_Con) is nonzero
        # (although not just a single level).
        # Then these lambda _must_ be taken into account elsewhere, _so that_ we can write
        loc_lambda <- 1
      }  else loc_lambda <- lambda_list[[old_rd]]
      if ( ! is.null(fixZAlist)) {
        Cno_InvCoo_Cof <- cov_newLv_oldv_list[[new_rd]] %gI*gI% invCov_oldLv_oldLv_list[[old_rd]] %gI*gI% t(cov_fixLv_oldv_list[[new_rd]]) # (invCoo.t() could have been precomputed.)
        Evar_C <- loc_lambda * (cov_newLv_fixLv_list[[new_rd]] - Cno_InvCoo_Cof)  
        if (inherits(Evar_C,"bigq")) Evar_C <- .mMatrix_bigq(Evar_C)
        terme <- newZAlist[[new_rd]] %id*id% Evar_C %id*id% t(fixZAlist[[new_rd]])
      } else {
        ## We're hoping for a simplification of the form Evar <- lambda_list[[old_rd]] * (1 - diag_Cno_InvCoo_Con)
        ## There are underlying assumptions:
        ## (1) The diagonal of W should be sufficient to have the diagonal of ZA % W % t(ZA);
        ##     that is, tcrossprod(ZA) is diagonal (in particular, if ZA is an incidence matrix : one 1 per row and 0 otherwise)   
        ##                       seems also OK if there are also rows of 0's (conditional ranef, (female|.))
        ##     In other cases we need cov_newLv_newLv_list (the full-cov-matrix code is OK) 
        ## (2) The '1': L is the chol factor of a correlation matrix itself with a diag of 1's
        ##     Otherwise, even if we do not need the full cov matrix, we need its diagonal.
        if (is.null(cov_newLv_newLv_list[[new_rd]])) { ## we compute only the var vector
          if (inherits(invCov_oldLv_oldLv_list[[old_rd]],"bigq")) {
            Cno_InvCoo_Con <- ( cov_newLv_oldv_list[[new_rd]]%gI*gI%  invCov_oldLv_oldLv_list[[old_rd]]) * cov_newLv_oldv_list[[new_rd]] 
            diag_Cno_InvCoo_Con <- rowSums(gmp::asNumeric(Cno_InvCoo_Con))
            # faster than:
            #nr <- nrow(tmp)
            #diag_Cno_InvCoo_Con <- numeric(nr)
            #for (it in seq_len(nr)) diag_Cno_InvCoo_Con[it] <- gmp::asNumeric(sum(tmp[it,]))
          } else diag_Cno_InvCoo_Con <- rowSums( ( cov_newLv_oldv_list[[new_rd]] %gI*gI%  invCov_oldLv_oldLv_list[[old_rd]]) * cov_newLv_oldv_list[[new_rd]] )
          Evar_C <- loc_lambda * (diag_cov_newLv_newLv_list[[new_rd]] - diag_Cno_InvCoo_Con)
          if (inherits(Evar_C,"bigq")) Evar_C <- .mMatrix_bigq(Evar_C)
          terme <- .calcZWZt_mat_or_diag(Z=newZAlist[[new_rd]],W=Evar_C,cholW=sqrt(Evar_C), returnMat=covMatrix) # READ the above comments to understand the complexity of this case
        } else { ## full cov matrix
          Cno_InvCoo_Con <- cov_newLv_oldv_list[[new_rd]] %gI*gI%  invCov_oldLv_oldLv_list[[old_rd]] %gI*gI%  t(cov_newLv_oldv_list[[new_rd]])
          Evar_C <- loc_lambda * (cov_newLv_newLv_list[[new_rd]] - Cno_InvCoo_Con)
          #cat(crayon::red("cov_newLv_newLv_list"));str(cov_newLv_newLv_list[[new_rd]])
          #cat(crayon::red("Cno_InvCoo_Con"));str(Cno_InvCoo_Con)
          if (inherits(Evar_C,"bigq")) Evar_C <- .mMatrix_bigq(Evar_C)
          if (as_tcrossfac_list) {
            eS <- eigen(Evar_C, symmetric = TRUE) # chol() typically fails and mat_sqrt() corrects more than below
            pos_ev <- (eS$values>0)
            tcrossfac <- .m_Matrix_times_Dvec(eS$vectors[,pos_ev,drop=FALSE], sqrt(eS$values[pos_ev]))
            terme <- newZAlist[[new_rd]] %*% tcrossfac
          } else terme <- .calcZWZt_mat_or_diag(Z=newZAlist[[new_rd]],W=Evar_C, cholW=try(chol(Evar_C), silent=TRUE), returnMat=covMatrix)
        }
      }
    }
    if (as_tcrossfac_list) {
      Evarlist[[new_rd]] <- terme
    } else if (covMatrix) {
      Evarlist[[new_rd]] <- as.matrix(terme)
    } else Evarlist[[new_rd]] <- terme
  }
  if ( as_tcrossfac_list) {
    return(Evarlist)
  } else {
    Evar <- Reduce("+",Evarlist)
    return(Evar)
  }
}

.calc_sliceVar <- function(slice, cov_newLv_oldv_list,
                           X.pv, newZAlist_slice, 
                           re_form_col_indices=NULL,
                           tcrossfac_beta_w_cov,
                           cov_newLv_newLv_list,
                           invCov_oldLv_oldLv_list, ## may be null if no newdata.
                           object, 
                           logdispObject, ## should remain NULL is disp not requested
                           newinold, 
                           covMatrix,blockSize,
                           diag_cov_newLv_newLv_list) {
  ## here  the problem is that newZA should map the new response levels 
  ## to the 'new' levels of random effects  
  newnrand <- length(newZAlist_slice) ## or of any other of the lists of matrices
  templateList <- vector("list", length = newnrand)
  ## subset new levels in all relevant matrices. The oldv levels must be untouched !
  locnewZA <- templateList
  locCov_newLv_oldv_list <- templateList
  if ( ! is.null(cov_newLv_newLv_list)) {locCov_newLv_newLv_list <- templateList} else locCov_newLv_newLv_list <- NULL
  if ( ! is.null(diag_cov_newLv_newLv_list)) {locDiag_cov_newLv_newLv_list <- templateList} else locDiag_cov_newLv_newLv_list <- NULL
  for (new_rd in seq_len(newnrand)) {
    req_levels <- which(colSums(newZAlist_slice[[new_rd]])>0L) 
    locnewZA[[new_rd]] <- newZAlist_slice[[new_rd]][ , req_levels, drop=FALSE] 
    locCov_newLv_oldv_list[[new_rd]] <- structure(
      cov_newLv_oldv_list[[new_rd]][req_levels, , drop=FALSE],
      isEachNewLevelInOld = attr(cov_newLv_oldv_list[[new_rd]],"isEachNewLevelInOld")[req_levels] ## for non-spatial effects; (qualifies sub cols of sub Cnewold)
    )
    #
    if ( ! is.null(locCov_newLv_newLv_list)) locCov_newLv_newLv_list[[new_rd]] <- cov_newLv_newLv_list[[new_rd]][req_levels, req_levels, drop=FALSE] 
    if ( ! is.null(locDiag_cov_newLv_newLv_list)) locDiag_cov_newLv_newLv_list[[new_rd]] <- diag_cov_newLv_newLv_list[[new_rd]][req_levels]
  }
  .calcPredVar(cov_newLv_oldv_list=locCov_newLv_oldv_list, ## needed only if newdata in predict() call
               X.pv=X.pv,## problem is that this creates the appearence of new data end more calculations
               newZAlist=locnewZA, ## either newZA or newZAC needed even if no newdata in predict() call
               re_form_col_indices=re_form_col_indices,
               tcrossfac_beta_w_cov=tcrossfac_beta_w_cov,
               cov_newLv_newLv_list=locCov_newLv_newLv_list,
               invCov_oldLv_oldLv_list=invCov_oldLv_oldLv_list, ## needed only if newdata in predict() call
               object=object,
               logdispObject=logdispObject,
               newinold=newinold,
               covMatrix=covMatrix,blockSize=blockSize,
               diag_cov_newLv_newLv_list=locDiag_cov_newLv_newLv_list,
               showpbar=FALSE)
}

.get_disp_effect_on_newZACw <- function(logdisp_cov, newZACw, fixZACw=NULL, covMatrix) {
  if (any(logdisp_cov>1e8)) {
    newZACw <- gmp::as.bigz(as.matrix(newZACw))
    logdisp_cov <- gmp::as.bigz(logdisp_cov)
    premul <- gmp::`%*%`(newZACw, logdisp_cov)
    if (covMatrix) {
      if (is.null(fixZACw)) {
        disp_effect_on_newZACw <- gmp::asNumeric(gmp::`%*%`(premul, t(newZACw)))   
      } else { # for get_predCov_var_fix(); then covMatrix is always TRUE
        fixZACw <- gmp::as.bigz(as.matrix(fixZACw))
        disp_effect_on_newZACw <- gmp::asNumeric(gmp::`%*%`(premul, t(fixZACw)))   
      }
    } else {
      disp_effect_on_newZACw <- rowSums(gmp::asNumeric(premul * newZACw))
    }
  } else {
    premul <- newZACw %*% logdisp_cov
    if (covMatrix) {
      if (is.null(fixZACw)) {
        disp_effect_on_newZACw <- premul %*% t(newZACw)  
      } else { # for get_predCov_var_fix(); then covMatrix is always TRUE
        disp_effect_on_newZACw <- premul %*% t(fixZACw)  
      }
    } else {
      disp_effect_on_newZACw <- rowSums(premul * newZACw)
    }
  }
  disp_effect_on_newZACw
}


.calcPredVar <- function(cov_newLv_oldv_list,
                         X.pv,newZAlist,
                         re_form_col_indices=NULL,
                         tcrossfac_beta_w_cov,
                         #newZACpplist=NULL,
                         cov_newLv_newLv_list=NULL,
                         invCov_oldLv_oldLv_list=NULL, ## may be null if no newdata.
                        object, 
                        logdispObject=NULL, ## should remain NULL is disp not requested
                        newinold, 
                        covMatrix=FALSE,blockSize,
                        diag_cov_newLv_newLv_list,
                        as_tcrossfac_list=FALSE,
                        showpbar
                        ) {
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
    progrbar_setup <- .set_progrbar(style = showpbar, char="s") # FIXME could implement parallee computation
    for (it in seq_len(nslices)) {
      slice <- (slices[it]+1L):slices[it+1L]
      X.pv_slice <- X.pv[slice,,drop=FALSE]
      for (rt in seq_along(newZAlist)) newZAlist_slice[[rt]] <- newZAlist[[rt]][slice,,drop=FALSE]
      predVar[[it]] <- .calc_sliceVar( ## -> recursive call to .calcPredVar 
                                      slice=slice,
                                      cov_newLv_oldv_list=cov_newLv_oldv_list,
                                      X.pv=X.pv_slice,
                                      newZAlist_slice=newZAlist_slice,
                                      re_form_col_indices=re_form_col_indices,
                                      tcrossfac_beta_w_cov=tcrossfac_beta_w_cov,
                                      cov_newLv_newLv_list=cov_newLv_newLv_list,
                                      invCov_oldLv_oldLv_list=invCov_oldLv_oldLv_list, ## may be null if no newdata.
                                      object=object, 
                                      logdispObject=logdispObject, ## should remain NULL is disp not requested
                                      newinold=newinold, 
                                      covMatrix=covMatrix,blockSize=blockSize,
                                      diag_cov_newLv_newLv_list=diag_cov_newLv_newLv_list) # (! covMatrix) hence as_tcrossfac_list must be FALSE
      if (showpbar) progrbar_setup$progress(slices[it+1L]/nrX) ## update progress bar
    }
    if (showpbar) close(progrbar_setup$pb)
    return(unlist(predVar))
  }
  ############################################################
  newZACvar <- .calc_newZACvar(newZAlist,cov_newLv_oldv_list)
  if (ncol(X.pv)>ncol(object$X.pv)) { # detection etaFix$beta case
    if (is.matrix(newZACvar)) {
      XZAC <- cbind2(as.matrix(X.pv)[,colnames(object$X.pv),drop=FALSE],newZACvar)
    } else XZAC <- cbind2(X.pv[,colnames(object$X.pv),drop=FALSE],newZACvar) ## no such thing in point pred as fix vs. estim does not matter there.
  } else if (is.matrix(newZACvar)) { # when C is dense, cbind2(X.pv,newZACvar) is quite dense and further operations from it may handle inefficient dge
    XZAC <- cbind2(as.matrix(X.pv),newZACvar)
  } else XZAC <- cbind2(X.pv,newZACvar) ## uses C (in C.w) rather than L (in L.v)
  ## First component of predVar
  # variance of expectation of Xbeta+Zb due to var of (hat(beta),old hat(v)) using E[b] as function of old hat(v)
  if (as_tcrossfac_list) { # covMatrix + is.null(invCov_oldLv_oldLv_list)
    if ( ! covMatrix ) stop("invalid call with covMatrix=FALSE and as_tcrossfac_list=TRUE")
    predVar <- list("1"= XZAC %*% tcrossfac_beta_w_cov)
  } else {
    predVar <- .calcZWZt_mat_or_diag(Z=XZAC, W=NULL, 
                                     tcrossfac_W=tcrossfac_beta_w_cov, returnMat=covMatrix) ## component for linPred=TRUE,disp=FALSE when newdata=ori data
  }
  #cat(crayon::red("predVar"));str(predVar)
  ## Second component of predVar: 
  # Evar: expect over distrib of (hat(beta),new hat(v)) of [variance of Xbeta+Zb given (hat(beta),old hat(v))]
  # Evar must be 0 when newdata=ori data
  if ( ! is.null(invCov_oldLv_oldLv_list)) { ## must be equivalent to the presence of newdata
    Evar <- .calc_Evar(newZAlist=newZAlist,newinold=newinold, cov_newLv_oldv_list=cov_newLv_oldv_list, 
               invCov_oldLv_oldLv_list=invCov_oldLv_oldLv_list, 
               cov_newLv_newLv_list=cov_newLv_newLv_list, covMatrix=covMatrix,
               diag_cov_newLv_newLv_list=diag_cov_newLv_newLv_list, object=object,
               as_tcrossfac_list=as_tcrossfac_list)
    if (as_tcrossfac_list) { # in that case Evar should be a list with elements for each ranef (except IMRF whose Evar contribution is 0)
      predVar <- c(predVar,Evar) # join lists
    } else predVar <- predVar + Evar ## component for linPred=TRUE,disp=FALSE whether newdata=ori data or not
  }
  #cat(crayon::red("Evar"));str(Evar)
  # If components for uncertainty in dispersion params were requested,
  #   logdispObject is not NULL
  # If some components are computable, $$dwdlogdisp should not be NULL
  # Former approach (changed 08/2016) was to test logdispObject and then 
  #   for any 'problem'. But there may some 'problem' and still a valid logdispObject
  # => In new version, dwdlogdisp should be either NULL or a conforming matrix;
  #  'problems" should not be tested.
  if ( ! is.null(logdispObject$dwdlogdisp) ) {
    dwdlogdisp <- logdispObject$dwdlogdisp ## original ranefs in original order (indep of re.form)
    logdisp_cov <- logdispObject$logdisp_cov ## idem
    phi_cols <- attr(dwdlogdisp,"col_info")$phi_cols ## make local copy before subsetting the matrix!
    if ( ! is.null(re_form_col_indices) ) { ## selection of blocks for re.form ranefs 
      whichcols <- c(re_form_col_indices$which_ranef_cols, phi_cols)
      dwdlogdisp <- dwdlogdisp[re_form_col_indices$subrange,whichcols] ## permuted ranefs => permuted rows and cols
      logdisp_cov <- logdisp_cov[whichcols,whichcols] ## permuted ranefs => permuted rows and cols
    } 
    #if new ranef are permuted wrt old ranefs, 
    # newZACvar rows are unchanged (they correspond to observations),
    # but its cols and the *local* dwdlogdisp's rows and cols have just been permuted 
    # hence newZACw rows will be unchanged bt its cols will be permuted
    if ( ! is.null(hyper_info <- .get_from_ranef_info(object,"hyper_info"))) {
      summingMat <- hyper_info$summingMat
      if ( ! is.null(re_form_col_indices) ) {
        summingMat <- summingMat[newinold,,drop=FALSE]
        colids <- numeric(nrow(summingMat))
        for (it in seq(nrow(summingMat))) colids[it] <- which(summingMat[it,]>0) # identifiy active cols for each row
        colids <- unique(colids) # but not sort()
        summingMat <- summingMat[,colids, drop=FALSE]
      }
      summingMat <- as.matrix(Matrix::bdiag(summingMat,rep(1,length(phi_cols))))
      dwdlogdisp <- dwdlogdisp %*% summingMat
      logdisp_cov <- t(summingMat) %*% logdisp_cov %*% summingMat
    }
    newZACw <- newZACvar %*% dwdlogdisp ## typically (nnew X (n_u_h*nrand)) %*% ((n_u_h*nrand) X (nrand+1)) = nnew * 2 hence small 
    if (as_tcrossfac_list) { 
      predVar[["tcross_disp_effect"]] <- newZACw %*% mat_sqrt(logdisp_cov)
    } else {
      disp_effect_on_newZACw <- .get_disp_effect_on_newZACw(logdisp_cov, newZACw, covMatrix=covMatrix)
      predVar <- predVar + disp_effect_on_newZACw
    }
  }
  return(predVar) ## may be a Matrix
}


.calcResidVar <- function(object,newdata=NULL, phi.object=object$phi.object, families=object$families, 
                          mv_it=NULL, #  used to pass non-default to .get_glm_phi()
                          family=object$family, fv,
                          nobs_info=nrow(object$X.pv)) {
  if ( ! is.null(families)) {
    residVars <- vector("list", length(families))
    cum_nobs <- attr(families,"cum_nobs")
    for (mv_it in seq_along(families)) {
      residVars[[mv_it]] <- .calcResidVar(object, newdata=newdata, phi.object=phi.object[[mv_it]], families=NULL, 
                                          mv_it=mv_it, 
                                          family=families[[mv_it]], fv=fv[.subrange(cumul=cum_nobs, it=mv_it)],
                                          nobs_info=object$vec_nobs[mv_it])
    }
    residVar <- unlist(residVars, recursive = FALSE, use.names = FALSE)
  } else {
    if (family$family %in% c("poisson","binomial","COMPoisson","negbin")) {
      residVar <- family$variance(fv)
    } else {
      if (is.null(phi_outer <- phi.object$phi_outer)) { ## valid whether newdata are NULL or not:      
        glm_phi <- .get_glm_phi(object, mv_it)
        residVar <- predict(glm_phi, newdata=newdata, type="response")
      } else { ## phi, but not glm_phi
        if (length(phi_outer)==1L) {
          if (is.null(newdata)) {
            nobs <- nobs_info            
          } else nobs <- nrow(newdata)
          residVar <- rep(phi_outer,nobs)    # __F I X M E   N O T__ This shortcut certainly does not handle prior.weights (but this is consistent with doc)        
        } else { # e.g. resid.model = list(formula=~0+offset(logphi)) example => there is a phi_outer vector,
          # so .get_glm_phi() was not previously called. But even for such offset .get_glm_phi() can return an object of class "glm"
          # so that one can predict() from it. This recycles existing code caring for links, prior weights etc...
          # Otherwise a more direct call to model.frame(.get_phiform(object, mv_it), data=newdata) might perhaps be used.
          glm_phi <- .get_glm_phi(object, mv_it) # look whether phi.Fix was set by a phiGLM with only an offset
          if (is.null(glm_phi)) {
            stop(paste0("Unable to compute 'residVar' given length(<fixed phi>)!=1L.\n", 
                        "A 'resid.model' argument, perhaps with only an offset term, might be the only way forward.")) 
          } else residVar <- predict(glm_phi, newdata=newdata, type="response")
        }
      }
    }
    if (family$family=="Gamma") residVar <- residVar * fv^2
  }
  residVar
} 


if (FALSE) {# This fn is documentation, 
  .calcResidVar_phiW <- function(object,newdata=NULL, phi.object=object$phi.object, families=object$families, 
                                 mv_it=NULL, #  used to pass non-default to .get_glm_phi()
                                 family=object$family, fv,
                                 nobs_info=nrow(object$X.pv),
                                 prior.weights=object$prior.weights,
                                 phi_mod_class=object$models[["phi"]]) {
    if ( ! is.null(families)) {
      residVars <- vector("list", length(families))
      cum_nobs <- attr(families,"cum_nobs")
      for (mv_it in seq_along(families)) {
        residVars[[mv_it]] <- .calcResidVar(object, newdata=newdata, phi.object=phi.object[[mv_it]], families=NULL, 
                                            mv_it=mv_it, 
                                            family=families[[mv_it]], fv=fv[.subrange(cumul=cum_nobs, it=mv_it)],
                                            nobs_info=object$vec_nobs[mv_it],
                                            prior.weights=prior.weights[[mv_it]],
                                            phi_mod_class=phi_mod_class[mv_it])
      }
      residVar <- unlist(residVars, recursive = FALSE, use.names = FALSE)
    } else {
      if (family$family %in% c("poisson","binomial","COMPoisson","negbin")) {
        residVar <- family$variance(fv)
      } else {
        if (is.null(phi_outer <- phi.object$phi_outer)) { ## valid whether newdata are NULL or not:      
          glm_phi <- .get_glm_phi(object, mv_it)
          residVar <- predict(glm_phi, newdata=newdata, type="response")
        } else { ## phi, but not glm_phi
          if (length(phi_outer)==1L) {
            if (is.null(newdata)) {
              nobs <- nobs_info            
            } else nobs <- nrow(newdata)
            residVar <- rep(phi_outer,nobs)    # __F I X M E   N O T__ This shortcut certainly does not handle prior.weights (but this is consistent with doc)        
          } else { # e.g. resid.model = list(formula=~0+offset(logphi)) example => there is a phi_outer vector,
            # so .get_glm_phi() was not previously called. But even for such offset .get_glm_phi() can return an object of class "glm"
            # so that one can predict() from it. This recycles existing code caring for links, prior weights etc...
            # Otherwise a more direct call to model.frame(.get_phiform(object, mv_it), data=newdata) might perhaps be used.
            glm_phi <- .get_glm_phi(object, mv_it) # look whether phi.Fix was set by a phiGLM with only an offset
            if (is.null(glm_phi)) {
              stop(paste0("Unable to compute 'residVar' given length(<fixed phi>)!=1L.\n", 
                          "A 'resid.model' argument, perhaps with only an offset term, might be the only way forward.")) 
            } else residVar <- predict(glm_phi, newdata=newdata, type="response")
          }
        }
        ## isoscape(raster = ElevRasterDE, isofit = GermanFit) is a good test. 26 points were fitted with prior weights
        ## prediction is by blocks of ~1000 points and one wants the response variance for a single obs. Prior weights are not useful here
        ## And the next call fails if the prior weights do not have a matching length ~1000.
        essai <- .get_phiW(object, newdata, dims=c(nobs,1L),
                           prior.weights=prior.weights, phi_type="predict", needed=1L,
                           phi_mod_class=phi_mod_class, mv_it)[,1]
        if (diff(range(essai-residVar))>1e-8) {browser()} else {cat(crayon::red("OK"),"\n")}
      }
      if (family$family=="Gamma") residVar <- residVar * fv^2
    }
    residVar
  } 
}


.get_oldlevels <- function(object, old_rd, fix_info) {
  if ( ! is.null(fix_info)) {
    oldlevels <- colnames(fix_info$newZAlist[[old_rd]]) ## old_rd despite the "newZAlist" name: specific to fix_info
  } else oldlevels <- colnames(object$ZAlist[[old_rd]])
  oldlevels
}

.assign_newLv_for_newlevels_corrMatrix <- function(newoldC, newlevels, newLv_env, new_rd, corr.model, which_mats, ranefs) {
  ## this is called if re.form... or pdep_effects() but need for newdata too (test by permuting data)
  # "test-predVar-Matern-corrMatrix" shows newdata working with corrMatrix, when newlevels are within oldlevels 
  oldlevels <- rownames(newoldC)
  if (is.null(oldlevels)) stop("is.null(rownames(newoldC))")
  if (length(setdiff(newlevels,oldlevels))) stop(paste0("Found new levels for a '",corr.model,"' random effect."))
  if (which_mats$no) newLv_env$cov_newLv_oldv_list[[new_rd]] <- structure(newoldC[newlevels, ,drop=FALSE],
                                                                          ranefs=ranefs[[new_rd]])
  if (which_mats$nn[new_rd]) {
    newLv_env$cov_newLv_newLv_list[[new_rd]] <- newoldC[newlevels, newlevels ,drop=FALSE]
  } else {
    # diag(x=newoldC,names=TRUE)[newlevels] # does not work: names=TRUE is ignored
    tmp <- diag(x=newoldC)
    names(tmp) <- oldlevels
    newLv_env$diag_cov_newLv_newLv_list[[new_rd]] <- tmp[newlevels]
  }
  # no reurn value, the newLv_env ahs been modified
}

.make_new_corr_lists <- function(object,
                                 locdata, ## a list of arrays for selected ranefs, with selected coordinates, if 'preprocessed' by get_predCov_var_fix(), 
                                          ## or a single data frame with all coordinates for all ranefs if from within .calc_new_X_ZAC()
                                 which_mats, newZAlist, newinold, fix_info=NULL,
                                 invCov_oldLv_oldLv_list) {
  newLv_env <- new.env(parent=emptyenv())
  if (which_mats$no) {
    newLv_env$cov_newLv_oldv_list <- vector("list",length(newZAlist)) ## declar-initialization, will be filled in the loop
  } else newLv_env$cov_newLv_oldv_list <- .make_corr_list(object,newZAlist=NULL) # we always need a non trivial value
  newLv_env$cov_newLv_newLv_list <- vector("list",length(newLv_env$cov_newLv_oldv_list))
  newLv_env$diag_cov_newLv_newLv_list <- vector("list",length(newLv_env$cov_newLv_oldv_list))
  ranefs <- attr(newZAlist,"exp_ranef_strings") 
  ## The following comment may be obsolete: ## only the matrices that get a "ranefs" attribute will be used in .compute_ZAXlist()
  if (any(unlist(which_mats))) {
    strucList <- object$strucList
    spatial_terms <- attr(object$ZAlist,"exp_spatial_terms")
    for (new_rd in seq_along(newLv_env$cov_newLv_oldv_list)) {
      old_rd <- newinold[new_rd]
      if ( which_mats$no || which_mats$nn[new_rd]) {
        corr.model <- attr(strucList[[old_rd]],"corr.model")
        if (is.null(corr.model)) {
          if ( ! is.null(fix_info)) {
            if (which_mats$no) newLv_env$cov_newLv_oldv_list[[new_rd]] <- 
              .calc_sub_diagmat_cov_newLv_oldv(oldZA=fix_info$newZAlist[[old_rd]], newZA=newZAlist[[new_rd]],
                                               namesTerms=attr(newZAlist,"namesTerms")[[new_rd]] ## should have length one, which is all that matters
              ) 
          } else { ## else the list elements remained NULL... until .calc_Var_given_fixef() needed them
            if (which_mats$nn[new_rd]) {
              newLv_env$cov_newLv_newLv_list[[new_rd]] <- 1 ## just 1 must suffice except if we subsetted (slice...) in which case we would need the names as in:
              # zut <- .symDiagonal(n=ncol(newZAlist[[old_rd]])) 
              # dimnames(zut) <- list(colnames(newZAlist[[old_rd]]),colnames(newZAlist[[old_rd]])) 
              # cov_newLv_newLv_list[[new_rd]] <- zut
              ## BUT slicing does not occur when any which_mats$nn[new_rd] is true. 
            } else { # occurs even in simulate if (new)ZAlist is an incidence matrix
              # and further tested in particular by 'HACorn' in test-predVar.
              newLv_env$diag_cov_newLv_newLv_list[[new_rd]] <- rep(1,ncol(newZAlist[[old_rd]])) # allowing subsetting
            }
          }
        } else if ( corr.model=="random-coef") {
          namesTerms <- attr(newZAlist,"namesTerms")[[new_rd]]
          if ( ! is.null(fix_info)) {
            oldlevels <- colnames(fix_info$newZAlist[[old_rd]]) ## old_rd despite the "newZAlist" name: specific to fix_info
          } else oldlevels <- colnames(object$ZAlist[[old_rd]])
          newlevels <- colnames(newZAlist[[new_rd]])
          oldornew <- unique(c(oldlevels,newlevels))
          numnamesTerms <- length(namesTerms) ## 2 for random-coef
          repnames <- rep(oldornew, numnamesTerms)
          newcols <- repnames %in% newlevels ## handle replicates (don't `[` newoldC using names !)
          oldcols <- repnames %in% oldlevels
          if ( ! is.null(kron_Y <- object$ranef_info$sub_corr_info$kron_Y_LMatrices[[old_rd]])) { # 
            namesTerms <- attr(newZAlist,"namesTerms")[[new_rd]]
            # if ( ! is.null(fix_info)) {
            #   oldlevels <- colnames(fix_info$newZAlist[[old_rd]]) ## old_rd despite the "newZAlist" name: specific to fix_info
            # } else oldlevels <- colnames(object$ZAlist[[old_rd]])
            # oldlevels <- unique(oldlevels)
            # newlevels <- unique(colnames(newZAlist[[new_rd]]))
            # # "test-predVar-Matern-corrMatrix" shows newdata working with corrMatrix, when newlevels are within oldlevels 
            # if (length(setdiff(newlevels,oldlevels))) stop(paste0("Found new levels for a '",corr.model,"' random effect."))
            # kron_Y <- .tcrossprod(kron_Y, perm=TRUE) ##  Can reconstruct permuted (consistent with perm of cols of Z) corrMatrix from its CHM factor
            # colnames(kron_Y) <- rownames(kron_Y) <- oldlevels
            kron_Y <- .tcrossprod(kron_Y, perm=TRUE) ##  Can reconstruct permuted (consistent with perm of cols of Z) corrMatrix from its CHM factor
            oldlevels <- colnames(kron_Y)
            newlevels <- unique(colnames(newZAlist[[new_rd]]))
            # "test-predVar-Matern-corrMatrix" shows newdata working with corrMatrix, when newlevels are within oldlevels 
            if (length(setdiff(newlevels,oldlevels))) stop(paste0("Found new levels for a '",corr.model,"' random effect."))
            # all newlevels must be in oldlevels hence the construction of the matrix is a little simplified
            compactcovmat <- attr(strucList[[old_rd]],"latentL_blob")$compactcovmat
            newoldC <- .makelong_kronprod(compactcovmat, kron_Y=kron_Y[newlevels,oldlevels]) 
            #colnames(newoldC) <- rep(oldlevels, numnamesTerms) ## these will be needed by .match_old_new_levels() or .assign_newLv_for_newlevels_corrMatrix()
            #rownames(newoldC) <- rep(newlevels, numnamesTerms)
            if (which_mats$no) {
              newLv_env$cov_newLv_oldv_list[[new_rd]] <- structure(newoldC, ranefs=ranefs[[new_rd]])
              # if (inherits(newoldC,"bigq")) attr(newLv_env$cov_newLv_oldv_list[[new_rd]], 
              #                                    "DIMNAMES") <- list(repnames[newcols],repnames[oldcols]) ## these will be needed by .match_old_new_levels()
            }
            newnewC <- .makelong(compactcovmat, longsize=length(newlevels)*numnamesTerms, kron_Y=kron_Y[newlevels,newlevels]) 
            if (which_mats$nn[new_rd]) {
              newLv_env$cov_newLv_newLv_list[[new_rd]] <- newnewC #names presumably not used
              #if (inherits(newoldC,"bigq")) attr(diag_cov_newLv_newLv_list[[new_rd]], "DIMNAMES") <- list(repnames[newcols],repnames[newcols])
            } else {
              nc <- ncol(newnewC)
              diagPos <- seq.int(1L,nc^2,nc+1L) ## it's not only an efficient syntax, that's the only way to extract the diag from bigq...
              #if (inherits(newoldC,"bigq")) {
              #  newLv_env$diag_cov_newLv_newLv_list[[new_rd]] <- newoldC[[]][diagPos[newcols]] ## with the [[]] to access the raw vector
              #} else 
              newLv_env$diag_cov_newLv_newLv_list[[new_rd]] <- newnewC[diagPos]
            }
          } else {
            if ( ! is.null(gmp_compactcovmat <- attr(strucList[[old_rd]],"latentL_blob")$gmp_compactcovmat)) {
              newoldC <- .makelong(gmp_compactcovmat, longsize=length(oldornew)*numnamesTerms, kron_Y=kron_Y) 
              #attr(newoldC,"DIMNAMES") <- list(repnames,repnames) # will be lost in next subsetting
            } else {
              compactcovmat <- attr(strucList[[old_rd]],"latentL_blob")$compactcovmat
              newoldC <- .makelong(compactcovmat, longsize=length(oldornew)*numnamesTerms, kron_Y=kron_Y) 
            }
            #
            if (inherits(newoldC,"bigq")) {
              attr(newoldC, "DIMNAMES") <- list(repnames,repnames) 
            } else colnames(newoldC) <- rownames(newoldC) <- repnames ## these will be needed by .match_old_new_levels() or .assign_newLv_for_newlevels_corrMatrix()
            if (which_mats$no) {
              # using which(newcols) etc suffice because there is so much redundant info in the newoldC matrix
              newLv_env$cov_newLv_oldv_list[[new_rd]] <- structure(newoldC[which(newcols),which(oldcols),drop=FALSE],
                                                                   ranefs=ranefs[[new_rd]])
              if (inherits(newoldC,"bigq")) attr(newLv_env$cov_newLv_oldv_list[[new_rd]], 
                                                 "DIMNAMES") <- list(repnames[newcols],repnames[oldcols]) ## these will be needed by .match_old_new_levels()
            }
            if (which_mats$nn[new_rd]) {
              newLv_env$cov_newLv_newLv_list[[new_rd]] <- newoldC[which(newcols),which(newcols),drop=FALSE] ## gmp bug report sent a long time ago
              #names presumably not used
              #if (inherits(newoldC,"bigq")) attr(diag_cov_newLv_newLv_list[[new_rd]], "DIMNAMES") <- list(repnames[newcols],repnames[newcols])
            } else {
              nc <- ncol(newoldC)
              diagPos <- seq.int(1L,nc^2,nc+1L) ## it's not only an efficient syntax, that's the only way to extract the diag from bigq...
              if (inherits(newoldC,"bigq")) {
                newLv_env$diag_cov_newLv_newLv_list[[new_rd]] <- newoldC[[]][diagPos[newcols]] ## with the [[]] to access the raw vector
              } else newLv_env$diag_cov_newLv_newLv_list[[new_rd]] <- newoldC[diagPos[newcols]]
            }
            #cat(crayon::red("newoldC"));str(newoldC)
          }
        } else if (corr.model %in% c("corrMatrix","adjacency")) {
          newoldC <- .tcrossprod(object$strucList[[old_rd]], perm=TRUE) ##  Can reconstruct permuted (consistent with perm of cols of Z) corrMatrix from its CHM factor
          colnames(newoldC) <- rownames(newoldC) <- .get_oldlevels(object, old_rd, fix_info)
          .assign_newLv_for_newlevels_corrMatrix(newoldC, newlevels=colnames(newZAlist[[new_rd]]), 
                                                 newLv_env, new_rd, corr.model, which_mats, ranefs)
          # assign to newLv_env$cov_newLv_oldv_list, cov_newLv_newLv_list, diag_cov_newLv_newLv_list
        } else if (corr.model == "IMRF") {
          ## in this case a new A matrix (by .get_new_AMatrices()) must be computed (elsewhere: it goes in newZAlist)
          ## but the corr matrix between the nodes is unchanged as node positions do not depend on response position => nn=no
          if (which_mats$no || which_mats$nn[new_rd]) uuColdold <- .tcrossprod(object$strucList[[old_rd]], perm=TRUE) # perm=TRUE means for ranefs permuted as in ZA
          if (which_mats$no) newLv_env$cov_newLv_oldv_list[[new_rd]] <- structure(uuColdold,ranefs=ranefs[[new_rd]]) ## always necess for .compute_ZAXlist(XMatrix = cov_newLv_oldv_list, ZAlist = newZAlist)
          ## correct code but should be useless:
          # warning("Suspect code: covariance matrix computation requested for IMRF term.")
          # if (which_mats$nn[new_rd]) {
          #   cov_newLv_newLv_list[[new_rd]] <- uuColdold
          # } else diag_cov_newLv_newLv_list[[new_rd]] <- diag(x=uuColdold)
        } else { ## all models where a correlation matrix must be computed from a distance matrix => $calc_corr_from_dist() needed
          old_char_rd <- as.character(old_rd)
          if ( ! is.null(fix_info)) {
            info_olduniqueGeo <- fix_info$newuniqueGeo
          } else {
            info_olduniqueGeo <- attr(object,"info.uniqueGeo") 
          }
          # olduniqueGeo needed in all cases
          if ( ! is.array(info_olduniqueGeo)) { ## test TRUE for attr(object,"info.uniqueGeo") for version > 2.3.18:
            olduniqueGeo <- info_olduniqueGeo[[old_char_rd]]
          } else olduniqueGeo <- info_olduniqueGeo 
          if ( is.data.frame(locdata)) {
            geonames <- colnames(olduniqueGeo) 
            newuniqueGeo <- locdata[,geonames,drop=FALSE] ## includes nesting factor
          } else { ## locdata is 'preprocessed' list of arrays (tested by get_predCov_var_fix() tests)
            newuniqueGeo <- locdata[[as.character(old_rd)]] ## preprocessed, [,geonames,drop=FALSE] not necess ## includes nesting factor 
            geonames <- colnames(newuniqueGeo)
          }
          ## distance matrix and then call to correl fn:
          control_dist_rd <- .get_control_dist(object, old_char_rd)
          if (corr.model=="AR1") {
            ### older, non nested code:
            # if (which$no) resu$uuCnewold <- proxy::dist(newuniqueGeo,olduniqueGeo) 
            # if (which$nn) resu$uuCnewnew <- proxy::dist(newuniqueGeo)
            ### new code recycling .get_dist_nested_or_not, but not fastest nor most transparent (fixme?)
            onGeo <- rbind(newuniqueGeo,olduniqueGeo) # includes nesting factor
            if (object$spaMM.version < "2.2.118") {
              blob <- .get_dist_nested_or_not(object$spatial_term, data=onGeo, distMatrix=NULL, uniqueGeo=NULL, 
                                              dist.method=control_dist_rd$dist.method, as_matrix=TRUE,
                                              needed=c(distMatrix=TRUE), geo_envir=NULL)
            } else {
              blob <- .get_dist_nested_or_not(spatial_terms[[old_rd]], data=onGeo, distMatrix=NULL, 
                                              uniqueGeo=NULL, ## FIXME provide uniqueGeo 'to save time'(??) ?
                                            dist.method=control_dist_rd$dist.method, as_matrix=TRUE,
                                            needed=c(distMatrix=TRUE), geo_envir=NULL) 
            }
            ## we merged old and new so need to get the respective cols (which may overlap) 
            uli_onGeo <- .ULI(onGeo) # this should give row and columns in the blob ## FIXME how to make sure of that? .get_dist_nested_or_not must use .ULI()
            uli_new <- uli_onGeo[seq(nrow(newuniqueGeo))]
            uli_old <- uli_onGeo[-seq(nrow(newuniqueGeo))]
            if (which_mats$no) uuCnewold <- blob$distMatrix[uli_new,uli_old,drop=FALSE] ## rows match the newZAlist, cols match th u_h 
            if (which_mats$nn[new_rd]) {
              uuCnewnew <- blob$distMatrix[uli_new,uli_new,drop=FALSE]
            }
          } else if (corr.model %in% c("Cauchy", "Matern")) {
            ### rho only used to compute scaled distances
            rho <- .get_cP_stuff(object$ranFix,"rho", which=old_char_rd)
            if ( ! is.null(rho_mapping <- control_dist_rd$rho.mapping) 
                 && length(rho)>1L ) rho <- .calc_fullrho(rho=rho,coordinates=geonames,rho_mapping=rho_mapping)
            # : if control_dist_rd comme from the call (vs from moreargs) rho_mapping may still be NULL 
            #     and then the code assumes that calling  .calc_fullrho() is not necessary (that seems OK) 
            ## rows from newuniqueGeo, cols from olduniqueGeo:
            txt <- paste(c(spatial_terms[[old_rd]][[2]][[3]])) ## the RHS of the ( . | . ) # c() to handle very long RHS
            if (length(rho)>1L && length(grep("%in%",txt))) { # nested geostatistical effect
              coord_within <- .extract_check_coords_within(spatial_term=spatial_terms[[old_rd]]) 
              msd.arglist <- list(uniqueGeo=newuniqueGeo[,coord_within,drop=FALSE],
                                  uniqueGeo2=olduniqueGeo[,coord_within,drop=FALSE],
                                  rho=rho,return_matrix=TRUE)
            } else msd.arglist <- list(uniqueGeo=newuniqueGeo,uniqueGeo2=olduniqueGeo,
                                rho=rho,return_matrix=TRUE)
            # If control_dist_rd$dist.method is NULL, do not over-write the non-NULL default of make_scaled_dist():
            if ( ! is.null(dist.method <- control_dist_rd$dist.method)) msd.arglist$dist.method <- dist.method 
            if (which_mats$no) uuCnewold <- do.call(make_scaled_dist,msd.arglist) ## ultimately allows products with Matrix ## '*cross*dist' has few methods, not even as.matrix
            if (which_mats$nn[new_rd])  {
              msd.arglist$uniqueGeo2 <- NULL
              if (nrow(msd.arglist$uniqueGeo)==1L) {
                uuCnewnew <- matrix(0) ## trivial _distance_ matrix for single point (converted to trivial cov below!)
              } else uuCnewnew <- do.call(make_scaled_dist,msd.arglist) 
            }
            if (length(grep("%in%",txt))) { # nested geostatistical effect
              onGeo <- rbind(newuniqueGeo,olduniqueGeo) # includes nesting factor
              isInf <- .get_dist_nested_or_not(spatial_term=spatial_terms[[old_rd]], 
                                               data=onGeo, distMatrix=NULL, 
                                               uniqueGeo=NULL, 
                                               dist.method = dist.method,needed=c(notSameGrp=TRUE),
                                               geo_envir=NULL)$notSameGrp
              ## we merged old and new so need to get the respective cols (which may overlap) 
              uli_onGeo <- .ULI(onGeo) # this should give row and columns in the blob ## FIXME how to make sure of that? .get_dist_nested_or_not must use .ULI()
              uli_new <- uli_onGeo[seq(nrow(newuniqueGeo))]
              uli_old <- uli_onGeo[-seq(nrow(newuniqueGeo))]
              if (which_mats$no) {
                isInfno <- isInf[uli_new,uli_old,drop=FALSE]
                uuCnewold[isInfno] <- Inf ## ultimately allows products with Matrix ## '*cross*dist' has few methods, not even as.matrix
              }
              if (which_mats$nn[new_rd]) {
                isInfnn <- isInf[uli_new,uli_new,drop=FALSE]
                uuCnewnew[isInfnn] <- Inf
              }
            }
          } else stop("Unhandled corr.model.")
          if (object$spaMM.version<"2.4.49") {
            if (which_mats$no) newLv_env$cov_newLv_oldv_list[[new_rd]] <- structure(.calc_corr_from_dist(uuCnewold, object, corr.model,char_rd=old_char_rd),
                                                                          corr.model=corr.model,
                                                                          ranefs=ranefs[[new_rd]])
            if (which_mats$nn[new_rd]) newLv_env$cov_newLv_newLv_list[[new_rd]] <- .calc_corr_from_dist(uuCnewnew, object, corr.model,char_rd=old_char_rd)
          } else {
            if (which_mats$no) {
              cov_newLv_oldv <- .get_from_ranef_info(object)$corr_families[[old_rd]]$calc_corr_from_dist(
                ranFix=object$ranFix, char_rd=old_char_rd, distmat=uuCnewold)
              # if (attr(strucList[[old_rd]],"need_gmp") && 
              if (inherits(invCov_oldLv_oldLv_list[[old_rd]],"bigq")) cov_newLv_oldv <- structure(gmp::as.bigq(cov_newLv_oldv),
                                                                                    DIMNAMES=dimnames(cov_newLv_oldv))
              newLv_env$cov_newLv_oldv_list[[new_rd]] <- structure(cov_newLv_oldv, corr.model=corr.model, ranefs=ranefs[[new_rd]])
            }
            if (which_mats$nn[new_rd]) {
              newLv_env$cov_newLv_newLv_list[[new_rd]] <- 
                .get_from_ranef_info(object)$corr_families[[old_rd]]$calc_corr_from_dist(
                  ranFix=object$ranFix, char_rd=old_char_rd, distmat=uuCnewnew)
            } else {
              newLv_env$diag_cov_newLv_newLv_list[[new_rd]] <- rep(1,nrow(newuniqueGeo)) # just 1 must suffice except when we subset (slice...)
              #   .get_from_ranef_info(object)$corr_families[[old_rd]]$calc_corr_from_dist(
              # ranFix=object$ranFix, char_rd=old_char_rd, distmat=diag(x=uuCnewnew))
            }
          }
        }
      }
    }
  }
  return(newLv_env)
}

.process_variances <- local({
  link_warned <- FALSE
  function(variances, object) {
    if ( (! link_warned) && identical(variances$predVar, TRUE)) {
      if (is.null(families <- object$families)) {
        any_nonid_link <- object$family$link!="identity"
      } else for (mv_it in seq_along(families)) if (any_nonid_link <- (families[[mv_it]]$link!="identity")) {break}
      if (any_nonid_link) {
        link_warned <<- TRUE
        message("Non-identity link: predVar is on linear-predictor scale.") #        [This message appears only once per session]")
      }
    }
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
    if (is.null(variances$as_tcrossfac_list)) variances$as_tcrossfac_list <- FALSE 
    if (is.null(variances$cov)) variances$cov <- variances$as_tcrossfac_list
    if (is.null(variances$naive)) variances$naive <- FALSE
    if (is.null(variances$cancel_X.pv)) variances$cancel_X.pv <- FALSE
    return(variances)
  }
})

.match_old_new_levels <- function(new_rd, new_X_ZACblob,
                                  old_cum_n_u_h, w_h_coeffs, rand.families) {
  old_rd <- new_X_ZACblob$newinold[new_rd]
  oldu.range <- (old_cum_n_u_h[old_rd]+1L):(old_cum_n_u_h[old_rd+1L])
  if (old_rd %in% new_X_ZACblob$spatial_old_rd) { 
    return(w_h_coeffs[oldu.range])          
  } else {
    oldlevels <- colnames(new_X_ZACblob$subZAlist[[new_rd]]) 
    newlevels <- colnames(new_X_ZACblob$newZACpplist[[new_rd]])
    interlevels <- intersect(oldlevels,newlevels)
    oldpos <- which(oldlevels %in% interlevels) ## positions: handle replicates for random-coef
    newpos <- which(newlevels %in% interlevels)
    oldv <- w_h_coeffs[oldu.range]
    names(oldv) <- oldlevels
    psi_M <- switch(attr(rand.families,"lcrandfamfam")[old_rd], # __F I X M E__ I could start to get rid of lcrandfamfam ?
                    gaussian = 0,
                    gamma = 1, 
                    beta = 1/2, 
                    "inverse.gamma" = 1
    )
    vpsi_M <- rand.families[[old_rd]]$linkfun(psi_M) 
    ## since vpsi_M can be non-zero, the expectation of the response can be modified in a re.form model compared to the original
    newv <- rep(vpsi_M,length(newlevels)) ## fills new levels with psi_M
    names(newv) <- newlevels
    newv[newpos] <- oldv[oldpos] 
    return(newv)
  }
}

.eta_linkfun <- function(mu_T, # must be vector, even in mv case 
                         family, families=NULL) { 
  if (! is.null(families)) {
    cum_nobs <- attr(families,"cum_nobs")
    eta <- vector("list", length(cum_nobs)-1L)
    for (mv_it in seq_along(eta)) {
      resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
      eta[[mv_it]] <- .eta_linkfun(mu_T[resp_range], family=families[[mv_it]]) 
    }
    return(unlist(eta, recursive = FALSE, use.names = TRUE))
    # : a redundant object with a list 'mv' of sub-muetablob's added to a synthetic muetablob with itself a list of sub-GLMweights'
  } else if ( ! is.null(zero_truncated <- family$zero_truncated)) {
    eta <- family$linkfun(mu_T, mu_truncated=zero_truncated) ## back to _U to add ranefs to the linear predictor.
  } else eta <- family$linkfun(mu_T)
  return(eta)
}

.subrange <- function(cumul, it) {cumul[it]+ seq_len(cumul[it+1L]-cumul[it])}

.fv_linkinv <- function(eta, # vector, not matrix : simulate(., nsim) -> .fv_linkinv(eta[, it], ...)
                        family, families=NULL, cum_nobs= attr(families,"cum_nobs")) {
  if (! is.null(families)) {
    fv <- p0 <- mu_U <- vector("list", length(cum_nobs)-1L)
    for (mv_it in seq_along(fv)) {
      resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
      if (length(resp_range)) fv[[mv_it]] <- .fv_linkinv(eta[resp_range], family=families[[mv_it]])
      # simulate(<mv>) with a Tnegbin shows the need to keep the attribute (-> .r_resid_var). 
      #     There is no code to compute it from mu_T without some attribute.
      p0[mv_it] <- list(attr(fv[[mv_it]], "p0"))
      mu_U[mv_it] <- list(attr(fv[[mv_it]], "mu_U"))
    }
    return(structure(unlist(fv, recursive = FALSE, use.names = TRUE), p0=p0, mu_U=mu_U))
  } else if ( ! is.null(zero_truncated <- family$zero_truncated)) {
    fv <- family$linkinv(eta,mu_truncated=zero_truncated)
  } else fv <- family$linkinv(eta) ## ! freqs for binomial, counts for poisson
  return(fv)
}

.calc_ZACw <- function(newZACpplist, augm_w_h_coeffs, map_rd_mv, newinold=NULL, cum_nobs) {
  if ( ! is.null(map_rd_mv)) {
    ZACwS <- vector("list", length(map_rd_mv))
    # newZACpplist should only contain the new ranefs
    for (mv_it in seq_along(map_rd_mv)) {
      # which(map_rd_mv[[mv_it]] %in% newinold) intersect two maps to old indices, and returns their position in the submodel
      # eg if newinold =2,3 and map_rd_mv ={{1,2},{3}} ....
      # pos_in_old <- which(map_rd_mv[[mv_it]] %in% newinold) # ... for 1st submodel, which= 2
      # newinold_it <- which(newinold %in% pos_in_old) # and this is the first ranef in newinold hence in the lists
      # which is simply...
      newinold_it <- which(newinold %in% map_rd_mv[[mv_it]])
      if (length(newinold_it)) {
        resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
        ZACwS[[mv_it]] <- .calc_ZACw(newZACpplist[newinold_it], 
                                     augm_w_h_coeffs[newinold_it], map_rd_mv=NULL)[resp_range]
      } else ZACwS[[mv_it]] <- numeric(cum_nobs[mv_it+1L]-cum_nobs[mv_it])
      # => final [resp_range] sacrificing mv performance to univariate one
    }
    return(unlist(ZACwS))
  }
  newnrand <- length(newZACpplist)
  if (newnrand>1L) {
    ZACw <- vector("list", length=newnrand )
    for (new_rd in seq_len(newnrand)) ZACw[[new_rd]] <- drop(newZACpplist[[new_rd]] %*% augm_w_h_coeffs[[new_rd]])
    ZACw <- do.call(cbind,ZACw)
    ZACw <- rowSums(ZACw)
  } else ZACw <- drop(newZACpplist[[1]] %*% augm_w_h_coeffs[[1]])
  return(ZACw)
}

.mvize <- function(fv, cum_nobs) {
  if (is.null(cum_nobs)) return(fv)
  # else:
  len <- length(cum_nobs)-1L
  fvs <- vector("list", len)  
  for (mv_it in seq_len(len)) {
    resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
    fvs[[mv_it]] <- fv[resp_range]
  }
  structure(fv, mv=fvs)
}

.point_predict <- function(object, newdata, new_X_ZACblob, variances, re.form, type, 
                           eta_fix=new_X_ZACblob$eta_fix ## may be NULL. addition of random-effect terms in the function
                           ) {
  if (.noRanef(re.form)) {
    if (type=="link") { # in those cases, we return eta but no useful $fv 
      fv <- NULL 
    } else {
      fv <- .fv_linkinv(eta_fix, object$family, object$families) ## ! freqs for binomial, counts for poisson
      fv <- .mvize(fv=fv, cum_nobs=new_X_ZACblob$cum_nobs)
    }
    return(list(fv=fv,eta=eta_fix))
  } else if ( is.null(newdata) && ! inherits(re.form,"formula")) {
    if (type=="link") {
      return(list(eta=object$eta)) 
    } else return(list(fv=object$fv)) ## eta will be reconstructed from fv on request
  } else { ## 
    newnrand <- length(new_X_ZACblob$subZAlist)
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
      augm_w_h_coeffs <- vector("list", newnrand)
      for (new_rd in seq_len(newnrand)) {
        augm_w_h_coeffs[[new_rd]] <- .match_old_new_levels(new_rd, new_X_ZACblob,
                                                           old_cum_n_u_h=attr(object$lambda,"cum_n_u_h"), 
                                                           w_h_coeffs=w_h_coeffs, 
                                                           rand.families=object$rand.families) # ~ solve(t(object$strucList[[1]]), object$v_h)
      }
      # Here we use a single rd-list of 'design' matrices new_X_ZACblob$newZACpplist for all submodels, 
      # which is thus redundant whenever there are redundant lines among different locdataS[[mv_it]]
      # .calc_ZACw() handles this by subsetting by resp_range's
      # The alternative would be to compute one new_X_ZACblob per submodel and call .point_predict() on each submodel 
      # This looks like asking for troubles.
      ZACw <- .calc_ZACw(newZACpplist=new_X_ZACblob$newZACpplist, 
                         augm_w_h_coeffs, map_rd_mv=attr(object$ZAlist, "map_rd_mv"), 
                         newinold=new_X_ZACblob$newinold,
                         cum_nobs=new_X_ZACblob$cum_nobs)
    }
    eta <- eta_fix + ZACw ## (length(eta)) col vector from coeffs = length(eta) row vector...
    # done with eta
    if (type=="link") {
      fv <- NULL # but eta will be returned
    } else {
      fv <- .fv_linkinv(eta, object$family, object$families, cum_nobs=new_X_ZACblob$cum_nobs) ## ! freqs for binomial, counts for poisson
      fv <- .mvize(fv=fv, cum_nobs=new_X_ZACblob$cum_nobs)
    }
    return(list(fv=fv,eta=eta))
  }
} 

.get_phiform <- function(object, mv_it=NULL) {
  if (object$spaMM.version < "3.5.23") {
    object$resid.predictor
  } else if (is.null(mv_it)) {
    object$residModel$formula
  } else object$residModels[[mv_it]]$formula
  
}

.get_phifam <- function(object, mv_it=NULL) { 
  if (object$spaMM.version < "3.5.23") {
    object$resid.family
  } else if (is.null(mv_it)) {
    object$residModel$family
  } else object$residModels[[mv_it]]$family
}

.to_respScale_var <- function(respVar, ppblob, object) {
  dmudeta <- NULL
  if ( ! is.null(families <- object$families)) {
    nonid_linkS <- logical(length(families))
    for (mv_it in seq_along(families)) nonid_linkS[mv_it] <- (families[[mv_it]]$link!="identity")
    if (any(nonid_linkS)) {
      cum_nobs <- attr(families,"cum_nobs")
      dmudetaS <- vector("list", length(nonid_linkS))
      if (is.null(eta <- ppblob$eta)) eta <- .eta_linkfun(ppblob$fv, family = NULL, families=families) 
      for (mv_it in seq_along(families)) {
        resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
        if (nonid_linkS[mv_it] && length(resp_range)) {
          dmudetaS[[mv_it]] <- families[[mv_it]]$mu.eta(eta[resp_range])
        } else dmudetaS[[mv_it]] <- rep(1, (cum_nobs[mv_it + 1L] - cum_nobs[mv_it]))
      }
      dmudeta <- unlist(dmudetaS, recursive = FALSE, use.names = TRUE)
    }
  } else if (object$family$link!="identity") {
    if (is.null(eta <- ppblob$eta)) eta <- object$family$linkfun(ppblob$fv) 
    dmudeta <- object$family$mu.eta(eta)
  }
  if ( ! is.null(dmudeta)) {
    if (is.list(respVar)) {  # (as_tcrossfac_list: repres of respVar as matrix factorization
      for (it in seq_len(length(respVar))) respVar[[it]] <- .Dvec_times_m_Matrix(dmudeta, respVar[[it]])
    } else if (!is.null(dim(respVar))) {
      respVar <-  .Dvec_times_m_Matrix(dmudeta, respVar) # sweep(respVar, MARGIN = 1, dmudeta, `*`)
      rownames(respVar) <- colnames(respVar) # _FIXME_ # Question is what is the logical behaviour here between the following two:
      # Rcpp_sweepZ1W loses names, not equivalent to sweepn and .Dvec_times_matrix() put backs the colnames only (by design)
      # sweep() keep all names
      ## Further _FIXME_ .matrix_times_Dvec(respVar, dmudeta) -> sweep(., MARGIN = 2...) keeps all names
      ##  (sparse matrix algo manipulating @x directly will also keep names)
      # so the reciprocal Dvec functions have inconsistent effects here
      # Anyway, for a symmetric respVar matrix, it makes sense to keep the rownames
      # and names are checked in a test-predVar item.
      respVar <- .m_Matrix_times_Dvec(respVar, dmudeta) # sweep(respVar, MARGIN = 2, dmudeta, `*`)
    } else respVar <- respVar*(dmudeta^2)
  }
  respVar
}

.warn_pw <- local({
  prior_weights_residvar_warned <- FALSE
  function(object) {
    pw <- object$prior.weights
    if (inherits(pw, "list")) {
      for (mv_it in seq_along(pw)) {
        if ((! prior_weights_residvar_warned ) && ! (attr(pw[[mv_it]],"unique") && pw[[mv_it]][1L]==1L)) {
          warning("Prior weights are not taken in account in residVar computation.")
          prior_weights_residvar_warned <<- TRUE
        }
      }
    } else if ((! prior_weights_residvar_warned) && ! (attr(pw,"unique") && pw[1L]==1L)) {
      warning("Prior weights are not taken in account in residVar computation.")
      prior_weights_residvar_warned <<- TRUE
    }
  }
})
 

.add_residVar <- function(object, resu, fv, locdata, respVar, variances) {
  .warn_pw(object) # residVar() wraps .get_phiW() that can handle prior weights 
  # .get_phiW() must be wrapped bc it is not OK for all objects, It can handle new data, but residVar() cannot =>
  # need other wrapper or new argument
  attr(resu,"residVar") <- .calcResidVar(object,newdata=locdata, fv=fv) 
  if (inherits(respVar,"matrix")) {
    nc <- ncol(respVar)
    diagPos <- seq.int(1L,nc^2,nc+1L)
    respVar[diagPos] <- respVar[diagPos] + attr(resu,"residVar")
  } else if ( is.list(respVar)) { # (as_tcrossfac_list)
    respVar[["residVar"]] <- Diagonal(x=sqrt(attr(resu,"residVar")))
  } else respVar <- respVar + attr(resu,"residVar")
  if (variances$respVar) attr(resu,"respVar") <- respVar
  return(resu)
}


.predict_body <- function(object, newdata, re.form, type,
                          variances, binding, intervals, level, blockSize, control, showpbar) {
  delayedAssign("invCov_oldLv_oldLv_list", .get_invColdoldList(object, control=control))
  if ( ! is.null(object$vec_nobs)) {
    new_X_ZACblob <- .calc_new_X_ZAC_mv(object=object, newdata=newdata, re.form = re.form,
                                     variances=variances, invCov_oldLv_oldLv_list=invCov_oldLv_oldLv_list)
  } else new_X_ZACblob <- .calc_new_X_ZAC(object=object, newdata=newdata, re.form = re.form,
                                   variances=variances, invCov_oldLv_oldLv_list=invCov_oldLv_oldLv_list) ## (fixme) still some unnecessary computation for predict(object)
  #
  ## (1) computes fv (2) compute predVar
  ##### fv
  ppblob <- .point_predict(object, newdata, new_X_ZACblob, variances, re.form, type) 
  if (type=="link" ) {
    resu <- ppblob$eta
  } else resu <- ppblob$fv ## mu(_T) by default
  # default binding appears to be FALSE ; map functions and blackbox::rbb use binding=<name>; binding=NA is used by simulate()
  if ( ! is.na(binding)){ ## suitable for objective function of optim() etc                ## ambiguous comment (about NA or !NA ?)    
    resu <- structure(as.matrix(resu),mu_U=attr(resu,"mu_U"),p0=attr(resu,"p0"),
                      mv=attr(resu,"mv") ## matrix ! maybe more suitable than data frame as objective function
    ) 
  }
  # if (identical(object$family$zero_truncated,TRUE)) {
  #   attr(resu,"p0") <- attr(ppblob$fv,"p0")
  #   attr(resu,"mu_U") <- attr(ppblob$fv,"mu_U")
  # }
  #\item 'make.names(.,unique= TRUE)' is applied to the input data, so the row names may be modified, if the input data did not contain unique, syntactically valid row names as defined by 'help(make.names)'.
  #rownames(resu) <- make.names(rownames(resu),unique = TRUE)
  locdata <- new_X_ZACblob$locdata 
  if ( ! is.logical(binding) ) { ## expecting a string
    if (inherits(locdata,"list")) { # mv case
      stop("'binding' operation not yet defined for multivariate-response models. Ask the maintainer.") 
      # wait for real-life example. (___FIXME___) and beware of code for intervals below
    } else {
      binding <- .makenewname(base=binding,varnames=colnames(locdata)) 
      resu <- structure(data.frame(resu), mu_U=attr(resu,"mu_U"),p0=attr(resu,"p0"))
      colnames(resu) <- binding
      resu <- cbind(locdata,resu) 
      attr(resu,"fittedName") <- binding
    }
  } else { ## alternative expecting binding= FALSE (but also handling binding = NA by doing nothing)
    if (! is.na(binding))  attr(resu,"frame") <- locdata 
  }
  if (inherits(locdata,"list")) attr(resu,"respnames") <- .get_from_terms_info(object=object, which="respnames")
  ##### (2) predVar
  if (variances$naive) { # the 'naive' estimate in BH98 (nu_i, without beta)
    naive <- .calc_Var_given_fixef(object, new_X_ZACblob=new_X_ZACblob, covMatrix=variances$cov, fix_X_ZAC.object=NULL)
    if (variances$cov) {
      rownames(naive) <- colnames(naive) <- rownames(locdata)
    } else names(naive) <- rownames(locdata)
    attr(resu,"naive") <- naive
  }
  #
  newnrand <- length(new_X_ZACblob$newZAlist) ## may be reduced if non trivial re.form
  newX.pv <- new_X_ZACblob$newX.pv
  #
  beta_w_cov_needed <- (# inherits(re.form,"formula") || 
                          (variances$linPred && newnrand)
                        )
  if (beta_w_cov_needed) loc_tcrossfac_beta_w_cov <- .get_tcrossfac_beta_w_cov(object) 
  if (beta_w_cov_needed && inherits(re.form,"formula")) {
    # identifies and selects columns for the [retained ranefs, which are given by newinold 
    re_form_col_indices <- .re_form_col_indices(new_X_ZACblob$newinold, cum_n_u_h=attr(object$lambda,"cum_n_u_h"), Xi_cols=attr(object$ZAlist, "Xi_cols"))
    Xncol <- ncol(object$X.pv)
    subrange <- c(seq_len(Xncol),re_form_col_indices$subrange + Xncol)
    loc_tcrossfac_beta_w_cov <- loc_tcrossfac_beta_w_cov[subrange,]
  } else re_form_col_indices <- NULL
  #
  if(variances$linPred) {
    if (newnrand) {
      if ( is.null(newdata) && ! inherits(re.form,"formula")) {
        newZAlist <- new_X_ZACblob$subZAlist ## (subset of) the original object$ZAlist 
      } else { ## 
        newZAlist <- new_X_ZACblob$newZAlist ## a new ZAlist built using .calc_Zlist(new formula...) 
      }
    }
    if (newnrand) { # if there are ranefs in 'new' formula (which may be orig formula)
      if (variances$cancel_X.pv) newX.pv[] <- 0 
      respVar <- .calcPredVar(X.pv=newX.pv,tcrossfac_beta_w_cov=loc_tcrossfac_beta_w_cov,
                              covMatrix=variances$cov,object=object,
                              newinold=new_X_ZACblob$newinold,blockSize=blockSize, 
                              newZAlist=newZAlist, re_form_col_indices=re_form_col_indices,
                              invCov_oldLv_oldLv_list=if (is.null(newdata)) {
                                NULL} else {invCov_oldLv_oldLv_list}, # promise for newdata
                              logdispObject=if (variances$disp) {
                                .get_logdispObject(object)} else {NULL}, # promise for variances$disp
                              cov_newLv_oldv_list=new_X_ZACblob$cov_newLv_oldv_list,
                              ## list for Cnewnew, which enters in  newZA %*% Cnewnew %*% tnewZA, hence should not represent newZA itself 
                              cov_newLv_newLv_list=new_X_ZACblob$cov_newLv_newLv_list,
                              diag_cov_newLv_newLv_list=new_X_ZACblob$diag_cov_newLv_newLv_list,
                              as_tcrossfac_list=variances$as_tcrossfac_list,
                              showpbar=showpbar) 
      if ( ! is.list(respVar)) {
        if (inherits(locdata,"list")) { # mv
          respVnames <- .unlist(lapply(locdata,rownames))
        } else respVnames <- rownames(locdata)
        if (variances$cov) {
          respVar <- as.matrix(respVar) ## matrix, not Matrix (assumed below)
          rownames(respVar) <- colnames(respVar) <- respVnames
        } else {
          names(respVar) <- respVnames
        }
      } 
    } else { # no ranefs
      if (inherits(locdata,"list")) { # mv
        respVnrow <- sum(.unlist(lapply(locdata,nrow)))
      } else respVnrow <- nrow(locdata)
      if (variances$cov) {
        respVar <- matrix(0,nrow=respVnrow,ncol=respVnrow)
      } else respVar <- rep(0,respVnrow)
    }
  } else if (any(unlist(variances))) {
    if (inherits(locdata,"list")) { # mv
      respVnrow <- sum(.unlist(lapply(locdata,nrow)))
    } else respVnrow <- nrow(locdata)
    respVar <- rep(0,respVnrow)
  } else respVar <- NULL 
  if ( variances$fixefVar || (newnrand==0L && variances$linPred) ) {
    beta_cov <- .get_beta_cov_any_version(object)
    if (! is.null(beta_cov)) {
      fixefcov <- newX.pv %*% beta_cov %*% t(newX.pv)
      if (variances$cov) {
        attr(resu,"fixefVar") <- fixefcov 
      } else attr(resu,"fixefVar") <- diag(fixefcov) # if beta_cov is a Matrix, fixefcov is too, and diag loses names -> detected by tests.
      if (newnrand==0L) { ## otherwise there is already such a term in predVar
        respVar <- respVar + attr(resu,"fixefVar") 
      }
    }
  }
  if ( is.list(respVar)) {
    attr(resu,"predVar") <- respVar 
  } else {
    if (variances$cov) respVar <- (respVar+t(respVar))/2 ## if numerically asym, rand_eta <- mvrnorm(.,Sigma=attr(point_pred_eta,"predVar")) fails
    attr(resu,"predVar") <- respVar ## vector or matrix
  }
  if ( ! is.null(respVar)) respVar <- .to_respScale_var(respVar, ppblob, object)
  if (variances$residVar) resu <- .add_residVar(object, resu, fv=ppblob$fv, locdata, respVar, variances)
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
  if(length(intervals)) { # for e.g. intervals="predVar"
    intervalresu <- NULL
    # need eta as vector in all cases (if eta is data.frame execution stops on cbind(<NULL>, interval))
    if (type=="link") {
      if (is.character(binding)) { # then resu is data.frame
        eta <- resu[[binding]] # data.frame -> vector
      } else {eta <- resu; dim(eta) <- NULL} # matrix -> vector
    } else {
      if (is.character(binding)) { # then resu is data.frame
        eta <- .eta_linkfun(resu[[binding]] , object$family, object$families) # with possible mu attribute
      } else eta <- .eta_linkfun(resu , object$family, object$families)  # with possible mu attribute
      ## mv: 'resu' ('center' of intervals) still a single column 
      eta <- eta[seq(nrow(resu))] ## [] with explicit indices, needed here to drop possible mu attribute otherwise linkinv(eta+/-sd) uses it! 
      #                           #   and so -> vector in all case  
    } 
    # 
    for (st in intervals) {
      varcomp <- attr(resu,st)
      if (is.null(varcomp)) warning(paste("Prediction variance component",st,"requested but not available: check input."))
      if (is.matrix(varcomp)) varcomp <- diag(varcomp)
      pv <- 1-(1-level)/2
      ## special case for simple LM
      if (length(object$rand.families)==0L && # not mixed
          object$family$family=="gaussian" && ## __FIXME__ so calling this on a mv-GLM will bug.
          .DEPARSE(.get_phiform(object))=="~1" # not heteroscedastic
      ) { 
        nobs <- length(object$y)
        resdf <- nobs - ncol(object$X.pv) ## don't use fixef here, that contains bot NAs and argument etaFix$beta! 
        is_REML <- ( .REMLmess(object,return_message=FALSE))
        if ( ! is_REML) {
          vart <- varcomp*nobs/resdf
        } else vart <- varcomp
        ## FR->FR (use a type attribute for fixef ?)
        sd <- stats::qt(pv,df=resdf)*sqrt(vart) # student for LM (not LMM)
      } else {
        sd <- qnorm(pv)*sqrt(varcomp) 
        # for GLMMs a better bound might be obtained by convolution of residual error distrib and of predVar  ~gaussian on linear scale
      }
      if (variances$respVar) { # sd is already on response scale
        mu <- .fv_linkinv(eta=eta, family=object$family, families=object$families)
        interval <- cbind(mu-sd,mu+sd)
      } else interval <- cbind(.fv_linkinv(eta=eta-sd, family=object$family, families=object$families),
                               .fv_linkinv(eta=eta+sd, family=object$family, families=object$families))
      colnames(interval) <- paste(st,c(signif(1-pv,4),signif(pv,4)),sep="_")
      intervalresu <- cbind(intervalresu,interval) # recursive cbind over types of intervals
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
                          variances=list(), binding = FALSE, intervals = NULL,
                          level = 0.95, blockSize = 2000L, type = "response", 
                          verbose=c(showpbar=eval(spaMM.getOption("barstyle"))), 
                          control=list(), ...) { ## but not new Y
  if (is.null(object$envir)) object$envir <- list2env(list(), ## back-compatibility fix for old objects
                                                     parent=environment(HLfit_body))
  ## the final components returned as attributes have names ...Var, other terms should be named differently
  #
  if ( ! is.null(intervals) && ! inherits(intervals,"character")) stop("'intervals' arguments should inherit from class 'character'.")
  checkIntervals <- (substr(x=intervals, nchar(intervals)-2, nchar(intervals))=="Var")
  if (any(!checkIntervals)) warning("Element(s)",intervals[!checkIntervals],"are suspect, not ending in 'Var'.")
  # possible elements in return value: fixefVar, predVar, residVar, respVar
  variances[intervals] <- TRUE 
  variances <- .process_variances(variances, object)
  # if (type=="rand_eta") { # called by simulate(., type="[predVar" || "(ranef|response)"] ))
  #   rand_eta <- list(...)$rand_eta # this may not even be useful except to avoid a CRAN check NOTE
  #   type <- "link"
  #   variances[names(rand_eta)] <- rand_eta
  # }
  nrX <-  NROW(newdata)
  if (!is.null(re.form) && inherits(re.form,"formula")) re.form <- .preprocess_formula(re.form)
  showpbar <- verbose[["showpbar"]]
  ############################## if (nrX>0L) newdata <- droplevels(newdata) FIXME perhaps here ? 
  if ( (! variances$cov) && nrX > blockSize) {
    ### this part of code is tested by the test-predVar code on Loaloa data
    # et par test geostat dans probitgem (iterateSEMSmooth -> .sampleNextPars -> .spaMM_rhullByEI)
    ## newdata <- droplevels(newdata) ## potential gain of time for droplevels(newdata_slice)
    slices <- unique(c(seq(0L,nrX,blockSize),nrX))
    nslices <- length(slices)-1L
    res <- vector("list",nslices)
    progrbar_setup <- .set_progrbar(style = showpbar, char="s") # FIXME could implement parallel computation
    for (it in seq_len(nslices)) {
      slice <- (slices[it]+1L):slices[it+1L]
      newdata_slice <- newdata[slice,,drop=FALSE]
      ## newdata_slice <- droplevels(newdata_slice) 
      res[[it]] <- .predict_body(object=object, newdata=newdata_slice, re.form = re.form, variances=variances, 
                                 binding=binding, type=type, intervals=intervals, level=level, blockSize=blockSize, ## blockSize should not be useful *here*
                                 control=control, showpbar=showpbar)
      if (showpbar) progrbar_setup$progress(slices[it+1L]/nrX)  ## update progress bar
    }
    if (showpbar) close(progrbar_setup$pb)
    if (is.character(binding)) {
      res <- do.call(rbind, res) # rbind data frames
    } else res <- .unlist_with_attributes(res) # vector for binding=NA, 1-col matrix for binding=FALSE
  } else if (type=="marginal") {
    res <- .predict_marg(object=object, newdata=newdata, re.form = re.form, control=control)
  } else res <- .predict_body(object=object, newdata=newdata, re.form = re.form,
                variances=variances, binding=binding, type=type,
                intervals=intervals, level=level, blockSize=blockSize, ## but blockSize could be useful *here* if newdata was NULL
                control=control, showpbar=showpbar)
  return(res)
}

print.vcov.HLfit <- function(x, expanded=FALSE, ...) {
  a <- attributes(x)
  print.default(x)
  cat("with additional attribute(s):")
  std.attr <- c("names","dim","dimnames","class") ## attributes not to be shown
  nam <- names(a)
  if (expanded) { # shows structure of attributes as in utils:::str.default. This was useful for beta_v_cov attribute, but now ?
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


`[.predictions` <- local({
  mv_warned <- FALSE # we need special handling for all attributes that become lists in mv code.
                    # => only 'frame' attribute ?
  function (x, i, j, 
            drop = TRUE ## by default, this function will return scalar/vector 
  ) {
    class(x) <- "matrix" ## avoids recursive call to `[.predictions` 
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
      if ( ! is.null(frame)) {
        if (inherits(frame,"list")) {
          # base unique() -> unique.matrix() -> do.call("[", c(list(x), args, list(drop = FALSE)))
          # -> bug on mv output where the frame attribute is a list of frames.
          # 'frame' might be needed only in the context of binding predictions, so best code not yet clear
          if ( ! mv_warned) {
            warning("'frame' attribute of predictions not yet handled for predictions from multivariate fits.\n Contact the maintaner if you encounter problems.")
            mv_warned <<- TRUE
          }
        } else frame <- frame[i,] ## dataframe => nodrop
      }
      residVar <- attr(x, "residVar")
      if ( ! is.null(frame)) residVar <- residVar[i,drop=FALSE]
      respVar <- attr(x, "respVar")
      if ( ! is.null(respVar)) {
        if (is.null(dim(respVar))) {
          respVar <- respVar[i]
        } else respVar <- respVar[i,i,drop=FALSE]
      }
      class(resu) <- c("predictions","matrix")
      structure(resu,fixefVar=fixefVar,predVar=predVar,residVar=residVar,frame=frame,fittedName=attr(x, "fittedName"))
    } else return(resu)
  } # Use unlist() to remove attributes from the return value
}) 

print.predictions <- function (x, expanded=FALSE, ...) {
  asvec <- as.vector(x) ## important to remove names and keep them separately
  rnames <- rownames(x)
  if (is.null(rnames)) rnames <- rownames(attr(x,"frame"))
  toolong <- nchar(rnames)>9
  rnames[toolong] <- paste0(substr(rnames[toolong],0,8),".")
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

# usage: 
# as.matrix(<1-col matrix inheriting from class "predictions" resulting from predict(, binding=FALSE)>, ...)
# or 
# as.matrix.predictions(<vector resulting from predict(, binding=NA)>, ...)
# as.matrix.predictions <- function(x, colnames.=attr(x, "respnames"), ...) {
#   #mv <- attr(x,"mv")
#   #vec_nobs <- sapply(mv, nrow)
#   resu <- matrix(x,ncol=length(attr(x,"mv"))) # or do.call(rbind, attr(x,"mv")) ?
#   colnames(resu) <- colnames.
#   resu
# } 
