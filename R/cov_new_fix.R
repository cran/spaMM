.calc_corr_from_dist <- function(distmat, object, corr.model) {
  if (corr.model=="AR1") {
    args <-object$ranFix[which(names(object$ranFix) %in% c("ARphi"))]
    corr_mat <- args$ARphi^distmat 
  } else {
    args <-object$ranFix[which(names(object$ranFix) %in% c("nu","Nugget"))] ## so that rho=1 in MaternCorr
    corr_mat <- do.call(MaternCorr,args=c(args,list(d=distmat)))  
  }
  return(corr_mat) 
}

# public wrapper for more transparent workflow
preprocess_fix_corr <- function(object, fixdata, re.form = NULL,
                                variances=list(residVar=FALSE)) {
  .calc_new_X_ZAC(object=object, newdata=fixdata, re.form = re.form,
                             variances=variances) 
}
###############################################

.make_corr_list <- function(strucList, newZAlist) {
  if ( ! is.null(newZAlist)) {
    locstrucList <- strucList
    for (it in seq_len(length(locstrucList))) {
      if (! is.null(locstrucList[[it]])) {
        if (attr(locstrucList[[it]],"corr.model")=="random-coef" &&
            ncol(newZAlist[[it]]) != nrow(strucList[[it]]) ## new groups for the random-coef term
        ) {
          locstrucList[[it]] <- .makelong(attr(locstrucList[[it]],"latentL_blob")$u,longsize = ncol(newZAlist[[it]]))
          rownames(locstrucList[[it]]) <- colnames(newZAlist[[it]]) 
          # .tcrossprod gets names from rownames
          ## names required for a predict(,newdata) on a random-coef model
        }
      }
    }
    corr_list <- lapply(locstrucList, .tcrossprod,y=NULL)
  } else corr_list <- lapply(strucList, .tcrossprod,y=NULL)
  for (it in seq_len(length(corr_list))) {
    if ( ! is.null(corr_list[[it]])) {
      attr(corr_list[[it]],"ranefs") <- attr(strucList[[it]],"ranefs")
      attr(corr_list[[it]],"corr.model") <- attr(strucList[[it]],"corr.model")
    }
  }
  return(corr_list)
}

.calc_sub_diagmat_cov_newLv_oldv <- function(oldZA, newZA,namesTerms) {
  ## provides a Diagonal cov matrix for un-autocorrelated ranefs. Ignores unrepresented old levels, 
  ##   while covs for autocorrelated ranefs keep unrepresented old levels because it takes time to remove them.
  oldlevels <- colnames(oldZA)
  newlevels <- colnames(newZA)
  if (identical(oldlevels,newlevels)) {
    newoldC <- Diagonal(n=length(oldlevels)) ## replaces old identityMatrix
  } else {
    oldornew <- unique(c(oldlevels,newlevels))
    numnamesTerms <- length(namesTerms) ## 2 for random-coef
    newoldC <- diag(length(oldornew)*numnamesTerms)
    colnames(newoldC) <- rownames(newoldC) <- rep(oldornew, numnamesTerms)
    newoldC <- newoldC[newlevels,oldlevels,drop=FALSE]
  }
  attr(newoldC,"isEachNewLevelInOld") <- newlevels %in% oldlevels  ## but this attr is unevaluated (-> NULL) for spatial models 
  return(newoldC)
}


.calc_new_X_ZAC <- function(object, newdata=NULL, re.form = NULL,
                         variances=list(residVar=FALSE)) {
  locform <- attr(object$predictor,"oriFormula")
  ## possible change of random effect terms
  if ( .noRanef(re.form)) { ## i.e. if re.form implies that there is no random effect
    locform <- .nobarsMM(locform)  ## 
  } else if (inherits(re.form,"formula")) { ## ie there is a nontrivial re.form
    lhs <- .DEPARSE(locform[[2L]]) ## response
    fixterms <- .DEPARSE(.nobarsMM(locform)[[3L]]) ## fixed effect terms
    ranterms <- .parseBars(re.form)  ## random effects  ## expand does not seem to be important here
    locform <- as.formula(paste(lhs,"~", paste(fixterms,"+",paste(ranterms,collapse="+")) )) 
  } ## else keep oriFormula
  # checking variables in the data
  if (length(locform)==3L) locform <- locform[-2] ## removes RHS for checking  vars on RHS etc
  allvars <- all.vars(locform) 
  if (variances$residVar) allvars <- unique(c(allvars,all.vars(attr(object$resid.predictor,"oriFormula")))) ## but the 'newdata' may not contain the resid.predictor vars. 
  if (is.vector(newdata)) { ## ## less well controlled case, but useful for maximization
    newdata <- data.frame(matrix(newdata,nrow=1))
    if (length(allvars)==ncol(newdata)) {
      names(newdata) <- allvars
    } else {
      stop(paste("(!) newdata has incorrect length. It should match the following variables:\n",paste(allvars,collapse=" ")))
    }
  } 
  ## it is important that newdata remains NULL if it was so initially because it is tested below. Hence copy in locdata
  if (is.null(newdata)) {
    locdata <- object$data[,allvars,drop=FALSE]
  } else {
    if( is.matrix(newdata) ) newdata <- as.data.frame(newdata)  
    # so that matrix 'newdata' arguments can be used as in some other predict methods.
    locdata <- newdata[,allvars,drop=FALSE] ## allvars checks only RHS variables
    checkNAs <- apply(locdata,1L,anyNA)
    if (any(checkNAs)) {
      message("NA's in required variables from 'newdata'. Prediction not always possible.")
      locdata <- locdata[!checkNAs, ]
    }
  }
  need_new_design <- ( ( ! is.null(newdata) ) || ! is.null(re.form)) ## newdata or new model
  if (need_new_design) {
    allFrames <- .HLframes(formula=locform,data=locdata,fitobject=object) ## also used for predVar computations
    ## preparation for fixed effects
    newX.pv <- allFrames$X ## contains columns for the offset and columns for the other variables
    # newX.pv must intersect non-NA elements of fixef; see comment and code in newetaFix
    validnames <- intersect(colnames(object$X.pv),colnames(newX.pv)) ## we don't want the etaFix cols (detected by bboptim)
    if (length(validnames)==0L) validnames <- c() ## without this, validnames could be character(0) and [,validnames,drop=FALSE] fails.
    RESU <- list(locdata=locdata,newX.pv=newX.pv[,validnames,drop=FALSE]) 
    RESU$etaFix <- .newetaFix(object,allFrames,validnames=validnames)  ## new fixed effects (or [part of] new linear predictor if re.form)      
  } else RESU <- list(locdata=locdata,newX.pv=object$X.pv) 
  ## preparation for random effects
  if ( ! is.null(newdata)) { 
    if (object$spaMM.version < "2.2.116") {
      bad_corr_types <- intersect(attr(attr(object$ZAlist,"ranefs"),"type"),c("adjacency","corrMatrix"))
    } else bad_corr_types <- intersect(attr(object$ZAlist,"exp_ranef_types"),c("adjacency","corrMatrix")) ## AR1 and Matern are OK
    if (length(bad_corr_types)) { 
      stop(paste("Prediction in 'newdata' not implemented or not possible for models including ", 
                 paste(bad_corr_types, collapse="and"), "term(s)."))
    } 
  }
  # newZAlist and subZAlist appear to have distinct usages since they are created under different conditions.
  # subZAlist is a subset of the old ZA, newZAlist contains new ZA
  # calling .make_corr_list(object$strucList,...) is always OK bc the fist argument may be a superset of the required list
  # all matching in .make_corr_list is through the ranef attributes.
  #
  ## matching ranef terms of re.form
  if (.noRanef(re.form)) {
    nrand <- 0L
  } else {
    if (object$spaMM.version < "2.2.116") {
      ori_exp_ranef_strings <- attr(object$ZAlist,"ranefs") 
      ori_exp_ranef_types <- attr(ori_exp_ranef_strings,"type") 
    } else {
      ori_exp_ranef_strings <- attr(object$ZAlist,"exp_ranef_strings")
      ori_exp_ranef_types <- attr(object$ZAlist,"exp_ranef_types") 
    }
    RESU$spatial_old_rd <- which(ori_exp_ranef_types != "(.|.)")   
    if (inherits(re.form,"formula")) {
      new_exp_ranef_strings <- .parseBars(locform,expand=TRUE) ## note importance for the computation of newZAlist below
      nrand <- length(new_exp_ranef_strings)
      newinold <- sapply(lapply(new_exp_ranef_strings, `==`, y= ori_exp_ranef_strings), which) ## e.g 1 4 5
      RESU$subZAlist <- object$ZAlist[newinold] ## and reordered
      sublambda <- object$lambda[newinold]
    } else {
      new_exp_ranef_strings <- ori_exp_ranef_strings
      nrand <- length(new_exp_ranef_strings)
      newinold <- seq(nrand) ## keep all ranefs
      RESU$subZAlist <- object$ZAlist
      sublambda <- object$lambda
    }    
    RESU$newinold <- newinold
  }
  #
  if (nrand) {
    strucList <- object$strucList
    if (object$spaMM.version<"1.11.57") stop("This fit object was created with spaMM version<1.11.57, and is no longer supported.\n Please recompute it.")
    which_mats <- list(no= need_new_design, 
                       ## cov_newLv_newLv_list used in .calc_Evar() whenever newdata, but elements may remain NULL if $cov not requested
                       ## However, for ranCoefs, we need Xi_cols rows for each response's predVar. (FIXME) we store the full matrix.
                       nn= ( attr(object$strucList,"isRandomSlope")[newinold] | identical(variances$cov,TRUE))) 
    if (need_new_design) {
      ## with newdata we need Evar and then we need nn... if newdata=ori data the Evar (computed with the proper nn) should be 0
      newZAlist <- .spMMFactorList(locform, allFrames$mf, 0L, drop=TRUE,sparse_precision=FALSE, type="seq_len") 
      ## must be ordered as parseBars result for the next line to be correct.
      attr(newZAlist,"exp_ranef_strings") <- new_exp_ranef_strings ## required pour .compute_ZAXlist to match the ranefs of LMatrix
    } else {
      newZAlist <- object$ZAlist
    }
    ## AT this point both newZAlist and subZAlist may be reduced to 'newnrand' elements relative to ori object$ZAlist.
    ## .match_old_new_levels will use it running over newnrand values
    ## The cov_ lists are reduced too. newinold should be used to construct them
    ## newZACpplist is reduced.
    blob <- .make_new_corr_lists(object=object,locdata=locdata, which_mats=which_mats, 
                                 newZAlist=newZAlist, newinold=newinold)
    # cov_newLv_oldv_list is always needed for cbind(X.pv,newZAC [which may be ori ZAC]); should ~corr_list when newdata=ori data
    cov_newLv_oldv_list <- blob$cov_newLv_oldv_list 
    cov_newLv_newLv_list <- blob$cov_newLv_newLv_list ## may be NULL
    # cov_newLv_oldv_list still contains NULL's, and we need something complete to compute newZAC without '%*% NULL'
    namesTerms <- attr(newZAlist,"namesTerms")
    for (new_rd in seq_along(cov_newLv_oldv_list)) { ## seems correct as everything in rhs is reduced
      if (is.null(cov_newLv_oldv_list[[new_rd]])) {
        cov_newLv_oldv_list[[new_rd]] <- .calc_sub_diagmat_cov_newLv_oldv(oldZA=RESU$subZAlist[[new_rd]], newZA=newZAlist[[new_rd]],
                                                                                         namesTerms=namesTerms[[new_rd]]) ## non-trivial value [I] if uuCnewold is NULL
        ## without a ranefs attribute, hence ignored by .compute_ZAXlist
        if (which_mats$nn[new_rd]) {
          newc <- Diagonal(n=ncol(newZAlist[[new_rd]]))
          rownames(newc) <- colnames(newc) <- colnames(newZAlist[[new_rd]])
          cov_newLv_newLv_list[[new_rd]] <- newc
        }
      }
    }
    if (any(which_mats$nn)) RESU$cov_newLv_newLv_list <- cov_newLv_newLv_list
    RESU$cov_newLv_oldv_list <- cov_newLv_oldv_list
    RESU$newZACpplist <- .compute_ZAXlist(XMatrix=cov_newLv_oldv_list, ZAlist=newZAlist) ## build from reduced list, returns a reduced list
    ## This $newZACpplist serves to compute new _point predictions_.
    #  .compute_ZAXlist affects elements of ZAlist that have a ranefs attribute. 
    #  It builds a design matrix to all oldv levels. It does not try to reduce levels. 
    #  But it use newZAlist <- .spMMFactorList(...) which may consider a reduced number of levels.
    #  The function .match_old_new_levels() will perform the match with 'w_h_coeffs' for point prediction.
    #  it assign values psi_M to new levels of ranefs.
    ## _Alternatively_ one may construct a newZACvar for _predVar_ 
    #  Here we have a slice mechanism (contrary for point pred) hence new ZA with different rows, and the columns of the 
    #  newZACvar constructed there must match those of beta_w_cov for  ZWZt_mat_or_diag( <cbind(newX,newZAC)> ,beta_w_cov)
    #  cov_newLv_oldv_list() provides the info for the expansion from the newZA cols to the oldZA cols.
    #  In that case one does not need to match levels. .calc_newZACvar() performs a simpler operation than .compute_ZAXlist.
    #############################
    RESU$spatial_term <- NaN
    if (object$spaMM.version < "2.2.20") {
      spatial.terms <- .findSpatial(locform) ## list of formula terms
      spatial_term <- spatial.terms[[1L]] ## one formula term, e.g Matern(1|...)
    } else if (object$spaMM.version < "2.2.118") {
      spatial_term <- object$spatial_term
    } else {
      spatial_terms <- attr(object$ZAlist,"exp_spatial_terms")
      if ( ! is.null(spatial_terms) && any( true_spatial_terms <- (! is.na(spatial_terms)))) {
        spatial_term <- (spatial_terms[true_spatial_terms])[[1L]] ## tempo befre spaMM 3.0
      } else spatial_term <- NULL
    }
    if ( ! is.null(spatial_term)) { 
      olduniqueGeo <- attr(object,"info.uniqueGeo")
      geonames <- colnames(olduniqueGeo)
      if (need_new_design) geonames <- intersect(geonames, colnames(locdata)) 
      RESU$newuniqueGeo <- locdata[,geonames,drop=FALSE]
    }
    # compatiblity code that should disappear:
    uuCnewnew <- NULL
    #
    RESU$newZAlist <- newZAlist
  } 
  return(RESU)
}

## get_predCov_var_fix: see example in predict.Rd (?get_predCov_var_fix), test in test-predVar 
get_predCov_var_fix <- function(object, newdata = NULL, fix_X_ZAC.object,fixdata, 
                                variances=list(disp=TRUE,residVar=FALSE),...) {
  newnrand <- length(fix_X_ZAC.object$newZAlist) 
  fixZACvar <- .calc_newZACvar(fix_X_ZAC.object$newZAlist,fix_X_ZAC.object$cov_newLv_oldv_list)
  new_X_ZACblob <- .calc_new_X_ZAC(object,newdata=newdata,variances=variances) ## called for a correlation block
  newZACvar <- .calc_newZACvar(new_X_ZACblob$newZAlist,new_X_ZACblob$cov_newLv_oldv_list)
  ## First component of predVar
  # covariance of expectation of Xbeta+Zb due to var of (hat(beta),hat(v)) using E[b] as function of hat(v)
  ## (X_n | C_no) %*% [ t(invL) %*% beta_v_cov[v.range,] %*% invL ] %*% t(X_f | C_fo)
  fix_X_ZAC <- cbind2(fix_X_ZAC.object$newX.pv, fixZACvar)
  new_X_ZAC <- cbind2(new_X_ZACblob$newX.pv, newZACvar)
  beta_w_cov <- .get_beta_w_cov(object)
  predVar <- new_X_ZAC[] %id*id% beta_w_cov[] %id*id% t(fix_X_ZAC)[]
  ## Second component of predVar:
  cov_newLv_fixLv_list <- .make_new_corr_lists(object,new_X_ZACblob$newuniqueGeo, 
                                               which_mats=list(no=TRUE,nn=rep(FALSE,newnrand)), 
                                               new_X_ZACblob$newZAlist, 
                                               newinold=fix_X_ZAC.object$newinold,
                                               fix_info=fix_X_ZAC.object)$cov_newLv_oldv_list
  if ( ! is.null(cov_newLv_fixLv_list) ) {
    # Evar: expect over distrib of (hat(beta),hat(v)) of [covvariance of Xbeta+Zb given (hat(beta),hat(v))]
    Evar <- .calc_Evar(newZAlist=new_X_ZACblob$newZAlist,newinold=fix_X_ZAC.object$newinold, 
                                    cov_newLv_oldv_list=new_X_ZACblob$cov_newLv_oldv_list, 
                                    lambda=object$lambda, invCov_oldLv_oldLv_list=.get_invColdoldList(object), 
                                    cov_newLv_fixLv_list=cov_newLv_fixLv_list, cov_fixLv_oldv_list=fix_X_ZAC.object$cov_newLv_oldv_list, 
                                    fixZAlist=fix_X_ZAC.object$newZAlist,covMatrix=TRUE)
    predVar <- predVar + Evar
  } 
  # If components for uncertainty in dispersion params were requested,
  #   logdispObject is not NULL
  # If some components ere computable, $$dwdlogdisp should not be NULL
  # Former approach (changed 08/2016) was to test logdispObject and then 
  #   for any 'problem'. But there may some 'problem' and still a valid logdispObject
  # => In new version, dwdlogdisp should be either NULL or a conforming matrix;
  #  'problems" should not be tested.
  if (variances$disp) logdispObject <- .get_logdispObject(object)
  if ( ! is.null(logdispObject$dwdlogdisp) ) {
    col_info <- attr(logdispObject$dwdlogdisp,"col_info") ## ranefs for which there is a col in dwdlogdisp
    if ( newnrand != col_info$nrand) { ## selection of blocks for re.form ranefs # fixme: this code is not (testthat-)checked
      submatrices <- .submatrices_for_disp_effect(logdispObject, col_info, fix_X_ZAC.object$newinold) 
    } else submatrices <- logdispObject
    newZACw <- newZACvar %*% submatrices$dwdlogdisp ## typically (nnew * n_u_h) %*% (n_u_h * 2) = nnew * 2 hence small 
    fixZACw <- fixZACvar %*% submatrices$dwdlogdisp ## typically (nnew * n_u_h) %*% (n_u_h * 2) = nnew * 2 hence small 
    disp_effect_on_newZACw <- newZACw %*% submatrices$logdisp_cov %*% t(fixZACw)      
    predVar <- predVar + disp_effect_on_newZACw
  }
  return(predVar)
}
