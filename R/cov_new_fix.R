.calc_corr_from_dist <- function(resu, object, spatial.model, which) {
  corr.model <- as.character(spatial.model[[1]])
  if (corr.model=="AR1") {
    args <-object$ranFix[which(names(object$ranFix) %in% c("ARphi"))]
    if (which$no) resu$uuCnewold <- args$ARphi^resu$uuCnewold  
    if (which$nn) resu$uuCnewnew <- args$ARphi^resu$uuCnewnew  
  } else {
    args <-object$ranFix[which(names(object$ranFix) %in% c("nu","Nugget"))] ## so that rho=1 in MaternCorr
    if (which$no) resu$uuCnewold <- do.call(MaternCorr,args=c(args,list(d=resu$uuCnewold)))  
    if (which$nn) resu$uuCnewnew <- do.call(MaternCorr,args=c(args,list(d=resu$uuCnewnew)))  
  }
  return(resu) 
}

.calc_corr_list <- function(uuCnewold, spatialOne, subZAlist, FL) {
  newZAlist <- FL$Design
  namesTerms <- FL$namesTerms
  nrand <- length(newZAlist)
  locfn <- function(it) {
    if (it %in% spatialOne) { ## avoids running the next algo which is slow on large matrices
      oldnewC <- t(uuCnewold)
    } else {
      oldlevels <- colnames(subZAlist[[it]])
      newlevels <- colnames(newZAlist[[it]])
      numnamesTerms <- length(namesTerms[[it]]) ## 2 for random-coef
      if (identical(oldlevels,newlevels)) {
        oldnewC <- Diagonal(length(oldlevels)*numnamesTerms) ## replaces old identityMatrix
      } else {
        oldornew <- unique(c(oldlevels,newlevels))
        oldnewC <- diag(length(oldornew)*numnamesTerms)
        colnames(oldnewC) <- rownames(oldnewC) <- rep(oldornew, numnamesTerms)
        oldnewC <- oldnewC[oldlevels,newlevels]
      }
      attr(oldnewC,"isEachNewLevelInOld") <- newlevels %in% oldlevels  ## but this attr is unevaluated (-> NULL) for spatial models 
    }
    return(oldnewC)
  }
  oldnewClist <- lapply(seq_len(nrand),locfn)
  return(oldnewClist)
}

## obsolete:
.calc_old_corr <- function(Lmatrix) {
  if(is.null(Lmatrix)) return(NULL)
  corr.model <- attr(Lmatrix,"corr.model")
  if (corr.model=="random-coef") {
    LforLv <- .calc_latentL(attr(Lmatrix,"cov.mat"))$u
    old_corr <- tcrossprodCpp(LforLv,NULL)
  } else old_corr <- tcrossprodCpp(Lmatrix,NULL)
  attr(old_corr,"ranefs") <- attr(Lmatrix,"ranefs")
  attr(old_corr,"corr.model") <- attr(Lmatrix,"corr.model")
  return(old_corr)
}

# public wrapper for more transparent workflow
preprocess_fix_corr <- function(object, fixdata, re.form = NULL,
                                variances=list(residVar=FALSE)) {
  .calc_new_X_ZAC(object=object, newdata=fixdata, re.form = re.form,
                             variances=variances, needNewNew=FALSE) 
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
          locstrucList[[it]] <- .makelong(attr(locstrucList[[it]],"Lcompact"),longsize = ncol(newZAlist[[it]]))
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


.calc_new_X_ZAC <- function(object, newdata=NULL, re.form = NULL,
                         variances=list(residVar=FALSE), needNewNew) {
  locform <- attr(object$predictor,"oriFormula")
  ## possible change of random effect terms
  if ( .noRanef(re.form)) { ## i.e. if re.form implies that there is no random effect
    locform <- nobarsMM(locform)  ## 
  } else if (inherits(re.form,"formula")) { ## ie there is a nontrivial re.form
    lhs <- .DEPARSE(locform[[2L]]) ## response
    fixterms <- .DEPARSE(nobarsMM(locform)[[3L]]) ## fixed effect terms
    ranterms <- parseBars(re.form)  ## random effects 
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
  ## preparation for fixed effects
  allFrames <- HLframes(formula=locform,data=locdata,fitobject=object)
  newX.pv <- allFrames$X ## contains columns for the offset and columns for the other variables
  neednewZAlist <- needNewEta <- ( ( ! is.null(newdata) ) || ! is.null(re.form))
  # newX.pv must intersect non-NA elements of fixef; see comment and code in newetaFix
  validnames <- intersect(colnames(object$X.pv),colnames(newX.pv)) ## we don't want the etaFix cols (detected by bboptim)
  if (length(validnames)==0L) validnames <- c() ## without this, validnames could be character(0) and [,validnames,drop=FALSE] fails.
  RESU <- list(locdata=locdata,newX.pv=newX.pv[,validnames,drop=FALSE]) 
  if (needNewEta) RESU$etaFix <- .newetaFix(object,allFrames,validnames=validnames)  ## new fixed effects (or [part of] new linear predictor if re.form)      
  ## preparation for random effects
  spatial.terms <- findSpatial(locform) ## list of formula terms
  RESU$spatial.model <- spatial.model <- spatial.terms[[1L]] ## one formula term, e.g Matern(1|...)
  if ( ! is.null(newdata)) {
    if ( ! is.null(spatial.model)) { 
      if ((asc <- as.character(spatial.model[[1L]])) %in% c("adjacency","ar1","corrMatrix")) { ## AR1  vs ar1 ?
        stop(paste("Prediction in 'newdata' not implemented or not possible for models including a", asc,"term."))
      } 
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
    ORIspatialOne <- which(attr(object$ZAlist,"ranefs") == spatial.model) ## strictly a spatial one, not other correlated ones
    if (inherits(re.form,"formula")) {
      ranefs <- parseBars(locform) ## note importance for the computation of newZAlist below
      nrand <- length(ranefs)
      newinold <- sapply(ranefs, function(v) {which(v==attr(object$ZAlist,"ranefs"))})
      subZAlist <- object$ZAlist[newinold] ## and reordered
      sublambda <- object$lambda[newinold]
    } else {
      ranefs <- attr(object$ZAlist,"ranefs")
      nrand <- length(ranefs)
      newinold <- seq(nrand) ## keep all ranefs
      subZAlist <- object$ZAlist
      sublambda <- object$lambda
    }    
    RESU$spatialOne <- spatialOne <- which(ranefs == spatial.model) ## strictly a spatial one, not other correlated ones    
    RESU$subZAlist <- subZAlist
    RESU$newinold <- newinold
  }
  #
  if (nrand>0L) {
    if (object$spaMM.version<"1.11.57") {
      strucList <- list(dummyid=attr(object$predictor,"LMatrix")) ## back compat
    } else strucList <- object$strucList
    if (neednewZAlist) {
      FL <- .spMMFactorList(locform, allFrames$mf, 0L, drop=TRUE) 
      newZAlist <- FL$Design ## must be ordered as parseBars result for the next line to be correct.
      attr(newZAlist,"ranefs") <- ranefs ## required pour .compute_ZAXlist to match the ranefs of LMatrix
      corr_list <- .make_corr_list(strucList,newZAlist=newZAlist)
    } else {
      corr_list <- .make_corr_list(strucList,newZAlist=NULL)
      FL <- list(Design=object$ZAlist,namesTerms=object$lambda.object$namesTerms)
      newZAlist <- object$ZAlist
    }
    #############################
    if ( is.null(spatial.model)) {
      uuCnewnew <- NULL
      uuCnewold <- NULL
    } else {
      if (object$spaMM.version<"1.11.57") {
        oldLMatrix <- attr(object$predictor,"LMatrix") ## back compat
      } else oldLMatrix <- object$strucList[[ORIspatialOne]]
      if (is.null(newdata)) {
        uuCnewold <- tcrossprodCpp(oldLMatrix,NULL) ## as in Corr()
        uuCnewnew <- NULL
      } else {
        which <- list(no= needNewEta, nn= needNewNew)
        uunewCorrs <- calcNewCorrs(object=object,locdata=locdata,
                                   which=which,
                                   spatial.model=spatial.model)
        ## matrices, not list of matrices which are constructed later
        uuCnewnew <- uunewCorrs$uuCnewnew
        uuCnewold <- uunewCorrs$uuCnewold
      }
      if ( ! is.null(uuCnewold)) attr(uuCnewold,"ranefs") <- attr(oldLMatrix,"ranefs")    
      corr_list[[ORIspatialOne]] <- uuCnewold
    }
    RESU$newZAClist <- .compute_ZAXlist(XMatrix=corr_list, ZAlist=newZAlist) 
    RESU$newZAlist <- newZAlist
    RESU$uuCnewold <- uuCnewold
    if (needNewNew) RESU$uuCnewnew <- uuCnewnew
    RESU$oldnewClist <- .calc_corr_list(uuCnewold, spatialOne, subZAlist, FL) ## non-trivial value [I] if uuCnewold is NULL
    olduniqueGeo <- attr(object,"info.uniqueGeo")
    geonames <- colnames(olduniqueGeo)
    RESU$newuniqueGeo <- calcUniqueGeo(data=locdata[,geonames,drop=FALSE])
  } 
  return(RESU)
}

## get_predCov_var_fix: see example in predict.Rd (?get_predCov_var_fix)
get_predCov_var_fix <- function(object, newdata = NULL, fix_X_ZAC.object,fixdata, 
                                variances=list(disp=TRUE,residVar=FALSE),...) {
  nrand <- length(fix_X_ZAC.object$newZAlist) 
  if (nrand>1L) {
    fixZAC <- do.call(cbind,fix_X_ZAC.object$newZAClist)
  } else {fixZAC <- fix_X_ZAC.object$newZAClist[[1]]}
  fix_X_ZAC <- cbind2(fix_X_ZAC.object$newX.pv, fixZAC)
  uuCfixold <- fix_X_ZAC.object$uuCnewold
  oldfixClist <- fix_X_ZAC.object$oldnewClist  # cf calcNewCorrs
  fixZAlist <- fix_X_ZAC.object$newZAlist
  spatialOne <- fix_X_ZAC.object$spatialOne
  spatial.model <- fix_X_ZAC.object$spatial.model
  #
  new_X_ZACblob <- .calc_new_X_ZAC(object,newdata=newdata,variances=variances,needNewNew=FALSE) ## called for a correlation block
  if (nrand>1L) {
    newZAC <- do.call(cbind,new_X_ZACblob$newZAClist)
  } else {newZAC <- new_X_ZACblob$newZAClist[[1]]}
  new_X_ZAC <- cbind2(new_X_ZACblob$newX.pv, newZAC)
  uuCnewold <- new_X_ZACblob$uuCnewold
  newuniqueGeo <- new_X_ZACblob$newuniqueGeo
  oldnewClist <- new_X_ZACblob$oldnewClist  # cf calcNewCorrs
  newZAlist <- new_X_ZACblob$newZAlist
  rho <- getPar(object$ranFix,"rho")
  if ( ! is.null(rho_mapping <- attr(object,"dist_info")$rho.mapping)
       && length(rho)>1L ) rho <- fullrho(rho=rho,coordinates=colnames(newuniqueGeo),rho_mapping=rho_mapping)
  msd.arglist <- list(uniqueGeo=newuniqueGeo,uniqueGeo2=fix_X_ZAC.object$newuniqueGeo,
                      rho=rho,return_matrix=TRUE)
  uuCnewfix <- do.call(make_scaled_dist,msd.arglist)
  uuCnewfix <- .calc_corr_from_dist(resu=list(uuCnewold=uuCnewfix), 
                                   object=object, spatial.model=spatial.model, 
                                   which=list(nn=FALSE,no=TRUE))$uuCnewold
  #
  CnewfixList <- lapply(seq_len(nrand),function(it) {
    if (it %in% spatialOne) { ## avoids running the next algo which is slow on large matrices
      return(uuCnewfix)
    } else return(NULL)
  }) ## should make sure it has the same structure as the other matrix lists
  #######################################################
  #
  ## First component of predCov
  # covariance of expectation of Xbeta+Zb due to var of (hat(beta),hat(v)) using E[b] as function of hat(v)
  ## (X_n | C_no) %*% [ t(invL) %*% beta_v_cov[v.range,] %*% invL ] %*% t(X_f | C_fo)
  beta_w_cov <- .get_beta_w_cov(object)
  predVar <- new_X_ZAC[] %id*id% beta_w_cov[] %id*id% t(fix_X_ZAC)[]
  ## Second component of predVar:
  # Evar: expect over distrib of (hat(beta),hat(v)) of [covvariance of Xbeta+Zb given (hat(beta),hat(v))]
  lambda <- object$lambda
  invColdoldList <- .get_invColdoldList(object)
  if (! is.null(CnewfixList) ) {
    Evarlist <- lapply(seq_len(nrand), function(it) {
      isEachNewLevelInOld <- attr(oldnewClist[[it]],"isEachNewLevelInOld")
      if ( ! is.null(isEachNewLevelInOld)) { ## non spatial effect: a vector of booleans indicating whether new is in old (qualifies cols of Cnewold)
        Evar <- Diagonal(x=lambda[it]*as.numeric(! isEachNewLevelInOld))
      } else { ## spatial effect
        Cno_InvCoo_Cof <- t(oldnewClist[[it]])[] %id*id% (invColdoldList[[it]])[] %id*id% oldfixClist[[it]][]
        Evar <- lambda[it] * (CnewfixList[[it]] - Cno_InvCoo_Cof)
      } 
      terme <- newZAlist[[it]][] %id*id% Evar[] %id*id% t(fixZAlist[[it]])[]
      return(as.matrix(terme))
    })
    if (nrand>1L) {Evar <- Reduce("+",Evarlist)} else {Evar <- Evarlist[[1]]}
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
  if ( ! is.null(dwdlogdisp <- logdispObject$dwdlogdisp) ) {
    if ( length(newinold <- fix_X_ZAC.object$newinold) != nrand) { ## ## selection of blocks for re.form ranefs
      cum_n_u_h <- attr(object$lambda.object$lambda,"cum_n_u_h")
      col_info <- attr(dwdlogdisp,"col_info") ## ranefs for which there is a col in dwdlogdisp
      which_ranef_cols <- intersect(col_info$ranef_ids,newinold)
      u_range <- unlist(lapply(which_ranef_cols,function(it) (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])))
      dwdlogdisp <- dwdlogdisp[u_range,c(which_ranef_cols,col_info$phi_cols)]
    }
    newZACw <- newZAC %*% dwdlogdisp ## typically (nnew * n_u_h) %*% (n_u_h * 2) = nnew * 2 hence small 
    fixZACw <- fixZAC %*% dwdlogdisp ## typically (nnew * n_u_h) %*% (n_u_h * 2) = nnew * 2 hence small 
    disp_effect_on_newZACw <- newZACw %*% logdispObject$logdisp_cov %*% t(fixZACw)      
    predVar <- predVar + disp_effect_on_newZACw
  }
  return(predVar)
}
