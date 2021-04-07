.calc_corr_from_dist <- function(distmat, object, corr.model,char_rd) { ## FIXME later we can pass only ranFix$corrPars[[char_rd]]
  if (corr.model=="AR1") {
    corr_mat <- .get_cP_stuff(object$ranFix,"ARphi",which=char_rd)^distmat 
  } else if (corr.model=="Cauchy") {
    corr_mat <- CauchyCorr(shape=.get_cP_stuff(object$ranFix,"shape",which=char_rd),
                           longdep=.get_cP_stuff(object$ranFix,"longdep",which=char_rd),
                           Nugget=.get_cP_stuff(object$ranFix,"Nugget",which=char_rd),
                           d=distmat)  ## so that rho=1 in CauchyCorr
  } else if (corr.model=="Matern") {
    corr_mat <- MaternCorr(nu=.get_cP_stuff(object$ranFix,"nu",which=char_rd),
                           Nugget=.get_cP_stuff(object$ranFix,"Nugget",which=char_rd),
                           d=distmat)  ## so that rho=1 in MaternCorr
  } else stop("Unhandled corr.model")
  return(corr_mat) 
} ## FIXME .calc_corr_from_dist() will ultimately become obsolete, being tied to object$spaMM.version<"2.4.49"


# public wrapper for more transparent workflow
preprocess_fix_corr <- function(object, fixdata, re.form = NULL,
                                variances=list(residVar=FALSE, cov=FALSE), control=list()) {
  delayedAssign("invCov_oldLv_oldLv_list", .get_invColdoldList(object, control=control))
  return(.calc_new_X_ZAC(object=object, newdata=fixdata, re.form = re.form,
                             variances=variances,invCov_oldLv_oldLv_list=invCov_oldLv_oldLv_list) )
}
###############################################

.make_corr_list <- function(strucList, newZAlist) {
  corr_list <- vector("list", length(strucList))
  for (it in seq_len(length(strucList))) {
    if (! is.null(strucList[[it]])) {
      if (attr(strucList[[it]],"corr.model")=="random-coef") {
        if ( ! is.null(newZAlist) || # banal
             ! .spaMM.data$options$replace_design_u # then strucList may not contain expanded-covmat ## in principle only devel case.
             ) {
          abyss <- is.null(longsize <- ncol(newZAlist[[it]])) && (longsize <- ncol(strucList[[it]])) # handling devel case
          corr_list[[it]] <- .makelong(attr(strucList[[it]],"latentL_blob")$compactcovmat,longsize = longsize)
          rownames(corr_list[[it]]) <- colnames(corr_list[[it]]) <- colnames(newZAlist[[it]]) 
          # .tcrossprod gets names from rownames
          ## names required for a predict(,newdata) on a random-coef model
        } else corr_list[[it]] <- .tcrossprod(strucList[[it]],y=NULL, perm=TRUE) # else we can use the original strucList[[it]] which is now and expanded covmat
      } else corr_list[[it]] <- .tcrossprod(strucList[[it]],y=NULL, perm=TRUE)
      attr(corr_list[[it]],"ranefs") <- attr(strucList[[it]],"ranefs")
      attr(corr_list[[it]],"corr.model") <- attr(strucList[[it]],"corr.model")
    }
  } 
  return(corr_list)
}

.calc_sub_diagmat_cov_newLv_oldv <- function(oldZA, newZA,namesTerms) {
  oldlevels <- colnames(oldZA)
  newlevels <- colnames(newZA)
  if (identical(oldlevels,newlevels)) {
    newoldC <- Diagonal(n=length(oldlevels)) ## replaces old identityMatrix
    colnames(newoldC) <- oldlevels ## provides colnames for some XMatrix'es -> some ZAX'es -> used by .match_old_new_levels
  } else {
    numnamesTerms <- length(namesTerms) 
    if (numnamesTerms==1L) {
      nold <- length(oldlevels)
      nnew <- length(newlevels)
      newinold <- match(newlevels,oldlevels)
      matched <- (!is.na(newinold))
      #newoldC <- matrix(0,nrow=nnew,ncol=nold)
      #newoldC[cbind(seq(nnew)[matched],newinold[matched])] <- 1 # cbind() ! don't fill a block !  
      newoldC <- as(sparseMatrix(i=seq(nnew)[matched], j=newinold[matched], dims=c(nnew,nold)),"dgCMatrix")
      colnames(newoldC) <- oldlevels
      rownames(newoldC) <- newlevels
    } else {## random-coef.. but .calc_sub_diagmat_cov_newLv_oldv appears never called in that case...
      nc <- length(oldlevels)
      nr <- length(newlevels)
      bloc_i <- which(newlevels %in% oldlevels)
      bloc_j <- which(oldlevels %in% newlevels)
      lb <- length(bloc_i)
      ii <- rep(bloc_i+(as.integer(gl(numnamesTerms,lb,lb*numnamesTerms))-1L)*nr, numnamesTerms)
      jj <- outer(rep(bloc_j,numnamesTerms), (seq(numnamesTerms)-1L)*nc,"+")
      dim(jj) <- NULL
      newoldC <- as(sparseMatrix(x=rep(1,length(ii)), i=ii, j=jj, dims=c(nr,nc)*numnamesTerms),
                    "dgCMatrix")
      colnames(newoldC) <- rep(oldlevels,numnamesTerms)
      rownames(newoldC) <- rep(newlevels,numnamesTerms)
      if (FALSE) { # NOT the same thing but could be useful elsewhere.
        oldornew <- unique(c(oldlevels,newlevels))
        x <- oldornew %in% intersect(oldlevels,newlevels)
        xx <- rep(x,numnamesTerms)
        nc <- length(x)*numnamesTerms
        newoldC <- as(sparseMatrix(x=rep(1,numnamesTerms^2), i=rep(which(xx),numnamesTerms), j=which(xx)[j],
                                   dims=c(nc,nc)),
                   "dgCMatrix")
        colnames(newoldC) <- rownames(newoldC) <- rep(oldornew, numnamesTerms)
      }
    }
  } 
  attr(newoldC,"isEachNewLevelInOld") <- newlevels %in% oldlevels  ## but this attr is unevaluated (-> NULL) for spatial models 
  return(newoldC)
}

.re_form_col_indices <- function(newinold, cum_n_u_h, Xi_cols) {
  cum_Xi_cols <- cumsum(c(0,Xi_cols))
  ranef_ids <- rep(seq_len(length(Xi_cols)),Xi_cols) ## (repeated for ranCoefs) indices of ranefs, not cols of ranefs
  subrange <- which_ranef_cols <- vector("list",length = length(newinold))
  for (it in seq_len(length(newinold))) {
    rd <- newinold[it] ## allowing possible permutation of ranefs
    subrange[[it]] <- (cum_n_u_h[rd]+1L):(cum_n_u_h[rd+1L])
    which_ranef_cols[[it]] <- which(ranef_ids==rd) ## for ranCoefs, several elements of ranef_ids can match one in newinold
  }
  subrange <- unlist(subrange)
  which_ranef_cols <- unlist(which_ranef_cols)
  return(list(subrange=subrange,which_ranef_cols=which_ranef_cols))
}

.calc_need_Cnn <- function(object, newinold, ori_exp_ranef_types, variances, newZAlist) {
  nrand <- length(newinold)
  need_Cnn <- rep(FALSE, nrand) # i.e. length of newZAlist
  for (new_rd in seq_len(nrand)) { # we don't vectorize to avoid evaluating isDiagonal when it's not needed.  
    old_rd <- newinold[new_rd]
    if ( ! ori_exp_ranef_types[old_rd] %in% c("IMRF","adjacency","corrMatrix")) { 
      need_Cnn[new_rd] <- (attr(object$strucList,"isRandomSlope")[old_rd] | variances$cov ) 
      # that is, we need them also for prediction variances (not cov) for random-slope.
      # and we check another condition where we may need them for prediction variances: 
      if ( ! need_Cnn[new_rd]) {
        if (is.null(is_incid <- attr(newZAlist[[new_rd]],"is_incid"))) {
          message("'is_incid' attribute missing, which suggests inefficient code in .calc_new_X_ZAC().") # should be corrected
          # typically occurred bc there was an A matrix (though not for IMRF) so that ZA differs from Z
          # or by subsetting (ZAlist[[it]][,.] is .preprocess) 
          # or when forgetting to copy the is_incid attribute on .Dvec_times_Matrix(newrd_in_obs, newZlist[[new_rd]])
          #####if (is.null(attr(newZAlist,"AMatrices"))) {
          # awful code potentially huge calculation to detect a nondiagonal element...
          # i could consider that the case is rare enough and not bother
          # or I could consider that the tcrossprod is unlikely to be diagonal...
          need_Cnn[new_rd] <- ( ! isDiagonal(.tcrossprod(newZAlist[[new_rd]]))) 
          # another approach requiring less computation is to check the sum of the values of non-diagonal elements of the tcrossp of abs(ZA), 
          # absZA <- abs(newZAlist[[new_rd]])
          # need_Cnn[new_rd] <- ((sum(rowSums(absZA)^2)-sum(absZA^2))>0)
          #####} else need_Cnn[new_rd] <- TRUE # assumin that the tcrossprod is unlikely to be diagonal...
        } else need_Cnn[new_rd] <- ( ! is_incid)
      }
    } # else for IMRF it remains FALSE bc for IMRF the Evar in nodes is 0 and so is ZA %*% Evar %*% t(ZA) )
  }
  return(need_Cnn)
}

.update_cov_no_nn <- function(RESU, blob, which_mats, newZAlist) {
  # These matrices do not include the lambda factor, as this will be provided in .calc_Evar()...
  cov_newLv_oldv_list <- blob$cov_newLv_oldv_list 
  cov_newLv_newLv_list <- blob$cov_newLv_newLv_list ## may be NULL
  # cov_newLv_oldv_list still contains NULL's, and we need something complete to compute newZAC without '%*% NULL'
  namesTerms <- attr(newZAlist,"namesTerms")
  for (new_rd in seq_along(cov_newLv_oldv_list)) { ## seems correct as everything in rhs is reduced
    if (is.null(cov_newLv_oldv_list[[new_rd]])) { ## excludes ranCoefs, corrMatrix, geostat...
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
  RESU$diag_cov_newLv_newLv_list <- blob$diag_cov_newLv_newLv_list
  return(RESU)
}

.update_newuniqueGeo <- function(info_olduniqueGeo, newinold, need_new_design, locdata) {
  if ( ! is.array(info_olduniqueGeo)) {
    newuniqueGeo <- list() ## to be indexed as the olduniqueGeo
    for (old_rd in newinold) {
      char_old_rd <- as.character(old_rd)
      olduniqueGeo <- info_olduniqueGeo[[char_old_rd]]
      if (!is.null(olduniqueGeo)) {
        coordinates <- colnames(olduniqueGeo)
        geonames <- colnames(olduniqueGeo)
        if (need_new_design) geonames <- intersect(geonames, colnames(locdata)) 
        newuniqueGeo[[char_old_rd]] <- locdata[ ,geonames,drop=FALSE]
      }
    }
  } else {
    olduniqueGeo <- info_olduniqueGeo
    coordinates <- colnames(olduniqueGeo)
    geonames <- colnames(olduniqueGeo)
    if (need_new_design) geonames <- intersect(geonames, colnames(locdata)) 
    newuniqueGeo <- locdata[ ,geonames,drop=FALSE]
  }
  return(newuniqueGeo)
}

# Return with original fixed-effect terms + only shared ranefs 
.update_formula_shared_ranefs <- function(locform, re.form, rm_RHS) {
  if ( .noRanef(re.form)) { ## i.e. if re.form implies that there is no random effect
    locform <- .stripRanefs(locform)  
  } else if (inherits(re.form,"formula")) { ## ie there is a nontrivial re.form
    ranterms <- intersect(.process_bars(re.form), .process_bars(locform)) # It is convenient to handle re.form with more effects than original formula in mv case  
    if (length(ranterms)) {
      lhs <- .DEPARSE(locform[[2L]]) ## response
      fixterms <- .DEPARSE(.stripRanefs(locform)[[3L]]) ## fixed effect terms
      locform <- as.formula(paste(lhs,"~", paste(fixterms,"+",paste(ranterms,collapse="+")) )) # orig fixterms + new ranterms
    } else locform <- .stripRanefs(locform)  
  } ## else keep original formula
  if (rm_RHS && length(locform)==3L) locform <- locform[-2] ## removes RHS for checking  vars on RHS etc
  return(locform)
}

.get_locdata <- function(newdata, locvars=NULL, locform, object, variances) {
  if (is.null(newdata)) {
    return(object$data)
  } 
  # ELSE  
  if (is.null(locvars)) {
    locvars <- all.vars(.strip_IMRF_args(locform)) ## strip to avoid e.g. 'stuff' being retained as a var from IMRF(..., model=stuff)
    if (variances$residVar) locvars <- unique(c(locvars,all.vars(.strip_IMRF_args(.get_phiform(object)))))  
  }
  #
  if (is.vector(newdata)) { ## ## less well controlled case, but useful for maximization
    locdata <- data.frame(matrix(newdata,nrow=1))
    if (length(locvars)==ncol(locdata)) {
      names(locdata) <- locvars
    } else {
      stop(paste("(!) newdata has incorrect length. It should match the following variables:\n",paste(locvars,collapse=" ")))
    }
  } else {
    if( is.matrix(newdata) ) newdata <- as.data.frame(newdata)  
    # so that matrix 'newdata' arguments can be used as in some other predict methods.
    locdata <- try(newdata[ , locvars,drop=FALSE]) ## allvars checks only RHS variables
    if (inherits(locdata,"try-error")) stop(paste0("Variable(s) ",
                                                   paste(setdiff(locvars,colnames(newdata)), collapse=","),
                                                   " appear to be missing from newdata."))
  }
  # => for any non-NULL newdata, locadta is a data.frame with valid variables. Check NAs:
  checkNAs <- apply(locdata,1L,anyNA)
  if (any(checkNAs)) {
    message("NA's in required variables from 'newdata'. Prediction not always possible.")
    locdata <- locdata[!checkNAs, ]
  }
  locdata
}

.calc_newFrames_fixed <- function (formula, data, fitobject, need_allFrames=TRUE, mv_it=NULL) {
  ## X may or may not contain offset info, which should not be used (see .newEtaFix()) 
  #  but fixef_mf should contain such info bc .newEtaFix calls off <- model.offset( newMeanFrames$mf)
  if (is.null(formula)) {
    X <- matrix(nrow=nrow(data),ncol=0L) ## model without fixed effects, not even an Intercept 
    fixef_mf <- NULL 
  } else { 
    fixef_off_form <- .stripRanefs_(formula) 
    if (inherits(fixef_off_form, "formula")) {
      if (is.character(formula[[2L]])) fixef_off_form <- fixef_off_form[-2L] ## something like ".phi" ....
      Terms <- terms(fixef_off_form)
      Terms <- stats::delete.response(Terms)
      attr(Terms,"predvars") <- .calc_newpredvars(fitobject, fixef_off_form) ## for poly()
      fixef_form <- .stripOffset_(fixef_off_form) # formula if something remains after the offset has been removed
      if ( ! inherits(fixef_form, "formula")) { ## only an offset in formula, not even an explicit 0: .stripOffset_(fixef_off_form) produced a 'name'
        attr(Terms,"intercept") <- 0L # removes offset that terms() assumes if there is no explicit '0'.
      }
      # handles offset:  (without the small shortcut used in .get_terms_info())
      fixef_mf <- model.frame(Terms, data, xlev = .get_from_terms_info(object=fitobject, which="fixef_levels", mv_it=mv_it)) ## xlev gives info about the original levels
      # :here for a poly(age,.) Terms and age=Inf in the 'data', fixef_mf had zero rows and linkinv will fail on numeric(0) 'eta'
      if (nrow(fixef_mf)!=nrow(data)) {
        if (any(pb <- which( ! sapply(lapply(data,is.finite), all)))) {
          stop(paste0("NA/NaN/Inf in 'data' for fixed-effects prediction: check variable(s) '", paste(names(pb), collapse="', '"),"'."))
        } else stop("nrow(fixef_mf)!=nrow(data) for undetermined reason") 
      }
      X <- model.matrix(Terms, fixef_mf, contrasts.arg=attr(fitobject$X.pv,"contrasts")) 
      ## : where original contrasts definition is used to define X cols that match those of original X, whatever was the contrast definition when the model was fitted
      if ( ! is.null(mv_it)) {
        cum_ncol_X <- attr(fitobject$X.pv,"cum_ncol")
        col_range <- cum_ncol_X[mv_it]+seq_len(cum_ncol_X[mv_it+1L]-cum_ncol_X[mv_it]) # avoid (n+1:n) problem 
        colnames(X) <- colnames(fitobject$X.pv)[col_range]
      }
    } else {
      fixef_mf <- NULL
      X <- matrix(nrow=nrow(data), ncol=0L) ## Not using NROW(Y) which is 0 if formula has no LHS
    }
  }
  storage.mode(X) <- "double" ## otherwise X may be logi[] rather than num[] in particular when ncol=0
  return(list(X = X, mf = fixef_mf)) 
}

if (FALSE) { # v3.5.121 managed to get rid of it 
  .calc_newFrames_ranef <- function (formula, data, fitobject) {
    formula <- .asNoCorrFormula(formula) ## strips out the correlation information, retaining the ranefs as (.|.)
    if (is.character(formula[[2L]])) formula <- formula[-2L] ## something like ".phi" ....
    plusForm <- .subbarsMM(formula) ## this comes from lme4 and converts (.|.) terms to (.+.) form 
    environment(plusForm) <- environment(formula)
    Terms <- terms(plusForm) ## assumes an Intercept implicitly
    Terms <- stats::delete.response(Terms)
    #attr(Terms,"predvars") <- .calc_newpredvars(fitobject$main_terms_info$all_terms, Terms) ## for poly in ranefs ? 
    mf <- model.frame(Terms, data, drop.unused.levels=TRUE) 
    return(list(mf = mf))
  }
  # } else {
  #   .calc_newFrames_ranef <- function (formula, data, fitobject) {list(mf=data)}
}


.get_newX_info <- function(locform, locdata, object, mv_it=NULL) {
  newFrames_fixed <- .calc_newFrames_fixed(formula=.stripRanefs(locform),data=locdata,fitobject=object, mv_it=mv_it) ## also used for predVar computations
  ## preparation for fixed effects
  newX.pv <- newFrames_fixed$X ## contains columns for the offset and columns for the other variables
  # newX.pv must intersect non-NA elements of fixef; see comment and code in newetaFix
  est_and_fix <- names(which(!is.na(object$fixef)))
  validnames <- intersect(est_and_fix,colnames(newX.pv)) ## we don't want the etaFix cols (detected by bboptim)
  if (length(validnames)==0L) validnames <- c() ## without this, validnames could be character(0) and [,validnames,drop=FALSE] fails.
  if (length(notfound <- setdiff(colnames(newX.pv), est_and_fix))) {
    # capture case where the newX.pv has colnames  not in object$X.pv (including weird case of mis-naming)
    stop(paste0("No fitted coefficient(s) for variables\n",paste(notfound,collapse=", "),"\nin the design matrix derived from 'newdata'."))
  }
  return(list(newX.pv=newX.pv[ , validnames,drop=FALSE],
              eta_fix=.newetaFix(object,newFrames_fixed,validnames=validnames)
              )
         ) 
}

.get_newinold <- function(re.form, locform, ori_exp_ranef_strings, rd_in_mv=NULL) {
  if (is.null(rd_in_mv) || ## assumping that in the univariate-response case, we call .get_newinold() only when we already tested inherits(re.form,"formula")
      inherits(re.form,"formula")) {
    new_exp_ranef_strings <- .process_bars(locform,expand=TRUE) ## to be added as attribute to the newZAlist created by .calc_normalized_ZAlist()
    newinold <- unlist(sapply(lapply(new_exp_ranef_strings, `==`, y= ori_exp_ranef_strings), which)) ## e.g 1 4 5
    if ( ! is.null(rd_in_mv)) newinold <- rd_in_mv[newinold]
    # : unlist bc sapply(list(F), which) return list(numeric(0)), etc.
  } else {
    new_exp_ranef_strings <- ori_exp_ranef_strings
    if ( ! is.null(rd_in_mv)) {
      newinold <- rd_in_mv
    } else {
      newinold <- seq_along(new_exp_ranef_strings) ## keep all ranefs
    }
  }    
  newinold
}

.get_newfixef_info <- function(newdata, locform, locdata, object, re.form) {
  if (  ! is.null(newdata) ) {
    RESU <- .get_newX_info(locform, locdata, object) # includes element 'eta_fix' which includes the offset
    RESU$locdata <- locdata
  } else if (! is.null(re.form)) { # then I still need eta_fix
    eta_fix <- drop(object$X.pv %*% fixef(object))
    off <- model.offset(model.frame(object)) # so if I later store the model.frame I would no longer need the raw data here.
    if (!is.null(off)) eta_fix <- eta_fix+off
    RESU <- list(locdata=locdata,newX.pv=object$X.pv, eta_fix=eta_fix)
  } else RESU <- list(locdata=locdata,newX.pv=object$X.pv) 
  RESU
}

# Currently never called for mv: cf .calc_new_X_ZAC_mv() instead
.calc_new_X_ZAC <- function(object, newdata=NULL, re.form = NULL,
                            variances=list(residVar=FALSE, cov=FALSE),invCov_oldLv_oldLv_list,
                            locform=formula.HLfit(object, which="")) {
  ## possible change of random effect terms
  locform <- .update_formula_shared_ranefs(locform, re.form, rm_RHS=TRUE)
  # checking variables in the data
  locdata <- .get_locdata(newdata, locform=locform, object=object, variances=variances) # always required (cf when newdata is not a data frame)
  #
  RESU <- .get_newfixef_info(newdata, locform, locdata, object, re.form)
  #
  ## subZAlist is a subset of the old ZA, newZAlist contains new ZA ; different uses=> computation required under disticnt conditions for each.
  ## calling .make_corr_list(object$strucList,...) is always OK bc the fist argument may be a superset of the required list
  ## all matching in .make_corr_list is through the ranef attributes.
  #
  ## matching ranef terms of re.form
  if ( ! .noRanef(re.form)) { # (which is by default TRUE, including for GLMs... but we donc care about optimizing code for GLMs) 
    if (object$spaMM.version < "2.2.116") {
      ori_exp_ranef_strings <- attr(object$ZAlist,"ranefs") 
      ori_exp_ranef_types <- attr(ori_exp_ranef_strings,"type")
      # next line is a long-after guess. We need 'ori_exp_ranef_terms' to simplify some code below.
      ori_exp_ranef_terms <- .process_bars(barlist=structure(ori_exp_ranef_strings, type=ori_exp_ranef_types),
                                           expand=FALSE, as_character=FALSE, which.="exp_ranef_terms")
    } else {
      ori_exp_ranef_terms <- attr(object$ZAlist,"exp_ranef_terms")
      ori_exp_ranef_strings <- attr(object$ZAlist,"exp_ranef_strings")
      ori_exp_ranef_types <- attr(object$ZAlist,"exp_ranef_types") 
    }
    RESU$spatial_old_rd <- which(ori_exp_ranef_types != "(.|.)")   
    #
    if (inherits(re.form,"formula")) {
      newinold <-.get_newinold(re.form, locform, ori_exp_ranef_strings, rd_in_mv=NULL)
      new_exp_ranef_strings <- ori_exp_ranef_strings[newinold]
      RESU$subZAlist <- object$ZAlist[newinold] ## and reordered
    } else {
      newinold <- seq_along(ori_exp_ranef_strings)
      new_exp_ranef_strings <- ori_exp_ranef_strings
      RESU$subZAlist <- object$ZAlist
    }
    RESU$newinold <- newinold
    #
    if (nrand <- length(newinold)) {
      strucList <- object$strucList
      if (object$spaMM.version<"1.11.57") stop("This fit object was created with spaMM version<1.11.57, and is no longer supported.\n Please recompute it.")
      need_new_design <- ( ( ! is.null(newdata) ) || ! is.null(re.form)) ## newdata or new model
      if (need_new_design) {
        ## with newdata we need Evar and then we need nn... if newdata=ori data the Evar (computed with the proper nn) should be 0
        #barlist <- .process_bars(locform,as_character=FALSE) 
        #new_mf_ranef <- .get_new_mf_ranef(barlist=barlist, locdata, object) 
        new_exp_ranef_terms <- structure(ori_exp_ranef_terms[newinold], type=attr(ori_exp_ranef_terms,"type")[newinold])
        ranef_form <- as.formula(paste("~",(paste(new_exp_ranef_strings,collapse="+")))) ## effective '.noFixef'
        newZlist <- .calc_Zlist(exp_ranef_terms=new_exp_ranef_terms, # .process_bars(barlist=barlist,as_character=FALSE, which.="exp_ranef_terms"), # != barlist, for IMRF notably
                                #locform, 
                                data=locdata, rmInt=0L, drop=TRUE,sparse_precision=FALSE, 
                                levels_type= "seq_len", ## superseded in specific cases: notably, 
                                ## the same type has to be used by .calc_AMatrix_IMRF() -> .as_factor() 
                                ##  as by .calc_Zmatrix() -> .as_factor() for IMRFs.
                                ## This is controlled by option uGeo_levels_type (default = "mf" as the most explicit; using ".ULI" appears OK).
                                ## The sames functions are called with the same arguments for predict with newdata.
                                corrMats_info=object$strucList, ## use is to provide info about levels in .calc_ZMatrix()
                                old_ZAlist=object$ZAlist, newinold=newinold, #barlist=barlist, 
                                lcrandfamfam=attr(object$rand.families,"lcrandfamfam")) 
        amatrices <- .get_new_AMatrices(object, newdata=locdata) # .calc_newFrames_ranef(formula=ranef_form,data=locdata,fitobject=object)$mf)
        ## ! complications:
        ## even if we used perm_Q for Matern, the permutation A matrix should not be necessary 
        ##  in building the new correlation matrix, although it night be used as well 
        ## explict colnames should handle both cases, so that
        ## newZAlist <- .calc_normalized_ZAlist( ignoring those A matrices)
        ## and 
        ## newZAlist <- object$ZAlist
        ## should be OK.
        ## But the other Amatrices should be processed before newZACpplist <- .compute_ZAXlist(.) is called
        requires_ZCpL <- (attr(newZlist,"exp_ranef_types") %in% c("Matern","Cauchy"))
        newZAlist <- .calc_normalized_ZAlist(Zlist=newZlist,
                                             # newZlist has names not necessarily starting at "1"
                                             AMatrices=amatrices[names(newZlist)[ ! requires_ZCpL]],
                                             vec_normIMRF=object$ranef_info$vec_normIMRF, 
                                             strucList=strucList)
        ## must be ordered as parseBars result for the next line to be correct.
        attr(newZAlist,"exp_ranef_strings") <- new_exp_ranef_strings ## required pour .compute_ZAXlist to match the ranefs of LMatrix
      } else {
        newZAlist <- object$ZAlist
      }
      RESU$newZAlist <- newZAlist
      # We determine which matrices we need for computation of Evar:
      need_Cnn <- .calc_need_Cnn(object, newinold, ori_exp_ranef_types, variances, newZAlist)
      which_mats <- list(no= need_new_design, 
                         ## cov_newLv_newLv_list used in .calc_Evar() whenever newdata, but elements may remain NULL if $cov not requested
                         ## However, for ranCoefs, we need Xi_cols rows for each response's predVar. (FIXME) we store the full matrix.
                         nn= need_Cnn ) 
      #nrand <- length(newinold)
      #
      ## AT this point both newZAlist and subZAlist may have been reduced to 'newnrand' elements relative to ori object$ZAlist.
      ## .match_old_new_levels will use it running over newnrand values
      ## The cov_ lists are reduced too. newinold should be used to construct them
      ## newZACpplist is reduced.
      if (any(unlist(which_mats)) 
          || any(unlist(variances)) # cov_newLv_oldv_list is always needed for cbind(X.pv,newZAC [which may be ori ZAC]); should ~corr_list when newdata=ori data
      ) {
        blob <- .make_new_corr_lists(object=object,locdata=locdata, which_mats=which_mats, 
                                     newZAlist=newZAlist, newinold=newinold,
                                     invCov_oldLv_oldLv_list=invCov_oldLv_oldLv_list)
        RESU <- .update_cov_no_nn(RESU, blob, which_mats, newZAlist)
        RESU$newZACpplist <- .compute_ZAXlist(XMatrix=RESU$cov_newLv_oldv_list, ZAlist=newZAlist) ## build from reduced list, returns a reduced list
        ## This $newZACpplist serves to compute new _point predictions_.
        #  # this comment may be obsolete : .compute_ZAXlist affects elements of ZAlist that have a ranefs attribute. 
        #  It builds a design matrix to all oldv levels. It does not try to reduce levels. 
        #  But it use newZlist <- .calc_Zlist(...) which may consider a reduced number of levels.
        #  The function .match_old_new_levels() will perform the match with 'w_h_coeffs' for point prediction.
        #  it assign values psi_M to new levels of ranefs.
        ## _Alternatively_ one may construct a newZACvar for _predVar_ 
        #  Here we have a slice mechanism (contrary for point pred) hence new ZA with different rows, and the columns of the 
        #  newZACvar constructed there must match those of beta_w_cov for  ZWZt_mat_or_diag( <cbind(newX,newZAC)> ,beta_w_cov)
        #  cov_newLv_oldv_list() provides the info for the expansion from the newZA cols to the oldZA cols.
        #  In that case one does not need to match levels. .calc_newZACvar() performs a simpler operation than .compute_ZAXlist.
        if (need_new_design) RESU$newZACpplist <- .calc_ZAlist(Zlist=RESU$newZACpplist, AMatrices=amatrices[names(newZAlist)[requires_ZCpL]])
      }
      #
      ## attribute added in version 2.3.18:
      info_olduniqueGeo <- attr(object,"info.uniqueGeo")
      if ( ! is.null(info_olduniqueGeo)) RESU$newuniqueGeo <- .update_newuniqueGeo(info_olduniqueGeo, 
                                                                                   newinold, need_new_design, locdata)
      # Despite the 'newuniqueGeo' name, it may be necess without newdata; preprocess_fix_corr() with NULL 'fixdata' providing info_olduniqueGeo <- fix_info$newuniqueGeo
    }
  }
  return(RESU)
}

## get_predCov_var_fix: see example in predict.Rd (?get_predCov_var_fix), test in test-predVar 
# get_predCov_var_fix -> .calc_new_X_ZAC -> evaluates 'which_mats' according to all relevant arguments
get_predCov_var_fix <- function(object, newdata = NULL, fix_X_ZAC.object,fixdata, re.form = NULL, 
                                variances=list(disp=TRUE,residVar=FALSE,cov=FALSE), control=list(), ...) {
  delayedAssign("invCov_oldLv_oldLv_list", .get_invColdoldList(object, control=control))
  variances <- .process_variances(variances, object)
  newnrand <- length(fix_X_ZAC.object$newZAlist) 
  fixZACvar <- .calc_newZACvar(fix_X_ZAC.object$newZAlist,fix_X_ZAC.object$cov_newLv_oldv_list)
  new_X_ZACblob <- .calc_new_X_ZAC(object,newdata=newdata,variances=variances,
                                   invCov_oldLv_oldLv_list=invCov_oldLv_oldLv_list) ## called for a correlation block
  newZACvar <- .calc_newZACvar(new_X_ZACblob$newZAlist,new_X_ZACblob$cov_newLv_oldv_list)
  ## First component of predVar
  # covariance of expectation of Xbeta+Zb due to var of (hat(beta),hat(v)) using E[b] as function of hat(v)
  ## (X_n | C_no) %*% [ t(invL) %*% beta_v_cov[v.range,] %*% invL ] %*% t(X_f | C_fo)
  fix_X_ZAC <- cbind2(fix_X_ZAC.object$newX.pv, fixZACvar)
  new_X_ZAC <- cbind2(new_X_ZACblob$newX.pv, newZACvar)
  loc_tcrossfac_beta_w_cov <- .get_tcrossfac_beta_w_cov(object)
  if (inherits(re.form,"formula")) {
    # identifies and selects columns for the [retained ranefs, which are given by newinold 
    ori_exp_ranef_strings <- attr(object$ZAlist,"exp_ranef_strings")
    new_exp_ranef_strings <- .process_bars(re.form,expand=TRUE) 
    newinold <- sapply(lapply(new_exp_ranef_strings, `==`, y= ori_exp_ranef_strings), which) ## e.g 1 4 5
    re_form_col_indices <- .re_form_col_indices(newinold, cum_n_u_h=attr(object$lambda,"cum_n_u_h"), Xi_cols=attr(object$ZAlist, "Xi_cols"))
    Xncol <- ncol(object$X.pv)
    subrange <- c(seq_len(Xncol),re_form_col_indices$subrange + Xncol)
    loc_tcrossfac_beta_w_cov <- loc_tcrossfac_beta_w_cov[subrange,]
  } else re_form_col_indices <- NULL
  if (variances$naive) {
    naive <- .calc_Var_given_fixef(object, new_X_ZACblob=new_X_ZACblob, covMatrix=variances$cov, fix_X_ZAC.object=fix_X_ZAC.object)
    return(naive)
  } 
  ## get_predCov_var_fix() is typically called once (if it is) so no use in saving beta_w_cov
  predVar <- new_X_ZAC %id*% .tcrossprod(loc_tcrossfac_beta_w_cov) %*id% t(fix_X_ZAC)
  ## Second component of predVar:
  cov_newLv_fixLv_list <- .make_new_corr_lists(object, locdata=new_X_ZACblob$newuniqueGeo, 
                                               which_mats=list(no=TRUE,nn=rep(FALSE,newnrand)), # all the covariance info being provided by the fix_X_ZAC.object 
                                               new_X_ZACblob$newZAlist, 
                                               newinold=fix_X_ZAC.object$newinold,
                                               fix_info=fix_X_ZAC.object,
                                               invCov_oldLv_oldLv_list=invCov_oldLv_oldLv_list
                                               )$cov_newLv_oldv_list
  if ( ! is.null(cov_newLv_fixLv_list) ) {
    # Evar: expect over distrib of (hat(beta),new hat(v)) of [covariance of Xbeta+Zb given (hat(beta),orig hat(v))]
    Evar <- .calc_Evar(newZAlist=new_X_ZACblob$newZAlist,newinold=fix_X_ZAC.object$newinold, 
                       cov_newLv_oldv_list=new_X_ZACblob$cov_newLv_oldv_list, 
                       invCov_oldLv_oldLv_list=invCov_oldLv_oldLv_list, 
                       cov_newLv_fixLv_list=cov_newLv_fixLv_list, cov_fixLv_oldv_list=fix_X_ZAC.object$cov_newLv_oldv_list, 
                       fixZAlist=fix_X_ZAC.object$newZAlist,covMatrix=TRUE,
                       diag_cov_newLv_newLv_list=fix_X_ZAC.object$diag_cov_newLv_newLv_list,
                       object=object) # we pass object to be able to assign it its environment
    predVar <- predVar + Evar
  } 
  # If components for uncertainty in dispersion params were requested,
  #   logdispObject is not NULL
  # If some components ere computable, $$dwdlogdisp should not be NULL
  # Former approach (changed 08/2016) was to test logdispObject and then 
  #   for any 'problem'. But there may some 'problem' and still a valid logdispObject
  # => In new version, dwdlogdisp should be either NULL or a conforming matrix;
  #  'problems" should not be tested.
  logdispObject <- .get_logdispObject(object)
  if (variances$disp && ! is.null(logdispObject$dwdlogdisp) ) {
    dwdlogdisp <- logdispObject$dwdlogdisp
    logdisp_cov <- logdispObject$logdisp_cov ## idem
    phi_cols <- attr(dwdlogdisp,"col_info")$phi_cols ## make local copy before subsetting the matrix!
    if ( ! is.null(re_form_col_indices) ) { ## selection of blocks for re.form ranefs 
      whichcols <- c(re_form_col_indices$which_ranef_cols, phi_cols)
      dwdlogdisp <- dwdlogdisp[re_form_col_indices$subrange,whichcols] ## permuted ranefs => permuted rows and cols
      logdisp_cov <- logdisp_cov[whichcols,whichcols] ## permuted ranefs => permuted rows and cols
    } 
    if (!is.null(hyper_info <- .get_from_ranef_info(object, "hyper_info"))) {
      summingMat <- hyper_info$summingMat
      if (!is.null(re_form_col_indices)) {
        summingMat <- summingMat[newinold, , drop = FALSE]
        colids <- numeric(nrow(summingMat))
        for (it in seq(nrow(summingMat))) colids[it] <- which(summingMat[it, ] > 0)
        colids <- unique(colids)
        summingMat <- summingMat[, colids, drop = FALSE]
      }
      summingMat <- as.matrix(Matrix::bdiag(summingMat, rep(1, length(phi_cols))))
      dwdlogdisp <- dwdlogdisp %*% summingMat
      logdisp_cov <- t(summingMat) %*% logdisp_cov %*% summingMat
    }
    # newZACvar = (ZAC_ranef1 | ZAC_ranef3... ) %*% dwdlogdisp which rows match the successive v_h (all ranefs) and cols match disp pars
    newZACw <- newZACvar %*% dwdlogdisp ## typically (nnew * n_u_h) %*% (n_u_h * 2) = nnew * 2 hence small 
    fixZACw <- fixZACvar %*% dwdlogdisp ## typically (nnew * n_u_h) %*% (n_u_h * 2) = nnew * 2 hence small 
    disp_effect_on_newZACw <- .get_disp_effect_on_newZACw(logdisp_cov, newZACw, fixZACw=fixZACw, covMatrix=TRUE) # newZACw %*% logdisp_cov %*% t(fixZACw)
    predVar <- predVar + disp_effect_on_newZACw
  }
  return(predVar)
}
