.calc_new_X_ZAC_mv <- function(object, newdata=NULL, re.form = NULL,
                            variances=list(residVar=FALSE, cov=FALSE),invCov_oldLv_oldLv_list,
                            control=list()) {
  locformS <- formula.HLfit(object, which="")
  if (inherits(re.form, "formula")) {
    re.formS <- vector("list", length(locformS))
    for (mv_it in seq_along(locformS)) re.formS[[mv_it]] <- .update_formula_shared_ranefs(re.form, locformS[[mv_it]], rm_RHS=TRUE)
  } else if (length(re.form)==1L && is.na(re.form)) {
    re.formS <- rep(NA, length(locformS))
  } else re.formS <- re.form # a list or NULL. 
  ## possible change of random effect terms
  for (mv_it in seq_along(locformS)) locformS[[mv_it]] <- .update_formula_shared_ranefs(locformS[[mv_it]], re.formS[[mv_it]], rm_RHS=TRUE)
  # checking variables in the data
  allvarsS <- vector("list", length(locformS))
  for (mv_it in seq_along(locformS)) {
    allvars_it <- all.vars(.strip_cF_args(locformS[[mv_it]])) ## strip to avoid e.g. 'stuff' being retained as a var from IMRF(..., model=stuff)
    if (variances$residVar || control$simulate) {
      allvars_it <- unique(c(allvars_it, all.vars(.strip_cF_args(.get_phiform(object, mv_it)))))
      allvars_it <- unique(c(allvars_it, all.vars(object$families[[mv_it]]$resid.model$resid.formula)))
    } 
    allvarsS[[mv_it]] <- allvars_it
  }
  #
  ## matching ranef terms of re.form
  if (length(object$ZAlist)) { 
    if (identical(control$keep_ranef_vars_for_simulate, TRUE) || # : condition for the case 
                      # where only eta_fixed is predicted for marginal simulation, hence re.form is NA ("no  prediction for ranef") BUT 
                      # we will need the locdataS with the variables for ranefs, to simulate these ranefs.
         ( re.form_allows_ranefs <- ( is.null(re.form) || # : this condition means that re.form_allows_ranefs may be TREU despite no ranef in model formula 
                                      any( ! sapply(re.formS, .noRanef))) )
       ) { 
      map_rd_mv <- attr(object$ZAlist, "map_rd_mv")
      ori_exp_ranef_strings <- attr(object$ZAlist,"exp_ranef_strings")
      #
      newinoldS <- vector("list", length(locformS))
      for (mv_it in seq_along(locformS)) newinoldS[[mv_it]] <- .get_newinold(re.formS[[mv_it]], locformS[[mv_it]], 
                                                                             ori_exp_ranef_strings, rd_in_mv=map_rd_mv[[mv_it]])
      newinold <- unique(.unlist(newinoldS))
    } else newinold <- NULL 
  } else newinold <- NULL
  #
  if (length(newinold)) {
    new_exp_ranef_strings <- ori_exp_ranef_strings[newinold]
    ranef_form <- as.formula(paste("~",(paste(new_exp_ranef_strings,collapse="+")))) ## effective '.noFixef'
    ranefvars <- all.vars(.strip_cF_args(ranef_form))
  } else ranefvars <- c()
  locdata <- .get_locdata(newdata=newdata, locvars=unique(c(.unlist(allvarsS),ranefvars)), 
                          object=object, variances=variances, na.rm=FALSE) # see comment on evaluation of RESU$newuniqueGeo
  #
  locdataS <- vector("list", length(allvarsS))
  need_new_design <- ( ( ! is.null(newdata) ) || ! is.null(re.form)) ## newdata or new model
  if (need_new_design) {
    newX.pv <- eta_fix <- NULL
    for (mv_it in seq_along(locformS)) {
      allvars_it <- unique(c(allvarsS[[mv_it]], ranefvars)) # all ranefvars + mv_it-specific fixefvars: note difference if ! need_new_design
      # presence of ranefvars for all ranefs in each allvars_it allows do.call(rbind, locdataS) later in this fn after selecting them
      locdataS[[mv_it]] <- ..get_locdata(locdata, allvars_it, na.rm=TRUE)
      newX_info <- .get_newX_info(locformS[[mv_it]], locdataS[[mv_it]], object, mv_it=mv_it)
      newX.pv <- .merge_Xs(newX.pv, newX_info$newX.pv, mv_it)
      eta_fix <- c(eta_fix, newX_info$eta_fix)
    }
    attr(locdataS,"allvarsS") <- allvarsS # allVarsS (used by map_ranef()) tells which variables were needed for which submodel predictions.
    RESU <- list(locdata=locdataS, # locdata in RESU for cbind()ing the predictions in .predict_body();
                 cum_nobs= cumsum(c(0L,lapply(locdataS, nrow))), # more widely used by .fv_linkinv()
                 newX.pv=newX.pv, eta_fix=eta_fix) 
  } else {
    for (mv_it in seq_along(allvarsS)) locdataS[[mv_it]] <-  ..get_locdata(locdata, allvarsS[[mv_it]], na.rm=TRUE)
    RESU <- list(locdata=locdataS, # locdata in RESU allowing (potential) cbind() with predictions in .predict_body(). (maybe not implemented for mv)
                 cum_nobs= cumsum(c(0L,lapply(locdataS, nrow))), # more widely used by .fv_linkinv()
                 newX.pv=model.matrix(object)) 
  } 
  ## so we save 'locdata=locdataS' in RESU, but we will locally modify 'locdataS' by selecting columns in locdataS[[mv_it]]
  ## and will create a local 'ranefdata' from this locally modified 'locdataS'.
  #
  ## newZAlist and subZAlist appear to have distinct usages since they are created under different conditions.
  ## subZAlist is a subset of the old ZA, newZAlist contains new ZA
  #
  if (nrand <- length(newinold)) {  
    if (length( info_olduniqueGeo <- .get_old_info_uniqueGeo(object) )) {
      RESU$newuniqueGeo <- .update_newuniqueGeo(info_olduniqueGeo, newinold, need_new_design, locdata)
      # Despite the $newuniqueGeo name, it may be necess without newdata; 
      # cf preprocess_fix_corr() with NULL 'fixdata' providing info_olduniqueGeo <- fix_info$newuniqueGeo (univariate case).
      # There may be a slight suboptimality as it uses 'locdata' produced with na.rm=FALSE. 
      # But na.rm=TRUE would remove rows used in a submodel, if they have NA's in variables not used in that submodel. 
    }
    RESU$subZAlist <- object$ZAlist[newinold] ## and reordered
    RESU$newinold <- newinold
    ori_exp_ranef_types <- attr(object$ZAlist,"exp_ranef_types") 
    RESU$spatial_old_rd <- which(ori_exp_ranef_types != "(.|.)") 
    
    strucList <- object$strucList
    if (need_new_design) {
      ## with newdata we need Evar and then we need nn... if newdata=ori data the Evar (computed with the proper nn) should be 0
      # barlist <- .process_bars_mv(predictors=formula.HLfit(object,which = ""),
      #                             map_rd_mv=map_rd_mv,
      #                             as_character=FALSE) ## but default expand =TRUE 
      # barlist <- structure(barlist[newinold], type=attr(barlist,"type")[newinold])
      
      for (mv_it in seq_along(locdataS)) locdataS[[mv_it]] <- locdataS[[mv_it]][,ranefvars, drop=FALSE]
      ranefdata <- do.call(rbind, locdataS)
      ori_exp_ranef_terms <- attr(object$ZAlist,"exp_ranef_terms")
      new_exp_ranef_terms <- structure(ori_exp_ranef_terms[newinold], type=attr(ori_exp_ranef_terms,"type")[newinold])
      newZlist <- .calc_Zlist(exp_ranef_terms=new_exp_ranef_terms, # .process_bars(barlist=barlist,as_character=FALSE, which.="exp_ranef_terms"), # != barlist, for IMRF notably
                              data=ranefdata, 
                              rmInt=0L, drop=TRUE,sparse_precision=FALSE, 
                              corr_info=.get_from_ranef_info(object),
                              levels_type= "seq_len", ## superseded in specific cases: notably, 
                              ## the same type has to be used by .calc_AMatrix_IMRF() -> .as_factor() 
                              ##  as by .calc_Zmatrix() -> .as_factor() for IMRFs.
                              ## This is controlled by package option 'uGeo_levels_type' (default = "data_order" as the most explicit).
                              ## The sames functions are called with the same arguments for predict with newdata.
                              ##
                              ## This means that if there are repeated geo positions in the newdata 
                              ## we save the time of trying to find them (which perhaps is less justifiable in mv case? __FIXME__)
                              sub_oldZAlist=object$ZAlist[newinold],  
                              lcrandfamfam=attr(object$rand.families,"lcrandfamfam"))
      # Each Z in newZlist then has non-zero rows even for response levels that are not affected by the ranef
      # bc it's built from 'ranefdata' built as if all submodels were affected by each (new) ranef
      # => need to correct this
      newrd_in_obsS <- vector("list", length(newinoldS))
      loc_cum_nobs <- RESU$cum_nobs
      for(new_rd in seq_along(newinold)) {
        for (mv_it in seq_along(newinoldS)) {
          nobs_it <- loc_cum_nobs[[mv_it+1L]] - loc_cum_nobs[[mv_it]]
          newrd_in_obsS[[mv_it]] <- rep((newinold[new_rd] %in% newinoldS[[mv_it]]), nobs_it)
        }
        newrd_in_obs <- 1L * .unlist(newrd_in_obsS)
        newZlist[[new_rd]] <- structure(.Dvec_times_Matrix(newrd_in_obs, newZlist[[new_rd]]),
                                        is_incid=attr(newZlist[[new_rd]],"is_incid"))
      }
      amatrices <- .get_new_AMatrices(object, newdata=ranefdata) # .calc_newFrames_ranef(formula=ranef_form,data=ranefdata,fitobject=object)$mf)
      ## ! complications:
      ## even if we used perm_Q for *fitting* Matern , the permutation A matrix should not be necessary in building the new correlation matrix, bc ./.
      ## explicit colnames should handle both cases, so that
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
      newZAlist <- .correct_newZAlist_mv_ranCoefs(ZAlist=newZAlist, 
                                                  cum_nobs=cumsum(c(0L,sapply(locdataS,nrow))))
      # : distinct .correct_...() function needed here bc [
      #   we cannot run an outer loop on mv_it, calling.merge_ZAlists(., ZAlist2 = newZAlist_i...) bc
      #     expected format of final newZAlist differs from that for the fit: for some ranef types
      #     we do not try to match the columns of ZA_i's different submodels, instead building larger but diagonal ZA's 
      #  ]
    } else {
      newZAlist <- object$ZAlist
      ranefdata <- object$data
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
      blob <- .make_new_corr_lists(object=object,locdata=ranefdata, which_mats=which_mats, 
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
  } 
  return(RESU)
}

.calc_logphiInfo <- function(invV_factors, which, cum_nobs, phi_est, w.resid, spprec) { ## called by .calc_logdisp_cov_mv(), using the result of .calc_invV.dV_info()
  n_phi <- length(which)
  logphiInfo <- matrix(ncol=n_phi,nrow=n_phi)
  lhs_invV.dVdphi_list <- wresid_list <- list()
  zerovec <- 0*w.resid
  for (colit in seq_along(which)) { 
    mv_it <- which[colit]
    iresp_range <- .subrange(cum_nobs, mv_it)
    indic_it <- zerovec
    indic_it[iresp_range] <- 1  
    i_lhs_invV.dVdphi <- .Dvec_times_Matrix(indic_it, invV_factors$n_x_r) # .fill_rhs_invV.dVdlam(template=zerotemplate, urange=uirange[,ilam], invV.dV_info)
    lhs_invV.dVdphi_list[[as.character(mv_it)]] <- i_lhs_invV.dVdphi
    wresid_list[[as.character(mv_it)]] <- wresid_it <- w.resid*indic_it
    if (spprec) {
      A <- invV_factors$r_x_n %*% invV_factors$n_x_r
      trAB <- sum(A^2)
    } else { # .traceAB() differs by checking which is the largest dimension:
      trAB <- .traceAB(lA=i_lhs_invV.dVdphi, rA=invV_factors$r_x_n, 
                       B_is_tA=TRUE # (lA is n x r, rA is r x n, so A is square as this cases requires)
      ) 
    }    
    # using the pattern (D-nXr.rXn)^2 = D^2 - 2 D nXr.rXn + (nXr.rXn)^2
    logphiInfo[colit,colit] <- sum(wresid_it^2) -2 * .traceDB(wresid_it,i_lhs_invV.dVdphi, invV_factors$r_x_n) + trAB 
    for (coljt in seq_len(colit-1L)) { 
      mv_jt <- which[coljt]
      wresid_jt <- wresid_list[[as.character(mv_jt)]]
      j_lhs_invV.dVdphi <-  lhs_invV.dVdphi_list[[as.character(mv_jt)]] 
      # trace[(D1-l*r) (D2-l*r)], with D1*D2=0:
      logphiInfo[coljt,colit] <- - .traceDB(wresid_it,i_lhs_invV.dVdphi, invV_factors$r_x_n) -
                                   .traceDB(wresid_jt,j_lhs_invV.dVdphi, invV_factors$r_x_n) +
                                   .traceAB(i_lhs_invV.dVdphi, invV_factors$r_x_n, t(invV_factors$r_x_n), t(j_lhs_invV.dVdphi))
      logphiInfo[colit,coljt] <- logphiInfo[coljt,colit]
    }
  }
  sub_phi_vec <- .unlist(phi_est[which]) 
  logphiInfo <- logphiInfo * (sub_phi_vec %*% t(sub_phi_vec)) # inspired from code for lambda 
  return(list(logphiInfo=logphiInfo,lhs_invV.dVdphi_list=lhs_invV.dVdphi_list, wresid_list=wresid_list))
}



.calc_logdisp_cov_mv <- function(object, dvdloglamMat=NULL, dvdlogphiMat=NULL, invV_factors=NULL, force_fixed=FALSE) { 
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
      if ( inherits(strucList[[randit]],"dCHMsimpl")) { # (not expected in default use in mv, for reasons explained in mv, sinc AUG_ZXy presumably FALSE )
        for_dw_i <- as(strucList[[randit]], "sparseMatrix") %*% dvdloglamMat[range_in_dw,  ] # i.e L_Q %*% lignes de (t(L_Q) %*% invG %*% L_Q %*% some rhs) 
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
        for (randjt in which(checklambda)) { ## NOT all ranefs!
          nblocks_randjt <- Xi_cols[randjt]
          cum_colrange_in_dw_i <- (cum_n_u_h[randjt]+1L):(cum_n_u_h[randjt+1L])
          cum_colranges_in_dw_i <- matrix(cum_colrange_in_dw_i,ncol=nblocks_randjt) ## this _splits_ seq(n_u_h[randit]) over two columns for a random-slope model
          for (col_in_colranges_dw_i in nblocks_randjt) { ## half-ranges for random-slope model
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
  is_phiScalS <- (phimodel=="phiScal" | force_fixed) # vector !
  if (any(is_phiScalS) ||  
      identical(object$envir$forcePhiComponent,TRUE) ## hack for code testing: force dispVar computation as if phi was not fixed.
  ) { ## semble impliquer pas outer phi.Fix... => no need to test object$phi.object$phi_outer,"type")
    phi_est <- object$phi ## no need to get object$phi.object$phi_outer
    n_phiScal <- length(which(is_phiScalS))
    for (mv_it in which(is_phiScalS)) if (length(phi_est[[mv_it]])!=1L) problems$stopphi <- warning("phimodel=\"phiScal\" but length(phi_est)!=1L.")
    if ( ! is.null(dvdlogphiMat)) {
      cum_nobs <- attr(object$vec_nobs, "cum_nobs")
      dvdlogphiS <- vector("list", length(is_phiScalS))
      for (mv_it in seq_along(dvdlogphiS)) { ## there are cols even for fixed phiS
        if (is_phiScalS[mv_it]) {
          resp_range <- .subrange(cum_nobs, mv_it)
          dvdlogphiS[[mv_it]] <- rowSums(dvdlogphiMat[,resp_range, drop=FALSE]) ## using each phi_i = phi # always a vector, even from 0-col matrix
        } else dvdlogphiS[[mv_it]] <- NULL # rep(0, nrow(dvdlogphiMat)) ## using each phi_i = phi # always a vector, even from 0-col matrix
      }
      dvdlogphi <- do.call(cbind, dvdlogphiS)
      dwdlogphi <- .calc_invL_coeffs(object,dvdlogphi) # input matrix -> output matrix
      col_info$phi_cols <- length(ranef_ids)+seq_len(n_phiScal) ## cols indices for phi 
      dispcolinfo$logphi <- rep("logphi", n_phiScal)
    } else if (object$models[["eta"]]=="etaHGLM") stop("phimodel=='phiScal' but is.null(dvdlogphiMat)")
  } else {  ## else phimodel="", e.g. binomial
    # if binomial or poisson, phimodel=""; warning for other phimodels
    if (any(phimodel[!is_phiScalS] !="")) {
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
    }) ## col ranges for lambda[checklambda], phi, etc (with wrong names)
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
      logphiInfo_blob <- .calc_logphiInfo(invV_factors, which=which(is_phiScalS), cum_nobs=cum_nobs, 
                                            phi_est=phi_est, w.resid=.get_H_w.resid(object),
                                            spprec=.is_spprec_fit(object))
      logdispInfo[dispcols$logphi,dispcols$logphi] <- logphiInfo_blob$logphiInfo 
      wresid_list <- logphiInfo_blob$wresid_list
      if ("loglambda" %in% dispnames) {
        sublambda <- .unlist(invV.dV_info$lambda)
        logdispInfoBlock <- matrix( nrow=nrand, ncol=length(wresid_list)) # a 
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
            char_mv_jtS <- names(wresid_list)
            for (block_jt in seq_along(char_mv_jtS)) { 
              char_mv_jt <-  char_mv_jtS[block_jt]
              j_lhs_invV.dVdphi <-  logphiInfo_blob$lhs_invV.dVdphi_list[[char_mv_jt]] 
              logdispInfoBlock[colit,block_jt] <- phi_est[[as.integer(char_mv_jt)]] * (
                .traceDB( wresid_list[[char_mv_jt]], invV.dV_info$lhs_invV.dVdlam, i_rhs_invV.dVdlam) -
                .traceAB(invV.dV_info$lhs_invV.dVdlam, i_rhs_invV.dVdlam, j_lhs_invV.dVdphi, invV_factors$r_x_n) )
            }
          }
        }
        logdispInfoBlock <- .Dvec_times_matrix(sublambda, logdispInfoBlock[which(checklambda),,drop=FALSE])  
        logdispInfo[dispcols$loglambda,dispcols$logphi] <- logdispInfoBlock
        logdispInfo[dispcols$logphi,dispcols$loglambda] <- t(logdispInfoBlock)
      }
      if ("rho" %in% dispnames) {
        char_mv_jtS <- names(wresid_list)
        for (jt in seq_along(char_mv_jtS)) { 
          char_mv_jt <-  char_mv_jtS[jt]
          j_lhs_invV.dVdphi <-  logphiInfo_blob$lhs_invV.dVdphi_list[[char_mv_jt]] 
          coljt <- dispcols$logphi[jt] # column of full logdispInfo matrix
          logdispInfo[dispcols$rho,coljt] <- 
            logdispInfo[coljt, dispcols$rho] <- phi_est[[as.integer(char_mv_jt)]] * .traceAB(lhs_invV.dVdrho,rhs_invV.dVdrho, 
                                                                            j_lhs_invV.dVdphi, invV_factors$r_x_n)  
        }
        # phi_est * sum(invV.dVdrho * invV)  
      }
    } 
    logdispInfo <- logdispInfo/2
    resu <- .wrap_solve_logdispinfo(logdispInfo, object)
    if (any( ! checklambda )) { ## if cols missing from logdisp_cov compared to dwdlogdisp
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
    resu$dwdlogdisp <- dwdlogdisp
    return(resu)
    ## more compact than storing ww %*% logdisp_cov %*% t(ww) which is nobs*nobs 
  }
}
