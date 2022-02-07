.calc_new_X_ZAC___2 <- function(object, newdata=NULL, re.form = NULL,
                            variances=list(residVar=FALSE, cov=FALSE),invCov_oldLv_oldLv_list,
                            locform=formula.HLfit(object, which="")) {
  ## possible change of random effect terms
  locform <- .update_formula_shared_ranefs(locform, re.form, rm_RHS=TRUE)
  # checking variables in the data
  if (object$spaMM.version > "3.6.19") {
    frame.form <- .subbarsMM(locform) ## converts (...|...) terms to some "+" form, from a similar concept in lme4 
    environment(frame.form) <- environment(locform)
    allframevars <- as.character(attr(terms(frame.form), "variables"))[-1L]
    # checking variables in the data
    
    { # possible replacement for .get_locdata()
      if (is.vector(newdata)) { ## ## less well controlled case, but useful for maximization
        newdata <- data.frame(matrix(newdata,nrow=1))
        allvarsori <- attr(object$data,"allvars")
        if (length(allvarsori)==ncol(newdata)) {
          names(newdata) <- allvarsori
        } else {
          allvars <- all.vars(.strip_IMRF_args(locform)) ## strip to avoid e.g. 'stuff' being retained as a var from IMRF(..., model=stuff)
          if (variances$residVar) allvars <- unique(c(allvars,all.vars(.strip_IMRF_args(.get_phiform(object)))))  
          stop(paste("(!) newdata has incorrect length. It should match the following variables:\n",paste(allvars,collapse=" ")))
        }
        locframe <- model.frame(data=newdata, formula=frame.form,  drop.unused.levels=TRUE)
      } else if (is.null(newdata)) {
        locframe <- object$data # [ , allframevars,drop=FALSE]  ### [...] a priori suspect 
      } else {
        if( is.matrix(newdata) ) newdata <- as.data.frame(newdata)  
        # so that matrix 'newdata' arguments can be used as in some other predict methods.
        locdata <- try(newdata[ , allvarsori,drop=FALSE]) ## allvars checks only RHS variables
        if (inherits(locdata,"try-error")) stop(paste0("Variable(s) ",
                                                       paste(setdiff(allvars,colnames(newdata)), collapse=","),
                                                       " appear to be missing from newdata."))
        checkNAs <- apply(locdata,1L,anyNA)
        if (any(checkNAs)) {
          message("NA's in required variables from 'newdata'. Prediction not always possible.")
          locdata <- locdata[!checkNAs, ]
        }
        locframe <- model.frame(data=newdata, formula=frame.form,  drop.unused.levels=TRUE)
      }
      # locframe
    }
    locdata <- locframe # for later code that use 'locdata' but slipt this later code in small fns ?
    need_new_design <- ( ( ! is.null(newdata) ) || ! is.null(re.form)) ## newdata or new model
    if (need_new_design) {
      { # replacement for .get_newX_info()
        
        # not clear what to do since it calls .calc_newFrames_fixed() -> model.frame() so it needs raw 'data"

      }
      RESU$locdata <- locframe
    } else RESU <- list(locdata=locframe,newX.pv=object$X.pv)
    
  } else {
    allvars <- all.vars(.strip_IMRF_args(locform)) ## strip to avoid e.g. 'stuff' being retained as a var from IMRF(..., model=stuff)
    if (variances$residVar) allvars <- unique(c(allvars,all.vars(.strip_IMRF_args(.get_phiform(object)))))  
    #
    locdata <- .get_locdata(newdata, locvars=allvars, object=object, variances=variances) # always required (cf when newdata is not a data frame)
    need_new_design <- ( ( ! is.null(newdata) ) || ! is.null(re.form)) ## newdata or new model
    if (need_new_design) {
      RESU <- .get_newX_info(locform, locdata, object)
      RESU$locdata <- locdata
    } else RESU <- list(locdata=locdata,newX.pv=object$X.pv) 
  }
  #
  #
  ## subZAlist is a subset of the old ZA, newZAlist contains new ZA ; different uses=> computation required under disticnt conditions for each.
  ## calling .make_corr_list(object,...) is always OK bc the fist argument may be a superset of the required list
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
      if (need_new_design) {
        ## with newdata we need Evar and then we need nn... if newdata=ori data the Evar (computed with the proper nn) should be 0
        #barlist <- .process_bars(locform,as_character=FALSE) 
        new_exp_ranef_terms <- structure(ori_exp_ranef_terms[newinold], type=attr(ori_exp_ranef_terms,"type")[newinold])
        ranef_form <- as.formula(paste("~",(paste(new_exp_ranef_strings,collapse="+")))) ## effective '.noFixef'
        newZlist <- .calc_Zlist(exp_ranef_terms=new_exp_ranef_terms, # .process_bars(barlist=barlist,as_character=FALSE, which.="exp_ranef_terms"), # != barlist, for IMRF notably
                                #locform, 
                                data=locdata, rmInt=0L, drop=TRUE,sparse_precision=FALSE, 
                                levels_type= "seq_len", ## superseded in specific cases: notably, 
                                ## the same type has to be used by .calc_AMatrix_IMRF() -> .as_factor() 
                                ##  as by .calc_Zmatrix() -> .as_factor() for IMRFs.
                                ## This is controlled by option uGeo_levels_type (default = "data_order" as the most explicit).
                                ## The sames functions are called with the same arguments for predict with newdata.
                                ZAlist_info=object$ZAlist[newinold],  
                                lcrandfamfam=attr(object$rand.families,"lcrandfamfam")) 
        amatrices <- .get_new_AMatrices(object, newdata=locdata) # .calc_newFrames_ranef(formula=ranef_form,data=locdata,fitobject=object)$mf)
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
