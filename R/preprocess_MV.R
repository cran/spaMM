## The file is in .Rbuildignore

# for structured list of parameters, fully named by type of param as assumed by as assumed by .rename_ranPars()
.merge_mv_parlist <- function(parlistS, merged, ZAlist=merged$ZAlist) {
  resu <- list()
  map_rd_mv <- attr(ZAlist, "map_rd_mv")
  nrand <- length(ZAlist)
  for (mv_it in seq_along(map_rd_mv)) {
    rd_in_mv <- map_rd_mv[[mv_it]]
    stuff_it <- .rename_ranPars(parlistS[[mv_it]], rd_in_mv, mv_it)
    resu <- .modify_list(resu,stuff_it) 
  }
  resu
}

# here mv_stuff is a list of the same length as map_rd_mv, whose elements have no particular structure and names
.merge_mv_list <- function(mv_stuff, merged, ZAlist=merged$ZAlist, full=FALSE) {
  nrand <- length(ZAlist)
  if (full) {
    resu <- structure(rep(list(NULL), nrand), names=seq_len(nrand)) 
  } else resu <- list() # and explicit NULLs in 'mv_stuff' are not passed to 'resu' 
  map_rd_mv <- attr(ZAlist, "map_rd_mv")
  for (mv_it in seq_along(map_rd_mv)) {
    rd_in_mv <- map_rd_mv[[mv_it]]
    if (length(stuff_it <- mv_stuff[[mv_it]])) {
      if ( is.null(nams <- names(stuff_it))) names(stuff_it) <- nams <- seq_along(stuff_it)
      names(stuff_it) <- rd_in_mv[nams]
      # resu <- .modify_list(resu,stuff_it) # that's not good bc explicit NULLs in 'stuff_it' replace corresponding elements in 'resu'
      # which may be expected behavious in other contexts of use of .modify_list() (incl. recursive ones), but not here.
      for (st in names(stuff_it)) if ( ! is.null(stuff_it[[st]])) resu[[st]] <- stuff_it[[st]]
    }
  }
  resu
}


.map_rd_mv <- function(ZAlist, unmerged) {
  nmodels <- length(unmerged)
  map_rd_mv <- vector("list", nmodels)
  for (mv_it in seq_len(nmodels)) {
    terms_it <- attr(unmerged[[mv_it]]$ZAlist,"exp_ranef_strings") 
    map_rd_mv[[mv_it]] <- structure(match(terms_it,attr(ZAlist,"exp_ranef_strings")), names=seq_along(terms_it))
  }
  return(map_rd_mv) # contains the rank of each rd of the total model in the rds of each submodel
}

.merge_rand_families <- function(unmerged, map_rd_mv=attr(ZAlist, "map_rd_mv"), nrand=length(ZAlist), ZAlist) {
  nmodels <- length(unmerged)
  rand_families <- vector("list", nrand)
  for (rd in seq_len(nrand)) {
    for (mv_it in seq_len(nmodels)) {
      if (length(rd_in_mv <- which(map_rd_mv[[mv_it]]==rd))) { # gets the position in map_rd_mv[[mv_it]] 
        # and the name associated to this position, which the position in unmerged[[mv_it]][["rand.families"]]:  
        #rd_in_submodel <- names(rd_in_mv)
        newfam <- unmerged[[mv_it]][["rand.families"]][[names(rd_in_mv)]]
        if ( is.null(rand_families[[rd]])) {
          rand_families[[rd]] <- newfam
        } else if ( ! identical(rand_families[[rd]],newfam, ignore.environment=TRUE)) { 
          stop("Conflicting families for random effect term shared accross models.")
        }
      }
    }
  }
  rand_families <- .checkRandLinkS(rand_families)
  return(rand_families)
}

.merge_Zlists <- function(ZAlist1, ZAlist2, nobs1, nobs2, ranefs1, ranefs2, mv_it) {
  if ( ! length(ZAlist1)) ranefs1 <- c() # here we need to handle case where first list is empty 
  if ( ! length(ZAlist2)) ranefs2 <- c() # here rather the case where second list is NULL
  Zlist <- .merge_ZZlists(ZAlist1, ZAlist2, nobs1, nobs2, 
                          ranefs1=ranefs1,
                          ranefs2=ranefs2,
                          mv_it, type="Zlist")  
  names(Zlist) <- unique(c(ranefs1,ranefs2))
  Zlist
}


.merge_ZAlists <- function(ZAlist1, ZAlist2, nobs1, nobs2, mv_it) {
  .merge_ZZlists(ZAlist1, ZAlist2, nobs1, nobs2, 
                 ranefs1=attr(ZAlist1,"exp_ranef_strings"),
                 ranefs2=attr(ZAlist2,"exp_ranef_strings"),
                 mv_it, type="ZAlist")  
}

.merge_ZZlists <- function(ZAlist1, ZAlist2, nobs1, nobs2, 
                           ranefs1,
                           ranefs2,
                           mv_it, type="ZAlist") {
  allranefs <- unique(c(ranefs1,ranefs2))
  ZAlist <- exp_ranef_terms <- namesTerms <- vector("list", length(allranefs))
  exp_ranef_types <- exp_ranef_strings <- type_attr <- character(length(allranefs))
  Xi_cols <- integer(length(allranefs))
  # exp_ranef_strings and exp_ranef_terms will have a global 'type' attibute
  for (ran_it in seq_along(allranefs)) {
    ranf <- allranefs[ran_it]
    in1 <- which(ranefs1==ranf)
    in2 <- which(ranefs2==ranf)
    if (length(in1) && length(in2)) { # true merge
      Zlistori <- ZAlist1
      ZA1 <- ZAlist1[[in1]]
      which_mv <- c(attr(ZA1,"which_mv"), mv_it)
      namesTerm <- attr(ZA1,"namesTerm")
      ZA2 <- ZAlist2[[in2]]
      if (is.null(ii1 <- attr(ZA1,"is_incid")) || is.null(ii2 <- attr(ZA2,"is_incid"))) {
        is_incid <- NULL # they can be NULL by design, cf .calc_ZAlist() when there is an A matrix
      } else is_incid <- ii1 && ii2
      # levels of the LHS of the term:
      LHS_levels <- unique(c(attr(ZA1,"LHS_levels"),attr(ZA2,"LHS_levels")))
      # levels of the RHS...
      ulevels1 <- unique(colnames(ZA1))
      ulevels2 <- unique(colnames(ZA2))
      if ( ! identical(ulevels1,ulevels2)) { # add columns (here) before rbinding rows (later)
        alllevels <- unique(c(ulevels1,ulevels2))
        extralevels1 <- setdiff(alllevels,ulevels1)
        #
        if ((Xi_ncol <- attr(ZAlist1,"Xi_cols")[in1])>1L) {
          #allsortedlevels <- as.character(sort(as.integer(alllevels))) # ugly but clear
          nlev1 <- length(ulevels1)
          # assuming some kind of ordering with each original matrix, we don't reorder after naming
          ZA1_colblocks <- vector("list",Xi_ncol)
          for (it in seq_len(Xi_ncol)) {
            ZA1_colblock <-  .adhoc_cbind_dgC_0(ZA1[,(it-1L)*nlev1+seq_len(nlev1),drop=FALSE],
                                                length(extralevels1))
            colnames(ZA1_colblock) <- alllevels
            ZA1_colblocks[[it]] <- ZA1_colblock #[,allsortedlevels]
          }
          ZA1 <- do.call(cbind,ZA1_colblocks)
        } else {
          ZA1 <- .adhoc_cbind_dgC_0(ZA1,length(extralevels1))
          colnames(ZA1) <- alllevels
        }
        #
        extralevels2 <- setdiff(alllevels,ulevels2)
        nextcolnames <- c(ulevels2,extralevels2) # so not the same order as alllevels: reordering after naming
        if ((Xi_ncol <- attr(ZAlist2,"Xi_cols")[in2])>1L) {
          nlev2 <- length(ulevels2)
          ZA2_colblocks <- vector("list",Xi_ncol)
          for (it in seq_len(Xi_ncol)) {
            ZA2_colblock <-  .adhoc_cbind_dgC_0(ZA2[,(it-1L)*nlev2+seq_len(nlev2),drop=FALSE],
                                                length(extralevels2))
            colnames(ZA2_colblock) <- nextcolnames
            #ZA2_colblocks[[it]] <- ZA2_colblock[,allsortedlevels]
            ZA2_colblocks[[it]] <- ZA2_colblock[,alllevels]
          }
          ZA2 <- do.call(cbind,ZA2_colblocks)
        } else {
          ZA2 <- .adhoc_cbind_dgC_0(ZA2,length(extralevels2))
          colnames(ZA2) <- nextcolnames
          ZA2 <- ZA2[,alllevels]
        }
      }
      ori <- in1
    } else  { 
      # Identify source of attributes:
      if (length(in1)) {
        Zlistori <- ZAlist1
        ori <- in1
        which_mv <- attr(Zlistori[[ori]],"which_mv")
        ZA1 <- ZAlist1[[in1]]
        namesTerm <- attr(ZA1,"namesTerm")
        LHS_levels <- attr(ZA1,"LHS_levels")
        is_incid <- attr(ZA1,"is_incid") ## will surely need to be tested ___FIXME___ vs <- FALSE
        ZA2 <- Matrix(0,ncol=ncol(ZA1), nrow=nobs2)
      } else {
        Zlistori <- ZAlist2
        ori <- in2
        which_mv <- mv_it
        ZA2 <- ZAlist2[[in2]]
        namesTerm <- attr(ZA2,"namesTerm")
        LHS_levels <- attr(ZA2,"LHS_levels")
        is_incid <- attr(ZA2,"is_incid") ## will surely need to be tested ___FIXME___ vs <- FALSE
        ZA1 <- Matrix(0,ncol=ncol(ZA2), nrow=nobs1)
      } 
    } 
    ZAlist[[ran_it]] <- structure(rbind(ZA1,ZA2), which_mv=which_mv, namesTerm=namesTerm, is_incid=is_incid, LHS_levels=LHS_levels)
    # ... leftOfBar_mf" "dataordered_unique_levels" "AR1_block_n_u_h_s" "uniqueGeo"
    namesTerms[ran_it] <- attr(Zlistori,"namesTerms")[ori] # namesTerms is a *named* list... but pathetically that does not copy the name
    names(namesTerms)[ran_it] <- names(attr(Zlistori,"namesTerms")[ori])
    type_attr[ran_it] <- attr(attr(Zlistori,"exp_ranef_terms"),"type")[ori]
    exp_ranef_terms[[ran_it]] <- attr(Zlistori,"exp_ranef_terms")[[ori]]
    exp_ranef_types[ran_it] <- attr(Zlistori,"exp_ranef_types")[ori]
    if (type=="ZAlist") exp_ranef_strings[ran_it] <- attr(Zlistori,"exp_ranef_strings")[ori]
    Xi_cols[ran_it] <- attr(Zlistori,"Xi_cols")[ori]
  }
  return(structure(ZAlist, 
                   namesTerms=namesTerms, 
                   exp_ranef_terms=structure(exp_ranef_terms,type=type_attr),
                   exp_ranef_types=exp_ranef_types,
                   exp_ranef_strings=structure(exp_ranef_strings,type=type_attr),
                   Xi_cols=Xi_cols
  ))
} 

.merge_Xs <- function(X1,X2, mv_it, REML=FALSE) {
  if ( ! is.null(colnames(X2))) colnames(X2) <- paste0(colnames(X2),"_",mv_it)
  if (is.null(X1)) { 
    X <- X2
  } else {
    XX1 <- cbind(X1,matrix(0, nrow=nrow(X1), ncol=ncol(X2)))
    colnames(XX1) <- c(colnames(X1),colnames(X2))
    XX2 <- cbind(matrix(0, nrow=nrow(X2), ncol=ncol(X1)),X2)
    X <- rbind(XX1,XX2)
    attr(X,"scaled:scale") <- c(attr(X1,"scaled:scale"),attr(X2,"scaled:scale"))
    # we ignore the 'assign' attribute which is used only for preprocessing => .determine_sparse_X_mv() will directly use those from 'unmerged'.
  }
  #
  if (REML) {
    if ( ! is.null(attr(X2, "extra_vars"))) {
      attr(X,"extra_vars") <- c(attr(X1, "extra_vars"), paste0(attr(X2, "extra_vars"),"_",mv_it))
    } else attr(X,"extra_vars") <- attr(X1, "extra_vars")
  }
  return(X)
}

.rename_ranPars <- function(parlist, rd_in_mv, mv_it, hy_in_mv=NULL) {
  if ( ! is.null(parlist$lambda)) names(parlist$lambda) <- rd_in_mv[names(parlist$lambda)]
  if ( ! is.null(parlist$trLambda)) names(parlist$trLambda) <- rd_in_mv[names(parlist$trLambda)]
  if ( ! is.null(parlist$corrPars)) names(parlist$corrPars) <- rd_in_mv[names(parlist$corrPars)]
  if ( ! is.null(parlist$ranCoefs)) names(parlist$ranCoefs) <- rd_in_mv[names(parlist$ranCoefs)]
  if ( ! is.null(parlist$trRanCoefs)) names(parlist$trRanCoefs) <- rd_in_mv[names(parlist$trRanCoefs)]
  
  if ( ! is.null(parlist$COMP_nu)) names(parlist$COMP_nu) <- mv_it
  if ( ! is.null(parlist$NB_shape)) names(parlist$NB_shape) <- mv_it
  if ( ! is.null(parlist$trNB_shape)) names(parlist$trNB_shape) <- mv_it
  if ( ! is.null(parlist$phi)) names(parlist$phi) <- mv_it #parlist$phi <- list(parlist$phi, names=mv_it)
  if ( ! is.null(parlist$trPhi)) names(parlist$trPhi) <- mv_it #parlist$trPhi <- list(parlist$trPhi, names=mv_it)
  
  if ( ! is.null(parlist$hyper)) names(parlist$hyper) <- hy_in_mv[names(parlist$hyper)]
  return(parlist)
} 

.merge_lambdas_mv <- function(lambda_merger, optim_blob) {
  for (char_rd in names(lambda_merger)) { 
    if ((len_lam <- length(lambdas <- lambda_merger[[char_rd]]))>1L) {
      lambda_merger[[char_rd]] <- exp(mean(log( lambda_merger[[char_rd]])))
    } # else if (len_lam==0L) lambda_merger[char_rd] <- NULL
  }
  lambda_merger <- unlist(lambda_merger)
  lambda_merger <- list(inits=list(init=list(lambda=lambda_merger),
                                   init.optim=list(trLambda=.dispFn(lambda_merger))))
  .modify_list(optim_blob,lambda_merger, obey_NULLs=FALSE) 
}

.calc_optim_args_mv <- function(processed, map_rd_mv, user_init_optim, fixedS, user.lower, user.upper, verbose, optim.scale) {
  unmerged <- processed$unmerged 
  optim_blob <-  NULL
  hyper_info <- processed$hyper_info
  if (has_hy <- (length(hyper_info$template)>0L)) {map_hy_mv <- hyper_info$map_hy_mv}  
  nrand <- length(processed$ZAlist)
  lambda_merger <- structure(vector("list", nrand), names=seq_len(nrand))
  for (mv_it in seq_along(unmerged)) {
    rd_in_mv <- map_rd_mv[[mv_it]]
    # operation requiring full model indices:
    rd_in_submv <- structure(seq_along(rd_in_mv), names=rd_in_mv) # reverse map (its names are full model indices)
    phistr <- structure(1,names=mv_it)
    skel <- list(phi=phistr, # ugly per se but fits in the general algo
                 ranCoefs=rd_in_submv,   
                 trRanCoefs=rd_in_submv,
                 lambda=rd_in_submv,   
                 trLambda=rd_in_submv,
                 corrPars=rd_in_submv,# fake skeleton (not sublist here) but appears sufficient # contains kappa...
                 NB_shape=phistr, 
                 COMP_nu=phistr)
    # 'hyper' elements are indexed in not by ranefs but only for reference by the map, which creates a new problem
    if (has_hy && length(hy_in_mv <- map_hy_mv[[mv_it]])) {
      hy_in_submv <- structure(names(hy_in_mv), names=hy_in_mv) # reverse map (its names are full model indices); names() and seq_along should be equiv ?
      hyper_skel <- list(hyper=hy_in_submv)
      skel <- c(skel, hyper_skel)
    } else hy_in_mv <- NULL
    # the global 'init' is used to complete the submodel inits
    # But no such thing for the lower, upper, which are used "globally" in an additional call to .makeLowerUpper()
    # This is a bit tortuous and the Q is whether we could perform only a final .makeLowerUpper() call  (__FIXME__)
    user_init_optim_it <- .subPars(user_init_optim,skeleton = skel)
    # : works on phi only if phi is a named list. result is a named sub*list*
    user_init_optim_it <- .rename_ranPars(user_init_optim_it, rd_in_submv, 1L, hy_in_mv=hy_in_mv)
    # we want user_init_optim_it[["phi"]] to be an unnamed scalar, or NULL, starting from any user input (input unnamed vector got names by .reformat_phi()) 
    # user_init_optim_it[["phi"]] <- user_init_optim[["phi"]][[as.character(mv_it)]] # OK for phi list or complete vector, but not incomplete vector 
    user_init_optim_it[["phi"]] <- as.vector(.unlist(user_init_optim_it[["phi"]])) 
    #
    user_lower_it <- .subPars(user.lower,skeleton = skel) # 
    user_upper_it <- .subPars(user.upper,skeleton = skel)
    fixed_it <- .subPars(fixedS,skeleton = skel) # for mv_it=2 we potentially have $phi =c("2"=.).... 
    fixed_it <- .rename_ranPars(fixed_it,rd_in_submv, 1L, hy_in_mv=hy_in_mv) # and then $phi =c("1"=.) bc of the 3rd argument of .rename_ranPars 
    # first argument list(phi=.) bc list(phi=NULL) is modified as expected while .modify_list(list(NULL), .) returns NULL.
    # Not so necessary for lambda as lambda.Fix is never NULL; 
    #   but the fact that it works for unnamed list list(unmerged[[mv_it]][["lambda.Fix"]]) is perhaps best not to be assumed.
    unmerged[[mv_it]][["phi.Fix"]] <- .modify_list(list(phi=unmerged[[mv_it]][["phi.Fix"]]), fixed_it["phi"])[[1]]
    unmerged[[mv_it]][["lambda.Fix"]] <- .modify_list(list(lambda=unmerged[[mv_it]][["lambda.Fix"]]), fixed_it["lambda"])[[1]] # bc : .../...
    # .calc_optim_args() -> .more_init_optim() distinguishes 'NA' in lambda.Fix (some provided by .preprocess()) versus NULL
    #          hence complex above call in order not to erase such NA by a NULL...
    # also used in HLfit_body (and equivalent fns) -> .calc_initial_init_lambda() ;
    # and .calc_optim_args() -> {init.optim <- .more_init_optim(proc1=proc1, corr_types=corr_types, init.optim=init.optim)}  ./.
    #  ./. -> .init_optim_lambda_ranCoefs() has to check not globally fixed lambdas ;
    # Maybe not optimal but need to distinguish globally fixed param at some point? (_F I X M E_) 
    optim_blob_it <- .calc_optim_args(proc_it=unmerged[[mv_it]], processed=processed,
                                      user_init_optim=user_init_optim_it, fixed=fixed_it, lower=user_lower_it, upper=user_upper_it, 
                                      verbose=verbose, optim.scale=optim.scale, For="fitmv") 
    # operation requiring full model indices:
    # .modify_list() requires that all vectors are named; and moreover, that ranef indices are those of the mv model, not of each sub-model
    # which(map_rd_mv[[mv_it]])
    optim_blob_it$inits[["init"]] <- .rename_ranPars(optim_blob_it$inits[["init"]], rd_in_mv, mv_it, hy_in_mv=hy_in_mv)
    optim_blob_it$inits[["init.optim"]] <- .rename_ranPars(optim_blob_it$inits[["init.optim"]], rd_in_mv, mv_it, hy_in_mv=hy_in_mv)
    optim_blob_it$inits[["init.HLfit"]] <- .rename_ranPars(optim_blob_it$inits[["init.HLfit"]], rd_in_mv, mv_it)
    optim_blob_it$inits[["ranFix"]] <- .rename_ranPars(optim_blob_it$inits[["ranFix"]], rd_in_mv, mv_it, hy_in_mv=hy_in_mv)
    optim_blob_it[["fixed"]] <- .rename_ranPars(optim_blob_it[["fixed"]], rd_in_mv, mv_it, hy_in_mv=hy_in_mv)
    if (length(rd_in_mv)) { # if ranef(s) in sub-model...
      names(optim_blob_it[["corr_types"]]) <- rd_in_mv # full length vector
      # corrMatrix argument no longer properly handled if this is removed...
    }
    if (FALSE) {
      optim_blob_it$LUarglist[["canon.init"]] <- .rename_ranPars(optim_blob_it$LUarglist[["canon.init"]], rd_in_mv, mv_it, hy_in_mv=hy_in_mv)
      optim_blob_it$LUarglist[["init.optim"]] <- .rename_ranPars(optim_blob_it$LUarglist[["init.optim"]], rd_in_mv, mv_it, hy_in_mv=hy_in_mv)
      optim_blob_it$LowUp[["lower"]] <- .rename_ranPars(optim_blob_it$LowUp[["lower"]], rd_in_mv, mv_it, hy_in_mv=hy_in_mv)
      optim_blob_it$LowUp[["upper"]] <- .rename_ranPars(optim_blob_it$LowUp[["upper"]], rd_in_mv, mv_it, hy_in_mv=hy_in_mv)
      if (length(rd_in_mv)) { # if ranef(s) in sub-model...
        names(optim_blob_it$LUarglist[["moreargs"]]) <- rd_in_mv[names(optim_blob_it$LUarglist[["moreargs"]])]
      }
    }
    lambdas <- optim_blob_it$inits$init$lambda
    for (char_rd in names(lambdas)) lambda_merger[[char_rd]] <- c(lambda_merger[[char_rd]], lambdas[[char_rd]]) 
    optim_blob <- .modify_list(optim_blob,optim_blob_it, obey_NULLs=FALSE) 
  } 
  if ( ! is.null(optim_blob$inits$init$lambda)) {
    optim_blob <- .merge_lambdas_mv(lambda_merger, optim_blob)
  }
  
  if ( ! is.null(optim_blob$inits$init.optim$hyper)) {
    attr(optim_blob$inits$init.optim$hyper,"hy_info") <- processed$hyper_info # *environment*: for .makeLowerUpper, with distinct name for easier tracking 
  }
  #optim_blob$LUarglist$corr_types <- optim_blob$corr_types # necessary since optim_blob_it$LUarglist$corr_types has no names to allow merging
  #                                                          in contrast to optim_blob_it$corr_types which just got names for that purpose.
  # (but $LUarglist$moreargs is still defective)
  optim_blob <- .makeLowUp_stuff_mv(optim_blob, user.lower=user.lower, user.upper=user.upper, optim.scale, processed, verbose)
  optim_blob
}

.process_bars_mv <- function(predictors, map_rd_mv, as_character=FALSE) {
  exp_barlist <- NULL
  for (mv_it in seq_along(predictors)) {
    exp_barlist_it <- .process_bars(predictors[[mv_it]],as_character=as_character)
    if (length(exp_barlist_it)) {
      type_it <- attr(exp_barlist_it, 'type')
      names(exp_barlist_it) <- map_rd_mv[[mv_it]][names(exp_barlist_it)]
      names(type_it) <- map_rd_mv[[mv_it]][names(type_it)]
      exp_barlist <- structure(.modify_list(exp_barlist, exp_barlist_it),
                               type=.modify_list(attr(exp_barlist, 'type'), type_it))
    }
  }
  exp_barlist
}

.determine_sparse_X_mv <- function(terms_info, X.pv, vec_nobs, unmerged) {
  sparse_X <- spaMM.getOption("sparse_X") 
  ## forcing sparse_X may (1) be slow for small problems 
  ## (2) entails the use of Matrix::Cholesky, which is less accurate => small bu visible effect on predVar in singular 'twolambda' case
  if (is.null(sparse_X)) {
    nmodels <- length(terms_info$Y)
    col_heuristic_densenesseS <- vector("list", nmodels)
    rel_nobs <- vec_nobs/sum(vec_nobs)
    for (mv_it in seq_along(nmodels)) {
      asgn <- attr(unmerged[[mv_it]][["AUGI0_ZX"]]$X.pv,"assign") ## "for each column in the matrix ... the term in the formula which gave rise to the column"
      col_heuristic_denseness_it <- rep(rel_nobs[mv_it],length(asgn))
      if ( length(fixef_levels <- .get_from_terms_info(terms_info=terms_info, which="fixef_levels", mv_it=mv_it)) ) {
        terms_densenesses_it <- rel_nobs[mv_it] * 
          .calc_terms_heuristic_denseness(terms_info,
                                          fixef_off_terms=.get_from_terms_info(terms_info=terms_info, which="fixef_off_terms", mv_it=mv_it),
                                          fixef_levels=fixef_levels) 
        for (jt in seq_along(asgn)) if (asgn[jt]>0L) col_heuristic_denseness_it[jt] <- terms_densenesses_it[asgn[jt]]
      }
      col_heuristic_densenesseS[[mv_it]] <- col_heuristic_denseness_it
    }
    col_heuristic_denseness <- .unlist(col_heuristic_densenesseS)
    sparse_X <- (mean(col_heuristic_denseness)<0.11) ## FIXME not enough tests of threshold; could use data.test in test-predVar which has mean density=0.19
  }
  return(sparse_X)
}

.correct_ZA_mv_ranCoefs <- function(ZAlist, mv_it) { # for the ZAlist of aubmodel mv_it...
  # the mv_it argument serves to avoid unnecessary operations, but it seems the block depending on it could be run in all cases
  # provided paste0(".mv",mv_it) is replaced by paste0(".mv",model_ids)
  # 
  # If, by oversight, the mv(...mv_it...) term is missing from submodel mv_it, 
  # there is no ZA matrix for this term in the ZAlist for this submodel, so there is a virtual matrix of zero
  # which does not need to be corrected here. So the actually fitted model is a meaningful interpretation of the model formulas.
  # But .check_mv_in_submodels() will warn about the issue.
  Xi_cols <- attr(ZAlist,'Xi_cols')
  exp_ranef_strings <- attr(ZAlist,"exp_ranef_strings")
  for(rd in seq_along(ZAlist)) {
    LHS_levels <- attr(ZAlist[[rd]],"LHS_levels") # "2", "3' for all relevant matrices if mv(2,3) or 0+(mv(2,3)) 
    if ( ! is.null(model_ids <- LHS_levels[[".mv"]])) { # "mv("-specific code 
      Xi_ncol <- Xi_cols[rd]
      if (mv_it %in% model_ids) { 
        ZA <- ZAlist[[rd]]
        n_levels <- ncol(ZA)/Xi_ncol
        ZAattr <- attributes(ZA)
        namesTerm <- attr(ZA,"namesTerm") # "(Intercept)" ".mv2" or ".mv1" ".mv2" for all relevant matrices depending on absence/presence of 0+
        ZA <- .Matrix_times_Dvec(ZA, rep(as.numeric(namesTerm %in% c("(Intercept)",paste0(".mv",mv_it))),
                                         rep(n_levels,Xi_ncol)))
        ZA <- drop0(ZA)
        attr(ZA,"is_incid") <- ! ("(Intercept)" %in% namesTerm) # while is_incid was FALSE for the template...
        names_lostattrs <- setdiff(names(ZAattr), names(attributes(ZA)))
        attributes(ZA)[names_lostattrs] <- ZAattr[names_lostattrs] 
        ZAlist[[rd]] <- ZA
      }
    }
  }
  return(ZAlist)
}

.correct_newZAlist_mv_ranCoefs <- function(ZAlist, Xi_cols=attr(ZAlist,'Xi_cols'), cum_nobs) {
  for(rd in seq_along(ZAlist)) {
    LHS_levels <- attr(ZAlist[[rd]],"LHS_levels") # "2", "3' for all relevant matrices if mv(2,3) or 0+(mv(2,3)) 
    if ( ! is.null(model_ids <- LHS_levels[[".mv"]])) { # "mv("-specific code 
      Xi_ncol <- Xi_cols[rd]
      ZA <- ZAlist[[rd]]
      n_levels <- ncol(ZA)/Xi_ncol
      ZAattr <- attributes(ZA)
      namesTerm <- attr(ZA,"namesTerm") # "(Intercept)" ".mv2" or ".mv1" ".mv2" for all relevant matrices depending on absence/presence of 0+
      for (mv_it in as.integer(model_ids)) { 
        obs.range <- (cum_nobs[mv_it]+1L):cum_nobs[mv_it+1L]
        ZA[obs.range,] <- .Matrix_times_Dvec(ZA[obs.range,], rep(as.numeric(namesTerm %in% c("(Intercept)",paste0(".mv",mv_it))),
                                                                 rep(n_levels,Xi_ncol)))
      }
      ZA <- drop0(ZA)
      attr(ZA,"is_incid") <- ! ("(Intercept)" %in% namesTerm) # while is_incid was FALSE for the template...
      names_lostattrs <- setdiff(names(ZAattr), names(attributes(ZA)))
      attributes(ZA)[names_lostattrs] <- ZAattr[names_lostattrs] 
      ZAlist[[rd]] <- ZA
    }
  }
  ZAlist
}


.check_identifiability_LMM_mv <- function(processed, vec_nobs, map_rd_mv, unmerged=processed$unmerged) {
  ## identifiability checks cf modular.R -> checkNlevels() in lmer:
  vec_n_u_h <- diff(processed$cum_n_u_h)
  if (any(vec_n_u_h<2L)) {
    problems <- which(vec_n_u_h<2L) 
    for (rd in problems) {
      in_submodel <- logical(length(map_rd_mv))
      for (mv_it in seq_along(unmerged)) in_submodel[mv_it] <- rd %in% map_rd_mv[[mv_it]]
      which_submodels <- which(in_submodel)
      LMMbools <- unlist(lapply(unmerged[which_submodels], function(v) attr(v[["models"]], "LMMbool")))
      if (all(LMMbools) && ! length(.unlist(processed$phi.Fix[which_submodels]))) {
        mess <- paste0("Only ",vec_n_u_h[rd]," level for random effect ",
                       attr(processed$ZAlist,"exp_ranef_strings")[rd],
                       ";\n   this model cannot be fitted unless phi is fixed in a submodel where the random effect appears.")
        warning(mess, immediate.=TRUE)
      }
    }
  }
  for (rd in seq_along(vec_n_u_h)) {
    in_submodel <- logical(length(map_rd_mv))
    for (mv_it in seq_along(unmerged)) in_submodel[mv_it] <- rd %in% map_rd_mv[[mv_it]]
    which_submodels <- which(in_submodel)
    if (length(which_submodels)==1L && attr(unmerged[[which_submodels]][["models"]],"LMMbool") ) {
      if (vec_n_u_h[rd]==vec_nobs[which_submodels] && processed$models[["phi"]][which_submodels] %in% c("phiScal","phiGLM")) { 
        if (attr(processed$residModels[[which_submodels]]$formula,"has_intercept")!=0L) { ## there is an intercept in the resid.model formula
          # cf comments in univariate version
          term_ranef <- attr(processed$ZAlist,"exp_ranef_strings")[rd]
          if (substr(term_ranef, 1, 1)=="(" ## excludes spatial (and more generally 'keyword') ranefs 
              && ! is.numeric(processed$lambda.Fix[rd])
          ) {
            mess <- paste0("Unable to ascertain full-model identifiablity from information for submodel ",which_submodels,
                           " alone,\n   where number of levels = number of observations for random effect ", term_ranef,
                           ";\n   Full model might not be identifiable unless phi is fixed",
                           ",\n   or the variance of this effect is fixed, or a non-trivial correlation matrix is given.") 
            message(mess) # warning bc not sure of correct detection; would be stop() otherwise
          }          
        }
      }
    }
  }
}

.check_mv_in_submodels <- function(ZAlist) {
  # Check that he mv() terms are where they are expected
  # each ZA matrix has the info attr(., "which_mv") which tells in which submodel a ranef was actually found;
  # The problem is recovering the mv() info ! The "LHS_levels" attribute is not unique to mv() factors
  # and the "namesTerm" attribute is ambiguous (cf "(Intercept)")
  exp_ranef_terms <- attr(ZAlist, "exp_ranef_terms")
  for (rd in seq_along(exp_ranef_terms)) {
    lhs <- .DEPARSE(exp_ranef_terms[[rd]][[2]])
    if (grepl("mv(",lhs, fixed=TRUE)) {
      model_ids <- sub("(mv)(\\([^|]+)","c\\2", lhs)
      model_ids <- eval(parse(text=model_ids))
      which_mv <- attr(ZAlist[[rd]], "which_mv")
      if ( ! setequal(model_ids, which_mv)) {
        warnmess <- paste("Random effect term", attr(ZAlist, "exp_ranef_strings")[[rd]], "expected in submodels",
                          paste(model_ids,collapse=","), "but instead found in submodel(s)", paste(which_mv,collapse=",") )
        warning(warnmess, immediate.=TRUE)
      }
    }
  }
}
#

.merge_hyper_infos <- function(ZAlist, unmerged) {
  nrand <- length(ZAlist)
  merged_map <- structure(rep(NA, nrand), names=seq_len(nrand)) 
  map_rd_mv <- attr(ZAlist, "map_rd_mv")
  idx <- 0L
  for (mv_it in seq_along(map_rd_mv)) {
    info_it <- unmerged[[mv_it]][["hyper_info"]]
    rd_in_mv <- map_rd_mv[[mv_it]]
    stuff_it <- na.omit(info_it$map)
    if (length(stuff_it)) {
      hy_idxes <- names(info_it$ranges)
      for (hy_it in hy_idxes) {
        idx <- idx+1L
        rd_in_submv <- names(info_it$ranges[[hy_it]])
        rd_in_mv[rd_in_submv]
        merged_map[rd_in_mv[rd_in_submv]] <- idx
      }
    }
  }
  umap <- unique(na.omit(merged_map))
  for (idx in seq_along(umap)) merged_map[merged_map==umap[idx]] <- idx
  umap <- unique(na.omit(merged_map))
  merged_ranges <- template <- structure(vector("list", length(umap)), names=seq_along(umap)) 
  for (idx in seq_along(umap)) {
    which_rds <- names(which(merged_map==idx))
    merged_ranges[[idx]] <- structure(as.integer(which_rds), names=which_rds)
  }
  summingMat <- .calc_summingMat_hyper(nrand, merged_map, merged_ranges)
  map_hy_mv <- vector("list", length(map_rd_mv))
  for (mv_it in seq_along(map_rd_mv)) {
    rd_in_mv <- map_rd_mv[[mv_it]]
    inverse_map <- structure(names(rd_in_mv), names=rd_in_mv)
    map_it <- unmerged[[mv_it]][["hyper_info"]]$map
    # merged_map[rd_in_mv] gives NA's or full-model hyper indices of [ranefs in submodel mv_it]
    # names(.) goes back to the names, ie the full-model indices of [ranefs in submodel mv_it] once the NA have been removed, 
    #                                     ie keeping only ranefs in hyper terms
    hy_full_indices <- na.omit(merged_map[rd_in_mv])
    # rd_in_mv[names(.)] thus gives the sub-model indices of ranefs in hyper terms
    # unique(map_it[.]) Then gives the submodel hyper indices of these ranefs
    hy_in_submv <- unique(map_it[inverse_map[names(hy_full_indices)]])
    map_hy_mv[[mv_it]] <- structure(unique(hy_full_indices), names=hy_in_submv) # values are full-model hyper-indices, names are submodel hyper indices
    # same value/names relationship as for map_rd_mv 
  }
  return(list2env(list(map=merged_map,ranges=merged_ranges, template= template, summingMat=summingMat, map_hy_mv=map_hy_mv),
                  parent=emptyenv()))
}

.merge_processed <- function(calls_W_processed, data, init=list(), control.HLfit=list(), method="ML", verbose=NULL, init.HLfit=list(),
                             covStruct=NULL, corrMatrix=NULL, adjMatrix=NULL, distMatrix=NULL, control.dist=list()) {
  # this fn passes no '...' so has no '...'
  nmodels <- length(calls_W_processed)
  namedlist <- structure(vector("list",nmodels), names=seq_len(nmodels))
  ### Fill lists for further processing:
  unmerged <- predictors <- families <- prior.weights <- clik_fns <- phiFixs <- Ys <- pS_fixef_phi <- namedlist
  AMatrices <- adjMatrices <- corrMatrices <- fixef_off_termsS <- fixef_termsS <- fixef_levelsS <- validrownames <- namedlist
  for (mv_it in seq_len(nmodels)) {
    unmerged[[mv_it]] <- calls_W_processed[[mv_it]][["processed"]]
    predictors[[mv_it]] <- unmerged[[mv_it]][["predictor"]]
    families[[mv_it]] <- unmerged[[mv_it]][["family"]]
    prior.weights[[mv_it]] <- unmerged[[mv_it]][["prior.weights"]] ## may be quote extression, etc.
    clik_fns[[mv_it]] <- unmerged[[mv_it]][["clik_fn"]]
    phiFixs[mv_it] <- list(unmerged[[mv_it]][["phi.Fix"]]) # ! syntax to allow explicit NULL's
    terms_info_it <- unmerged[[mv_it]][["main_terms_info"]]
    Ys[[mv_it]] <- terms_info_it[["Y"]] # one reason for parsing frames that way is to allow $Y as argument of .get_inits_from_glm()
    fixef_off_termsS[[mv_it]] <- terms_info_it[["fixef_off_terms"]]
    fixef_termsS[mv_it] <- list(terms_info_it[["fixef_terms"]]) 
    fixef_levelsS[mv_it] <- list(terms_info_it[["fixef_levels"]]) 
    corr_info_it <- unmerged[[mv_it]][["corr_info"]]     # already preprocessed info for:
    AMatrices[mv_it] <- list(corr_info_it$AMatrices)
    adjMatrices[mv_it] <- list(corr_info_it$adjMatrices)
    corrMatrices[mv_it] <- list(corr_info_it$corrMatrices)
    # corr_info_it also has ""corr_families" "corr_types""cov_info_mats" "G_diagnosis"
    pS_fixef_phi[[mv_it]] <- unmerged[[mv_it]][["p_fixef_phi"]] 
    #geo_infos[mv_it] <- list(unmerged[[mv_it]][["geo_info"]]) # each geo_info is a ranef-list of environments,  ou NULL; le list() est pourle second cas
  }
  #attr(phiFixs,"anyNULL") <- .anyNULL(phiFixs); attr(phiFixs,"allNULL") <- .allNULL(phiFixs) ## not a good idea bc its too easy to change the elements withouth changing the attributes
  #
  ### initialize 'merged' from unmerged[[1L]]:
  merged <- list2env(list(envir=list2env(list(), parent=environment(HLfit))))
  ## From unmerged[[1L]][[st]] to 'merged': that assignment should ultimately be only for elements not recursively updated:
  for (st in c(#"AUGI0_ZX",
    "REMLformula", # __FIXME__ this will really handle only standard ML (REMLformula has an isML attr) 
    #                                                      or standard REML (REMLformula is NULL). 
    # => no attempt to look REMLformula over models below (But  is built iteratively). 
    "verbose","control.glm","HL","p_v_obj",#"rand.families",
    "spaMM_tol",
    "break_conv_logL",
    "objective","port_env")
  ) assign(st,value=unmerged[[1L]][[st]],envir=merged)
  
  merged[["For"]] <- "fitme" # does not appear necessary except with adjmatrixas $For needed to determine inner_estim_adj_rho in .determine_spprec()
  merged[["phi.Fix"]] <- phiFixs
  merged[["p_fixef_phi"]] <- pS_fixef_phi
  merged$predictor <- predictors # that will be passed by HLfit_body to the result...
  merged$clik_fn <- clik_fns
  #
  ### Local values from unmerged[[1L]]:
  # That block should ultimately be for elements recursively updated
  for (st in c("off","y","BinomialDen","main_terms_info","iter_mean_dispVar","iter_mean_dispFix",
               "max.iter","maxLambda","models","vecdisneeded","bin_all_or_none")) assign(st,value=unmerged[[1L]][[st]])
  # Operatiosn on 'models' form the first submodel:
  LMMbool  <- attr(models,"LMMbool")
  GLMMbool  <- attr(models,"GLMMbool")
  LLM_const_w  <- attr(models,"LLM_const_w")
  GLGLLM_const_w  <- attr(models,"GLGLLM_const_w")
  GLMbool <- (models[["eta"]]=="etaGLM")
  LMbool <- GLMbool && unmerged[[1L]]$family$family == "gaussian" && attr(unmerged[[1L]]$family, "canonicalLink")
  const_wresid  <- attr(models,"const_wresid")
  phi_models <- character(nmodels)
  phi_models[[1L]] <- models[["phi"]]
  #
  vec_nobs <- integer(nmodels)
  vec_nobs[1L] <- length(y)  
  ZAlist <- unmerged[[1L]]$ZAlist
  ZAlist <- .correct_ZA_mv_ranCoefs(ZAlist, mv_it=1L)
  ZAlist <- .merge_ZAlists(list(), ZAlist, 0L, vec_nobs[1L], 1L)
  merged_X <- .merge_Xs(NULL, unmerged[[1L]][["AUGI0_ZX"]]$X.pv, mv_it=1L)
  merged_X.Re <- .merge_Xs(NULL, unmerged[[1L]][["X.Re"]], mv_it=1L, REML=TRUE)
  vec_ncol_X <- integer(nmodels)
  vec_ncol_X[1L] <- ncol(merged_X)
  validrownames[[1L]] <- rownames(unmerged[[1L]][["data"]])
  # Recursive updating:
  for (mv_it in (seq_len(nmodels-1L)+1L)) {
    p_i <- unmerged[[mv_it]]
    # merged$predictor... I should try to get rid of its use in HLfit_body, at least... FIXME
    ### random effects stuff
    validrownames[[mv_it]] <- rownames(p_i[["data"]])
    vec_nobs[mv_it] <- length(p_i[["y"]])  
    cum_nobs <- cumsum(c(0L, vec_nobs)) # quick & dirty rebuild cum_nobs from scratch in each iteration.
    ZAlist_i <- p_i$ZAlist
    ZAlist_i <- .correct_ZA_mv_ranCoefs(ZAlist_i, mv_it=mv_it)
    ZAlist <- .merge_ZAlists(ZAlist, ZAlist_i, nobs1=cum_nobs[mv_it], nobs2=vec_nobs[mv_it], mv_it)
    ### response stuff:
    y <- c(y,p_i[["y"]])
    BinomialDen <- c(BinomialDen,p_i[["BinomialDen"]])
    off <- c(off, p_i[["off"]])
    X_i <- p_i$AUGI0_ZX$X.pv # .unscale(p_i$AUGI0_ZX$X.pv)
    merged_X <- .merge_Xs(merged_X, X_i, mv_it=mv_it)  
    merged_X.Re <- .merge_Xs(merged_X.Re, p_i[["X.Re"]], mv_it=mv_it, REML=TRUE)  
    vec_ncol_X[mv_it] <- ncol(X_i)
    bin_all_or_none <- bin_all_or_none && p_i[["bin_all_or_none"]]
    # .setattr_G_LMMbool() examines single phi model so would need to be extended in order to replace the following lines
    GLMbool_it <- p_i[["models"]][["eta"]]=="etaGLM"
    LMbool_it <- GLMbool_it && p_i$family$family == "gaussian" && attr(p_i$family, "canonicalLink")
    modattrs_it <- attributes(p_i[["models"]])
    LLM_const_w_it <- modattrs_it[["LLM_const_w"]]
    GLGLLM_const_w_it <- modattrs_it[["GLGLLM_const_w"]] 
    LMMbool_it <- modattrs_it[["LMMbool"]]
    GLMMbool_it <- modattrs_it[["GLMMbool"]]
    const_wresid_it <- modattrs_it[["const_wresid"]]
    LLM_const_w <- ((LLM_const_w  || (LMbool && const_wresid) ) && LLM_const_w_it) ||
      ( LLM_const_w && LMbool_it && const_wresid_it )
    GLGLLM_const_w <- ((GLGLLM_const_w  || (GLMbool && const_wresid) ) && GLGLLM_const_w_it) ||
      (GLGLLM_const_w && GLMbool_it && const_wresid_it ) # (constant weights over iterations of IRLS)
    LMMbool <- ((LMMbool || LMbool) && LMMbool_it) ||
      (LMMbool && LMbool_it)
    GLMMbool <- ((GLMMbool || GLMbool) && GLMMbool_it) ||
      (GLMMbool && GLMbool_it)
    const_wresid  <- const_wresid && const_wresid_it
    GLMbool <- GLMbool && GLMbool_it
    LMbool <- LMbool && LMbool_it
    iter_mean_dispVar <- max(iter_mean_dispVar, p_i[["iter_mean_dispVar"]])
    iter_mean_dispFix <- max(iter_mean_dispFix, p_i[["iter_mean_dispFix"]])
    max.iter <- max(max.iter, p_i[["max.iter"]])
    maxLambda <- min(maxLambda, p_i[["maxLambda"]])
    vecdisneeded <- (vecdisneeded | p_i[["vecdisneeded"]])
    phi_models[[mv_it]] <- p_i[["models"]][["phi"]]
  }
  .check_mv_in_submodels(ZAlist)
  attr(data,"validrownames") <- validrownames
  merged[["data"]] <- data
  merged[["vec_nobs"]] <- structure(vec_nobs, cum_nobs=cum_nobs) # objetc will contain multiple copies of cum_nobs attribute, for conveniency
  merged[["prior.weights"]] <- prior.weights
  merged[["off"]] <- off
  merged[["bin_all_or_none"]] <- bin_all_or_none 
  merged[["iter_mean_dispVar"]] <- iter_mean_dispVar
  merged[["iter_mean_dispFix"]] <- iter_mean_dispFix
  merged[["max.iter"]] <- max.iter
  merged[["maxLambda"]] <- maxLambda
  merged[["vecdisneeded"]] <- vecdisneeded
  # 
  augZXy_cond_inner <- spaMM.getOption("allow_augZXy")
  if (is.null(augZXy_cond_inner)) augZXy_cond_inner <- TRUE ## for .makeCovEst1()
  if (augZXy_cond_inner) augZXy_cond_inner <-   LMMbool
  if (augZXy_cond_inner) augZXy_cond_inner <- ( is.null(merged_X.Re) || ! ncol(merged_X.Re)) ## exclude non-standard REML (avoiding NCOL(NULL)=1)
  merged$augZXy_cond <- structure(FALSE, inner=augZXy_cond_inner) # augZXy_cond would impose a unique phi accross submodels
  #
  attr(phi_models,"anyHGLM") <- any(phi_models=="phiHGLM")
  models[["phi"]] <- phi_models # som 'models' is a list whose element 'phi' is a vector
  attr(models, "LMMbool") <- LMMbool # add more attributes to avoid clumsy tests later
  attr(models,"GLMMbool") <- GLMMbool 
  attr(models,"LLM_const_w") <- LLM_const_w 
  attr(models,"GLGLLM_const_w") <- GLGLLM_const_w 
  merged[["models"]] <- models
  models <- NULL # make sure we work on only one 'models'
  residProcesseds <- residModels <- vector("list", length(merged[["models"]][["phi"]]))
  for (mv_it in seq_along(phi_models)) {
    residModels[mv_it] <- list(unmerged[[mv_it]]$residModel)
    residProcesseds[mv_it] <- list(unmerged[[mv_it]]$residProcessed)
  }
  merged$residProcesseds <- residProcesseds 
  merged$residModels <- residModels 
  #
  #
  attr(families, "cum_nobs") <- cum_nobs
  #
  has_estim_families_par <- FALSE
  for (mv_it in seq_along(unmerged)) {
    family_it <- families[[mv_it]]
    has_estim_families_par <- ((family_it$family=="negbin" && inherits(substitute(shape, env=environment(family_it$aic)),"call")) ||
                                 (family_it$family=="COMPoisson" && inherits(substitute(nu, env=environment(family_it$aic)),"call")))
    if (has_estim_families_par) break
  }
  attr(families,"has_estim_families_par") <- has_estim_families_par 
  #
  merged$families <- families
  # namestable <- table(colnames(merged_X)) 
  # if (length(namestable)<ncol(merged_X)) {
  #   vec_col_mod <- rep(seq(nmodels),vec_ncol_X)
  #   allnames <- names(namestable)
  #   for (namit in allnames) {
  #     if (namestable[namit]>1L) {
  #       whichcols <- which(colnames(merged_X)==namit)
  #       whichmods <- vec_col_mod[whichcols]
  #       colnames(merged_X)[whichcols] <- paste0(namit,"_",whichmods)
  #     }
  #   }
  # }
  merged_X <- .post_process_X(X.pv=merged_X, HL=merged$HL, rankinfo=NULL, 
                              sparse_X=.determine_sparse_X_mv(merged$main_terms_info, X.pv= merged_X, vec_nobs=vec_nobs, unmerged=unmerged) ) 
  # processing of merged_X and other elements of AUGI0_ZX:
  attr(merged_X,"cum_ncol") <- cumsum(c(0L,vec_ncol_X))
  attr(merged_X,"cum_nobs") <- cum_nobs
  #
  if ( ! is.null(merged_X.Re) && ncol(merged_X.Re)) { # non-standard REML... 
    # attr(., "extra_vars") has aleardy been updated by .merge_Xs(., REML=TRUE), but unrestricting_cols culd not as it refers to X.pv cols 
    unrestricting_cols <- which(colnames(merged_X) %in% setdiff(colnames(merged_X),colnames(merged_X.Re))) ## not in X.Re
    distinct.X.ReML <- c(length(unrestricting_cols), length(attr(merged_X.Re,"extra_vars"))) ## TWO integers often used as booleans 
    attr(merged_X.Re,"distinct.X.ReML") <- distinct.X.ReML 
    if ( distinct.X.ReML[1L]) attr(merged_X.Re,"unrestricting_cols") <- unrestricting_cols # cols of X.pv not in X.Re
  }
  merged[["X.Re"]] <- merged_X.Re
  #
  dim(y) <- c(length(y),1L) # cf comments in .preprocess(); but vector format worked in all tests for fitmv up to version 3.5.115; 
  merged[["y"]] <- y
  merged[["BinomialDen"]] <- BinomialDen
  main_terms_info[["Y"]] <- Ys 
  main_terms_info[["fixef_off_terms"]] <- fixef_off_termsS
  main_terms_info[["fixef_terms"]] <- fixef_termsS
  main_terms_info[["fixef_levels"]] <- fixef_levelsS
  merged[["main_terms_info"]] <- structure(main_terms_info, vec_nobs=vec_nobs)
  attr(ZAlist,"map_rd_mv") <- map_rd_mv <- .map_rd_mv(ZAlist, unmerged)
  # ZAlist including map_rd_mv available -> rand.families
  merged[["rand.families"]] <- rand.families <- .merge_rand_families(unmerged, ZAlist=ZAlist) 
  merged$lcrandfamfam <- attr(rand.families,"lcrandfamfam") ## else remains NULL
  #
  # map_rd_mv available -> exp_ranef_types -> soon used by .assign_corr_types_families()
  exp_barlist <- .process_bars_mv(predictors, map_rd_mv)
  exp_ranef_strings <- .process_bars(barlist=exp_barlist,expand=FALSE, as_character=TRUE) ## no need to expand again
  if (nrand <- length(exp_ranef_strings)) {
    names(ZAlist) <- seq_len(nrand) # not sure where it is used, but conform to fact that univar-resp ZAlist is named.
    exp_ranef_types <- attr(exp_ranef_strings,"type") ## expanded
    attr(ZAlist,"exp_ranef_strings") <- exp_ranef_strings ## expanded 
    attr(ZAlist,"exp_ranef_types") <- exp_ranef_types ## expanded
    #
    ####### Builds $corr_info (ASAP to assign_cov_matrices ASAP):
    merged$corr_info <- corr_info <- new.env() ## do not set parent=emptyenv() else with(corr_info,...) will not find trivial fns such as `[`
    .assign_corr_types_families(corr_info, exp_ranef_types) # provides 'corr_types', soon necessary, and 'corr_families'
    # follow the order of .preprocess():
    ## Assigns $corr_info list of matrices BEFORE determining sparse precision:
    corr_info$corrMatrices <- .merge_mv_list(corrMatrices, merged, ZAlist=ZAlist, full=TRUE) 
    corr_info$adjMatrices <- .merge_mv_list(adjMatrices, merged, ZAlist=ZAlist, full=TRUE)
    corr_info$AMatrices <- .merge_mv_list(AMatrices, merged, ZAlist=ZAlist, full=TRUE)
    ## for (mult)IMRF:
    ## HLCor_body -> .assign_geoinfo_and_LMatrices_but_ranCoefs() expects attr(processed$ZAlist, "AMatrices")
    ## Otherwise it may stop on  .............................................. -> .calc_IMRF_Qmat()
    ## In .preprocess they are added through ZAlist <- .calc_ZAlist(Zlist=Zlist, AMatrices=corr_info$AMatrices),
    ##    using corr_info$AMatrices previously assigned by
    ##    .assign_AMatrices_IMRF(corr_info, Zlist, exp_barlist=exp_barlist, processed$data, control_dist=processed$control_dist)
    #
    ## So we need first to replicate the effect on corr_info$AMatrices of 
    ##    .assign_AMatrices_IMRF(corr_info, Zlist, exp_barlist=exp_barlist, processed$data, control_dist=processed$control_dist)
    ## Then to replicate the effect of ZAlist <- .calc_ZAlist(Zlist=Zlist, AMatrices=corr_info$AMatrices) on attr(processed$ZAlist, "AMatrices")
    .add_cov_matrices__from_mv_global(corr_info, covStruct=covStruct, corrMatrix=corrMatrix, adjMatrix=adjMatrix) # using $corr_info$corr_type
    ## for adjacency:
    # filling corr_info$adjMatrices was delayed until the above .add_cov_matrices__from_mv_global. 
    # So they are not in the  corr_info in each unmerged[[mv_it]]. We might copy them by
    # for (mv_it in seq_along(unmerged)) {
    #   unmerged[[mv_it]]$corr_info$adjMatrices[names(map_rd_mv[[mv_it]])] <- corr_info$adjMatrices[map_rd_mv[[mv_it]]]
    #   # and let us do the same for the other auxiliary matricess ? 
    # }
    #
    # "cov_info_mats" will be provided by .init_assign_geoinfo(), and "G_diagnosis" is computed when needed
    ####### 
    #
    # Using corr_info:
    merged$control_dist <- .preprocess_control.dist(control.dist, corr_info$corr_types)
    #
    merged$init_HLfit <- .preprocess_init.HLfit(init.HLfit, corr_info)
  
    merged$is_spprec <- .wrap_determine_spprec(control.HLfit, ZAlist=ZAlist, processed=merged, X.pv=merged_X)
    ## post-processing of corr_info depending on sparse_precision
    .process_corr_info_spprec(corr_info=corr_info, For="fitme",sparse_precision=merged$is_spprec)
    # Heavily using corr_info, and spprec:
    ZAlist <- .init_assign_geoinfo(processed=merged, ZAlist=ZAlist, For="fitme", 
                                   exp_barlist=exp_barlist, distMatrix=distMatrix)  
    vec_normIMRF <- .calc_vec_normIMRF(exp_ranef_terms=attr(ZAlist, "exp_ranef_terms"), corr_types=corr_info$corr_types)   
    if (any(vec_normIMRF)) {
      Zlist <- .merge_Zlists(list(), attr(unmerged[[1L]]$ZAlist,"Zlist"), 0L, vec_nobs[1L], 
                             ranefs1=NULL,
                             ranefs2=map_rd_mv[[1]],
                             1L) # attribute present when IMRF present
      for (mv_it in (seq_len(nmodels-1L)+1L)) {
        ranefs1 <- names(Zlist)
        ranefs2 <- map_rd_mv[[mv_it]]
        Zlist <- .merge_Zlists(Zlist, attr(unmerged[[mv_it]]$ZAlist,"Zlist"), nobs1=cum_nobs[mv_it], nobs2=vec_nobs[mv_it], 
                               ranefs1=names(Zlist),
                               ranefs2=map_rd_mv[[mv_it]],
                               mv_it) # attribute present when IMRF present
      }
      Zlist <- Zlist[seq_along(Zlist)] # remove attributes to make clear they are not needed
      attr(ZAlist,"Zlist") <- Zlist  
    }
    attr(ZAlist,"AMatrices") <- corr_info$AMatrices 
    merged[["ZAlist"]] <- ZAlist
    merged$hyper_info <- .merge_hyper_infos(ZAlist, unmerged)
  } else {
    merged$is_spprec <- FALSE ## for .do_TRACE()
    merged$init_HLfit <- init.HLfit
  }
  merged$QRmethod <- .choose_QRmethod(ZAlist, corr_info=merged$corr_info, 
                                      is_spprec=merged$is_spprec, processed=merged, control.HLfit=control.HLfit)
  .check_time_G_diagnosis(.provide_G_diagnosis, processed=merged)
  #
  # merged_X <- .scale(merged_X) not necessary since the merged X's are already scaled
  if (nrand) {
    merged$models[["eta"]] <- "etaHGLM" 
    vec_n_u_h <- unlist(lapply(merged$ZAlist,ncol)) 
    merged[["cum_n_u_h"]] <- cumsum(c(0L, vec_n_u_h))
    nrd <- merged[["cum_n_u_h"]][nrand+1L]
    if (nrd==1L) {
      warning("Found a single random effect with a *single level*. Check formula?", immediate.=TRUE)
    }
    merged$models[["lambda"]] <- rep("lamScal",nrand) ## even for adjacency, random slope...
    merged[["reserve"]] <- .preprocess_arglists(merged)
    # Need this as long as we need  .calc_optim_args_mv(processed$unmerged....) ie as long as  .calc_optim_args() handles only single $family 
    for (mv_it in seq_along(unmerged)) {
      rd_in_mv <- map_rd_mv[[mv_it]]
      unmerged[[mv_it]]$geo_info <- merged$geo_info[rd_in_mv]
    }
  } else merged$models[["eta"]] <- "etaGLM" 
  merged[["AUGI0_ZX"]] <- .init_AUGI0_ZX(X.pv=merged_X, vec_normIMRF, ZAlist=merged$ZAlist, nrand, n_u_h=nrd, 
                                         sparse_precision=merged$is_spprec, 
                                         as_mat=.eval_as_mat_arg(merged))
  if (nrand) {
    merged[["X_lamres"]] <- .calc_X_lamres(merged, models=merged$models, ZAlist=ZAlist, nrand=nrand)
    # ranCoefs processing moved to after merging the ranCoefs in fitmv(); 
    # this does not seem to impact the final operation within the present fn; but still calls for some improvement
  }
  merged[["unmerged"]] <- unmerged # A list whose each element is the processed arg of each processed call
  .check_identifiability_LMM_mv(merged, vec_nobs=vec_nobs, map_rd_mv=map_rd_mv) # uses $unmerged
  merged[["LevenbergM"]] <- .preprocess_LevM(control.HLfit$LevenbergM, merged, nrand=length(ZAlist)) # uses $models, $HL & optionally $bin_all_or_none & $cum_n_u_h
  #
  merged$verbose <- .reformat_verbose(verbose,For="fitme") 
  .do_TRACE(merged)
  class(merged) <- c("mvarglist", class(merged))
  return(merged)
  # models <- NULL # prevents a R CMD check NOTE: .merge_processed: no visible binding for global variable 'models' (!)
  # Apparently the check is satisfied only if any variable is assigned  by recognized means, 
  # including '<-' but not assign(). Even this silly placement of '<-' works.
  # (For other variables such as 'off' created by assign(), off <- ....(off) plays this role)
}

# does not canonize (to canonical scale) but formats corrPars to structured list and lambda to standardly-named vector.
.reformat_parlist <- function(parlist, processed, ZAlist=processed$ZAlist, namesTerms=attr(ZAlist,"namesTerms")) { # for merging user's inputs mid-processing. Not for optimization hence no transformed params
  parlist$lambda <- .reformat_lambda(parlist$lambda, nrand=length(ZAlist), namesTerms=namesTerms, full_lambda=FALSE)
  parlist$phi <- .reformat_phi(parlist$phi,n_models=length(processed$unmerged), full_phi=FALSE)
  # This fn is only for fixed values and fixed COMP_nu is best given by the family argument...   
  # parlist$COMP_nu <- .reformat_phi(parlist$COMP_nu,n_models=length(processed$unmerged), full_phi=FALSE)
  # parlist$NB_shape <- .reformat_phi(parlist$NB_shape,n_models=length(processed$unmerged), full_phi=FALSE)
  parlist <- .reformat_corrPars(parlist,corr_families=processed$corr_info$corr_families)
}

fitmv <- function(submodels, data, fixed=NULL, init=list(), lower=list(), upper=list(),
                  control=list(), # needed to avoid partial matching of explicit 'control' argument with 'control.dist' one (bug when the latter is used) 
                  control.dist = list(), method="ML", init.HLfit=list(), ...) { # explicit arguments or dots depending on what requires specific documentation.
  assign("spaMM_glm_conv_crit",list(max=-Inf) , envir=environment(spaMM_glm.fit))
  time1 <- Sys.time()
  oricall <- match.call(expand.dots=TRUE) ## mc including dotlist
  oricall$control.HLfit <- eval(oricall$control.HLfit, parent.frame()) # to evaluate variables in the formula_env, otherwise there are bugs in waiting 
  n_models <- length(submodels) # so the promise is already evaluated here...
  calls_W_processed <- fixedS <- vector("list",n_models)
  for (mv_it in seq_along(calls_W_processed)) { # call .preprocess() on each submodel
    call_ <- oricall
    call_["submodels"] <- NULL # so that it remains in call_ the arguments others than mv.
    ## I need to match the names of mv[[mit]] to those of a fitme call to make sure that they all named...
    call_["fixed"] <- NULL ## so that the lambda fixing (in particular) is not the default value for each processed call
    ## *** global arguments => avoid to mix thenm with local arguments 
    ##     (although this is stricly necessary only for covStruct since...) ***  
    call_["corrMatrix"] <- NULL # not strictly necess since single matrix so never a problem of matching ranefs: The global corrMatrix is a locally usable corrMatrix, 
    call_["adjMatrix"] <- NULL # not strictly necess ... same comment...
    call_["covStruct"] <- NULL # => important to remove it since ranefs cannot be matched in .preprocess().
    #
    call_["init.HLfit"] <- NULL # We could leave it, that would be useless. OTOH, .merge_processed() will use it.
    ## *** ***  
    matched_args_it <- match.call(fitme, do.call("call",c(list(name="fitme"), submodels[[mv_it]]))) # match the elements of mv[[mv_it]] to those of a call to fitme
    #    => any explicit 'fixed' in the submodel will be in matched_args_it ; same for init but it is used by .preprocess() for something not relevant here (augZXy-related) 
    if ( ! is.null(matched_args_it[["init"]]) ) warning("'init' in sub-model is ignored. Use fitmv()'s 'init' argument instead.", immediate.=TRUE)
    # : I could implement a merging at a later step (it's not useful at .preprocess_fitme() step) but it does not seem worth the code.
    if ( ! is.null(matched_args_it[["distMatrix"]]) ) warning("'distMatrix' in sub-model is ignored. Use fitmv()'s 'distMatrix' argument instead.", immediate.=TRUE)
    # : merging distMatrices would be difficult since they are dispersed in element of geo_info, and locations therein may have been subsetted.
    #
    ## to make update(, formula.=<.>) work, fitmv handles a formula. argument through the '...'
    # we use it to update the 'formula' argument of each matched_args_it, and remove "formula." from the call_ to be evaluated.
    call_["formula."] <- NULL 
    if ( ! is.null(form._it <- oricall$formula.[[mv_it]])) matched_args_it[["formula"]] <- form._it
    #
    for (st in names(matched_args_it)[-1]) call_[[st]] <- matched_args_it[[st]] # so args within the mv[[mv_it]] list add to or replace those outside of mv
    # so if there was no explicit fixed in the submodel there is no fixed in matched_args_it nor in current call_. Hence...
    fixedS[[mv_it]] <- .modify_list(list(), eval(call_[["fixed"]], parent.frame())) # Ensures we have a list... but it's ugly.
    call_$formula <- .preprocess_formula(call_$formula)
    #
    call_[["what_checked"]] <- "arguments for .preprocess_fitme()" 
    call_[[1L]] <- get(".check_args_fitme", asNamespace("spaMM"), inherits=FALSE) 
    call_ <- eval(call_,parent.frame()) # 
    call_["what_checked"] <- NULL 
    #
    call_[["For"]] <- "fitmv"
    if ( ! is.null(main_terms_info <- attr(data,"updated_terms_info"))) { # from update_resp -> .update_main_terms_info() 
      main_terms_info_it <- list(mf=main_terms_info$mf[[mv_it]],fixef_off_terms=main_terms_info$fixef_off_terms[[mv_it]],
                                 fixef_terms=main_terms_info$fixef_terms[[mv_it]],
                                 fixef_levels=main_terms_info$fixef_levels[[mv_it]])
      #class(main_terms_info_it) <- "HLframes" # we tag the result again so that .preprocess() will recognize it as coming from .update_data()
      call_[["data"]] <- structure(data, updated_terms_info=main_terms_info_it)
    }
    call_[[1L]] <- get(".preprocess_fitme", asNamespace("spaMM"), inherits=FALSE) 
    calls_W_processed[[mv_it]] <- eval(call_,parent.frame()) # returns modified call including an element 'processed'
    residProcessed <- calls_W_processed[[mv_it]]$processed$residProcessed
    if ( ! is.null(validrownames <- attr(residProcessed$data, "validrownames"))) { # post-fit (confint...)
      residProcessed$data <- residProcessed$data[validrownames[[mv_it]],, drop=FALSE]
      attr(residProcessed$data, "validrownames") <- NULL
    }
    calls_W_processed[[mv_it]][["processed"]][["augZXy_cond"]] <- FALSE # not only to ensure the merged value but also for init.optim for each  
  }
  #
  ##### merge and finalize preprocessing
  mc <- oricall # with only the explicit arguments of the call (no 'fixed' if ...)
  mc["submodels"] <- NULL
  mc["formula."] <- NULL 
  mc[["what_checked"]] <- "fitmv() call" 
  mc[[1L]] <- get(".check_args_fitme", asNamespace("spaMM"), inherits=FALSE) 
  eval(mc,parent.frame()) # -> abyss 
  mc["what_checked"] <- NULL 
  mc["fixed"] <- NULL
  mc["upper"] <- NULL # to be used only in fitmv_body()
  mc["lower"] <- NULL
  mc["control"] <- NULL
  mc[["calls_W_processed"]] <- calls_W_processed
  # the fact that promises are evaluated within a call-execution is "local": they will appear not evaluated
  # when we reuse a call (here mc). E.g. corrMatrix=as_precision(.) would be evaluated twice 
  # => We need to put the evaluated value in the call list. 
  # Next line ad-hoc for corrMatrix (__F I X M E__?: What about other arguments ? Which would benefit from some preprocessing?)
  # if ("corrMatrix" %in% ...names()) # is an R >= 4.1.0 syntax 
  if ("corrMatrix" %in% names(mc)) mc["corrMatrix"] <- list(eval(mc[["corrMatrix"]])) 
  mc[[1L]] <-  get(".merge_processed", asNamespace("spaMM"), inherits=FALSE)
  merged <- eval(mc, parent.frame()) # means that arguments of *.merge_processed()* must have default values as mc does not contains defaults of fitmv()
  #
  fixed <- .reformat_parlist(fixed,processed = merged) # reformat user's global 'fixed' argument
  fixedS <- lapply(fixedS, .reformat_parlist, processed = merged)
  fixedS <- .merge_mv_parlist(fixedS, merged) # now fixedS is a single parlist from the  sub-models specifications
  fixedS <- .modify_list(fixedS,fixed) # now fixedS is a single parlist from both sub-model and global specifications
  #fixed <- .canonizeRanPars(ranPars=fixed,corr_info=merged$corr_info, checkComplete = FALSE, rC_transf=.spaMM.data$options$rC_transf)
  # These infos are ultimately used by summary() to distinguish "fix" from outer "var":
  merged[["lambda.Fix"]] <- .reformat_lambda(.getPar(fixed,"lambda"), nrand=length(merged$ZAlist), 
                                             namesTerms=attr(merged$ZAlist,"namesTerms"), full_lambda=TRUE)
  # HLfit_body() expects merged[["phi.Fix]] to be a full-length list, possibly with explicit NULLs.
  # merged[["phi.Fix"]] from .merge_processed() should be so, and .modify_list() should keep it so.
  merged[["phi.Fix"]] <- .modify_list(merged[["phi.Fix"]], fixedS$phi)
  #
  ranCoefs <- .getPar(fixedS,"ranCoefs") ## may be NULL
  merged$ranCoefs_blob <- .process_ranCoefs(merged, ranCoefs, use_tri_CORREL=TRUE) 
  merged$AUGI0_ZX$envir$finertypes[merged$ranCoefs_blob$isRandomSlope] <- "ranCoefs" 
  #
  mc["fixed"] <- oricall["fixed"]
  mc["upper"] <- oricall["upper"]
  mc["lower"] <- oricall["lower"]
  mc["control"] <- oricall["control"]
  mc["calls_W_processed"] <- NULL
  mc[["fixedS"]] <- fixedS # to build and merge the inits
  mc$processed <- merged
  mc[[1L]] <-  get("fitmv_body", asNamespace("spaMM"), inherits=FALSE)
  hlcor <- eval(mc,parent.frame()) 
  #
  for (mit in seq_along(calls_W_processed)) {
    if ("control.dist" %in% names(submodels[[mit]])) {
      oricall[["mv"]][["control.dist"]] <- calls_W_processed[[mit]][["control_dist"]]
    } else oricall$"control.dist" <- calls_W_processed[[mit]][["control_dist"]] ## maybe # [[]] <- does not work if [["control_dist"]] orginally absent
    hlcor$call <- oricall ## this is a call to fitmv()
  }
  lsv <- c("lsv",ls())
  if ( ! inherits(hlcor,"HLfitlist") && ! is.call(hlcor) ) {
    hlcor$how$fit_time <- .timerraw(time1)
    hlcor$how$fnname <- "fitmv"
    hlcor$fit_time <- structure(hlcor$how$fit_time,
                                message="Please use how(<fit object>)[['fit_time']] to extract this information cleanly.")
  }
  rm(list=setdiff(lsv,"hlcor")) ## empties the whole local envir except the return value
  return(hlcor)
}




