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
  return(map_rd_mv) # contains the rank of each rd of the total model (value) in the rds of each submodel (name)
  # e.g. *IF* there are 3 ranefs, one in each submodel,  map_rd_mv[[mv_it]] = c("1"=mv_it)
}
# Idiom:
# which_in_subm <- which(map_rd_mv[[mv_it]] %in% corr_info$is_cF_internally ) # so eg from (c("1"=1,"2"=3)) %in% 3), rd_in_mv is 2 
# pos_in_sub <- map_rd_mv[[mv_it]][which_in_subm]
# pos_in_merged <- as.integer(names(pos_in_sub)) 
# but pos_in_merged is more directly provided by map_rd_mv[[mv_it]] (_F I X M E__ rename?)

.merge_rand_families <- function(unmerged, map_rd_mv=attr(ZAlist, "map_rd_mv"), nrand=length(ZAlist), ZAlist) {
  n_submodels <- length(unmerged)
  rand_families <- vector("list", nrand)
  for (rd in seq_len(nrand)) { # indices for full model
    for (mv_it in seq_len(n_submodels)) {
      if (length(rd_in_mv <- which(map_rd_mv[[mv_it]]==rd))) { # gets the position in map_rd_mv[[mv_it]] 
        # and the name associated to this position, which the position in unmerged[[mv_it]][["rand.families"]]:  
        #rd_in_submodel <- names(rd_in_mv)
        newfam <- unmerged[[mv_it]][["rand.families"]][[names(rd_in_mv)]] # the names of the positions are the rd of the submodel
        if ( is.null(rand_families[[rd]])) {
          rand_families[[rd]] <- newfam
        } else if ( ! identical(rand_families[[rd]],newfam, ignore.environment=TRUE)) { 
          stop("Conflicting families for random effect term shared across models.")
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
      if (is_incid) attr(is_incid,"is01col") <- FALSE # conservative default (might be TRUE in some cases?)
      RHS_info1 <- attr(ZA1,"RHS_info")
      RHS_info2 <- attr(ZA2,"RHS_info")
      RHS_info <- list(splt=RHS_info1$splt,
                       dataordered_unique_levels=unique(RHS_info1$dataordered_unique_levels, 
                                                        RHS_info2$dataordered_unique_levels)
      )
      # levels of the LHS of the term:
      LHS_levels <- unique(c(attr(ZA1,"LHS_levels"),attr(ZA2,"LHS_levels")))
      # levels of the RHS...
      ulevels1 <- unique(colnames(ZA1))
      ulevels2 <- unique(colnames(ZA2))
      if ( ! identical(ulevels1,ulevels2)) { # add columns (here) before rbinding rows (later)
        #
        u_ranges_1 <- RHS_info1$AR1_block_u_h_ranges
        u_ranges_2 <- RHS_info2$AR1_block_u_h_ranges
        by_levels <- unique(names(u_ranges_1),names(u_ranges_2))
        AR1_block_u_h_ranges <- vector("list",length(by_levels))
        names(AR1_block_u_h_ranges) <- by_levels
        for (level_it in by_levels) AR1_block_u_h_ranges[[level_it]] <- range(c(u_ranges_1[[level_it]],u_ranges_2[[level_it]]))
        RHS_info$AR1_block_u_h_ranges <- AR1_block_u_h_ranges
        #
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
      } else RHS_info$AR1_block_u_h_ranges <- RHS_info1$AR1_block_u_h_ranges
      ori <- in1
    } else  { 
      # Identify source of attributes:
      if (length(in1)) {
        Zlistori <- ZAlist1
        ori <- in1
        which_mv <- attr(Zlistori[[ori]],"which_mv")
        ZA1 <- ZAlist1[[in1]]
        RHS_info <- attr(ZA1,"RHS_info")
        namesTerm <- attr(ZA1,"namesTerm")
        LHS_levels <- attr(ZA1,"LHS_levels")
        is_incid <- attr(ZA1,"is_incid") ## OK when adding rows of zeroes  
        ZA2 <- Matrix(0,ncol=ncol(ZA1), nrow=nobs2)
      } else {
        Zlistori <- ZAlist2
        ori <- in2
        which_mv <- mv_it
        ZA2 <- ZAlist2[[in2]]
        RHS_info <- attr(ZA2,"RHS_info")
        namesTerm <- attr(ZA2,"namesTerm")
        LHS_levels <- attr(ZA2,"LHS_levels")
        is_incid <- attr(ZA2,"is_incid") 
        ZA1 <- Matrix(0,ncol=ncol(ZA2), nrow=nobs1)
      } 
    } 
    ZAlist[[ran_it]] <- structure(rbind(ZA1,ZA2), 
                                  which_mv=which_mv, 
                                  namesTerm=namesTerm, 
                                  is_incid=is_incid,
                                  RHS_info=RHS_info, # RHS_info (here coming from one matrix) may not be appropriate for all models
                                  LHS_levels=LHS_levels)
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

# The matrices to be merged by .merge_Xs() have been individually rank-trimmed and scaled.
# .rankTrim() is also applied on the merged matrix.
# The final rankinfo attr stores rankinfo's of individual submodels in rankinfo$mvlist; this info is used 
# by anova.HLfit() -> contrast functions so the original rank-trimming is used. 
## Attributes have been added to each matrix, in the following order:
# .post_process_X() -> .rankTrim() on submodel matrices also computed {
#   namesOri attr, , keeps (modified) info about colnames of *submodel* matrices before *their* trimming 
#   assign, modified by .rankTrim() from the assign attr given by model.matrix(). 
#   assignOri, added by .rankTrim() copies the assign attr given by m odel.matrix().
# } 
# scaled:scale attr, added by .scale(),
# BUT: 
## Merging steps:
# As seen below, .merge_Xs() takes care of the namesOri and scaled:scale attrs
#    The info from assign and assignOri attrs appear to be lost in the merged return value, but
#      some attributes have been preprocessed and are put back after the .merge_Xs() call: cf
# attr(merged_X,"assign") <- assign_attr       # notably used by anova()
# attr(merged_X,"rankinfo") <- rankinfo_attr   # rankinfo$mvlist notably used by anova()
# attr(merged_X,"col_ranges") <- col_ranges
# attr(merged_X,"cum_nobs") <- cum_nobs
#
# The assignOri attr is apparently lost, (.anova.glm() will use it if present, and use "assign" otherwise)
#
## X2X product, if any, comes after merging in mv fits :.merge_processed() -> .designX_prod_mv(),
# which adds attr cols_lhs_X2X: the colnames of the merged_X before the X2X product. 
.merge_Xs <- function(X1,X2, mv_it, REML=FALSE) { 
  if ( ! is.null(colnames(X2))) colnames(X2) <- paste0(colnames(X2),"_",mv_it)
  if (is.null(X1)) { 
    X <- X2 # X can be NULL when merging REML matrices...
    if ( ! is.null(X)) {
      scale2 <- attr(X2,"scaled:scale")
      if ( length(scale2)) { # should be equiv to testing ncol(X2)
        names(scale2) <- paste0(names(scale2),"_",mv_it)
        attr(X,"scaled:scale") <- scale2
      } # else X1 was NULL (1st 'merge') and X2 had zero col (no fixef)
      if (length(attr(X2,"namesOri"))) {
        attr(X,"namesOri") <- paste0(attr(X2,"namesOri"),"_",mv_it)
      } else attr(X,"namesOri") <- attr(X1,"namesOri")
    }
  } else {
    XX1 <- cbind(X1,matrix(0, nrow=nrow(X1), ncol=ncol(X2)))
    colnames(XX1) <- c(colnames(X1),colnames(X2))
    XX2 <- cbind(matrix(0, nrow=nrow(X2), ncol=ncol(X1)),X2)
    X <- rbind(XX1,XX2)
    scale2 <- attr(X2,"scaled:scale")
    if ( length(scale2)) {
      names(scale2) <- paste0(names(scale2),"_",mv_it)
      attr(X,"scaled:scale") <- c(attr(X1,"scaled:scale"), scale2 )
    } else attr(X,"scaled:scale") <- attr(X1,"scaled:scale")
    if (length(attr(X2,"namesOri"))) {
      attr(X,"namesOri") <- c(attr(X1,"namesOri"),  paste0(attr(X2,"namesOri"),"_",mv_it))
    } else attr(X,"namesOri") <- attr(X1,"namesOri")
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
  if ( ! is.null(parlist$beta_prec)) names(parlist$beta_prec) <- mv_it
  if ( ! is.null(parlist$trbeta_prec)) names(parlist$trbeta_prec) <- mv_it
  if ( ! is.null(parlist$Tw_index)) names(parlist$Tw_index) <- mv_it
  if ( ! is.null(parlist$Tw_link)) names(parlist$Tw_link) <- mv_it
  # The input parlist$rdisPars may be a vector, or a list (even an empty one)
  # Notably: parlist may be provided by .calc_optim_args(): for each submodel that function provides 
  #           inits$init and other parlists as for a single resp fit -> parlist$rdisPars is a vector
  # but    : parlist may also be provided by .subPars() -> parlist$rdisPars is a list.
  # and... if a phimodel, its parameters can also end as a list in rdisPars... I should be careful thta it does not reaches here.
  # Changing the .calc_optim_args() return format looks bad: For clarity, it should be consistent with the single-response case.
  # Changing the .subPars() return format would make even less sense: it is a rather generic algorithm not distinguishing 
  # between different classes of parameters (or whatever its input means).
  if ( length(.unlist(parlist$rdisPars))) { # parlist has element "rdisPars"... not being an empty list...
    rdisPars <- parlist[["rdisPars"]] # either a vector or a single-element list
    if (is.numeric(rdisPars)) rdisPars <- list("1"=rdisPars) # -> single-element list, the element being the input rdisPars vector
    names(rdisPars) <- as.character(mv_it) # -> list with single element named <char_mv_it> being the input rdisPars vector
    parlist$rdisPars <- rdisPars # -> parlist has element "rdisPars" itself being a list with single element named <char_mv_it> being the input rdisPars vector
  } else parlist$rdisPars <- NULL
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
  if (length(lambda_merger)) {
    lambda_merger <- list(inits=list(init=list(lambda=lambda_merger),
                                     init.optim=list(trLambda=.dispFn(lambda_merger))))
    .modify_list(optim_blob,lambda_merger, obey_NULLs=FALSE) 
  } else optim_blob
}

# This is called by fitmv_body so not really preprocessing
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
                 rdisPars=structure(list(NULL),names=mv_it), # Another fake skeleton
                 beta_prec=phistr,
                 NB_shape=phistr, 
                 COMP_nu=phistr,
                 Tw_index=phistr,
                 Tw_link=phistr)
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
    #
    if (is.list(user_init_optim_it[["phi"]][[1L]])) {
      # at this point if we had a resid.model user_init_optim_it[["phi"]] 
      # contains phi$`1` (`1` for any mv_it) being a list of (fitted) parameters of the resid.model, 
      # which is sort-of not too bad but should not be used here
      user_init_optim_it[["phi"]] <- NULL # (__F I X M E___ perhaps we could better use the resid.model info?
      # how to pass info from here to init the resid.model ?
      # some test code would be test-mv -> confint(zut1,"(Intercept)_1") )
    } else {
      # we want user_init_optim_it[["phi"]] to be an unnamed scalar, or NULL, starting from any user input (input unnamed vector got names by .reformat_phi()) 
      # user_init_optim_it[["phi"]] <- user_init_optim[["phi"]][[as.character(mv_it)]] # OK for phi list or complete vector, but not incomplete vector
      user_init_optim_it[["phi"]] <- as.vector(.unlist(user_init_optim_it[["phi"]])) # replaces the trivial `1`
      # replaces
      # $ phi   :List of 1
      #  ..$ 1: num <num>
      # by num <num>
    }
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
    ## TRY:
    ##   fixed_it[["phi"]] <- as.vector(.unlist(fixed_it[["phi"]])) # __F I X M E___ why wasn't it necess (was it?) by compar with user_init_optim_it[["phi"]]...
    optim_blob_it <- .calc_optim_args(proc_it=unmerged[[mv_it]], processed=processed,
                                      # (__F I X M E___): intriguing heterogeneity of formats:
                                      user_init_optim=user_init_optim_it, # from .subPars + .rename_ranPars + ad hoc hack for phi... 
                                      fixed=fixed_it,  # from .subPars + .rename_ranPars; no ad hoc hack for phi
                                      lower=user_lower_it, upper=user_upper_it, # from .subPars
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
  # end loop on submodels
  if ( ! is.null(processed$init_HLfit)) {
    optim_blob$inits$init.HLfit <- .modify_list(optim_blob$inits$init.HLfit,processed$init_HLfit, obey_NULLs=FALSE)  
  }
  if ( ! is.null(optim_blob$inits$init$lambda)) {
    optim_blob <- .merge_lambdas_mv(lambda_merger, optim_blob)
  }
  if ( ! is.null(optim_blob$inits$init.optim$hyper)) {
    attr(optim_blob$inits$init.optim$hyper,"hy_info") <- processed$hyper_info # *environment*: for .makeLowerUpper, with distinct name for easier tracking 
  }
  
  if ( .has_X_off_betas(processed) &&
       length(user_init_optim$beta) ) { # outer beta... (not same workflow as in .preprocess(univar))
    X_off <- environment(processed$X_off_Xb_fn)$X_fixed
    optim_blob[["inits"]][["init"]]$beta <-  
      optim_blob[["inits"]][["init.optim"]]$beta <- user_init_optim$beta # outer beta. numeric values needed.s
    # fixef_lowup <- .calc_fixef_lowup(processed) # cannot work as such bc would need a nontrivial eta
    # (the .sanitize_eta call might to modified to avoid that).
    # Thus, currently, the user must have provided some bounds
    # which are already in user.lower and user.upper, and passed to .makeLowUp_stuff_mv()
    # by these arguments, but nevertheless currently have to be copied in:
    fixef_lowup <- list()
    # Keep in mind that for NULL beta, .scale() returns a scaled X matrix...
    if (length(user.lower$beta)) fixef_lowup$lower <- user.lower$beta 
    if (length(user.upper$beta)) fixef_lowup$upper <- user.upper$beta 
  } else fixef_lowup <- NULL
  
  optim_blob <- .makeLowUp_stuff_mv(optim_blob, user.lower=user.lower, user.upper=user.upper, 
                                    optim.scale, processed, verbose, fixef_lowup=fixef_lowup)
  
  # LUarglist contains any beta info (for outer beta optim) in $beta elements
  # this is converted to etaFix$beta by HLCor.obj() or HLfit.obj()
  # Note that .numInfo_objfn() expects a skeleton with $etaFix$beta, so 
  # elements of the optim_blob must be modified if to be used as arguments of .numInfo_objfn 
  # (____F I X M E____) is there a nicer way to control all these different computations?
  
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
      exp_barlist <- .modify_list(exp_barlist, exp_barlist_it)
      attr(exp_barlist,"type") <- .modify_list(attr(exp_barlist, 'type'), type_it)
    }
  }
  exp_barlist
}

.determine_sparse_X_mv <- function(terms_info, X.pv, 
                                   vec_nobs=attr(terms_info,"vec_nobs"), 
                                   assign.=attr(X.pv,"assign")) {
  sparse_X <- spaMM.getOption("sparse_X") 
  ## forcing sparse_X may (1) be slow for small problems 
  ## (2) entails the use of Matrix::Cholesky, which is less accurate => small bu visible effect on predVar in singular 'twolambda' case
  if (is.null(sparse_X)) {
    nmodels <- length(terms_info$Y)
    col_heuristic_densenesseS <- vector("list", nmodels)
    rel_nobs <- vec_nobs/sum(vec_nobs)
    for (mv_it in seq_len(nmodels)) {
      asgn <- assign.[[mv_it]] ## "for each column in the matrix ... the term in the formula which gave rise to the column" (O for the Intercept)
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

.designX_prod_mv <- function(merged_X, X2X, merged, assign_attr, rankinfo) {
  if (is.null(colnames(X2X))) stop("'X2X' *must* have column names")
  
  # determine sparse_X on X 'L'HS:
  sparse_X <- .determine_sparse_X_mv(merged$main_terms_info, X.pv= merged_X, assign.=assign_attr)
  
  rescale. <- ! is.null(attr(merged_X, "scaled:scale"))
  if (rescale.) merged_X <- .unscale(merged_X)
  cols_lhs_X2X <- colnames(merged_X)
  merged_X <- merged_X %*% X2X # loss of attributes, "assign" notably
  if (rescale.) merged_X <- .scale(merged_X)
  merged_X <- .post_process_X(X.pv=merged_X, HL=merged$HL, rankinfo=rankinfo, sparse_X = sparse_X) 
  attr(merged_X,"cols_lhs_X2X") <- cols_lhs_X2X # info used for predict(<mv with X2X factor>,newdata), about RHS of X2X product 
  merged_X
  
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
        is_incid <- ! ("(Intercept)" %in% namesTerm) # while is_incid was FALSE for the template...
        attr(is_incid,"is01col") <- FALSE
        attr(ZA,"is_incid") <- is_incid
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
        # Dvec is a rep(.,  rep(n_levels,Xi_ncol)) of a vector of booleans and distinguishes absence/presence of 0+
        ZA[obs.range,] <- .Matrix_times_Dvec(ZA[obs.range,], rep(as.numeric(namesTerm %in% c("(Intercept)",paste0(".mv",mv_it))),
                                                                 rep(n_levels,Xi_ncol)))
      }
      ZA <- drop0(ZA)
      is_incid <- ! ("(Intercept)" %in% namesTerm) # while is_incid was FALSE for the template...
      attr(is_incid,"is01col") <- FALSE
      attr(ZA,"is_incid") <- is_incid
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
      if (all(LMMbools) && ! length(unlist(processed$phi.Fix[which_submodels], use.names = FALSE))) {
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
    LHSstr <- .DEPARSE(exp_ranef_terms[[rd]][[2]])
    if (grepl("mv(",LHSstr, fixed=TRUE)) {
      model_ids <- sub("(mv)(\\([^|]+)",".mrdots\\2", LHSstr) # converts to an expr such as 0+c(1,2)
      model_ids <- eval(str2lang(model_ids)) # evaluates c(1,2)...
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
  if (length(umap)) {
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
  } else NULL # predVar computations check this and the alternative code would cause a bug.
}

.add_cov_matrices__from_mv_global <- function(corr_info, covStruct=NULL, corrMatrix=NULL, adjMatrix=NULL) {
  if ( ! is.null(covStruct)) covStruct <- .preprocess_covStruct(covStruct)
  if ( length(AMatrices <- attr(covStruct,"AMatrices"))) corr_info$AMatrices <- .modify_list(corr_info$AMatrices, AMatrices)
  corr_types <- corr_info$corr_types
  for (it in seq_along(corr_types)) {
    corr_type <- corr_types[[it]]
    if ( ! is.na(corr_type)) {
      if (corr_type=="adjacency" || corr_type=="SAR_WWt") {
        if ( is.null(adjMatrix) ) adjMatrix <- .get_adjMatrix_from_covStruct(covStruct,it)
        if (is.null(adjMatrix)) {
          # should probably check that the adjMatrix is already in. (_FIXME_)
        } else corr_info$adjMatrices[[it]] <-  .reformat_adjMatrix(adjMatrix)
      } else if (corr_type=="corrMatrix") {
        if (is.null(corrMatrix)) corrMatrix <- .get_corr_prec_from_covStruct(covStruct,it, required=FALSE) 
        if ( is.null(corrMatrix)) {
          # should probably check that the corrMatrix is already in. (_FIXME_)
        } else corr_info$corrMatrices[[it]] <- corrMatrix
        .check_corrMatrix(corr_info$corrMatrices[[it]], element=1) 
      } 
      # IMRF AMatrices are assigned later from Zlist info, not from covStruct info
    }
  }
}


.update_ZAlist <- function(ZAlist, AMatrices, which_ZA) { 
  if ( length(ZAlist)) {
    for (char_rd in names(ZAlist)[which_ZA]) { 
      Amatrix <- AMatrices[[char_rd]]
      if ( ! is.null(Amatrix)) ZAlist[[char_rd]] <- .ZxA_with_attrs(Z=ZAlist[[char_rd]], Amatrix) 
    }
  } 
  return(ZAlist)
}


.preprocess_betaFix <- function(X.pv, betaFix, merged) { 
  namesbetafix <- names(betaFix)
  if (is.null(namesbetafix)) {
    message("The elements of etaFix$beta should be named and the names should match the column names of the design matrix.")
  }
  if (length(setdiff(namesbetafix,colnames(X.pv)))==0L) { ## if no incorrect name
    betaFix <- .scale(X = X.pv, beta = betaFix)
    merged$X_off_Xb_fn <- .def_off_Xb_fn(X_fixed=.subcol_wAttr(X.pv, j=namesbetafix, drop=FALSE), 
                                       offsets=merged$off,
                                       sc_betaFix=betaFix)  # X.pv is scaled here so betaFix must be too
    merged$off <- merged$X_off_Xb_fn(NULL,NULL,NULL)
    namesbetavar <- setdiff(colnames(X.pv),namesbetafix)
    X.pv <- .subcol_wAttr(X.pv, j=namesbetavar, drop=FALSE)
  } else {
    stop("The names of elements of etaFix$beta should all match column names of the design matrix.")
  }
  X.pv
}

.merge_processed <- function(calls_W_processed, data, init=list(), control.HLfit=list(), method="ML", verbose=NULL, init.HLfit=list(),
                             covStruct=NULL, corrMatrix=NULL, adjMatrix=NULL, distMatrix=NULL, control.dist=list(), etaFix=list(),
                             X2X=NULL) {
  # this fn passes no '...' so has no '...'
  nmodels <- length(calls_W_processed)
  namedlist <- structure(vector("list",nmodels), names=seq_len(nmodels))
  ### Fill lists for further processing:
  unmerged <- predictors <- families <- prior.weights <- clik_fns <- phiFixs <- Ys <- pS_fixef_phi <- 
    assign_attr <- rankinfo_attr <- namedlist
  AMatrices <- adjMatrices <- corrMatrices <- fixef_off_termsS <- fixef_termsS <- fixef_levelsS <- 
    specials_levelsS <- validrownames <- namedlist
  vec_nobs <- integer(nmodels)

  # Two big loops: the second recursively updates in more or less complex ways 
  # elements initialized from unmerged[[1L]], 
  #
  # while the first, here, is a simpler collection of elements from all unmerged's: 
  for (mv_it in seq_len(nmodels)) {
    p_i <- unmerged[[mv_it]] <- calls_W_processed[[mv_it]][["processed"]]
    validrownames[[mv_it]] <- rownames(p_i[["data"]])
    vec_nobs[mv_it] <- length(p_i[["y"]])  
    predictors[[mv_it]] <- p_i[["predictor"]]
    families[[mv_it]] <- p_i[["family"]]
    prior.weights[[mv_it]] <- p_i[["prior.weights"]] ## may need to be quoted expression, etc.
    clik_fns[[mv_it]] <- p_i[["clik_fn"]]
    phiFixs[mv_it] <- list(p_i[["phi.Fix"]]) # ! syntax to allow explicit NULL's
    pS_fixef_phi[[mv_it]] <- p_i[["p_fixef_phi"]] 
    assign_attr[mv_it] <- list(attr(p_i[["AUGI0_ZX"]]$X.pv,"assign")) ## "for each column in the matrix ... the term in the formula which gave rise to the column"
    rankinfo_attr[mv_it] <- list(attr(p_i[["AUGI0_ZX"]]$X.pv,"rankinfo")) 
    
    terms_info_it <- p_i[["main_terms_info"]]
    Ys[[mv_it]] <- terms_info_it[["Y"]] # one reason for parsing frames that way is to allow $Y as argument of .get_inits_from_glm()
    fixef_off_termsS[[mv_it]] <- terms_info_it[["fixef_off_terms"]]
    fixef_termsS[mv_it] <- list(terms_info_it[["fixef_terms"]]) 
    fixef_levelsS[mv_it] <- list(terms_info_it[["fixef_levels"]]) 
    specials_levelsS[mv_it] <- list(terms_info_it[["specials_levels"]]) 
  }
  
  merged <- list2env(list(envir=list2env(list(), parent=environment(HLfit))))
  
  main_terms_info <- list(
    Y=Ys, 
    fixef_off_terms=fixef_off_termsS,
    fixef_terms=fixef_termsS,
    fixef_levels=fixef_levelsS,
    specials_level=specials_levelsS) # (not yet ready for inclusion in 'merged') 
  merged[["main_terms_info"]] <- structure(main_terms_info, vec_nobs=vec_nobs)
  cum_nobs <- cumsum(c(0L, vec_nobs)) 
  merged[["vec_nobs"]] <- structure(vec_nobs, cum_nobs=cum_nobs) # object will contain multiple copies of cum_nobs attribute, for conveniency
  merged[["For"]] <- "fitme" # does not appear necessary except with adjmatrixas $For needed to determine inner_estim_adj_rho in .determine_spprec()
  merged[["phi.Fix"]] <- phiFixs
  merged[["p_fixef_phi"]] <- pS_fixef_phi
  merged$predictor <- predictors # that will be passed by HLfit_body to the result...
  merged$clik_fn <- clik_fns
  
  ### initialize some 'merged' elements from unmerged[[1L]]:
  #
  # This small loop for elements *not* recursively updated later, so directly assigned to envir=merged:
  for (st in c(#"AUGI0_ZX",
    "REMLformula", # __FIXME__ this will really handle only standard ML (REMLformula has an isML attr) 
    #                                                      or standard REML (REMLformula is NULL). 
    # => no attempt to look REMLformula over models below (But  is built iteratively). 
    "verbose","control.glm","HL","p_v_obj",#"rand.families",
    "spaMM_tol",
    "break_conv_logL","intervalInfo",
    "objective","port_env","CONTROL")
  ) assign(st,value=unmerged[[1L]][[st]],envir=merged)
  #
  # This small loop for elements recursively updated later, so not yet assigned into 'merged':
  for (st in c("off","y","BinomialDen","iter_mean_dispVar","iter_mean_dispFix",
               "max.iter","models","vecdisneeded","bin_all_or_none")) assign(st,value=unmerged[[1L]][[st]])

  # Operations on 'models' form the first submodel:
  LMMbool  <- attr(models,"LMMbool")
  GLMMbool  <- attr(models,"GLMMbool")
  LLM_const_w  <- attr(models,"LLM_const_w")
  GLGLLM_const_w  <- attr(models,"GLGLLM_const_w")
  GLMbool <- (models[["eta"]]=="etaGLM")
  LMbool <- unmerged[[1L]]$family$flags$LMbool
  unit_GLMweights  <- attr(models,"unit_GLMweights")
  unit_Hobs_weights  <- attr(models,"unit_Hobs_weights")
  const_Hobs_wresid  <- attr(models,"const_Hobs_wresid")
  phi_models <- rdispar_models <- character(nmodels)
  phi_models[[1L]] <- models[["phi"]]
  rdispar_models[[1L]] <- models[["rdispar"]]
  #
  ZAlist <- unmerged[[1L]]$ZAlist # misnomer, as A woill be included much later (cf .update_ZAlist() call)
  ZAlist <- .correct_ZA_mv_ranCoefs(ZAlist, mv_it=1L)
  ZAlist <- .merge_ZAlists(list(), ZAlist, 0L, vec_nobs[1L], 1L)
  merged_X <- .merge_Xs(NULL, unmerged[[1L]][["AUGI0_ZX"]]$X.pv, mv_it=1L)
  merged_X.Re <- .merge_Xs(NULL, unmerged[[1L]][["X.Re"]], mv_it=1L, REML=TRUE)
  vec_ncol_X <- integer(nmodels)
  vec_ncol_X[1L] <- ncol(merged_X)
  validrownames[[1L]] <- rownames(unmerged[[1L]][["data"]])
  obsInfo <- unmerged[[1L]][["how"]][["obsInfo"]]
  # Recursive updating:
  for (mv_it in (seq_len(nmodels-1L)+1L)) {
    p_i <- unmerged[[mv_it]]
    # merged$predictor... I should try to get rid of its use in HLfit_body, at least... FIXME
    ### random effects stuff
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
    LMbool_it <- p_i$family$flags$LMbool
    modattrs_it <- attributes(p_i[["models"]])
    unit_GLMweights <- unit_GLMweights && modattrs_it[["unit_GLMweights"]]
    unit_Hobs_weights <- unit_Hobs_weights && modattrs_it[["unit_Hobs_weights"]]
    LLM_const_w_it <- modattrs_it[["LLM_const_w"]]
    GLGLLM_const_w_it <- modattrs_it[["GLGLLM_const_w"]] 
    LMMbool_it <- modattrs_it[["LMMbool"]]
    GLMMbool_it <- modattrs_it[["GLMMbool"]]
    const_Hobs_wresid_it <- modattrs_it[["const_Hobs_wresid"]]
    LLM_const_w <- ((LLM_const_w  || (LMbool && const_Hobs_wresid) ) && LLM_const_w_it) ||
      ( LLM_const_w && LMbool_it && const_Hobs_wresid_it )
    GLGLLM_const_w <- ((GLGLLM_const_w  || (GLMbool && const_Hobs_wresid) ) && GLGLLM_const_w_it) ||
      (GLGLLM_const_w && GLMbool_it && const_Hobs_wresid_it ) # (constant weights over iterations of IRLS)
    LMMbool <- ((LMMbool || LMbool) && LMMbool_it) ||
      (LMMbool && LMbool_it)
    GLMMbool <- ((GLMMbool || GLMbool) && GLMMbool_it) ||
      (GLMMbool && GLMbool_it)
    const_Hobs_wresid  <- const_Hobs_wresid && const_Hobs_wresid_it
    GLMbool <- GLMbool && GLMbool_it
    LMbool <- LMbool && LMbool_it
    # LLFbool <- LLFbool || modattrs_it[["LLFbool"]]
    obsInfo <- obsInfo || p_i[["how"]][["obsInfo"]] # it beign TRUE only for 'distinctly obsInfo' submodels
    iter_mean_dispVar <- max(iter_mean_dispVar, p_i[["iter_mean_dispVar"]])
    iter_mean_dispFix <- max(iter_mean_dispFix, p_i[["iter_mean_dispFix"]])
    max.iter <- max(max.iter, p_i[["max.iter"]])
    vecdisneeded <- (vecdisneeded | p_i[["vecdisneeded"]])
    phi_models[[mv_it]] <- p_i[["models"]][["phi"]]
    rdispar_models[[mv_it]] <- p_i[["models"]][["rdispar"]]
  }
  .check_mv_in_submodels(ZAlist)
  
  if (identical(control.HLfit$getFromTheDepths,"surrogate_info")) { # ____F I X M E____ assemble the info much earlier in this fn?
    has_dynoffset <- sapply(lapply(predictors,.DEPARSE),
                            grepl, pattern="offset(.dynoffset)", fixed=TRUE)
    surrogate_info <- list(
      has_dynoffset=has_dynoffset, 
      validPrownames=validrownames, # list of vectors of strings
      # is_fixefM=length(ZAlist)==0L,
      vec_nobs=vec_nobs # valid for poisson fit; cf .get_multinom_info() for use.
    ) 
    .sendFromTheDepths(surrogate_info=surrogate_info, class="spaMM.FTD.SI")
  }
  
  
  attr(data,"validrownames") <- validrownames
  merged[["data"]] <- data
  merged[["prior.weights"]] <- prior.weights
  merged[["off"]] <- off
  merged[["bin_all_or_none"]] <- bin_all_or_none 
  merged[["iter_mean_dispVar"]] <- iter_mean_dispVar
  merged[["iter_mean_dispFix"]] <- iter_mean_dispFix
  merged[["max.iter"]] <- max.iter
  merged[["vecdisneeded"]] <- vecdisneeded
  merged[["how"]] <- list(obsInfo=obsInfo)
  # 
  augZXy_cond_inner <- spaMM.getOption("allow_augZXy")
  if (is.null(augZXy_cond_inner)) augZXy_cond_inner <- TRUE ## for .makeCovEst1()
  if (augZXy_cond_inner) augZXy_cond_inner <-   LMMbool
  if (augZXy_cond_inner) augZXy_cond_inner <- ( is.null(merged_X.Re) || ! ncol(merged_X.Re)) ## exclude only non-standard REML 
  merged$augZXy_cond <- structure(FALSE, inner=augZXy_cond_inner) # augZXy_cond would impose a unique phi across submodels
  #
  # attr(phi_models,"anyHGLM") <- any(phi_models=="phiHGLM") # not used
  models[["phi"]] <- phi_models # 'models' is a list whose element 'phi' is a vector
  models[["rdispar"]] <- rdispar_models # 'models' is a list whose element 'rdispar' is a vector
  attr(models, "LMMbool") <- LMMbool # add more attributes to avoid clumsy tests later
  attr(models,"GLMMbool") <- GLMMbool 
  attr(models,"LLM_const_w") <- LLM_const_w 
  attr(models,"GLGLLM_const_w") <- GLGLLM_const_w 
  attr(models,"unit_GLMweights") <- unit_GLMweights
  attr(models,"unit_Hobs_weights") <- unit_Hobs_weights
  attr(models,"const_Hobs_wresid") <- const_Hobs_wresid
  merged[["models"]] <- models
  models <- NULL # make sure we work on only one 'models'
  if (merged$how$obsInfo) { # at least one 'distinctly obsInfo' submodel
    for (mv_it in seq_along(families)) if ( ! families[[mv_it]]$flags$LLgeneric) families[[mv_it]]  <- .statsfam2LLF(families[[mv_it]])
    merged$etaxLM_fn <- .calc_etaLLMblob
  } else merged$etaxLM_fn <- .calc_etaGLMblob
  
  residProcesseds <- residModels <- vector("list", nmodels)
  for (mv_it in seq_len(nmodels)) { 
    residModels[mv_it] <- list(unmerged[[mv_it]]$residModel) # (allows list(NULL) but generally non-NULL even if trivial;  
                                                             # *!* for non-GLM the resid.model is in the family closure environment)
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
    has_estim_families_par <- (
      (family_it$family %in% c("negbin1","negbin2") && 
         .is_fampar_missing(family=family_it,fampar="shape")) ||
        (family_it$family %in% c("beta_resp","betabin") && 
           .is_fampar_missing(family=family_it,fampar="prec"))||
        (family_it$family=="COMPoisson" && 
           .is_fampar_missing(family=family_it,fampar="nu"))||
        (family_it$family=="tweedie" && 
           (.is_fampar_missing(family=family_it,fampar="p") ||
              .is_fampar_missing(family=family_it,fampar="q")))
    )
    if (has_estim_families_par) break
  }
  attr(families,"has_estim_families_par") <- has_estim_families_par 
  attr(families,"famfams") <- sapply(families, `[[`, x="family")
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
  
  col_ranges <- vector("list", nmodels)
  cum_ncol_X <- cumsum(c(0L,vec_ncol_X))
  
  for (mv_it in seq_len(nmodels)) col_ranges[[mv_it]] <- cum_ncol_X[mv_it]+seq_len(cum_ncol_X[mv_it+1L]-cum_ncol_X[mv_it]) # avoid (n+1:n) problem

  if  (! is.null(X2X)) {
    if (inherits(X2X,"call")) {
      X2X[["names_ori"]] <- colnames(merged_X)
      X2X <- eval(X2X) 
    } # else if (identical(X2X,"coef_names")) .sendFromTheDepths(coef_names=colnames(merged_X))
    
    merged_X <- .designX_prod_mv( merged_X, X2X, merged, assign_attr, rankinfo=control.HLfit$rankinfo)
    abs_X_RHS <- abs(X2X)
    for (mv_it in seq_len(nmodels)) col_ranges[[mv_it]] <- which(colSums(abs_X_RHS[col_ranges[[mv_it]],,drop=FALSE])>0)
  } else {
    merged_X <- .post_process_X(X.pv=merged_X, HL=merged$HL, rankinfo=control.HLfit$rankinfo, 
                                sparse_X=.determine_sparse_X_mv(merged$main_terms_info, X.pv= merged_X, 
                                                                assign.=assign_attr) ) 
    # processing of merged_X and other elements of AUGI0_ZX:
  }
  attr(merged_X,"assign") <- assign_attr       # notably used by anova()
  attr(merged_X,"rankinfo")$mvlist <- rankinfo_attr   # notably used by anova()
  attr(merged_X,"col_ranges") <- col_ranges
  attr(merged_X,"cum_nobs") <- cum_nobs
  
  #### replacement for .preprocess_X_XRe_off():
  if ( length(betaFix <- etaFix$beta)>0 ) { # .preprocess_betaFix updates {any $off previously created by merging submodel $off's}
    if (is.null(merged_X.Re)) { # standard REML were it not for betaFix:
      # in standard REML the correction is based on the cols of X.pv (were it not for betaFix) 
      keepInREML <- attr(betaFix,"keepInREML") 
      if (is.null(keepInREML)) keepInREML <- FALSE
      if (keepInREML) merged_X.Re <- merged_X # keep cols of full X.pv in X.Re (as in std REML without betaFix)
      merged_X <- .preprocess_betaFix(X.pv=merged_X, betaFix=betaFix, merged)
      if ( ! keepInREML) {
        merged_X.Re <- merged_X # keeps only cols of sub X in X.pv
      } else { 
        # (keepInREML TRUE)=> merged_X has been modified by .preprocess_betaFix, but we don't update 'XReinput', so the two are now different, 
        # 'XReinput' being the merged_X before .preprocess_betaFix -> .subcol_wAttr()ing
        # 'XReinput' is the name in .preprocess() of what is here merged_X.Re. Maybe tidy .preprocess()? 
      }
    } else merged_X <- .preprocess_betaFix(X.pv=merged_X, betaFix=betaFix, merged)
    merged[["vecdisneeded"]] <- merged[["vecdisneeded"]] & ncol(merged_X)
  }
  #
  ### Following block replaces, at least partially,
  #   .assign_X.Re_objective(processed, XReinput=XReinput, processed$REMLformula, data, X.pv, objective) ##] assigns processed$X.Re and processed$objective
  ## in mv fit, the merged 'objective' is set by copying the one of the first submodel, and X.Re is first set by merging submodels' X.Re's.
  ## In principle this should be modified following the above modif of X.pv. HOWEVER:
  ## if standard ML: ... processed$X.Re is 0-col matrix: unchanged
  ## if standard REML: processed$X.Re is NULL: unchanged
  ## non standard REML: case where X.Re would have to be modified, per the following block:
  ##    which may handle correctly only a few subcases of non-std REML
  ##    I wrote it before trying to implement the keepInREML case in this fn, so it must havehad some previous usage...
  ##    In the the keepInREML case, merged.X.Re retains columns removed from merged_X.
  ##    The REMLformula, needed to handle more general cases, is not defined here.
  if ( ! is.null(merged_X.Re) && ncol(merged_X.Re)) { # non-standard REML... 
    # attr(., "extra_vars") has aleardy been updated by .merge_Xs(., REML=TRUE), but unrestricting_cols culd not as it refers to X.pv cols 
    unrestricting_cols <- which(colnames(merged_X) %in% setdiff(colnames(merged_X),colnames(merged_X.Re))) ## not in X.Re
    distinct.X.ReML <- c(length(unrestricting_cols), length(attr(merged_X.Re,"extra_vars"))) ## TWO integers often used as booleans 
    attr(merged_X.Re,"distinct.X.ReML") <- distinct.X.ReML 
    if ( distinct.X.ReML[1L]) attr(merged_X.Re,"unrestricting_cols") <- unrestricting_cols # cols of X.pv not in X.Re
  }
  merged[["X.Re"]] <- merged_X.Re
  ####
  
  if ( length(init_beta <- init[["beta"]])) { # outer beta. See above for fixed beta
    # The distinction by this test is between, say, 
    # * fitmv(., etaFix=list(beta=.)):    [ case if length(betaFix <- etaFix$beta) above ]
    #    an offset could be set at preprocessing time, there is no need for an X_off_Xb_fn and no final inner 'refit' of beta is considered.  
    # * fitmv(., init=list(beta=.)): outer optimization of beta
    #    an X_off_Xb_fn is used to compute  X_off %*% beta for each new beta, 
    #    and a final inner 'refit' of beta is possible => distinct 'vecdisneeded_ori' to be used then.  
    ## In both cases, cols of merged_X are removed, so the two are not compatible (at least for the same coefficients).
    betanames <- names(init_beta)
    if (length(intersect(colnames(merged_X),betanames))!=length(init_beta)) stop("init[['beta']] must have names matching those of the design matrix")
    X_off <-.subcol_wAttr(merged_X, j=betanames, drop=FALSE)
    merged_X <- .subcol_wAttr(merged_X, j=setdiff(colnames(merged_X),betanames), drop=FALSE)
    merged$X_off_Xb_fn <- .def_off_Xb_fn(X_fixed=X_off, offsets=merged$off)
    merged[["vecdisneeded_ori"]] <-  merged[["vecdisneeded"]]
    merged[["vecdisneeded"]] <- merged[["vecdisneeded"]] & ncol(merged_X)
  }
  
  dim(y) <- c(length(y),1L) # cf comments in .preprocess(); but vector format worked in all tests for fitmv up to version 3.5.115; 
  merged[["y"]] <- y
  merged[["BinomialDen"]] <- BinomialDen
  attr(ZAlist,"map_rd_mv") <- map_rd_mv <- .map_rd_mv(ZAlist, unmerged)
  #
  # map_rd_mv available -> exp_ranef_types -> soon used by .assign_corr_types_families()
  exp_barlist <- .process_bars_mv(predictors, map_rd_mv)
  exp_ranef_strings <- .process_bars(barlist=exp_barlist,expand=FALSE, as_character=TRUE) ## no need to expand again
  if (nrand <- length(exp_ranef_strings)) {
    # ZAlist including map_rd_mv available -> rand.families
    merged[["rand.families"]] <- rand.families <- .merge_rand_families(unmerged, ZAlist=ZAlist) 
    merged$lcrandfamfam <- attr(rand.families,"lcrandfamfam") ## else remains NULL
    
    names(ZAlist) <- seq_len(nrand) # not sure where it is used, but conform to fact that univar-resp ZAlist is named.
    exp_ranef_types <- attr(exp_ranef_strings,"type") ## expanded
    attr(ZAlist,"exp_ranef_strings") <- exp_ranef_strings ## expanded 
    attr(ZAlist,"exp_ranef_types") <- exp_ranef_types ## expanded
    #
    ####### Builds $corr_info (ASAP to assign_cov_matrices ASAP):
    merged$corr_info <- corr_info <- new.env() ## do not set parent=emptyenv() else with(corr_info,...) will not find trivial fns such as `[`
    covStruct <- .assign_corr_types_families(covStruct=covStruct, corr_info=corr_info, exp_ranef_types=exp_ranef_types, 
                                             exp_barlist=exp_barlist) # provides 'corr_types' and 'is_cF_internally', soon necessary, and 'corr_families'
    
    # either I assign the A matrices within each submodel and I merge them across submodels afterwards,
    # or I assign them on the merged "corr_info". In neither case .assign_AMatrices_corrFamily(corr_info, ...) is useful *here*.
    
    
    ## Some corrFamily stuff (next TWO loops)
    levels_types <- rep("", nrand)
    for (rd in which(corr_info$is_cF_internally)) { # (allows NA in $corr_types)
      # For the hard coded Matern(), AR1() etc. $corr_families[[it]] is already a list of functions $calc_moreargs, $canonize...
      # for corrFamily() by the next line it will be the corrFamily descriptor as an *environment* with $f, $tpar, $type, $template, $Af... $calc_moreargs, $canonize...
      corr_info$corr_families[[rd]] <- .preprocess_corrFamily(corrfamily=eval(covStruct[[rd]])) 
      levels_types[[rd]] <- corr_info$corr_families[[rd]]$levels_type 
      # For the hard coded Matern(), AR1() etc. $corr_families[[it]] is already a list of functions $calc_moreargs, $canonize...
      # for corrFamily() by the next line it will be the corrFamily descriptor as an *environment* with $Cf, $tpar, $type, $template, $Af... $calc_moreargs, $canonize...
      .initialize_corrFamily(corr_info$corr_families[[rd]], Zmatrix=ZAlist[[rd]])
    }
    corr_info$levels_types <- levels_types # so that this info can be looked up efficiently post-fit
    
    # copy some of the info in each submodel's corr_info
    for (mv_it in seq_along(unmerged)) {
      which_in_subm <- which(map_rd_mv[[mv_it]] %in% corr_info$is_cF_internally ) # so eg from (c("1"=1,"2"=3)) %in% 3), rd_in_mv is 2 
      pos_in_sub <- map_rd_mv[[mv_it]][which_in_subm]
      pos_in_merged <- as.integer(names(pos_in_sub)) 
      unmerged[[mv_it]]$corr_info$corr_families[pos_in_sub] <- corr_info$corr_families[pos_in_merged]
      # there are hidden subcases: if submodel is not MM, LHS's corr_families is NULL before and after the assignment (and RHS is 'list()')
      ##  Allow distinct maxLambda for each ranef according to the response links of the submodel(s) within which it is involved.
    }
    ## /corrFamily
    #
    for (mv_it in seq_len(nmodels)) {
      corr_info_it <- unmerged[[mv_it]][["corr_info"]]     # already preprocessed info for:
      AMatrices[mv_it] <- list(corr_info_it$AMatrices)
      adjMatrices[mv_it] <- list(corr_info_it$adjMatrices)
      corrMatrices[mv_it] <- list(corr_info_it$corrMatrices)
      # corr_info_it also has ""corr_families" "corr_types""cov_info_mats" "G_diagnosis"
    }
    # follow the order of .preprocess():
    ## Assigns $corr_info list of matrices BEFORE determining sparse precision:
    corr_info$corrMatrices <- .merge_mv_list(corrMatrices, merged, ZAlist=ZAlist, full=TRUE) 
    corr_info$adjMatrices <- .merge_mv_list(adjMatrices, merged, ZAlist=ZAlist, full=TRUE)
    corr_info$AMatrices <- .merge_mv_list(AMatrices, merged, ZAlist=ZAlist, full=TRUE)
    #
    for (rd in which(corr_info$is_cF_internally)) { 
      .assign_AMatrices_corrFamily(corr_info, ZAlist, exp_barlist=exp_barlist, merged$data, control_dist=merged$control_dist)
    }
    
    ## for (mult)IMRF:
    ## HLCor_body -> .assign_geoinfo_and_LMatrices_but_ranCoefs() expects corr_info$AMatrices
    ## Otherwise it may stop on  .............................................. -> .calc_IMRF_Qmat()
    ## In .preprocess they are added through ZAlist <- .calc_ZAlist(Zlist=Zlist, AMatrices=corr_info$AMatrices),
    ##    using corr_info$AMatrices previously assigned by
    ##    .assign_AMatrices_IMRF(corr_info, Zlist, exp_barlist=exp_barlist, processed$data, control_dist=processed$control_dist)
    ## (and later addition: .assign_AMatrices_corrFamily() )
    #
    ## So we need first to replicate the effect on corr_info$AMatrices of 
    ##    .assign_AMatrices_IMRF(corr_info, Zlist, exp_barlist=exp_barlist, processed$data, control_dist=processed$control_dist)
    ## Then later to replicate the effect of ZAlist <- .calc_ZAlist(Zlist=Zlist, AMatrices=corr_info$AMatrices) on corr_info$AMatrices
    .add_cov_matrices__from_mv_global(corr_info, covStruct=covStruct, corrMatrix=corrMatrix, adjMatrix=adjMatrix) # using $corr_info$corr_type
                                     #  which modifies corr_info$AMatrices if attr(covStruct,"AMatrices") is present
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
    # if (any(exp_ranef_types== "corrFamily") && is.null(control.HLfit$algebra)) {
    #   message("No control.HLfit$algebra specified: it is set to 'spcorr' by default.") 
    #   control.HLfit$algebra <- "spcorr"
    # } 
    #
    ## for corrFamily:
    # The AMatrices are deduced from the the covStruct argument, not from the formula terms of the submodels 
    # => they are not yet factored in ZAlist, although they are available from the above call to .assign_AMatrices_corrFamily
    ## Same problem for user-provided A matrix, whatever the corr_type (which may be NA).
    ## BUT "IMRF" is a "special ranef" with ad hoc code, not handled by the general corrFamily
    # interface, so it gets a special handling here (in contrast to even MAternIMRFa).
    # See related comment on .calc_ZAlist() call in .preprocess().
    is_IMRF_ss <- (attr(ZAlist, "exp_ranef_types") %in% c("IMRF")) 
    # test this rather than corr_info$corr_types which contains NA for (.|.) (incl. ranCoefs) 
    # and would then be ignored by which().
    # At this point, for MaternIMRFa() the ZA product has not yet been performed,
    # But for IMRF(), it has.
    ZAlist <- .update_ZAlist(ZAlist, AMatrices=corr_info$AMatrices, 
                             which_ZA=which( ! is_IMRF_ss )) # IMRFs are separately handled below
    # ____F I X M E____ it would be nice to be able to check that all A matrices are taken into account.
    # Using corr_info:
    merged$control_dist <- .preprocess_control.dist(control.dist, corr_info$corr_types)
    #
    .check_init.HLfit(init.HLfit)
    merged$init_HLfit <- .preprocess_init.HLfit(init.HLfit, corr_info)
  
    merged$is_spprec <- .wrap_determine_spprec(control.HLfit, ZAlist=ZAlist, processed=merged, X.pv=merged_X)
    ## post-processing of corr_info depending on sparse_precision
    .process_corr_info_spprec(corr_info=corr_info, For="fitme",sparse_precision=merged$is_spprec)
    # Heavily using corr_info, and spprec:
    ZAlist <- .init_assign_geoinfo(processed=merged, ZAlist=ZAlist, For="fitme", 
                                   exp_barlist=exp_barlist, distMatrix=distMatrix)  
    vec_normIMRF <- .calc_vec_normIMRF(exp_ranef_terms=attr(ZAlist, "exp_ranef_terms"), corr_info=corr_info)   
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
    merged[["ZAlist"]] <- ZAlist
    merged$hyper_info <- .merge_hyper_infos(ZAlist, unmerged)
  } else {
    merged$is_spprec <- FALSE ## for .do_TRACE()
    merged$init_HLfit <- init.HLfit
  }
  if (merged$is_spprec) {
    merged$solve_IRLS_fn <- .solve_IRLS_as_spprec
  } else  merged$solve_IRLS_fn <- .solve_IRLS_as_ZX
  
  merged$QRmethod <- .choose_QRmethod(ZAlist, corr_info=merged$corr_info, 
                                      is_spprec=merged$is_spprec, processed=merged, control.HLfit=control.HLfit)
  algebra <- .set_augX_methods(merged) # sets processed$sXaug_method nd $spprec method
  .check_time_G_diagnosis(.provide_G_diagnosis, processed=merged, algebra)
  #
  thread_nbr <- control.HLfit$NbThreads
  if (is.null(thread_nbr)) thread_nbr <- .spaMM.data$options$NbThreads # should be 1L by default
  
  if (nrand) {
    merged$models[["eta"]] <- "etaHGLM" 
    vec_n_u_h <- unlist(lapply(merged$ZAlist,ncol)) 
    merged[["psi_M"]] <- rep(attr(rand.families,"unique.psi_M"),vec_n_u_h)
    merged[["cum_n_u_h"]] <- cumsum(c(0L, vec_n_u_h))
    nrd <- merged[["cum_n_u_h"]][nrand+1L]
    if (algebra=="decorr" && nrd>900L && thread_nbr==1L) message(cli::format_message('Using parallelisation might be useful. See {.topic [setNbThreads](spaMM::setNbThreads)}'))
    if (nrd==1L) {
      warning("Found a single random effect with a *single level*. Check formula?", immediate.=TRUE)
    }
    merged$models[["lambda"]] <- rep("lamScal",nrand) ## even for adjacency, random slope...
    merged[["reserve"]] <- .preprocess_arglists(merged)
    #
    maxLambda <- rep(.spaMM.data$options$maxLambda, nrand) 
    for (mv_it in seq_along(unmerged)) {
      maxLambda[map_rd_mv[[mv_it]]] <- pmin(maxLambda[map_rd_mv[[mv_it]]], unmerged[[mv_it]][["maxLambda"]]) 
    }
    merged[["maxLambda"]] <- maxLambda
    #
    # Need this as long as we need  .calc_optim_args_mv(processed$unmerged....) ie as long as  .calc_optim_args() handles only single $family 
    for (mv_it in seq_along(unmerged)) {
      rd_in_mv <- map_rd_mv[[mv_it]]
      unmerged[[mv_it]]$geo_info <- merged$geo_info[rd_in_mv]
    }
  } else merged$models[["eta"]] <- "etaGLM" 
  .check_nb_cores(nb_cores = thread_nbr)
  .setNbThreads(thr=thread_nbr)
  
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
  merged$HLfit_body_fn <- .spaMM.data$options$HLfit_body
  .do_TRACE(merged)
  merged$HLfit_body_fn <- get(merged$HLfit_body_fn, asNamespace("spaMM"), inherits=FALSE) 
  delayedAssign("HLCor_body", get("HLCor_body", asNamespace("spaMM"), inherits=FALSE), assign.env = merged) 
  delayedAssign("HLCor", get("HLCor", asNamespace("spaMM"), inherits=FALSE), assign.env = merged) 
  merged$HLfit <- get("HLfit", asNamespace("spaMM"), inherits=FALSE) 
  merged$fitenv <- list2env(list(prevmsglength=0L))
  class(merged) <- c("mvarglist", class(merged))
  return(merged)
  # models <- NULL # prevents a R CMD check NOTE: .merge_processed: no visible binding for global variable 'models' (!)
  # Apparently the check is satisfied only if any variable is assigned  by recognized means, 
  # including '<-' but not assign(). Even this silly placement of '<-' works.
  # (For other variables such as 'off' created by assign(), off <- ....(off) plays this role)
}

# does not canonize (to canonical scale) but formats corrPars to structured list and lambda to standardly-named vector.
.reformat_parlist <- function(parlist, processed) { # for merging user's inputs mid-processing. Not for optimization hence no transformed params
  parlist$lambda <- .reformat_lambda(parlist$lambda, processed=processed, full_lambda=FALSE)
  parlist$phi <- .reformat_phi(parlist$phi,n_models=length(processed$unmerged), full_phi=FALSE)
  # This fn is only for fixed values and fixed COMP_nu is best given by the family argument...   
  # parlist$COMP_nu <- .reformat_phi(parlist$COMP_nu,n_models=length(processed$unmerged), full_phi=FALSE)
  # parlist$NB_shape <- .reformat_phi(parlist$NB_shape,n_models=length(processed$unmerged), full_phi=FALSE)
  parlist <- .reformat_corrPars(parlist,corr_families=processed$corr_info$corr_families)
}

.wrap_diff_str <- function(v, mismatch) {
  dif <- cli::diff_str(mismatch, v)
  # lcs <- dif$lcs
  # dif$lcs <- lcs[lcs$operation != "delete",]
  # dif
  format(dif)
}

# Functions derived from cli:::format_diff_str_color
.format_invalid <- function (x, ...) { # for 1st invalid name
  out <- lapply(seq_len(nrow(x$lcs)), function(i) {
    op <- x$lcs$operation[i]
    off <- x$lcs$offset[i]
    len <- x$lcs$length[i]
    if (op == "match") {
      paste0(x$old[off + 1:len], collapse = "")
    }
    else if (op == "delete") {
      cli::bg_red(cli::col_black(paste0(x$old[off + 1:len], collapse = "")))
    }
    else if (op == "insert") {
      cli::col_red(paste0(rep("\U00B7",len), collapse = ""))
    }
  })
  paste(out, collapse = "")
}

.format_valid <- function (x, ...) {
  out <- lapply(seq_len(nrow(x$lcs)), function(i) {
    op <- x$lcs$operation[i]
    off <- x$lcs$offset[i]
    len <- x$lcs$length[i]
    if (op == "match") {
      paste0(x$old[off + 1:len], collapse = "")
    }
    else if (op == "delete") {
      paste0(rep("\U00B7",len), collapse = "")
    }
    else if (op == "insert") {
      paste0(x$new[off + 1:len], collapse = "")
    }
  })
  paste(out, collapse = "")
}


.diagnose_string_mismatches <- function(extracols, obs_coln_LHS) {
  submodelstr <- function(v) gsub(pattern="(.+?)_([0-9]+)$",replacement = "\\2", x=v)
  # find first submodel with extracols
  extrasubmodels <- suppressWarnings(as.integer(sapply(extracols, submodelstr))) 
  # : suppresses ...NAs introduits lors de la conversion automatique
  if (anyNA(extrasubmodels)) {
    message(paste0("submodel indices appear to be missing from name(s): '",
                   paste(extracols[is.na(extrasubmodels)], collapse="', '"),"'"))
    return(NULL)
  }
  extra1stsub <- min(extrasubmodels)
  # find names implied by this submodel formula
  names1stsub <- obs_coln_LHS[as.integer(sapply(obs_coln_LHS, submodelstr))==extra1stsub]
  # first mismatches in first suspect submodel:
  mismatch1 <- extracols[extrasubmodels==extra1stsub][1]
  if (length(names1stsub)) {
    mess <- cli::format_message(paste0("Differences of first ",
                                       format(cli::diff_str("invalid",""))," name '", mismatch1, 
                                       "' with ",format(cli::diff_str("","valid"))," 
                                     names for submodel ", extra1stsub," are:"))
    message(mess)
    diffs <- lapply(names1stsub, function(v) cli::diff_str(mismatch1, v))
    # Sort by number of colored characters...
    diffs <- diffs[order(.unlist(lapply(diffs, function(v)  {
      lcs <- v$lcs
      sum(lcs[lcs$operation != "match", "length"])
    })), decreasing = FALSE)]
    diffs <- lapply(diffs, function(foo) {
      c(.format_invalid(foo),.format_valid(foo))
    })
    diffs <- do.call(rbind, diffs)
    writeLines(paste0(cli::bg_blue(cli::col_black("invalid")),": ", paste(diffs[,1], collapse="   ")))
    writeLines(paste0("  ",cli::bg_green(cli::col_black("valid")),": ", paste(diffs[,2], collapse="   ")))
  } else {
    mess <- cli::format_message(paste0("First invalid name '", mismatch1, 
                                       "' refers to submodel ",extra1stsub," which has no fixed effect."))
    message(mess)
  }
  
  # diffs <- paste(diffs, collapse="    ")
  # #  diffs <- gsub('(.{1,60})(\\s|$)', '\\1\n', diffs) # break lines (but collapses multiple spaces)
  # # format_message()collapse multiple spaces (bad) and provides line breaks (good).
  # mess <- cli::format_message(diffs)
  # message(mess)
  NULL
}

genX2X <- function(matches, names_ori) { 
  if (missing(names_ori)) return(match.call())
  # 'matches' is a list of the form list(<colname X2X>= c(<some colnames LHS>), ...)
  # 'names_ori' are the regressor variables implied by the model formulas; 
  # they are the colnames of the LHS 'X_ori' of the matrix product X_ori %*% X2X.
  matchable_coln_LHS <- .unlist(matches) # LHS regressors for which 'matches' give info  
  matching_coln_X2X <- names(matches) # to become some of the colnames of the matrix product.
  # By convention, all 'matchable_coln_LHS' should appear in the 'names_ori' information:
  # (this convention helps detecting potential errors in definition of 'matches')
  if (length(extracols <- setdiff(matchable_coln_LHS, names_ori))) { 
    .diagnose_string_mismatches(extracols, obs_coln_LHS=names_ori)
    stop(paste0("some submodel regressor name(s) are invalid: '", paste(extracols, collapse="', '"),
               "'. See additional message for hint."))
  }
  # LHS regressors for which 'matches' provide no info will the added to the final result: 
  unmatched_coln_LHS <- setdiff(names_ori, matchable_coln_LHS)
  colnames_X2X <- c(unmatched_coln_LHS, matching_coln_X2X)
  X2X <- matrix(0, ncol=length(colnames_X2X), nrow=length(names_ori),
                dimnames=list(names_ori, colnames_X2X))
  X2X[cbind(unmatched_coln_LHS,unmatched_coln_LHS)] <- 1
  for (st in names(matches)) X2X[matches[[st]],st] <- 1
  X2X
}
# genX2X(list(a1=c("aa","aaa")), c("a","b","c","aa","aaa"))
# genZX2X(list(a1=c("aa","aaa")), c("a","b","c")) # check -> stop

