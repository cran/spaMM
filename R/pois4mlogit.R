# quick base R reformatting
reshape2long  <- function(data, types) {
  template <- rep(0L,length(types))
  names(template) <- types
  vvv <- lapply(seq_len(nrow(data)),
                function(i) {
                  NA_template <- template
                  v  <- data[i,]
                  NA_template[is.na(v)] <- NA_integer_ 
                  vv <- lapply(types, function(typ) {
                    if (is.na(nr <- v[,typ])) {
                      NULL
                    } else {
                      NA_template[typ]  <- 1L
                      v[,types] <- NA_template[types]
                      v[rep(1,nr),]
                    }
                  })
                  do.call(rbind,vv)
                }
  )
  do.call(rbind,vvv)
}

.old_get_surrogate_info <- function(mc, data) {
  mc[["data"]] <- data
  mc[["verbose"]]["getCall"] <- TRUE
  mc[[1L]] <-  get("fitmv", asNamespace("spaMM"), inherits=FALSE)  
  mvp <- eval(mc,parent.frame()) # get processed call
  has_dynoffset <- sapply(lapply(mvp$processed$predictor,.DEPARSE),
                          grepl,pattern="offset(.dynoffset)", fixed=TRUE)
  cum_nobs <- attr(mvp$processed$families,"cum_nobs")
  list(has_dynoffset=has_dynoffset, 
       validPrownames=attr(mvp$processed$data,"validrownames"), # list of vectors of strings
       vec_nobs=diff(cum_nobs) # valid for poisson fit; cf .get_multinom_info() for use.
  ) 
}


# Routinely called by .p4m_by_iters()
.get_surrogate_info <- function(mc, data) {
  mc[["data"]] <- data
  mc[["control"]][["p4m_reactvt_warn"]] <- FALSE
  # mc[["verbose"]]["getCall"] <- TRUE
  mc[[1L]] <-  get("fitmv", asNamespace("spaMM"), inherits=FALSE)  
  mc[["control.HLfit"]] <- .modify_list(mc[["control.HLfit"]], list(getFromTheDepths="surrogate_info"))
  surrogate_info <- .getFromTheDepths(eval(mc,parent.frame()), what="surrogate_info", class="spaMM.FTD.SI")
  mc[["control.HLfit"]][["getFromTheDepths"]] <- NULL
  surrogate_info
}

.get_multinom_info <- function(data, surrogate_info,
                              has_dynoffset=surrogate_info$has_dynoffset,
                              validPrownames=surrogate_info$validPrownames,
                              types
                              ) {
  dim_data <- c(nrow(data), length(has_dynoffset))
  muP_template <- rep(NA_real_, prod(dim_data))
  mnpos_in_template <- rep(FALSE, prod(dim_data))
  dim(muP_template) <- dim(mnpos_in_template) <- dim_data
  rownames(muP_template) <- rownames(data)
  for (it in which(has_dynoffset)) muP_template[validPrownames[[it]],it] <- 0
  invalid_P_info <- is.na(muP_template[,has_dynoffset, drop=FALSE])
  validmncounts <- data[,types]
  validmncounts[invalid_P_info] <- NA_integer_ # 'remove' cases where predictors are missing
  mnsizes <- rowSums(validmncounts, na.rm = TRUE) # only for cases where all predictors are present 
  validmnrows <- 
    rowSums( ! invalid_P_info) > 1L & # predictors present for at least two cases 
    ## second condition required for multinomial fit (but not required for prediction from surrogate Poisson fit):
    mnsizes > 0L # some type(s) were observed for cases where predictors were present.
  mnsizes[ ! validmnrows] <- NA_real_     
  muP_template[ ! validmnrows, surrogate_info$has_dynoffset] <- NA_real_ # mnpos_y must take this into account
  mnpos_in_template[,has_dynoffset] <- TRUE
  mnpos_in_template <- ! is.na(muP_template) & mnpos_in_template
  
  # This code comes before preprocessing of the {fit with all the NA dynoffsets in the right place}
  # and aims to compute and use a vec_nobs identical to that produced by such preprocessing.
  # (distinct from the vec_nobs provided by surrogate$info which does not knew *all* the NA dynoffsets)
  vec_nobs <- colSums(mnpos_in_template) # colsum=0 when ! has dynoffset; hence correction:
  vec_nobs[ ! has_dynoffset]  <- surrogate_info$vec_nobs[ ! has_dynoffset] 
  cum_nobs <- cumsum(c(0,vec_nobs)) 
  mnpos_y <- .unlist(lapply(which(has_dynoffset), 
                            function(v) .subrange(cum_nobs, v)))
  
  validmnsizes <- mnsizes[validmnrows]
  # The logic of the P2M correction is explained in my technical doc.
  P2Mcorr <- sum( validmnsizes*(1-log(validmnsizes))+ lgamma(validmnsizes + 1L) ) 

  list(muP_template=muP_template, 
       mnsizes=mnsizes, # NA for *invalid rows* (incl those with a single type adjustable by a Poisson GLMM)
       MNSIZES=mnsizes, # NA for *invalid rows* (incl those with a single type adjustable by a Poisson GLMM)
       log_mnsizes=log(mnsizes), 
       valid_rows= ( ! is.na(mnsizes)), 
       has_dynoffset=has_dynoffset, 
       mnpos_in_template=mnpos_in_template,
       mnpos_y=mnpos_y,
       # is_fixefM=surrogate_info$is_fixefM,
       P2Mcorr=P2Mcorr
       # which_from_pois_pred= mnpos_in_template[! invalid_P_info] # 
       ##     Might be used to control which values to use from Poisson predict.
       ##     But a simpler approach is to constraint the Poisson fit by controlling NA's in dynoffset.
      )
}

# cf .check_identif_p4m_ranefs() to diagnose the ZAL matrix.
.diagnose_conv_pois4mlogit <- function(mvp, 
                                       processed_call, # for :
                                       processed=processed_call[["processed"]] , 
                                       X2X) {
  has_dynoffset <- processed$multinom_info$has_dynoffset
  if (problem_found <- ( ! all(sapply(processed$families[has_dynoffset],
                                      getElement,name="family")=="poisson"))) {
    warning("Something suspect. ARe submodel families for all multinomial types 'poisson'?",
            immediate. = TRUE)
  } 
  if ( (! problem_found) && ! is.null(X2X)) {
    X.pv <- model.matrix(mvp)
    Xunames <- unique(rownames(X.pv))
    cum_nobs <- attr(X.pv,"cum_nobs")
    col_ranges <- attr(X.pv,"col_ranges")
    n_submodels <- length(col_ranges) 
    for (jt in seq_len(ncol(X.pv))) {
      colj_match_mv_it <- lapply(col_ranges,intersect, y=jt)
      is_colj_in_mv_it <- sapply(colj_match_mv_it, length)>0L
      if (sum(is_colj_in_mv_it)==n_submodels) { # coeff shared among all submodels
        # NA + use rownames + na.omit() to deal with missing info for some submodels. 
        subXwide <- matrix(NA, ncol=n_submodels, nrow=length(Xunames), dimnames = list(Xunames,NULL))
        for (mv_it in seq_len(n_submodels)) {
          resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
          subX_it <- X.pv[resp_range, jt, drop=FALSE]
          subXwide[rownames(subX_it), mv_it] <- subX_it[,1]
        }
        subXwide <- na.omit(subXwide) # no clear alternative 
        subXwide <- round(subXwide,8) # duplicated() def is unclear, but even 1e-16 may matter...
        if (problem_found <- (all(duplicated(subXwide,MARGIN=2L)[-1]))) {
          if (all(subXwide[,1]==1)
              && all(sapply(processed$main_terms_info$fixef_off_terms,attr,which="intercept"))) {
            warning("All submodels seem to have a shared intercept, so model may not be identifiable.",
                    immediate. = TRUE)
            
          } else warning(paste(jt,.stndrdth(jt)," fixed-effect coefficient seems shared among all submodels\n", 
                               "  with identical predictor values, so model may not be identifiable."),
                  immediate. = TRUE)
        }
      }
    }
  }  
  if ( (! problem_found) &&
       (problem_found <- all(sapply(processed$main_terms_info$fixef_off_terms,attr,which="intercept")))) {
    # pb occurs even if intercepts are not shared: at least one model should not have an intercept.
    warning("Maybe unidentifiable model because all submodels have an intercept?",
            immediate. = TRUE)
  } 
  if ( ! problem_found) {
    mess <- paste("No obvious identifiability problem found. Maybe large lambda? lambda=",
                  paste(signif(mvp$lambda,4), collapse=", "))
    warning(mess, immediate. = TRUE)
  }
}

.check_identif_p4m_ranefs <- function(mvp) {
  ZAL <- get_ZALMatrix(mvp)
  if (is.null(ZAL)) return(NULL)
  X.pv <- model.matrix(mvp)
  Xunames <- unique(rownames(X.pv))
  cum_n_u_h <- attr(mvp$lambda,"cum_n_u_h")
  nrand <- length(cum_n_u_h)-1L
  ZALblocks <- vector("list", nrand)
  #
  cum_nobs <- attr(X.pv,"cum_nobs")
  n_submodels <- length(cum_nobs)-1L
  has_dynoffset <- mvp$p4m_info$multinom_info$has_dynoffset
  p4m_submodels <- seq(n_submodels)[has_dynoffset]
  warnings <- NULL
  for (rd in seq_len(nrand)) {
    ident_blocks <- FALSE
    u_h_range <- .subrange(cumul=cum_n_u_h, it=rd)
    resp_range <- .subrange(cumul=cum_nobs, it=p4m_submodels[1L])
    block1 <- ZAL[resp_range, u_h_range, drop=FALSE]
    rownames1 <- rownames(X.pv[resp_range, , drop=FALSE])
    rownames(block1) <- rownames1
    for (p4m_it in 2L:length(p4m_submodels)) {
      resp_range <- .subrange(cumul=cum_nobs, it=p4m_submodels[p4m_it])
      rownamesi <- rownames(X.pv[resp_range, , drop=FALSE])
      blocki <- ZAL[resp_range, u_h_range, drop=FALSE]
      rownames(blocki) <- rownamesi
      shared_rows <- intersect(rownames1, rownamesi)
      ident_blocks <- max(abs(blocki[shared_rows,,drop=FALSE]-block1[shared_rows,,drop=FALSE]))<1e-8
      if ( ! ident_blocks) break
    }
    if (ident_blocks) {
      warnings <- c(warnings,
                    paste0(" All submodels seem to identically share the ",
                           rd,.stndrdth(rd)," random effect."))
    }
  }
  warnings 
} # NULL or a vector of char strings

.rebuid_fixedS <- function(submodels, merged, user_fixed) {
  fixedS <- vector("list",length(submodels))
  # assemble submodels's 'fixed's and global 'fixed' first:
  for (mv_it in seq_along(fixedS)) {
    fixedS[[mv_it]] <- .modify_list(list(), eval(submodels[[mv_it]][["fixed"]], parent.frame())) 
  }
  fixedS <- lapply(fixedS, .reformat_parlist, processed = merged)
  fixedS <- .merge_mv_parlist(fixedS, merged) # now fixedS is a single parlist from the  sub-models specifications
  user_fixed <- .reformat_parlist(user_fixed,processed = merged) # reformat user's global 'fixed' argument
  fixedS <- .modify_list(fixedS, user_fixed) # now fixedS is a single parlist from both sub-model and global specifications
  .preprocess_fixed(fixedS)
}

.calc_p4mprobs <- function(muetablob, multinom_info) {
  has_dynoffset <- multinom_info$has_dynoffset
  if (all(has_dynoffset)) {
    muP_template_mncols <- multinom_info$muP_template
    muP_template_mncols[multinom_info$mnpos_in_template] <- muetablob$mu
  } else {
    muP_template_mncols <- multinom_info$muP_template[,has_dynoffset]
    mumnpos <- .unlist(lapply(muetablob$mv[has_dynoffset],getElement,name="mu"))
    muP_template_mncols[multinom_info$mnpos_in_template[,has_dynoffset]] <- mumnpos
  }
  denoms <- rowSums(muP_template_mncols,na.rm = TRUE)
  .Dvec_times_matrix(1/denoms,muP_template_mncols)
}


.get_new_dynoffset_from_fit <- function(mvp, 
                                        #
                                        dynoffset=mvp$data$.dynoffset,
                                        #
                                        multinom_info=mvp$p4m_info$multinom_info, 
                                        #
                                        etaP_template=multinom_info$muP_template,
                                        has_dynoffset=multinom_info$has_dynoffset,
                                        mnpos_in_template=multinom_info$mnpos_in_template,
                                        #
                                        log_mnsizes=multinom_info$log_mnsizes) {
  eta <- .mvize(mvp$eta,cum_nobs = attr(mvp$families,"cum_nobs"))
  mneta <- .unlist(attr(eta,"mv")[has_dynoffset])
  etaP_template[mnpos_in_template] <- mneta
  # forScrit computed before correction of muP_template by .dynoffset
  forScrit <- rowSums(exp(etaP_template), na.rm=TRUE)
  # From (log) poisson counts to (correctly normalized at convergence) (log) multinomial probabilities:
  etaP_template <- etaP_template - dynoffset 
  # exp(etaP_template) should be (correctly normalized at convergence) multinomial probabilities
  output.dynoffset <- log_mnsizes - matrixStats::rowLogSumExps(etaP_template, na.rm = TRUE)
  attr(output.dynoffset,"forScrit") <- forScrit
  output.dynoffset
}
# I have used a .sanitize_eta_log_link() correction before using rowLogSumExps
# and this was perhaps worse (in unidentif model bc all submodels had an intercept, though)


.init_dynoffset <- function(mvcall, 
                            data, # distinct from those in the call
                            valid_rows,
                            submodels # bc the one in the call is treated as 'symbol'
) {
  mvcall[["submodels"]] <- lapply( submodels, function(subm) {
    if (is.null(form <- subm$formula)) {
      subm[[1]] <- as.formula(sub("offset(.dynoffset)","(1|.id)", 
                                  .DEPARSE(.stripRanefs(subm[[1]])), fixed=TRUE))
    } else subm$formula <- as.formula(sub("offset(.dynoffset)","(1|.id)", 
                                          .DEPARSE(.stripRanefs(subm$formula)), fixed=TRUE))
    subm
  })
  mvcall["covStruct"] <- NULL # I had a conflict with an A matrix...
  data$.id <- seq(nrow(data))
  data$.id[ ! valid_rows] <- NA_real_ # tricky: otherwise pb detected only when fixing ranPars of a model.
  mvcall[["data"]] <- data
  mvcall["multinom_info"] <- NULL # otherwise is_p4m_H-specific code would be triggered, (there's not even P2Mcorr but we don't care).
  #   This would include modifs of design matrices using muetablob$mu values, which seem inapproriate here.
  mvcall$control$p4m <- "o"
  mvcall$control$dyndyn <- FALSE
  mvcall["init.HLfit"] <- NULL # remove v_h in particular which may have wrong dim
  mvcall["init"] <- NULL # remove any user-given value of this arg
  locfit <- eval(mvcall,parent.frame()) 
  ## {All v_h= all rows of the data} may not be in the first submodel 
  ## (even if all lines of the data are informative), 
  ## in which case v_h is not ordered as in the data. So we cannot simply return v_h.
  ## get_ZALMatrix(locfit) %*% locfit$v_h might be used to map to Poisson-valid positions 
  ## in muP_template, but more code needed. Presumably we can use names here:
  v_h <- ranef(locfit, type="uncorrelated")[[1]] # gets correct names(v_h) from the implied design matrix!
  v_h <- v_h[paste(data$.id)] # reorder according to .id \equiv rwos of the data, which 
  #    automatically insert NA's only in rows invalid for the Poisson fit, as opposed to multinom one. 
  #    So we make sure  that there are NA's in all multinom-invalid rows   
  v_h[ ! valid_rows] <- NA_real_ 
  names(v_h) <- NULL # (removes ugly NA names, maybe only cosmetic)
  v_h
}

.get_update_args_from_p4mPQLfit <- function(p4mPQLcall, curr_mvp, 
                                            multinom_info,
                                            log_mnsizes) {
  if (.REMLmess(curr_mvp, return_message=FALSE)) {
    p4mPQLcall$method <- "PQL" 
  } else p4mPQLcall$method <- "PQL/L" 
  p4mPQLcall$progress <- FALSE
  pqlfit <-  eval(p4mPQLcall,parent.frame()) # nested .p4m_by_iters() call !
  pql_dynoffset <- .get_new_dynoffset_from_fit(pqlfit, multinom_info = multinom_info,
                                               log_mnsizes = log_mnsizes)
  init.HLfit <- list(fixef=na.omit(fixef(pqlfit)),
                     v_h=ranef(pqlfit, type="bare.init"))
  list(dynoffset=pql_dynoffset, init.HLfit=init.HLfit)
}

.port_env_overcat <- function(locmess, port_env) {
  IT <- port_env$IT
  ch <- c("|","/","-","\\")[(IT %% 4L)+1L]
  port_env$IT <- IT+1L
  locmess <- paste0(ch,locmess)
  port_env$prevmsglength <- overcat(locmess,prevmsglength=port_env$prevmsglength)
}

# Use missing() systematically to avoid confusion because some of the reset values may be NULL.
#
# processed$port_env is what allows communication between HLfit_body() calls within a fitme call()
#
# (2) inits_by_xLM: not quite sure whether to reinit it if port_env info is missing 
#
# (3) Reinit processed$port_env$objective, otherwise the (irrelevant) final logL 
# of the previous fitme_body with different .dynoffset would control further updating:
#
# (4) init_HLfit is a 'fitmv' specificity to handle preprocessed init.HLfit accounting for mv specificities.
# fitmv uses   mc[["init.HLfit"]] <- merged$"init_HLfit" # we is a single copy in any usage
# Either the processed is not recycled or it is recycled and fitmv_body() is directly called.
# So this is created once and stays there indefinitely unless we reinit it.
#   Then is is read by fitmv_by -> .calc_optim_args_mv -> {init.HLfit <- optim_blob$inits$`init.HLfit`}
# which provides the actual init.HLfit of the HLCor.args within fitmv_body
#   In p4m calls it is useful to reinit it if the dynoffset is modified.
#
# control$port_env (seek also [["control"]]$port_env) 
# is used by pois4lmogit, among other to carry info about the bestfit between .p4m_by_iters 
# in the outer optim and, and to the final .p4m_by_iters refit,
# *through an init.HLfit argument in the .p4m_by_iters call.* This is captured by match.call()
# but not transparent. 
# For processed .p4m_by_iters() call (calling processed fitmv_body calls)
# It seems it is operative only if I use .reinit_processed(init_HLfit=init.HLfit)
.reinit_processed <- function(
    processed, 
    new_offsets, # when .dynoffset and off are updated ; Xb term must be added in HL...bofy
    inits_by_xLM, # reset between fitmv_body calls (2)
    objective, # reset between fitmv_body calls (3)
    init_HLfit # whenever there is some init.HLfit to take into account (4)
    # add 'data' case?
) {
  if ( ! missing(new_offsets)) {
    if ( ! is.null(X_off_Xb_fn <- processed$X_off_Xb_fn)) { 
      processed$off <- 
        X_off_Xb_fn(new_offsets=new_offsets,new_un_betaFix=NULL,new_sc_betaFix=NULL)
    } else { # vanilla fit with use_proc_call -> .reinit_processed() is called,
      # but no etaFix (not confint in particular) so there is no X_off_Xb_fn.
      processed$off <- new_offsets
    }
  }
  if ( ! missing(inits_by_xLM)) 
    processed$envir$inits_by_xLM <- inits_by_xLM
  if ( ! missing(objective))  
    processed$port_env$objective  <- objective
  if ( ! missing(init_HLfit))
    processed$init_HLfit <- init_HLfit
  invisible(NULL)
}

.add_selected_multinom_info <- function(multinom_info, is_p4m_H, is_4_proc_call_4_p4m) {
  if ( ! (is_p4m_H || is_4_proc_call_4_p4m)) multinom_info <-  multinom_info["P2Mcorr"] # fact: don't neet to include $has_dynoffset here.
    multinom_info
}

.cast_as_fitmv_call <- function(mc, submodels) {
  mc[[1L]] <- get("fitmv", asNamespace("spaMM"), inherits=FALSE)  
  mc[ c("to.long","types","n_iter","tol","progress","next_inits","initfn", 
        "fac")] <- NULL # .p4m_by_iters -> fitmv
  mc
}

.cast_as_fitmv_body_call <- function(mc, submodels) {
  # The workflow in fitmv() is .preprocess submodels, .merge_processed(), 
  # build a 'fixedS' argument passed to fitme_body():
  user_global_fixed <- mc$fixed # before it is overwritten
  mc["fixed"] <- NULL # ; and we rebuild a 'fixedS' argument:
  mc[["fixedS"]] <- .rebuid_fixedS(submodels, merged=mc$processed, 
                                   user_fixed = user_global_fixed)
  mc[[1L]] <- get("fitmv_body", asNamespace("spaMM"), inherits=FALSE)  
  mc[ c("to.long","types","n_iter","tol","progress","next_inits","initfn", "fac", # .p4m_by_iters -> fitmv
        "submodels","data","verbose") # fitmv -> fitmv_body
  ] <- NULL
  mc
}

# per se 1st step of pois4mlogit(); also 2nd step, as 'template4objfn' argument
# of .p4m_by_outer_optim() is a .p4m_by_iters() call;
# and formally run by ad-hoc LUarglist extractors.
.p4m_by_iters <- function(submodels, data, to.long=FALSE,
                        init=list(),  # names(init) possibly used in all iterations...
                        control=list(), verbose=c(),
                        ..., # eg, etaFix
                        processed=NULL,
                        next_inits=c("ranPars","v_h","fixef"), 
                        types, n_iter=1000L, tol=1e-5, 
                        max_succ_non_conv=10L,
                        initfn=get_inits_from_fit,
                        init.HLfit=list(), # 
                        multinom_info=NULL,
                        progress=FALSE) {
  mc <- p4mcall <- match.call(expand.dots = TRUE) 
  is_processed_call <- ! is.null(processed)
  if (is_processed_call) { # a 'processed' environment was provided
    mc <- match.call(expand.dots = TRUE) # p4m_by_iters that contains a processed_mvcall...
    mc <- .cast_as_fitmv_body_call(mc, submodels)
    processed <- mc$processed  
    multinom_info <- processed$multinom_info
    valid_rows <- multinom_info$valid_rows
    has_dynoffset <- multinom_info$has_dynoffset
    data  <- processed$data # the *processed* data
    .reinit_processed(processed=processed, 
                      # inits_by_xLM=NULL, # There one case below where we might not want this
                      objective= -Inf) # hopefully offsets has been left up to date by previous computations
  } else {
    if (missing(types)) stop("Argument 'types' is missing.")
    mc <- .cast_as_fitmv_call(mc)
    ## First data$.dynoffset needed before further preprocessing:  _TODO_
    null_init_dynoffset <- is.null(data$.dynoffset)
    if (null_init_dynoffset) data$.dynoffset <- 0 # only for .get_surrogate_info() call

    ## multinom_info, but data may be reshaped in to.long case.
    if (is.null(multinom_info)) { # This case is avoided when this code 
      # is reached through a .get_LUarglist_from_p4m_call() call,
      # BUT it does routinely occur otherwise.
      surrogate_info <- .get_surrogate_info(mc, data) 
      has_dynoffset <- surrogate_info$has_dynoffset
      if (length(types) != sum(has_dynoffset)) 
        stop("Length of 'types' does not match number of submodels with a .dynoffset.")
      if (to.long) {
        if  ( ! all(has_dynoffset)) 
          stop("'to.long=TRUE' feasible only when all submodels are components of multinomial model.")
        data <- reshape2long(data=data, types = types)
        surrogate_info <- .get_surrogate_info(mc, data) 
      } 
      multinom_info <- .get_multinom_info(data, surrogate_info, types=types) 
    } else has_dynoffset <- multinom_info$has_dynoffset
    
    ## mc[["multinom_info"]] 
    is_p4m_H <- control[["p4m"]]=="H"     
    is_4_proc_call_4_p4m <- .safe_true(control["get_proc_call_4_p4m"][[1L]])
    mc[["multinom_info"]] <- 
      .add_selected_multinom_info(multinom_info, is_p4m_H, is_4_proc_call_4_p4m)

    ## data$.dynoffset, mc[["data"]]
    valid_rows <- multinom_info$valid_rows
    if (null_init_dynoffset) { # : typically FALSE when is_p4m_byoo
      data$.dynoffset <- .init_dynoffset(mvcall=mc, data= data, valid_rows=valid_rows, 
                                         submodels=submodels) 
    }
    data$.dynoffset[ ! valid_rows] <- NA_real_ 
    #     (incl. those with a single type adjustable by a Poisson GLMM), this makes sure 
    #     that predict(<surrogate Poisson  fit>) generates only values matching 'mnpos_in_template'.
    mc[["data"]] <- data
    
    
    if (is_4_proc_call_4_p4m) {
      mc <- .get_HLCorcall_4_p4m(mc=mc) # result is call to inner-estimating fn with $processed in the call
      control["get_proc_call_4_p4m"] <- FALSE
      mc$control <- control
      return(mc)  
    } 
    
    if (.safe_true(verbose["getCall"][[1L]])) { # may not occur in tests
      mvp <- eval(mc,parent.frame()) 
      return(mvp) # getCall() a *fitmv* call, not necess desirable. ####
      # But note that getCall(<pois4mlogit>) returns a pois4mlogit() call. (in checks)
    }
  }
  
  muP_template <- multinom_info$muP_template
  mnpos_in_template <- multinom_info$mnpos_in_template
  mnsizes <- multinom_info[["MNSIZES"]] # presumably 1L for long data 
  log_mnsizes <- multinom_info$log_mnsizes
  
  # control$port_env is provided by .p4m_by_outer_optim() -> .get_template4objfn()
  
  if ( is_p4m_byoo <- (! is.null(p4m_port_env <- control$port_env))) {
    # pois4mlogit -> outer optim ('H' step) -> (NOT first) .p4m_by_inits here -> its own first fitmv_body 
    if ( length(init.HLfit) || # directly passed by .numInfo_objfn(). Maybe it shouldn't as it increases complexity ?
         ! is.null(init.HLfit <- p4m_port_env$"init.HLfit")) {
      if (is_processed_call) {
        .reinit_processed(processed=processed, 
                          inits_by_xLM = NULL,
                          init_HLfit= init.HLfit)
      } else  mc[["init.HLfit"]] <- init.HLfit   
    } # ELSE pois4mlogit -> outer optim ('H' step) -> first .p4m_by_inits here -> this first fitmv_body 
    # (after possibly many of 'o' step)
    # If we want an init.HLfit here we can pass it either through inits_by_xLM 
    # or through an explicit .p4m_by_inits(init.HLfit) (not the case AFAICS)
    # ** So these is one case where we might want to keep inits_by_xLM **
  } else if (is_meta_p4m <- ( ! is.null(meta_p4m_port_env <- control$meta_port_env))) {
    # confint -> (NOT first) .p4m_by_inits('H') here -> its own first fitmv_body 
    if ( ! is.null(init.HLfit <- meta_p4m_port_env$"init.HLfit")) { # speculative code
      # this is not the first .p4m_by_inits ('H')
      if (is_processed_call) {
        .reinit_processed(processed=processed, 
                          inits_by_xLM = NULL,
                          init_HLfit= init.HLfit)
      } else  mc[["init.HLfit"]] <- init.HLfit   
    } 
  } # ELSE {this is the first .p4m_by_inits ('H') in confint case} OR
    # {first fitmv_body 'o' of a model fit} (everything is presumably NULL in both cases)

  
  
  curr_mvp <- eval(mc,parent.frame()) # 1st curr_mvp (fitmv or fitmv_body call) ####
  # potential .sendFromTheDepths(LUarglist = LUarglist) from fitmv_body in this eval()

  PQL_already_run <- (curr_mvp$HL[1L]==0L)
  
  output.dynoffset <- .get_new_dynoffset_from_fit(curr_mvp, multinom_info = multinom_info,
                                                  log_mnsizes = log_mnsizes)
  names.init <- names(init) # possibly used in all iterations...
  prevmsglength <- 0L
  oldScrit <- oldOcrit <- Inf
  Scrit <- Ocrit <- NA
  oldlogL <- -Inf
  logL <- logLik(curr_mvp)
  w_opt_off <- 1
  notwarned <- progress>=0L
  control <- mc[["control"]]
  control.HLfit <- mc[["control.HLfit"]]
  control.HLfit$"algebra" <- curr_mvp$how$algebra
  latest.successful.input.dynoffset  <- data$".dynoffset"
  next.dynoffset <- w_opt_off*output.dynoffset +(1-w_opt_off)*latest.successful.input.dynoffset# \sum^J
  problem <- any(is.infinite(next.dynoffset)) || 
    anyNA((next.dynoffset[valid_rows]))
  if (problem) next.dynoffset <- .init_dynoffset(mvcall=mc, data= data, valid_rows=valid_rows, 
                                                 submodels=submodels) 
  data$".dynoffset" <- next.dynoffset

  successiveNotConv <- ! is.null(curr_mvp$warnings$innerNotConv)
  succInnerNotConv <- FALSE
  min_update_time <- Inf
  
  dyndyn <- control$dyndyn
  
  for (it in 1L+seq_len(n_iter-1L)) { # BEGIN MAIN LOOP ####
    if (progress > 2L) str(data$".dynoffset")
    newinits <- initfn(curr_mvp, to_fn="fitmv_body")[["init"]]
    # To remove values from the newinits, one should have  (! "ranPars" %in% next_inits) and explicit NA's 
    if ( ! "ranPars" %in% next_inits) { # non-default
      if ("init"  %in% next_inits) { # persistent use of the original user init in all iters: presumably very inefficient
        newinits <- .modify_list(newinits, init) 
      } else newinits <- newinits[names.init] # also presumably inefficient, but not tested recently.
    }  
    time1 <- Sys.time()
    if ( ! problem) {
      if ("v_h" %in% next_inits) init.HLfit$v_h <- ranef(curr_mvp, type="bare.init")
      if ("fixef" %in% next_inits) init.HLfit$fixef <- fixef(curr_mvp)
    }
    if (is_processed_call) { # direct update on fitmv_body call
      processed$data  <- data # with new .dynoffset
      # Initialize next fitmv_body call: 
      mc["init"] <- list(newinits) 
      .reinit_processed(processed=processed, 
                        new_offsets=model.offset.HLfit(curr_mvp, data=data),
                        inits_by_xLM=NULL, 
                        objective= -Inf,
                        init_HLfit= init.HLfit
      )
      # The current control$port_env (with scaled values) will be used to initialize the next iteration,
      # provided we remove any original init.HLfit ...
      curr_mvp <- eval(mc,parent.frame()) #  eval(fitmv_body,...)
    } else {
      # $processed is recreated in each iter so we use init.HLfit for the first "inner fit" of the call
      curr_mvp <- update(curr_mvp, data=data, init.HLfit=init.HLfit, # loop update (*fitmv* call)
                         init=newinits, control.HLfit=control.HLfit, control=control) 
    }
    update_time <- .timerraw(time1)
    min_update_time <- min(min_update_time, update_time) 
    # If convergence problem, get init dynoffset and init.HLfit from nested p4m_by_iters PQL fit:  
    # update_time>'10'... is only a guess but seems effective. 
    if (is_p4m_byoo && 
        ! PQL_already_run && 
        (successiveNotConv || update_time>10*min_update_time))  {
      info_from_pql <- .get_update_args_from_p4mPQLfit(
        p4mPQLcall=p4mcall, curr_mvp, multinom_info = multinom_info, log_mnsizes = log_mnsizes)
      PQL_already_run <- TRUE
      data$.dynoffset <- info_from_pql$dynoffset
      curr_mvp <- update(curr_mvp, data=data, # extra fitmv call to handle convergence problem
                         init.HLfit=info_from_pql$init.HLfit, 
                         init=newinits, control.HLfit=control.HLfit, control=control) 
    } 
    
    output.dynoffset <- .get_new_dynoffset_from_fit(
      curr_mvp,
      dynoffset=curr_mvp$data$.dynoffset, # has effectively used data$.dynoffset rather than mvp$... for a long time.
      # makes no difference as long as update() does not change data$.dynoffset
      multinom_info=multinom_info, 
      etaP_template=muP_template,
      has_dynoffset=has_dynoffset,
      mnpos_in_template=mnpos_in_template,
      log_mnsizes=log_mnsizes)
    # Ocrit use new offsets computed from (old-offset)-included muP's
    Ocrit <- mean(abs(output.dynoffset-data$.dynoffset), na.rm=TRUE)
    Scrit <- 2*(attr(output.dynoffset, "forScrit")/mnsizes-1) # '2*' adjusted for good results in tests
    Scrit <- mean(abs(Scrit), na.rm = TRUE)
    logL <- logLik(curr_mvp)
    dlogL <- logL-oldlogL # always Inf in first iteration
    
    mess <- paste0(it,": logL:",signif(logL,5), # note oldlogL can be restored before this is printed 
                   " Ocrit: ",signif(Ocrit,3L), # which is why this part of the message is built now.
                   " Scrit: ",signif(Scrit,3L),"         ")
    
    cond <- Ocrit<tol && Scrit<tol
    if (dyndyn) cond <- cond && abs(dlogL)<tol 
    if (cond) {
      case <- "+|"
      break
    }
    
    control.HLfit <- .modify_list(control.HLfit,
                                  list(spaMM_tol=list(logL_tol=5e-5*max(1,Scrit/tol))))
    control$bobyqa$rhoend <- substitute(
      max(1e-8, control$bobyqa$rhobeg*1e-8*max(1,val)),   list(val=Scrit/tol)
    )
    control$nloptr$xtol_rel <- substitute(
      4e-6*max(1,val),   list(val=Scrit/tol)
    )
    dScrit <- Scrit-oldScrit
    dOcrit <- Ocrit-oldOcrit
    if (notwarned && abs(dlogL) < tol) {
      # if a>0 in (Scrit=a+b lambda_b^t) Scrit_lam_b is O(lambda_b^t), t Scrit_lam_b should still vanish
      # Ideally, a=0, Scrit_lam_b is O(lambda_b), t Scrit_lam_b will diverge 
      Scrit_lam_b  <- - dScrit/oldScrit # ideally large. So we test if Scrit_lam_b is small, as this is suspect.
      Ocrit_lam_b  <- - dOcrit/oldOcrit # same.
      # Rkably one of the identifiable 'pollen' fits shows these two lam_b's 
      #   staying large and ~constant for some time. Only the logL improves.
      
      if ( abs(Ocrit_lam_b)<tol && abs(Scrit_lam_b)<tol) { 
        message("Something suspect. Diagnosing...")
        .diagnose_conv_pois4mlogit(mvp=curr_mvp, processed_call=.get_HLCorcall_4_p4m(mc=mc, data=data),
                                   X2X=eval(mc$X2X))
        notwarned <- FALSE
      }
    }
    
    var_dyn <- NA
    
    if (dlogL > 0) { # logL improves  
      # The dOcrit condition has an effect in the tests.
      if (dScrit>0 && dOcrit>0) { # "bad": try to exit this region of slow progress quickly
        fac <- (40*tol*Scrit/(tol/4 + abs(dScrit))+
        #      Penalize large relative increases in Scrit more than in "++<" case   
        #                                vv        
                  Scrit*Scrit/(Scrit/20 + 4*abs(dScrit)))/
          (40*tol+Scrit)
        if (successiveNotConv) {
          w_opt_off <- 1+ log(max(1, fac))
        } else w_opt_off <- max(1, fac) 
        case <- "+-"
      } else { # presumably good cases but some care needed
        var_dyn <- var(output.dynoffset*sign(data$.dynoffset)/(1e-6+abs(data$.dynoffset)), na.rm=TRUE)
        if (is.na(var_dyn)) break # notably, major convergence pb, muP is Inf, dynoffset is-Inf...
        if (var_dyn < 0.05) {
          ## when Scrit approaches tol, the first term dominates => higher factor
          fac <- (40*tol*Scrit/(tol/4 + abs(dScrit))+
                    Scrit*Scrit/(Scrit/20 + abs(dScrit)))/
            (40*tol+Scrit)
          w_opt_off <- 1+ log(max(1, fac)) # rather than max(1, fac): helps mmfit6var  

          case <- "++<"
        } else { # good case
          w_opt_off <- 1 
          case <- "++>"
        }
      }
    } else { # Three different cases adjusted on mmfit6var 
      if (dScrit<0 && dOcrit<0) {
        # logL decreases but dynoffset converges according to both crits, which is 
        # quite possibly good step since dynoffset does not maximizes logL.
        # => 'Accept" dynoffset (w_opt_off>1) 
        fac <- Scrit/(tol/20 + abs(dScrit))
        if (is_p4m_byoo) {
          w_opt_off <- max(1.5, log(fac)) # 1 + log(max(1, fac))
        } else w_opt_off <- max(1.5, fac)
        case <- "--"
      }  else if (dScrit>0 && dOcrit>0) { # {crits diverge = bad}, and logL decreases, 
        # following large 'w_opt_off'?, or rather following "--".
        w_opt_off <- (1+exp(-dOcrit))/2 # semble effectif sur mmfit6var
        logL <- oldlogL
        case <- "-+"
      } else  {
        # ambiguous case: dynoffset converges according to one crit but not other.
        w_opt_off <- 1
        case <- "-+"
      }
    }
    
    oldlogL <- logL
    if (case != "-+") latest.successful.input.dynoffset <- data$".dynoffset" 

    if (progress>1L) {
      if (progress>2L) mess <- paste(mess, case, signif(w_opt_off,3))
      prevmsglength <- overcat(mess, prevmsglength)
    }
    
    oldScrit <- Scrit
    oldOcrit <- Ocrit
    
    innerNotConv <- ! is.null(curr_mvp$warnings$innerNotConv)
    successiveNotConv <- innerNotConv*(successiveNotConv+1L)
    if (succInnerNotConv <- (successiveNotConv> max_succ_non_conv)) {
      break
    }
    
    next.dynoffset <- w_opt_off*output.dynoffset +(1-w_opt_off)*latest.successful.input.dynoffset# \sum^J
    problem <- any(is.infinite(next.dynoffset)) || 
      anyNA((next.dynoffset[valid_rows]))
    if (problem) {
      break 
    }
    data$".dynoffset" <- next.dynoffset
  } # END MAIN LOOP #######################################
  
  if (progress>1L) {
    if (progress>2L) mess <- paste(mess, case, signif(w_opt_off,3))
    prevmsglength <- overcat(mess, prevmsglength)
  }
  
  if (n_iter > 1L && progress>=0L) {
    if (! cond) {
      cat("\n")
      warnmess <- paste(".dynoffset did not converge in",n_iter,
                        "iterations (Ocrit: ",signif(Ocrit,3L),
                        ", Scrit: ",signif(Scrit,3L),")")
      curr_mvp$warnings$p4m_by_iters_not_conv <- warnmess
    } else if (progress>1L) { # case with overcat's in the inner loop above
      cat("\n")
    } else if (progress==1L) { 
      locmess <- paste(".dynoffset converged in", it,"iterations.   ")
      if (is_p4m_byoo) { # overcat over distinct .p4m_by_iter() calls
        .port_env_overcat(locmess, port_env=p4m_port_env)
      } else if (is_meta_p4m) { # overcat over distinct calls within .confint_LRT_single_par_p4m()
        .port_env_overcat(locmess, port_env=meta_p4m_port_env)
      } else cat(paste(locmess,"\n"))
    } else if (progress>0L) {
      cat(".")
    }
  }
  curr_mvp$warnings$succInnerNotConv <- succInnerNotConv
  if (is_p4m_byoo) {
    if (succInnerNotConv) p4m_port_env$any_succInnerNotConv <- succInnerNotConv
    if (logLik(curr_mvp) > p4m_port_env$logL+0.0001) {
      p4m_port_env$logL <- logLik(curr_mvp)    
      p4m_port_env$"init.HLfit" <- list(v_h=ranef(curr_mvp, type="bare.init"),
                                        fixef=fixef(curr_mvp))
      # print(head(p4m_port_env$"init.HLfit"$"v_h"), digits=6)
    }
  }
  if ( ! inherits(curr_mvp,"HLfitlist") && ! is.call(curr_mvp) ) {
    if (is_processed_call) p4mcall[["processed"]] <- NULL
    curr_mvp$call <- p4mcall
    curr_mvp$p4m_info <- list(Ocrit=Ocrit, Scrit=Scrit, it=it, mnsizes=mnsizes,
                              multinom_info=multinom_info) 
    curr_mvp$how$fnname <- "pois4mlogit"
  }
  class(curr_mvp) <- c("pois4mlogit",class(curr_mvp))
  curr_mvp 
}

# 2nd step of pois4mlogit()
.p4m_by_outer_optim <- function(fitobject, 
                                template4objfn,
                                objfn.extras,
                             skeleton, LowUp,
                             # from .numinfo:
                             transf, 
                             which=NULL,
                             verbose=FALSE,
                             refit_hacks=list(),
                             check_deriv,
                             ...) {
  objective <- .get_objective(fitobject)
  proc_info <- list(objective=objective) 
  
  if (TRUE) { # no processed call: nonstandard use of the hlcorcall argument
    # Here preprocessing will be called, 
    #   in each iteration within the call of .p4m_by_iters() in .numInfo_objfn(), 
    # 
    ## its time to define another objfn (____F I X M E____)

    # .get_HLCorcall_4_p4m(template4objfn) cannot be used without additional programming

    optr <- .safe_opt(init=unlist(skeleton), objfn = .numInfo_objfn, objfn.extras=objfn.extras,
                      LowUp = LowUp,
                      lower=unlist(LowUp$lower),upper=unlist(LowUp$upper),
                      verbose=FALSE,
                      # additional arguments for .numInfo_objfn:
                      skeleton=skeleton, 
                      hlcorcall=template4objfn, # .p4m_by_iters() call; cannot be a processed HLCorcall (HLCor.obj|HLfit.obj) as p4m iterations are needed
                      transf=transf, # signals transformed input (init, skeleton, LowUp). 
                      full_beta=fixef(fitobject), 
                      objective=proc_info$objective, 
                      moreargs=.get_moreargs(fitobject), ...)
    
  } else eval(.lot_of_crap_from_earlier_attempts)

  list(optr=optr,
       transf=transf, # signals transformed input. FALSE when called from numInfo(), 
       full_beta=fixef(fitobject), 
       objective=proc_info$objective, 
       moreargs=.get_moreargs(fitobject))
}

# p4m_reactvt_warn and ppc_reactvt_warn can both reactivate warnings from 
#  .warn_glm_poisson_rates_0_once_per_fit() when set to TRUE,
#  by setting the latter fn's 'warned_glm_poisson_rates_0' to FALSE.
# ppc_reactvt_warn is used in .preprocess. By default it is TRUE 
#  when set by .reformat.CONTROL() in .preprocess().
#  But this is overridden by being set to FALSE by pois4mlogit() by default. 
#  This inhibits multiple warnings for multiple fitmv calls within p4m fits.
# p4m_reactvt_warn is used at beginning of pois4mlogit(). 
#  By default it is set to TRUE here, and thus p4m reactivates warnings.
#  But confint for p4m fits sets it to FALSE, so the multiple pois4mlogit calls()
#  won't produce multiple warnings (although the first can, provided
#  confint is called while 'warned_glm_poisson_rates_0' is FALSE).
.reformat_p4m_controls <- function(control, has_bar, p4m=NULL, 
                                   p4m_reactvt_warn=TRUE, ppc_reactvt_warn=FALSE) {
  if ( ! length(control)) {
    control  <- list(wdfac=2, p4m=p4m,grad=FALSE, has_bar=has_bar, 
                     p4m_reactvt_warn=TRUE, ppc_reactvt_warn=FALSE)
  } else {
    if (is.null(control[["wdfac"]])) control[["wdfac"]] <- 2 
    if (is.null(control[["grad"]])) control[["grad"]] <- FALSE 
    if (is.null(control[["p4m_reactvt_warn"]])) control[["p4m_reactvt_warn"]] <- TRUE 
    if (is.null(control[["ppc_reactvt_warn"]])) control[["ppc_reactvt_warn"]] <- FALSE 
    if ( ! is.null(p4m)) control[["p4m"]] <- p4m 
    
  }
  if (is.null(control[["use_proc_call"]])) control[["use_proc_call"]] <- FALSE
  if (is.null(control[["p4m"]])) control[["p4m"]] <- ""
  control[["dyndyn"]] <- (control[["p4m"]]=="W")
  if (has_bar) {
    if (control[["p4m"]]=="o") warning('p4m="o" control is not recommended for mixed-effect models.',
                                       immediate. = TRUE, call. = FALSE)
    if (control[["p4m"]] %in% c("","W")) control[["p4m"]] <- "oH"
  } else {
    if (control[["p4m"]]=="oH") warning('p4m="oH" control is superfluous for fixed-effect models.',
                                       immediate. = TRUE, call. = FALSE)
    if (control[["p4m"]]=="") control[["p4m"]] <- "o"
    if (control[["p4m"]]=="W") control[["p4m"]] <- "H"
  }
  control
}

.grepBarsInSubmodel <- function(subm) {
  if (is.null(form <- subm$formula)) form <- subm[[1]] 
  length(grep(pattern = "|", .DEPARSE(form),fixed=TRUE,value=FALSE))
}

.check_gradient_p4m_outer_optim <- function(mc, optim_blob, solution, skeleton) {
  # gradient check. Note that this can return nonzero gradient if solution is at boundary
  cat("Gradient check:\n")
  side <- .calc_grad_side_arg(solution)
  uside <- unlist(side)
  gr_APHL <- - grad(func = .numInfo_objfn, x =optim_blob$optr$solution, side=uside, skeleton=skeleton, 
                    hlcorcall=mc, # cannot be a processed HLCorcall here
                    transf=optim_blob$transf, 
                    full_beta=optim_blob$full_beta, 
                    # no ... in grad call(): otherwise ... that go in the match.call() go there too.
                    # before, the ... ere here, but then I had to add a 'verbose' argument to pois4mlogit()
                    # so that this 'verbose' was not passed to grad().
                    objective=optim_blob$objective, moreargs=optim_blob$moreargs) 
  names(gr_APHL) <- c(names(unlist(skeleton[setdiff(names(skeleton), "etaFix")])), 
                      names(skeleton$etaFix$beta))
  print(gr_APHL)
}

# Build template call 'template4objfn' used by .p4m_by_outer_optim(). 
# The objective function in .p4m_by_outer_optim() itself calls *.p4m_by_iters()*, 
# allowing iterations of the dynoffset in the objective fn;
# so template4objfn cannot be an HLCorcall (which does not iterate). 
# 'template4objfn' has to be modified from the initial .p4m_by_iters() call, including
# .p4m_by_iters() local changes to it, such as the spaMM_tol$logL_tol),
# through to the control[["port_env"]] environment
.get_template4objfn <- function(mc, # .p4m_by_iters() call, with data from mvp$data, hence final .dynoffset of 1st step
                                mvp, # fit from 1st step
                                p4mcontrol, skeleton) {
  template4objfn <- mc # .p4m_by_iters()
  # mc itself will be used to get some of the outer optim controls, then replaced later.
  control <- template4objfn[["control"]]
  control[["port_env"]] <- list2env(list(logL=logLik(mvp),
                                         any_succInnerNotConv=FALSE,
                                         prevmsglength=0L, IT=0L), parent = emptyenv())
  template4objfn[["control"]] <- control
  template4objfn["initfn"] <- NULL 
  # THis template4objfn typically already has an init.HLfit.
  # The next line controls the init of the first .p4m_by_iters() within .p4m_by_outer_optim()
  template4objfn[["init.HLfit"]] <- list(fixef=na.omit(fixef(mvp)), v_h=mvp$v_h)
  if (length(template4objfn[["init"]])) { # tested by numInfo(BbyP) as numInfo() adds init values,
    # potentially creating a conflict between init and fixed value.
    template4objfn[["init"]] <- eval(template4objfn[["init"]])
    cskeleton <- .canonizeRanPars(skeleton, corr_info = mvp$ranef_info$sub_corr_info, 
                                  rC_transf=.spaMM.data$options$rC_transf)
    template4objfn[["init"]] <- .remove_from_cP(template4objfn[["init"]],u_names=names(unlist(cskeleton)))
  }
  template4objfn # .p4m_by_iters call
  # TRY template4objfn$processed <- # not the right return element
  #        .get_HLCorcall_4_p4m(template4objfn$processed
}
# .get_proc_p4m_bi_call() vs .get_proc_call_4_p4m():
# .get_proc_p4m_bi_call() returns the input .p4m_by_iters() call after
# adding $processed, while .get_proc_call_4_p4m()
# returns an HL(Cor)_body call.
.get_proc_p4m_bi_call <- function(mc) {
  mc$control$"get_proc_call_4_p4m" <- TRUE
  # eval(mc, parent.frame()) would be HLCor_body call,
  # here we keep only its $processed
  mc$processed <- eval(mc, parent.frame())$processed
  mc$control$"get_proc_call_4_p4m" <- FALSE
  mc
}

.adjust_multinom_info <- function(p4mbicall) {
  control <- p4mbicall[["control"]]
  processed <- p4mbicall[["processed"]]
  if (control$use_proc_call) {
    if (control[["p4m"]] =="o") {
      processed$multinom_info[["mnsizes"]] <- NULL
    } else if (is.null(processed$multinom_info[["mnsizes"]])) {
      processed$multinom_info[["mnsizes"]] <- processed$multinom_info$MNSIZES
      if (is.null(processed$multinom_info[["mnsizes"]])) stop("something suspect.") # Remove later ?
    }
  }
}


# 'control' is passed as fitmv's 'control', with eg $p4m used deep in fitting algos.
# So it's not just locally used and some elements have to be in this 'control' arg.
pois4mlogit <- function(submodels, data, to.long=FALSE,
                        init=list(), # it goes in all local calls except the final refit
                        control=list(wdfac=2, p4m="",grad=FALSE), 
                        ..., next_inits=c("ranPars","v_h","fixef"), 
                        types, n_iter=1000L, tol=c(1e-3,1e-5), initfn=get_inits_from_fit,
                        progress=FALSE) {
  time1 <- Sys.time()
  mc <- oricall <- match.call(expand.dots = TRUE)
  has_bar <- any(sapply(submodels, .grepBarsInSubmodel))
  mc[["control"]] <- control <- .reformat_p4m_controls(control, has_bar=has_bar) # adds ppc_reactvt_warn=FALSE
  p4mcontrol <- control[["p4m"]]
  ## whether is_p4m_H is TRUE in .solve_IRLS_as_ZX() in step _s_ or not depends on the _s_th character of control[["p4m"]]
  
  ##### 1ST STEP, using .p4m_by_iters()
  #
  if (nchar(p4mcontrol)>1L) { # then the tol of the first, iterative, step is made less stringent:
    if (missing(tol) || 
      length(tol)>1L) {
      mc[["tol"]] <- tol[[1]] 
    } else  { # explicit of length 1 but two steps
      mc[["tol"]] <- tol*100 
    } 
  } else if (missing(tol)) {
    if (p4mcontrol=="o") {
      mc[["tol"]] <- 1e-6 
    } else {
      mc[["tol"]] <- 1e-5 # p4mcontrol=="H". 
      #  Modified below if 'has_bar' but not 'has_ranPars2fit' (e.g., confint...)
    }
  } else mc[["tol"]] <- tol[[1]]
  mc[["n_iter"]] <- n_iter[1] 
  
  mc[["control"]][["p4m"]] <- substr(p4mcontrol,1,1) # 'o' or 'H' (befon .adjust_multinom_info)
  .warn_glm_poisson_rates_0_once_per_fit(reinit=control$p4m_reactvt_warn) # reactvt TRUE by default
  
  # If ranPars are fixed and a .dynoffset is available, 
  #    calling *pois4mlogit* with "H" directly may make sense. 

  # Otherwise, calling *.p4m_by_iters* with p4m='H' is not enough 
  # bc ranPars are not correctly estimated by maximizing logL for given .dynoffset. 
  # Therefore, the correct fit has to outer-optimize ranPars
  # using an internal objective that iterates the dynoffset.
  # However, the result of a such an 'H' step appears rather good. 
  
  mc[[1L]] <- get(".p4m_by_iters", asNamespace("spaMM"), inherits=FALSE)  
  
  if (control[["use_proc_call"]]) {
    mc <- .get_proc_p4m_bi_call(mc) # get processed .p4m_by_iters() call ####
    .adjust_multinom_info(p4mbicall = mc)
  }
  
  mvp <- eval(mc, parent.frame()) # "o" by .p4m_by_iters() or .p4m_by_iters(., processed) ####
  
  if (.safe_true(control["get_proc_call_4_p4m"][[1L]]))
    return(mvp) # return the processed inner fitmv call from .p4m_by_iters()
  # which can in particular be obtained by 
  # update(<p4m fit>, control=list("get_proc_call_4_p4m"=TRUE))
  
  if ( ! is.null(warnmess <- mvp$warnings$p4m_by_iters_not_conv)) 
    warning(warnmess, immediate. = TRUE)
  chk_ranefs_blob <- .check_identif_p4m_ranefs(mvp)
  if (length(chk_ranefs_blob)) sapply(chk_ranefs_blob,
                                      warning, immediate. = TRUE,call. = FALSE)
  dfs <- mvp$dfs 
  external_fix_in_out_info <-   mvp$ranef_info$internal_fix_in_out_info
  CorrEst_and_RanFix_type <- attr(mvp$CorrEst_and_RanFix,"type")
  lambda_type <- mvp$lambda.object$type
  
  # end 1st step
  
  ##### 2ND STEP
  if (nchar(p4mcontrol)>1L) {
    mc[["control"]][["p4m"]] <- substr(p4mcontrol,2,2) # "H" by default.
    if (control[["use_proc_call"]]) {
      .adjust_multinom_info(p4mbicall=mc) 
    }  
    control.HLfit <- mc[["control.HLfit"]]
    control.HLfit$"algebra" <- mvp$how$algebra
    mc[["control.HLfit"]] <- control.HLfit
    #
    if (missing(tol) || length(tol)>1L) { 
      mc[["tol"]] <- tol[[2]] 
    } else  {
      mc[["tol"]] <- tol
    } 
    mc[["init.HLfit"]] <- list(fixef=na.omit(fixef(mvp)),
                               v_h=ranef(mvp, type="bare.init"))
    has_ranPars2fit <- sum(.unlist(mvp$dfs[c("p_lambda","p_corrPars")]))
    if (has_ranPars2fit) { # has_ranPars2fit ####
      # there are ranPars to estimate
      #### full outer-optim step:
      ## get initial values:
      fullinit4LUarglist <- get_inits_from_fit(mvp,inner_lambdas = TRUE)
      fullinit4LUarglist$init.HLfit$fixef <- NULL # otherwise .get_LUarglist_from_p4m_call() 
      #    will include beta in its various returned elements.
      # # if (DEVELp4m) eval(.init_call_outer_from_binomial_fit)
      mc[["n_iter"]] <- n_iter[min(2L, length(n_iter))]
      mc[["tol"]] <- tol[min(2L, length(tol))] # restores value from initial call.
      mc[["control"]] <- .modify_list(mc[["control"]], list(refit=FALSE)) # overrides any 'refit'
      # v_h correspond to ranef(mvp,type = "uncorrelated"), to be used in product ZAL . v_h
      mc[["data"]] <- mvp$data # for the .dynoffset in particular
      
      mc[["multinom_info"]] <- mvp$p4m_info$multinom_info # 
      # : .get_LUarglist_from_p4m_call() expects this at this place 
      # (ugly duplicate info if use_proc_call, but not immediate tidy)
      LUarglist  <- .get_LUarglist_from_p4m_call(mc=mc, fullinit=fullinit4LUarglist) 
      
      use_transf <- TRUE
      LowUp <- optimBounds(mvp, LUarglist=LUarglist, 
                           transf=use_transf) # previous comment: "! user's lower/upper ranCoefs are handled by .objfn_locoptim()"
      # but here it's rather "! user's lower/upper ranCoefs are handled by .numInfo_objfn()",
      # which receives here the appropriate info by the 'objfn.extras' arg passed through
      # .p4m_by_outer_optim() -> .safe_opt(., objfn = .numInfo_objfn, objfn.extras=objfn.extras, ....)
      
      # 'skeleton' will provide the initial value of .safe_optim
      if (use_transf) {
        skeleton <- LUarglist$init.optim 
      } else {
        skeleton <- LUarglist$canon.init 
        stop("need to convert LowUp to canonical space !")
      }
      template4objfn <- .get_template4objfn(mc=mc, mvp=mvp, p4mcontrol=p4mcontrol, skeleton=skeleton)
      ### 'H' by .p4m_by_outer_optim() itself calls *.p4m_by_iters()*, ####
      optim_blob <- .p4m_by_outer_optim(mvp, # this provides optimization controls
                                        template4objfn=template4objfn, 
                                        objfn.extras=LUarglist,
                                        skeleton=skeleton, LowUp=LowUp, 
                                        check_deriv=FALSE,
                                        transf=use_transf)
      
      optr <- optim_blob$optr
      solution <- relist(optr$solution, skeleton)
      
      if (! is.null(trRanCoefs <- solution$trRanCoefs)) {
        # (1) partially fixed values
        constraints <- LUarglist$ranFix$ranCoefs 
        if (length(constraints)) {
          solution$trRanCoefs <- .partially_fix_trRancoefs(trRanCoefs, 
                                                           constraints=LUarglist$ranFix$ranCoefs)
        }
        # (2) box constraints
        if ((length(LUarglist[["user.lower"]]$ranCoefs) || 
             length(LUarglist[["user.upper"]]$ranCoefs)) ) solution <- 
            .apply_transformed_box_constr(solution, skeleton=NULL, 
                                          user.lower=LUarglist[["user.lower"]], 
                                          user.upper=LUarglist[["user.upper"]], transf=use_transf)
      }
      
      if (control$grad) .check_gradient_p4m_outer_optim(template4objfn, optim_blob, solution, skeleton) 
      
      ## AFTER THE OUTER OPTIM ####
      # with fixed ranPars to get the full fit object
      # Don't use trLambda directly: .preprocess() does not handle it correctly. See details in .canonizeRanPars() 
      outerP <- .canonizeRanPars(solution,corr_info=mvp$ranef_info$sub_corr_info, 
                                 checkComplete=FALSE, rC_transf=.spaMM.data$options$rC_transf)
      fixed_or_outer <- .modify_list(oricall$fixed, outerP) 
      template4objfn[["fixed"]] <- fixed_or_outer
      bestfit <- template4objfn$control$port_env$bestfit
      template4objfn[["data"]][[".dynoffset"]] <- 
        .get_new_dynoffset_from_fit(bestfit,multinom_info = bestfit$p4m_info$multinom_info,
                                    log_mnsizes = log(bestfit$p4m_info[["mnsizes"]]))
      if ( ! is.null(solution$etaFix)) template4objfn[["etaFix"]] <- 
        .modify_list(solution$etaFix, # typically NULL ("H" refitting using init.HLfit instead)    
                     eval(oricall$etaFix))
      template4objfn["init"] <- NULL 
      template4objfn[["init.HLfit"]] <- list(fixef=na.omit(fixef(bestfit)),
                                             v_h=ranef(bestfit, type="bare.init"))
      
      mvp <- eval(template4objfn,parent.frame()) ## final .p4m_by_iters 'H' refit ####  

      if (template4objfn$control$port_env$prevmsglength) cat("\n")
      if ( ! is.null(warnmess <- mvp$warnings$p4m_by_iters_not_conv)) 
        warning(warnmess, immediate. = TRUE)
      ## ad-hoc fixes, restoring features of original call without fixed ranPars.
      mvp$dfs <- dfs # counting the dfs of the estimated parameters not of the final refit...
      mvp$ranef_info$external_fix_in_out_info <- external_fix_in_out_info
      attr(mvp$CorrEst_and_RanFix,"type") <- CorrEst_and_RanFix_type
      mvp$lambda.object$type <- lambda_type
      mvp$call <- oricall
      # Final result should have same ZAL, v_h as an equivalent binomial fit and the .dynoffset implied by denom of NLPredictor.
      attr(mvp,"optimInfo") <- list(LUarglist=LUarglist, init.optim=skeleton,
                                    optim.pars=solution,
                                    objective=optim_blob$objective,
                                    rC_transf=.spaMM.data$options$rC_transf)
      if (mvp$warnings$succInnerNotConv) { # $warnings$ has kept info about final refit;
        if (template4objfn$control$port_env$any_succInnerNotConv) { # port_env has kept info 
          # about successive _by_iters (="Intermediate step(s)") fits of the H step;
          # ('Inner' refers to steps within a by_iters fit)
          mvp$warning$any_succInnerNotConv <- "Intermediate step(s) and last step of pois4mlogit fit interrupted."
        } else mvp$warnings$succInnerNotConv <- "Last step of pois4mlogit fit interrupted."
      } else {
        mvp$warnings$succInnerNotConv <- NULL
        if (template4objfn$control$port_env$any_succInnerNotConv) {
          mvp$warnings$any_succInnerNotConv <- "Intermediate step(s) of pois4mlogit fit interrupted."
        } 
      }
    } else if (has_bar) { # has_bar, here implying that ranPars are all fixed ####
      # if has_bar here, then by default p4mcontrol is "oH" and this block is run
      # (which is OK bc "o" and "H" fits differ when ranPars are all fixed).
      # Only if user asked for "o", this "H" step is not run (with a note)
      if (missing(tol)) mc[["tol"]] <- 1e-6
      # ideally the p4m should have been simply "H" in this case: has_bar but not has_ranPars2fit.
      # we would need to detect this case before the first step while currently has_ranPars2fit is evaluated after.
      mvp <- eval(mc, parent.frame()) #### final .p4m_by_iters # ____F I X M E____ not consistently spprec through all steps? (trace .solve_IRLS_as_ZX)
      if (mvp$warnings$succInnerNotConv) { # unique .p4m_by_iters of second (H) step has been interrupted.
        mvp$warnings$succInnerNotConv <- "Last step of pois4mlogit fit interrupted." 
      } else mvp$warnings$succInnerNotConv <- NULL
    }
  } else { # no second step
    if (mvp$warnings$succInnerNotConv) { # unique .p4m_by_iters of first step has been interrupted
      mvp$warnings$succInnerNotConv <- "Unique step of pois4mlogit fit interrupted." 
    } else mvp$warnings$succInnerNotConv <- NULL
  }
  
  if ( is.list(sXaug <- mvp$envir$sXaug) && # -> spprec
       ! is.null(dcdb_p4m <- sXaug$AUGI0_ZX$dcdb_p4m)) { 
    mvp$envir$sXaug$AUGI0_ZX$dcdb_p4m <- .unscale(dcdb_p4m, scale=attr(mvp$X.pv,"scale_info")) 
    # $dcdb_p4m used post-fit in spprec to compute beta SE
  }
  if (! is.call(mvp) ) {
    fit_time <- .timerraw(time1) 
    mvp$warnings$"chk_ranefs_blob" <- chk_ranefs_blob
    mvp$how$fit_time <- structure(fit_time,
                                  message="Please use how(<fit object>)[['fit_time']] to extract this information cleanly.")
  }
  
  mvp
}

predict.pois4mlogit <- function(object, newdata=NULL, verbose=NULL, 
                                na.action=na.omit, ...) {
  .predict.pois4mlogit(object=object, newdata=newdata, verbose=verbose, 
                       na.action=na.action, ...) # this hides the private 'check' argument
}

.predict.pois4mlogit.check <- function(object, mnsizes, ...) {
  ndata <- object$data
  ndata$".dynoffset" <- ndata$".dynoffset" - log(mnsizes) # quick correction without predict.HLfit call.
  # This code was first conceived to avoid the double predict call in .predict.pois4mlogit(), 
  # but may be confusing because the predictions then 'exactly' sum to 1 only if convergence was 'exact'.
  predict.HLfit(object, newdata=ndata, ...)
}

#tests: test-simulate.R but not enough bc no case with missing predictors in valid rows
# test-simulate-composite-antisym.R fills this gap.
.predict.pois4mlogit <- function(object, newdata=NULL, verbose, binding=FALSE, 
                                 na.action, ..., check=FALSE) {
  multinom_info <- object$p4m_info$multinom_info
  if ( is.null(has_dynoffset <- multinom_info$has_dynoffset)) 
    stop("Old object, info not found where expected.")
  if (check) {
    # This code may be useful to provide some measure of convergence inaccuracy on prediction.
    if ( ! is.null(newdata)) stop("check=TRUE does not handle 'newdata'." ) 
    if ( is.null(mnsizes <- object[["mnsizes"]])) stop("Old object, $mnsizes not found where expected.")
    .predict.pois4mlogit.check(object, mnsizes=mnsizes, verbose=verbose, ...)
  } else { # default case: double predict call, the first without the '...'
    # This double call makes sure that frequencies sum to 1, even if convergence was not 'exact'. 
    nc <- length(formula(object))
    if (is_null_input_data <- is.null(newdata)) {
      newdata <- object$data # bc we want to control the .dynoffset before the 1st predict.
      good_pred_resps <- ! is.na(multinom_info$muP_template) # ! (invalid rows or missing resp)
      valid_rows <- object$p4m_info$multinom_info$valid_rows
    } else {
      if (FALSE) {
        # This block implied that all rows become 'valid' 
        # if there are newdata without .dynoffset (OK), BUT more problematically 
        # that otherwise the .dynoffset stored in the object's $data was used: problem
        # discussed on an example in test-pois4mlogit-paternity-interaction.R
        if (is.null(newdata$".dynoffset")) newdata$".dynoffset" <- 0
        valid_rows <- ! is.na(newdata$".dynoffset")
      } else valid_rows <- rep(TRUE, nrow(newdata))
    }
    newdata$".dynoffset"[valid_rows] <- 0
    # newdata is never NULL at this point!
    pred <- predict.HLfit(object, newdata=newdata, na.action=na.exclude, verbose=verbose) 
                          # without the dots! type for this first call is response, not "link"... 
    # 'na.action=na.exclude' keeps the NA in predictions, in line with the concept in base R.
    # but this does not work with extra arguments below.
    if ( is.character(binding) ) return(pred) # *RETURN* experimental
    good_predictors <- ! is.na(pred)
    
    predmat <- matrix(pred, ncol=nc) # hence this may have rather arbitrary NA's even in ! valid_rows
    if ( ! is_null_input_data) {
      # then the 'predmat' may contain some values for lines with info for a single type 
      # what should be invalid rows: we correct 'valid_rows'.
      valid_rows <- rowSums( ! is.na(predmat[,has_dynoffset, drop=FALSE]))>1L
      newdata$".dynoffset"[ ! valid_rows] <- NA
      predmat[ ! valid_rows,] <- NA
      good_next_predicts <- ! is.na(predmat)
    } else {
      good_next_predicts <- ! is.na(pred) # accouts for good_predictors but not for missing response.
      predmat[ ! good_pred_resps] <- NA # for missing response... important for next rowsums
    }
    valid_rowsums <- rowSums(predmat[valid_rows,has_dynoffset, drop=FALSE], na.rm =TRUE)
    newdata$".dynoffset"[valid_rows] <- - log(valid_rowsums)
    # newdata is never NULL at this point!
    pred <- predict.HLfit(object, newdata=newdata, verbose=verbose, binding=binding, ...) 
    # When called by simulate(), binding is typically NA internally so pred is here (a vector?) NOT a 1-col matrix
    # When called for pdep_effects (single submodel) pred is here a 1-col matrix.
    mostAttrs <- attributes(pred)
    predmat[good_next_predicts] <- pred
    if ( is_null_input_data) predmat[ ! good_pred_resps] <- NA # for missing response... but the rowSums may no longer be 1...
    # the new rowSums(predmat[valid_rows,has_dynoffset, drop=FALSE], na.rm =TRUE) 
    #   should be 1...1 *IF* type is response
    rownames(predmat) <- rownames(newdata)
    mv <- apply(predmat,2L, na.omit, simplify=FALSE)
    pred <- unlist(apply(predmat,2L, identity, simplify=FALSE)) # aim is to keep names... AND NAs
    if ( ! identical(na.action,na.exclude)) pred <- na.action(pred)
    if ("dim" %in% names(mostAttrs)) pred <- 
      matrix(pred, dimnames = list(names(pred), NULL))
    attr(pred,"mv") <- mv 
    attr(pred,"nobs") <- sapply(mv,length) # predict does not return NA's: this does not count NA's in simulate()
    # attr present but never used before? only here it can be expected to match the mv attr 
    # simulate() returns another nobs value ! 
    mostAttrs <- mostAttrs[setdiff(names(mostAttrs), 
                                   names(attributes(pred)))]
    for (st in setdiff(names(mostAttrs),c("dim","dimnames"))) attr(pred,st) <- mostAttrs[[st]] 
    pred
  }
}
# returned pred must be full length including NA's if any.

