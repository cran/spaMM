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
  
  list(muP_template=muP_template, 
       mnsizes=mnsizes, # NA for *invalid rows* (incl those with a single type adjustable by a Poisson GLMM)
       valid_rows= ( ! is.na(mnsizes)), 
       has_dynoffset=has_dynoffset, 
       mnpos_in_template=mnpos_in_template,
       mnpos_y=mnpos_y
       # which_from_pois_pred= mnpos_in_template[! invalid_P_info] # 
       ##     Might be used to control which values to use from Poisson predict.
       ##     But a simpler approach is to constraint the Poisson fit by controlling NA's in dynoffset.
      )
}

.diagnose_conv_pois4mlogit <- function(mvp, 
                                       processed_call, # for :
                                       processed=processed_call[["processed"]] , 
                                       X2X) {
 
  if ( ! all(sapply(processed$families,getElement,name="family")=="poisson")) {
    warning("Something suspect. Are all submodel families 'poisson'?",
            immediate. = TRUE)
  # } else if ( ! all(sapply(lapply(processed$predictor,.DEPARSE),grepl,pattern="offset(.dynoffset)", 
  #                            fixed=TRUE))) {
  #   warning("Something suspect. Maybe 'offset(.dynoffset)' term missing from some submodels?",
  #           immediate. = TRUE)
  } else if ( ! is.null(X2X)) {
    X.pv <- model.matrix(mvp)
    cum_nobs <- attr(X.pv,"cum_nobs")
    col_ranges <- attr(X.pv,"col_ranges")
    n_submodels <- length(col_ranges) 
    for (jt in seq_len(ncol(X.pv))) {
      colj_match_mv_it <- lapply(col_ranges,intersect, y=jt)
      is_colj_in_mv_it <- sapply(colj_match_mv_it, length)>0L
      if (sum(is_colj_in_mv_it)==n_submodels) { # coeff shared among all submodels
        subXcols <- vector("list",length=n_submodels) 
        for (mv_it in seq_along(subXcols)) {
          resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
          subXcols[[mv_it]] <- X.pv[resp_range, jt, drop=FALSE]
        }
        subXwide <- do.call(cbind, subXcols)
        if (all(duplicated(as.data.frame(t(subXwide)))[-1])) {
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
  } else if (all(sapply(processed$main_terms_info$fixef_off_terms,attr,which="intercept"))) {
    # pb occurs even if intercepts are not shared: at least one model should not have an intercept.
    warning("Maybe unidentifiable model because all submodels have an intercept?",
            immediate. = TRUE)
  } else {
    warning("Something suspect. Maybe unidentifiable model?",
            immediate. = TRUE)
  }
}

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

.get_new_dynoffset_from_fit <- function(mvp, 
                                        #
                                        dynoffset=mvp$data$.dynoffset,
                                        #
                                        multinom_info, 
                                        #
                                        muP_template=multinom_info$muP_template,
                                        has_dynoffset=multinom_info$has_dynoffset,
                                        mnpos_in_template=multinom_info$mnpos_in_template,
                                        #
                                        log_mnsizes) {
  eta <- .mvize(mvp$eta,cum_nobs = attr(mvp$families,"cum_nobs"))
  mneta <- .unlist(attr(eta,"mv")[has_dynoffset])
  muP_template[mnpos_in_template] <- mneta
  # forScrit computed before correction of muP_template by .dynoffset
  forScrit <- rowSums(exp(muP_template), na.rm=TRUE)
  # From (log) poisson counts to (correctly normalized at convergence) (log) multinomial probabilities:
  muP_template <- muP_template - dynoffset 
  muP_template[mnpos_in_template] <- .sanitize_eta_log_link(muP_template[mnpos_in_template], 
                                                            max=40, y=mvp$y[multinom_info$mnpos_y])
  muP <- exp(muP_template) #  (correctly normalized at convergence) multinomial probabilities
  output.dynoffset <- log_mnsizes-log(rowSums(muP,na.rm = TRUE))
  attr(output.dynoffset,"forScrit") <- forScrit
  output.dynoffset
}


.init_dynoffset <- function(mc, data, submodels, valid_rows) {
  
  mc[["submodels"]] <- lapply( submodels, function(subm) {
    if (is.null(form <- subm$formula)) {
      subm[[1]] <- as.formula(sub("offset(.dynoffset)","(1|.id)", 
                                  .DEPARSE(.stripRanefs(subm[[1]])), fixed=TRUE))
    } else subm$formula <- as.formula(sub("offset(.dynoffset)","(1|.id)", 
                                          .DEPARSE(.stripRanefs(subm$formula)), fixed=TRUE))
    subm
  })
  data$.id <- seq(nrow(data))
  data$.id[ ! valid_rows] <- NA_real_ # tricky: otherwise pb detected only when fixing ranPars of a model.
  mc[["data"]] <- data
  mc["multinom_info"] <- NULL # otherwise 'H' p4m-specific (is_p4m_H) code would be triggered by presence of multinom_info$mnsizes.
  #   This would include modifs of design matrices using muetablob$mu values, which seem inapproriate here.
  mc$control$p4m <- "o"
  mc["init.HLfit"] <- NULL # remove v_h in particular which may have wrong dim
  locfit <- eval(mc,parent.frame()) # first fit# using input 'init' if any
  # .get_new_dynoffset_from_fit(locfit, dynoffset=0, 
  #                             multinom_info=multinom_info, log_mnsizes=log_mnsizes, ...)
  # .get_new_dynoffset_from_fit(locfit, dynoffset=0, 
  #                             multinom_info=multinom_info, log_mnsizes=log_mnsizes, ...)
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

# per se 1st step of pois4mlogit(); also 2nd step, as 'template4objfn' argument
# of .p4m_by_outer_optim() is a .p4m_by_iters() call;
# and formally run by ad-hoc LUarglist extractors.
.p4m_by_iters <- function(submodels, data, to.long=FALSE,
                        init=list(),  # names(init) possibly used in all iterations...
                        control=list(), verbose=c(),
                        ...,
                        processed_call=NULL,
                        next_inits=c("ranPars","v_h","fixef"), 
                        types, n_iter=1000L, tol=1e-5, 
                        max_succ_non_conv=10L,
                        initfn=get_inits_from_fit,
                        curr_mvp=NULL,
                        multinom_info=NULL,
                        progress=FALSE) {
  if (is.null(processed_call)) {
    if (missing(types)) stop("Argument 'types' is missing.")
    update_fitmv_body <- control$update_fitmv_body
    if (is.null(update_fitmv_body)) update_fitmv_body <- 
        .safe_true(control["get_processed_p4m_by_iters_call"][[1L]])
    # fac_is_1  <- abs(fac-1) < 100*.Machine$double.eps
    mc <- p4mcall <- match.call(expand.dots = TRUE) 
    mc["to.long"] <- NULL
    mc["types"] <- NULL
    mc["n_iter"] <- NULL
    mc["tol"] <- NULL
    mc["progress"] <- NULL
    mc["next_inits"] <- NULL
    mc["initfn"] <- NULL
    mc["fac"] <- NULL
    # if (missing(types)) { # formula is not equiv to submodel 
    #   types <- sapply(submodels, function(form) deparse(form[[1]][[2]]))
    # }
    # Need a first .dynoffset before preprocessing:  
    init_dynoffset_run <- FALSE
    null_init_dynoffset <- is.null(data$.dynoffset)
    if (null_init_dynoffset) data$.dynoffset <- 0 # only for .get_surrogate_info() call
    # Build template matrix with NAs for missing response or predictors per submodel.
    
    if (is.null(multinom_info)) { # This case is avoided when this code 
      # is reached through a .get_LUarglist_from_p4m_call() call,
      # BUT it does routinely occur otherwise.
      # .get_LUarglist_from_p4m_call() -> .get_surrogate_info() might work now (not recently checked).
      surrogate_info <- .get_surrogate_info(mc, data) 
      has_dynoffset <- surrogate_info$has_dynoffset
      if (length(types) != sum(has_dynoffset)) 
        stop("Length of 'types' does not match number of submodels with a .dynoffset.")
      if (to.long) {
        if  ( ! all(has_dynoffset)) 
          stop("'to.long=TRUE' feasible only when all submodels are components of multinomial model.")
        # if length(setdiff(types, colnames(data))) stop("responses are not variables in the data.frame: provide a 'types' argument")
        data <- reshape2long(data=data, types = types)
        surrogate_info <- .get_surrogate_info(mc, data) 
      } 
      multinom_info <- .get_multinom_info(data, surrogate_info, types=types) 
    } else has_dynoffset <- multinom_info$has_dynoffset
    if (control[["p4m"]]=="H")   mc[["multinom_info"]] <- multinom_info
    
    valid_rows <- multinom_info$valid_rows
    muP_template <- multinom_info$muP_template
    mnpos_in_template <- multinom_info$mnpos_in_template
    mnsizes <- multinom_info$mnsizes # presumably 1L for long data
    log_mnsizes <- log(mnsizes)
    mc[[1L]] <- get("fitmv", asNamespace("spaMM"), inherits=FALSE)  
    if (null_init_dynoffset) {
      init_dynoffset_run <- TRUE
      data$.dynoffset <- .init_dynoffset(mc, data, submodels, valid_rows) 
    }
    data$.dynoffset[ ! multinom_info$valid_rows] <- NA_real_ 
    #     (incl. those with a single type adjustable by a Poisson GLMM), this makes sure 
    #     that predict(<Poisson surrogate fit>) generates only values matching 'mnpos_in_template'.
    
    mc[["data"]] <- data
    if (update_fitmv_body) {
      #
      user_global_fixed <- mc$fixed # before it is overwritten
      mc <- .get_processed_call(mc=mc, data=data) 
      if (.safe_true(control["get_processed_p4m_by_iters_call"][[1L]])) return(mc)
      processed <- mc[["processed"]]
      data  <- processed$data # the *processed* data
      #
      # The workflow in fitmv() is .preprocess submodels, .merge_processed(), build a 'fixedS'
      # argument passed to fitme_body(), but the latter is replaced by a 'fixed' argument
      # in {the call, to the inner-estimating fn, returned by .get_processed_call()}.
      # This 'fixed' has additional vakues for the outer-optimized parameters, so we cannot use it:
      mc["fixed"] <- NULL # ; and we rebuild a 'fixedS' argument:
      mc[["fixedS"]] <- .rebuid_fixedS(submodels, merged=processed, user_fixed = user_global_fixed)
      mc[[1L]] <- get("fitmv_body", asNamespace("spaMM"), inherits=FALSE)  
    } 
  } else { # a processed call was provided
    user_global_fixed <- mc$fixed # before it is overwritten
    mc <- processed_call
    processed <- mc[["processed"]]
    data  <- processed$data # the *processed* data
    #
    # The workflow in fitmv() is .preprocess submodels, .merge_processed(), build a 'fixedS'
    # argument passed to fitme_body(), but the latter is replaced by a 'fixed' argument
    # in {the call, to the inner-estimating fn, returned by .get_processed_call()}.
    # This 'fixed' has additional vakues for the outer-optimized parameters, so we cannot use it:
    mc["fixed"] <- NULL # ; and we rebuild a 'fixedS' argument:
    mc[["fixedS"]] <- .rebuid_fixedS(submodels, merged=processed, user_fixed = user_global_fixed)
    mc[[1L]] <- get("fitmv_body", asNamespace("spaMM"), inherits=FALSE)  
    update_fitmv_body <- TRUE
    multinom_info <- processed$multinom_info
    muP_template <- multinom_info$muP_template
    mnpos_in_template <- multinom_info$mnpos_in_template
  }
  
  if (.safe_true(verbose["getCall"][[1L]])) { 
    mvp <- eval(mc,parent.frame()) # first fit# using input 'init' if any
    return(mvp) # a *fitmv* call, not necess desirable. 
    # But note that getCall(<pois4mlogit>) returns a pois4mlogit() call. (in checks)
  } else if (is.null(curr_mvp)) curr_mvp <- eval(mc,parent.frame()) # potential .sendFromTheDepths(LUarglist = LUarglist) in this eval()
  
  output.dynoffset <- .get_new_dynoffset_from_fit(curr_mvp, multinom_info = multinom_info,
                                                  log_mnsizes = log_mnsizes)
  names.init <- names(init) # possibly used in all iterations...
  prevmsglength <- 0L
  oldScrit <- oldOcrit <- Inf
  Scrit <- Ocrit <- NA
  oldlogL <- -Inf
  logL <- logLik(curr_mvp)
  d_off_fac <- 1
  notwarned <- progress>=0L
  control <- mc[["control"]]
  control.HLfit <- mc[["control.HLfit"]]
  successful.input.dynoffset  <- data$".dynoffset"
  next.dynoffset <- d_off_fac*output.dynoffset +(1-d_off_fac)*successful.input.dynoffset# \sum^J
  problem <- any(is.infinite(next.dynoffset)) || 
    anyNA((next.dynoffset[valid_rows]))
  if (problem) next.dynoffset <- .init_dynoffset(mc, data, submodels, valid_rows) 
  data$".dynoffset" <- next.dynoffset

  successiveNotConv <- ! is.null(curr_mvp$warnings$innerNotConv)

  for (it in 1L+seq_len(n_iter-1L)) {
    if (progress > 2L) str(data$".dynoffset")
    newinits <- initfn(curr_mvp, to_fn="fitmv_body")[["init"]]
    # To remove values from the newinits, one should have  (! "ranPars" %in% next_inits) and explicit NA's 
    if ( ! "ranPars" %in% next_inits) { # non-default
      if ("init"  %in% next_inits) { # persistent use of the original user init in all iters: presumably very inefficient
        newinits <- .modify_list(newinits, init) 
      } else newinits <- newinits[names.init] # also presumably inefficient, but not tested recently.
    }  
    if (update_fitmv_body) { # direct update on fitmv_body call
      processed$data  <- data
      # The processed offset is stored in processed$off, which has to be recomputed
      # (and its mv-length is distinct from that of data$".dynoffset" )
      processed$off  <- model.offset.HLfit(curr_mvp, data=data)
      mc[["processed"]]  <- processed
      # Initialize next outer optim: 
      mc["init"] <- list(newinits) 
      # The current port_env (with scaled values) will be used to initialize the next "inner fit",
      # provided we remove any original init.HLfit ...
      mc["init.HLfit"] <- list(NULL) 
      # ... but we also have to hack the $port_env$objective, otherwise the (irrelevant) final logL 
      # of the previous fit with different .dynoffset would control further updating:
      processed$port_env$objective  <- -Inf
      curr_mvp <- eval(mc,parent.frame())
    } else {
      init.HLfit <- list()
      # $processed is recreated in each iter so we use init.HLfit for the first "inner fit" of the call
      if ( ! problem) {
        if ("v_h" %in% next_inits) init.HLfit$v_h <- ranef(curr_mvp, type="bare.init")
        if ("fixef" %in% next_inits) init.HLfit$fixef <- fixef(curr_mvp)
      }
      # if (DEVELp4m) eval(.all_inits_from_DEBUGfit) # contains DEVELp4m_verif <<- TRUE
      # if (FALSE) { # to examine very slow iterations
      #   update_args <- list(object=curr_mvp, data=data, init.HLfit=init.HLfit, init=newinits,
      #                       control.HLfit=control.HLfit, control=control)
      #   save(update_args, file=paste0("update_args.",it,".rda"))
      # }
      curr_mvp <- update(curr_mvp, data=data, init.HLfit=init.HLfit, init=newinits,
                    control.HLfit=control.HLfit, control=control) 
    }
       
    #
    output.dynoffset <- .get_new_dynoffset_from_fit(
      curr_mvp,
      dynoffset=curr_mvp$data$.dynoffset, # has effectively used data$.dynoffset rather than mvp$... for a long time.
      # makes no difference as long as update() does not change data$.dynoffset
      multinom_info=multinom_info, 
      muP_template=muP_template,
      has_dynoffset=has_dynoffset,
      mnpos_in_template=mnpos_in_template,
      log_mnsizes=log(mnsizes))
    # Ocrit use new offsets computed from (old-offset)-included muP's
    Ocrit <- mean(abs(output.dynoffset-data$.dynoffset), na.rm=TRUE)
    Scrit <- attr(output.dynoffset, "forScrit")-mnsizes
    Scrit <- mean(abs(Scrit), na.rm = TRUE)
    logL <- logLik(curr_mvp)
    dlogL <- logL-oldlogL # always Inf in first iteration
    
    mess <- paste0(it,": logL:",signif(logL,5), # note oldlogL can be restored before this is printed 
                   " Ocrit: ",signif(Ocrit,3L), # which is why this part of the message is built now.
                   " Scrit: ",signif(Scrit,3L),"         ")
    
    cond <- Ocrit<tol && Scrit<tol
    if (cond) break
    
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
        .diagnose_conv_pois4mlogit(mvp=curr_mvp, processed_call=.get_processed_call(mc=mc, data=data),
                                   X2X=eval(mc$X2X))
        notwarned <- FALSE
      }
    }
    
    var_dyn <- NA
    
    if (dlogL > 0) { # logL improves   # but dynoffset does not maximizes logL 
      # The dOcrit condition has an effect in the tests.
      if (dScrit>0 && dOcrit>0) { # "bad": try to exit this region of slow progress quickly
        fac <- (40*tol*Scrit/(tol/4 + abs(dScrit))+
        #      Penalize large relative increases in Scrit more than in "++<" case   
        #                                vv        
                  Scrit*Scrit/(Scrit/20 + 4*abs(dScrit)))/
          (40*tol+Scrit)
        d_off_fac <- max(1, fac) 
        case <- "+-"
      } else {
        var_dyn <- var(output.dynoffset*sign(data$.dynoffset)/(1e-6+abs(data$.dynoffset)), na.rm=TRUE)
        if (is.na(var_dyn)) break # notably, major convergence pb, muP is Inf, dynoffset is-Inf...
        if (var_dyn < 0.05) {
          ## when Scrit approaches tol, the first term dominates => higher factor
          fac <- (40*tol*Scrit/(tol/4 + abs(dScrit))+
                    Scrit*Scrit/(Scrit/20 + abs(dScrit)))/
            (40*tol+Scrit)
          d_off_fac <- max(1, fac) 
          case <- "++<"
        } else {
          d_off_fac <- 1 
          case <- "++>"
        }
      }
    } else {
      if (dScrit>0 && dOcrit>0) { # logL decreases and crits diverge...
        # typically occurs for large 'd_off_fac'
        d_off_fac <- 1 
        logL <- oldlogL
        case <- "-+"
      } else { #  # logL decreases but dynoffset converges according to at least one crit: *accept* step
        # OK since dynoffset does not maximizes logL.
        d_off_fac <- max(1, Scrit/(tol/20 + abs(dScrit)))
        case <- "--"
      }
    }
    
    oldlogL <- logL
    if (case != "-+") successful.input.dynoffset <- data$".dynoffset" 

    if (progress>1L) {
      if (progress>2L) mess <- paste(mess, case, signif(d_off_fac,3))
      prevmsglength <- overcat(mess, prevmsglength)
    }
    
    oldScrit <- Scrit
    oldOcrit <- Ocrit
    
    innerNotConv <- ! is.null(curr_mvp$warnings$innerNotConv)
    successiveNotConv <- innerNotConv*(successiveNotConv+1L)
    if (successiveNotConv> max_succ_non_conv) break # _____F I X M E_____ create message
    
    next.dynoffset <- d_off_fac*output.dynoffset +(1-d_off_fac)*successful.input.dynoffset# \sum^J
    problem <- any(is.infinite(next.dynoffset)) || 
      anyNA((next.dynoffset[valid_rows]))
    if (problem) {
      break 
      # next.dynoffset <- .init_dynoffset(mc, data, submodels, valid_rows) 
    }
    
    data$".dynoffset" <- next.dynoffset
  } # end main loop
  
  if (n_iter > 1L && progress>=0L) {
    if (! cond) {
      cat("\n")
      warnmess <- paste("pois4mlogit() fit did not converge in",n_iter,
                        "iterations (Ocrit: ",signif(Ocrit,3L),
                        ", Scrit: ",signif(Scrit,3L),")")
      curr_mvp$warnings$p4m_by_iters_not_conv <- warnmess
    } else if (progress>1L) { # case with overcat's
      cat("\n")
    } else if (progress>0L) {
      print(paste("Fit converged in", it,"iterations."), quote=FALSE)
    }
  }
  if ( ! inherits(curr_mvp,"HLfitlist") && ! is.call(curr_mvp) ) {
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

    # .get_processed_call(template4objfn) cannot be used without additional programming

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

.reformat_p4m_controls <- function(control, has_bar, p4m=NULL) {
  if (is.null(control)) {
    control  <- list(wdfac=2, p4m=p4m,grad=FALSE, has_bar=has_bar)
  } else {
    if (is.null(control[["wdfac"]])) control[["wdfac"]] <- 2 
    if (is.null(control[["grad"]])) control[["grad"]] <- FALSE 
    if ( ! is.null(p4m)) control[["p4m"]] <- p4m 
    
  }
  if (is.null(control[["p4m"]])) control[["p4m"]] <- ""
  if (has_bar) {
    if (control[["p4m"]]=="o") warning('p4m="o" control is not recommended for mixed-effect models.',
                                       immediate. = TRUE, call. = FALSE)
    if (control[["p4m"]]=="") control[["p4m"]] <- "oH"
  } else {
    if (control[["p4m"]]=="oH") warning('p4m="oH" control is ignored for fixed-effect models.',
                                       immediate. = TRUE, call. = FALSE)
    control[["p4m"]] <- "o"
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
  mc[["control"]] <- control <- .reformat_p4m_controls(control, has_bar=has_bar)
  p4mcontrol <- control[["p4m"]]
  ## whether is_p4m_H is TRUE in .solve_IRLS_as_ZX() in step _s_ or not depends on the _s_th character of control[["p4m"]]
  
  ##### 1ST STEP, using .p4m_by_iters()
  #
  if (nchar(p4mcontrol)>1L) { # then the tol of the first, iterative, step is made less stringent:
    if (missing(tol) || # _____F I X M E_____ more adaptive control of the loop? 
      length(tol)>1L) {
      mc[["tol"]] <- tol[[1]] 
    } else  { # explicit of length 1 but two steps
      mc[["tol"]] <- tol*100 
    } 
  } else if (missing(tol)) {
    if (p4mcontrol=="o") {
      mc[["tol"]] <- 1e-6
    } else mc[["tol"]] <- 1e-5 # p4mcontrol=="H"
  } else mc[["tol"]] <- tol[[1]]
  mc[["n_iter"]] <- n_iter[1] 
  
  mc[["control"]][["p4m"]] <- substr(p4mcontrol,1,1) # 'o' or 'H'
  # : a single 'H' step is not enough bc ranPars are not correctly estimated 
  # by maximizing logL for given .dynoffset. 
  # Therefore, the correct fit has to outer-optimize ranPars
  # using an internal objective that iterates the dynoffset.
  # Hower, the result of a 1st 'H' step are rather good. 
  
  # if (DEVELp4m) eval(.init_call_iter_from_binomial_fit)
  mc[[1L]] <- get(".p4m_by_iters", asNamespace("spaMM"), inherits=FALSE)  
  mvp <- eval(mc, parent.frame()) 
  if ( ! is.null(warnmess <- mvp$warnings$p4m_by_iters_not_conv)) 
    warning(warnmess, immediate. = TRUE)
  
  dfs <- mvp$dfs 
  # end 1st step
  
  ##### 2ND STEP
  fullinit4LUarglist <- get_inits_from_fit(mvp,inner_lambdas = TRUE)
  has_ranPars2fit <- sum(.unlist(mvp$dfs[c("p_lambda","p_corrPars")]))
  
  if (has_ranPars2fit) { # there are ranPars to estimate
    if (nchar(p4mcontrol)>1L) {
      #### full outer-optim step:
      ## get initial values:
      fullinit4LUarglist$init.HLfit$fixef <- NULL # otherwise .get_LUarglist_from_p4m_call() 
      #    will include beta in its various returned elements.
      # # if (DEVELp4m) eval(.init_call_outer_from_binomial_fit)
      mc[["n_iter"]] <- n_iter[min(2L, length(n_iter))] # restores value from initial call.
      mc[["tol"]] <- tol[min(2L, length(tol))] # restores value from initial call.
      mc[["control"]] <- .modify_list(mc[["control"]], list(refit=FALSE)) # overrides any 'refit'
      # v_h correspond to ranef(mvp,type = "uncorrelated"), to be used in product ZAL . v_h
      mc[["data"]] <- mvp$data # for the .dynoffset in particular
      mc[["multinom_info"]] <- mvp$p4m_info$multinom_info # very useful for .get_LUarglist_from_p4m_call()
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
      
      # if (DEVELp4m) {
      #   if (substr(control$p4m,1,1)=="o") skeleton$etaFix$beta <- beta_eta # was TRUE in first version checking the v_h
      #   if (length(rancoef <- fullinit4LUarglist$init$ranCoefs[[1]])) skeleton$trRanCoefs=list("1"=attr(rancoef,"transf"))
      # }
      
      # Build template call 'template4objfn' used by .p4m_by_outer_optim(). 
      # The objective function in .p4m_by_outer_optim() itself calls *.p4m_by_iters()*, 
      # allowing iterations of the dynoffset in the objective fn;
      # so this cannot be an HLCorcall (which does not iterate). 
      # 'template4objfn' has to be modified from the initial .p4m_by_iters() call, including
      # .p4m_by_iters() local changes to it, such as the spaMM_tol$logL_tol)
      template4objfn <- mc # 
      # mc itself will be used to get some of the outer optim controls, then replaced later.
      template4objfn[["control"]][["p4m"]] <- substr(p4mcontrol,2,2) # 'H' expected
      if (missing(tol) || length(tol)>1L) {
        mc[["tol"]] <- tol[[2]] 
      } else  {
        mc[["tol"]] <- tol
      } 
      template4objfn["initfn"] <- NULL 
      # THis template4objfn typically already has an init.HLfit.
      template4objfn[["init.HLfit"]] <-  list(beta=fixef(mvp), v_h=mvp$v_h)
      
      ## 2ND STEP OPTIMISATION:
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
      
      ## FINAL REFIT
      # with fixed parameters to get the full fit object
      # Don't use trLambda directly: .preprocess() does nto handle it correctly. See details in .canonizeRanPars() 
      template4objfn[["fixed"]] <- .modify_list(oricall$fixed, 
                                                .canonizeRanPars(solution,corr_info=mvp$ranef_info$sub_corr_info, 
                                                                 checkComplete=FALSE, rC_transf=.spaMM.data$options$rC_transf)) 
      template4objfn["etaFix"] <- list(etaFix=solution$etaFix) # typically NULL ("H" refitting using init.HLfit instead)    
      template4objfn["init"] <- NULL 
      # Unfortunately I don't have a clean way of using dynoffset from .p4m_by_outer_optim() implied final fit
      mvp <- eval(template4objfn,parent.frame())  # still using $multinom_info in .makeMatp4m()
      if ( ! is.null(warnmess <- mvp$warnings$p4m_by_iters_not_conv)) 
        warning(warnmess, immediate. = TRUE)
      ## ad-hoc fixes, restoring features of original call
      mvp$dfs <- dfs # counting the dfs of the estimated parameters not of the final refit...
      mvp$call <- oricall
      # Final result should have same ZAL, v_h as an equivalent binomial fit and the .dynoffset implied by denom of NLPredictor.
      attr(mvp,"optimInfo") <- list(LUarglist=LUarglist, init.optim=skeleton,
                                    optim.pars=solution,
                                    objective=optim_blob$objective,
                                    rC_transf=.spaMM.data$options$rC_transf)
    } # else do nothing
  } else if (has_bar && # ranPars are all fixed
             p4mcontrol[1L]!="H") { 
    mc[["control"]][["p4m"]] <- 'H' # _____F I X M E_____ any way to automate and control this more finely ?
    # ideally the default p4m should be simply "H" when has_bar but not has_ranPars2fit
    # we would need to detect this case before the first step while currently has_ranPars2fit is evaluated after.
    mvp <- eval(mc, parent.frame()) # .p4m_by_iters # _____F I X M E____ not consistently spprec through all steps? (trace .solve_IRLS_as_ZX)
  }
  if ( is.list(sXaug <- mvp$envir$sXaug) && # -> spprec
       ! is.null(dcdb_p4m <- sXaug$AUGI0_ZX$dcdb_p4m)) { 
    mvp$envir$sXaug$AUGI0_ZX$dcdb_p4m <- .unscale(dcdb_p4m, scale=attr(mvp$X.pv,"scale_info")) 
    # $dcdb_p4m used post-fit in spprec to compute beta SE
  }
  if ( ! inherits(mvp,"HLfitlist") && ! is.call(mvp) ) {
    fit_time <- .timerraw(time1) 
    mvp$how$fit_time <- structure(fit_time,
                                  message="Please use how(<fit object>)[['fit_time']] to extract this information cleanly.")
  }
  
  mvp
}

predict.pois4mlogit <- function(object, newdata=NULL, ...) {
  .predict.pois4mlogit(object=object, newdata=newdata, ...) # this hides the private 'check' argument
}

.predict.pois4mlogit.check <- function(object, mnsizes, ...) {
  ndata <- object$data
  ndata$".dynoffset" <- ndata$".dynoffset" - log(mnsizes) # quick correction without predict.HLfit call.
  # This code was first conceived to avoid the double predict call in .predict.pois4mlogit(), 
  # but may be confusing because the predictions then 'exactly' sum to 1 only if convergence was 'exact'.
  predict.HLfit(object, newdata=ndata, ...)
}

.predict.pois4mlogit <- function(object, newdata=NULL, ..., check=FALSE) {
  if ( is.null(has_dynoffset <- object$p4m_info$multinom_info$has_dynoffset)) 
    stop("Old object, info not found where expected.")
  if (check) {
    # This code may be useful to provide some measure of convergence inaccuracy on prediction.
    if ( ! is.null(newdata)) stop("check=TRUE does not handle 'newdata'." ) 
    if ( is.null(mnsizes <- object$mnsizes)) stop("Old object, $mnsizes not found where expected.")
    .predict.pois4mlogit.check(object, mnsizes=mnsizes, ...)
  } else { # default case: double predict call, the first without the '...'
    # This double call makes sure that frequencies sum to 1, even if convergence was not 'exact'. 
    if (is.null(newdata)) {
      newdata <- object$data
      valid_rows <- object$p4m_info$multinom_info$valid_rows
    } else {
      if (is.null(newdata$".dynoffset")) newdata$".dynoffset" <- 0
      valid_rows <- ! is.na(newdata$".dynoffset")
    }
    newdata$".dynoffset"[valid_rows] <- 0
    pred <- predict.HLfit(object, newdata=newdata, na.action=na.exclude, ...) 
    pred <- matrix(pred, ncol=length(formula(object)))
    valid_rowsums <- rowSums(pred[valid_rows,has_dynoffset, drop=FALSE], na.rm =TRUE)
    newdata$".dynoffset"[valid_rows] <- - log(valid_rowsums)
    # Using newdata[valid_rows,,drop=FALSE] in next line may not be better (NA still possible in valid_rows) 
    predict.HLfit(object, newdata=newdata, ...)
  }
}

