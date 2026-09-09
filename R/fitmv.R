fitmv <- function(submodels, data, fixed=NULL, init=list(), lower=list(), upper=list(),
                  control=list(), # needed to avoid partial matching of explicit 'control' argument with 'control.dist' one (bug when the latter is used) 
                  control.dist = list(), method="ML", init.HLfit=list(), 
                  X2X=NULL, aliases=NULL, 
                  # multinom_info can be passed in the \dots
                  ...) { # \dots being arguments not requiring specific documentation.
  .spaMM.data$options$xLM_conv_crit <- list(max=-Inf)
  time1 <- Sys.time()
  oricall <- match.call(expand.dots=TRUE) ## mc including dotlist
  oricall <- ..n_names2expr(oricall) 
  oricall$"control.HLfit" <- eval(oricall$control.HLfit, parent.frame()) # to evaluate variables in the formula_env, otherwise there are bugs in waiting
  # where oricall[["control.HLfit"]] <- ... wouldn't work when 'control.HLfit' was absent. Same for 'fixed'
  oricall$"fixed" <- .preprocess_fixed(fixed)
  if (length(class(data))>1L) oricall$"data" <- as.data.frame(data) # such as tibble. 
  # See explanation in .preprocess() (which rechecks the effect on the data) 
  n_models <- length(submodels) # so the promise is already evaluated here...
  calls_W_processed <- fixedS <- vector("list",n_models)
  for (mv_it in seq_len(n_models)) { # call .preprocess() on each submodel
    call_ <- oricall
    if (length(aliases)) {
      for (varname in names(aliases)) {
        alias_it <- aliases[[varname]][mv_it]
        if ( ! is.na(alias_it)) {
          val <- data[[alias_it]]
          if (is.null(val)) {
            stop(paste0("Variable '",alias_it,"' not found in the data when processing 'aliases'."))
          } else data[[varname]] <- val
        }
      }
      call_[["data"]] <- data
    }
    call_["aliases"] <- NULL 
    call_["submodels"] <- NULL 
    ## I need to match the names of mv[[mit]] to those of a fitme call to make sure that they all named...
    call_["fixed"] <- NULL ## so that the lambda fixing (in particular) is not the default value for each processed call
    call_["X2X"] <- NULL ## otherwise detected as suspect arg by .preprocess_fitme()
    call_["control.HLfit"]$rankinfo <- call_["control.HLfit"]$rankinfo$mvlist[[mv_it]]
    call_["etaFix"] <- NULL ## not the right step for fixing coefficients.
    ## *** global arguments => avoid mixing them with local arguments 
    ##     (although this is stricly necessary only for covStruct since...) ***  
    call_["corrMatrix"] <- NULL # not strictly necess since single matrix so never a problem of matching ranefs: The global corrMatrix is a locally usable corrMatrix, 
    call_["adjMatrix"] <- NULL # not strictly necess ... same comment...
    call_["covStruct"] <- NULL # => important to remove it since ranefs cannot be matched in .preprocess().
    #
    call_["init.HLfit"] <- NULL # useless for submodel .preprocess()ing IF .merge_processed() 
    # calls .check_init.HLfit(init.HLfit) using the oricall's init.HLfit, and if
    # the fitmv_body() call receives its init.HLfit$ from merged$"init_HLfit" 
    ## *** ***  
    matched_args_it <- 
      match.call(fitme, 
                 do.call("call",c(list(name="fitme"), submodels[[mv_it]]), quote=TRUE)) # match the elements of mv[[mv_it]] to those of a call to fitme
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
    # call_["submodel"] <- mv_it
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
      call_[["data"]] <- structure(call_[["data"]], updated_terms_info=main_terms_info_it)
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
  mc <- oricall 
  mc[c("aliases","submodels","formula.")] <- NULL # so that it remains in call_ the arguments others than mv.
  mc[["what_checked"]] <- "fitmv() call" 
  mc[[1L]] <- get(".check_args_fitme", asNamespace("spaMM"), inherits=FALSE) 
  eval(mc,parent.frame()) # -> abyss 
  mc[c("what_checked", "fixed","upper","lower","control", "multinom_info")] <- NULL # but user, lower in fitme_body call
  mc[["calls_W_processed"]] <- calls_W_processed
  # the fact that promises are evaluated within a call-execution is "local": they will appear not evaluated
  # when we reuse a call (here mc). E.g. corrMatrix=as_precision(.) would be evaluated twice 
  # => We need to put the evaluated value in the call list. 
  # Next line ad-hoc for corrMatrix (_F I X M E__?: What about other arguments ? Which would benefit from some preprocessing?)
  if ("corrMatrix" %in% ...names()) mc["corrMatrix"] <- list(eval(mc[["corrMatrix"]])) 
  mc[[1L]] <-  get(".merge_processed", asNamespace("spaMM"), inherits=FALSE)
  merged <- eval(mc, parent.frame()) # means that arguments of *.merge_processed()* must have default values as mc does not contains defaults of fitmv()
  
  # In p4m code: multinom_info with/out $mnsizes argument provided to fitmv() depending on p4m='H'/'o'.
  # so, here fitmv code: either no input multinom_info() or two types of input multinom_info, w/o $mnsizes. 
  # There is *always* an *output* $multinom_info:
  if (is.null(merged$multinom_info <- eval(oricall$multinom_info))) merged$multinom_info <- 
      list(has_dynoffset=rep(FALSE, n_models))
  #  
  fixed <- .reformat_parlist(fixed,processed = merged) # reformat user's global 'fixed' argument
  fixedS <- lapply(fixedS, .reformat_parlist, processed = merged)
  fixedS <- .merge_mv_parlist(fixedS, merged) # now fixedS is a single parlist from the  sub-models specifications
  fixedS <- .modify_list(fixedS,fixed) # now fixedS is a single parlist from both sub-model and global specifications
  fixedS <- .preprocess_fixed(fixedS)
  
  # These infos are ultimately used by summary() to distinguish "fix" from outer "var":
  merged[["lambda.Fix"]] <- .reformat_lambda(.getPar(fixed,"lambda"), processed=merged, full_lambda=TRUE)
  # HLfit_body() expects merged[["phi.Fix]] to be a full-length list, possibly with explicit NULLs.
  # merged[["phi.Fix"]] from .merge_processed() should be so, and .modify_list() should keep it so.
  merged[["phi.Fix"]] <- .modify_list(merged[["phi.Fix"]], fixedS$phi)
  #
  ranCoefs <- .getPar(fixedS,"ranCoefs") ## may be NULL
  merged$ranCoefs_blob <- .process_ranCoefs(merged, ranCoefs, use_tri_CORREL=TRUE) 
  merged$AUGI0_ZX$envir$finertypes[merged$ranCoefs_blob$isRandomSlope] <- "ranCoefs" 
  #
  mc["upper"] <- oricall["upper"]
  mc["lower"] <- oricall["lower"]
  mc["control"] <- oricall["control"]
  mc[["fixedS"]] <- fixedS # to build and merge the inits
  mc$processed <- merged
  not_in_fitmv_body <- c("init.HLfit", # fitmv_body directly use the processed$init_HLfit version; 
              # it would be confusing to suggest otherwise by keeping the arg.
              "calls_W_processed","data","family","prior.weights", "weights.form", 
              "HLmethod","method","rand.family","control.glm","REMLformula",
              "resid.model", "verbose","distMatrix","adjMatrix", "control.dist", "corrMatrix","covStruct","X2X") 
  mc[not_in_fitmv_body] <- NULL 
  mc[[1L]] <-  get("fitmv_body", asNamespace("spaMM"), inherits=FALSE)
  hlcor <- eval(mc,parent.frame()) 
  # if (.safe_true(processed[["verbose"]]["get_LUarglist"][[1L]])) return(hlcor)
  oricall$"control.dist" <- merged[["control_dist"]] 
  hlcor$call <- oricall ## this is a call to fitmv()
  lsv <- c("lsv",ls())
  if (is.call(hlcor)) {
    # ...
  } else if ( ! inherits(hlcor,"HLfitlist")) {
    X2X <- eval(oricall[["X2X"]], parent.frame())
    if (inherits(X2X,"call")) { # genX2X call
      if (deparse(X2X[[1]])=="genX2X") { # ____F I X M E____ allow user-def'd function ?
        X2X[["names_ori"]] <- attr(hlcor$X.pv,"cols_lhs_X2X")
        X2X <- eval(X2X) 
      } else warning("Fit object's 'X2X' element remains a call: this may be a problem in post-fit operations such as predict().")
    }
    hlcor$X2X <- X2X
    hlcor$aliases <- eval(oricall[["aliases"]], parent.frame())
    hlcor$how$fit_time <- .timerraw(time1)
    hlcor$how$fnname <- "fitmv"
    hlcor$fit_time <- structure(hlcor$how$fit_time,
                                message="Please use how(<fit object>)[['fit_time']] to extract this information cleanly.")
    if ( ! is.null(mc$control.HLfit$NbThreads)) .setNbThreads(thr=.spaMM.data$options$NbThreads)
    class(hlcor) <- c("fitmv", class(hlcor))
  }
  rm(list=setdiff(lsv,"hlcor")) ## empties the whole local envir except the return value
  return(hlcor)
}

