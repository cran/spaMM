getCall.HLfit <- function(x, NbThreads=1L, ...) { ## FIXME ? getCall()$resid.model does not look like a language object (cf DHGLM)
  # The [stats::: !]getCall.default method cannot be called b/c it is not exported from stats.
  # stats::getCall() with only call getCall.HLfit in an infinite recursion.
  #
  if (x$spaMM.version> "2.4.62") {
    call <- x$call
    call[["control.HLfit"]]$NbThreads <- NbThreads
    call
  } else {
    # O L D E R version
    # only one of these call may be present in the object: HLCorcall is removed by fitme and corrHLfit
    if ( ! is.null(call <- attr(x,"fitmecall"))) return(call) 
    if ( ! is.null(call <- attr(x,"HLCorcall"))) return(call) ## eg confint on an HLCor object
    if ( ! is.null(call <- attr(x,"corrHLfitcall"))) return(call) 
    return(x$call) ## this one is the HLfit call
  }
}


##### OLD comment before I use fixed <- .modify_list(.)
## to get a call with the structure of the final HLCorcall in fitme or corrHLfit
## fixed is mandatory: Do not set a default value, so that one has to think about the correct value.
## Therefore, the original ranFix of the outer_object is replaced, unless it is explicitly set to getCall(object)$ranFix or $fixed... (in confint.HLfit)
## Parameters not in ranFix are set to the initial value of of the optimization call.
##   
## NB to get the $processed, it suffices to call a fitting function with verbose=c(getCall=TRUE)>...
#
get_HLCorcall <- function(outer_object, ## accepts fit object, or call, or list of call arguments
                          fixed, ## see comments above
                          ... # anything needed to overcome promises in the call
                          ) { 
  
  outer_call <- getCall(outer_object) ## corrHLfit/fitme/HLCor/HLfit/fitmv call
  outer_call$data <- outer_object$data ## removes dependence on promise
  outer_call$fixed <- .modify_list(outer_call$fixed, fixed)
  #
  outer_fn <-.get_bare_fnname.HLfit(outer_object, call.=outer_call)
  if (outer_fn=="corrHLfit" && length(init <- eval(outer_call[["init.corrHLfit"]]))) {
    init <- .reformat_ranPars(init, fitobject = outer_object)
    outer_call[["init.corrHLfit"]] <- remove_from_parlist(init, removand=fixed)
  } else if (length(init <- eval(outer_call[["init"]]))) { # user-level => Need to standardize it before trimming it.
    init <- .reformat_ranPars(init, fitobject = outer_object)
    outer_call[["init"]] <- remove_from_parlist(init, removand=fixed)
  }
  ## compare to update.default, commented in R language Definition.
  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(outer_call)))
    dotlist <- list(...)
    for (a in names(extras)[existing]) outer_call[[a]] <- dotlist[[a]]
    if (any(!existing)) {
      outer_call <- c(as.list(outer_call), dotlist[!existing])
    }
  }
  # *After* modifications by "extras" arguments, which may contain a user-provided 'verbose':
  verbose <- outer_call$verbose
  verbose["getCall"] <- TRUE # but this is automatically converted to an integer if there are integer elsewhere in the vector...
  outer_call$verbose <- verbose
  #
  HLCorcall <- eval(as.call(outer_call)) ## calls outer fn and bypasses any optimization to get the inner call HLCor/HLfit... / fitmv?
  HLCorcall$call <- NULL ## $call kept the outer call! 
  if (inherits(HLCorcall[[1]], "function")) { # if it is of class "name", no need for this block
    if (outer_fn=="HLfit") {
      HLCorcall[[1L]] <- quote(HLfit)
    } else if (outer_fn=="fitme" || outer_fn=="fitmv" ) {
      if ("control.dist" %in% names(formals(HLCorcall[[1]]))) {
        HLCorcall[[1L]] <- quote(HLCor)
      } else HLCorcall[[1L]] <- quote(HLfit)
    } else HLCorcall[[1L]] <- quote(HLCor)
  }
  if ("hyper" %in% names(HLCorcall$fixed)) {
    # then some of the steps performed by HLCor.obj in the original unconstrained fit 
    # are not performed when generating the HLCorcall with fixed params
    # notably .merge_fixed() -> .expand_hyper()
    HLCorcall$fixed <- .expand_hyper(HLCorcall$fixed, HLCorcall$processed$hyper_info,
                                     moreargs=.get_moreargs(outer_object)) 
  }
  .assignWrapper(HLCorcall$processed,"verbose['getCall'] <- NA")
  return(HLCorcall)
}

# update is a generic, update.formula is not a generic...
# stats::update.default is  function (object, formula., ..., evaluate = TRUE) 
# stats::update.formula is function (old, new, ...) 
# which  does not look like a method for # stats::update ! But it is. The src/library/stats/man/update.formula.Rd file  has
#  \method{update}{formula}(old, new, \dots)

#update <- function(object, formula., ..., evaluate = TRUE) UseMethod("update")
#update.default <- function(object, formula., ..., evaluate = TRUE) stats::update(object, formula., ..., evaluate = evaluate)

update.HLfit <- function(object, formula., ..., evaluate = TRUE) {
  if (is.null(call <- getCall(object))) 
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.)) {
    ## fixme A long time ago I wrote "does not handle etaFix$beta"... 
    if (is.null(data <- extras$data)) data <- object$data ## fortunately keeping more than the variables required in the original formula
    oriform <- formula.HLfit(object, which="hyper") ## with hyper-ranefs and offset
    if (inherits(oriform,"list")) { 
      if ( ! inherits(formula.,"list")) stop("Old formula is list (presumably from fitmv()) but new formula is not")
      for (mv_it in seq_along(formula.)) {
        formula.[[mv_it]] <- .update_formula(oriform[[mv_it]],formula.[[mv_it]], ...)
      }
      call[["formula."]] <- formula.
      # # the easy part would be updating the list of formulas
      # # the challenge would be to modify the formulas in the expression for the structured list of submodels...
      # if ( ! is.null(object$families) ) {
      #   stop("the 'formula.' argument cannot be used to modify a multivariate-response fit")
      # } else stop("Something wrong; the object is not a multivariate-response fit but its original formula is a list.")
    } else call$formula <- .update_formula(oriform,formula.) 
  }
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call))) ## which to replace and which to add to the call
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]] ## replace
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing]) ## add
      call <- as.call(call)
    }
  }
  if (evaluate) {
    updated <- eval(call, parent.frame())
    if (missing(formula.)) { # message misleading if we change between FixedM and MixedM
      spprec_ori <- .is_spprec_fit(object)
      spprec_upd <- .is_spprec_fit(updated)
      if (spprec_ori!=spprec_upd) message("Note: *one of* original or updated objects used sparse-precision methods, the other did not.")
      # this is useful to diagnose 'unreproducible' problems due to setting the spaMM option in one session but not in the other. 
    }
    updated
  } else call
}

.update_data <- function(object, re_data=object$data, mf=model.frame(object), newresp, respName=object$respName) {
  # mf <- model.frame(object)
  #form <- formula.HLfit(object, which="hyper")
  if (inherits(mf,"list")) { # list of data.frames (with attributes): mv case : then newresp is assumed to be the result of simulate, an overlong vector.
    attr(re_data,"validrownames") <- NULL # Will be presumably overwritten anyway, but let's play safe.
    nobs <- nrow(re_data)
    frst <- 0L
    respNames <- object$respNames
    for (mv_it in seq_along(mf)) {
      lst <- frst+nobs
      resp_i <- newresp[(frst+1L):lst]
      re_data <- .update_data(object, re_data=re_data, mf=mf[[mv_it]], newresp=newresp[(frst+1L):lst], respName=respNames[mv_it])
      frst <- lst
    }
    return(re_data)
  }
  Y <- model.response(mf) # its nrow cannot be used (pb in mv zero-truncated models) 
  if (NCOL(Y)==2L) { ## paste(form[[2L]])[[1L]]=="cbind"
    # model.frame is a data frame whose 1st element is a 2-col matrix with unnamed first column
    # model.response() extracts this matrix:
    # If formula is cbind(npos,nneg) ~... the two columns have names "npos", "nneg"
    # If formula is cbind(npos,ntot-npos) ~... the two columns have names "npos", ""
    # If formula is cbind(ntot-nneg,nneg) ~... the two columns have names "", "nneg"
    colY <- colnames(Y)
    if (colY[1L]!="") re_data[,colY[1L]] <- newresp
    if (colY[2L]!="") re_data[,colY[2L]] <- rowSums(Y)-newresp
    # any "ntot" col is left unchanged. In particular, from # cbind(ntot-nneg,nneg)
    #  thsi code changes nneg to ntot-newresp so that ntot-nneg will be newresp
  } else { # colnames(Y) is typically NULL
    ## : from a formula of the form formula I(<fn>(var...)) ~ ... colnames(mf)[1L] is "I(<fn>(var...))" 
    # for a variable of class 'AsIs' which is NOT used in the refit... 
    if (inherits(mf[[1L]],"AsIs")) {
      stop("the response of the original fit is as 'AsIs' term, I(.), which is not handled by code updating response.")
      # the alternative would be to change internally the lhs of the formula...
    } else {
      if (is.null(respName) || is.na(respName)) {
        resp_expr <- colnames(mf)[1L]
        if (resp_expr %in% colnames(re_data)) {
          re_data[resp_expr] <- newresp 
        } else stop(paste0("Response, ",resp_expr,", is not a variable in the 'data'.\n",
                           "See help('respName') for how to handle this case."))
      } else if (respName %in% colnames(re_data)) {
        re_data[respName] <- newresp 
      } else stop(paste0("Variable, ",respName,", is not a variable in the 'data'."))
    }
  }
  return(re_data)
}

.update_main_terms_info <- function(object, newresp) {
  mf <- model.frame(object)
  if (inherits(mf,"list")) { # list of data.frames (with attributes): mv case : then newresp is assumed to be the result of simulate, an overlong vector.
    vec_nobs <- object$vec_nobs
    frst <- 0L
    for (mv_it in seq_along(mf)) {
      lst <- frst+vec_nobs[mv_it]
      resp_range <- (frst+1L):lst
      resp_i <- newresp[resp_range]
      frst <- lst
      Y <- model.response(mf[[mv_it]])
      if (NCOL(Y) == 2L) {
        Y <- cbind(resp_i, object$BinomialDen[resp_range] - resp_i)
      } else Y <- resp_i
      mf[[mv_it]][1] <- Y
    }
  } else {
    Y <- model.response(mf)
    if (NCOL(Y)==2L) { ## paste(form[[2L]])[[1L]]=="cbind"
      # model.frame is a data frame whose 1st element is a 2-col matrix with unnamed first column
      Y <- cbind(newresp, object$BinomialDen-newresp)
    } else Y <- newresp # could use all.vars(formula(object)[-3]) to get a single var from an expression... but if expression has several variabls...
    mf[1] <- Y
  }
  main_terms_info <- object$main_terms_info # no longer true: 'of class "HLframes"'
  mf <- structure(mf, # model.frame, or list of them for mv
                  # fixef_terms=.get_from_terms_info(object, which="fixef_terms"), # already in main_terms_info
                  # fixef_levels=.get_from_terms_info(object, which="fixef_levels"),  # already in main_terms_info
                  fixefvarnames=.get_from_data_attrs(object, which="fixefvarnames"), 
                  fixefpredvars=.get_from_data_attrs(object, which="fixefpredvars"))
  main_terms_info$mf <- mf # adding $mf keeps the class (in contrast to subsetting)
  return(main_terms_info)    ############################# no longer true: 'of class "HLframes", the class set to object$HLframes' 
}


update_resp <- function(object, newresp, ...,  evaluate = TRUE) {
  if (is.null(re_call <- getCall(object))) stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if ("data" %in% names(extras)) stop("'data' not allowed in update_resp()'s arguments")
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(re_call))) ## which to replace and which to add to the call
    for (a in names(extras)[existing]) re_call[[a]] <- extras[[a]] ## replace
    if (any(!existing)) {
      re_call <- c(as.list(re_call), extras[!existing]) ## add
      re_call <- as.call(re_call)
    }
  }
  if (.spaMM.data$options$use_terms_info_attr && # trying to no longer use:
      (preprocess_handles_terms_info_attr <- TRUE) # handles it (except for silent mutations... )
     ) {
    # the data were not updated, to deal with case now handled by $respName[s] [see Details in help("update.HLfit")], 
    #      but other variables in the $data may still be used for resid model.
    re_call$data <- structure(object$data, updated_terms_info=.update_main_terms_info(object, newresp=newresp))
  } else re_call$data <- .update_data(object, newresp=newresp)
  if (evaluate) 
    eval(re_call, parent.frame())
  else re_call
}

# Fix intercept issue in local def of stats::terms.formula
.fixFormulaObject <- function (object) { # object is a formula object not a fit object
  Terms <- terms(object)
  tmp <- attr(Terms, "term.labels")
  ind <- grep("|", tmp, fixed = TRUE)
  if (length(ind)) 
    tmp[ind] <- paste("(", tmp[ind], ")")
  if (length(ind <- attr(Terms, "offset"))) {
    tmp2 <- as.character(attr(Terms, "variables"))[-1L]
    tmp <- c(tmp, tmp2[ind])
  }
  rhs <- if (length(tmp)) 
    paste(tmp, collapse = " + ")
  else "1"
  if (attr(Terms, "intercept")) rhs <- paste("1 +", rhs) ## opposite logic to R
  if (length(form <- formula(object)) > 2L) {
    res <- formula(paste("lhs ~", rhs))
    res[[2L]] <- form[[2L]]
    res
  }
  else formula(paste("~", rhs))
}


.update_formula <- function (old, new, ...) { 
  C_updateform <- get("C_updateform",asNamespace("stats"), inherits=FALSE) ## not kocher?
  tmp <- do.call(".Call",list(C_updateform, as.formula(old), as.formula(new))) # circumventing RcppExports' kind bureaucracy...  
  out <- formula(terms.formula(tmp, simplify = FALSE))
  out <- .fixFormulaObject(out)
  environment(out) <- environment(tmp)
  if (!inherits(out, "formula")) 
    class(out) <- c(oldClass(out), "formula")
  return(out)
}

#update.formula <- function(object, formula., ...) UseMethod("update.formula")
#update.formula.default <- function(object, formula., ...) stats::update.formula(old=object, new=formula., ...)

update_formulas <- function(object, formula., ...) {
  old <- formula(object)
  if (inherits(old,"list")) { # mv case
    if ( ! inherits(formula.,"list")) stop("Old formula is list (presumably from fitmv()) but new formula is not")
    for (mv_it in seq_along(old)) {
      old[[mv_it]] <- .update_formula(old[[mv_it]],formula.[[mv_it]], ...)
    }
    return(old)
  } else .update_formula(old, new=formula., ...)
}