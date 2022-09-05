.is_2_in_1 <- function(X1,X2, tol=1e-10) {
  qrX <- qr(X1)
  proj2in1 <- X1 %*% qr.solve(qrX,X2)
  return(all(abs(proj2in1 -X2)<tol))
}

# deprecated:
.process_ranef_case <- function(object, object2, nest="") {
  l1 <- logLik(object)
  l2 <- logLik(object2)
  testlik <- unique(c(names(l1),names(l2)))
  if (length(testlik)!=1L) stop(paste("Unable to determine a unique objective function from logLik() names.\n",
                                      "Check that both fits are fitted by the same 'method'."))
  if (nest=="2in1" || (nest!="1in2" && l1>l2)) {
    nullfit <- object2
    fullfit <- object
  } else { # nest=="1in2" || l2>l1
    nullfit <- object
    fullfit <- object2
  }
  return(list(fullfit=fullfit,nullfit=nullfit,test_obj=testlik,df=NA))
}

.guess_Rnest <- function(object, object2, Xnest) {
  dfR1 <- sum(unlist(object$dfs, use.names=FALSE)) - object$dfs$pforpv 
  dfR2 <- sum(unlist(object2$dfs, use.names=FALSE)) - object2$dfs$pforpv
  if (dfR1>dfR2) {
    Rnest <- "2in1"
    if (is.null(Xnest)) { # fixed effects are identical
      message(paste("spaMM could not ascertain whether the models are nested in their random-effect specifications.\n", 
                    "  You are responsible for that.")) 
    } else warning(paste("Fixed-effects specifications differ between the models, and spaMM could not ascertain\n",
                         "  whether the models are similarly nested in their random-effect specifications.\n",
                         "  You are responsible for that. spaMM will guess nestedness from number of parameters."), immediate. = TRUE) 
  } else if (dfR1<dfR2) { 
    Rnest <- "1in2"
    if (is.null(Xnest)) { # fixed effects are identical
      message(paste("spaMM could not ascertain whether the models are nested in their random-effect specifications.\n", 
                    "  You are responsible for that.")) 
    } else warning(paste("Fixed-effects specifications differ between the models, and spaMM could not ascertain\n",
                         "  whether the models are similarly nested in their random-effect specifications.\n",
                         "  You are responsible for that. spaMM will guess nestedness from number of parameters."), immediate. = TRUE) 
  } else if (is.null(Xnest)) { # fixed effects are identical, dfs are identical
    stop("The models compared appear to have the same number of parameters.") 
  } else {
    Rnest <- NULL
    warning(paste("The models have the same number of parameters except for fixed effects, but spaMM could not ascertain\n",
                  "  whether the random-effect specifications are identical. You are responsible for that."), immediate. = TRUE) 
  }
  Rnest
}




.compare_model_structures <- function(object,object2) {
  if (inherits(object,"HLfitlist") || inherits(object2,"HLfitlist")) {
    stop("This does not yet work on HLfitlist objects")
  }
  X1 <- attr(object$`X.pv`,"namesOri") ## need to track NA beta's
  X2 <- attr(object2$`X.pv`,"namesOri")
  if (length(X1)==0L) {
    REML1 <- NULL ## compatible with both ML or REML tests
  } else REML1 <- (object$APHLs$p_v != object$APHLs$p_bv)
  if (length(X2)==0L) {
    REML2 <- NULL ## idem
  } else REML2 <- (object2$APHLs$p_v != object2$APHLs$p_bv)
  REML <- unique(c(REML1,REML2))
  meth1 <- object$HL
  meth2 <- object2$HL
  if (! identical(object$family[c("family","link")],object2$family[c("family","link")] ) ) {
    stop("Models may not be nested (distinct families).") ## but COMPoisson vs poisson ?
  }
  if (! identical(meth1,meth2) || length(REML)>1 ) {
    stop("object fitted by different methods cannot be compared.")
  }
  if ( ! is.null(X1)) X1 <- sapply(strsplit(X1,':'), function(x) paste(sort(x),collapse=':')) ## JBF 2015/02/23: sort variables in interaction terms before comparison
  if ( ! is.null(X2)) X2 <- sapply(strsplit(X2,':'), function(x) paste(sort(x),collapse=':'))
  dX12 <- setdiff(X1,X2)
  dX21 <- setdiff(X2,X1)
  if (length(dX12) && length(dX21)) {
    if (.is_2_in_1(X1=object$X.pv,  X2=object2$X.pv)) {
      Xnest <- "2in1"
    } else if (.is_2_in_1(X1=object2$X.pv,  X2=object$X.pv)) {
      Xnest <- "1in2"
    } else stop("Fixed effects seem non-nested.") 
  } else if (length(dX12)) {
    Xnest <- "2in1"
  } else if (length(dX21)) {
    Xnest <- "1in2"
  } else {
    Xnest <- NULL
  }
  if (object$spaMM.version < "2.2.116") {
    ranterms1 <- attr(object$ZAlist,"ranefs")
  } else ranterms1 <- attr(object$ZAlist,"exp_ranef_strings")
  if (object2$spaMM.version < "2.2.116") {
    ranterms2 <- attr(object2$ZAlist,"ranefs")
  } else ranterms2 <- attr(object2$ZAlist,"exp_ranef_strings")
  randist1 <- lapply(object$rand.families, function(v) paste(paste(v)[1:2],collapse="")) ## makes a string from each $family and $link 
  randist2 <- lapply(object2$rand.families, function(v) paste(paste(v)[1:2],collapse="")) ## makes a string from each $family and $link 
  ranterms1 <- paste(ranterms1,randist1) ## joins each term and its distrib
  ranterms2 <- paste(ranterms2,randist2) ## joins each term and its distrib
  dR12 <- setdiff(ranterms1,ranterms2)
  dR21 <- setdiff(ranterms2,ranterms1)
  if (length(dR12) && length(dR21)) { # no obvious nested structure on ranefs. Trying to guess 'Rnest'
    Rnest <- .guess_Rnest(object, object2, Xnest)
    #.process_ranef_case(object, object2, Xnest=Xnest)
  } else if (length(dR12)) {
    Rnest <- "2in1"
  } else if (length(dR21)) {
    Rnest <- "1in2"
  } else Rnest <- NULL
  nest <- c(Xnest,Rnest)
  unest <- unique(nest)
  if (length(unest)==2L) stop("Models not nested (opposite nestings for fixed and random terms). ")
  # ELSE
  if (length(unest)==0L) stop(paste("The two models appear equivalent (except perhaps for residual dispersion models).\n", 
                                    "This case is not handled."))
  # if (length(Rnest)) return(.process_ranef_case(object, object2, nest=Rnest)) # Possibly nested models, differing at least by their random effects.
  ###############################
  # # ELSE nested fixef, identical ranefs
  # df1 <- length(X1[!is.na(fixef(object))])
  # df2 <- length(X2[!is.na(fixef(object2))])
  # if (!is.null(Rnest)) {
  #   lambda.object <- object$lambda.object
  #   if (!is.null(lambda.object)) df1 <- df1+length(unlist(lambda.object$coefficients_lambdaS))
  #   cov.mats <- .get_compact_cov_mats(object$strucList)
  #   if (length(cov.mats)) {
  #     nrows <- unlist(lapply(cov.mats,NROW))
  #     df1 <- df1+sum(nrows*(nrows-1)/2)
  #   }
  #   lambda.object <- object2$lambda.object
  #   if (!is.null(lambda.object)) df2 <- df2+length(unlist(lambda.object$coefficients_lambdaS))
  #   cov.mats <- .get_compact_cov_mats(object2$strucList)
  #   if ( length(cov.mats)) {
  #     nrows <- unlist(lapply(cov.mats,NROW))
  #     df2 <- df2+sum(nrows*(nrows-1)/2)
  #   }
  # }
  df1 <- sum(unlist(object$dfs, use.names=FALSE)) # but recursive=TRUE bc for fitmv the $dfs are a structured list
  df2 <- sum(unlist(object2$dfs, use.names=FALSE))
  if (unest=="1in2") {
    fullm <- object2
    nullm <- object
    df <- df2-df1
  } else {
    fullm <- object
    nullm <- object2
    df <- df1-df2
  }
  if (length(nest)==2) {
    message("Models differing both by in their fixed and in their random terms. ")
    message("Tentatively using marginal likelihood to compare them... ")
    testlik <- "p_v" 
  } else {
    if (is.null(Rnest)) { ## fixed effect test 
      if (REML) {
        ## checking the comparability of REML fits
        if ( ! is.null(fullm$distinctX.Re) ) {
          df.f.Re <- ncol(fullm$distinctX.Re)
        } else df.f.Re <- ncol(fullm$`X.pv`)
        if ( ! is.null(nullm$distinctX.Re) ) {
          df.n.Re <- ncol(nullm$distinctX.Re)
        } else df.n.Re <- ncol(nullm$`X.pv`)
        if ( df.f.Re !=  df.n.Re ) {
          warning("LRT comparing REML fits with different fixed-effect conditions is highly suspect")
        }
      }
      testlik <- "p_v"
    } else { ## random effect test
      if ( ! REML) warning("ML fits used to compare different random-effects models...")
      testlik <- "p_bv" ## used in both case, identical to p_v in the non-standard case
    }
  } 
  return(list(fullfit=fullm,nullfit=nullm,test_obj=testlik,df=df))
}

.get_outer_inits_from_fit <- function(fitobject, keep_canon_user_inits) {
  canon.init <- attr(fitobject,"optimInfo")$LUarglist$canon.init ## includes user init
  #
  outer_ests <- get_ranPars(fitobject) ## CorrEst_and_RanFix only
  names_u_o_ests <- names(unlist(outer_ests))
  names_u_c_inits <- names(unlist(canon.init))
  # For each of the following family param there may already be attr(outer_ests,"type")$tr<fam par> but let's not assume that messy thing
  if ("NB_shape" %in% names_u_c_inits) {
    outer_ests$NB_shape <- environment(fitobject$family$aic)$shape
    attr(outer_ests,"type")$NB_shape <- "outer" 
  }
  if ("COMP_nu" %in% names_u_c_inits) {
    outer_ests$COMP_nu <- environment(fitobject$family$aic)$nu
    attr(outer_ests,"type")$COMP_nu <- "outer" 
  }
  if ("beta_prec" %in% names_u_c_inits) {
    outer_ests$beta_prec <- environment(fitobject$family$aic)$prec
    attr(outer_ests,"type")$beta_prec <- "outer" 
  }
  if (keep_canon_user_inits) { # keep them (as interpreted in canon.init: minimum phi is 1e-4, etc) in return value
    # => remove the fitted values from the nullranPars used to modify_list
    # => keep them in 'removand' list of pars to remove from nullranPars
    # => exclude them from 'not_user_inits_names' to remove from 'removand' !
    user_inits <- .reformat_corrPars(getCall(fitobject)$init,corr_families=fitobject$corr_info$corr_families)
    names_u_u_inits <- names(unlist(user_inits))
    not_user_inits_names <- setdiff(names_u_c_inits, names_u_u_inits) # names, excluding those of parameters with user inits
    removand <- setdiff(names_u_o_ests, not_user_inits_names) ## removand: user_inits, fixed, or inner optimized corrPars
    ## removand: user_inits, fixed, or inner optimized corrPars
    # locinit will retain parameters that were outer optimized without an explicit user init
  } else removand <- setdiff(names_u_o_ests, names_u_c_inits)
  if ( is.null(removand)) {
    locinit <- .modify_list(canon.init, outer_ests)
  } else { ## leaves user_inits as there are in LUarglist$canon.init, and do not add any fixed or inner-optimized par
    locinit <- .modify_list(canon.init,
                            .remove_from_cP(outer_ests, u_names=removand)) ## loses attributes
  }
  return(locinit)
}

.get_rC_inits_from_hlfit <- function(hlfit, type) {
  reinit <- .get_compact_cov_mats(hlfit$strucList)
  seq_rd <- seq_along(reinit)
  if (length(seq_rd)) {
    for (rd in seq_rd) {
      if ( ! is.null(ini_mix <- reinit[[rd]])) {
        ini_corr <- cov2cor(ini_mix)
        lowerblocF <- lower.tri(ini_corr,diag=FALSE)
        ini_mix[lowerblocF] <-  ini_corr[lowerblocF] # replaces covby corr
        reinit[[rd]] <- ini_mix[lower.tri(ini_mix,diag=TRUE)] ## mix cov/corr in vector form
      } else reinit[rd] <- NA # the syntax understood by fitting functions
    }
    names(reinit) <- paste(seq_rd)
    if ( ! is.null(type)) {
      reinit[hlfit$lambda.object$type!=type] <- NA
      reinit <- reinit[ ! is.na(reinit)] # otherwise for "outer" type the NA ends in the optimizer's init...
    }
  }
  return(reinit)
}

# for inner CAR it returns the single lambda factor which is OK 
.get_lambdas_notrC_from_hlfit <- function(hlfit, type, keep_names=(type=="adhoc"), 
                                          isRandomSlope=attr(hlfit$strucList,"isRandomSlope")) {
  lambdas <- hlfit$lambda.object$lambda_list
  seq_rd <- seq_along(lambdas)
  if (length(seq_rd)) {
    any_ranCoef <- FALSE
    for (rd in seq_rd) {
      if (isRandomSlope[[rd]]) {
        lambdas[[rd]] <- NA # the syntax understood by fitting functions
        any_ranCoef <- TRUE  
      } else {
        # pretty renaming of simple lambda:
        if (type=="adhoc")  if (length(lambdas[[rd]])==1L && names(lambdas[[rd]])=="(Intercept)") names(lambdas[[rd]]) <- NULL
      }
    }
    lambdas <- unlist(lambdas)
    if ( ! keep_names) names(lambdas) <- paste(seq_rd)
    if (type !="adhoc") lambdas[hlfit$lambda.object$type!=type] <- NA
    lambdas <- lambdas[ ! is.na(lambdas)] # otherwise for "outer" type the NA ends in the optimizer's init...
    # for "inner" NA's are harmless but not necessary when names are paste(seq_rd)
    if (any_ranCoef && type=="adhoc") lambdas <- structure(lambdas,
                                           message="Random-coefficient variances removed. Use e.g. VarCorr() to get them.")
  }
  return(lambdas) # vector wwith NAs for everything not wanted
}

# to initiate a *f*ullfit                                     # , inner_lambdas= TRUE
get_inits_from_fit <- function(from, template=NULL, to_fn=NULL, inner_lambdas=FALSE) { # 'to_fn' may differ from that of 'from' and 'to'
  new_outer_inits <- .get_outer_inits_from_fit(fitobject=from, keep_canon_user_inits = FALSE)
  # check fromfn and to_fn
  if (is.null(to_fn)) {
    # recent objects should have how$fnname. Otherwise, .get_bare_fnname() may not return a valid name.
    fromfn <- .get_bare_fnname.HLfit(from)
    if (is.null(template)) {
      to_fn <- fromfn
    } else if (inherits(template,"HLfit")) {
      to_fn <- .get_bare_fnname.HLfit(template)
    } else stop("Invalid 'template' argument.")
  } else if (to_fn=="fitme_body") { ## using to_fn to modify fromfn...
    ## ad hoc fix for residModel: fitme_body is called directly so the final object's call is to HLCor of HLfit
    fromfn <- "fitme"
  } else fnname <- .get_bare_fnname.HLfit(from)
  # Inner-estimated lambda and ranCoefs (FIXME could add phi)
  init.HLfit <- NULL
  rC_inner_inits <- .get_rC_inits_from_hlfit(from, type="inner")
  if (length(rC_inner_inits) ) init.HLfit <- list(ranCoefs=rC_inner_inits)
  if (inner_lambdas) {
    lambda_inner_inits <- .get_lambdas_notrC_from_hlfit(from, type="inner")
    if (length(lambda_inner_inits)) init.HLfit <- c(init.HLfit, list(lambda=lambda_inner_inits))
  } 
  #
  if (to_fn %in% c("fitme", "fitmv", "fitme_body")) {
    new_inits <- list(init=new_outer_inits,init.HLfit=init.HLfit)
  } else if (to_fn=="corrHLfit") {
    new_inits <- list(init.corrHLfit=new_outer_inits,init.HLfit=init.HLfit)
  } else new_inits <- list(init.HLfit=init.HLfit)
  # Add initial value for fixed effects
  if (length(fixef_from <- na.omit(fixef(from)))) {
    # fixef() returns a vector with NA's ifor non-estimable parameters; 
    #     but these NA should not reach the code initializing eta as X.beta in HLfit_body, using beta from the inits... => na.omit
    if (inherits(template,"HLfit")) { # This was motivated by the Leucadendron_hard.R bootstrap replicates.
      new_HLfit_inits <- na.omit(fixef(template))
      new_HLfit_inits[names(new_HLfit_inits)] <- 0
      shared_names <- intersect(names(new_HLfit_inits), names(fixef_from))
      new_HLfit_inits[shared_names] <- fixef_from[shared_names]
      new_HLfit_inits <- list(fixef=new_HLfit_inits)
      new_inits[["init.HLfit"]] <- c(new_inits[["init.HLfit"]], new_HLfit_inits)
    } else new_inits[["init.HLfit"]] <- c(new_inits[["init.HLfit"]], list(fixef=fixef_from))
  }
  return(new_inits)
}

.update_control <- function(fit_call, optim_boot, from_fn=NULL) {
  if (is.null(from_fn)) from_fn <- paste(fit_call[[1]])
  if (from_fn=="fitme") {
    ctrl_opt <- fit_call[["control"]]
    if (is.null(ctrl_opt)) {
      ctrl_opt <- list(optimizer=optim_boot) 
    } else ctrl_opt[["optimizer"]] <- optim_boot
  } else if (from_fn=="corrHLfit") {
    ctrl_opt <- fit_call[["control.corrHLfit"]]
    if (is.null(ctrl_opt)) {
      ctrl_opt <- list(optimizer=optim_boot) 
    } else ctrl_opt[["optimizer"]] <- optim_boot
  } else ctrl_opt <- NULL
  return(ctrl_opt)
}

eval_replicate <- function(y) { # no additional arguments, to ease parallel programming => next lines instead
  # the function will be called within e.g. pbapply so it's useless to refer to parent.frame() here
  enclosing_env <- parent.env(environment()) ## this is not necess for the code to run, but for the CRAN checks not to complain
  nullfit <- get("nullfit", enclosing_env)
  fullfit <- get("fullfit", enclosing_env)
  #dotargs <- get("dotargs", enclosing_env)
  test_obj <- get("test_obj", enclosing_env)
  debug. <- get("debug.", enclosing_env)
  #  .condition <- get(".condition", enclosing_env)
  null_call <- getCall(nullfit)
  null_fit_fn <- .get_bare_fnname.HLfit(nullfit, call.=null_call)
  full_fit_fn <- .get_bare_fnname.HLfit(fullfit) 
  newinits <- .get_outer_inits_from_fit(fitobject=nullfit, keep_canon_user_inits = FALSE) # for next new_nullfit
  ctrl_opt <- .update_control(fit_call=null_call, optim_boot=.spaMM.data$options$optim_boot, from_fn=null_fit_fn) # need .safe_opt when newinits are at bound.
  if (null_fit_fn=="fitme") {
    new_args <- list(init=newinits, control=ctrl_opt)
  } else if (null_fit_fn=="corrHLfit") {
    new_args <- list(init.corrHLfit=newinits, control.corrHLfit=ctrl_opt)
  } else new_args <- NULL
  # pbbly never good not to use try(). It's the handling of try-error that may vary
  # debug(pbapply) (or pblapply ?) may be useful to 
  if (debug.==2) {# Shuld be prevented by spaMM's calling fns in a parallel session   
    re_nullfit <- do.call(update_resp, c(list(object=nullfit, newresp = y),new_args))
  } else {
    re_nullfit <- try(do.call(update_resp, c(list(object=nullfit, newresp = y),new_args)))
    if (inherits(re_nullfit,"try-error")) { ## (debug.= TRUE or 1L) to return error info in parallel mode: return the try-error object
      if (debug.) {
        utils::dump.frames(dumpto="dump_on_re_nullfit", to.file=TRUE) # but doesn't stop
        return(list(full=structure(NA,info="no fullfit"), 
                    null=structure(NA,re_nullfit=re_nullfit)))
      } else return(c(full=NA, null=NA)) # attributes would be lost at the level of the pbapply() closure. 
      # cf apply -> array -> as.vector -> loses all attributes as doc'ed in apply and as.vector
    } ## ELSE continue
  }
  logL_re_null <- logLik(re_nullfit,which=test_obj)
  # Allows the user to control the starting values of the re_fullfit:
  # edotargs <- dotargs
  # for (st in setdiff(names(edotargs),"prior.weights")) {
  #   edotargs[[st]] <- eval(edotargs[[st]],env=environment()) ## evaluate the promises in the current execution envir
  # }
  # #
  newinits <- get_inits_from_fit(from=re_nullfit, template=fullfit, to_fn=full_fit_fn )
  lens <- rep(NA, length(newinits))
  for (it in seq_along(newinits)) lens[it] <- length(newinits[[it]])
  newinits <- newinits[lens>0L]
  subsets_inits <- .all.subsets(newinits) ## limited scope: could effectively do this is a nested way (with a skeleton+relist technique)
  # subsets_inits is always a list of lists, with at least one element (minimally: list(list()))
  for (it in seq_along(subsets_inits)) {
    inits <- subsets_inits[[it]]
    if (full_fit_fn=="fitme") {
      new_args <- list(init=inits[["init"]], init.HLfit=inits[["init.HLfit"]], control=ctrl_opt)
    } else if (full_fit_fn=="corrHLfit") {
      new_args <- list(init.corrHLfit=inits[["init.corrHLfit"]], init.HLfit=inits[["init.HLfit"]], control.corrHLfit=ctrl_opt)
    } else new_args <- list(init.HLfit=inits[["init.HLfit"]])
    if (debug.==2) {
      re_fullfit <- do.call(update_resp, c(list(object=fullfit, newresp = y),new_args)) # may stop on error
    } else {
      re_fullfit <- try(do.call(update_resp, c(list(object=fullfit, newresp = y),new_args)))
      if (inherits(re_fullfit,"try-error")) { ## (debug.= TRUE or 1L) to return error info in parallel mode: return the try-error object
        if (debug.) {
          utils::dump.frames(dumpto="dump_on_re_fullfit", to.file=TRUE) # but doesn't stop
          return(list(full=structure(NA,re_fullfit=re_fullfit), # but the attributes seem lost at the level of the pbapply() closure.
                      null=structure(logL_re_null,re_nullfit=re_nullfit)))
        } else return(c(full=NA, null=logL_re_null)) # attributes would be lost at the level of the pbapply() closure. 
      } ## ELSE continue
    }
    logL_re_full <- logLik(re_fullfit,which=test_obj)
    if (logL_re_full>logL_re_null-1e-04) break # break the for loop
  }
  resu <- c(full=logL_re_full,null=logL_re_null)
  # if ( ! is.null(.condition)) {
  #   condition <- with(list(nullfit=re_nullfit, fullfit=re_fullfit), eval(.condition))
  #   resu <- c(resu, condition=condition)
  # }
  return(resu)
}


.eval_replicate2 <- function(y) { 
  enclosing_env <- parent.env(environment())
  nullfit <- get("nullfit",enclosing_env)
  fullfit <- get("fullfit",enclosing_env)
  # dotargs <- get("dotargs",enclosing_env)
  test_obj <- get("test_obj",enclosing_env)
  debug. <- get("debug.", enclosing_env)
#  .condition <- get(".condition", enclosing_env)
  null_call <- getCall(nullfit)
  null_fit_fn <- .get_bare_fnname.HLfit(nullfit, call.=null_call)
  full_fit_fn <- .get_bare_fnname.HLfit(fullfit)
  conv_full <- conv_null <- FALSE
  best_logL_full <- best_logL_null <- prev_logL_full <- prev_logL_null <- -Inf
  # Allows the user to control the starting values of the initial new_nullfit
  # edotargs <- dotargs
  # for (st in setdiff(names(edotargs),"prior.weights")) {
  #   edotargs[[st]] <- eval(edotargs[[st]],env=environment()) ## evaluate the promises in the current execution envir
  # }
  ctrl_opt <- .update_control(fit_call=null_call, optim_boot=.spaMM.data$options$optim_boot, from_fn=null_fit_fn) # need .safe_opt when newinits are at bound.
  newinits <- .get_outer_inits_from_fit(fitobject=nullfit, keep_canon_user_inits = FALSE) # for next new_nullfit
  if (null_fit_fn=="fitme") {
    new_args <- list(init=newinits, control=ctrl_opt)
  } else if (null_fit_fn=="corrHLfit") {
    new_args <- list(init.corrHLfit=newinits, control.corrHLfit=ctrl_opt)
  } else new_args <- NULL
  while ( TRUE ) {
    if (debug.==2) {# Shuld be prevented by spaMM's calling fns in a parallel session   
      new_nullfit <- do.call(update_resp, c(list(object=nullfit, newresp = y),new_args))
    } else {
      new_nullfit <- try(do.call(update_resp, c(list(object=nullfit, newresp = y),new_args)))
      if (inherits(new_nullfit,"try-error")) { ## (debug.= TRUE or 1L) to return error info in parallel mode: return the try-error object
        if (debug.) {
          utils::dump.frames(dumpto="dump_on_new_nullfit", to.file=TRUE) # but doesn't stop
          return(list(full=structure(NA,new_nullfit=new_nullfit), # but the attributes seem lost at the level of the pbapply() closure.
                      null=structure(NA,info="no nullfit")))
        } else return(c(full=NA, null=NA)) # attributes would be lost at the level of the pbapply() closure. 
      } ## ELSE continue
    }
    logL_new_null <- logLik(new_nullfit,which=test_obj)
    #cat(logL_new_null)
    conv_null <- (abs(logL_new_null - prev_logL_null)<1e-4)
    if (logL_new_null>best_logL_null) { # always true the first time
      best_logL_null <- logL_new_null
      best_nullfit <- new_nullfit
    }
    if (conv_null) break # no point in refitting the full model if the new inits from null fit don't change
    # ELSE
    prev_logL_null <- logL_new_null
    newinits <- get_inits_from_fit(from=best_nullfit, template=fullfit, to_fn=full_fit_fn) # using default 'inner_lambdas' argument
    lens <- rep(NA, length(newinits))
    for (it in seq_along(newinits)) lens[it] <- length(newinits[[it]])
    newinits <- newinits[lens>0L]
    subsets_inits <- .all.subsets(newinits) ## limited scope: could effectively do this is a nested way (with a skeleton+relist technique)
    # subsets_inits is always a list of lists, with at least one element (minimally: list(list()))
    for (it in seq_along(subsets_inits)) {
      inits <- subsets_inits[[it]]
      if (full_fit_fn=="fitme") {
        new_args <- list(init=inits[["init"]], init.HLfit=inits[["init.HLfit"]], control=ctrl_opt)
      } else if (full_fit_fn=="corrHLfit") {
        new_args <- list(init.corrHLfit=inits[["init.corrHLfit"]], init.HLfit=inits[["init.HLfit"]], control.corrHLfit=ctrl_opt)
      } else new_args <- list(init.HLfit=inits[["init.HLfit"]])
      if (debug.==2) {
        new_fullfit <- do.call(update_resp, c(list(object=fullfit, newresp = y),new_args)) # may stop on error
      } else {
        new_fullfit <- try(do.call(update_resp, c(list(object=fullfit, newresp = y),new_args)))
        if (inherits(new_fullfit,"try-error")) { ## (debug.= TRUE or 1L) to return error info in parallel mode: return the try-error object
          if (debug.) {
            utils::dump.frames(dumpto="dump_on_new_fullfit", to.file=TRUE) # but doesn't stop
            return(list(full=structure(NA,new_fullfit=new_fullfit), # but the attributes seem lost at the level of the pbapply() closure.
                        null=structure(logL_new_null,new_nullfit=new_nullfit)))
          } else return(c(full=NA, null=logL_new_null)) # attributes would be lost at the level of the pbapply() closure. 
        } ## ELSE continue
      }
      logL_new_full <- logLik(new_fullfit,which=test_obj)
      #cat(" ",logL_new_full,"\n")
      conv_full <- (abs(logL_new_full - prev_logL_full)<1e-4)
      if (logL_new_full>best_logL_full) { # always true the first time
        best_fullfit <- new_fullfit
        best_logL_full <- logL_new_full
        if (logL_new_full>best_logL_null) break # the for loop
      }
    }
    if (conv_full) break # no point in refitting the null model if the new inits from full fit don't change
    # ELSE
    prev_logL_full <- logL_new_full
    newinits <- .get_outer_inits_from_fit(fitobject=best_fullfit, keep_canon_user_inits = FALSE) # for next new_nullfit
    if (null_fit_fn=="fitme") {
      new_args <- list(init=newinits, control=ctrl_opt)
    } else if (null_fit_fn=="corrHLfit") {
      new_args <- list(init.corrHLfit=newinits, control.corrHLfit=ctrl_opt)
    } else new_args <- NULL
  } # end while()
  # print(logLik(new_fullfit,which=test_obj) - logLik(new_nullfit,which=test_obj)) 
  resu <- c(full=best_logL_full,null=best_logL_null)
  return(resu)
}


# (fixme?) : create as.lm method for HLfit object?
LRT <- function(object,object2,boot.repl=0,# nb_cores=NULL, 
                resp_testfn=NULL, simuland=eval_replicate, 
#                .condition = NULL, ## bc expected by simuland, but not operational,
                ...) { 
  if (nrow(object$data)!=nrow(object2$data)) {
    stop("models were not both fitted to the same size of dataset.")
  }
  #if (length(list(...))) warning("...' arguments are currently ignored in LRT()", immediate. = TRUE) 
  #  which is a bit unfortunate ( say ...control=list(optimizer="bobyqa")) but makes parallelisation so much more straightforward...
  info <- .compare_model_structures(object,object2)
  nullfit <- info$nullfit
  fullfit <- info$fullfit
  test_obj <- info$test_obj
  LRTori <- 2*(logLik(fullfit,which=test_obj)-logLik(nullfit,which=test_obj))
  if (is.na(df <- info$df)) {
    resu <- list(nullfit=nullfit,fullfit=fullfit,basicLRT = data.frame(chi2_LR=LRTori,df=NA,p_value=NA))
  } else {
    pvalue <- 1-pchisq(LRTori,df=df) ## but not valid for testing null components of variance
    resu <- list(nullfit=nullfit,fullfit=fullfit,basicLRT = data.frame(chi2_LR=LRTori,df=df,p_value=pvalue)) ## format appropriate for more tests  
  }
  if (boot.repl) {
    if (boot.repl<100L && ! is.na(df)) message("It is recommended to set boot.repl>=100 for Bartlett correction")
#    dotargs <- match.call(expand.dots = FALSE)$... ## produce a pairlist of (essentially) promises. No quote() needed
    isdebugd <- isdebugged(simuland) # bc the assignment of environment drops this status
    environment(simuland) <- environment() # enclosing env(simuland) <- evaluation env(LRT)
    if (isdebugd) debug(simuland)

    bootblob <- spaMM_boot(object=nullfit,nsim = boot.repl,
                           simuland=simuland, 
                           #nb_cores = nb_cores,
                           resp_testfn = resp_testfn,
                           #aslistfull=aslistfull, aslistnull=aslistnull#, simbData=simbData,
                           # debug.=debug., now in the ...
                           type="marginal", # mandatory arg of spaMM_boot()
                           #, control.foreach=list(.errorhandling="pass")
                           ...
    )
    resu <- .add_boot_results(bootblob, resu, LRTori, df, test_obj)
  }
  class(resu) <- c("fixedLRT",class(resu)) 
  return(resu)
}

## anova treated as alias for LRT
anova.HLfit <- function(object, object2=NULL, ..., method="") {
  # if (method=="anova.lm" && is.null(object2)) {
  #   #identical(fullm$models[c("eta","lambda","phi")],list(eta="etaGLM",lambda="",phi="phiScal"))
  #   .anova_HLfit_lm(object, ...) ## may now handle factors but not continuosu variance => source in 'ignored' directory
  # } else 
    LRT(object,object2, ...)
}

