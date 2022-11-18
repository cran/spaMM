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
  if (length(unest)==2L) {
    warning("Models not nested (opposite nestings for fixed and random terms). No test performed.")
    return(list(fullfit=NULL,nullfit=NULL,test_obj=NULL,df=NA))
  }
  # ELSE
  if (length(unest)==0L) {
    warning(paste("The two models appear equivalent (except perhaps for residual dispersion models).\n", 
                                    "No test performed."))
    return(list(fullfit=NULL,nullfit=NULL,test_obj=NULL,df=0))
  }
  # if (length(Rnest)) return(.process_ranef_case(object, object2, nest=Rnest)) # Possibly nested models, differing at least by their random effects.
  ###############################
  # # ELSE nested fixef, identical ranefs
  # df1 <- length(X1[!is.na(fixef(object))])
  # df2 <- length(X2[!is.na(fixef(object2))])
  # if (!is.null(Rnest)) {
  #   lambda.object <- object$lambda.object
  #   if (!is.null(lambda.object)) df1 <- df1+length(unlist(lambda.object$coefficients_lambdaS, use.names=FALSE))
  #   cov.mats <- .get_compact_cov_mats(object$strucList)
  #   if (length(cov.mats)) {
  #     nrows <- unlist(lapply(cov.mats,NROW))
  #     df1 <- df1+sum(nrows*(nrows-1)/2)
  #   }
  #   lambda.object <- object2$lambda.object
  #   if (!is.null(lambda.object)) df2 <- df2+length(unlist(lambda.object$coefficients_lambdaS, use.names=FALSE))
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
          warning("LRT comparing REML fits with different fixed-effect conditions is highly suspect", 
                  immediate.=TRUE)
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

.add_famPars_outer <- function(parlist, fitobject, names_u_c_inits=NULL, type_attr) { 
  if (is.null(names_u_c_inits)) {
    canon.init <- attr(fitobject,"optimInfo")$LUarglist$canon.init ## includes user init
    names_u_c_inits <- names(unlist(canon.init))
  }
  # For each of the following family param there may already be attr(outer_ests,"type")$tr<fam par> but let's not assume that messy thing
  if ( ! is.null(families <- fitobject$families)) {
    for (mv_it in seq_along(families)) {
      fam_it <- families[[mv_it]]
      char_mv_it <- as.character(mv_it)
      if ((parname <- paste("NB_shape",mv_it,sep=".")) %in% names_u_c_inits) {
        parlist[["NB_shape"]][char_mv_it] <- environment(fam_it)$shape # <vector element> <- 
        if (type_attr) attr(parlist,"type")[["NB_shape"]][char_mv_it] <- "outer" 
      } else if ((parname <- paste("COMP_nu",mv_it,sep=".")) %in% names_u_c_inits) {
        parlist[["COMP_nu"]][char_mv_it] <- environment(fam_it)$nu # <vector element> <- 
        if (type_attr) attr(parlist,"type")[["COMP_nu"]][char_mv_it] <- "outer" 
      } else if ((parname <- paste("beta_prec",mv_it,sep=".")) %in% names_u_c_inits) {
        parlist[["beta_prec"]][char_mv_it] <- environment(fam_it)$prec # <vector element> <- 
        if (type_attr) attr(parlist,"type")[["beta_prec"]][char_mv_it] <- "outer" 
      } 
    }
  }
  ##
  if ("NB_shape" %in% names_u_c_inits) {
    parlist$NB_shape <- environment(fitobject$family$aic)$shape
    if (type_attr) attr(parlist,"type")$NB_shape <- "outer" 
  }
  if ("COMP_nu" %in% names_u_c_inits) {
    parlist$COMP_nu <- environment(fitobject$family$aic)$nu
    if (type_attr) attr(parlist,"type")$COMP_nu <- "outer" 
  }
  if ("beta_prec" %in% names_u_c_inits) {
    parlist$beta_prec <- environment(fitobject$family$aic)$prec
    if (type_attr) attr(parlist,"type")$beta_prec <- "outer" 
  }
  parlist
}

# Construct new outer inits from outer fitted values and/or initial outer values of input fit.
.get_outer_inits_from_fit <- function(fitobject, keep_canon_user_inits) {
  canon.init <- attr(fitobject,"optimInfo")$LUarglist$canon.init ## includes sanitized user init
  if (FALSE) {
    outer_ests <- get_ranPars(fitobject,lambda_names = "") ## "CorrEst_and_RanFix only" 
    #   which means that only these params were conceived to be controlled in the new outer inits.
    #    But get_ranPars(, which=NULL) is not formally defined, so has been changing... the original comment is no longer true.
    ## If the alternative not valid in the long run, this get_ranPars(.) calls should be modified.
  } else {
    optimInfo <- attr(fitobject,"optimInfo")
    outer_ests <- optimInfo$optim.pars # transparent (but potentially more comprehensive set of params)
    if ( ! is.null(outer_ests)) {
      attr(outer_ests,"moreargs") <- optimInfo$LUarglist$moreargs # necess to canonize Matern params...
      outer_ests <- .canonizeRanPars(outer_ests, corr_info=fitobject$ranef_info$sub_corr_info, 
                                     checkComplete=FALSE,  rC_transf=.spaMM.data$options$rC_transf)
    }
  }
  if (keep_canon_user_inits &&
      length(user_inits <- .reformat_corrPars(getCall(fitobject)$init,corr_families=fitobject$corr_info$corr_families))
      ) { # keep them (as interpreted in canon.init: minimum phi is 1e-4, etc) in return value
    # => remove the fitted values from the nullranPars used to modify_list
    # => keep them in 'removand' list of pars to remove from nullranPars
    # => exclude them from 'not_user_inits_names' to remove from 'removand' !
    in_init_not_user_inits <- setdiff(names(unlist(canon.init)),  names(unlist(user_inits))) # names, excluding those of parameters with user inits
    in_user_inits <- setdiff(names(unlist(outer_ests)), in_init_not_user_inits) ## it rem
    ## removand: user_inits, fixed, or inner optimized corrPars
    # locinit will retain parameters that were outer optimized without an explicit user init
    if ( is.null(in_user_inits)) {
      locinit <- outer_ests 
    } else { # replace initial value by [fitted values, except those that had a user_init]
      # => the locinit retains the sanitized user init
      not_in_user_inits <- .remove_from_cP(outer_ests, u_names=in_user_inits)
      locinit <- .modify_list(canon.init, not_in_user_inits) ## loses attributes
    }
  } else locinit <- outer_ests
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
      reinit[ ! hlfit$lambda.object$type %in% type] <- NA
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
  # Inner-estimated lambda and ranCoefs (__F I X M E__ could add phi: amusing has this was never done... inner estimated mv phi exist, incidentally)
  init.HLfit <- NULL
  rC_inner_inits <- .get_rC_inits_from_hlfit(from, type="inner") # (yes, inner, not inner_ranCoefs)
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

.anova.lm <- function (object, ...) {
  if (!inherits(object, "HLfit")) 
    warning("calling anova.lm(<fake-HLfit-object>) ...")
  w <- weights(object, type="prior")
  ssr <- sum(if (is.null(w)) residuals.HLfit(object,type="response")^2 else w * residuals.HLfit(object,type="response")^2)
  mss <- sum(if (is.null(w)) fitted(object)^2 else w * 
               fitted(object)^2)
  if (ssr < 1e-10 * mss) 
    warning("ANOVA F-tests on an essentially perfect fit are unreliable")
  dfr <- df.residual(object)
  p <- object$dfs$pforpv
  if (p > 0L) {
    p1 <- 1L:p
    qr_sXaug <- object$envir$qr_X
    comp <- crossprod(qr.Q(qr_sXaug),object$y) # object$effects[p1] # the better tibco doc says it is t(Q) %*% y 
    #  but given it is the QR for scaled X varaibles, 
    if (is.null(qr_sXaug))  stop("HLfit object does not have a 'qr' factorization of the model matrix.")
    asgn <- attr(object$X.pv,"assign")[qr_sXaug$pivot][p1]   # ____F I X M E____ use model.matrix() extractor everywhere for object$X.pv  ? 
    nmeffects <- c("(Intercept)", attr(terms(object), "term.labels"))
    tlabels <- nmeffects[1 + unique(asgn)]
    ss <- c(vapply(split(comp^2, asgn), sum, 1), ssr)
    df <- c(lengths(split(asgn, asgn)), dfr)
  }
  else {
    ss <- ssr
    df <- dfr
    tlabels <- character()
  }
  ms <- ss/df
  f <- ms/(ssr/dfr)
  P <- pf(f, df, dfr, lower.tail = FALSE)
  
  table <- data.frame(df, ss, ms, f, P)
  table[length(P), 4:5] <- NA ## row of residual SS
  dimnames(table) <- list(c(tlabels, "Residuals"), c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)"))
  #if (attr(object$terms, "intercept")) table <- table[-1, ]
  table <- table[ ! rownames(table) == "(Intercept)", ]
  structure(table, heading = c("Analysis of Variance Table\n", 
                               paste("Response:", deparse((formula(object))[[2L]]))),  ## formula ## FIXME use formula.HLfit()
            class = c("anova", "data.frame"))
  
}

.null.deviance <- function(object, intercept=attr(terms(object), "intercept")) { # see glm.fit
  if (length(offset <- model.offset(model.frame(object))) && intercept > 0L) {
    # in that case it is necessary to fit the model with intercept and offset to obtain the sensible null deviance. 
    # This is what glm() does, calling [if (length(offset) && attr(mt, "intercept") > 0L) { fit2 <- ...]
    termsv <- terms(object) # sufficient for fixed-effect model
    offset_term <- (as.character(attr(termsv, "variables"))[-1L])[[attr(termsv, "offset")]]
    newform <- as.formula(paste(" . ~ 1 +",offset_term))
    null_dev_fit <- update(object, formula.= newform )
    deviance(null_dev_fit)
  } else {
    # that's what is returned as $null.deviance by glm.fit, 
    pw <- weights(object, type="prior")
    wtdmu <- if (intercept) {
      sum(pw * object$y)/sum(pw)
    } else object$family$linkinv(model.offset(model.frame(object)))
    sum(object$family$dev.resids(object$y, wtdmu, pw))
  }
}

.df.null <- function(object, intercept=attr(terms(object), "intercept")) {
  n.ok <- length(object$y) - sum(weights(object, type="prior") == 0)
  nulldf <- n.ok - as.integer(intercept)
}

.drop.terms <- function (termobj, dropx = NULL, keep.response = FALSE) {
  if (is.null(dropx)) {
    termobj
  } else {
    if (!inherits(termobj, "terms")) 
      stop(gettextf("'termobj' must be a object of class %s", 
                    dQuote("terms")), domain = NA)
    itcp <- attr(termobj, "intercept")
    if ( ( ! length(dropped <- attr(termobj, "term.labels")[-dropx])) && 
         itcp) dropped <- '1'
    newformula <- reformulate(dropped, response = if (keep.response) termobj[[2L]], 
                              intercept = itcp, env = environment(termobj))
    result <- terms(newformula, specials = names(attr(termobj, "specials")))
    response <- attr(termobj, "response")
    dropOpt <- if (response && !keep.response) 
      c(response, dropx + length(response))
    else dropx + max(response)
    if (!is.null(predvars <- attr(termobj, "predvars"))) {
      attr(result, "predvars") <- predvars[-(dropOpt + 
                                               1)]
    }
    if (!is.null(dataClasses <- attr(termobj, "dataClasses"))) {
      attr(result, "dataClasses") <- dataClasses[-dropOpt]
    }
    result
  }
}

.anova.glm <- function(object, ..., dispersion = NULL, test = NULL) {
  doscore <- !is.null(test) && test == "Rao"
  x <- model.matrix(object)
  varseq <- attr(object$X.pv,"assign")
  nvars <- max(0, varseq)
  resdev <- resdf <- NULL
  termsv <- terms(object)
  if (doscore) {
    score <- numeric(nvars) # misnomer: regressors, not variables
    subx <- x[, varseq == 0, drop = FALSE]
    # E.g., the design matrix 'x' may have a col for intercept, two cols for a first factor, two cols for a second factor
    # varseq is 0 1 1 2 2 and there are 3 terms in the euation. dropterns removes the terms from the equation 
    #   and x[, varseq <= i, drop = FALSE] removes the corresponding blocs of cols
    # 0 always stands for the Intercept, not for the first term.
    # This first doscore fit uses the original GLM family ect to generate residuals used in the next doscore fits which are LMs 
    # x = x[, varseq == 0, drop = FALSE], y = y, weights = object$prior.weights, 
    # start = object$start, offset = object$offset, family = object$family,     
    # drop.terms fails to remove all terms because because attr(,"term.labels") does not include the intercept
    # hence calling the patch fn .drop.terms()
    newform <- .drop.terms(termsv, seq_len(nvars), keep.response = TRUE) 
    refit <- update(object, formula.= newform)
    r <- residuals(refit, type="working")
    w <- weights(refit, type="working")
    icpt <- attr(termsv, "intercept")
  }
  if (nvars > 1 || doscore) {
    method <- object$method
    for (i in seq_len(max(nvars - 1L, 0))) {
      newform <- drop.terms(termsv, (i+1L):nvars, keep.response = TRUE) 
      refit <- update(object, formula.= newform)
      if (doscore) {
        subx <- x[, varseq <= i, drop = FALSE]
        zz <- glm.fit(x=subx, y = r, weights = w, intercept = icpt)
        score[i] <- zz$null.deviance - zz$deviance
        r <- residuals(refit, type="working")
        w <- weights(refit, type="working")
      }
      resdev <- c(resdev, deviance(refit))
      resdf <- c(resdf, df.residual(refit))
    }
    if (doscore) {
      zz <- glm.fit(x=x, y = r, weights = w, intercept = icpt)
      score[nvars] <-  zz$null.deviance - zz$deviance
    }
  }
  resdf <- c(.df.null(object), resdf, df.residual(object))
  resdev <- c(.null.deviance(object), resdev, deviance(object))
  table <- data.frame(c(NA, -diff(resdf)), c(NA, pmax(0, -diff(resdev))), 
                      resdf, resdev)
  tl <- attr(termsv, "term.labels")
  if (length(tl) == 0L) 
    table <- table[1, , drop = FALSE]
  dimnames(table) <- list(c("NULL", tl), c("Df", "Deviance", 
                                           "Resid. Df", "Resid. Dev"))
  if (doscore) 
    table <- cbind(table, Rao = c(NA, score))
  varlist <- attr(termsv, "variables")
  title <- paste0("Analysis of Deviance Table", "\n\nModel: ", 
                  object$family$family, ", link: ", object$family$link, 
                  "\n\nResponse: ", as.character(varlist[-1L])[1L], "\n\nTerms added sequentially (first to last)\n\n")
  df.dispersion <- Inf
  if (is.null(dispersion)) {
    dispersion <- residVar(object,"fit")
    # for Gamma GLM anova.glm uses the MME estimate based on Pearson residuals, dispersion <- summary(object, dispersion = dispersion)$dispersion
    df.dispersion <- if (object$dfs$p_fixef_phi==0L) {
      # we actually tested phi model = phiScal in the calling function
      Inf
    } else df.residual(object)
  }
  if (!is.null(test)) {
    if (test == "F" && df.dispersion == Inf) {
      fam <- object$family$family
      if (fam == "binomial" || fam == "poisson") 
        warning(gettextf("using F test with a '%s' family is inappropriate", 
                         fam), domain = NA)
      else warning("using F test with a fixed dispersion is inappropriate")
    }
    table <- stat.anova(table = table, test = test, scale = dispersion, 
                        df.scale = df.residual(object), n = NROW(x))
  }
  structure(table, heading = title, class = c("anova", "data.frame"))
}

# LU facto with 1 on diag of L (Doolittle's scheme)
.lu_doo <- function(x, eps = 1e-06) {
  stopifnot(ncol(x) > 1L)
  n <- nrow(x)
  L <- U <- matrix(0, nrow = n, ncol = n)
  diag(L) <- rep(1, n)
  for (i in 1:n) {
    ip1 <- i + 1
    im1 <- i - 1
    for (j in 1:n) {
      U[i, j] <- x[i, j]
      if (im1 > 0) {
        for (k in 1:im1) {
          U[i, j] <- U[i, j] - L[i, k] * U[k, j]
        }
      }
    }
    if (ip1 <= n) {
      for (j in ip1:n) {
        L[j, i] <- x[j, i]
        if (im1 > 0) {
          for (k in 1:im1) {
            L[j, i] <- L[j, i] - L[j, k] * U[k, i]
          }
        }
        L[j, i] <- if (abs(U[i, i]) < eps) 
          0
        else L[j, i]/U[i, i]
      }
    }
  }
  L[abs(L) < eps] <- 0
  U[abs(U) < eps] <- 0
  list(L = L, U = U)
}

.term2colX <- function (term_names, itcp, asgn) {
  col_terms <- if (itcp) 
    c("(Intercept)", term_names)[asgn + 1L]
  else term_names[asgn[asgn > 0L]]
  nm <- union(unique(col_terms), term_names)
  res <- lapply(setNames(as.list(nm), nm), function(x) numeric(0L))
  map <- split(seq_along(col_terms), col_terms)
  res[names(map)] <- map
  res[nm]
}

.get_type1_contrasts <- function (model, 
                                  termsv=terms(model), # ____F I X M E___ rework fn for mv fits
                                  X=model.matrix(model)) {
  p <- ncol(X)
  if (p == 0L) 
    return(list(matrix(numeric(0L), nrow = 0L)))
  itcp <- attr(termsv, "intercept")
  if (p == 1L && itcp) 
    return(list(matrix(numeric(0L), ncol = 1L)))
  L <- if (p == 1L) 
    matrix(1L)
  else t(.lu_doo(crossprod(X))$L)
  dimnames(L) <- list(colnames(X), colnames(X))
  term_names <- attr(termsv, "term.labels")
  ind.list <- .term2colX(term_names=term_names, itcp=itcp, asgn= attr(X, "assign"))[term_names]
  lapply(ind.list, function(rows) L[rows, , drop = FALSE])
}

.term_contain <- function (term, factors, dataClasses, term_names) 
{
  get_vars <- function(term) rownames(factors)[factors[, term] == 1L]
  contain <- function(F1, F2) {
    all(vars[[F1]] %in% vars[[F2]]) && 
      length(setdiff(vars[[F2]], vars[[F1]])) > 0L && setequal(numerics[[F1]], numerics[[F2]])
  }
  vars <- lapply(setNames(term_names, term_names), get_vars)
  numerics <- lapply(vars, function(varnms) varnms[which(dataClasses[varnms] == "numeric")])
  sapply(term_names, function(term_nm) contain(term, term_nm))
}

.containment <- function (termsv) {
  data_classes <- attr(termsv, "dataClasses")
  term_names <- attr(termsv, "term.labels")
  factor_mat <- attr(termsv, "factors")
  lapply(setNames(term_names, term_names), function(term) {
    term_names[.term_contain(term, factor_mat, data_classes, 
                            term_names)]
  })
}

.get_type2_contrasts <- function (model, 
                                  termsv=terms(model), # ____F I X M E___ rework fn for mv fits. Set assign attribute first...
                                  X=model.matrix(model)) {
  # if (inherits(termsv,"list")) {
  #   resu <- vector("list",length(termsv))
  #   for (mv_it in seq_along(termsv)) {
  #     resu[[mv_it]] <- .get_type2_contrasts(model, termsv=termsv[[mv_it]], X=X)
  #   }
  #   return(unlist(resu,use.names = FALSE, recursive = FALSE))
  # }
  data_classes <- attr(termsv, "dataClasses")
  asgn <- attr(X, "assign")
  term_names <- attr(termsv, "term.labels")
  #
  which <- term_names
  if (ncol(X) <= 1L || length(term_names) <= 1L) 
    return(.get_type1_contrasts(model))
  #
  else stopifnot(is.character(which), all(which %in% term_names))
  which <- setNames(as.list(which), which)
  is_contained <- .containment(model)
  itcp <- attr(termsv, "intercept") > 0
  col_terms <- if (itcp) 
    c("(Intercept)", term_names)[asgn + 1]
  else term_names[asgn[asgn > 0]]
  term2colX <- split(seq_along(col_terms), col_terms)[unique(col_terms)]
  lapply(which, function(term) {
    cols_term <- unlist(term2colX[c(term, is_contained[[term]])])
    Xnew <- cbind(X[, -cols_term, drop = FALSE], X[, cols_term, drop = FALSE])
    newXcol_terms <- c(col_terms[-cols_term], col_terms[cols_term])
    Lc <- t(.lu_doo(crossprod(Xnew))$L)
    dimnames(Lc) <- list(colnames(Xnew), colnames(Xnew))
    Lc[newXcol_terms == term, colnames(X), drop = FALSE]
  })
}

.anova_fallback <- function(fitobject, type="2", rhs=NULL, test="Chisq.", ...) {
  beta <- fixef(fitobject)
  beta_cov <- vcov(fitobject)
  Llist <- switch(type,
                  "I" = .get_type1_contrasts(fitobject),
                  "1" = .get_type1_contrasts(fitobject),
                  "II" = .get_type2_contrasts(fitobject),
                  "2" = .get_type2_contrasts(fitobject),
                  "III" = stop("type-3 contrasts not implemented (and may never be)"),
                  "3" = stop("type-3 contrasts not implemented (and may never be)"),
                  stop("Unknown ANOVA 'type' speification")
  )
  # no line for Intercept: as in lmerTest output
  numtable <- matrix(NA,nrow=length(Llist), ncol = 3)
  colnames(numtable) <- c("Df", test, paste("Pr(>", test, ")", sep = ""))
  for (it in seq_along(Llist)) {
    L <- Llist[[it]]
    if (is.null(dim(L))) L <- t(L)
    if (is.null(rhs)) rhs <- rep(0, nrow(L))
    q <- NROW(L)
    Lbeta <- L %*% beta - rhs
    vcov_Lbeta <- L %*% tcrossprod(beta_cov,L)
    waldX2 <- as.vector(t(Lbeta) %*% solve(vcov_Lbeta,Lbeta))
    p <- pchisq(waldX2, q, lower.tail = FALSE)
    numtable[it, ] <- c(q, waldX2, p)
  }
  resu <- as.data.frame(numtable)
  rownames(resu) <-  names(Llist)
  attr(resu, "heading") <- paste(test, "tests for each term (type", utils::as.roman(as.integer(type)), "contrasts)")
  attr(resu, "hypotheses") <- Llist
  class(resu) <-  c("anova", "data.frame")
  resu
}



anova.HLfit <- function(object, object2=NULL, type="2", method="", ...) {
  if (is.null(object2)) {
    if (method != "t.Chisq") {
      models <- object$models[c("eta","phi")]
      if (length(models$phi)==1L && models$phi %in% c("phiScal","")) {
        if (models$eta=="etaGLM") { 
          if (object$family$family=="gaussian" && object$family$link=="identity") {
            return(.anova.lm(object, ...))
          } else return(.anova.glm(object, ...))
        } else if (object$family$family=="gaussian" && object$family$link=="identity") { # LMM
          if (requireNamespace("lmerTest",quietly=TRUE)) {
            lmlt <-  as_LMLT(object, ...)
            return(anova(lmlt, type=type)) # other possible arguments not currently meaningful
          } else if ( ! identical(spaMM.getOption("lmerTest_warned"),TRUE)) {
            message("If the lmerTest package were installed, a traditional anova table could be computed.")
            .spaMM.data$options$lmerTest_warned <- TRUE
          } 
        } 
      }
    }
    return(.anova_fallback(fitobject=object, type=type, ...)) # dots ma contain 'rhs' and (later) 'test'
  } else {
    ## anova treated as alias for LRT()
    LRT(object,object2, ...)
  }
}

# anova is an S3 class, LMLT an S4 class
## Aim of this method is to ensure that the lmerTest package is loaded when anova() is called
# bc it may be called when object has class LMLT, but no def may yet be accessible for such objects.
## See the 'ZZZ' version for detailed explanations of the code.
anova.LMLT <- function(object, ...) { 
  if (is.null(getClassDef("LMLT", where = .spaMM.data$class_cache, inherits = FALSE))) {
    # is as_LMLT has not been called since restarting R session-> lmerTest presumably not loaded
    if (requireNamespace("lmerTest",quietly=TRUE)) {
      # Hack to define object not of class LMLT (-> infinite recursion) but with similar "contains", using only default coerce methods
      setClass(Class="LMLT", contains = c("LMLTslots","lmerModLmerTest"), where=.spaMM.data$class_cache) 
      setClass(Class="LMLTinternal", contains = c("LMLTslots","lmerModLmerTest"), where=.spaMM.data$class_cache) 
    } else message("If the lmerTest package were installed, a traditional anova table could be computed.")
  } 
  object <- as(as(object,"LMLTslots"),"LMLTinternal")
  anova(object, ...)
}

