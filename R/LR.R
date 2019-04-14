.is_2_in_1 <- function(X1,X2, tol=1e-10) {
  qrX <- qr(X1)
  proj2in1 <- X1 %*% qr.solve(qrX,X2)
  return(all(abs(proj2in1 -X2)<tol))
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
    stop("Models may not be nested (distinct families)") ## but COMpoisson vs poisson ?
  }
  if (! identical(meth1,meth2) || length(REML)>1 ) {
    stop("object fitted by different methods cannot be compared")
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
  if (length(dR12) && length(dR21)) { 
    stop(paste("Fixed-effects specifications from both models seem equivalent,\n",
               "andd random-effect specifications may be non-nested\n", 
               "(a better algorithm would be required to check this properly).\n",
               "This case is not handled.")) 
  } else if (length(dR12)) {
    Rnest <- "2in1"
    warning(paste("Random-effect specifications appear distinct.\n", 
                  "This procedure is unreliable for comparing models with different random effects."),
            immediate.=TRUE)
  } else if (length(dR21)) {
    Rnest <- "1in2"
    warning(paste("Random-effect specifications appear distinct.\n", 
                  "This procedure is unreliable for comparing models with different random effects."),
            immediate.=TRUE)
  } else {
    Rnest <- NULL
    if (is.null(Xnest)) stop(paste("The two models appear equivalent (except perhaps for residual dispersion models).\n", 
                                   "This case is not handled."))
  }
  nest <- c(Xnest,Rnest)
  unest <- unique(nest)
  if (length(unest)==2L) {
    stop("Models not nested (opposite nestings for fixed and random terms). ")
  } else {
    df1 <- length(X1[!is.na(fixef(object))])
    df2 <- length(X2[!is.na(fixef(object2))])
    if (!is.null(Rnest)) {
      lambda.object <- object$lambda.object
      if (!is.null(lambda.object)) df1 <- df1+length(unlist(lambda.object$coefficients_lambdaS))
      cov.mats <- .get_compact_cov_mats(object$strucList,later=FALSE)
      if (length(cov.mats)) {
        nrows <- unlist(lapply(cov.mats,NROW))
        df1 <- df1+sum(nrows*(nrows-1)/2)
      }
      lambda.object <- object2$lambda.object
      if (!is.null(lambda.object)) df2 <- df2+length(unlist(lambda.object$coefficients_lambdaS))
      cov.mats <- .get_compact_cov_mats(object2$strucList,later=FALSE)
      if ( length(cov.mats)) {
        nrows <- unlist(lapply(cov.mats,NROW))
        df2 <- df2+sum(nrows*(nrows-1)/2)
      }
    }
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
      message("Nested models differing both by in their fixed and in their random terms. ")
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
            warning("LRT comparing REML fits with different designs is highly suspect")
          }
        }
        testlik <- "p_v"
      } else { ## random effect test
        if ( ! REML) warning("ML fits used to compare different random-effects models...")
        testlik <- "p_bv" ## used in both case, identical to p_v in the non-standard case
        stop("The two models have identical fixed-effect formulas\n and cannot yet be compared properly by this function.")
        ## need to take into account correlations in random slope models for example
      }
    } 
  }
  return(list(fullfit=fullm,nullfit=nullm,test_obj=testlik,df=df))
}

.get_inits_from_fit <- function(fitobject, keep_canon_user_inits) {
  canon.init <- attr(fitobject,"optimInfo")$LUarglist$canon.init ## includes user init
  #
  nullranPars <- get_ranPars(fitobject)
  names_u_nullranPars <- names(unlist(nullranPars))
  names_u_c_inits <- names(unlist(canon.init))
  if (keep_canon_user_inits) { # keep them (as interpreted in canon.init: minimum phi is 1e-4, etc) in return value
    # => remove the fitted values from the nullranPars used to modify_list
    # => keep them in 'removand' list of pars to remove from nullranPars
    # => exclude them from 'not_user_inits_names' to remove from 'removand' !
    user_inits <- .post_process_parlist(getCall(fitobject)$init,corr_families=fitobject$corr_info$corr_families)
    names_u_u_inits <- names(unlist(user_inits))
    not_user_inits_names <- setdiff(names_u_c_inits, names_u_u_inits) # names, excluding those of parameters with user inits
    removand <- setdiff(names_u_nullranPars, not_user_inits_names) ## removand: user_inits, fixed, or inner optimized corrPars
    ## removand: user_inits, fixed, or inner optimized corrPars
    # locinit will retain parameters that were outer optimized without an explicit user init
  } else removand <- setdiff(names_u_nullranPars, names_u_c_inits)
  if ( is.null(removand)) {
    locinit <- .modify_list(canon.init,nullranPars)
  } else { ## leaves user_inits as there are in LUarglist$canon.init, and do not add any fixed or inner-optimized par
    locinit <- .modify_list(canon.init,
                            .remove_from_cP(nullranPars,u_names=removand)) ## loses attributes
  }
  return(locinit)
}

.eval_replicate <- function(y) { # no additional arguments, to ease parallel programming => next lines instead
  # the function will be called within e.g. pbapply so it's useless to refer to parent.frame() here
  enclosing_env <- parent.env(environment()) ## this is not necess for the code to run, but for the CRAN checks not to complain
  nullfit <- get("nullfit", enclosing_env)
  fullfit <- get("fullfit", enclosing_env)
  #dotargs <- get("dotargs", enclosing_env)
  test_obj <- get("test_obj", enclosing_env)
  debug. <- get("debug.", enclosing_env)
  #  .condition <- get(".condition", enclosing_env)
  fittingFunction <- paste(getCall(nullfit)[[1]])
  newinits <- .get_inits_from_fit(fitobject=nullfit, 
                                  keep_canon_user_inits = FALSE)                                     
  if (fittingFunction=="fitme") {
    new_args <- list(init=newinits)
  } else if (fittingFunction=="corrHLfit") {
    new_args <- list(init.corrHLfit=newinits)
  } else new_args <- NULL
  if (debug.==2) {
    re_nullfit <- do.call(update_resp, c(list(object=nullfit, newresp = y),new_args)) # may stop on error
  } else {
    re_nullfit <- try(do.call(update_resp, c(list(object=nullfit, newresp = y),new_args)))
    if (inherits(re_nullfit,"try-error")) {
      if (debug.) { ## (debug.= TRUE or 1L) to return error info in parallel mode: return the try-error object
        return(c(NA,re_nullfit))
      } else return(c(NA,NA))
    }
    ## ELSE return pair of likelihoods
  }
  # Allows the user to control the starting values of the re_fullfit:
  # edotargs <- dotargs
  # for (st in setdiff(names(edotargs),"prior.weights")) {
  #   edotargs[[st]] <- eval(edotargs[[st]],env=environment()) ## evaluate the promises in the current execution envir
  # }
  # #
  newinits <- .get_inits_from_fit(fitobject=re_nullfit, keep_canon_user_inits = FALSE)
  if (fittingFunction=="fitme") {
    new_args <- list(init=newinits)
  } else if (fittingFunction=="corrHLfit") {
    new_args <- list(init.corrHLfit=newinits)
  } else new_args <- NULL
  if (debug.==2) {
    re_fullfit <- do.call(update_resp, c(list(object=fullfit, newresp = y),new_args)) # may stop on error
  } else {
    re_fullfit <- try(do.call(update_resp, c(list(object=fullfit, newresp = y),new_args)))
    if (inherits(fullfit,"try-error")) {
      if (debug.) { 
        return(c(re_fullfit,re_nullfit))
      } else return(c(NA,re_nullfit))
    }
    ## ELSE:
  }
  LRstat <- 2*(logLik(re_fullfit,which=test_obj)-logLik(re_nullfit,which=test_obj))  
  if (1-pchisq(LRstat,df=sum(re_fullfit$dfs)-sum(re_nullfit$dfs))<0) { # ie, if (FALSE)
    newinits <- .get_inits_from_fit(fitobject=re_fullfit, keep_canon_user_inits = FALSE)
    if (fittingFunction=="fitme") {
      new_args <- list(init=newinits)
    } else if (fittingFunction=="corrHLfit") {
      new_args <- list(init.corrHLfit=newinits)
    }
    if (debug.==2) {
      new_nullfit <- do.call(update_resp, c(list(object=nullfit, newresp = y),new_args)) # may stop on error
    } else {
      new_nullfit <- try(do.call(update_resp, c(list(object=nullfit, newresp = y),new_args)))
      if (inherits(new_nullfit,"try-error")) {
        if (debug.) { ## (debug.= TRUE or 1L) to return error info in parallel mode: return the try-error object
          return(c(NA,new_nullfit))
        } else return(c(NA,NA))
      }
    }
    logL_new_null <- logLik(new_nullfit,which=test_obj)
    if (logLik(new_nullfit,which=test_obj)>logLik(re_nullfit,which=test_obj)) re_nullfit <- new_nullfit
  }
  resu <- c(full=logLik(re_fullfit,which=test_obj),null=logLik(re_nullfit,which=test_obj))
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
  fittingFunction <- paste(getCall(nullfit)[[1]])
  conv_full <- conv_null <- FALSE
  best_logL_full <- best_logL_null <- prev_logL_full <- prev_logL_null <- -Inf
  # Allows the user to control the starting values of the initial new_nullfit
  # edotargs <- dotargs
  # for (st in setdiff(names(edotargs),"prior.weights")) {
  #   edotargs[[st]] <- eval(edotargs[[st]],env=environment()) ## evaluate the promises in the current execution envir
  # }
  newinits <- .get_inits_from_fit(fitobject=nullfit, 
                                  keep_canon_user_inits = FALSE)                 
  if (fittingFunction=="fitme") {
    new_args <- list(init=newinits)
  } else if (fittingFunction=="corrHLfit") {
    new_args <- list(init.corrHLfit=newinits)
  } else new_args <- NULL
  while ( TRUE ) {
    if (debug.==2) {
      new_nullfit <- do.call(update_resp, c(list(object=nullfit, newresp = y),new_args)) # may stop on error
    } else {
      new_nullfit <- try(do.call(update_resp, c(list(object=nullfit, newresp = y),new_args)))
      if (inherits(new_nullfit,"try-error")) {
        if (debug.) { ## (debug.= TRUE or 1L) to return error info in parallel mode: return the try-error object
          return(c(NA,new_nullfit))
        } else return(c(NA,NA))
      }
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
    newinits <- .get_inits_from_fit(fitobject=best_nullfit, keep_canon_user_inits = FALSE)
    if (fittingFunction=="fitme") {
      new_args <- list(init=newinits)
    } else if (fittingFunction=="corrHLfit") {
      new_args <- list(init.corrHLfit=newinits)
    }
    if (debug.==2) {
      new_fullfit <- do.call(update_resp, c(list(object=fullfit, newresp = y),new_args)) # may stop on error
    } else {
      new_fullfit <- try(do.call(update_resp, c(list(object=fullfit, newresp = y),new_args)))
      if (inherits(new_fullfit,"try-error")) {
        if (debug.) { 
          return(c(new_fullfit,new_nullfit))
        } else return(c(NA,new_nullfit))
      }
    }
    logL_new_full <- logLik(new_fullfit,which=test_obj)
    #cat(" ",logL_new_full,"\n")
    conv_full <- (abs(logL_new_full - prev_logL_full)<1e-4)
    if (logL_new_full>best_logL_full) { # always true the first time
      best_fullfit <- new_fullfit
      best_logL_full <- logL_new_full
    }
    if (conv_full) break # no point in refitting the null model if the new inits from full fit don't change
    # ELSE
    prev_logL_full <- logL_new_full
    newinits <- .get_inits_from_fit(fitobject=best_fullfit, keep_canon_user_inits = FALSE)
    if (fittingFunction=="fitme") {
      new_args <- list(init=newinits)
    } else if (fittingFunction=="corrHLfit") {
      new_args <- list(init.corrHLfit=newinits)
    }
  } # end while()
  # print(logLik(new_fullfit,which=test_obj) - logLik(new_nullfit,which=test_obj)) 
  resu <- c(full=best_logL_full,null=best_logL_null)
  # if ( ! is.null(.condition)) {
  #   condition <- with(list(nullfit=best_nullfit, fullfit=best_fullfit), eval(.condition))
  #   resu <- c(resu, condition=condition)
  # }
  return(resu)
}


# (fixme?) : create as.lm method for HLfit object?
LRT <- function(object,object2,boot.repl=0,nb_cores=NULL, boot_fn="spaMM_boot", 
                resp_testfn=NULL, debug.=FALSE, type="marginal", simuland=.eval_replicate, 
#                .condition = NULL, ## bc expected by simuland, but not operational,
                ...) { 
  if (nrow(object$data)!=nrow(object2$data)) {
    stop("models were not both fitted to the same size of dataset.")
  }
  if (length(list(...))) warning("...' arguments are currently ignored in LRT()", immediate. = TRUE) 
  #  which is a bit unfortunate ( say ...control=list(optimizer="bobyqa")) but makes parallelisation so much more straightforward...
  info <- .compare_model_structures(object,object2)
  nullfit <- info$nullfit
  fullfit <- info$fullfit
  test_obj <- info$test_obj
  df <- info$df
  LRTori <- 2*(logLik(fullfit,which=test_obj)-logLik(nullfit,which=test_obj))
  pvalue <- 1-pchisq(LRTori,df=df) ## but not valid for testing null components of variance
  resu <- list(nullfit=nullfit,fullfit=fullfit,basicLRT = data.frame(chi2_LR=LRTori,df=df,p_value=pvalue)) ## format appropriate for more tests  
  if (boot.repl) {
    if (boot.repl<100L) message("It is recommended to set boot.repl>=100 for Bartlett correction")
    nb_cores <- .check_nb_cores(nb_cores=nb_cores)
#    dotargs <- match.call(expand.dots = FALSE)$... ## produce a pairlist of (essentially) promises. No quote() needed
    isdebugd <- isdebugged(simuland) # bc the assignment of environment drops this status
    environment(simuland) <- environment() # enclosing env(simuland) <- evaluation env(LRT)
    if (isdebugd) debug(simuland)

    bootblob <- spaMM_boot(object=nullfit,nsim = boot.repl,
                           simuland=simuland, 
                           nb_cores = nb_cores,
                           resp_testfn = resp_testfn,
                           #aslistfull=aslistfull, aslistnull=aslistnull#, simbData=simbData,
                           debug.=debug., type=type
                           #, control.foreach=list(.errorhandling="pass")
    )
    
    bootreps <- bootblob$bootreps
    if (is.matrix(bootreps)) {
      colnames(bootreps)[1:2] <- paste0(c("full.","null."),test_obj) # which may already be the case
      bootdL <- bootreps[,1]-bootreps[,2]
      meanbootLRT <- 2*mean(bootdL)
      resu <- c(resu,list(rawBootLRT = data.frame(chi2_LR=LRTori,df=df,p_value=(1+sum(bootdL>=LRTori/2))/(boot.repl+1)))) ## format appropriate for more tests  
      LRTcorr <- LRTori*df/meanbootLRT
      resu <- c(resu,list(BartBootLRT = data.frame(chi2_LR=LRTcorr,df=df,p_value=1-pchisq(LRTcorr,df=df)))) ## format appropriate for more tests  
      bootInfo <- list(meanbootLRT = meanbootLRT,bootreps = bootreps, RNGstates=bootblob$RNGstates)
      resu <- c(resu,list(bootInfo=bootInfo)) ## keeps the sublist structure, which is not compatible with hglmjob.R...  
    } else {## a list in debug. case; following code to avoid bug
      bootInfo <- list(meanbootLRT = data.frame(chi2_LR=NA,df=df,p_value=NA),
                       bootreps = bootreps, RNGstates=bootblob$RNGstates)
      resu <- c(resu,list(bootInfo=bootInfo,
                          rawBootLRT = data.frame(chi2_LR=LRTori,df=df,p_value=NA))) 
    }
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

