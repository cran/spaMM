.check_nb_cores <- function(nb_cores=NULL) {
  if (is.null(nb_cores)) nb_cores <- spaMM.getOption("nb_cores") ## may be NULL
  machine_cores <- parallel::detectCores()
  if (is.null(nb_cores)) {
    nb_cores <- 1L ## default
    if (machine_cores>1L && interactive()) {
      if (! identical(spaMM.getOption("cores_avail_warned"),TRUE)) {
        message(paste(machine_cores,
                      "cores are available for parallel computation\n(you may be allowed to fewer of them on a shared cluster).\nChange 'nb_cores' argument to use some of them.\nUse spaMM.options(nb_cores=<n>) to control nb_cores globally.\n"))
        .spaMM.data$options$cores_avail_warned <- TRUE
      } else if (nb_cores>machine_cores) {
        if (! identical(spaMM.getOption("nb_cores_warned"),TRUE)) {
          warning("More cores are requested than found by parallel::detecCores(). Check 'nb_cores' argument.")
          ## + reduce it ?
          .spaMM.data$options$nb_cores_warned <- TRUE
        }
      }
    }
  }
  return(nb_cores)
}

.check_binomial_formula <- function(nullfit, fullfit, data) { ## has become obsolete as update_resp can use original names
  res <- list()
  nform <- formula.HLfit(nullfit, which="hyper")
  if (paste(nform[[2L]])[[1L]]=="cbind") {
    ## We have different possible (exprL,exprR) arguments in cbind(exprL,exprR), 
    ## but in all case the predictor is that of exprL and exprR is $BinomialDen - exprL. We standardize: 
    res$nposname <- nposname <- .makenewname("npos",names(data))
    res$nnegname <- nnegname <- .makenewname("nneg",names(data))
    nform <- paste(nform)
    nform[2L] <- paste0("cbind(",nposname,",",nnegname,")")
    res$null_formula <- as.formula(paste(nform[c(2,1,3)],collapse=""))
    if ( ! is.null(fullfit)) {
      fform <- paste(formula.HLfit(fullfit, which="hyper"))
      fform[2L] <- nform[2L]
      res$full_formula <- as.formula(paste(fform[c(2,1,3)],collapse=""))
    }
    res$cbindTest <- TRUE
  } else res$cbindTest <- FALSE
  return(res)
}

.preprocess_data <- function(data, null.formula, formula, resid.formula=NULL, prior.weights, ...) {
  mc <- match.call()
  if ( ! "prior.weights" %in% names(mc)) mc["prior.weights"] <- list(NULL) # don't test the value of prior.weights!
  mc[[1L]] <- get(".GetValidData_info", asNamespace("spaMM"), inherits=FALSE) 
  if ( inherits(data,"list")) {
    for (lit in seq_along(data)) {
      mc[["data"]] <- data[[lit]]
      mc[["formula"]] <- null.formula[-2L]
      null.rownames <- eval(mc)$rownames ## will remove rows with NA's in required variables
      mc[["formula"]] <- formula[-2L]
      full.rownames <- eval(mc)$rownames
      data[[lit]] <- data[[lit]][intersect(null.rownames,full.rownames),,drop=FALSE]
    }
  } else {
    mc[["formula"]] <- null.formula[-2L]
    null.rownames <- eval(mc)$rownames ## will remove rows with NA's in required variables
    mc[["formula"]] <- formula[-2L]
    full.rownames <- eval(mc)$rownames
    data <- data[intersect(null.rownames,full.rownames),,drop=FALSE]
  }  
  return(data)
}

.eval_both_fits <- function(nullm_call, fullm_call, fittingFunction) {
  nullfit <- eval(nullm_call)
  locinit <- .get_outer_inits_from_fit(nullfit, keep_canon_user_inits=TRUE) # to initiate fullfit (no previous fullfit available)
  if (fittingFunction=="fitme") {
    fullm_call$init <- locinit
  } else if (fittingFunction=="corrHLfit") {
    fullm_call[["init.corrHLfit"]] <- locinit
  }
  fullfit <- eval(fullm_call)
  if (logLik(fullfit)<logLik(nullfit)) { ## evidence of fullfit being trapped in a local maximum (or so)
    ## We did not overwrite user inits: we do so (but perhaps not optimal if fitted phi at 1e-6
    re_locinit <- get_inits_from_fit(from=nullfit, template=fullfit, to_fn=fittingFunction )# to initiate next fullfit
    if ( ! identical(locinit,re_locinit)) {
      if (fittingFunction=="fitme") {
        fullm_call$init <- locinit$init
      } else if (fittingFunction=="corrHLfit") {
        fullm_call[["init.corrHLfit"]] <- locinit[["init.corrHLfit"]]
      }
      if ( ! is.null(locinit$init.HLfit)) fullm_call[["init.HLfit"]] <- locinit[["init.HLfit"]]
      fullfit <- eval(fullm_call)
    } ## else it's not clear what to do, short of reproducing the eval_replicate2() concept
  }
  # # No comparable evidence that nullfit is trapped in a local maximum: check with fullfit ranPars
  locinit <- .get_outer_inits_from_fit(nullfit, keep_canon_user_inits=FALSE) # for next renullfit
  if (fittingFunction=="fitme") {
    nullm_call$init <- locinit
  } else if (fittingFunction=="corrHLfit") {
    nullm_call[["init.corrHLfit"]] <- locinit ## but currently not used.
  }
  if (fittingFunction=="fitme") {
    renullfit <- eval(nullm_call)
    if (logLik(nullfit)<logLik(renullfit)) nullfit <- renullfit ## test may be FALSE ./. 
    # ./. (typically if fullfit yields a low lambda which is a local maximum of the nullfit)
  } ## ELSE currently not this recheck for corrHLfit
  return(list(nullfit=nullfit, fullfit=fullfit))
}


.add_boot_results <- function(bootblob, resu, LRTori, df, test_obj) {
  bootreps <- bootblob$bootreps
  colnames(bootreps)[1:2] <- paste0(c("full.","null."),test_obj) # which may already be the case
  if (is.matrix(bootreps)) {
    bootreps <- data.frame(bootreps)
    if (anyNA(bootreps)) {
      bootreps <- na.omit(bootreps)
      n_omitted <- length(attr(na.omit(bootreps),"na.action"))
      warnmess <- paste0(n_omitted," bootstrap replicate(s) apparently failed and are omitted for p-value computations.")
      #warning(warnmess, immediate.=TRUE)
      bootblob$warnlist$n_omitted <- warnmess
    } # but thse are retained in the bootblob that is returned
    bootdL <- bootreps[,1L]-bootreps[,2L]
    # if ( ! is.null( .condition )) {
    #   lrfit <- lm(bootdL ~ condition ,data=bootreps)
    #   meanbootLRT <- 2*predict(lrfit,newdata=data.frame(condition=eval(.condition))) ## conditional mean
    #   attr(meanbootLRT,"boot_type") <- "conditional"
    # } else {
    if (any(bootdL < -2e-04)) {
      neg_values <- bootdL[bootdL<0]
      diagn_test <- stats::ks.test(x= -neg_values,y=.neg_r_chinorm, alternative = "less")
      ##### .neg_r_chinorm was produced as follows for a quick and dirty test (null hypo: chi2+gaussian noise(SD=1e-4))
      # set.seed(123)
      # r_chinorm <- rnorm(1e6, mean=rchisq(n=1e6, df=1), sd = 1e-4)
      # .neg_r_chinorm <- -r_chinorm[r_chinorm<0]
      # save(.neg_r_chinorm, file="C:\\home\\francois\\travail\\stats\\spaMMplus\\spaMM\\package/R/sysdata.rda", 
      #    compress="bzip2", version = 2) # version control compared to:
      ## devtools::use_data(.neg_r_chinorm, internal=TRUE) # while directory is set to .../package/ 
      ##### Ideally we would like to have a one-sample test instead.
      # Also, The test ignores the frequency of negatives (length(.neg_r_chinorm)=3284).
      if (diagn_test$p.value<0.01) {
        warnmess <- paste0("Suspiciously large negative values in bootstrap distribution of likelihood ratio\n",
                           "  (up to ",signif(min(neg_values),3L),"). These are treated as 0.")
        #warning(warnmess, immediate.=TRUE)
        bootblob$warnlist$neg_values <- warnmess
      }
    }
    bootdL <- pmax(0,bootdL)
    meanbootLRT <- 2*mean(bootdL)
    attr(meanbootLRT,"boot_type") <- "marginal"
    # }
    LRTcorr <- LRTori*df/meanbootLRT
    resu$BartBootLRT <- data.frame(chi2_LR=LRTcorr,df=df,p_value=1-pchisq(LRTcorr,df=df))
    rawPvalue <- (1+sum(bootdL>=LRTori/2))/(nrow(bootreps)+1) ## DavisonH, p.141
    ## as documented in ?LRT
    resu$rawBootLRT <- data.frame(chi2_LR=LRTori,df=df,p_value=rawPvalue)
    resu$bootInfo <- c(bootblob,list(meanbootLRT=meanbootLRT)) 
  } else {
    warning(" spaMM_boot is not returning a *matrix* of logLik values, which points to a bug (or an explicit debug.)", 
            immediate.=TRUE)
    # trying to return the info without stopping:
    bootInfo <- c(bootblob,list(meanbootLRT=data.frame(chi2_LR=NA,df=df,p_value=NA)))
    resu <- c(resu,list(bootInfo=bootInfo,
                        rawBootLRT = data.frame(chi2_LR=LRTori,df=df,p_value=NA))) 
  }
  return(resu)
}



# fixedLRT -> .LRT -> spaMM_boot
.LRT <- function(null.formula=NULL,formula,
                 null.disp=list(), boot.repl=0,
                 ## currently trace always false; this is not an argument t be forwarded as is to corrHLfit! 
                 #trace=FALSE, ## T means lead to calls of corrHLfit(... trace=list(<file name>,<over/append>))
                 verbose=c(TRACE=FALSE),
                 fittingFunction="corrHLfit",  
                 simuland= eval_replicate,
                 data,
                 resp_testfn=NULL,
                 #                 .condition = NULL, ## only an argument of the internal function .LRT() so not visible at user level.
                 seed=NULL,
                 nb_cores=NULL, 
                 debug.=FALSE, 
                 #type="marginal", # explicit in call to spaMM_boot() below.
                 ... # I cannot use the dots to pass arguments to spaMM_boot bc they also contain arguments for the fitting functions
                 ) {
  # if (is.na(verbose["trace"])) verbose["inner"] <- FALSE
  if (is.na(verbose["all_objfn_calls"])) verbose["all_objfn_calls"] <- FALSE   
  oricall <- match.call(expand.dots = TRUE)
  #
  boot_call <- oricall[which( ! names(oricall) %in% names(formals(.LRT)))] ## includes the position of the called fn in oricall
  ## birth pangs :
  if ("predictor" %in% names(oricall)) {
    stop(".LRT() called with 'predictor' argument which should be 'formula'" )
  }
  if ("null.predictor" %in% names(oricall)) {
    stop(".LRT() called with 'null.predictor' argument which should be 'null.formula'" )
  }
  ## here we makes sure that *predictor variables* are available for all data to be used under both models
  resid.model <- .reformat_resid_model(oricall$resid.model) ## family not easily checked here; not important
  loccall <- oricall
  loccall[[1L]] <- get(".preprocess_data", asNamespace("spaMM"), inherits=FALSE)
  loccall$"resid.formula" <- resid.model$formula
  data <- eval(loccall,parent.frame()) 
  #
  boot_call$data <- data
  boot_call$verbose <- verbose 
  if (fittingFunction == "fitme") {
    if (is.null(boot_call$method)) stop("'method' argument is required when fittingFunction=\"fitme\".")
    locmethod <- boot_call$method
  } else {
    locmethod <- boot_call$HLmethod
  }
  fullm_call <- boot_call
  fullm_call[[1L]] <- as.name(fittingFunction) ## so that I can paste() it in eval_replicate
  nullm_call <- fullm_call
  fullm_call$formula <- formula
  
  #### "limited use, for initializing bootstrap replicates:"
  if ( ! is.null(null.formula)) { ## ie if test effet fixe
    testFix <- TRUE
    if (locmethod =="SEM") {
      test_obj <- "logLapp"
    } else test_obj <- "p_v"
    ## check fullm.list$REMLformula, which will be copied into nullm in all cases of fixed LRTs
    if (locmethod %in% c("ML","PQL/L","SEM") || substr(locmethod,0,2) == "ML") {
      fullm_call$REMLformula <- NULL
      nullm_call$REMLformula <- NULL
    } else { ## locmthod secifies standard or non-standard REML. 
      ## Two attempts to get an LRT from REML; in each case the same conditioning is used in both fits
      if (is.null(fullm_call$REMLformula)) { ## default call
        # keep fullm_call$REMformula NULL => standard REML with conditioning defined by fullm_call$formula with be run
        nullm_call$REMLformula <- fullm_call$formula # non-standard REML fit, ame conditioning as in fullm
      } else { ## non-default REMLformula => it's already integrated in the two calls
        # => same conditioning for both fits (nothing to do)
      }
      namesinit <- names(fullm_call$init.corrHLfit)
      namesinit <- setdiff(namesinit,c("rho","nu","Nugget","ARphi"))
      len <- length(namesinit)
      if ( len) {
        if (len > 1L) namesinit <- paste(c(paste(namesinit[-len],collapse=", "),namesinit[len]),collapse=" and ")
        message("Argument 'init.corrHLfit' is used in such a way that")
        message(paste0("  ",namesinit," will be estimated by maximization of p_v."))
        message("  'REMLformula' will be inoperative if all dispersion")
        message("  and correlation parameters are estimated in this way.")
      }
    }
    nullm_call$formula <- null.formula
  } else if ( length(null.disp)>0 ) { ## test disp/corr param
    testFix <- FALSE
    #
    stop("Models differ in their random effects or residual dispersion structure.")
    #
    namescheck <- names(fullm_call$lower)
    namescheck <- namescheck[ ! namescheck %in% names(null.disp)]  
    nullm_call$lower <- nullm_call$lower[namescheck]
    nullm_call$upper <- nullm_call$upper[namescheck]
    nullm_call$ranFix <- c(nullm_call$ranFix,null.disp) ## adds fixed values to any preexisting one
  } else testFix <- NA
  
  blob <- .eval_both_fits(nullm_call, fullm_call, fittingFunction)
  nullfit <- blob$nullfit
  fullfit <- blob$fullfit
  if (testFix) {df <- fullfit$dfs[["pforpv"]]-nullfit$dfs[["pforpv"]]} else {df <- length(null.disp)} 
  if (df < 0) {
    tmp <- fullfit
    fullfit <- nullfit
    nullfit <- tmp
    df <- - df
  }
  if (inherits(fullfit,"HLfitlist")) {
    fullL <- attr(fullfit,"APHLs")[[test_obj]]
  } else fullL <- fullfit$APHLs[[test_obj]]
  if (inherits(nullfit,"HLfitlist")) {
    nullL <- attr(nullfit,"APHLs")[[test_obj]]
  } else nullL <- nullfit$APHLs[[test_obj]]
  LRTori <- 2*(fullL-nullL)
  ## BOOTSTRAP
  if ( ! is.na(testFix)) {
    resu <- list(fullfit=fullfit,nullfit=nullfit)
    resu$basicLRT <- data.frame(chi2_LR=LRTori,df=df,p_value=1-pchisq(LRTori,df=df))
    if (boot.repl>0L) {
      environment(simuland) <- environment() # enclosing env(simuland) <- evaluation env(LRT)
      bootblob <- spaMM_boot(object=nullfit, nsim = boot.repl, simuland=simuland, 
                             nb_cores = nb_cores, #in the dots
                             resp_testfn = resp_testfn, ## for simulate
                             debug.=debug., #in the dots
                             type="marginal", # mandatory arg of spaMM_boot()
                             seed=seed
      )
      resu <- .add_boot_results(bootblob, resu, LRTori, df, test_obj)
    } ## end bootstrap
  } else { ## nothing operativ yet
    warning("code missing here")
    #bootreps <- matrix(NA,nrow=boot.repl,ncol=length(unlist(fullfit$APHLs))) 
    ## more needed here ?
    resu <- list(fullfit=fullfit)
    resu$bootInfo <- list()
    resu$basicLRT <- list()
  }
  class(resu) <- c("fixedLRT",class(resu)) 
  return(resu)
}
