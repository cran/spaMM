.check_nb_cores <- function(nb_cores=NULL) {
  if (is.null(nb_cores)) nb_cores <- spaMM.getOption("nb_cores") ## may be NULL
  machine_cores <- parallel::detectCores()
  if (is.null(nb_cores)) {
    nb_cores <- 1L ## default
    if (machine_cores>1L && interactive()) {
      if (! identical(spaMM.getOption("cores_avail_warned"),TRUE)) {
        message(paste(machine_cores,
                      "cores are available for parallel computation\n(you may be allowed to fewer of them on a shared cluster).\nChange 'nb_cores' argument to use some of them.\nUse spaMM.options(nb_cores=<n>) to control nb_cores globally.\n"))
        spaMM.options(cores_avail_warned=TRUE)
      } else if (nb_cores>machine_cores) {
        if (! identical(spaMM.getOption("nb_cores_warned"),TRUE)) {
          warning("More cores are requested than found by parallel::detecCores(). Check 'nb_cores' argument.")
          ## + reduce it ?
          spaMM.options(nb_cores_warned=TRUE)
        }
      }
    }
  }
  return(nb_cores)
}

.check_binomial_formula <- function(nullfit, fullfit, data) {
  res <- list()
  nform <- formula.HLfit(nullfit)
  if (paste(nform[[2L]])[[1L]]=="cbind") {
    ## We have different possible (exprL,exprR) arguments in cbind(exprL,exprR), 
    ## but in all case the predictor is that of exprL and exprR is $BinomialDen - exprL. We standardize: 
    res$nposname <- nposname <- .makenewname("npos",names(data))
    res$nnegname <- nnegname <- .makenewname("nneg",names(data))
    nform <- paste(nform)
    nform[2L] <- paste0("cbind(",nposname,",",nnegname,")")
    res$null_formula <- as.formula(paste(nform[c(2,1,3)],collapse=""))
    fform <- paste(formula.HLfit(fullfit))
    fform[2L] <- nform[2L]
    res$full_formula <- as.formula(paste(fform[c(2,1,3)],collapse=""))
    res$cbindTest <- TRUE
  } else res$cbindTest <- FALSE
  return(res)
}

.eval_boot_replicates <- local({
  doSNOW_warned <- FALSE
  function(eval_replicate, ## function to run on each replicate
                                    boot.repl, ## number of bootstrap replicates
                                    nullfit, ## the fitted model object for which bootstrap replicates are drawn 
                                    nb_cores, ## passing explicit value from user
                                    ...) { ## ... are arguments used by functions called by the eval_replicate function
    #
    msg <- "Bootstrap replicates:"
    msglength <- nchar(msg) ## not used in the parallel case
    cat(msg)
    time1 <- Sys.time() 
    bootreps <- matrix(nrow=0,ncol=2)
    cumul_nsim <- 0L
    RNGstateList <- vector("list")
    boot.repl <- as.integer(boot.repl) ## impacts type of RNGstates...
    ii <- 0L ## 'global definition' (!)
    #nb_cores <- .check_nb_cores(nb_cores=nb_cores) # checked by calling fn
    if (nb_cores > 1L) {
      #cl <- parallel::makeCluster(nb_cores,outfile="essai.txt") 
      cl <- parallel::makeCluster(nb_cores) 
      R.seed <- get(".Random.seed", envir = .GlobalEnv)
      if (has_doSNOW <- ("package:doSNOW" %in% search())) { ## allows progressbar but then requires foreach
        # loading (?) the namespace of 'snow' changes the global RNG state!
        assign(".Random.seed", R.seed, envir = .GlobalEnv)
        fn <- get("registerDoSNOW", asNamespace("doSNOW"))
        do.call(fn,list(cl=cl)) 
        pb <- txtProgressBar(max = boot.repl, style = 3, char="P")
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        parallel::clusterExport(cl=cl, list("progress"),envir=environment()) ## slow! why?
      } else {
        if ( ! doSNOW_warned) {
          message("If the 'doSNOW' package were attached, better load-balancing might be possible.")
          doSNOW_warned <<- TRUE
        } 
        parallel::clusterEvalQ(cl, library("spaMM"))
      }
      dotenv <- list2env(list(...))
      parallel::clusterExport(cl=cl, as.list(ls(dotenv)),envir=dotenv) ## much faster...
    } else cl <- NULL
    repeat { 
      block_nsim <- boot.repl-nrow(bootreps)
      cumul_nsim <- cumul_nsim+block_nsim
      RNGstate <- get(".Random.seed", envir = .GlobalEnv)
      newy_s <- simulate(nullfit,nsim = block_nsim,verbose=FALSE) ## some replicates may not be analyzable  !
      if (block_nsim==1L) dim(newy_s) <- c(length(newy_s),1)
      if (nb_cores > 1L && has_doSNOW) {
        foreach_args <- list(
          ii = 1:boot.repl,
          .combine = "rbind",
          .inorder = TRUE,
          .packages = "spaMM",
          #.export = c("progress","simbData","bootlist"),
          .errorhandling = "remove",
          .options.snow = opts
        )
        foreach_blob <- do.call(foreach::foreach,foreach_args)
        bootblock <- foreach::`%dopar%`(foreach_blob, eval_replicate(newy_s[,ii]))
        close(pb)
      } else {
        bootblock <- pbapply(X=newy_s,MARGIN = 2L,FUN = eval_replicate, cl=cl)
        bootblock <- t(bootblock)
        bootblock <- stats::na.omit(bootblock)
      }
      if (is.null(bootblock)) {
        stop("All bootstrap replicates failed. Maybe a programming issue in parallel computation?")
      }
      if (nrow(bootblock)<block_nsim ) { ## eg separation in binomial models... 
        warning(paste("Analysis of",block_nsim-nrow(bootblock),"bootstrap samples out of",block_nsim," failed. 
                      Maybe no statistical information in them ?"))
        if (cumul_nsim>10*boot.repl) { ## to avoid an infinite loop
          stop("Analysis of bootstrap samples fails repeatedly. Maybe no statistical information in them ?")
        } ## otherwise repeat!
      }
      RNGstateList[[length(RNGstateList)+1L]] <- c(block_nsim=block_nsim, RNGstate)
      bootreps <- rbind(bootreps,bootblock)
      if (nrow(bootreps)>=boot.repl) break
  }
    if (nb_cores > 1L) {
      if (has_doSNOW) foreach::registerDoSEQ() ## https://stackoverflow.com/questions/25097729/un-register-a-doparallel-cluster
      parallel::stopCluster(cl) 
    } 
    cat(paste(" bootstrap took",.timerraw(time1),"s.\n")) 
    return(list(bootreps=bootreps,RNGstates=unlist(RNGstateList)))
  } ## end bootstrap
})

.preprocess_data <- function(data, null.formula, formula, resid.model, prior.weights, callargs) {
  if ( inherits(data,"list")) {
    for (lit in seq_along(data)) {
      null.validdata <- .getValidData(formula=null.formula[-2],resid.formula=resid.model$formula,data=data[[lit]],
                                      callargs=callargs["prior.weights"]) ## will remove rows with NA's in required variables
      full.validdata <- .getValidData(formula=formula[-2],resid.formula=resid.model$formula,data=data[[lit]],
                                      callargs=callargs["prior.weights"]) ## will remove rows with NA's in required variables
      data[[lit]] <- data[[lit]][intersect(rownames(null.validdata),rownames(full.validdata)),,drop=FALSE]
    }
  } else {
    null.validdata <- .getValidData(formula=null.formula[-2],resid.formula=resid.model$formula,data=data,
                                    callargs=callargs["prior.weights"]) ## will remove rows with NA's in required variables
    full.validdata <- .getValidData(formula=formula[-2],resid.formula=resid.model$formula,data=data,
                                    callargs=callargs["prior.weights"]) ## will remove rows with NA's in required variables
    data <- data[intersect(rownames(null.validdata),rownames(full.validdata)),,drop=FALSE]     
  }  
  return(data)
}

.LRT <- function(null.formula=NULL,formula,
                 null.disp=list(),REMLformula=NULL,boot.repl=0,
                 ## currently trace always false; this is not an argument t be forwarded as is to corrHLfit! 
                 #trace=FALSE, ## T means lead to calls of corrHLfit(... trace=list(<file name>,<over/append>))
                 verbose=c(trace=FALSE),
                 fittingFunction="corrHLfit",  
                 nb_cores=NULL,
                 data,
                 boot_fn="spaMM_boot",
                 resp_testfn=NULL,
                 .condition = NULL, ## only an argument of the internal function .LRT() so not visible at user level.
                 ...) {
  if (is.na(verbose["trace"])) verbose["trace"] <- FALSE
  if (is.na(verbose["all_objfn_calls"])) verbose["all_objfn_calls"] <- FALSE   
  callargs <- match.call(expand.dots = TRUE)
  ## We extract relevant arguments as promises (prior.weights, in particular, are held unevaluated)
  boot_call <- callargs[setdiff(names(callargs),names(formals(.LRT)))] 
  boot_call[[1L]] <- callargs[[1L]] ## we put back the called function (here .LRT, not fittingFunction) that was lost by subsetting
  ## birth pangs :
  if ("predictor" %in% names(callargs)) {
    stop(".LRT() called with 'predictor' argument which should be 'formula'" )
  }
  if ("null.predictor" %in% names(callargs)) {
    stop(".LRT() called with 'null.predictor' argument which should be 'null.formula'" )
  }
  ## here we makes sure that *predictor variables* are available for all data to be used under both models
  resid.model <- .reformat_resid_model(callargs$resid.model) ## family not easily checked here; not important
  data <- .preprocess_data(data=data, null.formula=null.formula, formula=formula, resid.model=resid.model,
                           callargs=callargs)
  boot_call$data <- data
  if ( inherits(data,"list")) {
    if ( ! is.null(boot_call$"prior.weights" <- data[[1L]]$"(weights)")) {
      warning(paste0("'data' is a list of data frames, and there are 'prior.weights':\n", # F I X M E
                  "the same prior.weights are applied to all data\n",
                  "which may be wrong or even fail in some cases."))
    }
  } else boot_call$"prior.weights" <- data$"(weights)" ## this is where .preprocess_data() -> .getValidData() -> model.frame() puts the evaluated weights.
  predictor <- formula   
  if (! inherits(formula,"predictor")) predictor <- .stripTerms(formula,data=data)
  null.predictor <- null.formula   
  if (! inherits(null.formula,"predictor")) null.predictor <- .stripTerms(null.formula,data=data)
  form <- predictor
  if (!is.null(boot_call$LamFix)) stop("LamFix is obsolete")
  if (!is.null(boot_call$PhiFix)) stop("PhiFix is obsolete")  
  boot_call$formula <- predictor 
  boot_call$verbose <- verbose 
  if (fittingFunction == "fitme") {
    if (is.null(boot_call$method)) stop("'method' argument is required when fittingFunction=\"fitme\".")
    locmethod <- boot_call$method
  } else {
    if (is.null(boot_call$HLmethod)) {
      boot_call$HLmethod <- boot_call$method
      boot_call$method <- NULL
    }
    locmethod <- boot_call$HLmethod
    if (fittingFunction=="HLfit") {
      FHFnames <- intersect(names(formals(HLfit)),names(boot_call))
      boot_call <- boot_call[FHFnames]
      boot_call[[1L]] <- callargs[[1L]] ## we put back the called function [here .LRT] that was lost by subsetting
    } ## HLfit has no ... args: we need to remove all arguments from the ... (FIXME: avoidable complication ?)
  }
  boot_call$nb_cores <- nb_cores ## not yet in boot_call bc nb_cores argument is in formals(.LRT) 
  fullm_call <- boot_call
  fullm_call$"prior.weights" <- callargs$"prior.weights" ## put back any promise that should be evaluated on the boot resample
  fullm_call[[1L]] <- get(fittingFunction, asNamespace("spaMM"))
  nullm_call <- fullm_call
  fullm_call$formula <- predictor
  
  #### "limited use, for initializing bootstrap replicates:"
  if ( ! is.null(null.predictor)) { ## ie if test effet fixe
    testFix <- TRUE
    if (locmethod =="SEM") {
      test_obj <- "logLapp"
    } else test_obj <- "p_v"
    ## check fullm.list$REMLformula, which will be copied into nullm in all cases of fixed LRTs
    if (locmethod %in% c("ML","PQL/L","SEM") || substr(locmethod,0,2) == "ML") {
      fullm_call$REMLformula <- NULL
      nullm_call$REMLformula <- NULL
    } else { ## an REML variant
      if (is.null(REMLformula)) { ## default
        fullm_call$REMLformula <- NULL ## HLfit will reconstruct proper REML from this
        nullm_call$REMLformula <- fullm_call$formula ## 
      } else { ## allows alternative choice, but still the same *both* the full and the null fit
        fullm_call$REMLformula <- REMLformula
        nullm_call$REMLformula <- REMLformula
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
    nullm_call$formula <- null.predictor
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
  
  #  trace.info <- NULL
  nullfit <- eval(nullm_call)
  canon.init <- attr(nullfit,"optimInfo")$LUarglist$canon.init ## includes user init
  # Laborious but matching the code in preprocess
  # true_corr_types <- c("adjacency","Matern","AR1","corrMatrix", "Cauchy")
  # corr_types <- true_corr_types[match(attr(nullfit$ZAlist, "exp_ranef_types"), true_corr_types)] ## full length
  # corr_families <- vector('list',length(corr_types))
  # for (rd in which( ! is.na(corr_types))) corr_families[[rd]] <- do.call(corr_types[rd],list())
  user_inits <- .post_process_parlist(nullm_call$init,corr_families=nullfit$corr_info$corr_families)
  #
  names_u_u_inits <- names(unlist(user_inits))
  names_u_c_inits <- names(unlist(canon.init))
  not_user_inits_names <- setdiff(names_u_c_inits, names_u_u_inits) 
  nullranPars <- get_ranPars(nullfit)
  names_u_nullranPars <- names(unlist(nullranPars))
  
  # at this point we have names of parameters that were outer optimized without an explicit user init
  if ( ! is.null(not_user_inits_names)) {
    removand <- setdiff(names_u_nullranPars, not_user_inits_names) ## removand: user_inits, fixed, or inner optimized corrPars
    if ( is.null(removand)) {
      locinit <- .modify_list(canon.init,nullranPars)
    } else { ## leaves user_inits as there are in LUarglist$canon.init, and do not add any fixed or inner-optimized par
      locinit <- .modify_list(canon.init,
                                      .remove_from_cP(nullranPars,u_names=removand)) ## loses attributes
    }
    if (fittingFunction=="fitme") {
      fullm_call$init <- locinit
    } else {
      fullm_call[["init.corrHLfit"]] <- locinit
    }
  }  else removand <- names_u_nullranPars
  fullfit <- eval(fullm_call)
  if (logLik(fullfit)<logLik(nullfit)) { ## evidence of fullfit being trapped in a local maximum
    ## We did not overwrite user inits: we do so
    if ( ! is.null(names_u_u_inits)) {
      removand <- setdiff(names_u_nullranPars, names_u_c_inits) ## removand2: fixed, or inner optimized corrPars
      if ( is.null(removand)) { ## may be true when this was previously false for removand, as there were user_inits
        locinit <- .modify_list(canon.init,nullranPars)
      } else { ## overwrites user_inits, and do not add any fixed or inner-optimized par
        locinit <- .modify_list(canon.init,
                                        .remove_from_cP(nullranPars,u_names=removand)) ## loses attributes
      }
      if (fittingFunction=="fitme") {
        fullm_call$init <- locinit
      } else {
        fullm_call[["init.corrHLfit"]] <- locinit
      }
      fullfit <- eval(fullm_call)
    } ## else it's not clear what to do here since we preemptively set the initial values of the full fit to the null fit values.
  }
  # # No comparable evidence that nullfit is trapped in a local maximum: check with fullfit ranPars
  fullranPars <- get_ranPars(fullfit) ## latest fulfit...
  if ( is.null(removand)) { ## ...and latest removand
    locinit <- .modify_list(canon.init,fullranPars)
  } else { ## leaves user_inits as there are in LUarglist$canon.init, and do not add any fixed or inner-optimized par
    locinit <- .modify_list(canon.init,
                                    .remove_from_cP(fullranPars,u_names=removand)) ## loses attributes
  }
  if (fittingFunction=="fitme") {
    nullm_call$init <- locinit
  } else {
    nullm_call[["init.corrHLfit"]] <- locinit ## but currently not used.
  }
  if (fittingFunction=="fitme") {
    renullfit <- eval(nullm_call)
    if (logLik(nullfit)<logLik(renullfit)) nullfit <- renullfit ## test may be FALSE ./. 
    # ./. (typically if fullfit yields a low lambda which is a local maximum of the nullfit)
  } ## ELSE currently not this recheck for corrHLfit
  
  if (testFix) {df <- fullfit$dfs[["pforpv"]]-nullfit$dfs[["pforpv"]]} else {df <- length(null.disp)} 
  if (df<0) {
    tmp <- fullfit
    fullfit <- nullfit
    nullfit <- tmp
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
    if (boot.repl>0L) {
      ## boot_call already has (full)formula and (optionally) ranFix
      more_boot_args <- list(null.formula=null.predictor,null.disp=null.disp,REMLformula=REMLformula,fittingFunction=fittingFunction)
      boot_call[names(more_boot_args)] <- more_boot_args ## unchanged user REMLformula forwarded
      boot_call$verbose <- c(trace=FALSE)
      boot_call$boot.repl <- 0L ## avoids recursive call of bootstrap
      if ( ! is.null(not_user_inits_names)) {
        removand <- setdiff(names_u_nullranPars, not_user_inits_names) ## removand: user_inits, fixed, or inner optimized corrPars
        if ( is.null(removand)) {
          boot_init <- .modify_list(canon.init,nullranPars)
        } else { ## leaves user_inits as there are in LUarglist$canon.init, and do not add any fixed or inner-optimized par
          boot_init <- .modify_list(canon.init, .remove_from_cP(nullranPars,u_names=removand)) ## loses attributes
        }
      }
      if (fittingFunction=="corrHLfit") {
        boot_call$init.corrHLfit <- boot_init
      } else if (fittingFunction=="fitme") {
        boot_call$init <- boot_init
      } #else boot_call$init.HLfit[notfixed] <- FIXME ## might do something here
      if (tolower(nullfit$family$family)=="binomial") {
        cbf <- .check_binomial_formula(nullfit=nullfit, data=data, fullfit=fullfit)
        cbindTest <- cbf$cbindTest
        if (cbindTest) {
          boot_call$null.formula <- cbf$null_formula
          boot_call$formula <- cbf$full_formula
          nnegname <- cbf$nnegname
          nposname <- cbf$nposname
        }
      } else cbindTest <- FALSE
      nb_cores <- .check_nb_cores(nb_cores=nb_cores)
      ## FIXME?  probable problem with prior.weights? except that prior.weights need not be in boot_call!
      if (nb_cores>1) for(st in names(boot_call)) if (st != "") boot_call[[st]] <- eval(boot_call[[st]]) ## force evaluation before running in another R session
      ## the data contain any original variable not further used; e.g original random effect values in the simulation tests  
      locitError <- 0
      #
      simbData <- nullfit$data
      eval_replicate <- function(y) {
        if (cbindTest) {
          simbData[[nposname]] <- y
          simbData[[nnegname]] <- .get_BinomialDen(nullfit)  - y
        } else {simbData[[as.character(nullfit$predictor[[2L]])]] <- y} ## allows y~x syntax for binary response
        boot_call$data <- simbData
        bootrepl <- eval(boot_call) 
        if (inherits(bootrepl,"try-error")) {
          return(c(NA,NA))
        } else {
          resu <- c(logLik(bootrepl$fullfit,which=test_obj),
                        logLik(bootrepl$nullfit,which=test_obj),
                    condition = with(bootrepl,eval(.condition))
                    )
          return(resu)
        }
      }
      if (boot_fn==".eval_boot_replicates") {
        bootblob <- .eval_boot_replicates(eval_replicate=function(newy) {eval_replicate(y=newy)},boot.repl=boot.repl,nullfit=nullfit,nb_cores=nb_cores,
                                          boot_call=boot_call, simbData=simbData # possibly exported to child processes for use by used by eval_replicates()
        )
      } else {
        bootblob <- spaMM_boot(object=nullfit,nsim = boot.repl,simuland=eval_replicate,nb_cores = nb_cores,
                               resp_testfn = resp_testfn, ## for simulate
                               boot_call=boot_call, simbData=simbData, .condition = .condition ## for exporting
                               # .condition is a .LRT() argument, used by the hard-coded eval_replicate function
        )
      }
      bootreps <- bootblob$bootreps
      #print(paste(boot_fn,paste(dim(bootreps),collapse=",")))
      # 
      colnames(bootreps) <- c(paste0(c("full.","null."),test_obj),"condition")[seq_len(ncol(bootreps))]
    } ## end bootstrap
  } else { ## nothing operativ yet
    warning("code missing here")
    bootreps <- matrix(NA,nrow=boot.repl,ncol=length(unlist(fullfit$APHLs))) 
    ## more needed here ?
  }
  ## prepare output
  if ( ! is.na(testFix)) {
    if (testFix) {df <- fullfit$dfs[["pforpv"]]-nullfit$dfs[["pforpv"]]} else {df <- length(null.disp)} 
    resu <- list(fullfit=fullfit,nullfit=nullfit)
    resu$basicLRT <- data.frame(chi2_LR=LRTori,df=df,p_value=1-pchisq(LRTori,df=df))
    if (boot.repl>0L) {
      bootreps <- data.frame(bootreps)
      bootreps$bootdL <- bootreps[,1L]-bootreps[,2L]
      rawPvalue <- (1+sum(bootreps$bootdL>=LRTori/2))/(boot.repl+1) ## DavisonH, p.141
      if ( ! is.null( .condition )) {
        lrfit <- lm(bootdL ~ condition ,data=bootreps)
        meanbootLRT <- 2*predict(lrfit,newdata=data.frame(condition=eval(.condition))) ## conditional mean
        attr(meanbootLRT,"boot_type") <- "conditional"
      } else {
        meanbootLRT <- 2*mean(bootreps$bootdL)
        attr(meanbootLRT,"boot_type") <- "marginal"
      }
      LRTcorr <- LRTori*df/meanbootLRT
      ## as documented in ?LRT
      resu$rawBootLRT <- data.frame(chi2_LR=LRTori,df=df,p_value=rawPvalue)
      resu$BartBootLRT <- data.frame(chi2_LR=LRTcorr,df=df,p_value=1-pchisq(LRTcorr,df=df))
      resu$bootInfo <- c(bootblob,list(meanbootLRT=meanbootLRT)) 
    }
  } else {
    resu <- list(fullfit=fullfit)
    resu$bootInfo <- list()
    resu$basicLRT <- list()
  }
  class(resu) <- c("fixedLRT",class(resu)) 
  return(resu)
}

get_boot_response <- function(object, replicate) {
  RNGstates <- object$bootInfo$RNGstates
  R.seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  if (is.matrix(RNGstates)) {
    cum_nsim <- cumsum(RNGstates[,1L])
    block <- Position(function(x) {x>=replicate},cum_nsim)
    assign(".Random.seed", RNGstates[block,-1L], envir = .GlobalEnv)
    nsim <- replicate-cum_nsim[block-1L]
    sim <- simulate(object$nullfit,nsim=nsim)
  } else {
    assign(".Random.seed", RNGstates[-1L], envir = .GlobalEnv)
    nsim <- replicate
  }
  sim <- simulate(object$nullfit,nsim=nsim,verbose=FALSE) 
  if (nsim>1L) {sim <- sim[,nsim]} else return(sim)
}

