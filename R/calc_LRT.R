.check_nb_cores <- function(nb_cores=NULL) {
  if (is.null(nb_cores)) nb_cores <- spaMM.getOption("nb_cores") ## may be NULL
  machine_cores <- parallel::detectCores()
  if (is.null(nb_cores)) {
    nb_cores <- 1L ## default
    if (machine_cores>1L && interactive()) {
      if (! identical(spaMM.getOption("cores_avail_warned"),TRUE)) {
        message(paste(machine_cores,
                      "cores are available for parallel computation\n(you may be allowed to fewer of them on a shared cluster).\nChange 'nb_cores' argument to use some of them.\nUse spaMM.options(nb_cores=<n>) to control nb_cores globally."))
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
  nform <- attr(nullfit$predictor,"oriFormula")  
  if (is.null(nform)) stop("a 'predictor' object must have an 'oriFormula' member.")
  if (paste(nform[[2L]])[[1L]]=="cbind") {
    ## We have different possible (exprL,exprR) arguments in cbind(exprL,exprR), 
    ## but in all case the predictor is that of exprL and exprR is $BinomialDen - exprL. We standardize: 
    res$nposname <- nposname <- .makenewname("npos",names(data))
    res$nnegname <- nnegname <- .makenewname("nneg",names(data))
    nform <- paste(nform)
    nform[2L] <- paste("cbind(",nposname,",",nnegname,")",sep="")
    res$null_formula <- as.formula(paste(nform[c(2,1,3)],collapse=""))
    fform <- paste(attr(fullfit$predictor,"oriFormula"))
    fform[2L] <- nform[2L]
    res$full_formula <- as.formula(paste(fform[c(2,1,3)],collapse=""))
    res$cbindTest <- TRUE
  } else res$cbindTest <- FALSE
  return(res)
}

.eval_boot_replicates <- function(eval_replicate, ## function to run on each replicate
                                  boot.repl, ## number of bootstrap replicates
                                  nullfit, ## the fitted model ibject for which bootstrap replicates are drawn 
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
  ii <- 0 ## 'global definition' (!)
  #nb_cores <- .check_nb_cores(nb_cores=nb_cores) # checked by calling fn
  if (nb_cores > 1L) {
    #cl <- parallel::makeCluster(nb_cores,outfile="essai.txt") 
    cl <- parallel::makeCluster(nb_cores) 
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    if (has_doSNOW <- ("doSNOW" %in% .packages() )) { ## allows progressbar but then requires foreach
      # loading (?) the namespace of 'snow' changes the global RNG state!
      assign(".Random.seed", R.seed, envir = .GlobalEnv)
      eval(as.call(c(quote(registerDoSNOW),list(cl=cl)))) 
      `%foreachdopar%` <- foreach::`%dopar%`
      pb <- txtProgressBar(max = boot.repl, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      parallel::clusterExport(cl=cl, list("progress"),envir=environment()) ## slow! why?
    } else if ( ! identical(spaMM.getOption("doSNOW_warned"),TRUE)) {
      message("If the 'doSNOW' package were attached, the progress of the bootstrap computation could be reported.")
      spaMM.options(doSNOW_warned=TRUE)
    } 
    dotenv <- list2env(list(...))
    parallel::clusterExport(cl=cl, as.list(ls(dotenv)),envir=dotenv) ## much faster...
  } 
  repeat { 
    block_nsim <- boot.repl-nrow(bootreps)
    cumul_nsim <- cumul_nsim+block_nsim
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    newy_s <- simulate(nullfit,nsim = block_nsim,verbose=FALSE) ## some replicates may not be analyzable  !
    if (block_nsim==1L) dim(newy_s) <- c(length(newy_s),1)
    #print(head(newy_s))
    if (nb_cores > 1L) {
      if (has_doSNOW) {
        bootblock <- foreach::foreach(
          ii = 1:boot.repl,
          .combine = "rbind",
          .inorder = TRUE,
          .packages = "spaMM",
          #.export = c("progress","simbData","bootlist"),
          .errorhandling = "remove",
          .options.snow = opts
        ) %foreachdopar% {
          eval_replicate(newy_s[,ii])
        }
        close(pb)
      } else {
        bootblock <- parallel::parApply(cl,newy_s,MARGIN = 2L,FUN = eval_replicate)
        bootblock <- t(bootblock)
        bootblock <- stats::na.omit(bootblock)
      }
    } else {
      eval_wrap <- function(v) {
        res <- eval_replicate(v)
        #if (res[1]< res[2]) {save(v,file="zut.rda");stop("ICI")}
        ii <<- ii+1
        tused <- .timerraw(time1)
        ttotal <- tused* boot.repl/ii
        if (interactive()) {
          for (bidon in 1:msglength) cat("\b")
          msg <<- paste("Estimated time remaining for bootstrap: ",signif(ttotal-tused,2)," s.",sep="")
          msglength <<- nchar(msg)
          cat(msg)
        } else {
          cat(ii);cat(" ")
          if ((ii %% 40)==0L) cat("\n")
        }
        return(res)
      }
      bootblock <- apply(newy_s,MARGIN = 2L,FUN = eval_wrap)
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
  if (nb_cores > 1L) { parallel::stopCluster(cl) } 
  cat(paste(" bootstrap took",.timerraw(time1),"s.\n")) 
  return(list(bootreps=bootreps,RNGstates=unlist(RNGstateList)))
} ## end bootstrap





.LRT <- function(null.formula=NULL,formula,
                     null.disp=list(),REMLformula=NULL,boot.repl=0,
                     ## currently trace always false; this is not an argument t be forwarded as is to corrHLfit! 
                     #trace=FALSE, ## T means lead to calls of corrHLfit(... trace=list(<file name>,<over/append>))
                     verbose=c(trace=FALSE),
                     fittingFunction="corrHLfit",  
                     nb_cores=NULL,
                     ...) {
  if (is.na(verbose["trace"])) verbose["trace"] <- FALSE
  if (is.na(verbose["all_objfn_calls"])) verbose["all_objfn_calls"] <- FALSE   
  callargs <- match.call(expand.dots = TRUE)
  ## as list(...) but holding elements (prior.weights, in particular) unevaluated
  dotlist <- callargs[setdiff(names(callargs),names(formals(.LRT)))] 
  dotlist$prior.weights <- NULL
  dotlist <- lapply(dotlist[-1L],eval,envir=parent.frame(1L))
  ## birth pangs :
  if ("predictor" %in% names(callargs)) {
    stop(".LRT() called with 'predictor' argument which should be 'formula'" )
  }
  if ("null.predictor" %in% names(callargs)) {
    stop(".LRT() called with 'null.predictor' argument which should be 'null.formula'" )
  }
  ## here we makes sure that *predictor variables* are available for all data to be used under both models
  data <- dotlist$data
  if ( ! is.null(callargs$resid.formula)) callargs$resid.model <- callargs$resid.formula
  resid.model <- .reformat_resid_model(callargs$resid.model) ## family not easily checked here; not important
  if ( inherits(data,"list")) {
    data <- lapply(data, function(dt) {
      null.validdata <- .getValidData(formula=null.formula[-2],resid.formula=resid.model$formula,data=dt,
                                     callargs=callargs["prior.weights"]) ## will remove rows with NA's in required variables
      full.validdata <- .getValidData(formula=formula[-2],resid.formula=resid.model$formula,data=dt,
                                     callargs=callargs["prior.weights"]) ## will remove rows with NA's in required variables
      dt[intersect(rownames(null.validdata),rownames(full.validdata)),,drop=FALSE]     
    })
  } else {
    null.validdata <- .getValidData(formula=null.formula[-2],resid.formula=resid.model$formula,data=data,
                                   callargs=callargs["prior.weights"]) ## will remove rows with NA's in required variables
    full.validdata <- .getValidData(formula=formula[-2],resid.formula=resid.model$formula,data=data,
                                   callargs=callargs["prior.weights"]) ## will remove rows with NA's in required variables
    data <- data[intersect(rownames(null.validdata),rownames(full.validdata)),,drop=FALSE]     
  }  
  dotlist$data <- data
  predictor <- formula   
  if (! inherits(formula,"predictor")) predictor <- Predictor(formula)
  null.predictor <- null.formula   
  if (! inherits(null.formula,"predictor")) null.predictor <- Predictor(null.formula)
  form <- predictor
  if (!is.null(dotlist$LamFix)) {
    dotlist$ranFix$lambda <- dotlist$LamFix
    dotlist$LamFix <- NULL
  }
  if (!is.null(dotlist$PhiFix)) {
    dotlist$ranFix$phi <- dotlist$PhiFix
    dotlist$PhiFix <- NULL
  }  
  dotlist$formula <- predictor 
  dotlist$verbose <- verbose 
  if (fittingFunction == "fitme") {
    if (is.null(dotlist$method)) stop("'method' argument is required when fittingFunction=\"fitme\".")
    locmethod <- dotlist$method
  } else {
    if (is.null(dotlist$HLmethod)) {
      dotlist$HLmethod <- dotlist$method
      dotlist$method <- NULL
    }
    locmethod <- dotlist$HLmethod
    if (fittingFunction=="HLfit") {
      FHFnames <- intersect(names(formals(HLfit)),names(dotlist))
      dotlist <- dotlist[FHFnames]
    } ## HLfit has no ... args: we need to remove all arguments from the ... (FIXME: avoidable complication ?)
  }
  dotlist$nb_cores <- nb_cores ## not yet in dotlist bc nb_cores argument is explicit for bootstrap 
  ## do not change dotlist afterwards !
  fullm.list <- dotlist
  nullm.list <- dotlist
  fullm.list$formula <- predictor
  #### "limited use, for initializing bootstrap replicates:"
  if ( ! is.null(null.predictor)) { ## ie if test effet fixe
    testFix <- TRUE
    if (locmethod =="SEM") {
      test.obj <- "logLapp"
    } else test.obj <- "p_v"
    ## check fullm.list$REMLformula, which will be copied into nullm in all cases of fixed LRTs
    if (locmethod %in% c("ML","PQL/L","SEM") || substr(locmethod,0,2) == "ML") {
      fullm.list$REMLformula <- NULL
      nullm.list$REMLformula <- NULL
    } else { ## an REML variant
      if (is.null(REMLformula)) { ## default
        fullm.list$REMLformula <- NULL ## HLfit will reconstruct proper REML from this
        nullm.list$REMLformula <- fullm.list$formula ## 
      } else { ## allows alternative choice, but still the same *both* the full and the null fit
        fullm.list$REMLformula <- REMLformula
        nullm.list$REMLformula <- REMLformula
      }
      namesinit <- names(fullm.list$init.corrHLfit)
      namesinit <- setdiff(namesinit,c("rho","nu","Nugget","ARphi"))
      len <- length(namesinit)
      if ( len>0) {
        if (len > 1) namesinit <- paste(c(paste(namesinit[-len],collapse=", "),namesinit[len]),collapse=" and ")
        message("Argument 'init.corrHLfit' is used in such a way that")
        message(paste("  ",namesinit," will be estimated by maximization of p_v.",sep=""))
        message("  'REMLformula' will be inoperative if all dispersion")
        message("  and correlation parameters are estimated in this way.")
      }
    }
    nullm.list$formula <- null.predictor
  } else if ( length(null.disp)>0 ) { ## test disp/corr param
    testFix <- FALSE
    #
    stop("Models differ in their random effects or residual dispersion structure.")
    #
    namescheck <- names(fullm.list$lower)
    namescheck <- namescheck[ ! namescheck %in% names(null.disp)]  
    nullm.list$lower <- nullm.list$lower[namescheck]
    nullm.list$upper <- nullm.list$upper[namescheck]
    nullm.list$ranFix <- c(nullm.list$ranFix,null.disp) ## adds fixed values to any preexisting one
  } else testFix <- NA
  
  #  trace.info <- NULL
  nullfit <- do.call(fittingFunction,nullm.list)
  nullranPars <- nullfit$CorrEst_and_RanFix
  notfixed <- names(which(attr(nullranPars,"type")!="fix"))
  if (fittingFunction=="fitme") {
    notfixed_notinit <- setdiff(notfixed,names(fullm.list$init)) ## not overwrite a user-explicit init
    fullm.list$init[notfixed_notinit] <- nullranPars[notfixed_notinit] 
  }
  fullfit <- do.call(fittingFunction,fullm.list)
  if (fittingFunction=="fitme") {
    if (logLik(fullfit)<logLik(nullfit)) { ## evidence of fullfit being trapped in a local maximum
      nullranPars <- nullfit$CorrEst_and_RanFix
      fullm.list$init[notfixed_notinit] <- nullranPars[notfixed_notinit]
      fullfit <- do.call(fittingFunction,fullm.list)
    }
    # # No comparable evidence that nullfit is trapped in a local maximum: check with fullfit ranPars
    fullranPars <- fullfit$CorrEst_and_RanFix
    nullm.list$init[notfixed_notinit] <- fullranPars[notfixed_notinit]
    renullfit <- do.call(fittingFunction,nullm.list)
    if (logLik(nullfit)<logLik(renullfit)) nullfit <- renullfit ## test may be FALSE ./. 
    # ./. (typically if fullfit yields a low lambda which is a local maximum of the nullfit)
  }
  if (testFix) {df <- length(fullfit$fixef)-length(nullfit$fixef)} else {df <- length(null.disp)}
  if (df<0) {
    tmp <- fullfit
    fullfit <- nullfit
    nullfit <- tmp
  }
  if (inherits(fullfit,"HLfitlist")) {
    fullL <- attr(fullfit,"APHLs")[[test.obj]]
  } else fullL <- fullfit$APHLs[[test.obj]]
  if (inherits(nullfit,"HLfitlist")) {
    nullL <- attr(nullfit,"APHLs")[[test.obj]]
  } else nullL <- nullfit$APHLs[[test.obj]]
  LRTori <- 2*(fullL-nullL)
  ## BOOTSTRAP
  if ( ! is.na(testFix)) {
    if (boot.repl>0) {
      bootlist <- dotlist ## copies (full)formula and (optionally) ranFix
      bootlist <- c(bootlist,list(null.formula=null.predictor,null.disp=null.disp,REMLformula=REMLformula,fittingFunction=fittingFunction)) ## unchanged user REMLformula forwarded
      bootlist$verbose <- c(trace=FALSE)
      #bootlist$trace <- FALSE 
      bootlist$boot.repl <- 0 ## avoids recursive call of bootstrap
      if (fittingFunction=="corrHLfit") {
        bootlist$init.corrHLfit[notfixed] <- nullfit$CorrEst_and_RanFix[notfixed]
      } else if (fittingFunction=="fitme") {
        bootlist$init <- nullfit$CorrEst_and_RanFix[notfixed]
      } else bootlist$init.HLfit[notfixed] <- nullfit$CorrEst_and_RanFix[notfixed]
      if (tolower(nullfit$family$family)=="binomial") {
        cbf <- .check_binomial_formula(nullfit=nullfit, data=data, fullfit=fullfit)
        cbindTest <- cbf$cbindTest
        if (cbindTest) {
          bootlist$null.formula <- cbf$null_formula
          bootlist$formula <- cbf$full_formula
          nnegname <- cbf$nnegname
          nposname <- cbf$nposname
        }
      } else cbindTest <- FALSE
      nb_cores <- .check_nb_cores(nb_cores=nb_cores)
      if (nb_cores>1) for(st in names(bootlist)) bootlist[[st]] <- eval(bootlist[[st]]) ## force evaluation before running in another R session
      
      ## the data contain any original variable not further used; e.g original random effect values in the simulation tests  
      thisFnName <- as.character(sys.call()[[1]]) ## prevents a bug when we change "this" function name
      thisFnName <- strsplit(thisFnName,":")[[1]]
      thisFnName <- thisFnName[length(thisFnName)]
      locitError <- 0
      #
      simbData <- nullfit$data
      eval_replicate <- function(newy,only_vector=TRUE) {
        if (cbindTest) {
          simbData[[nposname]] <- newy
          simbData[[nnegname]] <- .get_BinomialDen(nullfit)  - newy
        } else {simbData[[as.character(nullfit$predictor[[2L]])]] <- newy} ## allows y~x syntax for binary response
        bootlist$data <- simbData
        bootrepl <- try(do.call(thisFnName,bootlist)) ## appears to fail IF thisFnName is internal (.do_LRT)
        #bootrepl <- try(eval(as.call(c(quote(thisFnName),bootlist)))) ## never worked (without the strsplit at least)
        if (inherits(bootrepl,"try-error")) {
          if (only_vector) {
            return(c(NA,NA))
          } else return(bootrepl)
        } else return(c(logLik(bootrepl$fullfit,which=test.obj),
                        logLik(bootrepl$nullfit,which=test.obj)))
      }
      bootblob <- .eval_boot_replicates(eval_replicate=eval_replicate,boot.repl=boot.repl,nullfit=nullfit,nb_cores=nb_cores,
                                        bootlist=bootlist, simbData=simbData)
      bootreps <- bootblob$bootreps
    } ## end bootstrap
  } else { ## nothing operativ yet
    warning("code missing here")
    bootreps <- matrix(NA,nrow=boot.repl,ncol=length(unlist(fullfit$APHLs))) 
    colnames(bootreps) <- names(unlist(fullfit$APHLs))
    ## more needed here ?
  }
  ## prepare output
  if ( ! is.na(testFix)) {
    if (testFix) {df <- length(fullfit$fixef)-length(nullfit$fixef)} else {df <- length(null.disp)}
    resu <- list(fullfit=fullfit,nullfit=nullfit)
    resu$basicLRT <- data.frame(chi2_LR=LRTori,df=df,p_value=1-pchisq(LRTori,df=df))
    if (boot.repl>0) {
      bootdL <- bootreps[,1]-bootreps[,2]
      meanbootLRT <- 2*mean(bootdL)  
      rawPvalue <- (1+sum(bootdL>=LRTori/2))/(boot.repl+1) ## DavisonH, p.141
      LRTcorr <- LRTori*df/meanbootLRT
      ## as docuumented in ?LRT
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

