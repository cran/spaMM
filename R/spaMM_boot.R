simulate4boot <- function(object, nsim, seed=NULL, resp_testfn=NULL, type="marginal", showpbar=eval(spaMM.getOption("barstyle")),
                          ...) {
  RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  msg <- "Bootstrap replicates: "
  msglength <- nchar(msg) ## not used in the parallel case
  cat(msg)
  cumul_nsim <- 0L
  nsim <- as.integer(nsim) 
  if (nsim<1L) {
    warning("'nsim' must be at least 1.")
    return(list())
  }
  ####
  newy_s <- simulate(object,nsim = nsim,verbose=c(type=FALSE,showpbar=showpbar), resp_testfn=resp_testfn, type=type, seed=seed, ...) 
  if (nsim==1L) dim(newy_s) <- c(length(newy_s),1L)
  list(bootreps=newy_s,RNGstates=RNGstate)
}

spaMM_boot <- function(object, simuland, nsim, nb_cores=NULL,
           seed=NULL,
           resp_testfn=NULL, control.foreach=list(),
           debug. = FALSE, 
           type,
           fit_env=NULL,
           cluster_args=NULL,
           showpbar= eval(spaMM.getOption("barstyle")),
           boot_samples=NULL,
           ...) {
  if (missing(type) && is.null(boot_samples)) {
    warning("'type' is now a mandatory argument of spaMM_boot() when 'boot_samples' remains NULL.\n Assuming type='marginal' for consistency with previous versions.",
            immediate. = TRUE)
    type <- "marginal"
  }
  if ("fn" %in% ...names()) warning("spaMM_boot() expects 'simuland', not 'fn'.", immediate. = TRUE)
  if (is.null(boot_samples)) {
    boot_samples <- simulate4boot(object=object, nsim=nsim, seed=seed, resp_testfn=resp_testfn, type=type, showpbar=showpbar)
  } else if ( ! (is.list(boot_samples) && ! is.null(boot_samples$bootreps))) { # if not like a simulate4boot() result
    boot_samples <- list(bootreps=boot_samples)
  }
  #
  # If the simuland has (say) arguments y, what=NULL, lrt, ...   , we should not have lrt in the dots. Since the dots are not directly manipulable
  # we have to convert them to a list, and ultimately to use do.call()
  control.foreach$.combine <- "rbind"
  wrap_parallel <- get(.spaMM.data$options$wrap_parallel, asNamespace("spaMM"), inherits=FALSE) # dopar
  boot_samples$bootreps <- wrap_parallel(newresp = boot_samples$bootreps, nb_cores = nb_cores, # wrap_parallel() is typically dopar()
                                         fn = simuland, fit_env = fit_env,   
                                         control=control.foreach, debug.=debug., pretest_cores = .pretest_fn_on_cores, 
                                         showpbar = showpbar, ...) 
  return(boot_samples)
}

spaMM2boot <- function(object, statFUN, nsim, nb_cores=NULL,
           seed=NULL,
           resp_testfn=NULL, control.foreach=list(),
           debug. = FALSE, 
           type="marginal",
           fit_env=NULL,
           cluster_args=NULL,
           showpbar= eval(spaMM.getOption("barstyle")),
           boot_samples=NULL,
           ...) {
  if ("fn" %in% ...names()) warning("spaMM2boot() expects 'simuland', not 'fn'.", immediate. = TRUE)
  if (is.null(boot_samples))  {
    boot_samples <- simulate4boot(object=object, nsim=nsim, seed=seed, resp_testfn=resp_testfn, type=type, showpbar=showpbar)
  } else if ( ! (is.list(boot_samples) && ! is.null(boot_samples$bootreps))) { # if not like a simulate4boot() result
    boot_samples <- list(bootreps=boot_samples)
  }
  #
  # If the statFUN has (say) arguments y, what=NULL, lrt, ...   , we should not have lrt in the dots. Since the dots are not directly manipulable
  # we have to convert them to a list, and ultimately to use do.call()
  control.foreach$.combine <- "rbind"
  wrap_parallel <- get(.spaMM.data$options$wrap_parallel, asNamespace("spaMM"), inherits=FALSE) # dopar
  simuland <- function(y, ...) {
    .refit <- update_resp(object, y) # pb if '.refit=' is also in the dots
    statFUN(.refit, ...)  
  }
  bootreps <- wrap_parallel(newresp = boot_samples$bootreps, nb_cores = nb_cores, # wrap_parallel() is typically dopar()
                                         fn = simuland, fit_env = fit_env, 
                                         control=control.foreach, debug.=debug., pretest_cores = .pretest_fn_on_cores, 
                                         showpbar = showpbar, ...) 
  return(list(t=t(bootreps), t0=statFUN(refit=object, ...), R=nsim, sim="parametric", .Random.seed=boot_samples$RNGstate))
}

.set_progrbar <- function(style, ...) {
  if (style==0L) {
    setup <- list(progress=NULL)
  } else {
    pb <- txtProgressBar(style=style, ...)
    progress <- function(n) setTxtProgressBar(pb, n)
    setup <- list(pb=pb, progress = progress)
  }
  setup
}

.warn_once_progressr <- local({
  progressr_warned <- FALSE
  function() {
    if ( ! progressr_warned) {
      message("If the 'progressr' package were attached, a progress bar would be available.")
      progressr_warned <<- TRUE
    } 
  }
})

.check_call_on_core <- function(fitobject) { # this is run on a node
  obj_call <- getCall(fitobject)
  expr_list <- as.list(obj_call)
  checkand <- setdiff(names(expr_list),c("","formula","data","prior.weights"))
  errors <- list()
  for (st in checkand) {
    arg <- expr_list[[st]] # an expression
    chk <- try(eval(arg))
    if (inherits(chk, "try-error")) errors[[st]] <- attr(chk,"condition")$message
  }
  errors  <- unlist(errors)
  return(errors)
}

.pretest_fn_on_cores <- function(fn, cluster) { # this is run on the parent process 
  errors_null <- errors_full <- NULL
  pbopt <- pboptions(type="none") 
  nullfit <- environment(fn)$nullfit
  if ( ! is.null(nullfit)) errors_null <- pbreplicate(1, .check_call_on_core(fitobject = nullfit), cl=cluster)
  #foreach_blob <- foreach::foreach(i=1)
  #if ( ! is.null(nullfit)) errors_null <- foreach::`%dopar%`(foreach_blob, .check_call_on_core(fitobject = nullfit))
  fullfit <- environment(fn)$fullfit
  if ( ! is.null(fullfit)) errors_full <- pbreplicate(1, .check_call_on_core(fitobject = fullfit), cl=cluster)
  #if ( ! is.null(fullfit)) errors_null <- foreach::`%dopar%`(foreach_blob, .check_call_on_core(fitobject = fullfit))
  errors <- unique(unlist(c(errors_null,errors_full)))
  if (length(errors)) { 
    errmess <- paste0("'",unlist(lapply(strsplit(errors, split="'"),`[`,i=2)),"'", collapse=", ")
    errmess <- paste("Object(s)",errmess,"not found on cluster node:\n add them to 'fit_env' argument? (see ?spaMM_boot for details).")
    warning(.spaMM.data$options$stylefns$hardwarn(errmess),immediate.=TRUE)
  }
  pboptions(pbopt) 
}  # warnings on the parent process, no return value

.set_cluster_type <- function(cluster_args, nb_cores) {
  if (is.null(cluster_args$spec)) cluster_args$spec <- .check_nb_cores(nb_cores=nb_cores)
  if (cluster_args$spec>1L) {
    type_user <- cluster_args$type
    if (.Platform$OS.type == "windows") {
      if (is.null(type_user)) {
        cluster_args$type <- "PSOCK" # default, but explicit => can be tested # On windows, or on linux if explicitly requested
      } else if (type_user=="FORK") {
        message('cluster_args$type=="FORK" not feasible under Windows')
        cluster_args$type <- "PSOCK"
      }
    } else { # linux alikes
      if (is.null(type_user)) {
        cluster_args$type <- "PSOCK" # default is socket cluster in all cases (the doc says so)
      } else if (type_user=="FORK" && .inRstudio(silent=TRUE)) {
        message('cluster_args$type=="FORK" not feasible when R is called from an Rstudio session.')
        cluster_args$type <- "PSOCK"
      }
    }
  }
  cluster_args
}

# fn more generic than spaMM_boot: there is no call to other spaMM fns such as simulate(object, .) so this acts as a general wrapper for 
# foreach or pbapply, and not specifically for bootstrap computations.
dopar <- local({
  doSNOW_warned <- FALSE
  function(newresp, fn, nb_cores=NULL, 
           fit_env, control=list(), cluster_args=NULL,
           debug.=FALSE, iseed=NULL, showpbar=eval(spaMM.getOption("barstyle")),
           pretest_cores=NULL,
           ... # passed to fn... unless captured by pbapply (in which case 'simplify' may have a distinct effect).
  ) {
    if (is.list(fit_env)) fit_env <- list2env(fit_env)
    cluster_args <- .set_cluster_type(cluster_args, nb_cores) # PSOCK vs FORK
    nb_cores <- cluster_args$spec
    if (debug. && nb_cores>1L ) debug. <- 1L 
    assign("debug.", debug., environment(fn))
    if (is.null(dim(newresp))) newresp <- matrix(seq(newresp),ncol=newresp,nrow=1) # assuming newresp is an integer
    nsim <- ncol(newresp)
    time1 <- Sys.time() 
    if (nb_cores>1L) {
      if ( ! is.null(iseed) ) {
        ori <- RNGkind("L'Ecuyer-CMRG")
        set.seed(iseed)
      }
      if (cluster_args$type=="FORK") {
        if (is.null(mc.silent <- control$mc.silent)) mc.silent <- TRUE 
        if (is.null(mc.preschedule <- control$mc.preschedule)) mc.preschedule <- TRUE 
        has_progressr <- ("progressr" %in% loadedNamespaces())
        seq_nr <- seq_len(nsim)
        if (has_progressr) {
          # progressor is the only progress function that 'works' with mclapply
          # although not with load-balancing (mc.preschedule=FALSE)
          # Here we use the default (no balancing), and it is the steps with max value shown below that are reported.  
          prog_fn <- get("progressor", asNamespace("progressr"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
          with_fn <- get("with_progress", asNamespace("progressr"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
          with_fn({
            p <- prog_fn(steps=ceiling(nsim/nb_cores))
            p_fn <- function(it, ...) { # it OK for mclapply... not for apply on a matrix
              res <- fn(newresp[,it], ...)
              p() # p() call necessary for actual progress report 
              res
            }
            bootreps <- try(
              parallel::mclapply(seq_nr, FUN = p_fn, mc.silent=mc.silent, mc.cores=nb_cores,
                                 mc.preschedule = mc.preschedule)
            )
          })
          
        } else {
          .warn_once_progressr()
          bootreps <- try(
            parallel::mclapply(seq_nr, FUN = function(it) fn(newresp[,it], ...), mc.silent=mc.silent, mc.cores=nb_cores)
          )
        }
        if (identical(control$.combine,"rbind")) {
          bootreps <- do.call(rbind,bootreps)
        } else bootreps <- do.call(cbind,bootreps)
      } else { # PSOCK
        cl <- do.call(parallel::makeCluster, cluster_args) # note that _this_ line would make sense for fork clusters too. BUT
        # ... the foreach = dot args combination may not work for FORK type. Only pbapply would work with makeCluster+FORK, 
        # but pbmcapply is a better way to get a pb one a fork cluster as [pb]mclapply have better load balancing than pbapply. 
        # has_doSNOW <- ("package:doSNOW" %in% search()) # result of library()
        has_doSNOW <- ("doSNOW" %in% loadedNamespaces()) # result of library() or requireNamespace()
        if (has_doSNOW) {
          # loading (?) the namespace of 'snow' changes the *parent* RNG state (as well as sons' ones)! so we save and restore it 
          R.seed <- get(".Random.seed", envir = .GlobalEnv) # save parent RNG state
          rdS_fn <- get("registerDoSNOW", asNamespace("doSNOW"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
          do.call(rdS_fn,list(cl=cl)) # this is what makes foreach see it and perform parallel computations
          assign(".Random.seed", R.seed, envir = .GlobalEnv) # restore parent RNG state
          if ( ! is.null(iseed) ) parallel::clusterSetRNGStream(cl = cl, iseed) 
          #
          # if (cluster_args$type == "PSOCK") {
            if (is.environment(fit_env)) parallel::clusterExport(cl=cl, varlist=ls(fit_env), envir=fit_env) 
            pb_char <- "P"
          # } else pb_char <- "F"
          # A first foreach_blob for a first dopar before defining the progress bar (otherwise we see a progress bar on this dopar)
          i <- NULL ## otherwise R CMD check complains that no visible binding for global variable 'i' (in expression newy_s[,i])
          foreach_blob <- foreach::foreach(i=1:nb_cores)
          #if (cluster_args$type == "PSOCK") {
            abyss <- foreach::`%dopar%`(foreach_blob, Sys.setenv(LANG = "en")) # before setting the progress bar...
            if (is.function(pretest_cores)) pretest_cores(fn, cl)
          #}
          # define the progress bar:
          barstyle <- eval(spaMM.getOption("barstyle"))
          progrbar_setup <- .set_progrbar(max = nsim, style = barstyle, char=pb_char)
          # :where opts are needed to define a second foreach_blob
          foreach_args <- list( 
            i = 1:nsim, 
            .combine = "cbind", 
            .inorder = TRUE, .packages = "spaMM", 
            .errorhandling = "remove", ## use "pass" to see problems
            .options.snow = progrbar_setup["progress"]
          )
          foreach_args[names(control)] <- control # replaces the above defaults by user controls
          foreach_blob <- do.call(foreach::foreach,foreach_args) 
          if (TRUE) {
            fn_dots <- list(...)
            for (st in names(fn_dots)) {
              # Add an enclosing quote():
              if ( is.language(fn_dots[[st]])) fn_dots[[st]] <- substitute(quote(what),list(what=fn_dots[[st]]))
            }
            bootreps <- try(foreach::`%dopar%`(foreach_blob, do.call(fn, c(list(newresp[, i]), fn_dots)))) 
          } else {
            # Standard passing of the dots with foreach does not seem to work. (good test is the doSNOW case nested within test-LRT-boot.R)
            # bootreps <- try(foreach::`%dopar%`(foreach_blob, fn(newresp[, i], ...)))
          }
          # the try() is useful if the user interrupts the dopar, in which case it allows close(pb) to be run. (? But doSNOW appear to close the nodes asynchronously?)
          foreach::registerDoSEQ() ## https://stackoverflow.com/questions/25097729/un-register-a-doparallel-cluster
          parallel::stopCluster(cl)
          #
          if (foreach_args[[".errorhandling"]]=="remove" && is.null(bootreps)) {
            cat(crayon::bold(paste0(
              "Hmmm. It looks like all parallel processes failed. Maybe rerun spaMM_boot() \n",
              "with  ' control.foreach=list(.errorhandling=\"stop\") '  to diagnose the problem.\n"
            )))
          } else if (foreach_args[[".errorhandling"]]=="stop" && inherits(bootreps,"try-error")) {            
            # foreach alters the condition message => seel '\"' after 'could not find'
            if (length(grep("could not find",(condmess <- conditionMessage(attr(bootreps,"condition")))))) {
              firstpb <- strsplit(condmess,"could not find")[[1]][2]
              firstpb <- strsplit(firstpb,"\"")[[1]][2]
              cat(crayon::bold(paste0(
                "Hmmm. It looks like some variables were not passed to the parallel processes.\n",
                "Maybe add   ' ",firstpb," = ",firstpb," '  to spaMM_boot()'s 'fit_env' argument?\n"
              )))
            } else cat(crayon::bold(condmess))
          }
          #
          if (showpbar) close(progrbar_setup$pb)
        } else { # no doSNOW
          if ( ! doSNOW_warned) {
            message("Note: If the 'doSNOW' package were attached, better load-balancing might be possible.")
            doSNOW_warned <<- TRUE
          } 
          pb_char <- "p"
          parallel::clusterCall(cl, Sys.setenv, LANG = "en")
          if ( ! is.null(iseed) ) parallel::clusterSetRNGStream(cl = cl, iseed) 
          packages2export <- control$.packages
          if (is.null(packages2export)) packages2export <- "spaMM"
          parallel::clusterCall(cl,
                                function(packages) {for (p in packages) library(p, character.only = TRUE)}, 
                                packages2export)
          if (is.environment(fit_env)) try(parallel::clusterExport(cl=cl, varlist=ls(fit_env), envir=fit_env)) 
          
          # in that case, ## We will use pbapply, with argument cl=cl; 
          # Given no doSNOW, a direct call to foreach would require doParallel::registerDoParallel(cl)
          # or doFuture::registerDoFuture(). 
          # There are fake solutions suggesting a progress bar can be set up with doParallel, 
          # but it actually progresses only after all the processes have been run
          # (proposed examples have too short processes for this to be apparent).
          # So ultimately... we need doFuture, ( => see distinct wrapper).
          # and current we use pbapply with an ad hoc treatment for combining the results 
          
          if (is.function(pretest_cores)) pretest_cores(fn, cl)
          if (showpbar) {
            pbopt <- pboptions(nout=min(100L,2L*nsim),type="timer",char=pb_char) 
          } else pbopt <- pboptions(type="none") 
          #try() so that an interrupt does not prevent running stopCluster():
          bootreps <- try(pbapply(X=newresp,MARGIN = 2L,FUN = fn, cl=cl, ...))
          parallel::stopCluster(cl)
          pboptions(pbopt)
          if (inherits(bootreps,"try-error")) {
            if (length(grep("could not find",(condmess <- conditionMessage(attr(bootreps,"condition")))))) {
              firstpb <- strsplit(condmess,"\"")[[1]][2]
              cat(crayon::bold(paste0(
                "Hmmm. It looks like some variables were not passed to the parallel processes.\n",
                "Maybe add   ' ",firstpb," = ",firstpb," '  to spaMM_boot()'s 'fit_env' argument?\n"
              )))
            } else cat(crayon::bold(condmess))
          }
          # LRT -> spaMM_boot -> eval_replicate with debug.=TRUE and not doSNOW can return more elaborate objects in case of error.
          # But these should not be diagnosed in this generic function.
          if (identical(control$.combine,"rbind")) bootreps <- t(bootreps) # this means the pbapply version handles cbind or rbind but not other 
        } # has_doSNOW ... else
      } # FORK ... else
      if ( ! is.null(iseed) ) do.call("RNGkind", as.list(ori)) # reste to state pre-parallel computation
    } else { ## nb_cores=1L
      pb_char <- "s"
      if (FALSE) {
        # :where opts are needed to define a second foreach_blob
        foreach_args <- list( 
          i = 1:ncol(newresp), 
          .combine = "cbind", 
          .inorder = TRUE, .packages = "spaMM", 
          .errorhandling = "remove" ## use "pass" to see problems
        )
        
        foreach_args[names(control)] <- control # replaces the above defaults by user controls
        
        if (showpbar) { # optionally wrap the combine function with progress bar code
          barstyle <- eval(spaMM.getOption("barstyle"))
          
          .combine <- foreach_args$.combine
          if (inherits(.combine,"character")) .combine <- get(.combine)
          
          progrbar_setup <- .set_progrbar(max = nsim, style = barstyle, char=pb_char)
          combine_with_pb <- function(nsim, pb){
            count <- 0
            force(pb)
            function(...) {
              count <<- count + length(list(...)) - 1L
              setTxtProgressBar(pb, count)
              flush.console()
              cbind(...) # this can feed into .combine option of foreach
            }
          } # returns a function that increments the bar then actually combines.
          foreach_args$.combine <- combine_with_pb(nsim, pb=progrbar_setup$pb)
        } 
        
        foreach_blob <- do.call(foreach::foreach,foreach_args) 
        
        fn_dots <- list(...)
        # for (st in names(fn_dots)) {
        #   # Add an enclosing quote():
        #   if ( is.language(fn_dots[[st]])) fn_dots[[st]] <- substitute(quote(what),list(what=fn_dots[[st]]))
        # }
        bootreps <- try(foreach::`%do%`(foreach_blob, do.call(fn, c(list(newresp[, i]), fn_dots)))) 
        # the try() is useful if the user interrupts the %do%, in which case it allows close(pb) to be run.
        
        if (showpbar) close(progrbar_setup$pb)
      } else { # older version using pbapply
        if (showpbar) {
          pbopt <- pboptions(nout=min(100L,2L*nsim),type="timer",char=pb_char) 
        } else pbopt <- pboptions(type="none") 
        bootreps <- pbapply(X=newresp,MARGIN = 2L,FUN = fn, cl=NULL, ...)
        pboptions(pbopt)
        if (identical(control$.combine,"rbind")) bootreps <- t(bootreps)
      }
    }
    cat(paste(" bootstrap took",.timerraw(time1),"s.\n")) 
    return(bootreps)
  }
})


.init_cores <- local({
  doSNOW_warned <- FALSE
  function(cluster_args=list()) { 
    cluster_args$spec <- .check_nb_cores(nb_cores=cluster_args$spec) # if cluster_args was NULL it is converted to list here => no need for special handling code.
    cores_info <- list(nb_cores=cluster_args$spec)
    #
    if (cluster_args$spec > 1L) {
      cores_info$cl <- do.call(parallel::makeCluster, cluster_args) 
      #dotenv <- list2env(list(...))
      #parallel::clusterExport(cl=cores_info$cl, as.list(ls(dotenv)),envir=dotenv) 
      ## foreach is NOT a parallel backend so there is no point using it if doSNOW is not available
      if (cluster_args$type!="FORK") {
        if (cores_info$has_doSNOW <- (isNamespaceLoaded("doSNOW"))) {
          R.seed <- get(".Random.seed", envir = .GlobalEnv)
          ## allows progressbar but then requires foreach
          assign(".Random.seed", R.seed, envir = .GlobalEnv) # loading (?) the namespace of 'snow' changes the global RNG state!
          fn <- get("registerDoSNOW", asNamespace("doSNOW"))
          do.call(fn,list(cl=cores_info$cl)) 
        } else {
          if ( ! doSNOW_warned) {
            message("If the 'doSNOW' package were attached, better load-balancing might be possible (at the expense of control of RNG).")
            doSNOW_warned <<- TRUE
          } 
        }
      } else cores_info$has_doSNOW <- FALSE
    }
    return(cores_info)
  }
})

.close_cores <- function(cores_info) {
  if ( cores_info$nb_cores > 1L) {
    if (cores_info$has_doSNOW) foreach::registerDoSEQ() ## https://stackoverflow.com/questions/25097729/un-register-a-doparallel-cluster
    parallel::stopCluster(cores_info$cl)
  }
}

# This one is derived from the Infusion code, not from spaMM::dopar(). 
# It has no serial subcase. (___F I X M E___)
# There is more control of the iterator, which needs not be a seq of integers
# and this this can apply easily to non-matrix objects.
# .dopar() needs more testing before going public ____F I X M E____
.dopar <- function(iterator, FUN, cluster_args = cluster_args, # required argS
                   # FUN args passed as the \dots, apart from its FUN argument (each element of 'iterator'):
                   ...) {
  
  cluster_args <- .set_cluster_type(cluster_args=cluster_args) # PSOCK vs FORK
  cores_info <- .init_cores(cluster_args=cluster_args)
  if (cluster_args$type=="FORK") {
    nb_cores <- cores_info$nb_cores
    # if (is.null(mc.silent <- control$mc.silent)) # conflict with 'control' arg of FUN=..slice_n_predict_body(). Other name needed
      mc.silent <- TRUE 
    # if (is.null(mc.preschedule <- control$mc.preschedule)) 
      mc.preschedule <- TRUE 
    has_progressr <- ("package:progressr" %in% search())
    if (has_progressr) {
      # progressor is the only progress function that 'works' with mclapply
      # although not with load-balancing (mc.preschedule=FALSE)
      # Here we use the default (no balancing), and it is the steps with max value shown below that are reported.  
      prog_fn <- get("progressor", asNamespace("progressr"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
      with_fn <- get("with_progress", asNamespace("progressr"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
      with_fn({
        p <- prog_fn(steps=length(iterator))
        p_FUN <- function(it, ...) {
          res <- FUN(it, ...)
          p() # p() call necessary for actual progress report 
          res
        }
        RESU <- try(
          parallel::mclapply(iterator, FUN = p_FUN, ...)
        )
      })
    } else {
      .warn_once_progressr()
      RESU <- try(
        parallel::mclapply(iterator, FUN = FUN, ...)
      )
    }
  } else { # PSOCK
    cl <- cores_info$cl
    packages <- "Infusion"
    parallel::clusterExport(cl, "packages",envir=environment()) ## passes the list of packages to load
    abyss <- parallel::clusterEvalQ(cl, {sapply(packages,library,character.only=TRUE)}) ## snif
    if (cores_info$has_doSNOW) {
      show_pb <- (# verbose$most && 
        ! isTRUE(getOption('knitr.in.progress')))
      if (show_pb) {
        pb <- txtProgressBar(max = length(iterator), style = 3, char="P")
        progress <- function(n) setTxtProgressBar(pb, n)
        parallel::clusterExport(cl=cl, "progress",envir=environment()) ## slow?
        .options.snow = list(progress = progress)
      } else .options.snow = NULL
      it <- NULL ## otherwise R CMD check complains that no visible binding for global variable 'st'
      foreach_args <- list(
        it = iterator, 
        .packages= packages,
        .options.snow = .options.snow,
        .inorder = TRUE, .errorhandling = "remove"
        #                                 "pass"## "pass" to see error messages
      )
      foreach_blob <- do.call(foreach::foreach,foreach_args)
      RESU <- foreach::`%dopar%`(
        foreach_blob,
        FUN(it, ...) )
      if (show_pb) close(pb)
    } else { # PSOCK without doSNOW
      pbopt <- pboptions(nout=min(100L,2L*length(iterator)),type="timer", char="p")
      RESU <- pblapply(X=iterator, FUN = FUN, cl= cl, ...)
      pboptions(pbopt)
    }
  }
  .close_cores(cores_info)
  RESU
} 



