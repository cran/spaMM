.diagnose_bootreps_pbs <- function(bootreps, foreach_args) {
  if (foreach_args[[".errorhandling"]]=="remove" && is.null(bootreps)) {
    cat(crayon::bold(paste0(
      "Hmmm. It looks like all parallel processes failed. Maybe rerun spaMM_boot() \n",
      "with  ' control.foreach=list(.errorhandling=\"stop\") '  to diagnose the problem.\n"
    )))
  } else if (foreach_args[[".errorhandling"]]=="stop" && inherits(bootreps,"try-error")) {            
    # foreach alters the condition message => seel '\"' after 'could not find'
    if (length(grep("could not find",(condmess <- conditionMessage(attr(bootreps,"condition")))))) {
      firstpb <- strsplit(strsplit(condmess,"could not find")[[1]][2],"\"")[[1]][2]
      cat(crayon::bold(paste0(
        "Hmmm. It looks like some variables were not passed to the parallel processes.\n",
        "Maybe add   ' ",firstpb," = ",firstpb," '  to spaMM_boot()'s 'fit_env' argument?\n"
      )))
    } else cat(crayon::bold(condmess))
  }
}

.wrap_register_doFuture <- function(cl, iseed, nb_cores, PSOCK) {
  # .find_socket_backend() succesfully requireN...() the required packages.
  # doFuture::registerDoFuture() # before any %dopar%, becomes:
  rdS_fn <- get("registerDoFuture", asNamespace("doFuture"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
  do.call(rdS_fn,list()) # this is what makes foreach see it and perform parallel computations
  # plan_ <- get("plan", asNamespace("future"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
  if (PSOCK) {
    # cluster_ <- get("cluster", asNamespace("future"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
    future::plan(future::cluster, workers=cl)
  } else {
    # multicore_ <- get("multicore", asNamespace("future"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
    future::plan(future::multicore)
  }
  if ( ! is.null(iseed) ) parallel::clusterSetRNGStream(cl = cl, iseed) 
  #
  # A first foreach_blob for a first dopar before defining the progress bar (otherwise we see a progress bar on this dopar)
  i <- NULL ## otherwise R CMD check complains that no visible binding for global variable 'i' (in expression newy_s[,i])
  foreach_blob <- foreach::foreach(i=1:nb_cores)
  if (PSOCK) {
    abyss <- foreach::`%dopar%`(foreach_blob, Sys.setenv(LANG = "en")) # before setting the progress bar...
  }
}

.wrap_registerDoSNOW <- function(cl, iseed, nb_cores) {
  # loading (?) the namespace of 'snow' changes the *parent* RNG state (as well as sons' ones)! so we save and restore it 
  R.seed <- get(".Random.seed", envir = .GlobalEnv) # save parent RNG state
  rdS_fn <- get("registerDoSNOW", asNamespace("doSNOW"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
  do.call(rdS_fn,list(cl=cl)) # this is what makes foreach see it and perform parallel computations
  assign(".Random.seed", R.seed, envir = .GlobalEnv) # restore parent RNG state
  if ( ! is.null(iseed) ) parallel::clusterSetRNGStream(cl = cl, iseed) 
  #
  # A first foreach_blob for a first dopar before defining the progress bar (otherwise we see a progress bar on this dopar)
  i <- NULL ## otherwise R CMD check complains that no visible binding for global variable 'i' (in expression newy_s[,i])
  foreach_blob <- foreach::foreach(i=1:nb_cores)
  #if (cluster_args$type == "PSOCK") {
  abyss <- foreach::`%dopar%`(foreach_blob, Sys.setenv(LANG = "en")) # before setting the progress bar...
}


.wrap_registerDoParallel <- function(cl, iseed) {
  rdS_fn <- get("registerDoParallel", asNamespace("doParallel"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
  do.call(rdS_fn,list(cl=cl)) # this is what makes foreach see it and perform parallel computations
  if ( ! is.null(iseed) ) parallel::clusterSetRNGStream(cl = cl, iseed) 
  
  parallel::clusterCall(cl, Sys.setenv, LANG = "en")
}

.foreach_try_progressr <- function(newresp, fn, control, cluster_args, ...) {
  # some example at     # https://stackoverflow.com/questions/75252629/error-when-making-a-fork-cluster-and-registerdosnow-in-r
  nsim <- ncol(newresp)

  i <- NULL ## otherwise R CMD check complains that no visible binding for global variable 'i' (in expression newy_s[,i])
  foreach_args <- list( 
    i = 1:nsim, 
    .combine = "cbind", # may be overwritten by control$.combine
    .inorder = TRUE, 
    # .packages = "spaMM", # added by default in control.packages 
    .errorhandling = "remove", ## use "pass" to see problems
    .options.future = list(scheduling = 1.0)
  )
  foreach_args[names(control)] <- control # replaces the above defaults by user controls
  foreach_blob <- do.call(foreach::foreach,foreach_args) 
  
  fn_dots <- list(...)
  for (st in names(fn_dots)) {
    # Add an enclosing quote():
    if ( is.language(fn_dots[[st]])) fn_dots[[st]] <- substitute(quote(what),list(what=fn_dots[[st]]))
  }

  if ( do.call("requireNamespace",list(package="progressr", quietly = TRUE))) {
    # progressor is the only progress function that 'works' with mclapply
    # although not with load-balancing (mc.preschedule=FALSE)
    # Here we use the default (no balancing), and it is the steps with max value shown below that are reported.  
    progressor_ <- get("progressor", asNamespace("progressr"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
    with_progress_ <- get("with_progress", asNamespace("progressr"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
    if (cluster_args$spec>1L) {
      if (cluster_args$type=="FORK") {
        pb_char <- "F" # (no particular backend needed)
      } else pb_char <- "p" # doFuture backend case 
    } else pb_char <- "S" # presumably fictitious case
    
    handlers_ <- get("handlers", asNamespace("progressr"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
    handlers_(global=TRUE) # hmf. passes global variables ...? (try to remove it?)
    handler_txtprogressbar_ <- get("handler_txtprogressbar", asNamespace("progressr"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
    # define the progress bar:
    barstyle <- eval(spaMM.getOption("barstyle"))
    # (1) do not duplicate progressor(steps) by setting 'max'; (2) clear =FALSE otherwise pb disappears immediately
    handlers_(handler_txtprogressbar_(style = barstyle, char=pb_char, clear=FALSE)) 
    
    with_progress_({
      p <- progressor_(steps = nsim) 
      bootreps <- try(foreach::`%dopar%`(foreach_blob, {
        ## Fail to capture the package loading messages (<=> in the foreach but outside this expression) 
        # prom <- future::future(do.call(fn, c(list(newresp[, i]), fn_dots)), stdout=FALSE)
        # res <- future::value(prom)
        res <- do.call(fn, c(list(newresp[, i]), fn_dots))
        p() # p(sprintf("i=%g", i)) might be nicer but the barstyle seems to override this
        res
      }))
    })
    
  } else {
    .warn_once_progressr()
    bootreps <- try(foreach::`%dopar%`(foreach_blob, {do.call(fn, c(list(newresp[, i]), fn_dots))}))
  }
  .diagnose_bootreps_pbs(bootreps, foreach_args)
  bootreps
}



.foreach_PSOCK_nofuture <- function(newresp=newresp, fn=fn, control=control, ...) { # either doSNOW or doParallel backend
  
  nsim <- ncol(newresp)

  i <- NULL ## otherwise R CMD check complains that no visible binding for global variable 'i' (in expression newy_s[,i])
  foreach_args <- list( 
    i = 1:nsim, 
    .combine = "cbind",  # may be overwritten by control$.combine
    .inorder = TRUE, 
    # .packages = "spaMM", # added by default in control.packages 
    .errorhandling = "remove" ## use "pass" to see problems
  )
  foreach_args[names(control)] <- control # replaces the above defaults by user controls
  foreach_blob <- do.call(foreach::foreach,foreach_args) 
  
  fn_dots <- list(...)
  for (st in names(fn_dots)) {
    # Add an enclosing quote():
    if ( is.language(fn_dots[[st]])) fn_dots[[st]] <- substitute(quote(what),list(what=fn_dots[[st]]))
  }
  bootreps <- try(foreach::`%dopar%`(foreach_blob, {do.call(fn, c(list(newresp[, i]), fn_dots))}))
  
  .diagnose_bootreps_pbs(bootreps, foreach_args)
  bootreps
}

.foreach_snow_bar <- function(newresp, fn, control, ...) {
  
  # define the progress bar:
  progrbar_setup <- .set_progrbar(max = ncol(newresp), style = eval(spaMM.getOption("barstyle")), char="P") # pb_char
  control$.options.snow <- progrbar_setup["progress"]
  
  bootreps <- .foreach_PSOCK_nofuture(newresp=newresp, fn=fn, control=control, ...)
  
  close(progrbar_setup$pb)
  bootreps
}

.foreach_serial_bar <- function(newresp, fn, control, ...) {
  nsim <- ncol(newresp)

  # :where opts are needed to define a second foreach_blob
  i <- NULL
  foreach_args <- list( 
    i = 1:nsim, 
    .combine = "cbind",  # may be overwritten by control$.combine
    .inorder = TRUE, 
    # .packages = "spaMM", # added by default in control.packages 
    .errorhandling = "remove" ## use "pass" to see problems
  )
  
  foreach_args[names(control)] <- control # replaces the above defaults by user controls
  
  # The serial code is distinct in wrapping the .combine function with progress bar code, instead of using progressr
  .combine <- foreach_args$.combine
  if (inherits(.combine,"character")) .combine <- get(.combine)
  pb <- txtProgressBar(max = nsim, style = eval(spaMM.getOption("barstyle")), char="s") # pb_char
  combine_with_pb <- function() {
    count <-0L
    function(...) {
      count <<- count + length(list(...))
      setTxtProgressBar(pb, count)
      .combine(...) 
    }
  }
  foreach_args$.combine <- combine_with_pb()
  
  foreach_blob <- do.call(foreach::foreach,foreach_args) 
  
  fn_dots <- list(...)
  bootreps <- try(foreach::`%do%`(foreach_blob, do.call(fn, c(list(newresp[, i]), fn_dots)))) 
  # the try() is useful if the user interrupts the %do%, in which case it allows close(pb) to be run.
  
  close(pb)
  bootreps
}


.find_socket_backend <- local({
  nobar_warned <- FALSE
  nofuture_warned <- FALSE
  function() {
    if (suppressWarnings(do.call("requireNamespace",list(package="doSNOW", quietly = TRUE)))) return("doSNOW")
    #
    if (suppressWarnings(do.call("requireNamespace",list(package="doFuture", quietly = TRUE)))) {
      if ( ! suppressWarnings(do.call("requireNamespace",list(package="future", quietly = TRUE)))) {
        if ( ! nofuture_warned) {
          message("'doFuture' installed but not 'future': trying 'doParallel'.")
          nofuture_warned <<- TRUE
        } 
      } else if ( ! suppressWarnings(do.call("requireNamespace",list(package="progressr", quietly = TRUE)))) {
        if ( ! nofuture_warned) {
          message("'doFuture' installed but not 'progressr': trying 'doParallel'.")
          nofuture_warned <<- TRUE
        } 
      } else return("doFuture")
    }
    #
    if ( ! nobar_warned) {
      message("Neither 'doSNOW' nor 'doFuture' (& co.) installed: no progress bar.")
      nobar_warned <<- TRUE
    } 
    if ( ! suppressWarnings(do.call("requireNamespace",list(package="doParallel", quietly = TRUE)))) {
      stop("You *need* to install some parallel backend (doSNOW | doFuture & co. | doParallel) to use a cluster with foreach().")
    } else return("doParallel")
  }
})

.setCluster <- function(nb_cores, cluster_args, iseed, fit_env=NULL) {
  
  cluster_args <- .set_cluster_type(cluster_args, nb_cores=cluster_args$spec) # If I extract this call from the .setCluster() 
  #   I must make sure that spaMM:::.set_cluster_type is accessible (VERSUS: here in the EXPORTED spaMM:::.setCluster())
  nb_cores <- cluster_args$spec
  
  if (nb_cores>1L) {
    if (cluster_args$type=="FORK") {
      cl <- parallel::makeForkCluster(nnodes = nb_cores)
      .wrap_register_doFuture(cl, iseed=iseed, nb_cores=nb_cores, PSOCK=FALSE)
    } else { # PSOCK
      cl <- do.call(parallel::makeCluster, cluster_args) # note that _this_ line would make sense for fork clusters too. BUT
      # ... the foreach = dot args combination may not work for FORK type. Only pbapply would work with makeCluster+FORK, 
      # but pbmcapply is a better way to get a pb one a fork cluster as [pb]mclapply have better load balancing than pbapply. 
      if (is.environment(fit_env)) parallel::clusterExport(cl=cl, varlist=ls(fit_env), envir=fit_env)
      backend <- .find_socket_backend()
      if (backend=="doSNOW") {
        .wrap_registerDoSNOW(cl, iseed, nb_cores)
      } else if (backend=="doFuture") { # implies that .find_socket_backend() succesfully requireN...() the required packages.
        .wrap_register_doFuture(cl, iseed=iseed, nb_cores=nb_cores, PSOCK=TRUE)
      } else { # neither doSNOW nor doFuture => doParallel but not bar
        .wrap_registerDoParallel(cl, iseed)
      } # has_doSNOW ... else has_doFuture ... else ...
    } # FORK ... else
    cl
  } else NULL
}


combinepar <- function(newresp, fn, nb_cores=NULL, cluster=NULL,
           fit_env, control=list(), cluster_args=NULL,
           debug.=FALSE, iseed=NULL, showpbar=eval(spaMM.getOption("barstyle")),
           pretest_cores=NULL, 
           ... # passed to fn... unless captured by pbapply (in which case 'simplify' may have a distinct effect)
) {
  if (is.list(fit_env)) fit_env <- list2env(fit_env)
  if ( ! "spaMM" %in% control$.packages) control$.packages <- c("spaMM", control$.packages) 

  if (cluster_is_local <- is.null(cluster)) {
    cluster_args <- .set_cluster_type(cluster_args, nb_cores) # PSOCK vs FORK
    nb_cores <- cluster_args$spec
  } else if (inherits(cluster,"cluster")) {
    nb_cores <- length(cluster)
    if (inherits(cluster[[1]],"SOCKnode")) {
      type <- "PSOCK"
    } else type <- "FORK" #note that 'cluster' itself inherits from type "SOCKcluster" also in that case (!)
    cluster_args <- list(spec=nb_cores, type=type)
    cl <- cluster
  } else stop("Unhandled type of 'cluster'")
  
  if (debug. && nb_cores>1L ) debug. <- 1L 
  assign("debug.", debug., environment(fn))
  if (is.null(dim(newresp))) newresp <- matrix(seq(newresp),ncol=newresp,nrow=1) # assuming newresp is an integer
  nsim <- ncol(newresp)
  if (nb_cores>1L) {
    if ( ! is.null(iseed) ) {
      ori <- RNGkind("L'Ecuyer-CMRG")
      set.seed(iseed)
    }
    if (cluster_args$type=="FORK") {
      if (cluster_is_local) {
        cl <- parallel::makeForkCluster(nnodes = nb_cores)
        .wrap_register_doFuture(cl, iseed=iseed, nb_cores=nb_cores, PSOCK=FALSE)
      } 
      bootreps <- .foreach_try_progressr(newresp=newresp, fn=fn, control=control, cluster_args=cluster_args, ...)
    } else { # PSOCK
      if (cluster_is_local) {
        cl <- do.call(parallel::makeCluster, cluster_args) # note that _this_ line would make sense for fork clusters too. BUT
        # ... the foreach = dot args combination may not work for FORK type. Only pbapply would work with makeCluster+FORK, 
        # but pbmcapply is a better way to get a pb one a fork cluster as [pb]mclapply have better load balancing than pbapply. 
        if (is.environment(fit_env)) parallel::clusterExport(cl=cl, varlist=ls(fit_env), envir=fit_env)
      } 
      backend <- .find_socket_backend()
      if (backend=="doSNOW") {
        if (cluster_is_local) .wrap_registerDoSNOW(cl, iseed, nb_cores)
        bootreps <- .foreach_snow_bar(newresp=newresp, fn=fn, control=control, ...) # wraps .foreach_PSOCK_nofuture() with the added bar
      } else if (backend=="doFuture") { 
        if (cluster_is_local) .wrap_register_doFuture(cl, iseed=iseed, nb_cores=nb_cores, PSOCK=TRUE)
        bootreps <- .foreach_try_progressr(newresp=newresp, fn=fn, control=control, cluster_args=cluster_args, ...)
      } else { # neither doSNOW nor doFuture => doParallel but not bar
        if (cluster_is_local) .wrap_registerDoParallel(cl, iseed)
        bootreps <- .foreach_PSOCK_nofuture(newresp=newresp, fn=fn, control=control, ...)
      } # has_doSNOW ... else has_doFuture ... else ...
      if (inherits(bootreps,"try-error") ) {
        if (length(grep("could not find",(condmess <- conditionMessage(attr(bootreps,"condition")))))) {
          firstpb <- strsplit(condmess,"\"")[[1]][2]
          cat(crayon::bold(paste0(
            "Hmmm. It looks like some variables were not passed to the parallel processes.\n",
            "Maybe add    ",firstpb," = ",firstpb,"   to spaMM_boot()'s 'fit_env' argument?\n"
          )))
        } else cat(crayon::bold(condmess))
      }
    } # FORK ... else
    if (cluster_is_local) {
      foreach::registerDoSEQ() ## https://stackoverflow.com/questions/25097729/un-register-a-doparallel-cluster    [not necessarily DoParallel]
      parallel::stopCluster(cl)
    }
    if ( ! is.null(iseed) ) do.call("RNGkind", as.list(ori)) # restore to state pre-parallel computation 
    # (makes sense if not preset cluster. If preset cluster, it may make more sense to control RNG once when creating and once when closing it)
  } else {
    bootreps <-  .foreach_serial_bar(newresp, fn, control, ...)
  }
  return(bootreps)
}
