spaMM_boot <- local({
  doSNOW_warned <- FALSE
  function(object, simuland, nsim, nb_cores=NULL,
           seed=NULL,
           resp_testfn=NULL, control.foreach=list(),
           debug. = FALSE, 
           type,
           fit_env=NULL,
           cluster_args=NULL,
           showpbar= eval(spaMM.getOption("barstyle")),
           ...) {
    if (missing(type)) {
      warning("'type' is now a mandatory argument of spaMM_boot().\n Assuming type='marginal' for consistency with previous versions.",
              immediate. = TRUE)
      type <- "marginal"
    }
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    msg <- "Bootstrap replicates:"
    msglength <- nchar(msg) ## not used in the parallel case
    cat(msg)
    cumul_nsim <- 0L
    nsim <- as.integer(nsim) 
    if (nsim<1L) {
      warning("'nsim' must be at least 1.")
      return(list())
    }
    ####
    newy_s <- simulate(object,nsim = nsim,verbose=c(type=FALSE,showpbar=showpbar), resp_testfn=resp_testfn, type=type, seed=seed) 
    if (nsim==1L) dim(newy_s) <- c(length(newy_s),1L)
    #
    # If the simuland has (say) arguments y, what=NULL, lrt, ...   , we should not have lrt in the dots. Since the dots are not directly manipulable
    # we have to convert them to a list, and ultimately to use do.call()
    control.foreach$.combine <- "rbind"
    wrap_parallel <- get(.spaMM.data$options$wrap_parallel, asNamespace("spaMM"), inherits=FALSE)
    bootreps <- wrap_parallel(newresp = newy_s, nb_cores = nb_cores,fn = simuland, fit_env = fit_env, 
                       control=control.foreach, debug.=debug., pretest_cores = .pretest_fn_on_cores, 
                      showpbar = showpbar, ...) 
    return(list(bootreps=bootreps,RNGstates=RNGstate))
  }
})

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
  if ( ! is.null(nullfit)) errors_null <- pbapply::pbreplicate(1, .check_call_on_core(fitobject = nullfit), cl=cluster)
  #foreach_blob <- foreach::foreach(i=1)
  #if ( ! is.null(nullfit)) errors_null <- foreach::`%dopar%`(foreach_blob, .check_call_on_core(fitobject = nullfit))
  fullfit <- environment(fn)$fullfit
  if ( ! is.null(fullfit)) errors_full <- pbapply::pbreplicate(1, .check_call_on_core(fitobject = fullfit), cl=cluster)
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
        if (.inRstudio(silent=TRUE)) {
          cluster_args$type <- "PSOCK"
        } else cluster_args$type <- "PSOCK" # default is socket cluster in all cases
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
           ... # passed to fn
  ) {
    if (is.list(fit_env)) fit_env <- list2env(fit_env)
    cluster_args <- .set_cluster_type(cluster_args, nb_cores) # PSOCK vs FORK
    nb_cores <- cluster_args$spec
    if (debug. && nb_cores>1L ) debug. <- 1L 
    assign("debug.", debug., environment(fn))
    nsim <- ncol(newresp)
    time1 <- Sys.time() 
    if (nb_cores>1L) {
      if (cluster_args$type=="FORK") {
        if ( ! is.null(iseed) ) {
          ori <- RNGkind("L'Ecuyer-CMRG")
          set.seed(iseed)
        }
        if (is.null(mc.silent <- control$mc.silent)) mc.silent <- TRUE 
        has_progressr <- ("package:progressr" %in% search())
        seq_nr <- seq_len(ncol(newresp))
        if (has_progressr) {
          # progressor is the only progress function that 'works' with mclapply
          # although not with load-balancing (mc.preschedule=FALSE)
          # Here we use the default (no balancing), and it is the steps with max value shown below that are reported.  
          prog_fn <- get("progressor", asNamespace("progressr"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
          with_fn <- get("with_progress", asNamespace("progressr"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
          with_fn({
            p <- prog_fn(steps=ceiling(ncol(newresp)/nb_cores))
            p_fn <- function(it, ...) { # it OK for mclapply... not for apply on a matrix
              res <- fn(newresp[,it], ...)
              p() # p() call necessary for actual progress report 
              res
            }
            bootreps <- try(
              parallel::mclapply(seq_nr, FUN = p_fn, mc.silent=mc.silent, mc.cores=nb_cores)
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
        if ( ! is.null(iseed) ) do.call("RNGkind", as.list(ori)) # reste to state pre-parallel computation
      } else { # PSOCK
        cl <- do.call(parallel::makeCluster, cluster_args) # note that _this_ line would make sense for fork clusters too. BUT
        # ... the foreach = dot args combination may not work for FORK type. Only pbapply would work with makeCluster+FORK, 
        # but pbmcapply is a better way to get a pb one a fork cluster as [pb]mclapply have better load balancing than pbapply. 
        has_doSNOW <- ("package:doSNOW" %in% search())
        if (has_doSNOW) {
          # loading (?) the namespace of 'snow' changes the *parent* RNG state (as well as sons' ones)! so we save and restore it 
          R.seed <- get(".Random.seed", envir = .GlobalEnv) # save parent RNG state
          rdS_fn <- get("registerDoSNOW", asNamespace("doSNOW"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
          do.call(rdS_fn,list(cl=cl)) # this is what makes foreach see it and perform parallel computations
          assign(".Random.seed", R.seed, envir = .GlobalEnv) # restore parent RNG state
          if ( ! is.null(iseed) ) parallel::clusterSetRNGStream(cl = cl, iseed) 
          #
          if (cluster_args$type == "PSOCK") {
            if (is.environment(fit_env)) parallel::clusterExport(cl=cl, varlist=ls(fit_env), envir=fit_env) 
            pb_char <- "P"
          } else pb_char <- "F"
          # A first foreach_blob for a first dopar before defining the progress bar (otherwise we see a progress bar on this dopar)
          i <- NULL ## otherwise R CMD check complains that no visible binding for global variable 'i' (in expression newy_s[,i])
          foreach_blob <- foreach::foreach(i=1:nb_cores)
          if (cluster_args$type == "PSOCK") {
            abyss <- foreach::`%dopar%`(foreach_blob, Sys.setenv(LANG = "en")) # before setting the progress bar...
            if (is.function(pretest_cores)) pretest_cores(fn, cl)
          }
          # define the progress bar:
          barstyle <- eval(spaMM.getOption("barstyle"))
          progrbar_setup <- .set_progrbar(max = nsim, style = barstyle, char=pb_char)
          # :where opts are needed to define a second foreach_blob
          foreach_args <- list( 
            i = 1:ncol(newresp), 
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
              # Add an enclusing quote():
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
          if (showpbar) close(progrbar_setup$pb)
        } else { # no doSNOW
          # in that case, ## We will use pbapply, with argument cl=cl; a direct call to foreach would require doParallel::registerDoParallel(cl)
          if ( ! doSNOW_warned) {
            message("If the 'doSNOW' package were attached, better load-balancing might be possible.")
            doSNOW_warned <<- TRUE
          } 
          pb_char <- "p"
          if ( ! is.null(iseed) ) parallel::clusterSetRNGStream(cl = cl, iseed) 
          parallel::clusterEvalQ(cl, library("spaMM")) 
          if (is.environment(fit_env)) try(parallel::clusterExport(cl=cl, varlist=ls(fit_env), envir=fit_env)) 
          if (is.function(pretest_cores)) pretest_cores(fn, cl)
          if (showpbar) {
            pbopt <- pboptions(nout=min(100L,2L*nsim),type="timer",char=pb_char) 
          } else pbopt <- pboptions(type="none") 
          #try() so that an interrupt does not prevent running stopCluster():
          bootreps <- try(pbapply(X=newresp,MARGIN = 2L,FUN = fn, cl=cl, ...))
          parallel::stopCluster(cl)
          pboptions(pbopt)
          if (identical(control$.combine,"rbind")) bootreps <- t(bootreps)
        } # has_doSNOW ... else
      } # FORK ... else
    } else { ## nb_cores=1L
      pb_char <- "s"
      if (showpbar) {
        pbopt <- pboptions(nout=min(100L,2L*nsim),type="timer",char=pb_char) 
      } else pbopt <- pboptions(type="none") 
      bootreps <- pbapply(X=newresp,MARGIN = 2L,FUN = fn, cl=NULL, ...)
      pboptions(pbopt)
      if (identical(control$.combine,"rbind")) bootreps <- t(bootreps)
    }
    cat(paste(" bootstrap took",.timerraw(time1),"s.\n")) 
    return(bootreps)
  }
})
