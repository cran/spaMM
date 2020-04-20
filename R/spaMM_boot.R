spaMM_boot <- local({
  doSNOW_warned <- FALSE
  function(object, simuland, nsim, nb_cores=NULL,
           resp_testfn=NULL, control.foreach=list(),
           debug. = FALSE, 
           type,
           fit_env=NULL,
           cluster_args=NULL,
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
    newy_s <- simulate(object,nsim = nsim,verbose=FALSE,resp_testfn=resp_testfn, type=type) 
    if (nsim==1L) dim(newy_s) <- c(length(newy_s),1L)
    #
    # If the simuland has (say) arguments y, what=NULL, lrt, ...   , we should not have lrt in the dots. Since the dots are not directly manipulable
    # we have to convert them to a list, and ultimately to use do.call()
    control.foreach$.combine <- "rbind"
    bootreps <- dopar(newresp = newy_s, nb_cores = nb_cores,fn = simuland, fit_env = fit_env, 
                       control=control.foreach, debug.=debug., ...) 
    return(list(bootreps=bootreps,RNGstates=RNGstate))
  }
})

# fn more generic than spaMM_boot: there is no call to other spaMM fns such as simulate(object, .) so this acts as a general wrapper for 
# foreach or pbapply, and not specifically for bootstrap computations.
dopar <- local({
  doSNOW_warned <- FALSE
  function(newresp, fn, nb_cores=NULL, 
           fit_env, control=list(), cluster_args=NULL,
           debug.=FALSE, iseed=NULL,
           ... # passed to fn
  ) {
    if (is.list(fit_env)) fit_env <- list2env(fit_env)
    nb_cores <- .check_nb_cores(nb_cores=nb_cores) 
    if (debug. && nb_cores>1L ) debug. <- 1L 
    assign("debug.", debug., environment(fn))
    nsim <- ncol(newresp)
    time1 <- Sys.time() 
    if (nb_cores>1L) {
      cl <- do.call(parallel::makeCluster, c(list(spec=nb_cores), cluster_args)) 
      has_doSNOW <- ("package:doSNOW" %in% search())
      if (has_doSNOW) {
        # loading (?) the namespace of 'snow' changes the *parent* RNG state (as well as sons' ones)! so we save and restore it 
        R.seed <- get(".Random.seed", envir = .GlobalEnv) # save parent RNG state
        rdS_fn <- get("registerDoSNOW", asNamespace("doSNOW")) # syntax for using an undeclared package...
        do.call(rdS_fn,list(cl=cl)) # this is what makes foreach see it and perform parallel computations
        assign(".Random.seed", R.seed, envir = .GlobalEnv) # restore parent RNG state
        if ( ! is.null(iseed) ) parallel::clusterSetRNGStream(cl = cl, iseed) 
        #
        if (is.environment(fit_env)) parallel::clusterExport(cl=cl, varlist=ls(fit_env), envir=fit_env) 
        # A first foreach_blob for a first dopar before defining the progress bar (otherwise we see a progress bar on this dopar)
        i <- NULL ## otherwise R CMD check complains that no visible binding for global variable 'i' (in expression newy_s[,i])
        foreach_blob <- foreach::foreach(i=1:nb_cores)
        abyss <- foreach::`%dopar%`(foreach_blob, Sys.setenv(LANG = "en")) # before setting the progress bar...
        # define the progress bar:
        pb <- txtProgressBar(max = nsim, style = 3, char="P")
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        # :where opts are needed to define a second foreach_blob
        foreach_args <- list( 
          i = 1:ncol(newresp), 
          .combine = "cbind", 
          .inorder = TRUE, .packages = "spaMM", 
          .errorhandling = "remove", ## use "pass" to see problems
          .options.snow = opts
        )
        foreach_args[names(control)] <- control # replaces the above defaults by user controls
        foreach_blob <- do.call(foreach::foreach,foreach_args) ## rbinds 1-row data frames (into a data frame) as well as vectors (into a matrix) (result has nsim rows)
        if (TRUE) {
          fn_dots <- list(...)
          for (st in names(fn_dots)) {
            # Add an enclusing quote():
            if ( is.language(fn_dots[[st]])) fn_dots[[st]] <- substitute(quote(what),list(what=fn_dots[[st]]))
          }
          bootreps <- try(foreach::`%dopar%`(foreach_blob, do.call(fn, c(list(newresp[, i]), fn_dots)))) 
        } else {
          # Standard passing of the dots with foreach does not seem to work. (good test is the doSNOW case nested within test-LRT-boot.R)
          #bootreps <- try(foreach::`%dopar%`(foreach_blob, fn(newresp[, i], ...)))
        }
        # the try() is useful if the user interrupts the dopar, in which case it allows close(pb) to be run. (? But doSNOW appear to close the nodes asynchronously?)
        foreach::registerDoSEQ() ## https://stackoverflow.com/questions/25097729/un-register-a-doparallel-cluster
        parallel::stopCluster(cl)
        close(pb)
      } else {
        # in that case, ## We will use pbapply, with argument cl=cl; a direct call to foreach would require doParallel::registerDoParallel(cl)
        if ( ! doSNOW_warned) {
          message("If the 'doSNOW' package were attached, better load-balancing might be possible.")
          doSNOW_warned <<- TRUE
        } 
        pb_char <- "p"
        parallel::clusterEvalQ(cl, library("spaMM")) 
        if (is.environment(fit_env)) try(parallel::clusterExport(cl=cl, varlist=ls(fit_env), envir=fit_env)) 
        pbopt <- pboptions(nout=min(100L,2L*nsim),type="timer",char=pb_char) 
        #try() so that an interrupt does not prevent running stopCluster():
        bootreps <- try(pbapply(X=newresp,MARGIN = 2L,FUN = fn, cl=cl, ...))
        parallel::stopCluster(cl)
        pboptions(pbopt)
        if (identical(control$.combine,"rbind")) bootreps <- t(bootreps)
      }
    } else { ## nb_cores=1L
      pb_char <- "s"
      pbopt <- pboptions(nout=min(100L,2L*nsim),type="timer",char=pb_char) 
      bootreps <- pbapply(X=newresp,MARGIN = 2L,FUN = fn, cl=NULL, ...)
      pboptions(pbopt)
      if (identical(control$.combine,"rbind")) bootreps <- t(bootreps)
    }
    cat(paste(" bootstrap took",.timerraw(time1),"s.\n")) 
    return(bootreps)
  }
})
