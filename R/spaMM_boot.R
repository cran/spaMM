spaMM_boot <- local({
  doSNOW_warned <- FALSE
  function(object, simuland, nsim, nb_cores=NULL,
           resp_testfn=NULL, control.foreach=list(),
           debug. = FALSE, type="marginal",
           ...) {
    msg <- "Bootstrap replicates:"
    msglength <- nchar(msg) ## not used in the parallel case
    cat(msg)
    time1 <- Sys.time() 
    cumul_nsim <- 0L
    RNGstateList <- vector("list")
    nsim <- as.integer(nsim) 
    if (nsim<1L) {
      warning("'nsim' must be at least 1.")
      return(list())
    }
    nb_cores <- .check_nb_cores(nb_cores=nb_cores) 
    if (nb_cores > 1L) {
      #cl <- parallel::makeCluster(nb_cores,outfile="essai.txt") 
      cl <- parallel::makeCluster(nb_cores) 
      R.seed <- get(".Random.seed", envir = .GlobalEnv)
      has_doSNOW <- ("package:doSNOW" %in% search())
      if (has_doSNOW) { ## allows progressbar but then requires foreach
        assign(".Random.seed", R.seed, envir = .GlobalEnv) # loading (?) the namespace of 'snow' changes the global RNG state!
        fn <- get("registerDoSNOW", asNamespace("doSNOW"))
        do.call(fn,list(cl=cl)) 
      } else {
        if ( ! doSNOW_warned) {
          message("If the 'doSNOW' package were attached, better load-balancing might be possible.")
          doSNOW_warned <<- TRUE
        } 
        parallel::clusterEvalQ(cl, library("spaMM")) ## for pbapply
      }
    } else cl <- NULL
    ####
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    newy_s <- simulate(object,nsim = nsim,verbose=FALSE,resp_testfn=resp_testfn, type=type) 
    dots <- list(...) ## We will use a list, not a call, for exports etc.
    nfs <- names(formals(simuland))
    simuland_dots <- intersect(names(dots),nfs)
    simuland_dots <- dots[simuland_dots]
    for (st in names(simuland_dots)) {
      # Add an enclusing quote():
      if ( is.language(simuland_dots[[st]])) simuland_dots[[st]] <- substitute(quote(what),list(what=simuland_dots[[st]]))
    }
    if (nsim==1L) dim(newy_s) <- c(length(newy_s),1)
    if (nb_cores > 1L && has_doSNOW) {
      other_dots <- setdiff(names(dots),nfs)
      other_dots <- list2env(dots[other_dots])
      parallel::clusterExport(cl=cl, as.list(ls(other_dots)),envir=other_dots) 
      foreach_blob <- foreach::foreach(i=1:nb_cores)
      abyss <- foreach::`%dopar%`(foreach_blob, Sys.setenv(LANG = "en")) # before setting the progress bar...
      pb <- txtProgressBar(max = nsim, style = 3, char="P")
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      # parallel::clusterExport(cl=cl, list("progress"),envir=environment()) ## slow! why? # and not useful ??
      i <- NULL ## otherwise R CMD check complains that no visible binding for global variable 'i' (in expression newy_s[,i])
      if (debug.) {
        # One may want to redefine simuland so that it can return anything, 
        # but if it is redefined in debugging session, it loses its original environment, in which case we need:
        if (debug_interactive_redef_simuland <- FALSE) { # block useful in debuguing session
          simfun <- eval(match.call()$simuland,-2)
          #ls(environment(simfun))
          environment(simuland) <- environment(simfun)
        }
        foreach_args <- list(
          i = 1:nsim, 
          .inorder = TRUE, .packages = "spaMM", .errorhandling = "pass", ## pass error messages
          .options.snow = opts
        )
      } else {
        foreach_args <- list(
          i = 1:nsim, 
          .combine = "rbind", ## rbinds 1-row data frames (into a data frame) as well as vectors (into a matrix) (result has nsim rows)
          .inorder = TRUE, .packages = "spaMM", .errorhandling = "remove", 
          .options.snow = opts
        )
      }
      foreach_args[names(control.foreach)] <- control.foreach
      foreach_blob <- do.call(foreach::foreach,foreach_args)
      bootreps <- foreach::`%dopar%`(foreach_blob,
                                     do.call(simuland,c(list(y=newy_s[,i]),simuland_dots)))
      close(pb)
    } else { ## serial or parallel pbapply
      if (nb_cores>1L) {pb_char <- "p"} else pb_char <- "s"
      pbopt <- pboptions(nout=min(100,2*nsim),type="timer",char=pb_char) # I tested ,use_lb=TRUE (again) on 
      #                              spaMM v2.6.75, pbapply 1.4.0 (4/2019) and this has no measurable effect.
      bootreps <- do.call("pbapply",c(list(X=newy_s,MARGIN = 2L,FUN = simuland, cl=cl), simuland_dots)) 
      if ( ! is.list(bootreps)) {
        # tries to return an nsim-rows matrix
        if ( ! is.matrix(bootreps)) {
          dim(bootreps) <- c(length(bootreps),1L)
        } else bootreps <- t(bootreps)
      } else if (is.list(bootreps) && (! debug.) && identical(unique(unlist(lapply(bootreps,nrow))),1L)) {
        bootreps <- do.call(rbind, bootreps) # rbinds 1-row data frames (into a data frame) to match the doSNOW-based result.
      } # in debug. case bootreps is a list not suitable for rbind()ing
      pboptions(pbopt)
    }
    ####
    if (nb_cores > 1L) {
      if (has_doSNOW) foreach::registerDoSEQ() ## https://stackoverflow.com/questions/25097729/un-register-a-doparallel-cluster
      parallel::stopCluster(cl) 
    } 
    cat(paste(" bootstrap took",.timerraw(time1),"s.\n")) 
    return(list(bootreps=bootreps,RNGstates=RNGstate))
  }
})

