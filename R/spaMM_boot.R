spaMM_boot <- local({
  doSNOW_warned <- FALSE
  function(object, simuland, nsim, nb_cores=NULL,
           resp_testfn=NULL, control.foreach=list(),
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
      if (has_doSNOW <- ("package:doSNOW" %in% search())) { ## allows progressbar but then requires foreach
        assign(".Random.seed", R.seed, envir = .GlobalEnv)
        fn <- get("registerDoSNOW", asNamespace("doSNOW"))
        do.call(fn,list(cl=cl)) 
        #`%foreachdopar%` <- foreach::`%dopar%`
        pb <- txtProgressBar(max = nsim, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        parallel::clusterExport(cl=cl, list("progress"),envir=environment()) ## slow! why?
      } else {
        if ( ! doSNOW_warned) {
          message("If the 'doSNOW' package were attached, the progress of the bootstrap computation could be reported.")
          doSNOW_warned <<- TRUE
        } 
        parallel::clusterEvalQ(cl, library("spaMM")) ## for pbapply
      }
      #dotenv <- list2env(list(...))
      #parallel::clusterExport(cl=cl, as.list(ls(dotenv)),envir=dotenv) ## much faster...
    } else cl <- NULL
    ####
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    newy_s <- simulate(object,nsim = nsim,verbose=FALSE,resp_testfn=resp_testfn) 
    if (nsim==1L) dim(newy_s) <- c(length(newy_s),1)
    if (nb_cores > 1L && has_doSNOW) {
      #    dotenv <- list2env(list(...))
      #  parallel::clusterExport(cl=cl, as.list(ls(dotenv)),envir=dotenv) ## much faster...
      dots <- list(...)
      i <- NULL ## otherwise R CMD check complains that no visible binding for global variable 'i'
      foreach_args <- list(
        i = 1:nsim, .combine = "rbind",
        .inorder = TRUE, .packages = "spaMM", .errorhandling = "remove", 
        #"pass", ## "pass" to see error messages
        .options.snow = opts
      )
      foreach_args[names(control.foreach)] <- control.foreach
      foreach_blob <- do.call(foreach::foreach,foreach_args)
      bootreps <- foreach::`%dopar%`(foreach_blob,
                                     do.call(simuland,c(list(y=newy_s[,i]), dots)))
      #bootreps <- sapply(bootreps,identity)
      close(pb)
    } else {
      pbopt <- pboptions(nout=min(100,2*nsim),type="timer")
      bootreps <- pbapply(X=newy_s,MARGIN = 2L,FUN = simuland, cl=cl, ...) 
      pboptions(pbopt)
    }
    ####
    if (nb_cores > 1L) { parallel::stopCluster(cl) } 
    cat(paste(" bootstrap took",.timerraw(time1),"s.\n")) 
    return(list(bootreps=bootreps,RNGstate=RNGstate))
  }
})

