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
        assign(".Random.seed", R.seed, envir = .GlobalEnv) # loading (?) the namespace of 'snow' changes the global RNG state!
        fn <- get("registerDoSNOW", asNamespace("doSNOW"))
        do.call(fn,list(cl=cl)) 
        #`%foreachdopar%` <- foreach::`%dopar%`
        pb <- txtProgressBar(max = nsim, style = 3, char="P")
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        parallel::clusterExport(cl=cl, list("progress"),envir=environment()) ## slow! why?
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
    newy_s <- simulate(object,nsim = nsim,verbose=FALSE,resp_testfn=resp_testfn) 
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
      i <- NULL ## otherwise R CMD check complains that no visible binding for global variable 'i'
      foreach_args <- list(
        i = 1:nsim, 
        .combine = "rbind", ## rbinds 1-row data frames (into a data frame) as well as vectors (into a matrix) (result has nsim rows)
        .inorder = TRUE, .packages = "spaMM", .errorhandling = "remove", 
        #"pass", ## "pass" to see error messages
        .options.snow = opts
      )
      foreach_args[names(control.foreach)] <- control.foreach
      foreach_blob <- do.call(foreach::foreach,foreach_args)
      bootreps <- foreach::`%dopar%`(foreach_blob,
                                     do.call(simuland,c(list(y=newy_s[,i]),simuland_dots)))
      #bootreps <- sapply(bootreps,identity)
      close(pb)
    } else { ## serial or parallel pbapply
      if (nb_cores>1L) {pb_char <- "p"} else pb_char <- "s"
      pbopt <- pboptions(nout=min(100,2*nsim),type="timer",char=pb_char)
      bootreps <- do.call("pbapply",c(list(X=newy_s,MARGIN = 2L,FUN = simuland, cl=cl), simuland_dots)) 
      if ( ! is.list(bootreps)) {
        # tries to return an nsim-rows matrix
        if ( ! is.matrix(bootreps)) {
          dim(bootreps) <- c(length(bootreps),1L)
        } else bootreps <- t(bootreps)
      } else if (is.list(bootreps) &&  identical(unique(unlist(lapply(bootreps,nrow))),1L)) {
        bootreps <- do.call(rbind, bootreps) # rbinds 1-row data frames (into a data frame) to match the doSNOW-based result.
      }
      pboptions(pbopt)
    }
    ####
    if (nb_cores > 1L) {
      if (has_doSNOW) foreach::registerDoSEQ() ## https://stackoverflow.com/questions/25097729/un-register-a-doparallel-cluster
      parallel::stopCluster(cl) 
    } 
    cat(paste(" bootstrap took",.timerraw(time1),"s.\n")) 
    return(list(bootreps=bootreps,RNGstate=RNGstate))
  }
})

