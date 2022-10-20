.future_apply <- function(newresp, cluster_args, nb_cores, fn, iseed, steps, ...) {
  has_progressr <- ("package:progressr" %in% search())
  if (has_progressr) {
    # progressor is the only progress function that 'works' with mclapply
    # although not with load-balancing (mc.preschedule=FALSE)
    # Here we use the default (no balancing), and it is the steps with max value shown below that are reported.  
    progressor_ <- get("progressor", asNamespace("progressr"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
    with_progress_ <- get("with_progress", asNamespace("progressr"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
    if (nb_cores>1L) {
      if (cluster_args$type=="FORK") {
        pb_char <- "F"
      } else pb_char <- "P"
    } else pb_char <- "S"
    handlers_ <- get("handlers", asNamespace("progressr"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
    handler_txtprogressbar_ <- get("handler_txtprogressbar", asNamespace("progressr"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
    handlers_(handler_txtprogressbar_(char = pb_char))
    with_progress_({
      p <- progressor_(steps=steps)
      p_fn <- function(newy, ...) {
        res <- fn(newy, ...)
        p() # p() call necessary for actual progress report 
        res
      }
      bootreps <- try(future.apply::future_apply(X=newresp, MARGIN = 2L, FUN=p_fn, ..., future.seed = iseed))
    })
  } else {
    .warn_once_progressr()
    bootreps <- try(future.apply::future_apply(X=newresp, MARGIN = 2L, FUN=fn, ..., future.seed = iseed))
  }
  bootreps
}


# fn more generic than spaMM_boot: there is no call to other spaMM fns such as simulate(object, .) so this acts as a general wrapper for 
# foreach or pbapply, and not specifically for bootstrap computations.
dofuture <- function(newresp, fn, nb_cores=NULL, 
           fit_env, control=list(), cluster_args=NULL,
           debug.=FALSE, iseed=NULL, showpbar="ignored",
           pretest_cores="ignored",
           ... # passed to future.apply::future_apply then possibly to fn
) {
  if ( requireNamespace("future", quietly = TRUE) &&  requireNamespace("future.apply", quietly = TRUE)) { # both in Suggests
    if (is.list(fit_env)) fit_env <- list2env(fit_env)
    cluster_args <- .set_cluster_type(cluster_args, nb_cores=nb_cores)
    nb_cores <- cluster_args$spec
    if (debug. && nb_cores>1L ) debug. <- 1L 
    assign("debug.", debug., environment(fn))
    if (is.null(dim(newresp))) newresp <- matrix(seq(newresp),ncol=newresp,nrow=1) # assuming newresp is an integer
    nsim <- ncol(newresp)
    time1 <- Sys.time() 
    if (nb_cores>1L) {
      if (cluster_args$type=="FORK") {
        cl <- parallel::makeForkCluster(nnodes = nb_cores) 
        future::plan(future::cluster, workers=cl)
        bootreps <- .future_apply(newresp, cluster_args, nb_cores, fn, iseed, steps=ceiling(ncol(newresp)/nb_cores), ...)
      } else { # PSOCK
        cl <- do.call(parallel::makeCluster, cluster_args) # create a *socket* cluster
        future::plan(future::cluster, workers=cl)
        packages2export <- control$.packages
        if (is.null(packages2export)) packages2export <- "spaMM"
        parallel::clusterCall(cl,
                              function(packages) {for (p in packages) library(p, character.only = TRUE)}, 
                              packages2export)
        if (is.environment(fit_env)) try(parallel::clusterExport(cl=cl, varlist=ls(fit_env), envir=fit_env)) 
        bootreps <- .future_apply(newresp, cluster_args, nb_cores, fn, iseed, steps=ceiling(ncol(newresp)/nb_cores), ...)
      } # FORK ... else
      parallel::stopCluster(cl)
    } else { ## nb_cores=1L
      future::plan("sequential")
      bootreps <- .future_apply(newresp, cluster_args, nb_cores, fn, iseed, steps=ncol(newresp), ...)
    }
    if (identical(control$.combine,"rbind")) bootreps <- t(bootreps)
    cat(paste(" bootstrap took",.timerraw(time1),"s.\n")) 
    return(bootreps)
  } else {stop(paste("'future.apply' required but not both available.",sep=""))}
}
