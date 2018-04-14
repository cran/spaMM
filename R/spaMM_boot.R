spaMM_boot <- function(object, simuland, nsim, nb_cores=NULL,
                       resp_testfn=NULL,
                       ...) { ## ... are arguments used by functions called by the simuland function
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
  ii <- 0 ## 'global definition' (!)
  nb_cores <- .check_nb_cores(nb_cores=nb_cores) 
  if (nb_cores > 1L) {
    #cl <- parallel::makeCluster(nb_cores,outfile="essai.txt") 
    cl <- parallel::makeCluster(nb_cores) 
    parallel::clusterEvalQ(cl, library("spaMM"))
    dotenv <- list2env(list(...))
    parallel::clusterExport(cl=cl, as.list(ls(dotenv)),envir=dotenv) ## much faster...
  } else cl <- NULL
  ####
  RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  newy_s <- simulate(object,nsim = nsim,verbose=FALSE,resp_testfn=resp_testfn) 
  if (nsim==1L) dim(newy_s) <- c(length(newy_s),1)
  pbopt <- pboptions(nout=min(100,2*nsim),type="timer")
  bootreps <- pbapply(X=newy_s,MARGIN = 2L,FUN = simuland, cl=cl, ...) ## F I X M E allow doSNOW again
  pboptions(pbopt)
  ####
  if (nb_cores > 1L) { parallel::stopCluster(cl) } 
  cat(paste(" bootstrap took",.timerraw(time1),"s.\n")) 
  return(list(bootreps=bootreps,RNGstate=RNGstate))
}

if (FALSE) {
  
  myfun <- function(y, what=NULL, ...) { ## only a named first argument and ...
    data <- br$nullfit$data
    data$succes <- y
    data$echec <- 2-y
    aslistnull <- as.list(getCall(br$fullfit)) ## to estimate the slope
    aslistnull$data <- data
    res <- eval(as.call(aslistnull))
    if (!is.null(what)) res <- eval(what)
    return(res)
  }
  
  spaMM_boot(br$nullfit,simuland = myfun,1L,br=br)[["bootreps"]] ## no 'what'
  
  big <- spaMM_boot(br$nullfit,simuland = myfun,7L,nb_cores=7,
                    ## ... args of 
                    what=quote(fixef(res)[2L]),br=br)[["bootreps"]]
  
  zut <- unlist(spaMM_boot(br$nullfit,simuland = myfun,nsim=7,
                             what=quote(fixef(res)[2L]),br=br)[["bootreps"]])
  
}
