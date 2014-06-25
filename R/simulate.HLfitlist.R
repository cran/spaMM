simulate.HLfitlist <-
function(object,nsim=1,seed=NULL,newX=object[[1]]$data,sizes=object[[1]]$weights,...) {
  ## RNG stuff copied from simulate.lm
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else { ## this makes changes to RNG local where 'seed' is used:
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  replicate(nsim, {
    allrownames <- unique(unlist(lapply(object,function(hl){rownames(hl$data)})))
    resu <- matrix(0,nrow=length(allrownames),ncol=length(object)) ## two cols if 3 types
    cumul <- 0
    if (length(sizes) != nrow(newX)) {
      mess <- pastefrom("length(sizes) != nrow(newX).",prefix="(!) From ")
      stop(mess)
    }
    for (it in seq(ncol(resu))) {
      ## it = 1 se ramène à simulate(object[[1]])
      resu[,it] <- simulate(object[[it]],newX=newX,sizes=sizes - cumul)
      cumul <- rowSums(resu)  
    }
    resu <- cbind(resu,sizes - cumul) ## now 3 cols if 3 types
    rownames(resu) <- allrownames
    colnames(resu) <- attr(object,"sortedTypes")
    as.data.frame(resu)
  },simplify=FALSE)
}
