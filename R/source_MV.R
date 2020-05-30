## The file is in .Rbuildignore

.fitmv <- function(...) {
  assign("spaMM_glm_conv_crit",list(max=-Inf) , envir=environment(spaMM_glm.fit))
  time1 <- Sys.time()
  oricall <- match.call(expand.dots=TRUE) ## mc including dotlist
  mc <- oricall
  for (nam in names(mc)[-1]) {
    call_it <- match.call(get("fitme", asNamespace("spaMM")), mc[[nam]]) # matches args, e.g. adding formula= if first arg was unnamed, necess for next line.
    oricall[[nam]]$formula <- .preprocess_formula(call_it$formula) # cf comments in fitme()
    call_it[[1L]] <- get(".preprocess_fitme", asNamespace("spaMM"))
    mc[[nam]] <- eval(call_it,parent.frame()) ## returns a call with changed arguments, typically with argument "processed"
    # could be evaluated by mc[[nam]][[1L]] <- get("fitmv_body", asNamespace("spaMM")) ; eval(mc[[nam]])
  }
  mc[[1L]] <- get("fitmv_body", asNamespace("spaMM")) 
  hlcor <- eval(mc,parent.frame()) 
  .check_conv_glm_reinit()
  if (inherits(hlcor,"HLfitlist")) {
    attr(hlcor,"call") <- oricall
  } else {
    for (nam in names(mc)[-1]) {
      oricall[[nam]]$control.dist <- mc[[nam]]$processed$control_dist ## but never in the fitme_body() call
      hlcor$call <- oricall ## this is a call to fitme()
    }
  }
  attr(hlcor,"HLCorcall") <- NULL # presumably no more needed
  lsv <- c("lsv",ls())
  if ( ! .is.multi(family) && ! is.call(hlcor) ) {
    hlcor$how$fit_time <- .timerraw(time1)
    hlcor$fit_time <- structure(hlcor$how$fit_time,
                                message="Please use how(<fit object>)[['fit_time']] to extract this information cleanly.")
  }
  rm(list=setdiff(lsv,"hlcor")) ## empties the whole local envir except the return value
  return(hlcor)
}

