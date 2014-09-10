eval.update.call <-
function(mc,...) {
  mc <- as.list(mc)
  dotlist <- list(...)
  # mc[names(dotlist)] <- dotlist
  eval(as.call(mc))  
}
