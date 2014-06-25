HL_process.args <-
function (...) { ## from gsl package
  a <- list(...)
  attr <- attributes(a[[which.max(unlist(lapply(a, length)))]])
  a <- lapply(a, as.vector)
  out <- do.call("rbind", a)
  out <- list(`1`=out[1,],`2`=out[2,]) ## == out <- split(out, row(out)) but profiling shows the latter is slow
  names(out) <- paste("arg", 1:length(a), sep = "")
  return(c(out, attr = list(attr)))
}
