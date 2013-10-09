summary.HLCor <-
function(object, ...) {
  summary(object$hlfit) ## currently not more than the corresponding HLfit summary
  invisible(object)
}
