summary.corrHLfit <-
function(object, ...) {
  summobject <- object$hlfit
  summobject$ranFixNames <- object$ranFixNames
  summobject$objective <- object$objective
  summary(summobject) ## 
  invisible(object)
}
