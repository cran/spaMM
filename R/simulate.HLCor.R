simulate.HLCor <-
function(object, ...) { 
  object <- object$hlfit
  simulate.HLfit(object, ...)
}
