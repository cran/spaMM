getCallHL <-
function(object) {
  if (is.null(call <- attr(object,"corrHLfitcall"))) {
    if (is.null(call <- attr(object,"HLCorcall"))) call <- getCall(object)
  }
  call
}
