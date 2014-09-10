getProcessed <-
function(object,element,idx=1) {
  if ( ! is.null(attr(object,"multiple"))) object <- object[[idx]]  
  eval(parse(text=paste("object$",element,sep="")))
}
