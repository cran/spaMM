setProcessed <-
function(object,element,value=1) {
  if ( ! is.null(attr(object,"multiple"))) {
    for (nam in names(object)) {
      eval(parse(text=paste("object[[",nam,"]]$",element," <- ",value,sep="")))
    }
  } else eval(parse(text=paste("object$",element," <- ",value,sep="")))
  return(object)
}
