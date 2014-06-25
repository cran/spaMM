getgroups <-
function(formlist,data) { ## reduction extreme de nlme:::getGroups.data.frame
  vlist <- lapply(formlist, function(x) { 
    val <- eval(x[[length(x)]], data)
    if (length(val) == 1) {
      return(as.factor(rep(val, nrow(data))))
    }
    else {
      return(as.factor(val)[drop = TRUE])
    }
  })
  if (length(vlist) == 1) 
    return(vlist[[1]])
  value <- do.call("data.frame", vlist)
  return(value) 
}
