findSpatial <-
function (term) { ## derived from findbars
    if (is.name(term) || !is.language(term)) 
        return(NULL)
    if (term[[1]] == as.name("(")) ## i.e. ( | ) expression
        return(NULL) 
    if (!is.call(term)) 
        stop("term must be of class call")
    if (term[[1]] == as.name("|")) 
        return(term)
    if (length(term) == 2) {
      term1 <- as.character(term[[1]])
      if (term1 %in% c("adjacency","Matern","corMatern")) {
        return(term) 
      } else return(NULL) ## or findSpatial[[2]] ?
     }
    c(findSpatial(term[[2]]), findSpatial(term[[3]]))
}
