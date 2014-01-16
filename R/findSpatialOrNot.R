findSpatialOrNot <-
function (term) { ## derived from findbars
  if (is.name(term) || !is.language(term)) 
    return(NULL)
  if (term[[1]] == as.name("(")) ## i.e. (blob) but not ( | ) 
    return(findSpatialOrNot(term[[2]])) 
  if (!is.call(term)) 
    stop("term must be of class call")
  if (term[[1]] == as.name("|"))  ## i.e. ( | ) expression
    return(term)
  if (length(term) == 2) {
    term1 <- as.character(term[[1]])
    if (term1 %in% c("adjacency","Matern","corMatern")) {
      return(term) 
    } else return(NULL) 
  }
  c(findSpatialOrNot(term[[2]]), findSpatialOrNot(term[[3]]))
}
