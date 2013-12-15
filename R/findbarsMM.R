findbarsMM <-
function (term) 
{
    if (is.name(term) || !is.language(term)) 
        return(NULL)
    if (term[[1]] == as.name("(")) 
        return(findbarsMM(term[[2]]))
    if (!is.call(term)) 
        stop("term must be of class call")
    if (term[[1]] == as.name("|")) 
        return(term)
    if (length(term) == 2) 
        return(findbarsMM(term[[2]]))
    c(findbarsMM(term[[2]]), findbarsMM(term[[3]]))
}
