spaMM.options <-
function (...) {
    if (nargs() == 0) return(.SpaMM)
    current <- .SpaMM
    temp <- list(...)
    if (length(temp) == 1L && is.null(names(temp))) {
        arg <- temp[[1]]
        switch(mode(arg), list = temp <- arg, character = return(.SpaMM[arg]), 
            stop("invalid argument: ", sQuote(arg)))
    }
    if (length(temp) == 0) 
        return(current)
    n <- names(temp)
    if (is.null(n)) 
        stop("options must be given by name")
    current[n] <- temp
    if (sys.parent() == 0) {
      env <- asNamespace("spaMM")
    } else env <- parent.frame()
    assign(".SpaMM", current, envir = env) 
    invisible(current)
}
