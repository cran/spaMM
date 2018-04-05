.probitgemWrap <- function(chr_fnname,arglist, pack="probitgem") {
  if (length(grep(pack,packageDescription("spaMM")$Imports))) {
    ## then the necessary functions must be imported-from in the NAMESPACE  
    do.call(chr_fnname,arglist) 
  } else if (length(grep(pack,packageDescription("spaMM")$Suggests))) {
    ## then the necessary functions cannot be imported-from in the NAMESPACE  (and the package must be written in an appropriate way)
    if ( requireNamespace(pack, quietly = TRUE)) {
      myfun <- get(chr_fnname, asNamespace(pack)) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
      do.call(myfun,arglist) 
    } else {stop(paste("'",pack,"' required but not available.",sep=""))}
  } else { ## package not declared in DESCRIPTION; private for example
    if (do.call("require",list(package=pack, quietly = TRUE))) {
      do.call(chr_fnname,arglist) 
    } else {stop(paste("'",pack,"' required but not available.",sep=""))}
  }
}

