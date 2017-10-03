.modify_list <- function (x, val) { ## derived from utils::modifyList
  stopifnot(is.list(x), is.list(val))
  xnames <- names(x)
  vnames <- names(val)
  vnames <- vnames[nzchar(vnames)]
  for (v in vnames) {
    if (v %in% xnames) {
      if ( is.list(x[[v]]) && is.list(val[[v]])) {
        x[[v]] <- .modify_list(x[[v]], val[[v]])
      } else if ( ! is.null(nam <- names(val[[v]]))) {
        x[[v]][nam] <- val[[v]]
      } else x[[v]] <- val[[v]]
    } else x[[v]] <- val[[v]] 
  }
  x
}


.merge_parlist <- function(parlist=NULL, new, types) {
  if (is.null(parlist)) parlist <- structure(list(),types=list()) ## implicit or explicit NULL may occur
  if (is.null(new)) return(parlist)
  parlist <- .modify_list(parlist,new)
  if (is.null(types)) { ## then new ms already have a types attribute
    template <- attr(new,"types")
  } else {
    template <- unlist(new)
    template[] <- types
    template <- relist(template,new)
  }
  attr(parlist,"types") <- .modify_list(attr(parlist,"types"),template)
  return(parlist)
} 

if (FALSE) {
  blalist <- .merge_parlist(,new=list(phi=1,lambda=c(1,`2`=2)),types="bla")
  .merge_parlist(blalist,new=list(rho=2,lambda=c(`3`=3)),types="hop")
}