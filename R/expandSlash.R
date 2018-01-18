.slashTerms <- function (x) {
  if (!("/" %in% all.names(x))) 
    return(x)
  if (x[[1]] != as.name("/")) 
    stop("unparseable formula for grouping factor")
  list(.slashTerms(x[[2]]), .slashTerms(x[[3]]))
}

## current fns from lme4 for comparison are equivalent 
##  https://github.com/lme4/lme4/blob/master/R/utilities.R


.makeInteraction <- function(x) {
  if (length(x) < 2) return(x)
  trm1 <- .makeInteraction(x[[1]])
  trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
  list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
}

## Create the interaction terms for nested effects
.spMMexpandSlash <- function (bb) {
  if (!is.list(bb)) return(.spMMexpandSlash(list(bb)))
  #### ELSE :
  unlist(lapply(bb, function(x) {
    if (length(x) > 2 && is.list(trms <- .slashTerms(x[[3]])))
      ## lapply(unlist(...)) - unlist returns a flattened list
      lapply(unlist(.makeInteraction(trms)),
             function(trm) structure(substitute(foo|bar, list(foo = x[[2]], bar = trm)),
                                     type=attr(x,"type") ## replicates the type '(.|.)' of (1|./.) in each of the terms created 
             ))
    else x
  }))
}


