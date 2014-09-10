spMMexpandSlash <-
function (bb) {
  ## Create the interaction terms for nested effects
  makeInteraction <- function(x)
  {
    if (length(x) < 2) return(x)
    trm1 <- makeInteraction(x[[1]])
    trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
    list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
  }
  ##
  if (!is.list(bb)) 
    return(spMMexpandSlash(list(bb)))
  #### ELSE :
    unlist(lapply(bb, function(x) {
      if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
        ## lapply(unlist(...)) - unlist returns a flattened list
        lapply(unlist(makeInteraction(trms)),
               function(trm) substitute(foo|bar, list(foo = x[[2]], bar = trm)))
      else x
    }))
}
