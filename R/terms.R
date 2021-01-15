.get_from_HLframes <- function(object=NULL, HLframes=object$HLframes, which="fixef_off_terms",mv_it=NULL) {
  if ( ! is.null(object) && object$spaMM.version < "2.5.9") {
    return(HLframes$fixef_terms)
  } 
  if (which=="respname" || which=="respnames") {
    resu <- HLframes[['Y']]
  } else resu <- HLframes[[which]]
  if (which=="respnames") {
    resu <- sapply(resu, attr, which="respname")
  } else {
    if ( ! is.null(mv_it)) resu <- resu[[mv_it]]
    if (which=="respname") resu <- attr(resu,"respname")
  }
  resu
}

# Initially [for for MSFDR -> stats::step(); not directly called in spaMM code]
terms.HLfit <- function(x, ...) { ## the full formula with the attributes for the fixed effects only (OK for MSFDR -> stats::step())
  # distinct attributes for ranefs would surely work.
  form <- formula.HLfit(x, which="") ## hyper does not seem necessary (nor offset, probably but the attribute will keep offset info bc it's the info available)
  if (inherits(form,"list")) { 
    for (mv_it in seq_along(form)) {
      attributes(form[[mv_it]]) <- attributes(.get_from_HLframes(object=x, mv_it=mv_it))
    }
  } else attributes(form) <- attributes(.get_from_HLframes(object=x))
  return(form)
}
