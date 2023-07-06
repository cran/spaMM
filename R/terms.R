.get_from_data_attrs <- function(object=NULL, which, mv_it=NULL) {
  if (which=="fixefpredvars") {
    # object must not be NULL
    if (object$spaMM.version > "3.6.20") {
      resu <- attr(object$data,"fixefpredvars") # notably different from the fixefvarnames for poly()
      if ( ! is.null(mv_it)) resu <- resu[[mv_it]]
    } else {
      resu <- object$HLframes$fixef_off_terms # < 3.6.27
      if ( ! is.null(mv_it)) resu <- resu[[mv_it]]
      resu <- attr(resu,"predvars") # old fitobject structure
    }
    return(resu)
  }
  if (which=="fixefvarnames") {
    # object must not be NULL
    if (object$spaMM.version > "3.6.20") {
      resu <- attr(object$data,"fixefvarnames")
      if ( ! is.null(mv_it)) resu <- resu[[mv_it]]
    } else {
      resu <- object$HLframes$fixef_off_terms # < 3.6.27
      if ( ! is.null(mv_it)) resu <- resu[[mv_it]]
      resu <- rownames(attr(resu,"factors")) # old fitobject structure
    }
    return(resu)
  }
}

.get_from_terms_info <- function(object=NULL, terms_info, which="fixef_off_terms",mv_it=NULL) {
  # Either we have input from a "processed" object: we can assume the most recent 'processed' structure and terms_info can (must) be directly specified
  # Or we have only a fitobject available, and the following should be back-compat within limits of features of older spaMM
  if ( ! is.null(object)) {
    if (object$spaMM.version < "2.5.9") {
      return(object$HLframes$fixef_terms) # irrespective of 'which': should be compatible which features of spaMM < 2.5.9
    } else if (object$spaMM.version < "3.6.27") {
      terms_info <- object$HLframes
    } else terms_info <- object$main_terms_info 
  } 
  if (which=="respname" || which=="respnames") { # not used through API
    resu <- terms_info[['Y']] # Y from .preprocess -> main_terms_info$Y <- .get_Y(full_frame=main_terms_info$mf, famfam=family$family)
  } else resu <- terms_info[[which]]
  if (which=="respnames") {
    resu <- sapply(resu, attr, which="respname")
  } else {
    if ( ! is.null(mv_it)) resu <- resu[[mv_it]]
    if (which=="respname") resu <- attr(resu,"respname")
  }
  resu
}

# ___F I X M E___ can only return the fixef term. Makes it the default of a more general extractor.
# Initially [for for MSFDR -> stats::step(); not directly called in spaMM code]
terms.HLfit <- function(x, ...) { ## the full formula with the attributes for the fixed effects only (OK for MSFDR -> stats::step())
  # distinct attributes for ranefs would surely work.
  form <- formula.HLfit(x, which="") ## hyper does not seem necessary (nor offset, probably but the attribute will keep offset info bc it's the info available)
  if (inherits(form,"list")) { 
    for (mv_it in seq_along(form)) {
      attributes(form[[mv_it]]) <- attributes(.get_from_terms_info(object=x, mv_it=mv_it))
    }
  } else attributes(form) <- attributes(.get_from_terms_info(object=x))
  return(form)
}
