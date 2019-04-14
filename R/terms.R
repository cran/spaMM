.get_fixef_off_terms <- function(object) {
  if (object$spaMM.version < "2.5.9") {
    return(object$HLframes$fixef_terms)
  } else return(object$HLframes$fixef_off_terms)
}

# Initially [for for MSFDR -> stats::step(); not directly called in spaMM code]
terms.HLfit <- function(x, ...) { ## the full formula with the attributes for the fixed effects only (OK for MSFDR -> stats::step())
  # distinct attributes for ranefs wold surely work.
  form <- formula.HLfit(x, which="") ## hper does not seem necessary (nor offset, probably but the attribute will keep offset info bc it's the info available)
  attributes(form) <- attributes(.get_fixef_off_terms(x))
  return(form)
}

# model.terms <- function(object, which="fixef_off", ...) { 
#   if (which=="fixef_off") {
#     .get_fixef_off_terms(object)
#   } else stop("'which' value not handled")
# }
# 
