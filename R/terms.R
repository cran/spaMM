.get_fixef_off_terms <- function(object) {
  if (object$spaMM.version < "2.5.9") {
    return(object$HLframes$fixef_terms)
  } else return(object$HLframes$fixef_off_terms)
}

terms.HLfit <- function(x, ...) { ## the full formula with the attributes for the fixed effects only (OK for MSFDR -> stats::step())
  # distinct attributes for ranefs wold surely work.
  form <- x$predictor
  attributes(form) <- attributes(.get_fixef_off_terms(x))
  return(form)
}

# model.terms <- function(object, which="fixef_off", ...) { 
#   if (which=="fixef_off") {
#     .get_fixef_off_terms(object)
#   } else stop("'which' value not handled")
# }
# 
