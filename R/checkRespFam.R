checkRespFam <-
function(family) {
  ## four lines from glm()
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  family
}
