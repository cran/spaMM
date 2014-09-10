checkRandLinkS <-
function(rand.families) {
  unlist(lapply(rand.families,checkRandLink)) ## a vector of lcrandfamfam := tolower(rand.family$family)
}
