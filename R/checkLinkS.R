checkLinkS <-
function(rand.families) {
  unlist(lapply(rand.families,checkLink)) ## a vector of lcrandfamfam := tolower(rand.family$family)
}
