summary.HLfitlist <-
function(object, ...) {
  sapply(object,summary.HLfit) ## because summary(list object) displays nothing (contrary to print(list.object)) => rather call summary(each HLfit object)
  cat(" -------- Global likelihood values  --------\n")    
  zut <- attr(object,"APHLs")
  cat(paste(names(zut),": ",signif(unlist(zut),6), sep="",collapse="; "),"\n")
  invisible(object)
}
