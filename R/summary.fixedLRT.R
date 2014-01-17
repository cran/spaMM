summary.fixedLRT <-
function(object,verbose=TRUE,...) {
  if (verbose) {
    cat(" ========      'full' model:     ========\n")    
    summary(object$fullfit,...) 
    cat(" ========      'null' model:     ========\n")    
    summary(object$nullfit,...) 
    cat(" ======== Likelihood ratio test: ========\n")    
  }
  outst <- paste(" LR statistic (",object$df," df): ",signif(object$LRTori,3),sep="")    
  cat(outst)
  if (!is.null(object$meanbootLRT)) {
    X2 <- object$LRTori*object$df/object$meanbootLRT
    outst <- paste("\n Bartlett-corrected LR statistic (",object$df," df): ",signif(X2,3),sep="")    
    cat(outst)
  } else X2 <- object$LRTori
  outst <- paste(" (p = ",signif(1-pchisq(X2,df=object$df),3),")\n",sep="")    
  cat(outst)    
}
