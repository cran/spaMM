summary.fixedLRT <-
function(object,verbose=TRUE,...) {
  print(object$basicLRT)
  #   if (verbose) {
  #     cat(" ========      'full' model:     ========\n")    
  #     summary(object$fullfit,...) 
  #     cat(" ========      'null' model:     ========\n")    
  #     summary(object$nullfit,...) 
  #     cat(" ======== Likelihood ratio test: ========\n")    
  #   }
  #   outst <- paste(" LR statistic (",object$df," df): ",signif(object$LRTori,3),sep="")    
  #   cat(outst)
  bootInfo <- object$bootInfo
  if (!is.null(bootInfo)) {
    cat(" ======== Bootstrap: ========\n")    
    outst <- paste("Raw simulated p-value: ",signif(bootInfo$rawPvalue,3),sep="")    
    cat(outst)
    X2 <- bootInfo$LRTcorr
    outst <- paste("\nBartlett-corrected LR statistic (",bootInfo$df," df): ",signif(X2,3),sep="")    
    cat(outst)
  } 
}
