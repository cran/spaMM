.reset_processed_sparse_precision <- function(processed,sparseBool) {
  if (  is.list(processed) ) {
    for (lit in seq_along(processed)) {
      processed[[lit]]$sparsePrecisionBOOL <- sparseBool
      if (sparseBool) {
        ZAfix <- as(processed[[lit]]$AUGI0_ZX$ZAfix,"sparseMatrix")
        processed[[lit]]$AUGI0_ZX$ZAfix <- ZAfix
        rsZA <- rowSums(ZAfix) ## test that there a '1' per row and 'O's otherwise:  
        processed[[lit]]$AUGI0_ZX$is_unitary_ZAfix <- (unique(rsZA)==1 && all(rowSums(ZAfix^2)==rsZA)) ## $ rather than attribute to S4 ZAfix
      }
    } 
  } else {
    processed$sparsePrecisionBOOL <- sparseBool
    if (sparseBool) {
      ZAfix <- as(processed$AUGI0_ZX$ZAfix,"sparseMatrix")
      processed$AUGI0_ZX$ZAfix <- ZAfix
      rsZA <- rowSums(ZAfix) ## test that there a '1' per row and 'O's otherwise:  
      processed$AUGI0_ZX$is_unitary_ZAfix <- (unique(rsZA)==1 && all(rowSums(ZAfix^2)==rsZA)) ## $ rather than attribute to S4 ZAfix
    }
  }
  return(processed) ## practically a reference or list of references
}
