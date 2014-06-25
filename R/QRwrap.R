QRwrap <-
function(mat) {
  if (.SpaMM$USEEIGEN) {
    QR <- Rcpp_QR(mat) ## the main benefit is that the solve method for Rcpp-QR objects is numerically more stable 
#    if(any(is.nan(QR$R))) browser() 
  } else QR <-qr(mat)
  return(QR)
}
