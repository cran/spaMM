predictionCoeffs <-
function(object) {
  v_h_coeffs <- object$v_h ## v_h for non spatial and will later contain coeffs for spatial
  LMatrix <- attr(object$predictor,"LMatrix") 
  if ( ! is.null(LMatrix)) {
    vec_n_u_h <- attr(object$lambda,"n_u_h")
    cum_n_u_h <- cumsum(c(0,vec_n_u_h))  
    ZAlist <- object$ZAlist 
    if ( ! is.list(LMatrix)) LMatrix <- list(LMatrix)
    for (Lit in seq_len(length(LMatrix))) {
      lmatrix <- LMatrix[[Lit]]
      affecteds <- which(attr(ZAlist,"ranefs") %in% attr(lmatrix,"ranefs"))
      for (it in affecteds) {
        u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
        v_h_coeffs[u.range] <- solve(t(lmatrix),v_h_coeffs[u.range])   
      }  
    }
  }
  return(v_h_coeffs)
}
