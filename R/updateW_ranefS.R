updateW_ranefS <-
function(cum_n_u_h,rand.families,lambda,u_h,v_h) {
  nrand <- length(rand.families)
  blob <- lapply(seq(nrand), function(it) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    updateWranef(rand.family=rand.families[[it]],lambda[u.range],u_h[u.range],v_h[u.range])
  })
  w.ranef <- unlist(lapply(blob,function(b) {b$w.ranef}))
  w.ranef <- pmin(w.ranef,1e10) ## patch useful to avoid singular d2hdv2 in PLoG model
  dlogWran_dv_h <- unlist(lapply(blob,function(b) {b$dlogWran_dv_h}))
  dvdu <- unlist(lapply(blob,function(b) {b$dvdu}))
  return(list(w.ranef=w.ranef,dlogWran_dv_h=dlogWran_dv_h,dvdu=dvdu))
}
