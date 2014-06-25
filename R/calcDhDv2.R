calcDhDv2 <-
function(ZAL,w.resid,w.ranef) {
  d2hdv2 <- - ZtWZ(ZAL,w.resid) 
  diag(d2hdv2) <- diag(d2hdv2) - w.ranef ## small benefit that diag(w.ranef) not called on length-1 w.ranef which may occasionally be meaningful
  d2hdv2
}
