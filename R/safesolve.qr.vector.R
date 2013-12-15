safesolve.qr.vector <-
function(qr.a,b,silent=T,stop.on.error=T) { ## solve.qr with fall-back; qr.a should be a qr object, b must be a vector
  ## there was a 'Matrix' subcode prior to 10/03/2013
    res <- try(solve.qr(qr.a,b),silent=silent)
    if (class(res)=="try-error") {   ## then some weird code, but...
      ## we try to solve(<original 'a' matrix>) using the QR decomp... this may work when solve.qr fails !
      ## The following is equivalent to solveA <- try(solve(qr.R(qr.a)[,pivI])%*%t(qr.Q(qr.a)),silent=silent) then return solveA %*% b 
      ## but uses only backward mutliplication of vectors and transpose of vectors  
      res <- t(t(b) %*% (qr.Q(qr.a))) ## not yet res
      pivI <- sort.list(qr.a$pivot)
      res <- try(solve(qr.R(qr.a)[, pivI]) %*% res, silent=silent)
      if (class(res)=="try-error") {
        if (stop.on.error) {
          mess <- pastefrom("class(solveA)='try-error'.",prefix="(!) From ")
          message(mess)
          stop("More code is needed to handle this... I exit.") ## perhaps recover A by qr.X and solve(A) ?
        } else return(res) ## passes control to calling function
      } 
    }  
  return(res)
}
