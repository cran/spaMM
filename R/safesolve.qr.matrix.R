safesolve.qr.matrix <-
function(qr.a,B,silent=T,stop.on.error=T) { ## solve.qr with fall-back; qr.a should be a qr object, B a matrix
  ## there was a 'Matrix' subcode prior to 10/03/2013
    res <- try(solve.qr(qr.a,B),silent=silent)
    if (class(res)=="try-error") { ##FR->FR omprendre pourquoi c'est utile (systematique pour blackcap... Re: phi -> 0 ? )
      pivI <- sort.list(qr.a$pivot)
      solveA <- try(solve(qr.R(qr.a)[,pivI])%*%t(qr.Q(qr.a)),silent=silent)  
      ## solveA is solve(<original 'a' matrix>) using the QR decomp... but this may work when solve.qr fails !
      if (class(solveA)=="try-error") {
        if (stop.on.error) {
          mess <- pastefrom("class(solveA)='try-error'.",prefix="(!) From ")
          message(mess)
          stop("More code is needed to handle this... I exit.") ## perhaps recover A by qr.X and solve(A) ?
        } else return(solveA) ## passes control to calling function
      } else res <- solveA %*% B  
    }
  return(res)
}
