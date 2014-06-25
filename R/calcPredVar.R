calcPredVar <-
function(phi,lambda,Corr,Coldnew,Cnewnew,oldX.pv,X.pv,ZA) {
  if (attr(ZA,"identityMatrix")) {
    V <- phi * diag(nrow(Corr))+lambda * Corr ## r*r
    invV.Corrnew <- solve(V,lambda * Coldnew) ## r*1
    predVar <- lambda * (Cnewnew - t(Coldnew) %*% invV.Corrnew)
    bigRHS <- solve(t(oldX.pv) %*% solve(V,oldX.pv),t(X.pv)- t(oldX.pv) %*% invV.Corrnew)
    predVar <- predVar + (X.pv %*% bigRHS- lambda* t(Coldnew) %*% solve(V,oldX.pv %*% bigRHS))
  } else {
    # Cnn -c'.inv(Coo).c + (x-X'.inv(V).c)'.inv(X'inv(V)X). (x-X'.inv(V).c)
    V <- phi * diag(nrow(ZA))+lambda * ZA %*% Corr %*% t(ZA)  ## obs*r.r*r.r*obs =obs*obs
    invV.ZA.Corrnew <- solve(V,lambda * ZA %*% Coldnew) ## obs*obs.obs*r.r*new = obs*new
    Corrnew.tZA <- t(Coldnew) %*% t(ZA) ## new*r.r*obs = new*obs
    predVar <- lambda * (Cnewnew - Corrnew.tZA %*% invV.ZA.Corrnew) ## new*new- new*obs.obs*new = new*new
    bigRHS <- solve(t(oldX.pv) %*% solve(V,oldX.pv),t(X.pv)-t(oldX.pv) %*% invV.ZA.Corrnew) ## (p*p).(p*new) = p*new
    predVar <- predVar + (X.pv %*% bigRHS - lambda * Corrnew.tZA %*% solve(V,oldX.pv %*% bigRHS)) ## new*obs . obs*obs . obs*p . p*new
  }
  return(predVar)
}
