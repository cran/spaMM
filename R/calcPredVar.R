calcPredVar <-
function(phi,lambda,Coldnew,Cnewnew,oldX.pv,X.pv,ZA,newZA,beta_cov,Sig) { 
## the diag of predVar requires only the diag of Cnewnew but this leads nowhere
  if (attr(ZA,"identityMatrix")) {tZA <- ZA} else tZA <- t(ZA) ## to keep the attribute...
  if (attr(newZA,"identityMatrix")) {tnewZA <- newZA} else tnewZA <- t(newZA)
  #   if (FALSE && attr(ZA,"identityMatrix")) { ## FR->FR if FALSE
  #     V <- phi * diag(nrow(oldCorr))+lambda * oldCorr ## r*r ##FR->FR code pas general si  phi est modélisé !!! 
  #     invV.Corrnew <- solve(V,lambda * Coldnew %*% t(newZA)) ## r* (# of points to predict) 
  #     predVar <- lambda * newZA %*% (Cnewnew %*% t(newZA)- t(Coldnew) %*% invV.Corrnew)
  #     bigRHS <- solve(t(oldX.pv) %*% solve(V,oldX.pv),t(X.pv)- t(oldX.pv) %*% invV.Corrnew)
  #     predVar <- predVar + (X.pv %*% bigRHS- lambda* newZA %*% t(Coldnew) %*% solve(V,oldX.pv %*% bigRHS))
  #   } else if(FALSE) { ## old version pas lisible
  #     # Cnn -c'.inv(Coo).c + (x-X'.inv(V).c)'.inv(X'inv(V)X). (x-X'.inv(V).c)
  #     V <- phi * diag(nrow(ZA))+lambda * ZA %id*id% oldCorr %id*id% tZA  ## obs*r.r*r.r*obs =obs*obs
  #     invV.ZA.Corrnew <- solve(V,lambda * ZA %id*id% Coldnew %id*id% tnewZA) ## obs*obs.obs*r.r*new = obs*new
  #     Corrnew.tZA <-  t(Coldnew) %id*id% tZA ## new*r.r*obs = new*obs
  #     predVar <- lambda * newZA %id*id% ( Cnewnew %id*id% tnewZA- Corrnew.tZA %*% invV.ZA.Corrnew) ## new*new- new*obs.obs*new = new*new
  #     bigRHS <- solve(t(oldX.pv) %*% solve(V,oldX.pv),t(X.pv)-t(oldX.pv) %*% invV.ZA.Corrnew) ## (p*p).(p*new) = p*new
  #     predVar <- predVar + (X.pv %*% bigRHS - lambda *  newZA %id*id% Corrnew.tZA %*% solve(V,oldX.pv %*% bigRHS)) ## new*obs . obs*obs . obs*p . p*new
  #   } else {
    # Cnn -c'.inv(Coo).c + (x-X'.inv(V).c)'.inv(X'inv(V)X). (x-X'.inv(V).c)
  ZCZt_on <- ZA %id*id% Coldnew %id*id% tnewZA ## obs*new
  invSig.ZCZt_on <- solve(Sig,lambda * ZCZt_on) ## obs*obs.obs*new = obs*new
  if (attr(newZA,"identityMatrix")) {ZCZt_nn <- Cnewnew} else ZCZt_nn <- newZA %*% Cnewnew %*% tnewZA
  ## for known beta:
  predVar <- lambda * (ZCZt_nn- t(ZCZt_on) %*% invSig.ZCZt_on) ## new*new- new*obs.obs*new = new*new
  ## correction for estimated beta:
  K <- t(X.pv)-t(oldX.pv) %*% invSig.ZCZt_on ## p*new -(p*obs).(obs*new) = p*new
  predVar <- predVar + t(K) %*% beta_cov %*% K 
  # }
  return(predVar)
}
