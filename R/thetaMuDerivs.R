thetaMuDerivs <-
function(mu,BinomialDen,familyfam) {
  if (familyfam=="binomial") muFREQS <- mu/BinomialDen
  ## these definitions depend only on the canonical link
  Dtheta.Dmu <- switch(tolower(familyfam),
    gaussian = rep(1,length(mu)) ,
    poisson = 1/mu ,
    binomial = 1/(muFREQS*(1-muFREQS)),
    gamma = 1/mu^2
  ) ## values for given mu
  if (familyfam=="binomial") Dtheta.Dmu <- Dtheta.Dmu/BinomialDen
  D2theta.Dmu2 <- switch(tolower(familyfam),
    gaussian = rep(0,length(mu)) ,
    poisson = -1/mu^2 ,
    binomial = -(1-2*muFREQS)/(muFREQS*(1-muFREQS))^2,
    gamma = -2/mu^3
  ) ## values for given mu
  if (familyfam=="binomial") D2theta.Dmu2 <- D2theta.Dmu2/(BinomialDen^2)
  return(list(Dtheta.Dmu=Dtheta.Dmu,D2theta.Dmu2=D2theta.Dmu2))
}
