d2mudeta2fn <-
function(link,mu=NULL,eta=NULL,BinomialDen=NULL) { ## d2 MuCOUNTS d etaFREQS^2
  switch(link,
      identity = 0,
      log = mu, 
      inverse = 2 * mu^3 , ## canonical for Gamma()
      ## next three make sense for Binomial response data
      logit = {muFREQS <- mu/BinomialDen;
               d2muFREQS <- muFREQS*(1-muFREQS)*(1-2*muFREQS);
               d2muFREQS * BinomialDen
               },
      probit = -eta*dnorm(eta) * BinomialDen,
      cloglog = exp(eta-exp(eta))*(1-exp(eta)) * BinomialDen ## D[1 - E^-E^\[Eta], {\[Eta], 2}]
  )
}
