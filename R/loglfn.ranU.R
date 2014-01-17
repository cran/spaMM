loglfn.ranU <-
function(RandDist,y,nu) { ## functions with standardized mean and only a dispersion param
   switch(RandDist,
      gaussian = {- ((y^2)*nu+log(2*pi/nu))/2}, 
      gamma = {-nu*y+nu*(log(nu*y))-lgamma(nu)-log(y)}, ## p. 180 with psi=1 gives log pdf ranV assuming V=logU
      beta = {(nu/2-1)*log(y*(1-y))-lbeta(nu/2,nu/2)}, ## version explained p. 181 LeeNP
      ## Log[PDF[InverseGammaDistribution[1 + \[Nu], \[Nu] \[Mu]], uh]] with Mu=1 + |du/dv|
      "inverse.gamma" = {-nu/y - (2+nu)* log(y) + (1+nu)*log(nu) - lgamma(1+nu)} ## p. 181 with psi=1 gives log pdf ranV assuming V=-1/U, not log pdf ranU
   )
}
