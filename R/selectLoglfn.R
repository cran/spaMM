selectLoglfn <-
function(family) {
   family <- tolower(family)
   switch(family,
      gaussian = function(theta,y,nu) {nu*(theta*y-(theta^2)/2)- ((y^2)*nu+log(2*pi/nu))/2}, 
      poisson = function(theta,y,nu) {nu*(theta*y-exp(theta))   -  lfactorial(y)},
      binomial = function(theta,freqs,sizes,nu) {nu*sizes*(freqs*theta-log(1+exp(theta))) +lchoose(sizes,round(sizes*freqs))},
      # gamma = function(theta,y,nu) {nu*(y*theta+log(-theta))+nu*(log(nu*y))-lgamma(nu)-log(y)} ## mean mu=-1/th, **** var = mu^2 / vu ****
      # same bu using ad hoc C library...
      gamma = function(theta,y,nu) {
        disp <- 1/nu
        mu <- -1/theta
        dgamma(y, shape=nu , scale = -1/(nu*theta), log = TRUE) ## from Gamma(log)$aic
      }
    )
}
