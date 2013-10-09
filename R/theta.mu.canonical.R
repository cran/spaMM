theta.mu.canonical <-
function(mu,family) { 
   ## the (fixed) canonical link between theta and mu, not the family link between eta and mu 
   switch(tolower(family),
      gaussian = mu ,
      poisson = log(mu) ,
      binomial = make.link("logit")$linkfun(mu),  # correct syntax, does not use 'non-public API' such as .Call to access code from dlls from the base packages...
## if this does no work, use 
#                 { 
#                    theta <- logit(mu)
#                    theta[theta>27.6310] <- 27.6310 ## mu>1-1e-12
#                    theta[theta < -27.6310] <- -27.6310 
#                    theta
#                 },
      gamma = -1/mu, ## "-": McC&N p. 290
   )
}
