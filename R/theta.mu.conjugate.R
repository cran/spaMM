theta.mu.conjugate <-
function(mu,family) { 
   ## theta(u) in LeeN01... this is more pedagogy than efficient code
   switch(tolower(family),
      gaussian = theta.mu.canonical(mu,"gaussian") , ## mu 
      gamma = theta.mu.canonical(mu,"poisson"), ## log(mu)
      beta = theta.mu.canonical(mu,"binomial"), ## improved logit(mu)      
      "inverse.gamma" = theta.mu.canonical(mu,"gamma") ## -1/mu
   )
}
