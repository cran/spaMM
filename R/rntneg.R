rntneg <-
function(n,mu,sigma2)
{
  # produce n samples from the
  # specified rigth-truncated to 0 gaussian
  # distribution
  pn <- runif(n)*pnorm(0,mu,sqrt(sigma2))
  pn[pn==0] <-  .Machine$double.eps ## because if mu is large -> qnorm(0) is -Inf which later cause NaN's
  pn[pn==1] <-  1-.Machine$double.eps 
  qnorm(pn,mu,sqrt(sigma2))
  ## alternatively use ... qn[pn==0] <- ... 
}
