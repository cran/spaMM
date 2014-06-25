rntpos <-
function(n,mu,sigma2)
{
  # produce n samples from the
  # specified left-truncated to 0 gaussian
  # distribution
  -rntneg(n,-mu,sigma2)
}
