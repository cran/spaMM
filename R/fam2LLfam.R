
# The <stats:: family>2LLF functions set LLgeneric to TRUE, so they must add the member functions required by .add_Md_logcLdeta_terms() 

.Gamma2LLF <- function(family) {
  # clik is dgamma(y, shape=nu , scale = attr(theta,"mu")/nu, log = TRUE)
  DlogLDmu <- function(mu, y, wt, n, phi, nu=1/phi) { drop(wt*nu*(y/mu-1)/mu)} # or (y-mu)/V(mu), general result for GLMs
  D2logLDmu2 <- function(mu, y, wt, n, phi, nu=1/phi) { drop(wt*nu*(-2*y/mu+1)/(mu^2)) }
  D3logLDmu3 <- function(mu, y, wt, n, phi, nu=1/phi) { drop(wt*nu*(6*y/mu-2)/(mu^3))}
  D2muDeta2 <- .D2muDeta2(family$link)
  D3muDeta3 <- .D3muDeta3(family$link)
  #
  environment(DlogLDmu) <- environment(D2logLDmu2) <- environment(D3logLDmu3) <- environment(D2muDeta2) <- environment(D3muDeta3) <- environment(family$aic)
  family$DlogLDmu <- DlogLDmu
  family$D2logLDmu2 <- D2logLDmu2
  family$D3logLDmu3 <- D3logLDmu3
  family$D2muDeta2 <- D2muDeta2
  family$D3muDeta3 <- D3muDeta3
  family$flags <- list(obs=TRUE, # has info to fit by obs Info matrix
                       exp=TRUE, # has info to fit by exp Info matrix
                       LLgeneric=TRUE)
  class(family) <- c("LLF","family")
  return(family)
}

.gaussian2LLF <- function(family) {
  DlogLDmu <- function(mu, y, wt, n, phi, nu=1/phi) { drop(wt*nu*(y-mu))}
  D2logLDmu2 <- function(mu, y, wt, n, phi, nu=1/phi) { res <- drop(-wt*nu); if(length(res)==1L) {res <- rep(res,length(y))}; res}
  D3logLDmu3 <- function(mu, y, wt, n, phi, nu=1/phi) { rep(0,length(y))}
  D2muDeta2 <- .D2muDeta2(family$link)
  D3muDeta3 <- .D3muDeta3(family$link)
  #
  environment(DlogLDmu) <- environment(D2logLDmu2) <- environment(D3logLDmu3) <- environment(D2muDeta2) <- environment(D3muDeta3) <- environment(family$aic)
  family$DlogLDmu <- DlogLDmu
  family$D2logLDmu2 <- D2logLDmu2
  family$D3logLDmu3 <- D3logLDmu3
  family$D2muDeta2 <- D2muDeta2
  family$D3muDeta3 <- D3muDeta3
  family$flags <- list(obs=TRUE, exp=TRUE, LLgeneric=TRUE)
  class(family) <- c("LLF","family")
  return(family)
}



.binomial2LLF <- function(family) {
  DlogLDmu <- function(muCOUNT, muFREQS, y, BinomialDen) {(y-muCOUNT)/(muFREQS*(1-muFREQS)) }#  DlogLDmuFREQS ... 
  D2logLDmu2 <- function(muFREQS, y, BinomialDen) { (-BinomialDen + y)/(1 - muFREQS)^2  -  y/muFREQS^2 }
  D3logLDmu3 <- function(muFREQS, y, BinomialDen) {2 *( -(BinomialDen - y)/(1 - muFREQS)^3  +  y/muFREQS^3 ) }
  D2muDeta2 <- .D2muDeta2(family$link)
  D3muDeta3 <- .D3muDeta3(family$link)
  #
  environment(DlogLDmu) <- environment(D2logLDmu2) <- environment(D3logLDmu3) <- environment(D2muDeta2) <- environment(D3muDeta3) <- environment(family$aic)
  family$DlogLDmu <- DlogLDmu
  family$D2logLDmu2 <- D2logLDmu2
  family$D3logLDmu3 <- D3logLDmu3
  family$D2muDeta2 <- D2muDeta2
  family$D3muDeta3 <- D3muDeta3
  family$flags <- list(obs=TRUE, exp=TRUE, LLgeneric=TRUE)
  class(family) <- c("LLF","family")
  return(family)
}

.statsfam2LLF <- function(family) {
  famfam <- family$family
  flags <- family$flags
  # Conversions from from stats:: to LLF
  if (famfam=="gaussian") {
    family <- .gaussian2LLF(family) 
  } else if (famfam=="Gamma") {
    family <- .Gamma2LLF(family)  # from stats:: to LLF
  } else if (famfam=="binomial") { 
    family <- .binomial2LLF(family)
  } else if (famfam=="poisson") { # stats::poisson here ([T]Poisson would already be spaMM's *P*oisson) =>
    family <- Poisson(link=family$link) # untruncated here since was stats::poisson. And for LLgeneric=FALSE, the user should have called *P*oisson
  }
  # for .merge_processed() -> .statsfam2LLF() in particular, keep these flags from submodel processing:
  family$flags$LMbool <- flags$LMbool
  family$flags$canonicalLink <- flags$canonicalLink
  family
}

