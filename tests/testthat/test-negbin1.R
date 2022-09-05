cat(crayon::yellow("\ntest neg.bin.1 by negbin[2]:\n"))

# Here 'negbin1' denotes the case where the variance of the negbin is a linear function of the mean.
# This can be represented as a Poisson-Gamma mixture model with heteroscedastic independent Gamma random effects.
# Let mu_negbin|u := mu_pois * u for u ~Gamma(mean 1 and variance vg)
# Then the marginal mu_negbin = mu_pois, and the marginal variance of the negbin is mu_pois*(1+mu_pois*vg).
# We let vg=disp/mu_pois so that the variance of the negbin is mu_pois*(1+disp)
# We fit this through a random effect (wei-1|id) with 'prior weights' wei=1/sqrt(mu_pois). 
# The estimated lambda parameter of the Gamma random effect is then disp.

# The main limitation of this approach may be that the computation burden of the ad-hoc random effect increases marked with nobs.

if (spaMM.getOption("example_maxtime")>1) {
  set.seed(123)
  nobs <- 400
  X <- cbind(1,env=runif(nobs)-0.5)
  etafix <- X %*% c(2,2)
  mufix <- exp(etafix)
  disp <- 2
  vargamma <- disp/mufix
  mu <- mufix*rgamma(nobs,shape=1/vargamma,scale=vargamma)
  wei <- rep(1,nobs)
  dat <- data.frame(y=rpois(nobs,lambda=mu),env=X[,"env"],
                    id=seq(nobs), # defines the independent Gamma random effects for the Poisson-Gamma mixture
                    wei=wei)
  
  {
    # Initialization
    mfit <- fitme(y~env+(wei-1|id),family=poisson(), 
                  rand.family=list(Gamma(log)), data=dat)
    sqmu <- sqrt(predict(mfit,re.form=NA)[,1]) ## re.form should include all random effects except the ad-hoc Gamma one
    wei <- 1/sqmu # it might be useful to normalize weights to normalize variance, 
                  # but the interpretation of the lambda estimate would then bemodified
    dat$wei <- wei
    it <- 0L
    # iterations 
    while (it <100) {
      mfit <- fitme(y~env+(wei-1|id),family=poisson(), 
                    rand.family=list(Gamma(log)), data=dat)
      sqmu <- sqrt(predict(mfit,re.form=NA)[,1]) ## re.form as above
      wei <- 1/sqmu 
      convcrit <- max(abs(range(wei-dat$wei))) ## convergence criterion is convergence of predictions
      print(paste("iter",it,":",convcrit," ",logLik(mfit)))
      if (convcrit <1e-3) break()
      it <- it+1L
      dat$wei <- wei 
    }
    testthat::expect_equal(mfit$lambda[1],c(id.wei=1.4710451324)) 
    ## estimation of the Gamma variance is not precise, but approaches disp=2 in larger samples
  }
  simulate(mfit)
}

