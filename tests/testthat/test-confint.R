cat(crayon::yellow("\ntest confint() for HLfit, corrHLfit, and fitme (with sparse_precision=",
                   capture.output(spaMM.getOption("sparse_precision")),"):\n"))
data("wafers")
wfit <- HLfit(y ~X1+(1|batch),family=Gamma(log),data=wafers,HLmethod="ML")
ci <- confint(wfit,"X1",verbose=FALSE)
testthat::expect_equal(ci$interval[[1]],0.01828659,tolerance=1e-4)
testthat::expect_equal(ci$interval[[2]],0.17271333,tolerance=1e-4)

# HGLM... DHGLM!
wfit <- HLfit(y ~ X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch), family=Gamma(log),HLmethod="ML", 
              rand.family=inverse.Gamma(log),
              resid.model = ~ X3+I(X3^2) , data=wafers)
ci <- confint(wfit,"X1",verbose=FALSE)
testthat::expect_equal(ci$interval[[1]],0.0361157,tolerance=1e-4)
testthat::expect_equal(ci$interval[[2]],0.1313484,tolerance=1e-4)


#### Checks of consistency of procedures for profiling out one or two parameters (with fixed phi and lambda only for a faster test).
## The CI's of the more constrained model should be within the other (even if both fits coincide at the ML, the additional constraint may matter at the bounds)

## with corrHLfit
data("blackcap")
fitobject <- corrHLfit(migStatus ~ means + Matern(1|longitude+latitude),data=blackcap,HLmethod="ML",
                       init.corrHLfit=list(nu=0.535929513,rho=0.007485725),ranFix=list(phi=0.05,lambda=2))
ci <- confint(fitobject,"means",verbose=FALSE)
testthat::expect_equal(ci$interval[[1]],0.1606797,tolerance=1e-4) 
testthat::expect_equal(ci$interval[[2]],0.9558826,tolerance=1e-4)
refitobject <- corrHLfit(migStatus ~ means + Matern(1|longitude+latitude),data=blackcap,HLmethod="ML",
                         init.corrHLfit=list(rho=0.007485725),ranFix=list(nu=0.535929513,phi=0.05,lambda=2))
ci <- confint(refitobject,"means",verbose=FALSE)
testthat::expect_equal(ci$interval[[1]],0.1586536,tolerance=1e-4)
testthat::expect_equal(ci$interval[[2]],0.9570798,tolerance=1e-4)

## with fitme (derived from gentle intro)
rSample <- function(nb,rho,sigma2_u,resid,intercept,slope,pairs=TRUE) {
  ## sample pairs of adjacent locations
  if (pairs) {
    x <- rnorm(nb/2); x <- c(x,x+0.001)
    y <- rnorm(nb/2); y <- c(y,y+0.001)
  } else {x <- rnorm(nb);y <- rnorm(nb)}
  dist <- dist(cbind(x,y)) ## distance matrix between locations
  m <- exp(-rho*as.matrix(dist)) ## correlation matrix
  b <- MASS::mvrnorm(1,rep(0,nb),m*sigma2_u) ## correlated random effects
  pred <- sample(nb) ##  some predictor variable
  obs <- intercept+slope*pred + b +rnorm(nb,0,sqrt(resid)) ## response
  data.frame(obs=obs,x,y,pred=pred)
}

rngcheck <- ("sample.kind" %in% names (formals(RNGkind)))
if (rngcheck) suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding"  )) 
set.seed(123)
d1 <- rSample(nb=40,rho=3,sigma2_u=0.5,resid=0.5,intercept=-1,slope=0.1)
if (rngcheck) RNGkind("Mersenne-Twister", "Inversion", "Rejection"  )
# Rnews: The included LAPACK sources have been updated to 3.10.1. ... 
# ... Using 3.10.x may give some different signs from earlier versions in SVDs or eigendecompositions... 
# => MASS::mvrnorm affected.
d1$obs <- c(1.56404438213562,0.838761745570047,1.18539152814119,2.77609094294087,-0.362117933355623,1.05566922884628,1.04891131933423,-1.36129005102859,
            3.14126410723178,-0.267510545755898,1.68242492955138,-0.294886401760225,1.70572380683716,2.25124059611453,-0.850815921177352,0.439615106586318,
            2.0392252364533,2.42453461515297,-0.53309120933375,-0.554803498441966,1.71344222448971,1.43407925526473,1.27070104772916,-0.557869677646131,
            1.29371981088591,1.98292535259867,-0.440298987269969,0.493786122749392,0.0552051443241781,0.943700535986965,-0.0148351370346089,0.260693893323206,
            2.45333962421848,0.276596082599056,-1.53290736864249,1.88815293467526,-1.90295848044005,2.84485580265504,-0.514747120202024,-2.44239024049091)

HLMf <- fitme(obs~pred+Matern(1|x+y),init=list(rho=59.11287),fixed=list(nu=48.96201,phi=0.447761,lambda=0.3697),data=d1,method="ML")
ci <- confint(HLMf,"pred",verbose=FALSE) 
testthat::expect_equal(ci$interval[[1]],0.06483438,tolerance=1e-4)  
testthat::expect_equal(ci$interval[[2]],0.10567544,tolerance=1e-4)  
HLM <- fitme(obs~pred+Matern(1|x+y),init=list(nu=48.96201,rho=59.11287),fixed=list(phi=0.447761,lambda=0.3697),data=d1,method="ML")
ci <- confint(HLM,"pred",verbose=FALSE) ## practically identical in the two fits (+ nu drifts to higher values whatever the initial one)
testthat::expect_equal(ci$interval[[1]],0.06483437,tolerance=1e-4)  
testthat::expect_equal(ci$interval[[2]],0.10567543,tolerance=1e-4)

# compar to lme4
if(requireNamespace("lme4", quietly = TRUE)) {
  data("sleepstudy",package = "lme4")
  mlfit <- fitme(Reaction ~ Days + (1|Subject), data = sleepstudy)
  fitci <- confint(mlfit, parm = "Days",verbose=FALSE)$interval
  rl_mer <- lme4::lmer(Reaction ~ Days + (1|Subject), data = sleepstudy)
  ci_mer <- suppressMessages(suppressWarnings(confint(rl_mer)))[4, ] ## suppress bobyqa conv warning in lme4:::optwrap, and message "Computing profile confidence intervals ..."
  testthat::expect_true(diff(range(ci_mer-fitci))<1e-5) 
}

# GLM with a 'method'...
library(spaMM)
set.seed(123)
iris$bidon <- rbinom(nrow(iris), 1, 0.5)
fML <- fitme(bidon ~ Species, data = iris, family = binomial(link = "logit"))
fPQL <- fitme(bidon ~ Species, data = iris, family = binomial(link = "logit"), method = "PQL/L")
testthat::expect_true(identical(confint(fML, "Speciesversicolor", verbose = FALSE)$interval, 
                                confint(fPQL, "Speciesversicolor", verbose = FALSE)$interval))
  
