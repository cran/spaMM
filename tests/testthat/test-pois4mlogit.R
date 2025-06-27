cat(cli::col_yellow("\ntest multinomial logit:\n"))

#### Fitting a binomial(logit) model by a bivariate poisson(log) surrogate:
## Toy data 
set.seed(123)
ssize <- 10L
shape <- 0.35
toydata <- data.frame(
  yellow=rbinom(ssize, 16, prob=rbeta(ssize,shape,shape)),
  purple=rbinom(ssize, 16, prob=rbeta(ssize,shape,shape)), # (purple ignored below)
  blue=rbinom(ssize, 16, prob=rbeta(ssize,shape,shape)),
  phenotype=rnorm(ssize)
)
## Standard binomial fit (purple flowers are ignored here)
(byB <- fitme(cbind(yellow,blue) ~ phenotype, family = binomial(), 
      data=toydata))


{
  ## Surrogate fit
  (byP2 <- pois4mlogit(submodels = list(
    list(yellow ~ offset(.dynoffset) + phenotype, family = poisson()),
    list(blue ~ offset(.dynoffset) + phenotype, family = poisson())),
    data = toydata, types=c("yellow","blue")))
  
  #### Recover summaries of the binomial fit from the surrogate fit:
  
  ## Coefficients of binomial model recovered as:
  fixef(byP2)[1:2]-fixef(byP2)[3:4]
  
  ## SEs of coefficients of binomial model recovered as:  
  {
    P2B <- rbind(c(1,0,-1,0),c(0,1,0,-1))
    P2B %*% vcov(byP2) %*% t(P2B)
  } 
  # practically equivalent to 
  vcov(byB)
  
  ## logLiks
  # The logLiks differ between the binomial and poisson surrogate fits because the combinatorial coefficients differ.
  # The relationship between these logLiks is simple when the long form of the data is used:
  str(long2 <- reshape2long(toydata, c("yellow","blue")))
  
  # Fits on long data:
  (byBlong  <- fitme(cbind(yellow,blue) ~ phenotype, family = binomial(), 
                     data=long2))
  (byPlong <- pois4mlogit(submodels = list(
    list(yellow ~ offset(.dynoffset) + phenotype, family = poisson()),
    list(blue ~ offset(.dynoffset) + phenotype, family = poisson())),
    data = toydata, to.long=TRUE, types=c("yellow","blue")))
  
  # The logLik of the surrogate model is
  sum(byPlong$eta*byPlong$y) - sum(exp(byPlong$eta)) # = logLik(byPlong)
  # while le logLik of the binomial one can be obtained as 
  sum(byPlong$eta*byPlong$y) # = logLik(byBlong)
  # the difference is 156  __= total multinomial sample size__
  
  # This illustrates that the logLiks of different surrogate models differ 
  # only by a constant independent of the fitted model. 
  # Likelihood ratios between surrogate models are then 
  # identical to likelihood ratios between corresponding binomial models,
  # provided a single data format is used throughout the comparisons.
  
}   
     
{ ## NA handling
  toydataNA <- toydata
  toydataNA$pheno1 <- toydataNA$pheno2 <- toydataNA$pheno3 <- toydata$phenotype
  toydataNA$yellow[1] <- NA 
  toydataNA$pheno1[2] <- toydataNA$pheno2[2] <- NA 
  (byP3 <- pois4mlogit(submodels = list(
      list(yellow ~ offset(.dynoffset) + pheno1, family = poisson()),
      list(blue ~ offset(.dynoffset) + pheno2, family = poisson()),
      list(purple ~ offset(.dynoffset) + pheno3, family = poisson())),
      data = toydataNA, types=c("yellow","blue","purple")))
  length(predict(byP3)) # 27 bc only second draw cannot be predicted
  rowSums(matrix(predict(byP3), ncol=3))
}  

if (FALSE) { ## trying to test identifiability... but this converges immediately...
  X_4to3 <- 
    matrix(c(1,0,0,
             0,1,0,
             1,0,0,
             0,0,1), nrow=4, ncol=3, byrow=TRUE,
           dimnames=list(NULL, c("(Intercept)","phenotype_1","phenotype_2")))
  
  (unidentif <- pois4mlogit(submodels = list(
    list(yellow ~ offset(.dynoffset) + phenotype, family = poisson()),
    list(blue ~ offset(.dynoffset) + phenotype, family = poisson())),
    data = toydata, X2X=X_4to3, types=c("yellow","blue"), progress=2))
}
