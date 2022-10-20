cat(crayon::yellow("\ntest DHGLM:"))

data("crack") # crack data, LeeNP06 chapter 11 etc
hlfit <- HLfit(y~crack0+(1|specimen),family=Gamma(log),data=crack, HLmethod="REML", 
               rand.family=inverse.Gamma(log), 
               resid.model=list(formula=~cycle+(1|specimen),fixed=list(phi=NA))   ) 
# where phi=NA is a way to force estimatation of a scalar phi OF the residual model.
#  This replicates what Lee et al. did, but is not spaMM's default.

# testthat::expect_equal(hlfit$APHLs$p_v,789.60762,tolerance=1e-4)
# testthat::expect_equal(hlfit$APHLs$p_bv,785.36745,tolerance=1e-4)
if (how(hlfit,verbose=FALSE)$obsInfo) {
  testthat::expect_equal(AIC(hlfit,short.names = TRUE)[2],c(cAIC=-1606.466),tolerance=1e-4)
} else testthat::expect_equal(AIC(hlfit,short.names = TRUE)[2],c(cAIC=-1607.22365),tolerance=1e-4) 

set.seed(123)
simulate(hlfit)

if (FALSE) {  
  hlfit <- HLfit(y~crack0+(1|specimen),family=Gamma(log),data=crack, HLmethod="REML", 
                 rand.family=inverse.Gamma(log), resid.model=list(formula=~cycle+(1|specimen),fixed=list(lambda=0.666))   )
  # and with partially-fixed ranCoefs in the residual dispersion model:
  hlfit <- HLfit(y~crack0+(1|specimen),family=Gamma(log),data=crack, HLmethod="REML", 
                 rand.family=inverse.Gamma(log), resid.model=list(formula=~cycle+(cycle|specimen),fixed=list(ranCoefs=list("1"=c(NA,-0.5,NA))))   )
}

