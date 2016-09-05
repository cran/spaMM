cat("\ntest of DHGLM:")

data("crack") # crack data, LeeNP06 chapter 11 etc
hlfit <- HLfit(y~crack0+(1|specimen),family=Gamma(log),data=crack, HLmethod="REML", 
               rand.family=inverse.Gamma(log), resid.model=list(formula=~cycle+(1|specimen))   )
expect_equal(hlfit$APHLs$p_v,789.60762,tolerance=1e-4)
expect_equal(hlfit$APHLs$p_bv,785.36745,tolerance=1e-4)
# definitions of the APHLs seem variable among DHGLM publications (?) 
# The cAIC compute here further differs from that reported by dhglm version 1.5
# as the latter does not count the parameter for lambda, while spaMM counts it
# Cf Vaida & Blanchard, Biometrika 2005, theorem 2 vs theorem 1.