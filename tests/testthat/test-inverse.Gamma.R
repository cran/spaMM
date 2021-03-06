cat(crayon::yellow("test of inverse Gamma:\n"))

data("wafers")
HLig <- HLfit( y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),
       family=Gamma(log),rand.family=inverse.Gamma(log),
       resid.model= ~ X3+I(X3^2) ,data=wafers)
testthat::expect_equal(HLig$APHLs$p_v,-1157.523,tolerance=1e-3)
