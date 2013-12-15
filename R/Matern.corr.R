Matern.corr <-
function (d, rho=1, smoothness, nu=smoothness, Nugget=0L) { ## rho is alpha in fields
  ## ideally (but not necess) on a 'dist' so the diagonal is not  manipulated 
    if (any(d < 0)) 
        stop("distance argument must be nonnegative")
    dscal <- d * rho
    dscal[d == 0L] <- 1e-10 ## avoids errors on distance =0; but the resulting corrvals can be visibly < 1 for small nu 
    #con <- (2^(nu - 1)) * gamma(nu) ## overflow for nu > 151
    #con <- 1/con
    #corrvals <- con * (dscal^nu) * besselK(dscal, nu) ## besselK overflows for small x
    logcon <- (nu - 1)*log(2)+ lgamma(nu) 
    corrvals <- - logcon + nu*log(dscal)+ bessel_lnKnu(x=dscal, nu=nu) ## 
##    corrvals <- - logcon + nu*log(dscal)+ log(besselK(x=dscal, nu=nu)) ## function from package gsl
    corrvals <- exp(corrvals) 
    corrvals[d != 0L] <- (1-Nugget)* corrvals[d != 0L]
    corrvals[d == 0L] <- 1 ## 
    corrvals[corrvals < 1e-16] <- 0L ## an attempt to deal with problem in chol/ldl/svd which don't like 'nearly-identity' matrices
    return(corrvals)
}
