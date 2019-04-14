.probitgemWrap <- function(chr_fnname,arglist, pack="probitgem") {
  if (length(grep(pack,packageDescription("spaMM")$Imports))) {
    ## then the necessary functions must be imported-from in the NAMESPACE  
    do.call(chr_fnname,arglist) 
  } else if (length(grep(pack,packageDescription("spaMM")$Suggests))) {
    ## then the necessary functions cannot be imported-from in the NAMESPACE  (and the package must be written in an appropriate way)
    if ( requireNamespace(pack, quietly = TRUE)) {
      myfun <- get(chr_fnname, asNamespace(pack)) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
      do.call(myfun,arglist) 
    } else {stop(paste0("'",pack,"' required but not available."))}
  } else { ## package not declared in DESCRIPTION; private for example
    if (do.call("require",list(package=pack, quietly = TRUE))) {
      do.call(chr_fnname,arglist) 
    } else {stop(paste0("'",pack,"' required but not available."))}
  }
}

.get_rho_from_SEM_output <- function(SEMblob, lambda.Fix) {
  if (is.null(glm_lambda <- SEMblob$glm_lambda)) {
    rho <- SEMblob$corr_est["rho"] ## may again be NULL
  } else {
    # parallels a block of code from .calcRanefPars(), which is not run in SEM case 
    coeffs <- coefficients(glm_lambda)
    if ( is.na(coeffs["adjd"])) {
      rho <- SEMblob$corr_est["rho"] ## may again be NULL
    } else {
      if (is.na(lambda.Fix[1L])) { # [1L] is ad hoc not spaMM 3.0
        rho <- - coeffs[["adjd"]]/ coeffs[["(Intercept)"]]
      } else {
        rho <- - coeffs[1]*lambda.Fix
      }
    }
  }
  return(rho) 
}

