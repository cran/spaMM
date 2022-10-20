
get_infomat <- function(object, # currenty assumed to be fitme() result
                         verbose=TRUE,
                         ...
                         ) {
  iff <- get_inits_from_fit(object, inner_lambdas = TRUE)
  skeleton <- list(fixef=na.omit(fixef(object)))
  # Avoid explicit NULLs:
  skeleton$loglambda <- log(c(iff[["init.HLfit"]]$lambda,iff[["init"]]$lambda))
  skeleton$corrPars <- iff[["init"]]$corrPars
  
  
  if (verbose) print(skeleton)
  update_wrap <- function(x) {
    pt <- relist(x, skeleton)
    if (verbose) cat(".")
    logLik(update(object, 
                  fixed=list(lambda=exp(pt$loglambda),
                             corrPars=pt$corrPars), # => *fitme* results
                  etaFix=list(beta=pt$fixef)))
  }
  infomat <- - hessian(func=update_wrap, x=unlist(skeleton), ...) # ___F I X M E___ rethink this...
  if (verbose) cat("\n")
  rownames(infomat) <- colnames(infomat) <- names(unlist(skeleton)) 
  infomat # a matrix
}

