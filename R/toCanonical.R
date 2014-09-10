toCanonical <-
function(ranPars,corr.model,checkComplete=TRUE) {
  trueCorrpars <- list()
  if (corr.model %in% c("Matern","corMatern")) {
    if (!is.null(ranPars$trNu)) { ## either we have nu,rho or trNu,trRho 
      ranPars$nu <- nuInv(ranPars$trNu,ranPars$trRho) ## before trRho is removed...
      ranPars$trNu <- NULL
      attr(ranPars,"type")$nu <- attr(ranPars,"type")$trNu
      attr(ranPars,"type")$trNu <- NULL
    } 
    nu <- ranPars$nu
    if (is.null(nu) && checkComplete) {
      mess <- pastefrom("nu missing from ranPars (or correlation model mis-identified).",prefix="(!) From ")
      stop(mess)
    }
    trueCorrpars$nu <- nu 
  } 
  if (corr.model=="AR1") {
    ARphi <- ranPars$ARphi
    if (is.null(ARphi) && checkComplete) {
      mess <- pastefrom("ARphi missing from ranPars.",prefix="(!) From ")
      stop(mess)
    }
    trueCorrpars$ARphi <- ARphi    
  } else if (corr.model != "corrMatrix") { ## all models with a 'rho' parameter
    if (!is.null(ranPars$trRho)) {
      ranPars$rho <- rhoInv(ranPars$trRho)
      ranPars$trRho <- NULL
      attr(ranPars,"type")$rho <- attr(ranPars,"type")$trRho
      attr(ranPars,"type")$trRho <- NULL
    } ## else there may simply be rho rather than trRho (including for adjacency model through optim procedure !)
    rho <- ranPars$rho
    if (is.null(rho) && checkComplete) {
      mess <- pastefrom("rho missing from ranPars.",prefix="(!) From ")
      stop(mess)
    }
    trueCorrpars$rho <- rho
  }
  Nugget <- ranPars$Nugget
  if (! is.null(Nugget)) trueCorrpars$Nugget <- Nugget 
  if (!is.null(ranPars$trPhi)) {
    ranPars$phi <- dispInv(ranPars$trPhi)
    ranPars$trPhi <- NULL
    attr(ranPars,"type")$phi <- attr(ranPars,"type")$trPhi
    attr(ranPars,"type")$trPhi <- NULL
  } else if (!is.null(ranPars$logphi)) { ## debug code
    ## HL.info$ranFix$phi <- exp(ranPars$logphi)
    stop("logphi in HLCor...")
  } #####################  else HL.info$ranFix$phi <- ranPars$phi ## y st deja !?
  if (!is.null(ranPars$trLambda)) {## 
    ranPars$lambda <- dispInv(ranPars$trLambda)
    ranPars$trLambda <- NULL
    attr(ranPars,"type")$lambda <- attr(ranPars,"type")$trLambda
    attr(ranPars,"type")$trLambda <- NULL
  } else if (!is.null(ranPars$loglambda)) { ## debug code
    stop("loglambda in HLCor...")
  } ##################### else HL.info$ranFix$lambda <- ranPars$lambda
  return(list(trueCorrpars=trueCorrpars,ranPars=ranPars))
}
