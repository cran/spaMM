## better for development to avoid name conflicts with OKsmooth :toCanonical and :canonize
.canonizeRanPars <- function(ranPars, ## should have a RHOMAX attribute when trRho in input
                            corr_types,checkComplete=TRUE) {
  trueCorrpars <- list()
  for (it in seq_along(corr_types)) {
    corr_type <- corr_types[it]
    if (! is.na(corr_type)) {
      if (corr_type == "AR1") {
        ARphi <- ranPars$ARphi
        if (is.null(ARphi) && checkComplete) {
          stop("ARphi missing from ranPars.")
        }
        trueCorrpars$ARphi <- ARphi    
      } else if (corr_type != "corrMatrix") { ## all models with a 'rho' parameter
        ## nu, before trRho is removed; and Nugget: 
        if (corr_type == "Matern") {
          if (!is.null(ranPars$trNu)) { ## either we have nu,rho or trNu,trRho 
            ranPars$nu <- .nuInv(ranPars$trNu,ranPars$trRho,NUMAX=attr(ranPars,"NUMAX")) ## before trRho is removed...
            ranPars$trNu <- NULL
            attr(ranPars,"type")$nu <- attr(ranPars,"type")$trNu
            attr(ranPars,"type")$trNu <- NULL
          } 
          nu <- ranPars$nu
          if (is.null(nu) && checkComplete) {
            stop("nu missing from ranPars (or correlation model mis-identified).")
          }
          trueCorrpars$nu <- nu 
          Nugget <- ranPars$Nugget
          if (! is.null(Nugget)) trueCorrpars$Nugget <- Nugget 
        } 
        ## rho, Matern or not Matern:
        if ( ! is.null(ranPars$trRho)) { ## assuming a single trRho with possibly several elements
          ranPars$rho <- .rhoInv(ranPars$trRho,RHOMAX=attr(ranPars,"RHOMAX"))  
          ranPars$trRho <- NULL
          attr(ranPars,"type")$rho <- attr(ranPars,"type")$trRho
          attr(ranPars,"type")$trRho <- NULL
        } ## else there may simply be rho rather than trRho (including for adjacency model through optim procedure !)
        trueCorrpars$rho <- rho <- ranPars$rho
        if (is.null(rho)) {
          if(corr_type=="adjacency") { ## then provide initial rho to allow a direct call through HLCor 
            ranPars$rho <- 0  ## Should spaMM 3.0 correct this ? Direct call through HLCor may exclude other rho terms
            attr(ranPars,"type")$rho <- "var"
          } else if (checkComplete) {
            stop("rho missing from ranPars.")
          }
        }
      }
    }
  }
  if (!is.null(ranPars$trPhi)) {
    ranPars$phi <- .dispInv(ranPars$trPhi)
    ranPars$trPhi <- NULL
    attr(ranPars,"type")$phi <- attr(ranPars,"type")$trPhi
    attr(ranPars,"type")$trPhi <- NULL
    if (spaMM.getOption("wDEVEL2")) {
      parlist <- attr(ranPars,"parlist")
      parlist$phi <- .dispInv(parlist$trPhi)
      parlist$trPhi <- NULL
      attr(parlist,"types")$phi <- attr(parlist,"types")$trPhi
      attr(parlist,"types")$trPhi <- NULL
      attr(ranPars,"parlist") <- parlist
    }
  } # else ranPars$phi unchanged
  # ranPars may have trLambda and lambda (eg fitme(...lambda=c(<value>,NA)))  
  # It may have $trLambda (from notlambda) for what is optimized,
  #              and $lambda (from ranPars$lambda) for what was fixed in the whole outer fit, and also ini.value  
  if ( ! is.null(ranPars$trLambda)) {## 
    lambda <- ranPars$lambda
    ## At this point Fix and init.HLfit are merged 
    #  and the type info will be used by HLCor_body to separate what comes from ranPars ("fix") and what comes from init.HLfit ("var")
    #  Currently e.g. lambda can be of only one type
    if (is.null(lambda)) { ## only trLambda, not lambda
      ranPars$lambda <- .dispInv(ranPars$trLambda)
      type <- attr(ranPars,"type")$trLambda 
    } else { ## merge lambda and trLambda
      len_lam <- seq(length(lambda))
      type <- attr(ranPars,"type")$lambda # presumably "fix", or "var"<=> from init.HLfit 
      if (is.null(names(lambda))) names(lambda) <- len_lam ## but do not try to assign a vector of names for a single 'type' value
      fromTr <- .dispInv(ranPars$trLambda)
      lambda[names(fromTr)] <- fromTr  
      #type[names(fromTr)] <- attr(ranPars,"type")$trLambda ## presumably "fix" (for fully fix or outer estimated)
      ranPars$lambda <- lambda
    }
    attr(ranPars,"type")$lambda <- type
    ranPars$trLambda <- NULL
    attr(ranPars,"type")$trLambda <- NULL
    if (spaMM.getOption("wDEVEL2")) {
      parlist <- attr(ranPars,"parlist")
      lambda <- parlist$lambda
      if (is.null(lambda)) { ## only trLambda, not lambda
        types <- attr(parlist,"types")$trLambda 
        parlist$lambda <- .dispInv(parlist$trLambda)
      } else { ## merge lambda and trLambda
        len_lam <- seq(length(lambda))
        types <- attr(parlist,"types")$lambda # presumably "fix", or "var"<=> from init.HLfit 
        if (is.null(names(lambda))) names(lambdas) <- len_lam 
        fromTr <- .dispInv(parlist$trLambda)
        lambda[names(fromTr)] <- fromTr  
        types[names(fromTr)] <- attr(parlist,"types")$trLambda ## presumably "fix" (for fully fix or outer estimated)
        parlist$lambda <- lambda
      }
      parlist$trLambda <- NULL
      attr(parlist,"types")$lambda <- types
      attr(parlist,"types")$trLambda <- NULL
      attr(ranPars,"parlist") <- parlist
    }
  } ## else ranPars$lambda unchanged  
  if ( ! is.null(ranPars$trRanCoefs)) {
    ranPars$ranCoefs <- lapply(ranPars$trRanCoefs,.ranCoefsInv)
    ranPars$trRanCoefs <- NULL
  }
  if ( ! is.null(ranPars$trNB_shape)) {
    ranPars$NB_shape <- .NB_shapeInv(ranPars$trNB_shape)
    ranPars$trNB_shape <- NULL
  }
  return(list(trueCorrpars=trueCorrpars,ranPars=ranPars))
}

