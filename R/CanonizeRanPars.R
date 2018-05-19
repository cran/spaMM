## the "type" attribute is important, see comments explaining the computation of $CorrEst_and_RanFix
## There is also an "init.HLfit" attribute which usage is explained below for $rho
.canonizeRanPars <- function(ranPars, ## should have a RHOMAX attribute when trRho in input
                             corr_types,checkComplete=TRUE) {
  init.HLfit <- list() 
  for (it in seq_along(corr_types)) {
    corr_type <- corr_types[it]
    if (! is.na(corr_type)) {
      char_rd <- as.character(it)
      corrPars_rd <- ranPars$corrPars[[char_rd]]
      cP_type_rd <- attr(ranPars,"type")$corrPars[[char_rd]]
      if (corr_type == "AR1") {
        ARphi <- corrPars_rd$ARphi
        if (is.null(ARphi) && checkComplete) {
          stop("ARphi missing from ranPars.")
        }
      } else if (corr_type != "corrMatrix") { ## all models with a 'rho' parameter
        moreargs_rd <- attr(ranPars,"moreargs")[[char_rd]]
        ## nu, before trRho is removed; and Nugget: 
        if (corr_type == "Matern") {
          if (!is.null(corrPars_rd$trNu)) { ## either we have nu,rho or trNu,trRho 
            corrPars_rd$nu <- .nuInv(corrPars_rd$trNu,NUMAX=moreargs_rd$NUMAX) ## before trRho is removed...
            corrPars_rd$trNu <- NULL
            cP_type_rd$nu <- cP_type_rd$trNu 
            cP_type_rd$trNu <- NULL
          } 
          nu <- corrPars_rd$nu
          if (is.null(nu) && checkComplete) {
            stop("nu missing from ranPars (or correlation model mis-identified).")
          }
        } else if (corr_type == "Cauchy") {
          if (!is.null(corrPars_rd$trLongdep)) { ## either we have longdep,rho or trLongdep,trRho 
            corrPars_rd$longdep <- .longdepInv(corrPars_rd$trLongdep,LDMAX=moreargs_rd$LDMAX)
            corrPars_rd$trLongdep <- NULL
            cP_type_rd$longdep <- cP_type_rd$trLongdep 
            cP_type_rd$trLongdep <- NULL
          }
          longdep <- corrPars_rd$longdep
          if (is.null(longdep) && checkComplete) {
            stop("longdep missing from ranPars (or correlation model mis-identified).")
          }
        }
        ## rho, Matern or not Matern:
        if ( ! is.null(corrPars_rd$trRho)) { ## assuming a single trRho with possibly several elements
          corrPars_rd$rho <- .rhoInv(corrPars_rd$trRho,RHOMAX=moreargs_rd$RHOMAX)  
          corrPars_rd$trRho <- NULL
          cP_type_rd$rho <- cP_type_rd$trRho
          cP_type_rd$trRho <- NULL
        } ## else there may simply be rho rather than trRho (including for adjacency model through optim procedure !)
        rho <- corrPars_rd$rho
        if (is.null(rho)) {
          if(corr_type=="adjacency") { ## then provide initial rho to allow a direct call through HLCor 
            init.HLfit$corrPars[[char_rd]] <- list(rho=0)
          } else if (checkComplete) {
            stop("rho missing from ranPars.")
          }
        }
      }
      ranPars$corrPars[[char_rd]] <- corrPars_rd
      attr(ranPars,"type")$corrPars[[char_rd]] <- cP_type_rd
    }
  }
  if (!is.null(ranPars$trPhi)) {
    ranPars$phi <- .dispInv(ranPars$trPhi)
    ranPars$trPhi <- NULL
    attr(ranPars,"type")$phi <- attr(ranPars,"type")$trPhi
    attr(ranPars,"type")$trPhi <- NULL
  } # else ranPars$phi unchanged
  # ranPars may have trLambda and lambda (eg fitme(...lambda=c(<value>,NA)))  
  # It may have $trLambda (from notlambda) for what is optimized,
  #              and $lambda (from ranPars$lambda) for what was fixed in the whole outer fit, and also ini.value  
  if ( ! is.null(ranPars$trLambda)) {## 
    lambda <- ranPars$lambda
    if (is.null(lambda)) { ## only trLambda, not lambda
      ranPars$lambda <- .dispInv(ranPars$trLambda)
      type <- attr(ranPars,"type")$trLambda 
    } else { ## merge lambda and trLambda
      len_lam <- seq(length(lambda))
      type <- attr(ranPars,"type")$lambda
      if (is.null(names(lambda))) names(lambda) <- len_lam ## but do not try to assign a vector of names for a single 'type' value
      fromTr <- .dispInv(ranPars$trLambda)
      lambda[names(fromTr)] <- fromTr  
    }
    attr(ranPars,"type")$lambda <- type
    ranPars$trLambda <- NULL
    attr(ranPars,"type")$trLambda <- NULL
  } ## else ranPars$lambda unchanged  
  if ( ! is.null(ranPars$trRanCoefs)) {
    ranPars$ranCoefs <- lapply(ranPars$trRanCoefs,.ranCoefsInv)
    ranPars$trRanCoefs <- NULL
  }
  if ( ! is.null(ranPars$trNB_shape)) {
    ranPars$NB_shape <- .NB_shapeInv(ranPars$trNB_shape)
    ranPars$trNB_shape <- NULL
  }
  attr(ranPars,"init.HLfit") <- init.HLfit  
  return(ranPars)
}
