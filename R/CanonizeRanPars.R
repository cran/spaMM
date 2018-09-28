## the "type" attribute is important, see comments explaining the computation of $CorrEst_and_RanFix
## There is also an "init.HLfit" attribute which usage is explained below for $rho
.canonizeRanPars <- function(ranPars, ## should have a RHOMAX attribute when trRho in input
                             corr_info, ## NULL in HLfit fns assuming no corr_types is processed there. Else processed$corr_info or <HLfit object>$sub_corr_info
                             checkComplete=TRUE) {
  init.HLfit <- list() 
  corr_types <- corr_info$corr_types
  for (rd in seq_along(corr_types)) {
    corr_type <- corr_types[rd]
    if (! is.na(corr_type)) {
      char_rd <- as.character(rd)
      canonizeblob <- corr_info$corr_families[[rd]]$canonize(corrPars_rd=ranPars$corrPars[[char_rd]],
                                                             cP_type_rd=attr(ranPars,"type")$corrPars[[char_rd]], 
                                                             checkComplete=checkComplete,
                                                             moreargs_rd = attr(ranPars,"moreargs")[[char_rd]])  
      ranPars$corrPars[[char_rd]] <- canonizeblob$corrPars_rd
      if (!is.null(canonizeblob$cP_type_rd)) attr(ranPars,"type")$corrPars[[char_rd]] <- canonizeblob$cP_type_rd
      if(corr_type=="adjacency" && is.null(canonizeblob$corrPars_rd$rho)) { ## then provide initial rho to allow a direct call through HLCor 
        init.HLfit$corrPars[[char_rd]] <- list(rho=0)
      }
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
  if ( ! is.null(trRanCoefs <- ranPars$trRanCoefs)) {
    ranPars$ranCoefs <- trRanCoefs ## copies the non-trivial names
    for (rd in seq_along(trRanCoefs)) ranPars$ranCoefs[[rd]] <- .ranCoefsInv(trRanCoefs[[rd]])
    ranPars$trRanCoefs <- NULL
  }
  if ( ! is.null(ranPars$trNB_shape)) {
    ranPars$NB_shape <- .NB_shapeInv(ranPars$trNB_shape)
    ranPars$trNB_shape <- NULL
  }
  attr(ranPars,"init.HLfit") <- init.HLfit  
  return(ranPars)
}
