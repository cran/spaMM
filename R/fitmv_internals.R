# presumably not used:
.fix_LowUp <- function(LowUp, user.lower,user.upper, moreargs, corr_types, init.optim) {
  lower <- LowUp$lower ## list ! which elements may have length >1 !
  upper <- LowUp$upper ## list !
  for (it in seq_along(corr_types)) {
    corr_type <- corr_types[[it]]
    if (! is.na(corr_type)) {
      char_rd <- as.character(it)
      if (corr_type %in% c("SAR_WWt","adjacency") && 
          ! is.null(.get_cP_stuff(init.optim,"rho",char_rd)) ## to exclude inner optimization (was not previously necessary bc 
          # sequence was calc.inits / add rhorange conditionally on inits$init.optim / compute LowUp)
          # while now this this is calc rhorange before calc.inits )
      ) { ## adjacency model
        rhorange <- moreargs[[char_rd]]$rhorange  
        eps <- (rhorange[2L]-rhorange[1L])/(2e6)  
        ## no transfo for adjacency model
        if (is.null(lower$corrPars[[char_rd]][["rho"]] <- .get_cP_stuff(user.lower,"rho",char_rd))) lower$corrPars[[char_rd]][["rho"]] <- rhorange[1L]+eps ## may remain NULL  
        ## no transfo again:
        if (is.null(upper$corrPars[[char_rd]][["rho"]] <- .get_cP_stuff(user.upper,"rho",char_rd))) upper$corrPars[[char_rd]][["rho"]] <- rhorange[2L]-eps
      }
    }
  }
  LowUp$lower <- lower
  LowUp$upper <- upper
  return(LowUp)
}


.makeLowUp_stuff_mv <- function(optim_blob, user.lower, user.upper, optim.scale, processed, verbose) {
  init.optim <- optim_blob$inits$`init.optim` ## list; subset of all estimands, as name implies, and in transformed scale
  LUarglist <- list(canon.init=optim_blob$inits$`init`, # this has template values for all required residual dispersion parameters (phi's, negbin scale's, COMPnu's):
                    # one or more values depending on submodel families.  
                    init.optim=init.optim,
                    user.lower=user.lower,user.upper=user.upper, # those directly in the fitmv_body() call
                    corr_types=processed$corr_info$corr_types,
                    ranFix=optim_blob$fixed, # inits$ranFix, # Any change in $ranFix would be ignored 
                    optim.scale=optim.scale) 
  LUarglist$moreargs <- .calc_moreargs(processed=processed, 
                                       corr_types=processed$corr_info$corr_types, fixed=optim_blob$fixed, init.optim=init.optim, 
                                       control_dist=processed$control_dist, NUMAX=50, LDMAX=50, 
                                       KAPPAMAX=100.000001, # so that users can set it to 100...
                                       init.HLfit=optim_blob$inits$`init.HLfit`, corr_info=processed$corr_info, verbose, 
                                       lower=user.lower, upper=user.upper) 
  optim_blob$LUarglist <- LUarglist
  optim_blob$LowUp <- do.call(".makeLowerUpper",LUarglist)
  optim_blob
}

