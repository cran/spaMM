.get_refit_args <- function(fixed, optPars, processed, moreargs, proc1, refit_info, HLCor.args, augZXy_phi_est) {
  ufixed <- unlist(fixed)
  fixtyp <- rep("fix",length(ufixed)) 
  # special handling for partially-fixed ranCoefs:
  fixtyp[is.na(ufixed)] <- "outer"
  # so that ranPars_in_refit will have apparently conflicting info such as "outer" "fix"   "outer" on $ranCoefs$`1` and "outer" "outer" "outer" on $trRanCoefs$`1`
  # so we use the presence of both ranCoefs and trRanCoefs to distinguish the user setting and the internal optim over 3 params
  ranPars_in_refit <- structure(.modify_list(fixed,optPars),
                                type=.modify_list(relist(fixtyp,fixed), #attr(fixed,"type"),
                                                  relist(rep("outer",length(unlist(optPars))),optPars)))
  ranPars_in_refit <- .expand_hyper(ranPars_in_refit, processed$hyper_info,moreargs=moreargs)
  
  #any_nearly_singular_covmat <- FALSE
  if (! is.null(optPars$trRanCoefs)) {
    constraints <- fixed$ranCoefs # hack that uses the presence of the constraint here, jointly with $trRanCoefs
    ranCoefs <- optPars$trRanCoefs # copy names...
    for (char_rd in names(ranCoefs)) {
      ranCoefs[[char_rd]] <- 
        .constr_ranCoefsInv(trRanCoef=ranCoefs[[char_rd]], constraint=constraints[[char_rd]], rC_transf=.spaMM.data$options$rC_transf)
      covmat <- .C_calc_cov_from_ranCoef(ranCoefs[[char_rd]])
      #if (kappa(covmat)>1e14 || min(eigen(covmat,only.values = TRUE)$values)<1e-6) any_nearly_singular_covmat <- TRUE # use of this removed 2019/12/16
    }
  }
  init_refit <- list()
  if ( is.null(refit_info) ) { ## not the case with SEM
    refit_info <- list(phi=FALSE,lambda=(! is.null(optPars$trRanCoefs)),ranCoefs=FALSE) # defines default for refit_info
  } else if ( is.list(refit_info) ) { ## never the default
    refit_info <- .modify_list( list(phi=FALSE, lambda=(! is.null(optPars$trRanCoefs)), ranCoefs=FALSE), 
                                refit_info) 
    ## lambda=TRUE becomes the default (test random-slope 'ares' shows the effect)
  } else {
    refit_info <- list(phi=refit_info,lambda=refit_info,ranCoefs=refit_info) # default with SEM is all FALSE
  }
  # refit: (change ranPars_in_refit depending on what's already in and on refit_info)
  if (identical(refit_info$phi,TRUE)) {
    if ( ! is.null(optPars$trPhi)) {
      init_refit$phi <- .dispInv(optPars$trPhi)
      ranPars_in_refit$trPhi <- NULL
      attr(ranPars_in_refit,"type")$trPhi <- NULL 
    } else if (proc1$augZXy_cond ) {
      init_refit$phi <- augZXy_phi_est # might still be NULL...
    }
  } else if (proc1$augZXy_cond &&  ! is.null(augZXy_phi_est)) ranPars_in_refit$trPhi <- .dispFn(augZXy_phi_est)
  # refit, or rescale by augZXy_phi_est even if no refit:
  if ( ! is.null(augZXy_phi_est) || identical(refit_info$lambda,TRUE)) {
    # For outer optim with augZXy algo, we have fitted for the model (lambda, 1/prior_weights) and deduced the optimal (lamphifac*lambda, lamphifac/prior_weights)
    # this lamphifac being stored as augZXy_phi_est. We then correct the lambda (only now in all cases: this was outer-optim, the objfn returned almost only likelihoods)
    # We will further correct the ranCoefs and hy_lam below.
    if (length(optPars$trLambda)) { # does NOT include the lambdas of the ranCoefs
      lambdapos <- as.integer(names(optPars$trLambda))
      lambda <- rep(NA,max(lambdapos))
      names(lambda) <- paste(seq_len(length(lambda)))
      lambda[lambdapos] <- .dispInv(optPars$trLambda)
      if ( ! is.null(augZXy_phi_est)) { lambda <- lambda * augZXy_phi_est }
      if (identical(refit_info$lambda,TRUE)) {
        init_refit$lambda <- lambda
        ranPars_in_refit$trLambda <- NULL
        attr(ranPars_in_refit,"type")$trLambda <- NULL
      } else ranPars_in_refit$trLambda <- .dispFn(lambda) ## rescale, but no refit
    }
  }
  # rescaling hyper-controlled lambdas
  if ( ! is.null(augZXy_phi_est) && ! is.null(optPars$hyper)) {
    ranges <- processed$hyper_info$ranges #  attr(optPars$hyper,"hy_info")$ranges : same environment
    for (rg_it in names(ranges)) {
      hyper_el <- ranPars_in_refit$hyper[[rg_it]]
      if ( ! is.null(hy_trL <- hyper_el$hy_trL)) { # estimated ones
        ranPars_in_refit$hyper[[rg_it]]$hy_trL <- .dispFn(.dispInv(hy_trL)*augZXy_phi_est) # no direct effect on refit
        char_rd_range <- as.character(ranges[[rg_it]]) 
        trL <- ranPars_in_refit$trLambda[char_rd_range]
        ranPars_in_refit$trLambda[char_rd_range] <- .dispFn(.dispInv(trL)*augZXy_phi_est) # the effective rescaling
      }
    } 
  }
  ## refit, or rescale by augZXy_phi_est even if no refit:
  if ( ! is.null(augZXy_phi_est)  || identical(refit_info$ranCoefs,TRUE)) {
    if (! is.null(optPars$trRanCoefs)) {
      for (char_rd in names(ranCoefs)) {
        rancoef <- ranCoefs[[char_rd]] # full size even for DiagFamily fit; at this point the vector no longer has NA so isDiagFamily has been set to FALSE
        if ( ! is.null(augZXy_phi_est)) {
          Xi_ncol <- attr(rancoef,"Xi_ncol")
          lampos <- rev(length(rancoef) -cumsum(seq(Xi_ncol))+1L)  ## NOT cumsum(seq(Xi_cols))
          rancoef[lampos] <- rancoef[lampos] *augZXy_phi_est
          rancoef <- as.vector(rancoef) ## as.vector drops attributes
          attr(rancoef,"isDiagFamily") <- FALSE
        } 
        #  Note that both ranCefs and trRanCoefs may be used to count dfs so they are not quite redundant
        if (identical(refit_info$ranCoefs,TRUE)) { 
          init_refit$ranCoefs[[char_rd]] <- rancoef
          ranPars_in_refit$trRanCoefs[char_rd] <- NULL
          attr(ranPars_in_refit,"type")$trRanCoefs[char_rd] <- NULL
        } else { # ranCoefs are fixed in the refit.
          ranPars_in_refit$ranCoefs[[char_rd]] <- rancoef
          ranPars_in_refit$trRanCoefs[[char_rd]] <- .ranCoefsFn(rancoef, rC_transf=.spaMM.data$options$rC_transf_inner)  
        }
      }
    }      
  }
  if (length(init_refit)) HLCor.args$init.HLfit <- .modify_list(HLCor.args$init.HLfit, init_refit) 
  return(list(HLCor.args=HLCor.args, ranPars_in_refit=ranPars_in_refit))
} 


