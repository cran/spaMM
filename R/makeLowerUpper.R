.get_rdisPars_LowUp <- function(disp_env, 
                                X={ if (is.null(XX <- disp_env$X)) {XX <- disp_env$scaled_X} ; XX }, 
                                off=disp_env$off,lo, hi) { # lo and hi on 'response' scale for dispersion param
  if (length(off)==1L) off <- rep(off,nrow(X))
  # makeH represent the constraints a1 %*% beta <= b1
  # here X beta < log(hi) -off and [X beta >= log(lo) -off => - X beta <= off- log(lo) ]
  if (requireNamespace("rcdd",quietly=TRUE)) {
    vrepr <- rcdd::scdd(rcdd::makeH(a1=rbind(X,-X),b1=c(log(hi)-off,off-log(lo))))$output ## convex hull of feasible coefs
    if (nrow(vrepr)==0L) {
      warning("The dispersion model appear to impose extreme dispersion values so may be difficult to fit.")
      list(lower=rep(log(.Machine$double.xmin),ncol(X)), upper=rep(log(.Machine$double.xmax),ncol(X)))
    } else {
      verticesRows <- (vrepr[,2]==1)
      betaS <- vrepr[verticesRows, -c(1:2), drop = FALSE] ## row vectors of boundary beta values
      # The convex hull cannot be directly used. Hopefully the enclosing rectangle is good enough:
      list(lower=apply(betaS, 2L, min),
           upper=apply(betaS, 2L, max))
    }
  } else if ( ! identical(spaMM.getOption("rcdd_warned"),TRUE)) {
    message("If the 'rcdd' package were installed, spaMM could find a good range for the family dispersion parameters.")
    .spaMM.data$options$rcdd_warned <- TRUE
    list(lower=rep(log(.Machine$double.xmin),ncol(X)), upper=rep(log(.Machine$double.xmax),ncol(X)))
  }
}


.wrap_calc_famdisp_lowup <- function(processed, family=processed$family, prior.weights=processed$prior.weights) {
  
  # provide famdisp_lowup bc LUarglist, which is returned in the fit object, should not include 'processed'
  lo <- switch(family$family,
               "COMPoisson" = 0.05, # no prior.weights handling for COMPoisson 
               "beta_resp" = 1e-6/prior.weights, # '/'pw bc the disp param is here a prec param: precision =prec*pw must be within 1e-6, 1e6
               1e-6) #.Machine$double.xmin) # 1e-6 is typical lower bound for NB_shape
  hi <- switch(family$family,
               "COMPoisson" = 10, # no prior.weights handling for COMPoisson 
               "beta_resp" = 1e6/prior.weights,
               1e6) # .Machine$double.xmax) # 1e6 is typical upper bound for NB_shape
  .get_rdisPars_LowUp(disp_env=family$resid.model, lo=lo, hi=hi)
}



.calc_fixef_lowup <- function(processed) {
  eta_range <- .sanitize_eta(eta=c(-Inf,Inf),family = processed$family)
  lo <- eta_range[1]
  hi <- eta_range[2]
  X <- environment(processed$X_off_fn)$X_off # should be the scaled version. rcdd is not robust to extreme values
  off <- environment(processed$X_off_fn)$ori_off
  # makeH represent the constraints a1 %*% beta <= b1
  # here X beta < (hi) -off and [X beta >= (lo) -off => - X beta <= off- (lo) ]
  if (requireNamespace("rcdd",quietly=TRUE)) {
    vrepr <- rcdd::scdd(rcdd::makeH(a1=rbind(X,-X),b1=c((hi)-off,off-(lo))))$output ## convex hull of feasible coefs
    if (nrow(vrepr)==0L) {
      warning("The fixed-effect model appears to impose extreme 'eta' values so may be difficult to fit.")
      NULL
    } else {
      verticesRows <- (vrepr[,2]==1)
      betaS <- vrepr[verticesRows, -c(1:2), drop = FALSE] ## row vectors of boundary beta values
      colnames(betaS) <- colnames(X)
      # The convex hull cannot be directly used. Hopefully the enclosing rectangle is good enough:
      list(lower=apply(betaS, 2L, min),
           upper=apply(betaS, 2L, max))
    }
  } else NULL
}


.makeLowerUpper <- function(canon.init, ## cf calls: ~ in user scale, must be a full list of relevant params
                            init.optim, ## ~in transformed scale : is has all pars to be optimized
                            user.lower=list(),user.upper=list(),
                            corr_types=NULL, ranFix=list(),
                            optim.scale, moreargs=NULL, rC_transf=.spaMM.data$options$rC_transf,
                            famdisp_lowup=NULL,
                            fixef_lowup=NULL) {
  lower <- upper <- init.optim   
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
      } else if (corr_type =="AR1") {
        if ( ! is.null(.get_cP_stuff(canon.init,"ARphi",char_rd))) {
          ARphi <- .get_cP_stuff(user.lower,"ARphi",char_rd)
          if (is.null(ARphi)) ARphi <- -1 + 1e-6
          lower$corrPars[[char_rd]][["ARphi"]] <- ARphi
          ARphi <- .get_cP_stuff(user.upper,"ARphi",char_rd)
          if (is.null(ARphi)) ARphi <- 1 - 1e-6
          upper$corrPars[[char_rd]][["ARphi"]] <- ARphi
        }
      } else if (corr_type =="IMRF") {
        KAPPAMAX <- moreargs[[char_rd]]$KAPPAMAX
        hyper <- init.optim$hyper
        hyper_map <- attr(hyper,"hy_info")$map # hy_info is environment, same as in 'processed', the later which is not itself locally accessible
        trKappa <- hyper[[hyper_map[it]]]$hy_trK 
        if (is.null(trKappa)) { ## NOT hyper (independent IMRF terms)    or   fixed hyper
          if ( ! is.null(.get_cP_stuff(canon.init,"kappa",char_rd))) { ## optimized, independent IMRF term
            kappa <- .get_cP_stuff(user.lower,"kappa",char_rd)
            if (is.null(kappa)) {
              kappa <- moreargs[[char_rd]]$minKappa
              if (is.null(kappa)) kappa <- 1e-4
            } else {
              minKappa <- moreargs[[char_rd]]$minKappa
              if ( ! is.null(minKappa) && kappa<minKappa) {
                warning("User-provided minimum value of kappa is low. Numerical errors may result.",
                        immediate.=TRUE)
              }
            }
            if (optim.scale=="transformed") {
              lower$corrPars[[char_rd]][["trKappa"]] <- .kappaFn(kappa,KAPPAMAX=KAPPAMAX)
            } else lower$corrPars[[char_rd]][["kappa"]] <- kappa
            kappa <- .get_cP_stuff(user.upper,"kappa",char_rd)
            if (is.null(kappa)) kappa <- KAPPAMAX
            kappa <- min(KAPPAMAX-1e-6, kappa) # optimization should not try to approach infinity
            if (optim.scale=="transformed") {
              upper$corrPars[[char_rd]][["trKappa"]] <- .kappaFn(kappa,KAPPAMAX=KAPPAMAX)
            } else upper$corrPars[[char_rd]][["kappa"]] <- kappa
          } ## else fixed : do nothing
        } else { ## multIMRF hyper, optimized
          # slightly inelegant as repeated for several it... (same for lambda)
          lower$hyper[[hyper_map[it]]]$hy_trK <- .kappaFn(1e-4,KAPPAMAX=KAPPAMAX)
          upper$hyper[[hyper_map[it]]]$hy_trK <- .kappaFn(KAPPAMAX-1e-6,KAPPAMAX=KAPPAMAX) # .kappaFn(KAPPAMAX,.) is Inf => optimize() stops
        }
        trLambda <- hyper[[hyper_map[it]]]$hy_trL
        if ( ! is.null(trLambda)) { ## multIMRF hyper, optimized
          hy_lam <- user.lower$hyper[[hyper_map[it]]]$hy_lam
          if (is.null(hy_lam)) hy_lam <- 1e-6
          lower$hyper[[hyper_map[it]]]$hy_trL <- .dispFn(max(.dispInv(trLambda)/1e4,hy_lam))
          #
          hy_lam <- user.upper$hyper[[hyper_map[it]]]$hy_lam
          if (is.null(hy_lam)) hy_lam <- 1e7
          upper$hyper[[hyper_map[it]]]$hy_trL <- .dispFn(min(.dispInv(trLambda)*1e4,hy_lam))
        } # else see general code for lambda or all corr_type's 
        ## (which means that there are canon.init$lambda slots in and only in case of independent optimized IMRFs)
      } else if (corr_type %in% c("Matern","Cauchy")) { 
        lower_cP <- lower$corrPars[[char_rd]]
        if (is.null(lower_cP)) lower_cP <- list()
        upper_cP <- upper$corrPars[[char_rd]]
        if (is.null(upper_cP)) upper_cP <- list()
        if (! is.null(canon_rho <- .get_cP_stuff(canon.init,"rho",char_rd))) { # where canon.inits rho typically comes from .calc_inits_geostat_rho()
          RHOMAX <- moreargs[[char_rd]]$RHOMAX
          rho <- .get_cP_stuff(user.lower,"rho",char_rd)
          if (is.null(rho)) rho <- canon_rho/150
          if (optim.scale=="transformed") {
            lower_cP$trRho <- .rhoFn(rho,RHOMAX=RHOMAX)
          } else lower_cP$rho <- rho
          rho <- .get_cP_stuff(user.upper,"rho",char_rd)
          if (is.null(rho)) {
            rho <- canon_rho*2*moreargs[[as.character(it)]]$nbUnique ## The following was a bit too low for experiments with nu=0.5 : 1/(maxrange/(2*nbUnique)) ## nb => unique rows !
            ## *modify* upper rho so that it does not exceed RHOMAX => /($RHOMAX+...)
            if (optim.scale=="transformed") rho <- 2*rho * RHOMAX/(RHOMAX+2*rho)
          }
          if (optim.scale=="transformed") {
            upper_cP$trRho <- .rhoFn(rho,RHOMAX=RHOMAX)
          } else upper_cP$rho <- rho
          rhoForNu <- canon_rho
        } else rhoForNu <- .getPar(ranFix,"rho")
        if (corr_type == "Matern") {
          if (! is.null(canon_nu <- .get_cP_stuff(canon.init,"nu",char_rd))) {
            NUMAX <- moreargs[[char_rd]]$NUMAX
            nu <- .get_cP_stuff(user.lower,"nu",char_rd)
            if (is.null(nu)) nu <- canon_nu/100
            if (optim.scale=="transformed") {
              lower_cP$trNu <- .nuFn(nu,rhoForNu,NUMAX)
            } else lower_cP$nu <- nu
            nu <- .get_cP_stuff(user.upper,"nu",char_rd)
            if (is.null(nu)) {
              control_dist_rd <- moreargs[[char_rd]]$control.dist
              if ( ! is.null(dm <- control_dist_rd$`dist.method`) && dm %in% c("Geodesic","Earth")) {
                nu <- 0.5
              } else {
                ## constructs upper nu from NUMAX => /(1+...)
                ## nu should not diverge otherwise it will diverge in Bessel_lnKnu, whatever the transformation used
                nu <- NUMAX * canon_nu/(1+canon_nu) 
                ## FR->FR hmmm. If canon.init$nu= NUMAX-1 then 
                ##  (upper) nu= canon.init$nu and possibly < canon.init$nu by numerical accuracy issues => nloptr stops
              }
            }
            if (optim.scale=="transformed") {
              upper_cP$trNu <- .nuFn(nu,rhoForNu,NUMAX)
            } else upper_cP$nu <- nu
            #print(c(rhoForNu,nu,upper$trNu))
          }
        } else { ## ,"Cauchy"
          if (! is.null(canon_longdep <- .get_cP_stuff(canon.init,"longdep",char_rd))) {
            LDMAX <- moreargs[[char_rd]]$LDMAX
            longdep <- .get_cP_stuff(user.lower,"longdep",char_rd)
            if (is.null(longdep)) longdep <- canon_longdep/100
            if (optim.scale=="transformed") {
              lower_cP$trLongdep <- .longdepFn(longdep,LDMAX)
            } else lower_cP$longdep <- longdep
            longdep <- .get_cP_stuff(user.upper,"longdep",char_rd)
            if (is.null(longdep)) longdep <- LDMAX * canon_longdep/(1+canon_longdep) 
            if (optim.scale=="transformed") {
              upper_cP$trLongdep <- .longdepFn(longdep,LDMAX)
            } else upper_cP$longdep <- longdep
          }
          if (! is.null(canon_shape <- .get_cP_stuff(canon.init,"shape",char_rd))) {
            shape <- .get_cP_stuff(user.lower,"shape",char_rd)
            if (is.null(shape)) shape <- canon_shape/100
            lower_cP$shape <- shape
            shape <- .get_cP_stuff(user.upper,"shape",char_rd)
            if (is.null(shape)) {
              control_dist_rd <- moreargs[[char_rd]]$control.dist
              if ( ! is.null(dm <- control_dist_rd$`dist.method`) && dm %in% c("Geodesic","Earth")) {
                shape <- 1
              } else {
                shape <- 2
              }
            }
            upper_cP$shape <- shape
          }
        }
        if ( ! is.null(.get_cP_stuff(canon.init,"Nugget",char_rd))) {
          Nugget <- .get_cP_stuff(user.lower,"Nugget",char_rd)
          if (is.null(Nugget)) Nugget <- 0 
          lower_cP$Nugget <- Nugget 
          Nugget <- .get_cP_stuff(user.upper,"Nugget",char_rd)
          if (is.null(Nugget)) Nugget <- 0.999999 
          upper_cP$Nugget <- Nugget 
        }
        lower$corrPars[[char_rd]] <- lower_cP
        upper$corrPars[[char_rd]] <- upper_cP
        ## end else if Matern/Cauchy case
      } else if (corr_type =="corrFamily") { # that looks like a quick patch but part of the work is done elsewhere; .calc_inits_corrFamily() checks the user input, which is required here. 
        parnames <- names(init.optim$corrPars[[char_rd]]) 
        if (length(parnames)) {
          user_lo <- user.lower$corrPars[[char_rd]]
          if (is.null(user_lo)) {
            def_lo <- moreargs[[char_rd]]$lower
            if (is.null(def_lo)) stop("lower bounds provided neither by 'lower' nor by corrFamily's 'calc_moreargs' element.")
            lower$corrPars[[char_rd]] <- def_lo[parnames]
          } else {
            if (is.null(names(user_lo))) names(user_lo) <- parnames
            lower$corrPars[[char_rd]] <- user_lo
          }
          user_up <- user.upper$corrPars[[char_rd]]
          if (is.null(user_up)) {
            def_up <- moreargs[[char_rd]]$upper
            if (is.null(def_up)) stop("upper bounds provided neither by 'upper' nor by corrFamily's 'calc_moreargs' element.")
            upper$corrPars[[char_rd]] <- def_up[parnames]
          } else {
            if (is.null(names(user_up))) names(user_up) <- parnames
            upper$corrPars[[char_rd]] <- user_up
          }
        }
      } 
    }
  }
  
  if (! is.null(canon.init$phi)) { # *p*min, *p*max introduced for vector phi's for fitmv : might affect univariate-response fits
    phi <- user.lower$phi
    if (is.null(phi)) phi <- pmax(pmin(1e-6,canon.init$phi/1.01),canon.init$phi/1e5) # >=1e-6 if canon.init$phi>1e-6
    names(phi) <- names(canon.init$phi) # late addition for mv code (merging inits...)
    lower$trPhi <- .dispFn(phi)
    phi <- user.upper$phi
    if (is.null(phi)) phi <- pmin(pmax(1e8,canon.init$phi*1.01),canon.init$phi*1e7) 
    names(phi) <- names(canon.init$phi) # late addition for mv code
    ## if phi is badly initialized then it gets a default which may cause hard to catch problems in the bootstrap...
    upper$trPhi <- .dispFn(phi)
  }
  if (! is.null(canon.init$lambda)) {
    lambda <- user.lower$lambda
    if (is.null(lambda)) lambda <- pmax(1e-6,canon.init$lambda/1e5)
    names(lambda) <- names(canon.init$lambda) # late addition for mv code
    lower$trLambda <- .dispFn(lambda)
    lambda <- user.upper$lambda
    if (is.null(lambda)) lambda <- pmin(1e8,canon.init$lambda*1e7)
    names(lambda) <- names(canon.init$lambda) # late addition for mv code
    upper$trLambda <- .dispFn(lambda)
  }
  if (! is.null(canon.init$COMP_nu)) {
    if (length(canon.init$COMP_nu)>1L) { # mv case with >1 COMPoisson submodels
      # then all vectors mus be named and canon.init must have values for all COMpoisson submodels
      lower$COMP_nu <- .modify_list(pmin(canon.init$COMP_nu/2,0.05), user.lower$COMP_nu)
      upper$COMP_nu <- .modify_list(pmax(canon.init$COMP_nu*10,10), user.upper$COMP_nu)
    } else {
      if (is.null(COMP_nu <- user.lower$COMP_nu)) COMP_nu <- pmin(canon.init$COMP_nu/2,0.05)
      lower$COMP_nu <- COMP_nu
      if (is.null(COMP_nu <- user.upper$COMP_nu)) COMP_nu <- pmax(canon.init$COMP_nu*10,10)
      upper$COMP_nu <- COMP_nu
    }
  }
  if (! is.null(canon.init$NB_shape)) { ## shape param of latent [Gamma(1:sh,sh) with mean 1 and variance sh]
    if (length(canon.init$NB_shape)>1L) { # mv case with >1 COMPoisson submodels
      # then all vectors mus be named and canon.init must have values for all COMpoisson submodels
      NB_shape <- .modify_list(rep(1e-6, length(canon.init$NB_shape)), user.lower$NB_shape)
      lower$trNB_shape <- .NB_shapeFn(NB_shape)
      NB_shape <- .modify_list(pmax(100*canon.init$NB_shape,1e6), user.upper$NB_shape)
      upper$trNB_shape <- .NB_shapeFn(NB_shape)
    } else {
      if (is.null(NB_shape <- user.lower$NB_shape)) NB_shape <- 1e-6
      lower$trNB_shape <- .NB_shapeFn(NB_shape)
      if (is.null(NB_shape <- user.upper$NB_shape)) NB_shape <- max(100*canon.init$NB_shape,1e6)
      upper$trNB_shape <- .NB_shapeFn(NB_shape)
    }
  }
  if (! is.null(canon.init$beta_prec)) { 
    if (length(canon.init$beta_prec)>1L) { # mv case with >1 beta_resp submodels
      # then all vectors mus be named and canon.init must have values for all beta_resp submodels
      beta_prec <- .modify_list(rep(1e-6, length(canon.init$beta_prec)), user.lower$beta_prec)
      lower$trbeta_prec <- .beta_precFn(beta_prec)
      beta_prec <- .modify_list(pmax(100*canon.init$beta_prec,1e6), user.upper$beta_prec)
      upper$trbeta_prec <- .beta_precFn(beta_prec)
    } else {
      if (is.null(beta_prec <- user.lower$beta_prec)) beta_prec <- 1e-6
      lower$trbeta_prec <- .beta_precFn(beta_prec)
      if (is.null(beta_prec <- user.upper$beta_prec)) beta_prec <- max(100*canon.init$beta_prec,1e6)
      upper$trbeta_prec <- .beta_precFn(beta_prec)
    }
  } 
  if (length(fd <- canon.init$rdisPars)) { 
    if (is.list(fd)) { # mv case with >1 submodels
      for (char_mv_it in names(fd)) {
        lower$rdisPars[[char_mv_it]] <-  .modify_list(famdisp_lowup[[char_mv_it]]$lower, user.lower$rdisPars[[char_mv_it]])
        upper$rdisPars[[char_mv_it]] <-  .modify_list(famdisp_lowup[[char_mv_it]]$upper, user.upper$rdisPars[[char_mv_it]])
      }
    } else { # single vector...
      lower$rdisPars <-  .modify_list(famdisp_lowup$lower, user.lower$rdisPars)
      upper$rdisPars <-  .modify_list(famdisp_lowup$upper, user.upper$rdisPars)
    } 
  } 
  if ( ! is.null( ranCoefs <- canon.init$ranCoefs)) { ## whenever there are ranCoefs to outer-optimize 
    upper$trRanCoefs <- lower$trRanCoefs <- ranCoefs # so that assignments such as lower$trRanCoefs[[it]] <- ... will not fail
    for (char_rd in names(ranCoefs)) {
      constraint <- ranFix$ranCoefs[[char_rd]]
      init_trRancoef <- .constr_ranCoefsFn(ranCoefs[[char_rd]], constraint=constraint, rC_transf=rC_transf)
      if ( ! is.null(constraint) && attr(constraint,"isDiagFamily")) {
        if (! is.null(rC <- canon.init$ranCoefs[[char_rd]])) { # then same algo as for lambda; rC is vector
          ranCoef <- na.omit(user.lower$ranCoefs[[char_rd]]) # na.omit() means that users _may_ include NA in bounds for partially-fixed ranCoefs. 
          if (is.null(ranCoef)) ranCoef <- pmax(1e-6,rC/1e5)
          # names(ranCoef) <- names(rC) # hmm no names...
          lower$trRanCoefs[[char_rd]] <- .dispFn(ranCoef)
          ranCoef <- na.omit(user.upper$ranCoefs[[char_rd]])
          if (is.null(ranCoef)) ranCoef <- pmin(1e8,rC*1e7)
          # names(ranCoef) <- names(rC) 
          upper$trRanCoefs[[char_rd]] <- .dispFn(ranCoef)
        }
      } else {
        trRancoef_LowUp <- .calc_LowUp_trRancoef(init_trRancoef,Xi_ncol=attr(init_trRancoef,"Xi_ncol"),
                                                 tol_ranCoefs=.spaMM.data$options$tol_ranCoefs_outer,
                                                 rC_transf=rC_transf)
        lower$trRanCoefs[[char_rd]] <- trRancoef_LowUp$lower
        upper$trRanCoefs[[char_rd]] <- trRancoef_LowUp$upper
      }
    }
  }
  if ( ! is.null( beta <- canon.init$beta)) { # outer beta
    if (.spaMM.data$options$tr_beta) {
      lower$trBeta[names(beta)] <- -Inf
      upper$trBeta[names(beta)] <- Inf
    } else { # using template presumably defined by the explicit init that the user gave to select outer beta estimation
      lower$beta[names(beta)] <- -Inf
      upper$beta[names(beta)] <- Inf
      lower$beta[names(fixef_lowup$lower)] <- fixef_lowup$lower # (all wrt scaled X)
      upper$beta[names(fixef_lowup$upper)] <- fixef_lowup$upper
    }
  }
  ## names() to make sure the order of elements match; remove any extra stuff (which?... hmmm erroneous inclusion of some pars...) 
  return(list(lower=lower[names(init.optim)],upper=upper[names(init.optim)])) 
}
