.makeLowerUpper <- function(canon.init, ## cf calls: ~ in user scale, must be a full list of relevant params
                                init.optim, ## ~in transformed scale : is has all pars to be optimized
                                user.lower=list(),user.upper=list(),
                                corr_types=NULL, ranFix=list(),
                                optim.scale, moreargs=NULL) {
  lower <- upper <- init.optim   ## init.optim not further used...
  for (it in seq_along(corr_types)) {
    corr_type <- corr_types[it]
    if (! is.na(corr_type)) {
      char_rd <- as.character(it)
      if (corr_type %in% c("SAR_WWt","adjacency") && 
          ! is.null(.get_cP_stuff(init.optim,"rho",char_rd)) ## to exclude inner optimization (was not previously necessary bc 
          # sequence was calc.inits / add rhorange conditionnally on inits$init.optim / compute LowUp)
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
        hyper_map <- attr(hyper,"map")
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
        if (! is.null(canon_rho <- .get_cP_stuff(canon.init,"rho",char_rd))) {
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
      } ## end else if Matern/Cauchy case
    }
  }
  
  if (! is.null(canon.init$phi)) {
    phi <- user.lower$phi
    if (is.null(phi)) phi <- max(min(1e-6,canon.init$phi/1.01),canon.init$phi/1e5) # >=1e-6 if canon.init$phi>1e-6
    lower$trPhi <- .dispFn(phi)
    phi <- user.upper$phi
    if (is.null(phi)) phi <-  min(max(1e8,canon.init$phi*1.01),canon.init$phi*1e7) 
    ## if phi is badly initialized then it gets a default which may cause hard to catch problems in the bootstrap...
    upper$trPhi <- .dispFn(phi)
  }
  if (! is.null(canon.init$lambda)) {
    lambda <- user.lower$lambda
    if (is.null(lambda)) lambda <- pmax(1e-6,canon.init$lambda/1e5)
    lower$trLambda <- .dispFn(lambda)
    lambda <- user.upper$lambda
    if (is.null(lambda)) lambda <- pmin(1e8,canon.init$lambda*1e7)
    upper$trLambda <- .dispFn(lambda)
  }
  if (! is.null(canon.init$COMP_nu)) {
    COMP_nu <- user.lower$COMP_nu
    if (is.null(COMP_nu)) COMP_nu <- min(canon.init$COMP_nu/2,0.05)
    lower$COMP_nu <- COMP_nu
    COMP_nu <- user.upper$COMP_nu
    if (is.null(COMP_nu)) COMP_nu <- max(canon.init$COMP_nu*10,10)
    upper$COMP_nu <- COMP_nu
  } else if (! is.null(canon.init$NB_shape)) { ## for Gamma(1:sh,sh) with mean 1 and variance sh
    NB_shape <- user.lower$NB_shape
    if (is.null(NB_shape)) NB_shape <- 1e-6 
    lower$trNB_shape <- .NB_shapeFn(NB_shape)
    NB_shape <- user.upper$NB_shape
    if (is.null(NB_shape)) NB_shape <- max(100*canon.init$NB_shape,1e6)
    upper$trNB_shape <- .NB_shapeFn(NB_shape)
  }
  if ( ! is.null( ranCoefs <- canon.init$ranCoefs)) { ## whenever there are ranCoefs to outer-optimize (FIXME? no user control on canon.init?)
    upper$trRanCoefs <- lower$trRanCoefs <- ranCoefs
    for (it in seq_along(ranCoefs)) {
      init_trRancoef <- .ranCoefsFn(ranCoefs[[it]])
      trRancoef_LowUp <- .calc_LowUp_trRancoef(init_trRancoef,Xi_ncol=attr(init_trRancoef,"Xi_ncol"),
                                               tol_ranCoefs=.spaMM.data$options$tol_rel_ranCoefs)
      lower$trRanCoefs[[it]] <- trRancoef_LowUp$lower
      upper$trRanCoefs[[it]] <- trRancoef_LowUp$upper
    }
  }
  ## names() to make sure the order of elements match; remove any extra stuff (which?... hmmm erroneous inclusion of some pars...) 
  return(list(lower=lower[names(init.optim)],upper=upper[names(init.optim)])) 
}
