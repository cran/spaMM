if (FALSE) {
  .dispFn <- function(x) log(0.9+x)
  .dispInv <- function(x) exp(x)-0.9
} else {
  ## see control_nloptr.nb for alternatives
  ## generic .dispInv for devel
  # .dispInv <- function(x,xref=1e-4,xreff=1e-2,xm=1e-3) {
  #   uniroot(function(v) dispFn(v)-x,lower=1e-6,upper=1e6,tol=1e-7)$root
  # }
  ## complex but still useful in
  .dispFn <- function(x,xref=5e-5,xreff=1e-1,xm=1e-3) { 
    num <- x
    x_lo_xm <- (x<=xm)
    num[which(x_lo_xm)] <- log(xref+x[which(x_lo_xm)])
    num[which( ! x_lo_xm)] <- log(xref+xm)+log((xreff+x[which( ! x_lo_xm)])/(xreff+xm))*(xreff+xm)/(xref+xm)
    #print(num)
    num
  } 
  .dispInv <- function(x,xref=5e-5,xreff=1e-1,xm=1e-3) {
    num <- x
    x_lo_xm <- (x<log(xref+xm))
    num[which(x_lo_xm)] <- exp(x[which(x_lo_xm)])-xref
    fac <- exp((x[which( ! x_lo_xm)]-log(xref+xm))*(xm+xref)/(xm+xreff))
    num[which( ! x_lo_xm)] <-  xm*fac + xreff*(fac-1)
    #print(num)
    num
  }
}
#sapply(sapply(10^(seq(-6,6)),dispFn),dispInv)
## These must be directly applicable to say lambda=c(0.1,0.1) hence need to vectorize the test on x
# but Vectorize is not quite efficient
#.dispFn <- Vectorize(.dispFn,"x")
#.dispInv <- Vectorize(.dispInv,"x")
## for rho/trRho: rho=0 exactly is meaningful in adjacency model. rhoFn and rhoInv should be used only for other models.
.rhoFn <- function(x,RHOMAX) {log(x/(RHOMAX-x))} ## rho should be constrained to <RHOMAX and trRho should diverge as rho approaches RHOMAX
.rhoInv <- function(trRho,RHOMAX) {
  RHOMAX*exp(trRho)/(1+exp(trRho))
} 
.nuFn <- function(nu,rho,NUMAX) {log(nu/(NUMAX-nu))} ## nu should be constrained to <NUMAX and trNu should diverge as nu approaches NUMAX
.nuInv <- function(trNu,trRho,NUMAX) {NUMAX*exp(trNu)/(1+exp(trNu))}
.longdepFn <- function(longdep,LDMAX) {log(longdep/(LDMAX-longdep))} 
.longdepInv <- function(trLongdep,LDMAX) {LDMAX*exp(trLongdep)/(1+exp(trLongdep))}
.ranCoefsFn <- function(vec) {
  if ( ! is.null(vec)) {
    Xi_cols <- floor(sqrt(2*length(vec)))
    lampos <- cumsum(seq(Xi_cols))
    vec[lampos] <- .dispFn(vec[lampos])
  }
  return(vec)
}
.ranCoefsInv <- function(vec) {
  if ( ! is.null(vec)) {
    Xi_cols <- floor(sqrt(2*length(vec)))
    lampos <- cumsum(seq(Xi_cols))
    vec[lampos] <- .dispInv(vec[lampos])
  }
  return(vec)
}

.NB_shapeFn <- function(x) {log(log(1+x))} ## drastic handling of flat likelihoods for high shape= LOW variance. .../...
# .../... negbin example in gentle intro is a test (using optimize() -> initial value cannot be controlled)
.NB_shapeInv <- function(x) {exp(exp(x))-1}


## FR->FR rho/sqrt(nu) scaling => should be fixed nu and transformed rho to handle vectorial rho
## thus nu fns should be indep of nu and rho fns should be functions of nu ! :
if (FALSE) {
  futurerhoFn <- function(rho,nu) {
    rhosc <- rho/sqrt(nu)
    log(rhosc/(.spaMM.data$options$RHOMAX-rhosc))
  }
  futurerhoInv <- function(trRhosc,trNu) {
    nu <- .spaMM.data$options$NUMAX*exp(trNu)/(1+exp(trNu))
    rhosc <- .spaMM.data$options$RHOMAX*exp(trRhosc)/(1+exp(trRhosc))
    rhosc*sqrt(nu)
  }
}

.calc_canon_LowUp_ranCoef <- function(ranCoefTerm,tol_ranCoefs=1e-08) { ## FIXME not yet user control
  lower <- upper <- numeric(length(ranCoefTerm))
  Xi_ncol <- floor(sqrt(2*length(ranCoefTerm)))
  varpos <- cumsum(seq(Xi_ncol))
  lower[varpos] <- tol_ranCoefs
  upper[varpos] <- 1/tol_ranCoefs
  lower[-varpos] <-   -(1-tol_ranCoefs)
  upper[-varpos] <-   (1-tol_ranCoefs)
  lower <- pmin(lower, ranCoefTerm)
  upper <- pmax(upper, ranCoefTerm)
  return(list(lower=lower,upper=upper)) ## CANONICAL scale
}

.match_coords_in_tUniqueGeo <- function(coords,tUniqueGeo) {any(apply(tUniqueGeo,1L,identical,y=coords))}

.checkDistMatrix <- function(distMatrix,data,coordinates) {
  ## HLCor + non-null distMatrix (cf example ?HLCor MLdistMat) => plantage un peu bizarre dans .checkDistMatrix, 
  #   problème de transposition de matrice uniqueGeo. A priori le code actuel est correct, mais logique pas claire.
  if (inherits(distMatrix,"dist")) {
    usernames <- labels(distMatrix)
  } else if (inherits(distMatrix,"matrix")) {
    usernames <- rownames(distMatrix)
  } else message(paste("(!) 'distMatrix' is neither a 'matrix' or 'dist' object. Check the input. I exit."))
  ## chol() fails on distances matrices with repeated locations (which are pos SD)... but chol() not used by default
  ## the following code assumes that distMatrix deals only with unique locations, and checks this
  ## HENCE ******* distMatrix must refer to unique values of a grouping variable *********
  checknames <- all(sapply(usernames,`%in%`, table=rownames(data))) ## 
  if (!checknames) {
    warning("The rownames of 'distMatrix' are not rownames of the 'data'. Further checking of 'distMatrix' is not possible.")
    nbUnique <- NA
  } else {
    uniqueGeo <- .calcUniqueGeo(data=data[usernames,coordinates,drop=FALSE]) ## check that this corresponds to unique locations
    nbUnique <- nrow(uniqueGeo)
    if (nbUnique != nrow(distMatrix)) {
      stop("The dimension of 'distMatrix' does not match the number of levels of the grouping variable")
    } else { ## check order
      redondGeo <- data[,coordinates,drop=F]
      tUniqueGeo <- t(uniqueGeo)
      designRU <- apply(redondGeo, 1L, .match_coords_in_tUniqueGeo, tUniqueGeo=tUniqueGeo) ## has no names
      # designRU <- apply(redondGeo,1, function(v) {which(apply(v==t(uniqueGeo),2,all))}) ## has no names
      ## eg 1 1 2 2 3 2 3 4 is valid for 8 obs, 4 unique locations
      designRU <- unique(as.vector(designRU)) ## should then be 1 2 3 4
      ## but if distMatrix in reverse order, the first row of redondGeo would match the 4th of uniqueGeo and then the following test is FALSE:
      if ( ! all (designRU==seq_len(length(designRU))) ) {
        stop("The rows of 'distMatrix' are not ordered as rows of the 'data'.")
      }
    } 
  }
  nbUnique ## if stop() did not occur
}

.expand_GeoMatrices <- function(w_uniqueGeo, e_uniqueGeo, coords_nesting, coord_within, dist.method) {
  rownames(e_uniqueGeo) <- seq(nrow(e_uniqueGeo)) ## local & unique rownames
  ## Remarkably the next line works only if the cols are not factors ! Unless we have a fix for this,
  #  uniqueGeo classes should be integer not factor: see instances of as.numeric(levels(fac)) in the sources.
  rows_bynesting <- by(e_uniqueGeo ,e_uniqueGeo[,coords_nesting],rownames) 
  distMatrix <- matrix(Inf,ncol=nrow(e_uniqueGeo),nrow =nrow(e_uniqueGeo))
  rownames(distMatrix) <- colnames(distMatrix) <- rownames(e_uniqueGeo) ## same trivial rownames
  for (lit in seq_len(length(rows_bynesting))) {
    blockrows <- rows_bynesting[[lit]]
    within_values <- e_uniqueGeo[blockrows,coord_within]
    w_dist <- proxy::dist(within_values,method=dist.method)
    distMatrix[blockrows,blockrows] <- as.matrix(w_dist)
  }
  return(list(distMatrix=distMatrix))
}

.calc_fullrho <- function(rho, coordinates,rho_mapping) {
  if ( is.null(rho_mapping) ) {
    if (length(rho)==1L ) rho <- rep(rho,length(coordinates))
    names(rho) <- coordinates
  } else {
    rho <- rho[rho_mapping]
    names(rho) <- names(rho_mapping)
  }  
  rho
}


.reformat_verbose <- function(verbose,For) {
  if (is.null(verbose)) verbose <- logical(0)
  if (is.list(verbose)) stop("is.list(verbose)")
  if (is.na(verbose["TRACE"])) verbose["TRACE"] <- FALSE ## for fitme
  if (is.na(verbose["trace"])) verbose["trace"] <- FALSE
  if (is.na(verbose["SEM"])) verbose["SEM"] <- FALSE
  if (is.na(verbose["iterateSEM"])) verbose["iterateSEM"] <- TRUE ## summary info and plots for each iteration
  ## when a fit is by HLCor or HLfit, the serious message at the end of HLfit_body are by default displayed,
  ## but not in the other cases:
  if (is.na(verbose["all_objfn_calls"])) verbose["all_objfn_calls"] <- switch(For, "HLCor" = TRUE, "HLfit" = TRUE, 
                                                        "corrHLfit" = FALSE, "fitme" = FALSE)
  return(verbose)
}

.reformat_LevenbergM <- function(LevenbergM) {
  if (is.list(LevenbergM)) stop("is.list(LevenbergM)")
  ## allows user to control the default globally by spaMM.options(): BUT this should typically be NULL
  default_LM_start <- .spaMM.data$options$LevenbergM
  if ( ! is.logical(default_LM_start)) default_LM_start <- NA ## converts NULL to NA (default resolved further in preprocess)
  default_prefit <- .spaMM.data$options$PQL_prefit
  if ( ! is.logical(default_prefit)) default_prefit <- NA ## converts NULL to NA (default resolved further in preprocess)
  resu <- c(user_LM=LevenbergM["LevenbergM"][[1L]], user_prefit=LevenbergM["PQL_prefit"][[1L]], 
            # [.][[1L]] important to get rid of name yet be valid for implicit NAs
            default_LM_start=default_LM_start, default_prefit=default_prefit, force=NA) 
  if (length(LevenbergM)==1L && is.logical(LevenbergM)) resu["user_LM"] <- LevenbergM ## T, F, or NA : handles non-vector input
  return(resu) 
} 

.determine_LevenbergM <- function(LevenbergM, which="user_LM") {
  if (is.na(info <- LevenbergM["force"])) {
    if (is.na(info <- LevenbergM[which])) {
      if (which=="user_prefit") {
        info <- LevenbergM["default_prefit"]
      } else info <- LevenbergM["default_LM_start"]
    }  
  }
  if (! is.logical(info) || is.na(info)) stop("invalid input: check input (invalid 'which' or invalid user input?)")
  return(info)
} 

.calc_inits_ranCoefs <- function(init,init.optim,init.HLfit,ranFix,user.lower,user.upper) {
  # fixme c/should allow a user-provided init.HLfit$ranCoefs to control what remains in init.optim relative to init.HLfit 
  if ( ! is.null(init.optim$ranCoefs)) { ## should always be a complete ordered list in the case for random-coefficient models
    if (! is.null(init$ranCoefs)) for (st in names(init$ranCoefs)) init.optim$ranCoefs[[st]] <- init$ranCoefs[[st]]
    # fixme ideally the code whould check that st corresponds to a random-coefficient term so that the user does not mess everything....
    init$ranCoefs <- init.optim$ranCoefs
    trRanCoefs <- lapply(init.optim$ranCoefs,.ranCoefsFn)
    init.optim$trRanCoefs <- trRanCoefs
    init.optim$ranCoefs <- NULL
  }
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

.calc_inits_dispPars <- function(init,init.optim,init.HLfit,ranFix,user.lower,user.upper) {
  ## does not modify init.HLfit, but keeps its original value. Also useful to keep ranFix for simple and safe coding
  init$lambda <- init.optim$lambda 
  fixedlambda <- .getPar(ranFix,"lambda") ## FIXME getPar not useful ?
  if (!is.null(fixedlambda)) {
    whichNA <- which(is.na(fixedlambda))
    if (length(init$lambda)==length(whichNA)) {
      if (length(whichNA)) names(init$lambda) <- whichNA
    } else if (length(init$lambda)<length(whichNA)) {
      ## init$lambda should already have names...
    } else stop("'fixed lambda' and 'init lambda' arguments conflict with each other.")
  }
  if (length(init$lambda)) { ## length rather than !is.null() bc .dispFn(numeric(0)) returns list() bc of the Vectorize code:
    ## zut <- function(x) return(666); zut <- Vectorize(zut); zut(numeric(0)) ## gives list()
    init$lambda <- pmax(init$lambda, 1e-4) ## see remark on init$phi below. Lower user- or bootstrap- inits are not heeded.
    init.optim$trLambda <- .dispFn(init$lambda) 
  }
  init.optim$lambda <- NULL
  if (is.null(.getPar(ranFix,"phi"))) { ## FIXME getPar not useful ?
    init$phi <- init.optim$phi 
    if (!is.null(init$phi)) {
      ## pmax(), really ? 
      ## Only for generic optimization
      ##    which then means that lower user- or bootstrap- inits are not heeded. But test-fixedLRT supports this.
      init$phi <- pmax(init$phi, 1e-4) 
      init.optim$trPhi <- .dispFn(init$phi)
      init.optim$phi <- NULL
    }
  } else {
    if (!is.null(init.optim$phi)) stop("(!) Arguments 'ranFix$phi' and 'init$phi' conflict with each other.")  
  } 
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

.calc_inits_ARphi <- function(init,init.optim,init.HLfit,ranFix,user.lower=list(),user.upper=list(), char_rd) {
  if (is.null(ranFix$corrPars[[char_rd]][["ARphi"]])) {
    if (! is.numeric(HLfit_ARphi <- init.HLfit$corrPars[[char_rd]][["ARphi"]])) { ## if condition is TRUE then HLfit_ARphi should be NA or NULL
      init_ARphi <- init.optim$corrPars[[char_rd]][["ARphi"]] 
      if (is.null(init_ARphi)) init_ARphi <- (.get_cP_stuff(user.lower,"ARphi",char_rd)+.get_cP_stuff(user.upper,"ARphi",char_rd))/2 
      if ( ! length(init_ARphi)) init_ARphi <- 0. ## may be numeric(0) from previous line
      if (! is.null(HLfit_ARphi)) { ## then typically NA, given that the original init.HLfit$ was not numeric... 
        init.HLfit$corrPars[[char_rd]] <- list(ARphi=init_ARphi) ## numeric
        init.optim$corrPars[[char_rd]] <- NULL ## 'should' not be necessary
      } else { ## then, HLfit_ARphi was NULL and stays NULL
        init.optim$corrPars[[char_rd]] <- list(ARphi=init_ARphi) ## numeric
      } ## now HLfit_ARphi is NULL or numeric, and optim_ARphi is numeric 
      init$corrPars[[char_rd]] <- list(ARphi=init_ARphi)
    } else { ## else init.HLfit$ -> HLfit_ARphi is numeric
      init.optim$corrPars[[char_rd]] <- NULL ## 'should' not be necessary
      init$corrPars[[char_rd]] <- list(ARphi=init_ARphi)
    }
  }
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

.calc_inits_nu <- function(init,init.optim,init.HLfit,ranFix,control_dist_rd,optim.scale,NUMAX,user.lower,user.upper,char_rd) {
  if (is.null(.get_cP_stuff(ranFix,"nu",which=char_rd))) { 
    nu <- .get_cP_stuff(init.optim,"nu",which=char_rd)
    if (is.null(nu)) {
      if ( ! is.null(dm <- control_dist_rd$`dist.method`) && dm %in% c("Geodesic","Earth")) {
        nu <- 0.25  
      } else nu <- 0.5 
    }
    init$corrPars[[char_rd]] <- .modify_list(init$corrPars[[char_rd]], list(nu=nu)) ## canonical scale
    optim_cP <- init.optim$corrPars[[char_rd]]
    if (optim.scale=="transformed") {
      rho <- .get_cP_stuff(ranFix,"rho",which=char_rd)
      if (is.null(rho)) rho <- .get_cP_stuff(init,"rho",which=char_rd)
      optim_cP <- .modify_list(optim_cP, list(trNu=.nuFn(nu,rho,NUMAX=NUMAX))) 
      optim_cP$nu <- NULL
    } else optim_cP <- .modify_list(optim_cP, list(nu=nu)) 
    init.optim$corrPars[[char_rd]] <- optim_cP
  } 
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

.calc_inits_cauchy <- function(init,init.optim,init.HLfit,ranFix,control_dist_rd,optim.scale,LDMAX,user.lower,user.upper,char_rd) {
  if (is.null(.get_cP_stuff(ranFix,"shape",which=char_rd))) { 
    shape <- .get_cP_stuff(init.optim,"shape",which=char_rd)
    if (is.null(shape)) shape <- 1
    init$corrPars[[char_rd]] <- .modify_list(init$corrPars[[char_rd]], list(shape=shape)) ## canonical scale
    optim_cP <- init.optim$corrPars[[char_rd]]
    optim_cP <- .modify_list(optim_cP, list(shape=shape)) 
    init.optim$corrPars[[char_rd]] <- optim_cP
  }
  if (is.null(.get_cP_stuff(ranFix,"longdep",which=char_rd))) { 
    longdep <- .get_cP_stuff(init.optim,"longdep",which=char_rd)
    if (is.null(longdep)) longdep <- 0.5
    init$corrPars[[char_rd]] <- .modify_list(init$corrPars[[char_rd]], list(longdep=longdep)) ## canonical scale
    optim_cP <- init.optim$corrPars[[char_rd]]
    if (optim.scale=="transformed") {
      optim_cP <- .modify_list(optim_cP, list(trLongdep=.longdepFn(longdep,LDMAX=LDMAX))) 
      optim_cP$longdep <- NULL
    } else optim_cP <- .modify_list(optim_cP, list(longdep=longdep)) 
    init.optim$corrPars[[char_rd]] <- optim_cP
  }
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}



.calc_inits_Auto_rho <- function(init,init.optim,init.HLfit,ranFix,rhorange,For,user.lower,user.upper,char_rd) {
  if (is.null(.get_cP_stuff(ranFix,"rho",which=char_rd))) { ## $rho NULL or NA
    init_rho <- .get_cP_stuff(init.HLfit,"rho",which=char_rd)
    if (For=="corrHLfit") { ## with corrHLfit, defaul is outer optim. We provide init.optim
      if (is.null(init_rho)) { ## default init.HLfit
        init_rho <- .get_cP_stuff(init.optim,"rho",which=char_rd)
        if (! is.numeric(init_rho)) init.optim$corrPars[[char_rd]] <- init_rho <- list(rho=mean(rhorange))
      } else if ( ! is.numeric(init_rho) ) { ## non-default: inner optim, but no numeric init  
        init.HLfit$corrPars[[char_rd]] <- init_rho <- list(rho=mean(rhorange))
        init.optim$corrPars[[char_rd]] <- NULL 
      } else { ## non-default: inner optim, already with init value
        init.optim$corrPars[[char_rd]] <- NULL 
      }
    } else if (For=="fitme") { ## with fitme, default is inner optim. We provide init.HLfit
      if (is.null(init_rho)) {
        init_rho <- .get_cP_stuff(init.optim,"rho",which=char_rd)
        if ( ! is.numeric(init_rho) ) { ## default: outer optim, but no numeric init
          init.optim$corrPars[[char_rd]] <- init_rho <- list(rho=mean(rhorange))
          init.HLfit$corrPars[[char_rd]] <- NULL
        } else { ## non-default: inner optim, already with init value
          init.HLfit$corrPars[[char_rd]] <- NULL
        }
      } else if (! is.numeric(init_rho)) init.HLfit$corrPars[[char_rd]] <- init_rho <- list(rho=mean(rhorange))
    }
    init$corrPars[[char_rd]] <- init_rho ## assignment absent (still in p2 version)  prior to v2.3.30
  }
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

.calc_inits_geostat_rho <- function(init,init.optim,init.HLfit,ranFix,maxrange,optim.scale,RHOMAX,user.lower,user.upper,char_rd) {
  if (is.null(.get_cP_stuff(ranFix,"rho",which=char_rd))) {
    rho <- .get_cP_stuff(init.optim,"rho",which=char_rd)
    if (is.null(rho)) {
      rho <- 30/(2*maxrange)
    } else if (any( narho <- is.na(rho))) rho[narho] <- 30/(2*maxrange[narho]) ## 05/2015 allows init.corrHLfit=list(rho=rep(NA,...
    optim_cP <- init.optim$corrPars[[char_rd]]
    init$corrPars[[char_rd]] <- .modify_list(init$corrPars[[char_rd]], list(rho=rho))
    if (optim.scale=="transformed") {
      optim_cP <- .modify_list(optim_cP, list(trRho=.rhoFn(rho,RHOMAX=RHOMAX)))
      optim_cP$rho <- NULL
    } else optim_cP <- .modify_list(optim_cP, list(rho=rho))
    init.optim$corrPars[[char_rd]] <- optim_cP
  } 
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

.calc_inits <- function(init.optim,init.HLfit,ranFix,corr_types,
                        moreargs,
                        user.lower,user.upper,
                       optim.scale,For) { 
  inits <- list(init=list(corrPars=list()), ## minimal structure for assignments $corrPars[[char_rd]][[<name>]] to work.
                init.optim=init.optim, init.HLfit=init.HLfit, ranFix=ranFix,
                user.lower=user.lower, user.upper=user.upper, init=list()) 
  for (rd in seq_along(corr_types)) {
    corr_type <- corr_types[rd]
    if (! is.na(corr_type)) {
      char_rd <- as.character(rd)
      if (corr_type =="Matern") {
        moreargs_rd <- moreargs[[char_rd]]
        inits <- .calc_inits_geostat_rho(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                                         user.lower=user.lower,user.upper=user.upper,
                                         maxrange=moreargs_rd$maxrange,optim.scale=optim.scale,RHOMAX=moreargs_rd$RHOMAX,char_rd=char_rd)
        inits <- .calc_inits_nu(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                                control_dist_rd=moreargs_rd$control.dist, optim.scale=optim.scale, NUMAX=moreargs_rd$NUMAX,char_rd=char_rd)
        # Nugget: remains NULL through all computations if NULL in init.optim
        if (is.null(.get_cP_stuff(ranFix,"Nugget",which=char_rd))) { 
          if (is.null(init.optim$corrPars)) { ## old sp2 code 
            if (is.null(.getPar(ranFix,"Nugget"))) { inits$init["Nugget"] <- init.optim$Nugget }
          } else { ## new spaMM3.0 code
            inits$init$corrPars[[char_rd]] <- .modify_list(inits$init$corrPars[[char_rd]],
                                                           list(Nugget=.get_cP_stuff(init.optim,"Nugget",which=char_rd)))
          }  
        }
      } else if (corr_type =="Cauchy") {
        moreargs_rd <- moreargs[[char_rd]]
        inits <- .calc_inits_geostat_rho(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                                         user.lower=user.lower,user.upper=user.upper,
                                         maxrange=moreargs_rd$maxrange,optim.scale=optim.scale,RHOMAX=moreargs_rd$RHOMAX,char_rd=char_rd)
        inits <- .calc_inits_cauchy(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                                control_dist_rd=moreargs_rd$control.dist, optim.scale=optim.scale, LDMAX=moreargs_rd$LDMAX,char_rd=char_rd)
        # Nugget: remains NULL through all computations if NULL in init.optim
        if (is.null(.get_cP_stuff(ranFix,"Nugget",which=char_rd))) { 
          if (is.null(init.optim$corrPars)) { ## old sp2 code 
            if (is.null(.getPar(ranFix,"Nugget"))) { inits$init["Nugget"] <- init.optim$Nugget }
          } else { ## new spaMM3.0 code
            inits$init$corrPars[[char_rd]] <- .modify_list(inits$init$corrPars[[char_rd]],
                                                           list(Nugget=.get_cP_stuff(init.optim,"Nugget",which=char_rd)))
          }  
        }
      } else if ( corr_type  %in% c("SAR_WWt","adjacency") ) { 
        inits <- .calc_inits_Auto_rho(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                                      user.lower=user.lower,user.upper=user.upper,
                                      rhorange=moreargs[[char_rd]]$rhorange,For=For,char_rd=char_rd)
      } else if (corr_type =="AR1") {
        inits <- .calc_inits_ARphi(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                                   user.lower=user.lower,user.upper=user.upper,char_rd=char_rd)
      }  
    }
  }
  # phi, lambda
  inits <- .calc_inits_dispPars(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                                user.lower=user.lower,user.upper=user.upper)
  # random coefficients
  inits <- .calc_inits_ranCoefs(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                                user.lower=user.lower,user.upper=user.upper)
  # GLM family parameters
  inits$init$COMP_nu <- inits$init.optim$COMP_nu ## may be NULL. No checks needed
  inits$init$NB_shape <- inits$init.optim$NB_shape ## may be NULL. No checks needed
  if (! is.null(inits$init.optim$NB_shape)) {
    inits$init.optim$trNB_shape <- .NB_shapeFn(inits$init.optim$NB_shape)
    inits$init.optim$NB_shape <- NULL
  }
  return(eval(inits))
}
