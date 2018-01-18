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
    if (x<=xm) { 
      num <- log(xref+x) # -log(xref+1e-6)+1
    } else {
      num <- log(xref+xm)+log((xreff+x)/(xreff+xm))*(xreff+xm)/(xref+xm) # -log(xref+1e-6)+1
    }
    #num/160 ## scaling/shifting have no effect on nloptr path nor clear effect on timings
    num
  } 
  .dispInv <- function(x,xref=5e-5,xreff=1e-1,xm=1e-3) {
    # x <- x*160 +log(xref+1e-6)-1
    if (x<log(xref+xm)) {
      exp(x)-xref
    } else {
      fac <- exp((x-log(xref+xm))*(xm+xref)/(xm+xreff))
      xm*fac + xreff*(fac-1)
    }
  }
}
#sapply(sapply(10^(seq(-6,6)),dispFn),dispInv)
## These must be directly applicable to say lambda=c(0.1,0.1) hence need to vectorize the test on x
.dispFn <- Vectorize(.dispFn,"x")
.dispInv <- Vectorize(.dispInv,"x")
## for rho/trRho: rho=0 exactly is meaningful in adjacency model. rhoFn and rhoInv should be used only for other models.
.rhoFn <- function(x,RHOMAX) {log(x/(RHOMAX-x))} ## rho should be constrained to <RHOMAX and trRho should diverge as rho approaches RHOMAX
.rhoInv <- function(trRho,RHOMAX) {
  RHOMAX*exp(trRho)/(1+exp(trRho))
} 
.nuFn <- function(nu,rho,NUMAX) {log(nu/(NUMAX-nu))} ## nu should be constrained to <NUMAX and trNu should diverge as nu approaches NUMAX
.nuInv <- function(trNu,trRho,NUMAX) {NUMAX*exp(trNu)/(1+exp(trNu))}
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

.makeLowerUpper <- function(canon.init, ## cf calls: ~ in user scale, must be a full list of relevant params
                           init.optim, ## ~in transformed scale
                           user.lower=list(),user.upper=list(),
                           corr_types=NULL,nbUnique,ranFix=list(),
                           control.dist=list(),
                           optim.scale, 
                           rhorange=NULL,
                           RHOMAX,NUMAX) {
  lower <- upper <- init.optim   ## init.optim not further used...
  for (it in seq_along(corr_types)) {
    corr_type <- corr_types[it]
    if (! is.na(corr_type)) {
      if (corr_type %in% c("SAR_WWt","adjacency")) { ## adjacency model
        eps <- (rhorange[2L]-rhorange[1L])/(2e6)  
        lower$rho <- user.lower$rho ## no transfo for adjacency model
        if (is.null(lower$rho)) lower$rho <- rhorange[1L]+eps ## may remain NULL  
        upper$rho <- user.upper$rho ## no transfo again
        if (is.null(upper$rho)) upper$rho <- rhorange[2L]-eps
      } else if (corr_type =="AR1") {
        if ( ! is.null(canon.init$ARphi)) {
          ARphi <- user.lower$ARphi
          if (is.null(ARphi)) ARphi <- -1 + 1e-6
          lower$ARphi <- ARphi
          ARphi <- user.upper$ARphi
          if (is.null(ARphi)) ARphi <- 1 - 1e-6
          upper$ARphi <- ARphi
        }
      } else if (corr_type =="Matern") { 
        if (! is.null(canon.init$rho)) {
          rho <- user.lower$rho
          if (is.null(rho)) rho <- canon.init$rho/150
          if (optim.scale=="transformed") {
            lower$trRho <- .rhoFn(rho,RHOMAX=RHOMAX)
          } else lower$rho <- rho
          rho <- user.upper$rho
          if (is.null(rho)) {
            if (inherits(nbUnique,"list")) nbUnique <- mean(unlist(nbUnique))
            rho <- canon.init$rho*2*nbUnique ## The following was a bit too low for experiments with nu=0.5 : 1/(maxrange/(2*nbUnique)) ## nb => unique rows !
            ## *modify* upper rho so that it does not exceed RHOMAX => /($RHOMAX+...)
            if (optim.scale=="transformed") rho <- 2*rho * RHOMAX/(RHOMAX+2*rho)
          }
          if (optim.scale=="transformed") {
            upper$trRho <- .rhoFn(rho,RHOMAX=RHOMAX) 
          } else upper$rho <- rho
          rhoForNu <- canon.init$rho
        } else rhoForNu <- .getPar(ranFix,"rho")
        if (! is.null(canon.init$nu)) {
          nu <- user.lower$nu
          if (is.null(nu)) nu <- canon.init$nu/100
          if (optim.scale=="transformed") {
            lower$trNu <- .nuFn(nu,rhoForNu,NUMAX)
            #print(c(rhoForNu,nu,lower$trNu))
          } else lower$nu <-nu
          nu <- user.upper$nu
          if (is.null(nu)) {
            if ( ! is.null(dm <- control.dist$`dist.method`) && dm %in% c("Geodesic","Earth")) {
              nu <- 0.5
            } else {
              ## constructs upper nu from NUMAX => /(1+...)
              ## nu should not diverge otherwise it will diverge in Bessel_lnKnu, whatever the transformation used
              nu <- NUMAX * canon.init$nu/(1+canon.init$nu) 
              ## FR->FR hmmm. If canon.init$nu= NUMAX-1 then 
              ##  (upper) nu= canon.init$nu and possibly < canon.init$nu by numerical accuracy issues => nloptr stops
            }
          }
          if (optim.scale=="transformed") {
            upper$trNu <- .nuFn(nu,rhoForNu,NUMAX)
          } else upper$nu <- nu
          #print(c(rhoForNu,nu,upper$trNu))
        }
        if ( ! is.null(canon.init$Nugget)) {
          lower$Nugget <- 0
          upper$Nugget <- 0.999999
        }
      } ## end else if Matern case
    }
  }

  if (! is.null(canon.init$phi)) {
    phi <- user.lower$phi
    if (is.null(phi)) phi <- max(1e-6,canon.init$phi/1e5)
    lower$trPhi <- .dispFn(phi)
    phi <- user.upper$phi
    if (is.null(phi)) phi <-  min(1e8,canon.init$phi*1e7)
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
  if ( ! is.null( ranCoefs <- canon.init$ranCoefs)) { ## not tested
    canon_LowUp_ranCoef <- lapply(ranCoefs, .calc_canon_LowUp_ranCoef) ## canonical scale: compare with init before log transfo
    canon_lower <- lapply(canon_LowUp_ranCoef, getElement,"lower")
    lower$trRanCoefs <- lapply(canon_lower, .ranCoefsFn)
    canon_upper <- lapply(canon_LowUp_ranCoef, getElement,"upper")
    upper$trRanCoefs <- lapply(canon_upper, .ranCoefsFn)
  }
  ## names() to make sure the order of elements match; remove any extra stuff (which?)
  return(list(lower=lower[names(init.optim)],upper=upper[names(init.optim)])) 
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
  fixedlambda <- .getPar(ranFix,"lambda") ## FR->FR his comes from user arg, not from preprocessing
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
    init$lambda[init$lambda<1e-4] <- 1e-4
    init.optim$trLambda <- .dispFn(init$lambda) 
  }
  init.optim$lambda <- NULL
  if (is.null(.getPar(ranFix,"phi"))) {
    init$phi <- init.optim$phi 
    if (!is.null(init$phi)) {
      init$phi[init$phi<1e-4] <- 1e-4
      init.optim$trPhi <- .dispFn(init$phi)
      init.optim$phi <- NULL
    }
  } else {
    if (!is.null(init.optim$phi)) stop("(!) Arguments 'ranFix$phi' and 'init$phi' conflict with each other.")  
  } 
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

.calc_inits_ARphi <- function(init,init.optim,init.HLfit,ranFix,user.lower=list(),user.upper=list()) {
  if (is.null(.getPar(ranFix,"ARphi")) && (! is.numeric(init.HLfit$ARphi))) { 
    init$ARphi <- init.optim$ARphi 
    if (is.null(init$ARphi)) init$ARphi <- (user.lower$ARphi+user.upper$ARphi)/2 ## only a quick debugging aide, but FIXME:
    # .... set this working more generally in .calc_inits_.... ? 
    if ( ! length(init$ARphi)==1L) init$ARphi <- 0. 
    if (! is.null(init.HLfit$ARphi)) {
      init.HLfit$ARphi <- init$ARphi 
    } else {
      init.optim$ARphi <- init$ARphi
    }
  } 
  if ( (! is.null(init.HLfit$ARphi)) && (! is.numeric(init.HLfit$ARphi))) {
    init.HLfit$ARphi <- init.optim$ARphi
    init.optim$ARphi <- NULL
  }      
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

.calc_inits_nu <- function(init,init.optim,init.HLfit,ranFix,control.dist,optim.scale,NUMAX,user.lower,user.upper) {
  if (is.null(.getPar(ranFix,"nu")) && (! is.numeric(init.HLfit$nu))) { 
    init$nu <- init.optim$nu 
    if (is.null(init$nu)) {
      if ( ! is.null(dm <- control.dist$`dist.method`) && dm %in% c("Geodesic","Earth")) {
        init$nu <- 0.25  
      } else init$nu <- 0.5 
    }
    if (! is.null(init.HLfit$nu)) {
      init.HLfit$nu <- init$nu ## avant transformation
    } else {
      if (optim.scale=="transformed") {
        Fixrho <- .getPar(ranFix,"rho")
        if (is.null(Fixrho)) { 
          init.optim$trNu <- .nuFn(init$nu,init$rho,NUMAX=NUMAX) 
        } else init.optim$trNu <- .nuFn(init$nu,Fixrho,NUMAX=NUMAX)
        init.optim$nu <- NULL
      } else init.optim$nu <- init$nu
    }
  } 
  if ( (! is.null(init.HLfit$nu)) && (! is.numeric(init.HLfit$nu))) {
    init.HLfit$nu <- init.optim$nu
    init.optim$nu <- NULL
  }
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

.calc_inits_Auto_rho <- function(init,init.optim,init.HLfit,ranFix,rhorange,For,user.lower,user.upper) {
  if (is.null(.getPar(ranFix,"rho"))) { ## $rho NULL or NA
    if (For=="corrHLfit") { ## with corrHLfit, defaul is outer optim. We provide init.optim
      if (is.null(init.HLfit$rho)) { ## default init.HLfit
        if (! is.numeric(init.optim$rho)) init.optim$rho <- mean(rhorange)
      } else if ( ! is.numeric(init.HLfit$rho) ) { ## non-default: inner optim, but no numeric init  
        init.HLfit$rho <- mean(rhorange)
        init.optim$rho <- NULL 
      } else { ## non-default: inner optim, already with init value
        init.optim$rho <- NULL 
      }
    } else if (For=="fitme") { ## with fitme, default is inner optim. We provide init.HLfit
      if ( ! is.null(init.HLfit$rho)) { ## non-default, forces inner estim
        if (! is.numeric(init.HLfit$rho)) init.HLfit$rho <- mean(rhorange)
      } else if ( ! is.numeric(init.optim$rho) ) { ## default: outer optim, but no numeric init
        init.optim$rho <- mean(rhorange)
        init.HLfit$rho <- NULL
      } else { ## non-default: inner optim, already with init value
        init.HLfit$rho <- NULL
      }
      # if (is.null(init.optim$rho)) { ## default : inner estim
      #   if (! is.numeric(init.HLfit$rho)) init.HLfit$rho <- mean(rhorange)
      # } else if ( ! is.numeric(init.optim$rho) ) { ## non-default: outer optim, but no numeric init  
      #   init.optim$rho <- mean(rhorange) 
      #   init.HLfit$rho <- NULL
      # } else { ## non-default: inner optim, already with init value
      #   init.HLfit$rho <- NULL 
      # }
    }
  }
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}


.calc_inits_geostat_rho <- function(init,init.optim,init.HLfit,ranFix,maxrange,optim.scale,RHOMAX,user.lower,user.upper) {
  if (is.null(.getPar(ranFix,"rho")) && (! is.numeric(init.HLfit$rho))) {
    init$rho <- init.optim$rho 
    if (is.null(init$rho)) {
      init$rho <- 30/(2*maxrange)
    } else if (any( narho <- is.na(init$rho))) init$rho[narho] <- 30/(2*maxrange[narho]) ## 05/2015 allows init.corrHLfit=list(rho=rep(NA,...
    if (! is.null(init.HLfit$rho)) { ## ! is.null, but ! is.numeric as tested above; but now it becomes numeric
      init.HLfit$rho <- init$rho ## avant transformation
    } else {
      if (optim.scale=="transformed") {
        init.optim$trRho <- .rhoFn(init$rho,RHOMAX=RHOMAX) ## we're in Matern model here
        init.optim$rho <- NULL
      } else init.optim$rho <- init$rho
    }
  } 
  if ( (! is.null(init.HLfit$rho)) && (! is.numeric(init.HLfit$rho))) {
    init.HLfit$rho <- init.optim$rho
    init.optim$rho <- NULL
  }
  return(list(init=init,init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix))
}

.calc_inits <- function(init.optim,init.HLfit,ranFix,corr_types,rhorange=NULL,maxrange=NULL,
                        user.lower,user.upper,
                       optim.scale,control.dist,RHOMAX,NUMAX,For) { 
  inits <- as.list(match.call()[c("init.optim","init.HLfit","ranFix","user.lower","user.upper")])
  inits <- c(inits,list(init=list()))
  for (it in seq_along(corr_types)) {
    corr_type <- corr_types[it]
    if (! is.na(corr_type)) {
      if (corr_type =="Matern") {
        arglist <- c(inits,list(maxrange=maxrange,optim.scale=optim.scale,RHOMAX=RHOMAX))
        inits <- do.call(.calc_inits_geostat_rho,arglist)
        arglist <- c(inits,list(control.dist=control.dist, optim.scale=optim.scale, NUMAX=NUMAX))
        inits <- do.call(.calc_inits_nu,arglist)
        # Nugget: remains NULL through all computations if init.optim$Nugget is NULL
        if (is.null(.getPar(ranFix,"Nugget"))) { inits$init["Nugget"] <- init.optim$Nugget }  
      } else if ( corr_type  %in% c("SAR_WWt","adjacency") ) { 
        arglist <- c(inits,list(rhorange=rhorange,For=For))
        inits <- do.call(.calc_inits_Auto_rho,arglist)
      } else if (corr_type =="AR1") {
        inits <- do.call(.calc_inits_ARphi,inits)
      }  
    }
  }
  # phi, lambda
  inits <- do.call(.calc_inits_dispPars,inits)
  # random coefficients
  inits <- do.call(.calc_inits_ranCoefs,inits)
  # GLM family parameters
  inits$init$COMP_nu <- inits$init.optim$COMP_nu ## may be NULL. No checks needed
  inits$init$NB_shape <- inits$init.optim$NB_shape ## may be NULL. No checks needed
  if (! is.null(inits$init.optim$NB_shape)) {
    inits$init.optim$trNB_shape <- .NB_shapeFn(inits$init.optim$NB_shape)
    inits$init.optim$NB_shape <- NULL
  }
  return(inits)
}


