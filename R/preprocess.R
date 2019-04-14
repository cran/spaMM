.checkRandLink <- function(rand.family) {
  lcrandfamfam <- tolower(rand.family$family) ## tolower once and for all
  oklink <- F
  ## cases where g(u)=th(u)
  if (lcrandfamfam=="gaussian" && rand.family$link=="identity") oklink <- TRUE          
  if (lcrandfamfam=="gamma" && rand.family$link=="log") oklink <- TRUE          
  if (lcrandfamfam=="inverse.gamma" && rand.family$link=="-1/mu") oklink <- TRUE
  if (lcrandfamfam=="beta" && rand.family$link=="logit") oklink <- TRUE
  ## cases where g(u)!=th(u)
  if (lcrandfamfam=="inverse.gamma" && rand.family$link=="log") oklink <- TRUE 
  if (lcrandfamfam=="gamma" && rand.family$link=="identity") oklink <- TRUE ## gamma(identity)
  if ( ! oklink) {
    allowed <- switch(lcrandfamfam,
                      gaussian= "is 'identity'",
                      gamma= "is 'log' (explicit 'identity' may be *tried*)", ## gamma(identity) not yet working
                      #  and the above code protext against Gamma(inverse); it could perhaps be allowed for trying ? 
                      beta= "is 'logit'",
                      "inverse.gamma" = "are '-1/mu' and 'log'"
    )
    mess <- paste0("(!) rand.family/link combination not handled;\nallowed link(s) for rand.family '",
                   rand.family$family,"' ",allowed)
    stop(mess)
  }
  lcrandfamfam
}


.checkRandLinkS <- function(rand.families) {
  for (it in seq_len(length(rand.families))) {
    rf <- rand.families[[it]]
    if (is.character(rf)) {
      rand.families[[it]] <- switch(tolower(rf),
                  gaussian = gaussian(),
                  gamma = Gamma(link=log), ## NOT the default link
                  beta = Beta(), 
                  "inverse.gamma" = inverse.Gamma(),
                  stop("rand.family argument not valid"))
    }
  }
  lcrandfamfam <- unlist(lapply(rand.families, .checkRandLink)) ## a _vector_ of lcrandfamfam := tolower(rand.family$family)
  unique.psi_M <- sapply(lcrandfamfam, switch,
                         gaussian = 0,
                         gamma = 1, 
                         beta = 1/2, 
                         "inverse.gamma" = 1
  )
  return(structure(rand.families,lcrandfamfam=lcrandfamfam,unique.psi_M=unique.psi_M))
}

.checkRespFam <- function(family, spaMM.=TRUE) {
  envname <- environmentName(environment(family))
  if ( ! envname %in% c("stats","spaMM","")) { ## " typically occurs when family has laready been checked...
    message(paste("family from namespace environment",envname,"possibly not correctly handled"))
    if (envname=="mgcv" && identical(paste(formals(family)$theta)[1],"stop")) { ## a call via family=family prevents testing the name
      mess <- "spaMM::negbin is masked by mgcv::negbin. Use explicitly 'family=spaMM::negbin', or assign 'negbin <- spaMM::negbin'"
      stop(mess)
    }
  }
  ## four lines from glm(), which should have no effect on processed$family, spec. those with a param.
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.language(family)) { # to bypass spaMM_Gamma by using family=quote(stats::Gamma(log))
    spaMM. <- (! length(grep("stats::",paste(family)[1L]))) ## spaMM. becomes FALSE if explicit quote(stats::...)
    family <- eval(family) 
  } 
  if (family$family=="Gamma" && spaMM.) {
    family <- spaMM_Gamma(link=family$link) 
  }
  return(family) ## input negbin(...) or COMPoisson(...) are returned as is => nu/shape unaffected
}

.eval_v_h_bounds <- function(cum_n_u_h, rand.families) {
  ## builds box constraints either NULL or non-trivial, of length n_u_h
  n_u_h <- cum_n_u_h[length(cum_n_u_h)]
  lower.v_h <- rep(-Inf,n_u_h)
  boxConstraintsBool <- FALSE
  gaussian_u_ranges <- integer(0)  
  for (it in seq_along(rand.families)) {
    if (rand.families[[it]]$family=="Gamma" && rand.families[[it]]$link=="identity") { ## gamma(identity)
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      lower.v_h[u.range] <- 1e-6 ## 1e-12 is disastrous
      boxConstraintsBool <- TRUE
    }
    if (rand.families[[it]]$family=="gaussian") gaussian_u_ranges <- c(gaussian_u_ranges,(cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])) 
  }
  if ( ! boxConstraintsBool ) lower.v_h <- NULL 
  boxConstraintsBool <- FALSE
  upper.v_h <- rep(Inf,n_u_h)
  for (it in seq_along(rand.families)) {
    if (rand.families[[it]]$family=="inverse.Gamma" && rand.families[[it]]$link=="-1/mu") { ## v ~ -Gamma
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      upper.v_h[u.range] <- -1e-06
      boxConstraintsBool <- TRUE
    } else if (rand.families[[it]]$family=="Gamma" && rand.families[[it]]$link=="log") {
      ## Gamma(log)$linkinv is pmax(exp(eta), .Machine$double.eps), ensuring that all gamma deviates are >= .Machine$double.eps
      ## we ensure that log(u_h) has symmetric bounds on log scale (redefine Gamma()$linkfun ?)
      upper.v_h <- Gamma(log)$linkfun(1/.Machine$double.eps)
    } 
  }
  if ( ! boxConstraintsBool ) upper.v_h <- NULL
  return(list(lower.v_h=lower.v_h,upper.v_h=upper.v_h))
}


# to convert back a family object to a call, EXCEPT for "multi". 
# get_param should be TRUE only in cases we are sure that the param has a value
.as_call_family <- function(family,get_param=FALSE) {
  # checkRespFam should have been run. 
  # Then we have either a call or an evaluated family object
  if (identical(family[[1L]],"multi")) {
    return(family) ## FR->FR tempo fix
  } else {
    as_call <- environment(family$aic)$mc ## shouldn't be called on multi...
    if (is.null(as_call)) { ## then we have a stats::  family 
      as_call <- call(family$family, link=family$link) ## need a substitute() ? 
    } else if (get_param) { ## param values are constant or dynamically assigned to processed$family
      famfam <- family$family 
      if (famfam == "negbin") {
        as_call$shape <- environment(family$aic)$shape
      } else if (famfam=="COMPoisson") {
        as_call$nu <- environment(family$aic)$nu
      }
    }
    return(as_call)
    ## avoids display of family def in fit results
  }
}

.calc_ZAlist <- function(Zlist, AMatrices) { 
  if ( length(Zlist) > 0L ) {
    if (length(AMatrices)) {
      Znames <- names(Zlist) ## explicit Znames necessary for prediction with re.form: Otherwise this should be OK:
      if (is.null(names(Zlist))) names(Zlist) <- Znames <- paste(seq_len(length(Zlist)))
      ## logic is Z[nresp,nUniqueRespLoc].A[nUniqueRespLoc,nHiddenv].L[nHiddenv,nHiddenv]
      for (char_rd in Znames) {
        if ( ! is.null(AMatrices[[char_rd]])) {
          colnams <- colnames(Zlist[[char_rd]])
          if ( ! setequal(rownames(AMatrices[[char_rd]]),colnams)) {
            stop(paste0("Any 'A' matrix must have row names that match the levels of the random effects\n",
                        "(i.e. the colnames of the 'Z' design matrix)"))
          }
          Zlist[[char_rd]] <- Zlist[[char_rd]] %*% AMatrices[[char_rd]][colnams,] # for the IMRF model, Z must be identity
        }
      }
      ## why this: ?? => removed
      # Anames <- names(AMatrices)
      # intnames <- as.integer(Anames)
      # intnames <- intnames[seq_len(length(AMatrices))]
      # intnames[is.na(intnames)] <- (seq_len(length(AMatrices)))[is.na(intnames)]
      # names(AMatrices) <- as.character(intnames)
      attr(Zlist,"AMatrices") <- AMatrices 
    }
  } 
  return(Zlist)
}

.evalWrapper <- function(object,
                         element, ## using try() has a visible impact on speed. Hence, the element must exist in the envir.
                         from=NULL) { ## exported for programming purposes
  if (  is.list(object)) {
    if (is.null(from)) from <- seq_len(length(object))
    if (length(from)>1L) {
      resu <- vector("list",length(from))
      for (id in seq_len(length(from))) {
        resu[[id]] <- eval(parse(text=element),envir=object[[id]])
      }
      names(resu) <- names(object[from])
      return(resu) ## a list, e.g a data list
    } else {
      object <- object[[from]]
    } 
  } 
  resu <- eval(parse(text=element),envir=object) ## try(eval(parse(text=element),envir=object),silent = TRUE)
  ## if (inherits(resu,"try-error")) resu <- NULL
  return(resu)
}

# .assignWrapper <- function(object,assignment) {
#   if (  is.list(object)) {
#     #    for (it in seq_len(length((object))) {eval(parse(text=paste0("object[[",it,"]]$",element," <- ",value)))} 
#     ## fails bc nam loses its enclosing "s :
#     #for (nam in names(object)) {eval(parse(text=paste0("object[[",nam,"]]$",element," <- ",value)))} 
#     for (nam in names(object)) {eval(parse(text=paste0("object[[\"",nam,"\"]]$",assignment)))} 
#   } else eval(parse(text=paste0("object$",assignment)))
#   ## no need to return the modified environment
# }
.assignWrapper <- function(object,assignment) { ## not using assign(), to allow eg. verbose['warn'] <- ... instead of simple variable name    # older version ".setProcessed" (<2.2.119) did not use de envir argument; this was tricky
  if (  is.list(object)) {
    for (it in seq_along(object)) eval(parse(text=assignment),envir=object[[it]]) 
  } else eval(parse(text=assignment),envir=object)
  ## no need to return the modified environment
}

.generateInitPhi <- function(formula,data,family,weights=NULL) {
  ## look for replicates to estimate phi
  form <- .subbarsMM(.asNoCorrFormula(formula)) ## removes all random effect decorum (converts (...|...) terms to some "+" form)
  # lhs
  mf <- model.frame(form,data=data)
  Y <- model.response(mf) ## evaluated rhs (e.g. log(y) rather than y...)
  # rhs
  off <- model.offset(mf)
  rhs_terms <- terms(.stripOffset(form),data=data) ## data argument allows e.g. ~ (.)^2; "Intercept" implicitly assumed by terms()
  rhs_terms <- stats::delete.response(rhs_terms) 
  mf <- model.frame(rhs_terms, data)
  # builds a model which indexes responses with identical predictor [ULI(mf)] 
  # and retains data that are replicates for each level of this index [uli]
  # but (g)lm complains about a model where uli has a single level [though this is meaningful]
  uli <- .ULI(mf) 
  # selection of data for replicates
  mf <- data.frame(y=Y,uli=as.factor(uli)) ## what's needed for both sides
  tuli <- table(uli)
  phiinfo <- which(tuli>1L) ## which levels of uli are informative
  ############################################
  # -> may generate a large # of levels -> of parameters to estimate -> glm takes time ! ## FR->FR may take time => quick patch    
  if (length(phiinfo)>30) return(0.1)  ## RETURN !
  ############################################
  whichRows <- uli %in% phiinfo
  mf <- mf[whichRows,] ## keep lines with replicates of all predictor variables
  off <- off[whichRows]
  if (!is.null(weights)) weights <- weights[whichRows]
  if (length(phiinfo)>1L) {
    locform <- y ~ uli 
  } else locform <- y ~ 1 ## only one level of uli had replicates
  # estimation of residual var
  if (NROW(mf)>1L) { ## if replicate observations available
    # formula must retain any operation on the lhs
    locglm <- spaMM_glm(formula=locform,data=mf,family=family,offset=off,weights=weights)  
    phi_est <- as.numeric(deviance(locglm)/locglm$df.residual)
  } else phi_est <- NULL
  return(phi_est)
}

## function for init.optim hence for fitme.
.get_init_phi <- function(processed,weights=NULL) {
  if (is.null(processed$envir$init_phi)) { 
    processed$envir$init_phi <- .generateInitPhi(formula=attr(processed$predictor,"no_offset"),data=processed$data,
                                family=processed$family,weights=weights) 
  } 
  return(processed$envir$init_phi)
}

.calc_fam_corrected_guess <- function(guess, For, processed, link_=NULL, nrand=NULL) {
  if (is.null(link_)) link_ <- processed$family$link
  if (is.null(nrand)) nrand <- length(processed$ZAlist)
  if (link_ != "identity") {
    if (For=="optim" || ## to avoid high initial values of lambda with spuriously high logLik by Laplace approx.
        processed$HL[1L]=="SEM") { 
      if (link_=="log") { ## test of family, not rand.family...
        fam_corrected_guess <- log(1.00001+guess/nrand) 
      } else {
        if (processed$bin_all_or_none) {
          maxinit <- 0.1 ##  a low init value is better even if final lambda estimate is high.
        } else maxinit <- 0.2 
        fam_corrected_guess <- min(guess, maxinit/nrand)
      }
    } else { ## iterative: allow larger values, but within some limits
      if (link_=="log") { ## test of family, not rand.family...
        fam_corrected_guess <- log(1.00001+guess/nrand) 
      } else {
        maxinit <- 2 ## even for processed$bin_all_or_none, test_all does not support lower value; but this is all based on ~ nothing
        fam_corrected_guess <- min(guess, maxinit/nrand)
      }
    }
  } else { ## identity link
    if (For=="optim") {
      maxinit <- 2
    } else {
      maxinit <- Inf
    }
    fam_corrected_guess <- min(guess, maxinit/nrand) 
  }
  return(fam_corrected_guess)
}

.eval_init_lambda_guess <- function(processed, stillNAs, ZAL=NULL, cum_n_u_h,For) {
  nrand <-  length(processed$ZAlist)
  if (is.null(processed$HLframes$Y)) { ## for resid model
    if (For=="optim") { 
      guess_from_glm_lambda <- 0.1 ## (FIXME ad hoc) Cannot let it NA as this will be ignored by further code, namely
    #                               optim_lambda_with_NAs[which_NA_simplelambda] <- init_lambda[which_NA_simplelambda]
    #                               init.optim$lambda <- optim_lambda_with_NAs[!is.na(optim_lambda_with_NAs)]
    } else guess_from_glm_lambda <- NA
  } else {
    inits_by_glm <- .get_inits_by_glm(processed) 
    if (For=="optim") { 
      guess_from_glm_lambda <- inits_by_glm$lambda*(3L*nrand)/((nrand+1L)) # +1 for residual
    } else guess_from_glm_lambda <- inits_by_glm$lambda*(3L*nrand+2L)/((nrand+1L)) # +1 for residual  ## old code was *5/(nr+1) ## f i x m e super ad hoc
  }
  init_lambda <- rep(NA,nrand)
  rand.families <- processed$rand.families
  lcrandfamfam <- attr(rand.families,"lcrandfamfam")
  corr_types <- processed$corr_info$corr_types
  link_ <- processed$family$link
  for (it in stillNAs) { ## fam_corrected_guess for each ranef in stillNAs
    if (is.null(ZAL)) {
      ZA <- processed$ZAlist[[it]]
    } else if (inherits(ZAL,"ZAXlist")) {
      ZA <- ZAL@LIST[[it]]
      if (inherits(ZA,"ZA_QCHM")) ZA <- ZA$ZA
    } else {
      u.range <- (cum_n_u_h[it]+1L):cum_n_u_h[it+1L]
      ZA <- ZAL[,u.range,drop=FALSE]
      
    }
    denom <- colSums(ZA*ZA) # colSums(ZA*ZA) are diag crossprod(ZA) (memo: length=ncol(ZA)). 
    denom <- denom[denom!=0] ## so that same result in dense and sparse if ZA has empty cols in sparse
    ZA_corrected_guess <- guess_from_glm_lambda/sqrt(mean(denom)) 
    #if (corr_types[it]=="AR1") ZA_corrected_guess <- log(1.00001+ZA_corrected_guess) ## ad hoc fix but a transformation for ARphi could be better FIXME
    fam_corrected_guess <- .calc_fam_corrected_guess(guess=ZA_corrected_guess, link_=link_, For=For, processed=processed, nrand=nrand)
    init_lambda[it] <- .preprocess_valuesforNAs(it, lcrandfamfam=lcrandfamfam, 
                                                rand.families=rand.families, init.lambda=fam_corrected_guess)
  }
  if (For != "optim") {   ## If called by HLfit: the present pmax() matters.
    init_lambda[stillNAs] <- pmax(init_lambda[stillNAs],1e-4) 
  } ## ELSE If called by fitme: calc_init_dispPars runs pmax on the whole vector
  return(init_lambda)
}

.preprocess_valuesforNAs <-  function(it, lcrandfamfam, rand.families, init.lambda) {
  if(lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="identity" && init.lambda==1) {
    adhoc <- 0.9999
  } else if(lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="log") {
    objfn <- function(lambda) {psigamma(1/lambda,1)-init.lambda}
    adhoc <- uniroot(objfn,interval=c(1e-8,1e8))$root
  } else if(lcrandfamfam[it]=="beta" && rand.families[[it]]$link=="logit") {
    #ad hoc approximation which should be quite sufficient; otherwise hypergeometric fns.
    objfn <- function(lambda) {8* lambda^2+3.2898*lambda/(1+lambda)-init.lambda}
    adhoc <- uniroot(objfn,interval=c(2.5e-6,1e8))$root
  } else if(lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="log") {
    ## this tries to controle var(v)<-init.lambda of v=log(u) by defining the right lambda for var(u)=lambda/(1-lambda)
    ## log (X~inverse.G) = - log (~Gamma), ie
    ##  log (X~inverse.G(shape=1+1/init.lambda,scale=1/init.lambda))~ - log (rgamma(shape=1+1/init.lambda,scale=1/init.lambda)
    ## where + log (~Gamma) has known mean (=...) and variance=trigamma(shape)
    ## (pi^2)/6=1.644934... is upper bound for trigamma(x->1) ie for lambda ->infty.
    ## however, there is something wrong with setting very large initial lambda values. Hence
    # if (init.lambda > 1.64491 ) { ## trigamma(1+1/100000) 
    #   adhoc <- 100000 
    if (init.lambda > 0.3949341 ) { ## trigamma(3) [for lambda=0.5]
      adhoc <- 0.5 
    } else { ## loer init.lambda -> lower lambda
      objfn <- function(lambda) {psigamma(1+1/lambda,1)-init.lambda}
      adhoc <- uniroot(objfn,interval=c(1e-8,1e8))$root
    }
    ## but the mean of v  is function of lambda and I should rather try to control the second moment of v
  } else if(lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="-1/mu") {
    adhoc <- (sqrt(1+4*init.lambda)-1)/2 # simple exact solution
  } else adhoc <- init.lambda
  adhoc
}

.reformat_lFix <- function(user_lFix, nrand) {
  lFix <- rep(NA,nrand)
  names(lFix) <- seq(nrand) # (character)
  if ( ! is.null(user_lFix)) {
    user_names <- names(user_lFix)
    # 'user'_names are not necessarily 1,2... (fitted values have names from rhs of ranef terms) 
    complex_names <- setdiff(user_names, paste(seq(nrand)))
    if (length(complex_names)) { ## cannot be used
      if (length(user_lFix)!=nrand) {
        stop("'fixed lambda' vector cannot be matched to random effects (different lengths).") ## needs to include NA's
      } else lFix[seq(nrand)] <- user_lFix ## keeping lFix default names
    } else if (length(user_names)) {
      lFix[user_names] <- user_lFix
    } else lFix[] <- user_lFix ## keep lFix names
  }
  return(lFix)
}

# function called within HLfit for missing inits... 
# but also sometimes in fitme_body prior to optimization (cf complex condition for call of .eval_init_lambda_guess())
.get_inits_by_glm <- function(processed, family=processed$family, reset=FALSE) {
  if (reset || is.null(processed$envir$inits_by_glm)) { 
    ## if .get_inits_by_glm is called prior to optimization the family parameters may not be assigned, so: 
    if (family$family=="negbin") {
      checktheta <- suppressWarnings(try(environment(family$aic)$shape,silent=TRUE))
      if (inherits(checktheta,"try-error")) family <- Poisson(family$link, trunc=environment(family$aic)$trunc) 
    } else if (family$family=="COMPoisson") {
      checknu <- suppressWarnings(try(environment(family$aic)$nu,silent=TRUE))
      if (inherits(checknu,"try-error")) family <- poisson() ## do not use the loglambda link of the COMPoisson!
    }
    y <- processed$y ## requested by the formula
    if (family$family=="binomial" && NCOL(y)==1L) { 
      BinomialDen <- processed$BinomialDen ## requested by the formula
      begform <-"cbind(y,BinomialDen-y)~"  
    } else {begform <-"y~"}
    ###################################################if (pforpv==0) {endform <-"0"} else 
    X.pv <- processed$AUGI0_ZX$X.pv ## possibly requested by the formula
    pforpv <- ncol(X.pv)
    if(pforpv) {
      endform <-"X.pv-1" ## pas besoin de rajouter une constante vue qu'elle est deja dans X
    } else {
      if (family$family %in% c("binomial","poisson")) {
        endform <- "1" ## no meaningful glm without fixed effect in this case !
      } else {endform <- "0"}
    }
    locform <- as.formula(paste(begform, endform))
    prior.weights <- processed$prior.weights   ## do not try to eval() it outside of the .wfit function call; else nasty crashes may occur.
    resu <- list() 
    if (family$family=="gaussian" && family$link=="identity") {
      if (inherits(X.pv,"sparseMatrix")) {
        resglm <- .spaMM_lm.wfit(x=X.pv,y=y,offset=processed$off,w=eval(prior.weights))
      } else resglm <- lm.wfit(x=X.pv,y=y,offset=processed$off,w=eval(prior.weights))
      fv <- fitted(resglm)
      dev <- resid(resglm)^2
      if (! is.null(resglm$weights)) dev <- dev/resglm$weights
      resu$lambda <- resu$phi_est <- sum(dev)/resglm$df.residual
    } else { ## GLM
      #
      resglm <- .tryCatch_W_E(glm.fit(x=X.pv, 
                                      y=processed$HLframes$Y, 
                                      weights = eval(prior.weights), 
                                      offset = processed$off, family = family, 
                                      control = processed[["control.glm"]]))$value 
      if (inherits(resglm,"error") || 
          ( ! resglm$converged && any(fitted(resglm)>1e30)) # this occurred in Gamma(log) models
          ) {
        if (TRUE) {
          #control.glm <- processed[["control.glm"]]
          #control.glm$maxit <- control.glm$maxit*2
          resglm <- spaMM_glm.fit(x=X.pv, 
                        y=processed$HLframes$Y, 
                        weights = eval(prior.weights), 
                        offset = processed$off, family = family, 
                        control = processed[["control.glm"]])
        } else {
          resglm <- withCallingHandlers(
            {
              #control.glm <- processed[["control.glm"]]
              #control.glm$maxit <- control.glm$maxit*2
              spaMM_glm.fit(x=X.pv, 
                            y=processed$HLframes$Y, 
                            weights = eval(prior.weights), 
                            offset = processed$off, family = family, 
                            control = processed[["control.glm"]])
            },
            warning = function(w){
              ## Ignore convergence diagnostics from *this* spaMM_glm.fit() call (only, as it is the first call within HLfit)
              #if(w$message != "glm.fit: algorithm did not converge" ) # the message is never this !
              warning(w$message)
              invokeRestart("muffleWarning")
            } 
          )
        }
        if ( ( ! identical(spaMM.getOption("spaMM_glm_conv_silent"),TRUE))
             && (conv_crit <- environment(spaMM_glm.fit)$spaMM_glm_conv_crit$max>0)) {
          resu$conv_info <- paste(".get_inits_by_glm() -> spaMM_glm.fit did not yet converge at iteration",resglm$iter,"(criterion:",
                                paste(names(conv_crit),"=",signif(unlist(conv_crit),3),collapse=", "),")")
          assign("spaMM_glm_conv_crit",list(max=-Inf) , envir=environment(spaMM_glm.fit)) # So that no-convergence in these glms will not be warned about
        }
      }
      #
      resu$phi_est <- as.numeric(deviance(resglm)/resglm$df.residual) 
      if (family$family=="binomial" && max(resglm$prior.weights)==1L) { ## binary response
        resu$lambda <- 1
      } else {
        fv <- fitted(resglm)
        # resu$lambda <- sum(resid(resglm)^2/(resglm$prior.weights*family$variance(fv)))/resglm$df.residual # was a birth pang...
        # resid(resglm) gives resglm$residuals, which are the "working residuals" ie (y-mu)/mu.eta one the RHS of the IRLS equation
        # hence they are not an appropriate way to reconstruct the deviance residuals
        # which are given by with(resglm, family$dev.resids(y, mu=fitted.values, weights)), already taking prior weights into account
        # but this leads to 
        resu$lambda <- resu$phi_est
      }
    } ## end else ('GLM')
    if (pforpv) {
      ## Two potential problems (1) NA's pour param non estimables (cas normal); 
      ## (2) "glm.fit: fitted probabilities numerically 0 or 1 occurred" which implies separation or large offset
      if (max(abs(c(coefficients(resglm))),na.rm=TRUE)>1e10) { ## na.rm v1.2 
        warning(paste0("(!) Apparent divergence of estimates in a *GLM* analysis of the data.\n",
                       "    Check your data for extreme values, separation or bad offset values."))
      } 
      beta_eta <- c(coefficients(resglm)) ## this may include NA's. Testcase: HLfit(Strength ~ Material*Preheating+Method,data=weld)
      if (all(names(beta_eta)=="X.pv")) { ## si la formula etait y ~X.pv-1
        names(beta_eta) <- colnames(resglm$model$X.pv)
      } else names(beta_eta) <- gsub("^X\\.pv","",names(beta_eta)) ## removes initial "X.pv" without guessing any order or length
      resu$beta_eta <- beta_eta
    } 
    #parent.env(environment(processed$get_inits_by_glm)) <- environment(stats::glm)
    processed$envir$inits_by_glm <- resu
  } 
  return(processed$envir$inits_by_glm)
} 

.preprocess_HLmethod <- function(HLmethod, family, lcrandfamfam) {
  ## conversion from user-friendly format to standard 'XX(...)' format
  ## first index is for (0) h / (1) p_v(h) / (2) p^s_v(h) ie whether h lik or marginal lik is used for fixed effect estimation
  ## Other components determine three options wrt to leverages, some for ReML correction, some for notEQL.
  ## ML() vs HL() determines whether the hatval leverages are computed ie whether some basic ReML correction in applied
  ## second index is for further ReML correction: no (0) or yes (1) D log Det / d log lambda correction (2) further p^s_bv correction
  ## third index is for use of (0) EQL deviance residuals (this affects the leverages) or not (1) (this is not an ReML correction... but impatcs only dispersion estimation..)
  ## thus overall we have <ReML/not>( <h/l> , <more ReML/not> , <not/EQL> )
  ## NohL07 table 1 has interesting terminology and further tables show even more cases
  # comme l'EQL est dans un monde quasi gaussien, il se ramene aux LMM et utilise les leverage standard,
  # Pour les GLMM non LMM, Dans la mesure ou HL(0,.) utilise h et non q+, ce n'est pas l'EQL pour les fixed params
  # Dans la mesure ou on utilise les leverages standard pour les dispersion param, c'est de l'EQL
  if (HLmethod=="ML") {
    HLmethod <- "ML(1,1,1)"  
  } else if (HLmethod=="SEM") {
    #if ( ! requireNamespace("probitgem",quietly = TRUE)) {## will pass CHECK if CRAN knows probitgem
    # if ( ! ("probitgem" %in% .packages()) ) { ## passed CRAN checks; appropriate sif probitgem explicitly loaded
    if ( ! ("probitgem" %in% .packages(all.available = TRUE)) ) { ## 
      stop("Package 'probitgem' not available for fitting by SEM.")
    } else do.call("require", list(package="probitgem", quietly = TRUE)) ## passes CRAN checks and transparent...
    #      eval(as.call(c(quote(require),list(package="probitgem", quietly = TRUE)))) ## passes CRAN checks (but temporary)
    #      eval(as.call(c(quote(requireNamespace),list(package="probitgem", quietly = TRUE)))) ## is not sufficient
    HLmethod <- "ML('SEM',NA,NA)" 
    if ( ! (family$family=="binomial" && family$link=="probit")) {
      stop("SEM is applicable only to binomial(probit) models.")
    }
  } else if (HLmethod=="REML") {
    HLmethod <- "HL(1,1,1)"  
  } else if (HLmethod=="REPQL" || HLmethod=="PQL") {
    if (any(lcrandfamfam!="gaussian"))  stop("PQL is not defined for HGLMs in general. Do you mean 'EQL-'?") 
    HLmethod <- "HL(0,0,1)" ## (=REPQL, equivalent to HL(0,0,0) ie EQL- for GLMMs )
  } else if (HLmethod=="PQL/L") { ## again no D log Det / d log lambda correction
    if (any(lcrandfamfam!="gaussian"))  stop("PQL is not defined for HGLMs in general. Do you mean 'EQL-'?") 
    HLmethod <- "ML(0,0,1)" ## (equivalent to ML(0,0,0) for GLMMs)
  } else if (HLmethod=="EQL-") { ## version LeeNP06 p.212 incomplete 1st order ## probably hglm package
    ## thus overall we have <ReML->HL >( <h->0> , <not more ReML->0> , <EQL -> 0> )
    HLmethod <- "HL(0,0,0)" ## 
  } else if (HLmethod=="EQL+") { ## version LeeN01 complete 1st order
    ## thus overall we have <ReML->HL >( <h->0> , <more ReML->1> , <EQL -> 0> )
    HLmethod <- "HL(0,1,0)" ## (0,...): gaussianise everything, hence no a(1) correction ## there is no HL(1) a(1) correction in GLM.MME
  }
  return(HLmethod)
}


.preprocess_SEMargs <- function(BinomialDen, nobs, control.HLfit, y) {
  if (sum(BinomialDen) != nobs) {
    stop("(!) SEM procedure: the data do not seem binary; non-binary data are not handled.")
    
    # repsucces <- rep(seq(nrow(currentSample)),currentSample$succes)
    # repechec <- rep(seq(nrow(currentSample)),currentSample$echec)
    # currentSample <- currentSample[c(repsucces,repechec),]
    # currentSample$echec <- rep(c(0,1),c(length(repsucces),length(repechec)))
    # currentSample$succes <- 1L-currentSample$echec
    
  }
  SEMseed <- control.HLfit$SEMseed ##  
  # if (is.null(SEMseed)) SEMseed <- 123 ## OK pour SEM *unique* mais remplace par NULL dans optimthroughSmooth
  SEMargs <- list(SEMseed=SEMseed)
  SEMargs$nMCint <- control.HLfit$nMCint ##  as is => SEM procedure must handle null value
  SEMargs$get_SEM_info <- control.HLfit$get_SEM_info ##  as is => SEM procedure must handle null value
  SEMargs$control_pmvnorm$maxpts <- control.HLfit$pmvnorm_maxpts ##  as is => SEM procedure must handle null value
  #SEMlogL <- control.HLfit$SEMlogL
  #if (is.null(SEMlogL)) SEMlogL <- "pMVN" ## FR->FR 1.8.25 !
  #SEMargs$SEMlogL <- SEMlogL
  SEMargs$SEMlogL <- control.HLfit$SEMlogL
  nSEMiter <- control.HLfit$nSEMiter ##  
  if ( (! is.null(nSEMiter)) && nSEMiter < 10) {
    stop(" 'nSEMiter' should be >9")
  } else SEMargs$nSEMiter <- nSEMiter ## default given by probitgem, typic 200 
  SEMargs$ngibbs <- control.HLfit$ngibbs ## default given by probitgem, typic 20 
  SEMargs$SEMsample <- control.HLfit$SEMsample ## stays NULL if NULL
  SEMargs$whichy1 <- (y==1) 
  return(SEMargs)
}

.is_link_canonical <- function(family) {
  canonicalLink <- FALSE
  if (family$family=="gaussian" && family$link=="identity") {
    canonicalLink <- TRUE
  } else if (family$family=="poisson" && family$link=="log") {
    canonicalLink <- TRUE
  } else if (family$family=="binomial" && family$link=="logit") {
    canonicalLink <- TRUE
  } else if (family$family=="Gamma" && family$link=="inverse") {
    canonicalLink <- TRUE
  } else if (family$family=="COMPoisson" && family$link=="loglambda") {
    canonicalLink <- TRUE
  } ## no implemented canonical link case for negbin
  return(canonicalLink)
}

.assign_canonLink_G_LMMbool <- function(family, processed) {
  processed$canonicalLink <- .is_link_canonical(family)  
  #
  GLMMbool <- (length(processed$lcrandfamfam) && all(processed$lcrandfamfam=="gaussian") ) ## only allowed gaussian rand.family is gaussian(identity) 
  if (GLMMbool) {
    LMMbool <- (family$family=="gaussian" && processed$canonicalLink) 
  } else LMMbool <- FALSE
  processed$LMMbool <- LMMbool
  processed$GLMMbool <- GLMMbool
  # no return value
}

.calc_rankinfo <- function(X.pv, verbose=FALSE, tol) {
  nc <- ncol(X.pv)
  nobs <- nrow(X.pv)
  tol <- eval(tol)
  if (inherits(X.pv,"sparseMatrix")) {
    if (nc>nobs) { ## qr does not handle "wide " matrices
      zerofill <- sparseMatrix(dims = c(nc-nobs,nc), i={}, j={})
      if (verbose) print("qr(rbind2(<sparseMatrix>,zerofill))")
      rankinfo <- qr(x=rbind2(X.pv,zerofill),tol=tol)
    } else {
      if (verbose) print("qr(<sparseMatrix>)")
      rankinfo <- qr(x=X.pv,tol=tol)
    }
  } else  {
    if ((nc^2)*(max(nc,nobs)-nc/3)>1e10) {
      message("The design matrix for fixed effects is quite large. Checking its rank will take time. 
              See help('rankinfo') for possible ways to circumvent this.")
    }
    ## .rankinfo request much memory, but is clearly faster for matrices of size ~ 5000*2900
    ## HOWEVER it does not always yield *p_bv* results consistent with qr !!!!!
    ## +
    ## FIXME this may evolve when  ColPivHouseholderQR implements "blocking strategies".
    ## Waiting for such improvements, spaMM.getOption("rankMethod") is "qr" by default.
    rankmethod <- spaMM.getOption("rankMethod")
    if (is.null(rankmethod)) {
      .rankinfo_usable <- (prod(dim(X.pv))<1e8)
      rankmethod <- c("qr",".rankinfo")[.rankinfo_usable+1L] 
    }
    if (verbose)  print(paste0(rankmethod,"(<matrix>)"))
    rankinfo <- do.call(rankmethod,list(x=X.pv,tol=tol))
  }
  if (inherits(rankinfo,"qr")) {
    rankinfo <- list(rank=rankinfo$rank, 
                     whichcols=sort(rankinfo$pivot[1:rankinfo$rank]), ## cf lm.fit or base::qr.coef. sort() is cosmetic, \neq sort.list() !
                     method="qr") 
  } else if (inherits(rankinfo,"sparseQR")) {
    checkd <- (abs(diag(rankinfo@R))>tol) ## 1:rankinfo$rank in base::qr; but here the zero-filled rows are not the last of R
    whichcols <- (rankinfo@q+1L)[checkd] ## $pivot[ <checkd> ] in base::qr
    rankinfo <- list(rank=sum(checkd), whichcols=sort((rankinfo@q+1L)[checkd]), method="sparseQR")
  } else if (rankinfo$method==".rankinfo") { ## .rankinfo:
    stop("FIXME possibly missing the 'checkd' line required in sparseQR method.")
    checkd <- NaN ## suspect that not 1:rankinfo$rank
    piv <- (rankinfo$pivot+1L)[checkd] ## +1L bc .rankinfo()$pivot indexes columns from 0
    rankinfo <- list(rank=rankinfo$rank, whichcols=sort((rankinfo$pivot+1L)[checkd]), method=".rankinfo")
  } else stop("Unknown rankinfo$method.")
  return(rankinfo)
}

as_precision <- function(corrMatrix) {
  if (inherits(corrMatrix,"dist")) {
    corrMatrix <- as.matrix(corrMatrix) ## leaves 0 on the diagonal! 
    diag(corrMatrix) <- 1 ## so that m is true correlation matrix 
  }
  #precmat <- try(chol2inv(chol(corrMatrix)),silent=TRUE) # this can succeed and produce an inverse non positive def...
  #if (inherits(precmat,"try-error")) { # Don' try solve() without regularization as this may produce a matrix with quite large negative eigenvalues.
    esys <- eigen(corrMatrix, only.values = TRUE)
    evalues <- esys$values
    min_d <- evalues[1L]/1e14 ## so that corrected condition number is at most the 1e14
    diagcorr <- max(c(0,min_d-evalues)) # SINGLE SCALAR
    if (diagcorr>0) {
      corrMatrix <- (1-diagcorr)*corrMatrix
      diag(corrMatrix) <- diag(corrMatrix) + diagcorr ## # all diag is corrected => added a constant diagonal matrix 
    }
    # problem: forceSymmetric() after the solve can strongly perturb the eigenvalues; solve(forceSymmetric()) is much better than forceSymmetric(solve())
    precmat <- solve(forceSymmetric(corrMatrix))
  #} else precmat <- forceSymmetric(precmat)
  colnames(precmat) <- rownames(precmat) <- colnames(corrMatrix) 
  # drop0 useful to convert to sparseMarix even if no 0 to drop
  return(structure(list(matrix=drop0(precmat)),class=c("list","precision"))) # return value must be sparse, not simply Matrix. 
}

.get_corr_prec_from_covStruct <- function(covStruct,it) { ## compatible with spaMM3.0 extended syntax
  if (length(covStruct)>1L) {
    corrMatrix <- covStruct[it][["corrMatrix"]]
  } else {
    corrMatrix <- covStruct[["corrMatrix"]]
  }
  if (is.null(corrMatrix)) {
    matrix_ <- covStruct[["precision"]]
    if (is.null(matrix_)) stop("missing covariance structure for corrMatrix model") ## no info in any form
    ## else the user provided covStruct$precision =>  reformattng 
    corrMatrix <- structure(list(matrix=matrix_),class=c("list","precision")) ## so that further code can detect its a precision matrix
  } else  corrMatrix <- corrMatrix
  return(corrMatrix)
}

.get_adjMatrix_from_covStruct <- function(covStruct,it) { ## compatible with spaMM3.0 extended syntax
  if (length(covStruct)>1L) {
    adjMatrix <- covStruct[it][["adjMatrix"]]
  } else {
    adjMatrix <- covStruct[["adjMatrix"]]
  }
  if (is.null(adjMatrix)) stop("missing 'adjMatrix' for adjacency model") ## or SAR_WWt model...
  return(adjMatrix)
}

.assign_cov_matrices__from_covStruct <- function(corr_info, covStruct=NULL, corrMatrix=NULL, adjMatrix=NULL) {
  if ( ! is.null(covStruct)) covStruct <- .preprocess_covStruct(covStruct)
  corr_info$AMatrices <- attr(covStruct,"AMatrices")
  corr_types <- corr_info$corr_types
  corr_info$adjMatrices <- vector("list",length(corr_types))
  corr_info$corrMatrices <- vector("list",length(corr_types))
  for (it in seq_along(corr_types)) {
    corr_type <- corr_types[it]
    if ( ! is.na(corr_type)) {
      if (corr_type=="adjacency" || corr_type=="SAR_WWt") {
        if ( is.null(adjMatrix) ) {
          corr_info$adjMatrices[[it]] <- .get_adjMatrix_from_covStruct(covStruct,it)
        } else corr_info$adjMatrices[[it]] <- adjMatrix   
      } else if (corr_type=="corrMatrix") {
        if (is.null(corrMatrix)) {
          corr_info$corrMatrices[[it]] <- .get_corr_prec_from_covStruct(covStruct,it) 
        } else corr_info$corrMatrices[[it]] <- corrMatrix
        .check_corrMatrix(corr_info$corrMatrices[[it]]) 
      } 
      # IMRF AMatrices are assigned later from Zlist info, not from covStruct info
    }
  }
}

.determine_sparse_X <- function(HLframes, X.pv) {
  sparse_X <- spaMM.getOption("sparse_X") 
  ## forcing sparse_X may (1) be slow for small problems 
  ## (2) entails the use of Matrix::Cholesky, which is less accurate => small bu visible effect on predVar in singular 'twolambda' case
  if (is.null(sparse_X)) {
    # les terms et cols sont reord selon attr(fixef_off_terms, "term.labels")
    # (1) identify terms that involve factors: 
    if ( length(HLframes$fixef_levels) ) {
      vars_terms_table <- attr(HLframes$fixef_off_terms,"factors") ## ""factors"" is a misleading name as table includes quantitative predictors
      n_levels <- sapply(HLframes$fixef_levels,length)
      col_heuristic_denseness <- rep(1,ncol(X.pv))
      term_heuristic_denseness <- rep(1,ncol(vars_terms_table)) 
      names(term_heuristic_denseness) <- colnames(vars_terms_table)
      for (term in colnames(vars_terms_table)) { ## the cols of vars_terms_table should match the terms in the order used by attr(X.pv,"assign")...
        vars_in_term <- names(which(vars_terms_table[,term]>0))
        factors_in_terms <- intersect(names(HLframes$fixef_levels),vars_in_term)
        if (any(vars_terms_table[,term]>0)) term_heuristic_denseness[term] <- 1/prod(n_levels[factors_in_terms])
      }
      asgn <- attr(X.pv,"assign") ## "for each column in the matrix ... the term in the formula which gave rise to the column"
      for (it in seq_along(asgn)) if (asgn[it]>0L) col_heuristic_denseness[it] <- term_heuristic_denseness[asgn[it]]
      sparse_X <- (mean(col_heuristic_denseness)<0.11) ## FIXME not enough tests of threshold; could use data.test in test-predVar which has mean density=0.19
    } else sparse_X <- FALSE
  }
  return(sparse_X)
}


.process_corr_info_spprec <- function(corr_types, corr_info, For, sparse_precision) {
  for (it in seq_along(corr_types)) {
    corr_type <- corr_types[it]
    if ( ! is.na(corr_type)) {
      if (corr_type=="corrMatrix" && sparse_precision) {
        if (! inherits(corr_info$corrMatrices[[it]],"precision")) {
          message("Conversion of corrMatrix to precision matrix")
          corr_info$corrMatrices[[it]] <- as_precision(corr_info$corrMatrices[[it]]) ## test by test-Matern-spprec
        }
      } else if (corr_type=="adjacency" || corr_type=="SAR_WWt") {
        if (For %in% c("fitme","corrHLfit")) {
          ## This block cannot be in .assign_cov_matrices__from_covStruct bc it requires sparse_precision value
          # We will need the $d to determine rhorange.
          # given the decomp is computed, we try to make it available for other computations.
          # But others computs should not assume it is available.
          # HLCor_body can also compute "symSVD" attribute and add in each processed environment
          decomp <- .provide_AR_factorization(corr_info$adjMatrices[[it]], sparse_precision, corr_type)
          if (corr_type=="SAR_WWt") attr(corr_info$adjMatrices[[it]],"UDU.") <- decomp
          if (corr_type=="adjacency") {
            corr_info$adjMatrices[[it]] <- forceSymmetric(drop0(corr_info$adjMatrices[[it]]))
            attr(corr_info$adjMatrices[[it]],"symSVD") <- decomp
          } 
        }
      }
    }
  }
}


.preprocess_control.dist <- function(control.dist,corr_types) {
  # handles either a list with possible elements "dist.method","rho.mapping", or a nested list with elements for the different ranefs
  if ( ! is.null(control.dist)) {
    if (length(control.dist)) { # list with elements (or NULL_s_)
      matchnames <- intersect (names(control.dist), c("dist.method","rho.mapping"))
      if (length(matchnames)) {
        control_dist <- vector("list",length(corr_types))
        for (rd in seq_along(corr_types)) {
          if (corr_types[[rd]] %in% c("Matern","Cauchy","AR1", "IMRF")) { ## not (yet) clearly useful for AR1 ?
            control_dist[[as.character(rd)]] <- control.dist[matchnames]
          }
        }
        return(control_dist) # return nested list with elements for the different ranefs built from the 'matchnames'
      } else return(control.dist) ## No 'matchnames': control.dist should then be a nested list with elements for the different ranefs 
    } else return(vector("list", length(corr_types))) # Return "empty nested list" built from empty list() (could return NULL ?)
  } else return(NULL)
}


.determine_augZXy <- function(processed, control.HLfit) {
  augZXy_cond_inner <- augZXy_cond <- spaMM.getOption("allow_augZXy")
  if (is.null(augZXy_cond)) {
    augZXy_cond <- (processed$models[["phi"]]=="phiScal") ## allows prior weights, but there must be a single phi scaling factor
    # => default is FALSE if phi is fixed but can be overcome by allow_augZXy=TRUE for devel purpose at least
    augZXy_cond_inner <- TRUE ## for .makeCovEst1() ## I have not thought about ths case with [phi fixed but overcome by allow_augZXy=TRUE]
  }
  # Conditions common to outer and inner optim
  if (augZXy_cond_inner) augZXy_cond_inner <-  processed$LMMbool
  if (augZXy_cond_inner) augZXy_cond_inner <- ( is.null(processed$X.Re) || ! ncol(processed$X.Re)) ## exclude non-standard REML (avoiding NCOL(NULL)=1)
  augZXy_cond <- (augZXy_cond_inner && augZXy_cond) 
  # Conditions specific to outer optim
  if (augZXy_cond) augZXy_cond <- (processed$For=="fitme")
  if (augZXy_cond) augZXy_cond <- (is.null(control.HLfit$intervalInfo)) 
  if (augZXy_cond) augZXy_cond <- (length(processed$init_HLfit)==0L) ## excludes (among others) internal rho estimation
  ## we need to detect numbers in lambda, and NaN in lambda (forced inner optimization: in the latter case,
  #                                               we could imagine using the augZXy code for given phi within HLfit)
  if (augZXy_cond) augZXy_cond <- ! any( !is.na(c(processed$lambda.Fix)) | is.nan(c(processed$lambda.Fix)))
  if (augZXy_cond) augZXy_cond <- (! any(processed$ranCoefs_blob$is_set)) ## at preprocessing stage this excludes cases with ranCoefs sets by users
  if (augZXy_cond) augZXy_cond <- ( ! is.call(processed$prior.weights) && attr(processed$prior.weights,"unique")) ## cf APHLs_by_augZXy code
  return(structure(augZXy_cond, inner=augZXy_cond_inner))
}

# fails on priorw in Infusion (long tests)
# .preprocess_pweights <- function(prior.weights, validdata) {
#   subs_p_weights <- substitute(prior.weights)
#   if ( ! (inherits(subs_p_weights,"call") && subs_p_weights[[1L]] == "quote") )  {
#     if (is.null(prior.weights)) {
#       prior.weights <- structure(rep(1L,nrow(validdata)),unique=TRUE) ## <- 1L prevented by glm -> model.frame(... prior.weights)
#     } else {
#       prior.weights <- as.vector(stats::model.weights(validdata)) ## as.vector as in say lm() protects against array1d
#       if ( ! is.numeric(prior.weights)) 
#         stop("'weights' must be a numeric vector")
#       if (any(prior.weights < 0)) 
#         stop("negative weights not allowed")
#       attr(prior.weights,"unique") <- (length(unique(prior.weights))==1L) 
#       #attr(prior.weights,"only1") <- all(upw==1L)
#     }
#   } else attr(prior.weights,"unique") <- FALSE ## when 'prior.weights' is a quoted expression   
#   return(prior.weights)
# }


.reformat_resid_model <- function(resid.model
                                  #  check_old_syntax for back compatibility,
                                  # to check whether control.HLfit$resid.family was used, when resid.model is only a formula
                                  # This use of control.HLfit is no longer documented
                                  ,check_old_syntax=NULL) { 
  fixed <- as.list(resid.model$fixed) ## converts NULL to list() as exp'd for 'fixed' in fitme_body()
  if ( ! is.null(resid.model)) { 
    if ( ! is.null(form <- resid.model$formula)) {
      resid.model$formula <- .preprocess_formula(form)
      # control of phi matters for phiHGLM but the phi model has not yet been determined
      if (is.null(fixed[["phi"]])) {
        fixed[["phi"]] <- c(default=1)
        #        message("'phi' of residual dispersion model set to 1 by default") ## inappropriate when resid.model=~1
      } else if (is.na(fixed[["phi"]])) fixed[["phi"]] <- NULL ## to force estimation of this phi; 
    } else { ## including resid.model=~1
      resid.model <- list(formula=.preprocess_formula(resid.model), ## otherwise it retains the local environment of the fn in which it is match.call()ed!
                          family=check_old_syntax) # may be NULL, in which case it is later processed as spaMM_Gamma(log)
    }
    resid.model$fixed <- fixed
    if (is.null(resid.model$resid.model)) resid.model$resid.model <- ~1
  } 
  return(resid.model)
}

.rankTrim <- function(X.pv, rankinfo, verbose=FALSE) {
  if (verbose)  str(rankinfo)
  if (rankinfo$rank < ncol(X.pv)) {   
    X.pv <- structure(X.pv[,rankinfo$whichcols,drop=FALSE], namesOri=colnames(X.pv)) 
    # etaFix$beta |         variables 
    #             |  valid vars | invalid vars
    #     (1)           (2)           (3)
    # (2): colnames(<HLfit>$envir$beta_cov_info$beta_cov) = colnames (<HLfit>$X.pv)
    # (1+2+3): namesOri, names(<HLfit>$fixef)
  } else attr(X.pv,"namesOri") <- colnames(X.pv)  
  attr(X.pv,"rankinfo") <- rankinfo
  return(X.pv)
}

.post_process_X <- function(X.pv, HL, HLframes, control.HLfit) {
  Xattr <- attributes(X.pv)
  if ( ncol(X.pv) && HL[1L] != "SEM") { ## gaussian bc .get_inits_by_glm -> spaMM_glm -> stats::glm fails
    sparse_X <- .determine_sparse_X(HLframes, X.pv)   ## sparse_X is useful for rankinfo bc Matrix::qr can be much faster
    if (sparse_X) {
      X.pv <- as(X.pv,"dgCMatrix") # .Rcpp_as_dgCMatrix(X.pv) # 
    }
  } 
  # determine true # cols before determining sparse_precision
  colnames(X.pv) <- colnames(HLframes$X)
  if (ncol(X.pv)) {
    rankinfo <- control.HLfit$rankinfo
    if (is.null(rankinfo)) rankinfo <- .calc_rankinfo(X.pv, tol=spaMM.getOption("rankTolerance"))
    if (is.list(rankinfo)) X.pv <- .rankTrim(X.pv,rankinfo = rankinfo)
  }
  names_lostattrs <- setdiff(names(Xattr), c(names(attributes(X.pv)),"dim","dimnames"))
  attributes(X.pv)[names_lostattrs] <- Xattr[names_lostattrs] # as in .subcol_wAttr(). 
  return(X.pv)
}

.preprocess_formula <- function(formula) {
  if (inherits(formula,"predictor")) { 
    return(formula) ## happens eg in confint
    # stop("Do not call '.preprocess_formula' on a predictor object.")
  }
  formlen <- length(formula)
  formula[[formlen]] <- .expandDoubleVerts(formula[[formlen]])
  rhs <- .expand_multIMRFs(formula[[formlen]]) 
  hyper_info <- attr(rhs,"hyper_info")
  if ( !is.null(hyper_info)) {
    hyper_form <- formula 
    environment(hyper_form) <- new.env()
    hyper_info$formula <- hyper_form
    attr(rhs,"hyper_info") <- NULL
  }
  formula[[formlen]] <- rhs
  environment(formula) <- new.env() # zut <- y~x; rezut <- spaMM:::.preprocess_formula(zut); ls(environment(zut)) confirms that the original 'formula's environment is unaffected
  predictor <- structure(formula,hyper_info=hyper_info) 
  class(predictor) <- c("predictor",class(formula))
  return(predictor)
}

.calc_vec_normIMRF <- function(exp_ranef_terms, corr_types) {
  nrand <- length(exp_ranef_terms)
  vec_normIMRF <- rep(FALSE, nrand)
  for (rd in seq_along(exp_ranef_terms)) { 
    corr_type <- corr_types[rd]
    if (! is.na(corr_type) && corr_type== "IMRF" ) {
      useNorm <- attr(attr(exp_ranef_terms[[rd]],"type"),"pars")$no 
      if (is.null(useNorm)) {
        stop("is.null(useNorm)")  # DEBUGGING
      } else vec_normIMRF[rd] <- useNorm
    }
  }
  return(vec_normIMRF)
}

.add_ZAfix_info <- function(AUGI0_ZX, ZAlist, sparse_precision, corr_types, processed) {
  if ( ! any(AUGI0_ZX$vec_normIMRF)) {
    ZAfix <- .ad_hoc_cbind(ZAlist, as_matrix=FALSE)  
    if (sparse_precision) {
      if ( ! inherits(ZAfix,"sparseMatrix")) ZAfix <- as(ZAfix,"dgCMatrix") # .Rcpp_as_dgCMatrix(ZAfix) ## 
      rsZA <- rowSums(ZAfix) ## test that there a '1' per row and '0's otherwise:  
      AUGI0_ZX$is_unitary_ZAfix <- (unique(rsZA)==1 && all(rowSums(ZAfix^2)==rsZA)) ## $ rather than attribute to S4 ZAfix
    } else {
      if (.eval_as_mat_arg(processed)) ZAfix <- as.matrix(ZAfix)  
    }
    AUGI0_ZX$ZAfix <- ZAfix # Used in code for ZAL in HLfit_body(); and extensively to fit  by spprec. 
                            # Special care is required when A is later modified (as by .calc_normalized_ZAlist())
  }
  return(AUGI0_ZX)
}


.preprocess <- function(control.HLfit, ranFix=NULL, HLmethod, 
                       predictor, resid.model,
                       REMLformula, data, family,
                       BinomialDen, rand.families, etaFix, prior.weights,
                       objective=NULL, 
                       control.glm, adjMatrix=NULL, verbose=NULL, For,
                       init.HLfit=list(),corrMatrix=NULL,covStruct=NULL,
                       uniqueGeo=NULL,distMatrix=NULL, 
                       control.dist=NULL 
                       ) {
  callargs <- match.call() 
  #
  ################ handling list of data #######################
  if (inherits(data,"list")) {
    locargs <- as.list(callargs)
    famfam <- family$family
    processed <- lapply(data, function(dd) {
      locargs$data <- dd
      if ( ! is.null(famfam) && famfam=="multi") locargs$family <- family$binfamily  
      eval(as.call(locargs)) ## call("preprocess",...) on each data set
    })
    return(processed) ## a list of environments
  }
  ###############################################################
  resid.model <- .reformat_resid_model(resid.model,check_old_syntax=control.HLfit$resid.family) ## calls .preprocess_formula() ## the list(...) is used even for poisson, binomial...
  # remove rows with NA's in required variables:
  validdata <- .getValidData(formula=predictor, resid.formula=resid.model$formula,
                            data=data, callargs=callargs["prior.weights"]) ## OK for all prior weights
  if ( inherits(data,"data.frame")) {
    data <- data[rownames(validdata),,drop=FALSE]  
  } else if  (inherits(data,"environment")) {
    data <- validdata
  } else {
    stop("'data' is not a data.frame not an environment.")
  }
  # add easily testable family name
  famfam <- family$family
  processed <- list2env(list(data=data, family=family, For=For,
                             envir=list2env(list(),parent=environment(HLfit)))) 
  processed$clik_fn <- .get_clik_fn(family)
  #
  stop.on.error <- control.HLfit$stop.on.error ##  
  if (is.null(stop.on.error)) stop.on.error <- FALSE
  processed$stop.on.error <- stop.on.error ##  
  break_conv_logL <- control.HLfit$break_conv_logL ##whether to stop if logL (p_v) appears to have converged  
  if (is.null(break_conv_logL)) break_conv_logL <- FALSE
  processed$break_conv_logL <- break_conv_logL ##  
  ## numerical control parameters 
  spaMM_tol <- spaMM.getOption("spaMM_tol") 
  if ( ! is.list(spaMM_tol)) stop("spaMM_tol must be a list")
  # then user's explicit conv.threshold controls
  if ( ! is.null(conv.threshold <- control.HLfit$conv.threshold)) spaMM_tol[["Xtol_rel"]] <- conv.threshold
  # then user's explicit spaMM_tol controls
  if ( ! is.null(user_spaMM_tol <- control.HLfit$spaMM_tol)) spaMM_tol[names(user_spaMM_tol)] <- user_spaMM_tol
  processed$spaMM_tol <- spaMM_tol
  #
  iter.mean.dispFix <- control.HLfit$iter.mean.dispFix ## private control
  if (is.null(iter.mean.dispFix)) iter.mean.dispFix <- control.HLfit$max.iter.mean ## public control
  if (is.null(iter.mean.dispFix)) iter.mean.dispFix <- 200 ## control of inner loop when no disp param is estimated ## was 40, 06/2014
  processed$iter.mean.dispFix <- iter.mean.dispFix  
  #
  iter.mean.dispVar <- control.HLfit$iter.mean.dispVar ## private control
  if (is.null(iter.mean.dispVar)) iter.mean.dispVar <- control.HLfit$max.iter.mean ## public control ## control of inner loop when some disp param is estimated
  if (is.null(iter.mean.dispVar)) iter.mean.dispVar <- 50 ## control of inner loop when some disp param is estimated  ## was 20, 06/2014
  processed$iter.mean.dispVar <- iter.mean.dispVar  
  #
  max.iter <- control.HLfit$max.iter  ## control of outer loop 
  if (is.null(max.iter)) max.iter <- 200
  processed$max.iter <- max.iter  
  #
  resid.family <- resid.model$family
  if (is.null(resid.family)) {
    resid.family <- spaMM_Gamma(log)
  } else {
    if (resid.family$family!= "Gamma") stop("resid.family must be Gamma.")## spaMM_Gamma also valid by this test
    resid.family <- spaMM_Gamma(resid.family$link) ## we will need the returned link, not the promise 
  } ## and "quoted" suitable for saving in HLfit return object
  attr(resid.family,"quoted") <- substitute(spaMM_Gamma(link),
                                            list(link=resid.family$link))
  resid.model$family <- resid.family ## resid.predictor will also be returned
  #
  HLframes <- .HLframes(formula=predictor,data=data) ## design matrix X, Y... 
  Y <- HLframes$Y
  nobs <- NROW(HLframes$X) ## not using Y which may be NULL
  if (nobs==0L) stop("No line in the data have all the variables required to fit the model.")
  if (family$family=="binomial" && NCOL(Y)>1) {
    BinomialDen <- rowSums(Y)
    if (any(BinomialDen == 0)) {
      stop("please remove missing data (i.e. for which binomial sample size is 0).")
    }
    ## It's not really possible to remove data at this stage as this may not match the dimension of the distance matrices
    ## moreover one cannot simply remove rows of a matrix "root"...
  } else {
    BinomialDen <- rep(1,nobs)
  }
  processed$BinomialDen <- BinomialDen
  y <- Y[,1L,drop=FALSE] # may work as vector, but seems marginally faster as matrix.
  if (family$family=="Gamma") {
    Gamma_min_y <- .spaMM.data$options$Gamma_min_y
    is_low_y <- (y < Gamma_min_y)
    if (any(is_low_y)) {
      #y[which(is_low_y)] <- Gamma_min_y
      warning(paste0("Found Gamma response < (Gamma_min_y=",Gamma_min_y,") . Troubles may happen."))
    }
  }
  processed$y <- y
  if (family$family=="binomial") {
    if (var(y)==0 && var(BinomialDen)==0) { warning("var(response) = 0, which may cause errors.") }  
    processed$bin_all_or_none <- all(pmin(y,BinomialDen-y)==0L)
  } else { 
    processed$bin_all_or_none <- FALSE
    if ( ! is.null(y)) { ## y may be NULL in evaluation of residProcessed
      if ( var(y)==0 ) { 
        warning("var(response) = 0, which may cause errors.") 
      } else if (var(y)<1e-3 && family$family=="gaussian") {
        warning("The variance of the response is low, which may lead to imperfect estimation of variance parameters.\n Perhaps rescale the response?")
      }
    }  
  } # (e1071::svm should fail when var response=0)
  #
  X.pv <- HLframes$X
  exp_barlist <- .process_bars(predictor,as_character=FALSE) ## but default expand =TRUE
  exp_ranef_strings <- .process_bars(barlist=exp_barlist,expand=FALSE, as_character=TRUE) ## no need to expand again
  #
  if (nrand <- length(exp_ranef_strings)) {
    ## Initialize $corr_info (ASAP to assign_cov_matrices ASAP):
    processed$corr_info <- corr_info <- new.env() ## do not set parent=emptyenv() else with(corr_info,...) will not find trivial fns such as `[`
    true_corr_types <- c("adjacency","Matern","Cauchy","AR1","corrMatrix", "IMRF")
    exp_ranef_types <- attr(exp_ranef_strings,"type") ## expanded
    corr_info$corr_types <- corr_types <- true_corr_types[match(exp_ranef_types, true_corr_types)] ## full length
    processed$control_dist <- .preprocess_control.dist(control.dist,corr_types)
    corr_families <- vector('list',length(corr_types))
    for (rd in which( ! is.na(corr_types))) corr_families[[rd]] <- do.call(corr_types[rd],list()) # corr_types' elements are function names!
    corr_info$corr_families <- corr_families
    ## need to process lcrandfamfam, a convenient object for immediate use (for the vector class more than 'lc')
    if (inherits(rand.families,"family")) rand.families <- list(rand.families) 
    if (nrand != 1L && length(rand.families)==1L) rand.families <- rep(rand.families,nrand) 
    # (fixme): this will make two copies of lcrandfamfam in processed
    rand.families <- .checkRandLinkS(rand.families)  
    processed$lcrandfamfam <- attr(rand.families,"lcrandfamfam") ## else remains NULL
    ## Assigns $corr_info$corrMatrices, $adjMatrices, $AMatrices using $corr_info$corr_type: BEFORE determining sparse precision:
    .assign_cov_matrices__from_covStruct(corr_info, covStruct=covStruct, corrMatrix=corrMatrix, adjMatrix=adjMatrix)
    Zlist <- .calc_Zlist(formula=predictor, mf=HLframes$mf, rmInt=0L, drop=TRUE, 
                         corrMats_info=corr_info$corrMatrices, barlist=exp_barlist, 
                         lcrandfamfam=processed$lcrandfamfam) 
    for (rd in seq_len(nrand)) {
      rand.families[[rd]]$prior_lam_fac <- attr(Zlist[[rd]],"prior_lam_fac")
    }
    ## Assigns $corr_info$AMatrices for IMRFs:
    .assign_AMatrices_IMRF(corr_info, Zlist, exp_barlist=exp_barlist, processed$data, control_dist=processed$control_dist) 
    ZAlist <- .calc_ZAlist(Zlist=Zlist, AMatrices=corr_info$AMatrices)
    attr(ZAlist,"exp_ranef_strings") <- exp_ranef_strings ## expanded 
    attr(ZAlist,"exp_ranef_types") <- exp_ranef_types ## expanded
    
  } else rand.families <- list() ## ## corrects the default [gaussian()] when nrand=0
  .assign_canonLink_G_LMMbool(family, processed) ## assigns processed$canonicalLink, $LMMbool, $GLMMbool using $lcrandfamfam
  #
  HLmethod <- .preprocess_HLmethod(HLmethod, family, processed$lcrandfamfam) ## not a member of the return object
  HL <- eval(parse(text=paste0("c",substr(HLmethod,3,100)))) ## extracts the (...) part into a vector
  if (length(HL)==2) HL <- c(HL,1)
  processed$HL <- HL
  ####
  ####
  ## X.pv post-processing; later extract column for fixed oefficients, but we need the offset for that
  X.pv <- .post_process_X(X.pv, HL, HLframes, control.HLfit)
  # .post_process_X() calls .rankTrim() which adds attributes "namesOri" and "rankinfo"
  #
  if (nrand) {
    # Standardize init.HLfit in the one case where it may have corrPars (simplified version of .post_process_parlist)
    if ( ! is.null(rho <- init.HLfit$rho)) {
      init.HLfit$corrPars <- list()
      if (length(adj_rd <- which(corr_types=="adjacency"))) {
        init.HLfit$corrPars[[as.character(adj_rd)]][["rho"]] <- rho
        init.HLfit$rho <- NULL
      } else stop("Invalid ambiguous 'init.HLfit' argument: single 'rho' but not single adjacency random-effect term.")
    }
    #
    ## assign sparse_precision (.determine_spprec uses $For, $corr_info, and $LMMbool:)
    if (inherits(X.pv,"sparseMatrix") && ncol(X.pv)> sum(unlist(lapply(ZAlist,ncol)))) { ## f i x m e heuristic rule
      sparse_precision <- FALSE ## presumably efficient use of Matrix::qr by .sXaug_Matrix_QRP_CHM_scaled algo
    } else { # general case
      sparse_precision <- .determine_spprec(ZAlist, processed, init.HLfit=init.HLfit) 
    }
    processed$sparsePrecisionBOOL <- sparse_precision
    #
    ## post-processing of corr_info depending on sparse_precision
    .process_corr_info_spprec(corr_types=corr_types,corr_info=corr_info, For=processed$For,sparse_precision=sparse_precision)
  } else processed$sparsePrecisionBOOL <- FALSE ## for .do_TRACE()
  processed$init_HLfit <- init.HLfit ## this can be modified by dhglm-specific code.
  ###   OFFSET
  off <- model.offset(HLframes$mf) ## look for offset from (ori)Formula 
  processed$predictor <- predictor ## save copy with offset
  if ( ! is.null(off) ) { ## offset (ori)Formula term. Check attribute (Predictor ensures offset is not both in formula and attr)
    predictor <- .stripOffset(predictor) ## note that .stripOffset()ing is not performed in .preprocess for the resid.formula 
    off <- pmax(log(.Machine$double.xmin),off) ## handling log(0) ## but if input off were NULL, output off would be is numeric(0) where it should remain NULL
  }
  attr(processed$predictor,"no_offset") <- predictor
  ## Extract columns for fixed oefficients (involves offset; ! can conflict with .rankTrim() results)
  X.Re <- X.pv ## may be updated from etaFix$beta or REMLformula, but this line makes it possible to avoid this updating when we need to avoid it
  ## reimplementation of etaFix$beta (2015/03)
  if ( length(betaFix <- etaFix$beta)>0 ) {
    namesbetafix <- names(betaFix)
    if (is.null(namesbetafix)) {
      message("The elements of etaFix$beta should be named and the names should match the column names of the design matrix.")
    }
    if (length(setdiff(namesbetafix,colnames(X.pv)))==0L) { ## if no incorrect name
      offFromEtaFix <- (X.pv[,namesbetafix,drop=FALSE] %*% betaFix)[,1] # must be vector not matrix
      namesbetavar <- setdiff(colnames(X.pv),namesbetafix)
      X.pv <- .subcol_wAttr(X.pv, j=namesbetavar, drop=FALSE)
      if (is.null(off)) {
        off <- offFromEtaFix
      } else off <- off + offFromEtaFix
      ## TRUE by default:
      if ( is.null( keepInREML <- attr(betaFix,"keepInREML") ) ||  ( ! keepInREML) ) X.Re <- X.pv ## can be overwritten according to REMLformula
    } else {
      stop("The names of elements of etaFix$beta should all match column names of the design matrix.")
    }
  } # X.Re may again be modified later
  if (is.null(off)) { ## model.frame.default(formula = locform, offset = off,...) expects a vector....
    processed$off <- rep(0,nobs) ## long form expected by spaMM_glm.fit() [as by glm.fit()] and maby by fitme_body()
  } else {
    processed$off <- off
  }
  #
  pforpv <- ncol(X.pv)
  #
  if(family$family == "binomial" && processed$bin_all_or_none && pforpv) {
    isSeparated <- is_separated(X.pv, as.numeric(y))
  }
  if (HL[1L]=="SEM") processed$SEMargs <- .preprocess_SEMargs(BinomialDen, nobs, control.HLfit, y)
  if (HL[1L]==0L) {processed$p_v_obj <-"hlik"} else processed$p_v_obj <-"p_v" ## objective for beta(_v) estim only: != outer obj 
  if (substr(HLmethod,0,2)=="ML") { # && HL[1]!="SEM") { ## FR->FR c'est bizarre d'exclure le SEM là... p_bv est il vraiment utilisé ?
    if (family$family == "binomial" && processed$bin_all_or_none && HL[1]==1L) {
      if ( ! identical(spaMM.getOption("PQL_warned"),TRUE)) {
        message("Fits using Laplace approximation may diverge for all-or-none binomial data:\n check PQL or PQL/L methods in that case.")
        spaMM.options(PQL_warned=TRUE)
      }
    }
    if ( ! is.null(REMLformula)) {
      message("Confusing combination of arguments: 'HLmethod=ML(...)' with non-null 'REMLformula'.")
      stop("  Make sure what you mean and simplify the arguments.")
    }
    if (length(predictor)==3) {
      lhs <- paste(predictor)[[2]] ## extract response, either cbind or not
    } else lhs <- "" ## occurs for predictor of phi if dhglm ML fit.
    
    # if ( nrand > 0L ) {
    #   ## build formula with only random effects  ##FR->FR  why have I done that ??
    #   exp_ranef_strings <- attr(ZAlist,"exp_ranef_strings") ## expand maybe not important here
    #   REMLformula <- as.formula(paste(lhs,"~",paste(exp_ranef_strings,collapse="+")))
    # } else 
    REMLformula <- as.formula(paste(lhs,"~ 0")) 
    attr(REMLformula,"isML") <- TRUE
  } ## else do nothing: keeps input REMLformula, which may be NULL or a non-trivial formula
  # REMLformula <- .preprocess_formula(REMLformula)
  processed$REMLformula <- REMLformula  
  if ( ! is.null(REMLformula) ) { ## differences affects only REML estimation of dispersion params, ie which p_bv is computed
    REMLFrames <- .HLframes(formula=REMLformula,data=data) ## design matrix X, Y...
    X.Re <- REMLFrames$X
    # wAugX will have lost its colnames...
    unrestricting_cols <- which(colnames(X.pv) %in% setdiff(colnames(X.pv),colnames(X.Re))) ## not in X.Re
    extra_vars <- setdiff(colnames(X.Re),colnames(X.pv)) ## example("update") tests this.
    distinct.X.ReML <- c(length(unrestricting_cols), length(extra_vars)) ## TWO booleans
    if (any(distinct.X.ReML)){
      attr(X.Re,"distinct.X.ReML") <- distinct.X.ReML 
      if (attr(X.Re,"distinct.X.ReML")[1L]) attr(X.Re,"unrestricting_cols") <- unrestricting_cols
      if (attr(X.Re,"distinct.X.ReML")[2L])attr(X.Re,"extra_vars") <- extra_vars ## example("update") tests this.
      processed$X.Re <- X.Re
    } ## else processed$X.Re  is NULL
  } ## else keep previously computed X.Re = X.pv
  # if standard ML: there is an REMLformula ~ 0 (or with ranefs ?); local X.Re and processed$X.Re is 0-col matrix
  # if standard REML: REMLformula is NULL: $X.Re is X.pv, processed$X.Re is NULL
  # non standard REML: other REMLformula: $X.Re and processed$X.Re identical, and may take essentially any value
  if (is.null(objective)) {
    if (ncol(X.Re)) { ## standard or non-standard REML
      processed$objective <- "p_bv"  ## info for fitme_body and corrHLfit_body, while HLfit instead may use return_only="p_bvAPHLs"
    } else processed$objective <- "p_v"
  } else processed$objective <- objective
  #
  models <- list(eta="",lambda="",phi="")
  if ( nrand) {
    models[["eta"]] <- "etaHGLM" 
    ##
    if ( length(corr_types[ ! is.na(corr_types)])) {
      if (For=="HLfit") {
        ranef_string <- attr(ZAlist,"exp_ranef_strings")[ ! is.na(corr_types)][1L]
        stop(paste("Term",ranef_string,"not allowed in HLfit(). Try another fitting function such as fitme()."))
      } else {
        ## Cannot be unlisted bc cf ?unlist " Non-vector elements ... are not coerced, and so a list containing ... remains a list":
        attr(ZAlist,"exp_spatial_terms") <- exp_spatial_terms <- .findSpatial(predictor,barlist=exp_barlist, nomatch=NA, expand=TRUE) ## match ZAlist for predict()
        geo_info <- vector("list",nrand) ## each element will be an environment with typical elements $distMatrix, $uniqueGeo, $nbUnique
        cov_info_mats <- vector("list",nrand)
        for (it in seq_along(corr_types)) {
          corr_type <- corr_types[it]
          if ( ! is.na(corr_type)) {
            if (corr_type== "corrMatrix") { ## could cover all specifications of a constant 'cov' Matrix without correlation parameter
              ## .check_subset_corrMatrix has the effect that the corr or prec mat used in later computation is a permutation of the inputed one
              #  according to the order of columns of ZAlist[[it]] 
              ## In the spprec case corr_info$corrMatrices[[it]] already contains a precision matrix
              cov_info_mats[[it]] <- .check_subset_corrMatrix(corrMatrix=corr_info$corrMatrices[[it]],ZA=ZAlist[[it]]) ## correlation or precision...
            } else { ## all cases where geo_info (even empty) is needed 
              geo_info[[it]] <- new.env(parent=emptyenv())
              if (corr_type == "AR1") {
                ugeo <-  attr(ZAlist[[it]],"uniqueGeo")
                if ( ! sparse_precision &&
                     ! is.null(dataordered_unique_levels <- attr(ZAlist[[it]],"dataordered_unique_levels"))) {
                  # The above test implies "if not already ! sparse when ZAlist was first evaluated"
                  # Subsetting ZAlist drops useful attributes => "uniqueGeo" must be secured in a more consistent place
                  
                  rownames(ugeo) <- apply(ugeo,1L,paste0,collapse=":")
                  geo_info[[it]]$uniqueGeo <- ugeo[(dataordered_unique_levels),,drop=FALSE]
                  #
                  ZAlist[[it]] <- ZAlist[[it]][,(dataordered_unique_levels)]
                } else {
                  for (nam in names(ugeo)) if (is.factor(fac <- ugeo[[nam]])) ugeo[[nam]] <- as.character(levels(fac))[fac]
                  geo_info[[it]]$uniqueGeo <- ugeo
                }
              # } else if (corr_type %in% c("Matern","Cauchy") ) {
                # in that case we either
                # (1) have columns in ZAlist[[it]] not reordered: no further issues. This is the current design
                #  or
                # (2) have columns in ZAlist[[it]] not reordered in .calc_Zmatrix Then
                #     the cov_info_mats[[it]] info is not yet available at preprocessing time...
                #     See further comment in .assign_geoinfo_and_LMatrices_but_ranCoefs()
              } else if (corr_type %in% c("Matern","Cauchy")) { 
                if (attr(ZAlist[[it]],"namesTerm") != "(Intercept)") {
                  geo_info[[it]]$activelevels <- activelevels <- which(colSums(ZAlist[[it]]!=0L)>0L)
                  ZAlist[[it]] <- ZAlist[[it]][,activelevels,drop=FALSE]
                } 
              } else if ( ! is.null(uniqueGeo)) { ## user-provided uniqueGeo (no example anywhere! :-) )
                if (is.list(uniqueGeo)) { ## spaMM3.0 extended syntax
                  geo_info[[it]]$uniqueGeo <- uniqueGeo[[it]] 
                } else geo_info[[it]]$uniqueGeo <- uniqueGeo #
              } else if ( ! is.null(distMatrix)) { ## user-provided distMatrix )
                if (is.list(distMatrix)) { ## spaMM3.0 extended syntax
                  geo_info[[it]]$distMatrix <- distMatrix[[it]] 
                } else geo_info[[it]]$distMatrix <- distMatrix 
              }
            }
          }
        }
        corr_info$cov_info_mats <- cov_info_mats
        processed$geo_info <- geo_info
      }
    } else { ## no true_corr_types
      if ( For %in% c("HLCor","corrHLfit")) {
        stop("Correlation model not specified in 'formula': was valid in version 1.0 but not later.")
        ## Matern or corrMatrix were allowed without a tag then
      }
    }
    # This, together with two commented lines in .is_identity(), is not clearly useful:
    #for (rd in seq_along(ZAlist)) attr(ZAlist[[rd]], "is_identity") <- .is_identity(ZAlist[[rd]], matrixcheck=TRUE)
    vec_n_u_h <- unlist(lapply(ZAlist,ncol)) ## nb cols each design matrix = nb realizations each ranef
    processed$cum_n_u_h <- cum_n_u_h <- cumsum(c(0L, vec_n_u_h)) ## if two ranef,  with q=(3,3), this is 0,3,6 ;
    vec_normIMRF <- .calc_vec_normIMRF(exp_ranef_terms=attr(ZAlist, "exp_ranef_terms"), corr_types=corr_types)   
    if (any(vec_normIMRF)) attr(ZAlist,"Zlist") <- Zlist
    processed$ZAlist <- ZAlist 
    processed$hyper_info <- .preprocess_hyper(processed=processed) # uses$ZAlist and $predictor
    #
    processed$QRmethod <- .choose_QRmethod(ZAlist, predictor) ## (fixme: rethink) typically sparse for large Matern model ?
    # QRmethod may well be dense for adjacency and then ZAfix will be dense. 
    nrd <- cum_n_u_h[nrand+1L]
    if ( ! .eval_as_mat_arg(processed)) {
      AUGI0_ZX <- list2env( list(I=.trDiagonal(n=nrd), ## avoids repeated calls to as() through rbind2...
                                 ZeroBlock= Matrix(0,nrow=nrd,ncol=pforpv), X.pv=X.pv) )
    } else {
      ## here in version up to 2.4.100 I had code trying to guess available memory by memory.limit() or /proc/meminfo
      AUGI0_ZX <- list2env( list(I=diag(nrow=nrd),ZeroBlock= matrix(0,nrow=nrd,ncol=pforpv), X.pv=X.pv) )
    } ## $ZAfix added later   and   X.pv scaled below  !!
    AUGI0_ZX$vec_normIMRF <- vec_normIMRF
    AUGI0_ZX$envir <- list2env(list(finertypes=attr(ZAlist,"exp_ranef_types"), ## to be modified later
                                    LMatrices=structure(vector("list",nrand),
                                                        is_given_by=rep("",nrand))),    
                               parent=environment(.preprocess))
    #
    if (sparse_precision) AUGI0_ZX$envir$method <- .spaMM.data$options$spprec_method  
    AUGI0_ZX <- .add_ZAfix_info(AUGI0_ZX, ZAlist, sparse_precision, corr_types, processed)
    
    # processed info for u_h inference
    processed$u_h_info <- .eval_v_h_bounds(cum_n_u_h, rand.families) 
  } else {
    models[["eta"]] <- "etaGLM" 
    AUGI0_ZX <- list(X.pv=X.pv)
  }
  if (.spaMM.data$options$X_scaling) { ## use scaling by default v.2.4.83
    ##       standard REML    ||      ML 
    if ( is.null(REMLformula) || ncol(X.Re)==0L) AUGI0_ZX$X.pv <- .scale(AUGI0_ZX$X.pv)
  }
  #if (nrand &&  ! sparse_precision ) { AUGI0_ZX$template_Xscal <- .make_Xscal(ZAL=NULL, ZAL_scaling = NULL, AUGI0_ZX=AUGI0_ZX) } # ## failure to use it efficiently
  if (nrand &&  ! sparse_precision ) { 
    if (inherits(AUGI0_ZX$ZeroBlock,"sparseMatrix")) {
      AUGI0_ZX$Zero_sparseX <- rbind2(AUGI0_ZX$ZeroBlock, as(AUGI0_ZX$X.pv,"CsparseMatrix"))
    } else AUGI0_ZX$Zero_sparseX <- rbind2(AUGI0_ZX$ZeroBlock, AUGI0_ZX$X.pv)
  }
  processed$AUGI0_ZX <- AUGI0_ZX
  #
  processed$coef12needed <- ! ((family$family=="gaussian" && family$link=="identity")
                              || (family$family=="Gamma" && family$link=="log")  ) ## two ad hoc cases
  #
  if (is.null(subs_p_weights <- substitute(prior.weights))) {
    prior.weights <- structure(rep(1L,nobs),unique=TRUE) ## <- 1L prevented by glm -> model.frame(... prior.weights)
  } else if ( ! (inherits(subs_p_weights,"call") && subs_p_weights[[1L]] == "quote") )  {
    prior.weights <- as.vector(stats::model.weights(validdata)) ## as.vector as in say lm() protects against array1d
    if ( ! is.numeric(prior.weights)) 
      stop("'weights' must be a numeric vector")
    if (any(prior.weights < 0)) 
      stop("negative weights not allowed")
    attr(prior.weights,"unique") <- (length(unique(prior.weights))==1L) 
    #attr(prior.weights,"only1") <- all(upw==1L)
  } else attr(prior.weights,"unique") <- FALSE ## when 'prior.weights' is a quoted expression   
  processed$prior.weights <- prior.weights
  #
  ## algorithms (control of defaults remains in the HLfit code)
  betaFirst <- control.HLfit$betaFirst ##  
  if (is.null(betaFirst)) {
    betaFirst <- FALSE
  } else if (betaFirst && HL[1]=="SEM") {
    message("betaFirst && HLmethod= SEM: betaFirst turned to FALSE")
    betaFirst <- FALSE
  }
  processed$betaFirst <- betaFirst ##
  LevenbergM <- .reformat_LevenbergM(control.HLfit$LevenbergM)
  if (is.na(LevenbergM["default_LM_start"])) { ## ie, if not set globally by spaMM.options()
    if (family$family == "binomial" && processed$bin_all_or_none ) { # both PQL and ML
      if (HL[1L]==0L) {
        LevenbergM["default_LM_start"] <- TRUE ## PQL/L + LevenbergM combine safely and relatively fast.
      } else  LevenbergM["default_LM_start"] <- TRUE ## main reason for TRUE is to reduce variance in fit time: cf BINARYboot test
    } else LevenbergM["default_LM_start"] <- FALSE
  }
  ## prefit less useful since Levenberg works better
  if (is.na(LevenbergM["default_prefit"])) { ## ie, if not set globally by spaMM.options()
    LevenbergM["default_prefit"] <- FALSE ## and should not be looked for anyway
  }
  processed$LevenbergM <- LevenbergM
  #
  if (nrand) {   
    processed$lambda.Fix <- .reformat_lFix(.getPar(ranFix,"lambda"), nrand)
    models[["lambda"]] <- rep("lamScal",nrand) ## even for adjacnency, random slope...
    ################################################################################
    # for a random slope term, ie v= v_1+x v_2 , the x went into the general ZAL matrix 
    # (construction of Zlist by .calc_Zlist(), and
    # we are still estimating the lambda's using a X_lamres with 0/1's only
    # unless there is a non-trivial model for the lambdas
    ################################################################################
    if (all(models[["lambda"]]=="lamScal")) { ## all mixed models handled in 06/2015 (random slope, adjacency...) hence currently always TRUE
      Xi_cols <- attr(ZAlist,"Xi_cols")
      if (any(Xi_cols>1 & ! processed$lcrandfamfam=="gaussian")) {
        stop("(!) random slope models implying correlated non-gaussian random effects are not fitted.")
      }
      cum_Xi_cols <- cumsum(c(0, Xi_cols)) ## if two ranef,  with n_u_h=(3,3), this is 0,3,6. cum_h_u_h[nrand+1] is then 6, the total # of realizations
      n_u_h <- rep(0, sum(Xi_cols))
      #for (i in 1:nrand) n_u_h[(cum_Xi_cols[i]+1L):cum_Xi_cols[i+1L]] <- ncol(ZAlist[[i]]) ##  nlevels(Subject[[i]])
      # if 18 responses in a random slope model ncol(ZAlist[[i]]) is 36 while nlevels(Subject[[i]]) was 18
      for (i in 1:nrand) n_u_h[cum_Xi_cols[i]+seq(Xi_cols[i])] <- ncol(ZAlist[[i]])/Xi_cols[i]
      # h_u_h not n_u_h ...
      cum_h_u_h <- cumsum(c(0, n_u_h)) ## if two "Intercept" ranefs,  with n_u_h=(3,3), this is 0,3,6. cum_h_u_h[nrand+1] is then 6, the total # of realizations
        ## if (1+X|...) +(1|...),  with n_u_h=(3,4), this is 0,3,6,10. cum_h_u_h[sum(Xi_cols)+1] is then 10, the total # of realizations
      X_lamres <- matrix(0,cum_h_u_h[sum(Xi_cols)+1L],sum(Xi_cols))
      colnames(X_lamres) <- unlist(attr(ZAlist,"namesTerms"))
      for (i in seq(nrand)) {
        for (j in (cum_Xi_cols[i]+1L):cum_Xi_cols[i+1L]) {
          X_lamres[(cum_h_u_h[j]+1L):cum_h_u_h[j+1L],j] <- 1L ## this maps the deviance residuals to the lambda's to be estimated from them. None of the random-slope columns is a constant full column because each dev res is used for estimating only one lambda. Nevertheless, the first col will be called "(Intercept)", and this makes a valid output.
        }
      } 
    } else {  ## linear predictor for variance of random effects (lambda) (lamGLM or lamHGLM) 
      if (any(models[["lambda"]]=="lamHGLM")) { ##need distinct calls... to fit each lambda model  
        if (length(formulaLambda)==2) formulaLambda <- as.formula(paste('"lambda"',paste(formulaLambda,collapse=" ")))
        if (!is.null(.parseBars(formulaLambda))) {  ## lamHGLM
          models[["lambda"]] <- list("lamHGLM")
        } else models[["lambda"]] <- list("lamGLM")  
        colnames(X_lambda) <- colnames(fr_lambda$X) ## but code not effective, fr_lambda not computed
      } else { ## can use a single design matrix for all random effects, which is convenient.
        stop("LIKELY missing code to handle linear predictor for lambda.")
        # la suite c'est dexu residus de code a assembler: il faut une liste de HLframes du type
        fr_lambda <- .HLframes(formula=formulaLambda,data=data) ## but the "data" should probably be distinct data here, with nrow=number of reals of ranefs 
        # (pobablement calculee en amont pour determiner lamScal aussi...) ensuite extraire les design matrices
        #X_lamres ? Xi_cols ?
      }
    } 
    processed$X_lamres <- X_lamres ## for glm for lambda, and SEMbetalambda
    ranCoefs <- .getPar(ranFix,"ranCoefs") ## may be NULL
    processed$ranCoefs_blob <- .process_ranCoefs(processed, ranCoefs, use_tri=TRUE) 
    processed$AUGI0_ZX$envir$finertypes[processed$ranCoefs_blob$isRandomSlope] <- "ranCoefs" ## (*creates* the variable in the *environment* so that .evalWrapper() finds it)
  } #else processed$ranCoefs_blob <- NULL
  #
  phi.Fix <- .getPar(ranFix,"phi")
  if (is.null(phi.Fix)) {
    if (family$family %in% c("poisson","binomial","COMPoisson","negbin")) {
      phi.Fix <- 1 
    } # else if (var(y)==0) phi.Fix <- .spaMM.data$options$min_disp
  } else if (any(phi.Fix==0)) stop("phi cannot be fixed to 0.")
  processed$phi.Fix <- phi.Fix
  #
  resid.formula <- resid.model$formula # .reformat_resid_model() has called .preprocess_formula()  but we have to call it again after pasting ".phi"
  if ( is.null(phi.Fix)) {
    if ( ! is.null(.parseBars(resid.formula))) {
      if (is.null(resid.model$rand.family)) resid.model$rand.family <- gaussian() # avoids rand.families being NULL in call below.
      preprocess_arglist <- list(control.HLfit=control.HLfit, ## constrained
                         ranFix=resid.model$fixed, 
                         HLmethod=HLmethod, ## constrained
                         predictor=resid.formula, ## obvious
                         resid.model=resid.model$resid.model, # potentially allows nested resid.model's... 
                         REMLformula=NULL, # constrained
                         data=data, # obvious (?) 
                         family=resid.family, # obvious
                         BinomialDen=NULL, # obviously no binomial response
                         rand.families=resid.model$rand.family, # (NULL not handled by preprocess); 
                         #   outer preprocess calls *receive* a default value from formals(HLfit)
                         etaFix=resid.model$etaFix, ## not constrained, but should rather use 'resid.model$fixed'
                         prior.weights=NULL, ## currently defined  dynamically using lev_phi...
                         control.glm=control.glm, ## constrained
                         verbose=NULL, ## TRACE would be overriden by the final do_TRACE call of the parent .preprocess()
                         For="fitme", ## constrained: preprocess must allow spatial and non-spatial models
                         init.HLfit=as.list(resid.model$init.HLfit) ## converts NULL to list() as exp'd by .preprocess()
      )
      ## preprocess formal arguments that were ignored up to v.2.4.30 14/05/2018:
      other_preprocess_args <- setdiff(names(formals(.preprocess)),names(preprocess_arglist))
      preprocess_arglist[other_preprocess_args] <- resid.model[other_preprocess_args]
      residProcessed <- do.call(.preprocess,preprocess_arglist) ## cf verbose explicitly set to NULL 
      # preprocess here plays the role of fitme as wrapper bringing the following info to fitme_body:
      #
      # we add ".phi" to attr(residProcessed$predictor - for summary() only ? But then same operation on version with hyper-ranefs
      fullform <-  .preprocess_formula(as.formula(paste(".phi",.DEPARSE(residProcessed$predictor))))
      mostattributes(fullform) <- attributes(residProcessed$predictor)
      if ( ! is.null(hy_form <- attr(fullform,"hyper_info")$formula)) attr(fullform,"hyper_info")$formula <- as.formula(paste(".phi",.DEPARSE(hy_form)))
      residProcessed$predictor <- fullform
      # residProcessed$y <- residProcessed$HLframes$Y <- "removed for safety" # must be NULL
      if (identical(names(resid.model$fixed$phi),"default")) message("'phi' of residual dispersion model set to 1 by default")
      processed$residProcessed <- residProcessed
      models[["phi"]] <- "phiHGLM" 
      p_phi <- NA
    } else {
      residFrames <- .HLframes(formula=resid.formula,data=data)
      attr(resid.formula,"off") <- model.offset(residFrames$mf) ## only for summary.HLfit()
      ## if formula= ~1 and data is an environment, there is no info about nobs, => fr_disp$X has zero rows, which is a problem later 
      p_phi <- NCOL(residFrames$X)
      namesX_disp <- colnames(residFrames$X)
      if (p_phi==1 && namesX_disp[1]=="(Intercept)"
          && is.null(attr(resid.formula,"off")) ## added 06/2016 (bc phiScal does not handle offset in a phi formula) 
      ) {
        models[["phi"]] <- "phiScal"
      } else { 
        models[["phi"]] <- "phiGLM"
      }
      resid.model$formula <- resid.formula  ## absent  if phiHGLM has been detected
    } 
  } else p_phi <- 0
  processed$models <- models
  processed$p_fixef_phi <- p_phi # no X_disp is saved in processed
  if ( family$family %in% c("binomial","poisson","COMPoisson","negbin")) {
    ## the response variable should always be Counts
    if (max(abs(y-as.integer(y)))>1e-05) {
      stop("response variable should be integral values.")
    }
    if ( .DEPARSE(resid.formula) != "~1") {
      warning(paste0("resid.model is ignored in ",family$family,"-response models"))
    }
  }
  #
  processed$augZXy_cond <- .determine_augZXy(processed,control.HLfit=control.HLfit)
  if (processed$augZXy_cond) {
    processed$augZXy_env <- new.env(parent=emptyenv()) 
    processed$augZXy_env$objective <- -Inf 
  }
  #
  if (processed$LMMbool) {
    ## identifiability checks cf modular.R -> checkNlevels() in lmer:
    vec_n_u_h <- diff(processed$cum_n_u_h)
    if (any(vec_n_u_h<2L) && is.null(phi.Fix)) {
      problems <- which(vec_n_u_h<2L) 
      for (iMat in problems) {
        mess <- paste0("Only ",vec_n_u_h[iMat]," level for random effect ",
                      attr(ZAlist,"exp_ranef_strings")[iMat],
                      ";\n   this model cannot be fitted unless phi is fixed.")
        warning(mess)
      }
    }
    if (any(vec_n_u_h==nobs) && models[["phi"]] %in% c("phiScal","phiGLM")) { ## tests (mean)ranefs if intercept in resid formula
      resid.mf <- residFrames$mf 
      if (attr(attr(resid.mf, "terms"),"intercept")!=0L) { ## there is an intercept in the resid.model formula
        # ideally for random-coefficients models we should compare the design columns... 
        ## FR->FR cf isNested check as in https://github.com/lme4/lme4/blob/master/R/utilities.R, 
        problems <- which(vec_n_u_h==nobs) 
        for (iMat in problems) {
          term_ranef <- attr(ZAlist,"exp_ranef_strings")[iMat]
          if (# is.null(LMatrix) && ## does not seem useful
            substr(term_ranef, 1, 1)=="(" ## excludes spatial ranefs 
            && ! is.numeric(ranFix$lambda[iMat])
          ) {
            mess <- paste0("Number of levels = number of observations for random effect ", term_ranef,
                          ";\n   this model cannot be fitted unless phi is fixed, or the variance",
                          "\n   of this effect is fixed, or a non-trivial correlation matrix is given.")
            stop(mess)
          }          
        }
      }
    }
  } 
  #
  processed$HLframes <- HLframes[c("Y","fixef_off_terms","fixef_levels")] ## Y for family$initialize() (binomial spec.), others for predict()
  #processed$HLframes$all_terms <- attr(HLframes$mf, "terms") ## also for predict, with poly(.,raw=FALSE) term
  processed$residModel <- resid.model 
  processed$rand.families <- rand.families
  processed[["control.glm"]] <- do.call("glm.control", control.glm) ## added 04/2016 (LHS <- RHS list)
  processed$port_env <- new.env(parent=emptyenv()) ## added 09/2016
  processed$verbose <- .reformat_verbose(verbose,For=For) ## added 09/2017*
  .do_TRACE(processed)
  class(processed) <- c("arglist",class(processed))
  spaMM.options(COMP_maxn_warned=FALSE,COMP_geom_approx_warned=FALSE)
  return(processed)
}

.eval.update.call <- function(mc,...) { # not currently used
  mc <- as.list(mc)
  dotlist <- list(...)
  mc[names(dotlist)] <- dotlist ## a un moment j'ai mis cette ligne en commentaire, ce qui rend la fonction ineffective !
  eval(as.call(mc))  
}