## 'has-no-ranef' can be tested by is.null(.parseBars(re.form))
# devel version of lme4 has noReForm() that assumes that there in no fixed effect in re.form, with result
# TRUE if argument is NA or ~0, 
# FALSE if argument is NULL or a non-trivial formula (re.form=NULL in particular means that prediction assumes the original random-effect terms).
# We additionally want TRUE if argument is ~<only fixef> (fixed effect will be ignored anyway)
# The following is equivalent to noReForm() except for the last test, is.null(.parseBars(re.form)) 
# ~0 returns TRUE, ~Days returns TRUE [as the last test is TRUE in both cases], ~<including ranefs> returns FALSE 
.noRanef <- function(re.form) { # Stupidly returns FALSE for explicit formula with LHS. So re.form better not have a LHS.
  (!is.null(re.form) && !inherits(re.form,"formula") && is.na(re.form)) ||
    (inherits(re.form,"formula") && length(re.form)==2 && is.null(.parseBars(re.form)))
}

.newetaFix <- function(object, newMeanFrames,validnames=NULL) { # newMeanFrames must have cols for the fixed beta's, absent from the object's $X.pv 
  ## newdata -> offset must be recomputed. 
  ## dans l'état actuel $fixef et complet, incluant les etaFix$beta: pas besoin de les separer
  ## mais il peut contenir des NA ! à enlever
  # le newMeanFrames$X contient a priori les cols des etaFix$beta, contrairement à object$X.pv => don't use the latter cols 
  if (is.null(validnames)) {
    est_and_fix <- names(which(!is.na(object$fixef))) ## estimated + etaFix$beta
    validnames <- intersect(colnames(newMeanFrames$X) ,est_and_fix) # would be colnames(newMeanFrames$X) when there is no NA if object$fixef. would validnames=est_and_fix suffice ?
  }
  if (length(validnames)) {
    etaFix <-  drop(newMeanFrames$X[,validnames,drop=FALSE] %*% object$fixef[validnames]) ## valid even if ncol(newMeanFrames$X) = 0
  } else etaFix <- rep(0, nrow(newMeanFrames$X))
  off <- model.offset( newMeanFrames$mf) ### look for offset from (ori)Formula 
  if ( ! is.null(off)) etaFix <- etaFix + off   
  return(etaFix)
}

.r_resid_var <- function(mu,phiW,sizes,family,
                         COMP_nu,NB_shape,zero_truncated, famfam, nsim=1L) { 
  # we cannot use family()$simulate bc it assumes a fit object as input
  resu <- switch(famfam,
                    gaussian = rnorm(nsim*length(mu),mean=mu,sd=sqrt(phiW)),
                    poisson = .rpois(nsim*length(mu),mu,zero_truncated=zero_truncated), 
                    binomial = rbinom(nsim*length(mu),size=sizes,prob=mu),
                    Gamma = {
                      y <- rgamma(nsim*length(mu), shape= 1 / phiW, scale=mu*phiW) # mean=sh*sc=mu, var=sh*sc^2 = mu^2 phiW
                      Gamma_min_y <- .spaMM.data$options$Gamma_min_y
                      is_low_y <- (y < Gamma_min_y)
                      if (any(is_low_y)) y[which(is_low_y)] <- Gamma_min_y 
                      y
                    }, ## ie shape increase with prior weights, consistent with Gamma()$simulate / spaMM_Gamma()$simulate
                    COMPoisson = {
                      lambdas <- attr(mu,"lambda") # F I X M E an environment would keep values ?
                      if (is.null(lambdas)) {
                        sapply(mu, function(muv) {
                          lambda <- family$mu2lambda(muv)
                          .COMP_simulate(lambda=lambda,nu=COMP_nu, nsim=nsim)
                        })
                      } else sapply(lambdas,.COMP_simulate,nu=COMP_nu, nsim=nsim)
                    },
                    negbin = .rnbinom(nsim*length(mu), size=NB_shape, mu_str=mu, zero_truncated=zero_truncated), 
                    stop("(!) random sample from given family not yet implemented")
  )
  if (nsim>1L) dim(resu) <- c(length(mu),nsim)
  resu
} ## vector-valued function from vector input

.r_resid_var_over_cols <- function(needed, 
                                   mu,  # a list over nsim !
                                   phiW, sizes, family, families, resp_range=NULL, is_mu_fix_btwn_sims,
                                   is_phiW_fix_btwn_sims=attr(phiW,"is_phiW_fix_btwn_sims"), nsim=1L, as_matrix=FALSE, mv_it=NULL) {
  if ( ! is.null(families)) {
    cum_nobs <- attr(families,"cum_nobs")
    rowS <- vector("list", length(families))
    for (mv_it in seq_along(families)) {
      resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
      rowS[[mv_it]] <- .r_resid_var_over_cols(needed, mu, phiW=phiW[resp_range,,drop=FALSE], sizes=sizes[resp_range], 
                                              family=families[[mv_it]], families=NULL, resp_range=resp_range,
                                              is_mu_fix_btwn_sims=is_mu_fix_btwn_sims,
                                              is_phiW_fix_btwn_sims=is_phiW_fix_btwn_sims[mv_it], nsim=nsim, as_matrix=TRUE, mv_it=mv_it)
    }
    return(do.call(rbind, rowS))
  }
  block <- NA*phiW
  famfam <- family$family
  if (is_mu_fix_btwn_sims && is_phiW_fix_btwn_sims) { # 2nd condition should be trivially true for binomial, poisson, COMPoisson (phi=1, no prior.weights)
    if (is.null(resp_range)) {
      mu_all <- mu[[1L]] 
    } else mu_all <- structure(mu[[1L]][resp_range], p0=attr(mu[[1L]],"p0")[[mv_it]], mu_U=attr(mu[[1L]],"mu_U")[[mv_it]])
    # : subsetting otherwise loses crucial p0 and mu_U attributes for Tnegbin
    # These attributes also exists for Tpoisson
    block <- .r_resid_var(mu_all, phiW=phiW[,1L],sizes=sizes, COMP_nu=environment(family$aic)$nu,
                          NB_shape=environment(family$aic)$shape, 
                          zero_truncated=identical(family$zero_truncated,TRUE), 
                          famfam=famfam, family=family, nsim=nsim) # vector or nsim-col matrix
    if (as_matrix && is.null(dim(block))) dim(block) <- c(length(block),1L) # forces matrix for mv code using rbind()
  } else {
    # get these variables once for all cols:
    if (famfam =="COMPoisson") {
      COMP_nu <- environment(family$aic)$nu
    } else if (famfam == "negbin") {
      NB_shape <- environment(family$aic)$shape
      zero_truncated <- identical(family$zero_truncated,TRUE) 
    } else if (famfam == "poisson") {
      zero_truncated <- identical(family$zero_truncated,TRUE) 
    }
    for (sim_it in seq_len(needed)) {
      if (is.null(resp_range)) {
        mu_it <- mu[[sim_it]]
      } else mu_it <- structure(mu[[sim_it]][resp_range], p0=attr(mu[[sim_it]],"p0")[[mv_it]], mu_U=attr(mu[[sim_it]],"mu_U")[[mv_it]])
      block[ ,sim_it] <- .r_resid_var(mu_it, phiW=phiW[ ,sim_it], sizes=sizes, COMP_nu=COMP_nu, NB_shape=NB_shape, 
                                  zero_truncated=zero_truncated, famfam=famfam, family=family, nsim=1L)
    }
  }
  return(block)
}

.calc_ZAlist_newdata <- function(object, newdata) {
  # we simulate with all ranefs (treated conditionnally|ranef or marginally) hence 
  # * we need design matrices for all ranefs
  # * we need values of all the original variables
  # hence we use an "## effective '.noFixef'" : formula with only ranefs of the fit. But 1st version fails for mv; second may be more straightforward anyway
  #old_ranef_form <- formula.HLfit(object, which="")[-2L] 
  #old_ranef_form <- as.formula(paste("~",(paste(.process_bars(old_ranef_form),collapse="+")))) ## 
  #frame_ranefs <- .calc_newFrames_ranef(formula=old_ranef_form,data=newdata, fitobject=object)$mf 
  #exp_ranef_terms <- .process_bars(old_ranef_form[[2L]], # barlist must remain missing here
  #                                 expand=TRUE, which. = "exp_ranef_terms")
  if (is.null(vec_nobs <- object$vec_nobs)) {
    old_ranef_form <- as.formula(paste("~",(paste(attr(object$ZAlist,"exp_ranef_strings"),collapse="+")))) 
    exp_ranef_terms <- attr(object$ZAlist, "exp_ranef_terms")
    Zlist <- .calc_Zlist(exp_ranef_terms=exp_ranef_terms, data=newdata, rmInt=0L, drop=TRUE,sparse_precision=FALSE,
                         corrMats_info=object$strucList,
                         lcrandfamfam=attr(object$rand.families,"lcrandfamfam"))
    amatrices <- .get_new_AMatrices(object,newdata=newdata) # .calc_newFrames_ranef(formula=old_ranef_form,data=newdata, fitobject=object)$mf)
    newZAlist <- .calc_normalized_ZAlist(Zlist=Zlist,
                                         AMatrices=amatrices,
                                         vec_normIMRF=object$ranef_info$vec_normIMRF, 
                                         strucList=object$strucList)
  } else { 
    map_rd_mv <- attr(object$ZAlist, "map_rd_mv")
    ori_exp_ranef_terms <- attr(object$ZAlist, "exp_ranef_terms")
    ori_exp_ranef_strings <- attr(object$ZAlist,"exp_ranef_strings")
    for (mv_it in seq_along(vec_nobs)) {
      rd_in_mv <- map_rd_mv[[mv_it]]
      exp_ranef_strings_it <- ori_exp_ranef_strings[rd_in_mv]
      old_ranef_form <- as.formula(paste("~",(paste(exp_ranef_strings_it,collapse="+")))) 
      exp_ranef_terms_it <- structure(ori_exp_ranef_terms[rd_in_mv], type=attr(ori_exp_ranef_terms,"type")[rd_in_mv])
      Zlist <- .calc_Zlist(exp_ranef_terms=exp_ranef_terms_it, data=newdata, rmInt=0L, drop=TRUE,sparse_precision=FALSE,
                           corrMats_info=object$strucList[rd_in_mv],
                           lcrandfamfam=attr(object$rand.families,"lcrandfamfam")[rd_in_mv])
      amatrices <- .get_new_AMatrices(object,newdata=newdata) # .calc_newFrames_ranef(formula=old_ranef_form,data=newdata, fitobject=object)$mf)[rd_in_mv]
      ZAlist_it <- .calc_normalized_ZAlist(Zlist=Zlist,
                                           AMatrices=amatrices,
                                           vec_normIMRF=object$ranef_info$vec_normIMRF,
                                           strucList=object$strucList[rd_in_mv])
      names(ZAlist_it) <- rd_in_mv
      attr(ZAlist_it,"exp_ranef_strings") <- exp_ranef_strings_it
      if (mv_it>1L) {
        newZAlist <- .merge_ZAlists(newZAlist, ZAlist_it, nobs1=nrow(newdata), nobs2=nrow(newdata), mv_it)
      } else newZAlist <- ZAlist_it
    }  
  }
  return(newZAlist)
}


# simulate.HLfit(fullm[[2]],newdata=fullm[[1]]$data,size=fullm[[1]]$data$total) for multinomial avec binomial nichées de dimension différentes
# FR->FR misses the computation of random effects for new spatial positions: cf comments in the code below
simulate.HLfit <- function(object, nsim = 1, seed = NULL, newdata=NULL,
                           type = "marginal", re.form, conditional=NULL, 
                           verbose=c(type=TRUE, showpbar= eval(spaMM.getOption("barstyle"))),
                           sizes=NULL , resp_testfn=NULL, phi_type="predict", prior.weights=object$prior.weights, 
                           variances=list(), ...) { ## object must have class HLfit; corr pars are not used, but the ZAL matrix is.
  ## RNG stuff copied from simulate.lm
  was_invColdoldList_NULL <- is.null(object$envir$invColdoldList) # to be able to restore initial state 
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else { ## this makes changes to RNG local where 'seed' is used:
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  if (inherits(object,"HLfitlist")) { 
    message("simulate does not yet work on list of fits as returned by multinomial fit:")
    message(" run simulate on each of the individual fits in the list")
    stop() ## FR->FR also some basic changes in fixedLRT but more would be needed 
  }  
  if ( ! is.null(conditional)) {
    warning("argument 'conditional' is obsolete and will be deprecated. Use 'type' instead.")
    if (conditional) {type <- "residual"} else type <- "marginal"
  }
  if (is.na(verbose["showpbar"])) { # e.g. verbose =TRUE or verbose=c(type=TRUE)
    if (is.na(verbose["type"])) verbose["type"] <- verbose # need at least a boolean argument here
    verbose["showpbar"] <- eval(.spaMM.data$options$barstyle)
  } else if (is.na(verbose["type"])) verbose["type"] <- TRUE # user set verbose=c(showpbar=.) but not type
  if (type=="predVar") {
    pred_type <- "predVar_s.lato" 
    if ( ! length(variances)) stop("A 'variances' argument must be specified (e.g., variances=list(predVar=TRUE))")
    variances$cov <- TRUE # mandatory, overriding any user's variances$cov argument
  } else if (type=="(ranef|response)") {
    pred_type <- "predVar_s.lato" 
    variances <- list(linPred=TRUE, disp=FALSE, cancel_X.pv=TRUE, cov=TRUE) # mandatory, overriding any user's variance argument
  } else pred_type <- ""
  if (is.null(sizes)) sizes <- .get_BinomialDen(object)
  nrand <- length(object$ZAlist)
  if (nrand>0L) {
    if ( missing(re.form)) {
      if (type=="marginal") {
        re.form <- NA
      } else if (type=="residual") re.form <- NULL
      # type "predVar" leaves 're.form' missing. as it is not used (which should be equivalent to re.form=NULL)
    } else if (inherits(re.form,"formula")) {
      if (pred_type=="predVar_s.lato") warning("Non-default 're.form' is *currently* ignored when type='",type,"'.")
      re.form <- .preprocess_formula(re.form)
      ori_exp_ranef_strings <- attr(object$ZAlist,"exp_ranef_strings")
      new_exp_ranef_strings <- .process_bars(re.form,expand=TRUE)
      newinold <- unlist(sapply(lapply(new_exp_ranef_strings, `==`, y= ori_exp_ranef_strings), which)) ## unlist() bc empty list() can otherwise occur  
    } else if (anyNA(re.form)) { # anyNA rather than is.na, to handle NULL
      if (pred_type=="predVar_s.lato") warning("Non-default 're.form' is *currently* ignored when type='",type,"'.")
    }
  }
  resu <- NULL
  done <- 0L
  verbtype <- verbose[["type"]]
  is_mu_fix_btwn_sims <- FALSE
  while((needed <- nsim-done)) { ## loop operates only for resp_testfn
    if (nrand==0L) { ## note that replicate mu's can still be variable for non-standard pred_type
      if (pred_type=="predVar_s.lato") { ## re.form ignored so de facto NULL
        if (type=="(ranef|response)") {
          stop("meaningless argument type='(ranef|response)' for a fixed-effect model")
        } else if (verbtype) cat("Simulation from linear predictor variance | observed response:\n") 
        variances$cov <- (NROW(newdata)!=1L)
        point_pred_eta <- predict(object,newdata=newdata, type="link", control=list(fix_predVar=NA),
                                  variances=variances, verbose=verbose, ...) 
        predVar <- attr(point_pred_eta,"predVar")
        if (is.null(predVar)) stop("A 'variances' argument should be provided so that prediction variances are computed.") 
        rand_eta <- mvrnorm(n=needed,mu=point_pred_eta[,1L], predVar)
        if (needed>1L) rand_eta <- t(rand_eta) ## else mvrnorn value is a vector
        mu <- .fv_linkinv(eta=rand_eta, family=object$family, families=object$families)
      } else {
        mu <- predict(object,newdata=newdata,binding=NA,control=list(fix_predVar=FALSE), verbose=verbose)
        is_mu_fix_btwn_sims <- TRUE
        mu <- replicate(needed,mu,simplify=FALSE) # always a list at this stage
      }
    } else { ## MIXED MODEL
      if (pred_type=="predVar_s.lato") { ## re.form ignored so de facto NULL
        if (verbtype) {
          if (type=="(ranef|response)") {
            cat("Simulation from random-effects variance | observed response:\n")
          } else cat("Simulation from linear predictor variance | observed response:\n") 
        } 
        if (all(attr(object$rand.families,"lcrandfamfam")=="gaussian")){
          variances$cov <- (NROW(newdata)!=1L) # simulate() always need a covmatrix but for a signle response it is trivial
          # detect all cases where the cov mat is large:
          if (is.null(variances$as_tcrossfac_list)) variances$as_tcrossfac_list <- ( (is.null(newdata) && length(object$y)>200L) ||
                                                                                      NROW(newdata)>200L)
          point_pred_eta <- predict(object,newdata=newdata, type="link", control=list(fix_predVar=NA),
                                    variances=variances, verbose=verbose, ...) 
          predVar <- attr(point_pred_eta,"predVar")
          if (is.null(predVar)) stop("A 'variances' argument should be provided so that prediction variances are computed.") 
          if (is.list(predVar)) {
            rand_eta <- vector('list',length(predVar))
            for (it in seq_len(length(predVar))) {
              if (it==1L) {mu <- point_pred_eta[,1L]} else {mu <- rep(0,length(point_pred_eta[,1L]))}
              if ( ! is.null(predVar[[it]])) rand_eta[[it]] <- .mvrnorm(n=needed,mu=mu, tcross_Sigma = predVar[[it]])
              # else predVar[[it]] remains NULL [cf IMRF terms for newdata] and lengths() is used to remove them:
            }
            rand_eta <- Reduce("+",rand_eta[lengths(rand_eta)>0L])
          } else rand_eta <- mvrnorm(n=needed,mu=point_pred_eta[,1L], predVar)
          if (needed>1L) rand_eta <- t(rand_eta) ## else mvrnorn value is a vector
          mu <- .fv_linkinv(eta=rand_eta, family=object$family, families=object$families) ## ! freqs for binomial, counts for poisson: suitable for final code
        } else stop("This conditional simulation is not implemented for non-gaussian random-effects")
      } else if ( is.null(re.form)) { ## conditional on (all) predicted ranefs, type = "residual"
        if (verbtype) cat("simulation of residuals, conditional on point predictions (hence on random effects):\n") 
        mu <- predict(object,newdata=newdata,binding=NA,control=list(fix_predVar=FALSE), verbose=verbose)
        is_mu_fix_btwn_sims <- TRUE
        mu <- replicate(needed,mu,simplify=FALSE) #matrix(rep(mu,nsim),ncol=nsim)
      } else if ( inherits(re.form,"formula") || is.na(re.form) ){ ## explicit re.form; or unconditional MIXED MODEL, type= "marginal"
        if (verbtype) {
          if (inherits(re.form,"formula")) {
            cat("Simulation conditional on random effect(s) retained in 're.form':\n")
          } else cat("Unconditional simulation:\n") 
        }
        # mu_fixed <- predict(object, newdata=newdata, re.form=re.form,binding=NA,control=list(fix_predVar=FALSE), verbose=verbose) ## mu_T
        # eta_fixed <- .eta_linkfun(mu_fixed, object$family, object$families)
        eta_fixed <- predict(object, newdata=newdata, type="link",
                             re.form=re.form,binding=NA,control=list(fix_predVar=FALSE), verbose=verbose)
        if (is.null(newdata)) { ## we simulate with all ranefs (treated conditionnally|ranef or marginally) hence no selection of matrix
          ZAL <- get_ZALMatrix(object)
          cum_n_u_h <- attr(object$lambda,"cum_n_u_h")
          vec_n_u_h <- diff(cum_n_u_h)
        } else {
          newZAlist <-  .calc_ZAlist_newdata(object, newdata)
          ZALlist <- .compute_ZAXlist(object$strucList,newZAlist)
          ZAL <- .ad_hoc_cbind(ZALlist, as_matrix=FALSE ) 
          vec_n_u_h <- unlist(lapply(ZALlist,ncol)) ## nb cols each design matrix = nb realizations each ranef
          cum_n_u_h <- cumsum(c(0,vec_n_u_h))
        }
        lcrandfamfam <- attr(object$rand.families,"lcrandfamfam") ## unlist(lapply(object$rand.families, function(rf) {tolower(rf$family)}))
        if (inherits(re.form,"formula")) lcrandfamfam[newinold] <- "conditional" ## if is.na(re.form), lcrandfamfam is unchanged
        fittedLambda <- object$lambda.object$lambda_est
        locfn <- function(it) {
          nr <- vec_n_u_h[it]
          u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
          if ( ! is.null(object$rand.families[[it]]$prior_lam_fac) && ! is.null(newdata)) { # prior_lam_fac is the 'design' for non-ranCoef (wei-1|.)
            leftOfBar_terms <- attr(object$ZAlist,"exp_ranef_terms")[[it]][[2L]]
            leftOfBar_mf <- model.frame(as.formula(paste("~",leftOfBar_terms)), newdata, xlev = NULL) 
            prior_lam_fac <- leftOfBar_mf[,1L]^2 ## assumes simple syntax (wei-1|.)
            loclambda <- object$lambda.object$lambda_list[[it]]* prior_lam_fac
          } else loclambda <- fittedLambda[u.range] ## includes prior_lam_fac
          newU <- replicate(needed,{
            switch(lcrandfamfam[it], ## remainder of code should be OK for rand.families
                   gaussian = rnorm(nr,sd=sqrt(loclambda)),
                   gamma = rgamma(nr,shape=1/loclambda,scale=loclambda),
                   beta = rbeta(nr,1/(2*loclambda),1/(2*loclambda)),
                   "inverse.gamma" = 1/rgamma(nr,shape=1+1/loclambda,scale=loclambda), ## yields inverse gamma (1+1/object$lambda,1/object$lambda)
                   conditional= rep(0,length(loclambda)), ## conditional random effects already in predictor
                   stop("(!) random sample from given rand.family not yet implemented")
            )},simplify=TRUE) ## should have nsim columns
          object$rand.families[[it]]$linkfun(newU) 
        }
        newV <- lapply(seq(length(vec_n_u_h)), locfn ) ## one multi-rand.family simulation
        newV <- do.call(rbind,newV) ## each column a simulation
        eta <-  matrix(rep(eta_fixed,needed),ncol=needed) + as.matrix(ZAL %id*% newV) ## nobs rows, nsim col
        mu <- vector("list",needed)
        for (it in seq_len(needed)) mu[[it]] <- .fv_linkinv(eta[,it], object$family, object$families)
      } else stop("Unknown simulate 'type' value.")
    }
    ## ! mu := freqs for binomial, counts for poisson ## vector or matrix
    # phiW is always a matrix but mu cannot bc its elements may have attributes =>
    if ( ! inherits(mu,"list")) mu <- data.frame(mu) ## makes it always a list
    if (length(mu) != needed) stop("Programming error in simulate.HLfit() (ncol(mu) != needed).")
    #
    phiW <- .get_phiW(object=object, newdata=newdata, 
                      mu=mu, ## must be a list
                      phi_type=phi_type, needed=needed, 
                      prior.weights=prior.weights) # phiW is always a matrix
    # up to version 3.5.49, there was ad hoc code calling COMPoisson()$simulate(object, nsim=nsim) when 
    #    the mu and phi vectors are constant accross the nsim simulations, 
    #    to avoid computing the distribution (cumprodfacs in .COMP_simulate()) nsim times.
    # The updated .r_resid_var() -> appears to manage that too (without calling COMPoisson()$simulate()), 
    #  .r_resid_var_over_cols() too (a distinct issue is that multivariate mu may loose the (COMP)-lambda attribute: __FIXME__?)
    block <- .r_resid_var_over_cols(needed, mu, phiW, sizes, family=object$family, families=object$families, is_mu_fix_btwn_sims=is_mu_fix_btwn_sims,
                                    nsim=nsim)
    if (is.null(resp_testfn)) {
      if (nsim==1L) block <- drop(block)
      return(block) 
    } else {
      check_cond <- apply(block,2L, resp_testfn)
      if (is.null(resu)) {resu <- block[,check_cond,drop=FALSE]} else resu <- cbind(resu,block[,check_cond,drop=FALSE])
      done <- done+length(which(check_cond))
    }
  }
  ## we reach this point only if there is a resp_testfn, and then 'resu'
  if (nsim==1L) resu <- drop(resu)
  if (was_invColdoldList_NULL) object$envir$invColdoldList <- NULL
  return(resu)    
}


simulate.HLfitlist <- function(object,nsim=1,seed=NULL,newdata=object[[1]]$data,sizes=NULL,...) {
  ## RNG stuff copied from simulate.lm
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else { ## this makes changes to RNG local where 'seed' is used:
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  replicate(nsim, {
    allrownames <- unique(unlist(lapply(object, function(hl){rownames(hl$data)})))
    resu <- matrix(0,nrow=length(allrownames),ncol=length(object)) ## two cols if 3 types
    cumul <- 0
    if (length(sizes) != nrow(newdata)) stop("length(sizes) != nrow(newdata).")
    for (it in seq(ncol(resu))) {
      ## it = 1 se ramène à simulate(object[[1]])
      if (is.null(sizes)) sizes <- .get_BinomialDen(object[[it]])
      resu[,it] <- simulate(object[[it]],newdata=newdata,sizes=sizes - cumul,verbose=FALSE) ## FIXME use of resp_testfn would require it to be a list
      cumul <- rowSums(resu)  
    }
    resu <- cbind(resu,sizes - cumul) ## now 3 cols if 3 types
    rownames(resu) <- allrownames
    colnames(resu) <- attr(object,"sortedTypes")
    as.data.frame(resu)
  },simplify=FALSE)
}
