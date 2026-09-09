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

# newMeanFrames must have cols for the fixed beta's, absent from the object's $X.pv 
## newdata -> offset must be recomputed. 
## dans l'état actuel $fixef et complet, incluant les etaFix$beta: pas besoin de les separer
## mais il peut contenir des NA ! à enlever
# The newMeanFrames$X has columns for these fixed values from etaFix$beta, 
# contrarily to object$X.pv => don't use the cols of the latter matrix.
.newetaFix <- function(object, newMeanFrames,validnames=NULL,
                       X=newMeanFrames$X,
                       mf=newMeanFrames$mf) { 
  if (is.null(validnames)) {
    est_and_fix <- names(which(!is.na(object$fixef))) ## estimated + etaFix$beta
    validnames <- intersect(colnames(X) ,est_and_fix) # would be colnames(newMeanFrames$X) when there is no NA if object$fixef. would validnames=est_and_fix suffice ?
  }
  if (length(validnames)) {
    etaFix <-  drop(X[,validnames,drop=FALSE] %*% object$fixef[validnames]) ## valid even if ncol(newMeanFrames$X) = 0
  } else etaFix <- rep(0, nrow(X))
  off <- model.offset( mf) ### look for offset from (ori)Formula 
  if ( ! is.null(off)) etaFix <- etaFix + off   
  return(etaFix)
}

# up to version 3.5.49, there was ad hoc code calling COMPoisson()$simulate(object, nsim=nsim) when 
#    the mu and phi vectors are constant across the nsim simulations, 
#    to avoid computing the distribution (cumprodfacs in .COMP_simulate()) nsim times.
# The updated .r_resid_var() -> appears to manage that too (without calling COMPoisson()$simulate()), 
#  .r_resid_var_over_cols() too (a distinct issue is that multivariate mu may have lost 
#     the (COMP)-lambda attribute in some cases, but this should no longer occur or else would cause an error).
.r_resid_var <- function(mu, phiW, sizes, family,
                         family_par=.get_family_par(family=family),
                         zero_truncated=identical(family$zero_truncated,TRUE), 
                         famfam, nsim=1L) { 
  # we cannot use family()$simulate bc it assumes a fit object as input
  if (length(mu)) {
    resu <- switch(famfam,
                   "phi=0" = rep(mu,nsim),
                   "gaussian" = rnorm(nsim*length(mu),mean=mu,sd=sqrt(phiW)),
                   "poisson" = .rpois(nsim*length(mu),mu,zero_truncated=zero_truncated), 
                   "binomial" = rbinom(nsim*length(mu),size=sizes,prob=mu),
                   "Gamma" = {
                     y <- rgamma(nsim*length(mu), shape= 1 / phiW, scale=mu*phiW) # mean=sh*sc=mu, var=sh*sc^2 = mu^2 phiW
                     Gamma_min_y <- .spaMM.data$options$Gamma_min_y
                     is_low_y <- (y < Gamma_min_y)
                     if (any(is_low_y)) y[which(is_low_y)] <- Gamma_min_y 
                     y
                   }, ## ie shape increase with prior weights, consistent with Gamma()$simulate / spaMM_Gamma()$simulate
                   "COMPoisson" = {
                     lambdas <- attr(mu,"lambda") # F I X M E an environment would keep values ?
                     if (is.null(lambdas)) {
                       sapply(mu, function(muv) {
                         lambda <- family$mu2lambda(muv)
                         .COMP_simulate(lambda=lambda,nu=family_par, #= COMP_nu 
                                        nsim=nsim)
                       })
                     } else sapply(lambdas,.COMP_simulate,nu=family_par, #= COMP_nu 
                                   nsim=nsim)
                   },
                   "negbin" = .rnbinom(nsim*length(mu), size=family_par, #= NB_shape
                                       mu_str=mu, zero_truncated=zero_truncated),
                   "negbin1" = .rnbinom(nsim*length(mu), size=mu*family_par, #= NB_shape
                                        mu_str=mu, zero_truncated=zero_truncated),
                   "negbin2" = .rnbinom(nsim*length(mu), size=family_par, #= NB_shape
                                        mu_str=mu, zero_truncated=zero_truncated),
                   "tweedie" = .rtweedie(nsim*length(mu), p = family_par, #= Tw_index
                                                 mu=mu, phi=phiW),
                   "beta_resp" = {
                     # spaMM's phi is here 1, so phiW must be 1/prior.weights and so for the precision parameter:
                     Wfamily_par <- family_par/phiW
                     y <- rbeta(n=nsim * length(mu), 
                                shape1=mu*Wfamily_par,shape2=(1-mu)*Wfamily_par)
                     beta_min_y <- .spaMM.data$options$beta_min_y
                     is_low_y <- (y < beta_min_y)
                     if (any(is_low_y)) y[which(is_low_y)] <- beta_min_y 
                     beta_max_y <- 1- beta_min_y
                     is_high_y <- (y > beta_max_y)
                     if (any(is_high_y)) y[which(is_high_y)] <- beta_max_y 
                     y
                   }, # fam_par being the precision parameter in the Cribari-Neto parametrisation 
                   "betabin" = {
                     # phiW: same logic as for beta_resp
                     ntot <- nsim*length(mu) 
                     Wfamily_par <- family_par/phiW
                     rmu <- rbeta(n=ntot, shape1=mu*Wfamily_par,shape2=(1-mu)*Wfamily_par)
                     rbinom(ntot, size=sizes, prob=rmu)
                   },
                   stop("(!) random sampling from given family not yet implemented")
    )
  } else resu <- mu # i.e. numeric(0); such mu is possible in mv-fit if predictor variables are missing for one submodel.
  if (nsim>1L) dim(resu) <- c(length(mu),nsim)
  resu
} ## vector-valued function from vector input

.get_family_parlist <- function(object, families=object$families, family=object$family, newdata) {
  if ( ! is.null(families)) {
    family_parlist <- vector("list", length(families))
    for (mv_it in seq_along(families)) {
      family_parlist[mv_it] <- list(.get_family_par(family=families[[mv_it]], newdata=newdata)) 
    }
    return(family_parlist)
  } else .get_family_par(family=family, famfam=family$family, newdata=newdata)
}

.get_family_par <- function(family, famfam=family$family, newdata=NULL, mv_it=NULL) {
  if ( ! is.null(newdata) &&  ! is.null(resid.formula <- (disp_env <- family$resid.model)$resid.formula)) {
    residFrames <- .get_terms_info(formula=resid.formula, data=newdata, famfam="")
    # handles both "rdiOff" and "rdiForm" cases:
    if ( is.null(off <- model.offset(residFrames$mf)) ) {family_par <- 0} else family_par <- off
    if ( ! is.null(disp_env$beta)) family_par <- family_par + residFrames$X %*% disp_env$beta
  } else if (famfam =="COMPoisson") {   ##  all the next cases are OK whether there is no new data or whether there is no resid.formula
    family_par <- environment(family$aic)$nu
  } else if (famfam %in% c("beta_resp","betabin")) {
    family_par <- environment(family$aic)$prec
  } else if (famfam  %in% c("negbin","negbin1","negbin2")) {
    family_par <- environment(family$aic)$shape
  } else if (famfam  == "tweedie") {
    family_par <- environment(family$aic)$"p"
  } else family_par <- NULL
  family_par
}

# ((For nsim>1 at least)) the mu has been expanded as a list, 
# each element of which is the mu for a simulation replicate.
# *Each such element* is expected by .r_resid_var_over_cols() to bear attributes, 
# incl. for mv fits an "mv" attribute that stores a list of mu's per submodel, 
#  each possibly with ZT attributes:
# List of 1          <=  example of mu list 
# $ : Named num [1:3057] 0.494 0.514 0.535 0.555 0.575 ...     <= a simulation replicate
# ..- attr(*, "names")= chr [1:3057] "1.ld02" "2.ld02" "3.ld02" "4.ld02" ...
# ..- attr(*, "mv")=List of 3        <= list of mu's per submodel
# .. ..$ : Named num [1:1374] 0.494 0.514 0.535 0.555 0.575 ...
# .. .. ..- attr(*, "names")= chr [1:1374] "1.ld02" "2.ld02" "3.ld02" "4.ld02" ...
# .. ..$ : Named num [1:1182] 0.185 0.226 0.248 0.26 0.325 ...
# .. .. ..- attr(*, "names")= chr [1:1182] "2.fl02" "6.fl02" "8.fl02" "9.fl02" ...
# .. ..$ : Named num [1:501] 1.73 1.72 1.74 1.59 1.86 ...
# .. .. ..- attr(*, "mu_U")= num [1:501] 1.22 1.2 1.23 1.01 1.4 ...     <= third sumbodel is ZT
# .. .. ..- attr(*, "p0")= num [1:501] 0.295 0.301 0.292 0.362 0.247 ...
# .. .. ..- attr(*, "names")= chr [1:501] "16.hdct02" "18.hdct02" "19.hdct02" "24.hdct02" ...
.r_resid_var_over_cols <- function(mu,  # a list over nsim !
                                   phiW, 
                                   family_par,
                                   sizes, family, families, resp_range=NULL, is_mu_fix_btwn_sims,
                                   is_phiW_fix_btwn_sims=attr(phiW,"is_phiW_fix_btwn_sims"),
                                   nsim=1L, as_matrix=FALSE, mv_it=NULL,
                                   zero_truncated=identical(family$zero_truncated,TRUE),
                                   cum_nobs=attr(families,"cum_nobs"),
                                   phi_type) {
  if ( ! is.null(families)) {
    rowS <- vector("list", length(families))
    for (mv_it in seq_along(families)) {
      resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
      family <- families[[mv_it]] # copy needed for zero_truncated to get the correct family...
      rowS[[mv_it]] <- .r_resid_var_over_cols(mu, 
                                              phiW=phiW[resp_range,,drop=FALSE], 
                                              family_par=family_par[[mv_it]],
                                              sizes=sizes[resp_range], 
                                              family=family, families=NULL, resp_range=resp_range,
                                              is_mu_fix_btwn_sims=is_mu_fix_btwn_sims,
                                              is_phiW_fix_btwn_sims=is_phiW_fix_btwn_sims[mv_it], nsim=nsim, as_matrix=TRUE, mv_it=mv_it,
                                              zero_truncated=identical(family$zero_truncated,TRUE),
                                              phi_type=phi_type)
    }
    return(do.call(rbind, rowS))
  }
  block <- NA*phiW
  if (phi_type=="phi=0") {famfam <- phi_type} else famfam <- family$family
  if (is_mu_fix_btwn_sims && is_phiW_fix_btwn_sims) { # 2nd condition should be trivially true for count models without prior.weights,
                                                      # even those with a resid.model as it has fixed effects only.  
                                                      # and "presumably" also even for count models with prior weights (in beta_resp, at least):
                                                      # only variable prior weights btwn_sims (how?) should make it FALSE but the code ignores that possibility.
    mu_all <- mu[[1L]] 
    if ( ! is.null(mv_it)) { # multivariate
      mu_all <- structure(attr(mu_all,"mv")[[mv_it]], 
                          p0=attr(mu_all,"p0")[[mv_it]], 
                          mu_U=attr(mu_all,"mu_U")[[mv_it]])
      if (is.null(mu_all)) {
        # this should no longer occur. See comments above the .r_resid_var_over_cols() call.
        warning("Suspect structure of .r_resid_var_over_cols()'s mu argument. Zero-truncation info may be lost.")
        mu_all <- structure(mu[[1L]][resp_range], p0=attr(mu[[1L]],"p0")[[mv_it]], mu_U=attr(mu[[1L]],"mu_U")[[mv_it]])
      }
    }
    block <- .r_resid_var(mu_all, phiW=phiW[,1L],sizes=sizes, 
                          zero_truncated=zero_truncated, 
                          # cases where family_par is not needed but the is a promise available... => it may be possible to merge the codes?
                          famfam=famfam, family=family, nsim=nsim) # vector or nsim-col matrix
    if (as_matrix && is.null(dim(block))) dim(block) <- c(length(block),1L) # forces matrix for mv code using rbind()
  } else {
    for (sim_jt in seq_len(nsim)) {
      mu_jt <- mu[[sim_jt]]
      if ( ! is.null(mv_it)) {
        mu_jt <- structure(attr(mu_jt, "mv")[[mv_it]], 
                           p0=attr(mu_jt,"p0")[[mv_it]], 
                           mu_U=attr(mu_jt,"mu_U")[[mv_it]])
      } 
      block[ ,sim_jt] <- .r_resid_var(mu_jt, 
                                      phiW=phiW[ ,sim_jt], # rlevant rowas have been selected in the mv case 
                                      sizes=sizes, 
                                      family_par= family_par,  
                                      zero_truncated=zero_truncated, famfam=famfam, family=family, nsim=1L)
    }
  }
  return(block)
}

# tests with NA's: simulate(byP3) 
# 'valid_rows' do not depend on NA's in response values,
# which makes sense in some cases, BUT if using only the 'valid_rows' info,
# simulated bootstrap replicates would have more info than the original data.
# 
.r_resid_var_p4m <- function(mu=mu, object, sizes, newdata, cum_nobs, template, ...) {
  multinom_info <- object$p4m_info$multinom_info
  has_dynoffset <- multinom_info$has_dynoffset
  good_positions <- ! is.na(template) # ! (invalid rows or missing resp)
  if (is.null(newdata)) {
    valid_rows <- multinom_info$valid_rows
  } else {
    # if (is.null(newdata$".dynoffset")) stop("suspect null .dynoffset in simulate procedure") 
    valid_rows <- rowSums(good_positions)>1L #  ! is.na(newdata$".dynoffset")
    # good_positions[valid_rows,] <- TRUE # No. Missing predictors are always missing =>
    # => the template must already match the 'mu'
  }
  ntypes <- length(cum_nobs)-1L
  for (it in seq_along(mu)) {
    if (is.null(template)) {
      pred <- matrix(mu[[it]], ncol=ntypes)
    } else {
      mv_i <- attr(mu[[it]],"mv") 
      for (jt in seq_along(mv_i)) template[,jt][good_positions[,jt]] <- mv_i[[jt]]
      pred <- matrix(template, ncol=ntypes)
    }
    valid_rowsums <- rowSums(pred[valid_rows, has_dynoffset, drop=FALSE], na.rm =TRUE)
    # changing newdata$".dynoffset" here would be useless, as predict has already been called.
    pred[valid_rows,] <- pred[valid_rows,]/valid_rowsums
    mu[[it]] <- pred # that's full size including invalid rows.
  } # now mu is a list of matrices, which is the case only for p4m
  # cannot use the object$envir$missingRespInfo$missingResp here as it contains only valid_rows
  block <- lapply(mu, # designed to match number of simulation replicates (nsim) (is_mu_fix_btwn_sims?)
                  function(mu_simrep) {
                    resu <- as.vector(
                      t(sapply(seq_len(nrow(mu_simrep)), function(resp_line) {
                        if (valid_rows[resp_line]) {
                          lineprob <- mu_simrep[resp_line,]
                          NApos <- is.na(lineprob)
                          lineprob[is.na(lineprob)] <- 0 # (*) locally replaces missing-resp NAs
                          rnd <- rmultinom(1, size=sizes[resp_line],prob = lineprob)
                          rnd[NApos] <- NA_integer_
                          rnd
                        } else rep(NA_integer_, ntypes)
                      })))
                    if (is.null(newdata)) resu[ ! good_positions] <- NA_integer_ # corrects (*)
                    resu
                  }
  ) # without as.vector(), that would be an nsim-list of matrices, each of dim #responses X #types
  # Now this is is a list of vectors, so
  block <- do.call(cbind, block)
}

#################### becoming obsolete
# .calc_ZAlist_newdata_mv <- function(object, new_X_ZACblob=NULL) {
#   map_rd_mv <- attr(object$ZAlist, "map_rd_mv")
#   ori_exp_ranef_terms <- attr(object$ZAlist, "exp_ranef_terms")
#   ori_exp_ranef_strings <- attr(object$ZAlist,"exp_ranef_strings")
#   locdataS <- new_X_ZACblob$locdata # needed because it result from check of all variables needed for prediction (such as residVar predictors)
#   # This currently do not bear the variable of the ranefs that were not "conditioned upon", but these variables are needed here
#   # since this function construct all design matrices (out of which those for "conditioned upon" ranefs will be ...%*% [ranV=0])
#   loc_cum_nobs <- c(0L,cumsum(lapply(locdataS,nrow)))
#   newZAlist <- list()
#   for (mv_it in seq_along(map_rd_mv)) {
#     rd_in_mv <- map_rd_mv[[mv_it]]
#     exp_ranef_strings_it <- ori_exp_ranef_strings[rd_in_mv]
#     newdata_it <- locdataS[[mv_it]]
#     if (length(exp_ranef_strings_it)) {
#       Zlist <- .calc_Zlist(exp_ranef_terms=ori_exp_ranef_terms, data=newdata_it, 
#                            For="simulate",
#                            rmInt=0L, sparse_precision=FALSE,
#                            corr_info=.get_from_ranef_info(object), 
#                            rd_in_mv=rd_in_mv,
#                            sub_oldZAlist=object$ZAlist, # OK if we use only colnames, not attributes of the list...
#                            lcrandfamfam=attr(object$rand.families,"lcrandfamfam"))
#       amatrices <- .get_new_AMatrices(object,newdata=newdata_it, newZlist=Zlist) # 
#       ZAlist_it <- .calc_normalized_newZAlist(Zlist=Zlist,
#                                            AMatrices=amatrices,
#                                            vec_normIMRF=object$ranef_info$vec_normIMRF,
#                                            strucList=object$strucList[rd_in_mv])
#       # In the *fit* preprocessing, .merge_ZAlists is called on ZA lists for each submodel, named in ref to the submodels only;
#       # .merge_ZAlists uses "exp_ranef_strings" or similar info to match the lists, not list names. 
#       names(ZAlist_it) <- rd_in_mv
#     } else ZAlist_it <- list()
#     attr(ZAlist_it,"exp_ranef_strings") <- exp_ranef_strings_it
#     newZAlist <- .merge_ZAlists(
#       newZAlist, ZAlist_it, 
#       nobs1=loc_cum_nobs[mv_it], # will create a 0-block with nobs1 row before 
#                                  #    the nonzero block in ZAlist_it for new ranef 
#       nobs2=nrow(newdata_it), mv_it)
#   }  
#   newZAlist
# }

#################### becoming obsolete
#  Build full Zlist for simulation; ultimately 
# only elements for "marginalized upon" ranefs will be used to *simulate* ranefs
# (Z's for "conditioned upon" ranefs were already provided by .calc_new_X_ZAC[_mv] 
# and already used in the .point_predict step, without random draws).
# See comment on simulate.HLfit() further explaining what this function does.
# .calc_ZAlist_newdata <- function(object, new_X_ZACblob) {
#   # we simulate with all ranefs (treated conditionally|ranef or marginally) hence 
#   # * we need design matrices for all ranefs
#   # * we need values of all the original variables
#   # hence we use an "## effective '.noFixef'" : formula with only ranefs of the fit. But 1st version fails for mv; second may be more straightforward anyway
#   if (is.null(vec_nobs <- object$vec_nobs)) { #  *univariate*-resp model 
#     locdata <- new_X_ZACblob$locdata
#     exp_ranef_terms <- attr(object$ZAlist, "exp_ranef_terms")
#     Zlist <- .calc_Zlist(exp_ranef_terms=exp_ranef_terms, data=locdata, rmInt=0L, sparse_precision=FALSE,
#                          corr_info=.get_from_ranef_info(object),
#                          sub_oldZAlist=object$ZAlist,
#                          For="simulate",
#                          lcrandfamfam=attr(object$rand.families,"lcrandfamfam"))
#     amatrices <- .get_new_AMatrices(object,newdata=locdata, newZlist=Zlist) 
#     newZAlist <- .calc_normalized_newZAlist(Zlist=Zlist,
#                                          AMatrices=amatrices,
#                                          vec_normIMRF=object$ranef_info$vec_normIMRF, 
#                                          strucList=object$strucList)
#   } else {
#     newZAlist <- .calc_ZAlist_newdata_mv(object, new_X_ZACblob = new_X_ZACblob)
# 
#   }
#   return(newZAlist)
# }

.simulate_ranef <- function(object, rd, newdata, 
                           cum_n_u_h=attr(object$lambda,"cum_n_u_h"), 
                           vec_n_u_h=diff(cum_n_u_h), 
                           fittedLambda=object$lambda.object$lambda_est, 
                           nsim, 
                           lcrandfamfam=attr(object$rand.families,"lcrandfamfam")) {
  
  if (is.null(newdata)) { # assume no new ranef levels, 
    # & u.range must refer to fitted object's u (otherwise lambda foranother ranef might be used)
    u.range <- (cum_n_u_h[rd]+1L):(cum_n_u_h[rd+1L])
    loclambda <- fittedLambda[u.range] ## includes prior_lam_fac
  } else { # rebuild lambda, handled new levels at least in last case. 
    if ( ! is.null(object$rand.families[[rd]]$prior_lam_fac)) { 
      # prior_lam_fac is the 'design' for non-ranCoef (wei-1|.)
      leftOfBar_terms <- attr(object$ZAlist,"exp_ranef_terms")[[rd]][[2L]]
      leftOfBar_mf <- model.frame(as.formula(paste("~",leftOfBar_terms)), newdata, xlev = NULL) 
      prior_lam_fac <- leftOfBar_mf[,1L]^2 ## assumes simple syntax (wei-1|.)
      loclambda <- object$lambda.object$lambda_list[[rd]]* prior_lam_fac
    } else loclambda <- object$lambda.object$lambda_list[[rd]] # scalar
  }
  
  nr <- vec_n_u_h[rd]
  newU <- replicate(nsim, {
    switch(lcrandfamfam[rd], ## remainder of code should be OK for rand.families
           "gaussian" = rnorm(nr,sd=sqrt(loclambda)),
           "gamma" = rgamma(nr,shape=1/loclambda,scale=loclambda),
           "beta" = rbeta(nr,1/(2*loclambda),1/(2*loclambda)),
           "inverse.gamma" = 1/rgamma(nr,shape=1+1/loclambda,scale=loclambda), ## yields inverse gamma (1+1/object$lambda,1/object$lambda)
           "conditional"= rep(0, nr), ## conditional random effects already in predictor
           stop("(!) random sample from given rand.family not yet implemented")
    )},simplify=TRUE) ## should have nsim columns
  object$rand.families[[rd]]$linkfun(newU) 
}

simulate_ranef <- function(object, which=NULL, newdata=NULL, nsim=1L) {
  cum_n_u_h <- attr(object$lambda,"cum_n_u_h")
  vec_n_u_h <- diff(cum_n_u_h)
  if (is.null(which)) which <- seq_along(vec_n_u_h)
  nwhich <- length(which)
  newb <- vector("list", nwhich)
  for (rd in seq_along(which)) {
    newb[[rd]] <- .simulate_ranef(object=object, rd=which[rd], newdata=newdata, 
                                  cum_n_u_h=cum_n_u_h, 
                                  vec_n_u_h=vec_n_u_h, 
                                  fittedLambda=object$lambda.object$lambda_est, 
                                  nsim=nsim, 
                                  lcrandfamfam=attr(object$rand.families,"lcrandfamfam"))
  }
  if (nsim>1L) {
    newb <- do.call(rbind, newb)
  } else newb <- unlist(newb)
  newb
}

.warn_size_mismatch <- function(expected, isNullnewData, orisizes, vec_nobs) {
  # Providing a right-sized value facilitates automated tests => warning rather than stop
  # But __F I X M E___ make this despendent on _LOCAL_TESTS_ ?
  if (length(vec_nobs)>1L) {
    warnmess <- paste0("Specifying 'sizes' of length ", expected, " (",
                       paste(vec_nobs,collapse="+")," for respective submodels) is necessary")
  } else warnmess <- paste0("Specifying 'sizes' of length ", expected, " is necessary")
  if ( ! isNullnewData) warnmess <- paste0(warnmess, " for these 'newdata'.\n")
  if (is.null(orisizes))  {
    warnmess <- paste0(warnmess, " NULL")
  } else warnmess <- paste0(warnmess, " Wrong-sized")
  warning(paste0(warnmess, " value is replaced by unit sizes"),
          immediate. = TRUE) 
}

# the previous version made more effort to re-use fit values using a vec_nobs[fam_it] %% ori_vec_nobs[fam_it] test
.guess_new_BinomialDen <- function(sizes, mu, cum_nobs, isNullnewData,
                                   famfams,
                                   size_control=famfams %in% c("binomial","betabin")) {
  if (length(famfams)>1L) { # mvfit
    expected <- cum_nobs[length(cum_nobs)]
    size_mismatch <- length(sizes) != expected
    if (size_mismatch) {
      if (any(size_control)) .warn_size_mismatch(expected=expected, isNullnewData, 
                                                 orisizes=sizes, vec_nobs=diff(cum_nobs))
      sizes <- rep(1L, expected) 
    } 
  } else {
    expected <- length(mu[[1]])
    size_mismatch <-  ! length(sizes) %in% c(expected,1L) # looser non-API control for univariate (unify?)
    if (size_mismatch) {
      if (any(size_control)) .warn_size_mismatch(expected=expected, isNullnewData, 
                                                 orisizes=sizes, vec_nobs=expected)
      sizes <- rep(1L, expected) 
    } 
  }
  sizes
}

.warn_pw_mismatch <- function(chr_expected, isNullnewData, default_pw) {
  # Providing a right-sized value facilitates automated tests => warning rather than stop
  # But __F I X M E___ make this despendent on _LOCAL_TESTS_ ?
  warnmess <- paste0("Specifying 'prior.weights' list (lengths=", chr_expected, " is necessary")
  if ( ! isNullnewData) warnmess <- paste0(warnmess, " for these 'newdata'.\n")
  warnmess <- paste0(warnmess, "Missing of invalid values are replaced by ")
  if (length(default_pw)>1L) {
    warnmess <- paste0(warnmess,"(",paste(default_pw,collapse=","), ") for respective submodels.")
  } else warnmess <- paste0(warnmess, default_pw,".")
  warning(warnmess, immediate. = TRUE) 
}

.check_simulate_pw <- function(prior.weights, mu, cum_nobs, famfams,
                               fit_pw, isNullnewData, 
                               pw_control= ! famfams %in% c("binomial", "poisson", "COMPoisson", 
                                                            "negbin1", "negbin2") # the latter ones have no pw but still a resid disp param
) {
  # Providing a right-sized value facilitates automated tests
  n_mv <- length(famfams)
  if (n_mv>1L) { # mvfit
    if ( ! is.list(prior.weights)) prior.weights <- vector("list", length(famfams))
    vec_nobs <- diff(cum_nobs)
    pw_mismatch <- sapply(prior.weights, length)!=vec_nobs
    if (any(pw_mismatch)) {
      ambiguous_pw <- pw_mismatch & pw_control
      if ( any(ambiguous_pw)) {
        are_units <- sapply(fit_pw, attr, which="is_unit")
        if (is.list(are_units)) are_units <- sapply(are_units, identical, y=TRUE)
        defaults <- rep(NA_real_, n_mv)
        defaults[are_units] <- 1
        if (anyNA(defaults)) {
          are_unique <- sapply(fit_pw, attr, which="is_unique")
          if (is.list(are_unique)) are_unique <- sapply(are_unique, identical, y=TRUE)
          for (mv_it in seq_along(n_mv)) {
            if (is.na(defaults[mv_it]) && are_unique[mv_it]) defaults[mv_it] <- fit_pw[[mv_it]][1]
          }
        }
        if (anyNA(defaults)) {
          defaults[is.na(defaults)] <- 1
          .warn_pw_mismatch(chr_expected=paste(vec_nobs, collapse=","), isNullnewData, 
                            default_pw=defaults)
        }
      } else defaults <- rep(1, n_mv)
      for (mv_it in which(pw_mismatch)) {
        prior.weights[[mv_it]] <- structure(rep(defaults[mv_it],vec_nobs[mv_it]), unique=TRUE)
      }
    }
  } else { # univariate
    if ( length(prior.weights) != length(mu[[1]])) {
      if (pw_control &&
          ! identical(attr(fit_pw,"is_unit"),TRUE) ) { # -> warn
        if (identical(attr(fit_pw,"unique"),TRUE)) {
          default_pw <- fit_pw[1]
        } else default_pw <- 1
        .warn_pw_mismatch(chr_expected=length(mu[[1]]), isNullnewData, default_pw=default_pw)
      } else default_pw <- 1 # no pw_control or pw is unit
      prior.weights <- structure(rep(default_pw,length(mu[[1]])), unique=TRUE)
    }
  }
  prior.weights
}

.wrap_compute_ZALlist4simulate <- function(new_X_ZACblob, newZAlist, strucList) {
  L_newLv_newLv_list <- new_X_ZACblob$L_newLv_newLv_list
  newinold <- new_X_ZACblob$newinold
  for (new_rd in seq_along(newinold)) {
    if (is.null(L_newLv_newLv_list[[new_rd]])) {
      if (! is.null(new_X_ZACblob$cov_newLv_newLv_list[[new_rd]])) {
        # One would wish to check whether
        # corr.model <- attr(object$strucList[[old_rd]],"corr.model")
        # is a corrFamily, but old_rf info is not available. Wait for pb to occur...
        warning("Contact the maintainer about presumably inefficient code for computation of L_newLv_newLv_list...")
        L_newLv_newLv_list[[new_rd]] <- mat_sqrt(new_X_ZACblob$cov_newLv_newLv_list[[new_rd]])
      } else { # e.g. for IMRF *no Cnn*, no Lnn
        old_rd <- newinold[new_rd]
        L_newLv_newLv_list[[new_rd]] <- strucList[[old_rd]]
      }
    }
  }
  .compute_ZAXlist(ZAlist=newZAlist, XMatrix = L_newLv_newLv_list, 
                   cols_from_RHS=FALSE) # may be "notBindable"
}

# is.null(newdata), & re.form=NA (typical marginal-type simulate() case)
# eta_fixed_cond <- predict(..., re.form=[NA]) 
# built a new X including rows for NAs in response values,
# but ZAL is the one from the fitted object, excluding such rows.
# The new X has the extra row at least bc 
# .calc_new_X_ZAC[_mv]() emphatically removes the resp variable from the checked ones.
# For this specific combination of arguments, we might consider not removing it...?

# In p4m fits, distinct info is available, and there are possibly distinct constraints:
# rows with one missing response 
# contain info from other responses and remain in the object's $data.
# Further the $cum_nobs associated to the new X in new_X_ZACblob 
# also counts the extra rows.
# [ The putative 'mv' attribute does not contain NAs but is not always present.
# attr(.,"mv") is added by .fv_linkinv() and (provisorily at least)
# by .predict.pois4mlogit(). Only in the p4m case it may be present here. ]
.provide_oriLinesInfo <- function(object, # always needed;
                                  newX_oldZACblob # needed only to build the info.
                                  ) {
  if ( is.null(oriLinesInfo <- object$envir$oriLinesInfo)) {
    newXnames <- rownames(newX_oldZACblob$newX.pv)
    if (inherits(object,"fitmv")) {
      newX_cum_nobs <- newX_oldZACblob$cum_nobs
      data <- object$data
      oldXnamelist <- attr(data,"validrownames")
      n_submodels <- length(oldXnamelist)
      template_NAall <- matrix(NA_integer_, ncol=n_submodels, nrow=nrow(data),
                              dimnames = list(rownames(data), NULL))
      # template_NApred <- template_NAall
      oriLines <- vector("list", n_submodels)
      for (subm in seq_len(n_submodels)) {
        newXsubrnge <- .subrange(newX_cum_nobs, subm)
        oriLines[[subm]] <- ! is.na(match(newXnames[newXsubrnge], oldXnamelist[[subm]])) 
        # template_NApred[newXnames[newXsubrnge],subm] <- 0L
        template_NAall[oldXnamelist[[subm]],subm] <- 0L
      } 
      # if ( ! anyNA(template_NApred)) template_NApred <- NULL
      if ( ! anyNA(template_NAall)) template_NAall <- NULL
      oriLinesInfo <- list(oriLines=.unlist(oriLines), # vector of T/F (this is used)
                           # template_NApred=template_NApred,
                           template_NAall=template_NAall ) 
    } else { # case presumably never used
      oldXnames <- names(object$fv)
      oriLines <- ! is.na(match(newXnames,oldXnames))
      oriLinesInfo <- list(oriLines=oriLines)
    }
    object$envir$oriLinesInfo <- oriLinesInfo
  }
  oriLinesInfo
}

.provide_newLinesInfo <- function(object, # always needed;
                                  new_X_ZACblob, # needed only to build the info.
                                  newdata
) {
  newXnames <- rownames(new_X_ZACblob$newX.pv)
  newX_cum_nobs <- new_X_ZACblob$cum_nobs
  n_submodels <- length(formula(object))
  template_NApred <- matrix(NA_integer_, ncol=n_submodels, nrow=nrow(newdata),
                            dimnames = list(rownames(newdata), NULL))
  oriLines <- vector("list", n_submodels)
  for (subm in seq_len(n_submodels)) {
    newXsubrnge <- .subrange(newX_cum_nobs, subm)
    template_NApred[newXnames[newXsubrnge],subm] <- 0L
  } 
  if ( ! anyNA(template_NApred)) template_NApred <- NULL
  list(template_NAall=template_NApred) 
}

.conditioned_upon <- function(object, re.form, pred_type, type, nrand) {
  if (inherits(re.form,"formula")) {
    if (pred_type=="predVar_s.lato") warning("Non-default 're.form' is *currently* ignored when type='",type,"'.")
    re.form <- .preprocess_formula(re.form)
    ori_exp_ranef_strings <- attr(object$ZAlist,"exp_ranef_strings")
    new_exp_ranef_strings <- .process_bars(re.form,expand=TRUE)
    conditioned_upon <- .unlist(lapply(new_exp_ranef_strings, `==`, y= ori_exp_ranef_strings)) 
  } else if (is.null(re.form)) {
    conditioned_upon <- rep(TRUE,nrand)
  } else if (is.na(re.form)) { 
    if (pred_type=="predVar_s.lato") warning("Non-default 're.form' is *currently* ignored when type='",type,"'.")
    if (is.na(re.form)) {
      conditioned_upon <- rep(FALSE, nrand)
    } else conditioned_upon <- rep(TRUE,nrand)
  } 
  conditioned_upon
}

# simulate.HLfit(fullm[[2]],newdata=fullm[[1]]$data,size=fullm[[1]]$data$total) for multinomial avec binomial nichées de dimension différentes

simulate.HLfit <- function(object, nsim = 1, seed = NULL, newdata=NULL,
                           type = "marginal", re.form, conditional=NULL, 
                           verbose=c(type=TRUE, showpbar= eval(spaMM.getOption("barstyle"))),
                           sizes=if (is.null(newdata)) get_drawSizes(object, p4m="M"), 
                           resp_testfn=NULL, phi_type="predict", 
                           prior.weights= if (is.null(newdata)) object$prior.weights, 
                           variances=list(), ...) { ## object must have class HLfit; corr pars are not used, but the ZAL matrix is.
  
  if (inherits(newdata,"tibble")) newdata <- as.data.frame(newdata) 
  ## RNG stuff copied from simulate.lm
  control <- list(simulate=TRUE,
                  keep_ranef_covs_for_simulate=FALSE) # modified in one case below
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
    message("simulate does not yet work on list of fits. Either run simulate on each")
    message("  of the individual fits, or use fitmv() or pois4mlogit() instead.")
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
  variances$residVar <- TRUE # so that new_X_ZACblob will provide variables 
  # for variances of drawn residuals, when passed to .get_phiW().   
  if (isNullUserSizes <- is.null(sizes)) sizes <- .get_BinomialDen(object)
  nrand <- length(object$ZAlist)
  if (nrand>0L) {
    if ( missing(re.form)) {
      if (type=="marginal") {
        re.form <- NA # Does not mean that ranefs are entirely ignored as uuCnewnew is computed when control$keep_ranef_covs_for_simulate is TRUE.
      } else if (type=="residual") {
        re.form <- NULL
      } else if (pred_type=="predVar_s.lato") {
        # type "predVar" leaves 're.form' missing. as it is not used (which should be equivalent to re.form=NULL)
      } else if (type=="conditional") {
        stop("'conditional' is note a valid simulate type. Consider type='residual', or some re.form value?")
      } else stop("Unhandled value of 'type' argument in simulate.HLfit().")
    } 
    if ( ! missing(re.form)) {
      conditioned_upon <- .conditioned_upon(object, re.form, pred_type, type, nrand)
      control$marginalized <- ! conditioned_upon
    }
  }
  resu <- NULL
  done <- 0L
  verbtype <- verbose[["type"]]
  is_mu_fix_btwn_sims <- FALSE
  
  if (length(object$families)) {
    famfams <- sapply(object$families, `[[`, x="family")
  } else famfams <- object$family$family
  
  while((needed <- nsim-done)) { ## loop operates only for resp_testfn
    if (nrand==0L) { ## note that replicate mu's can still be variable for non-standard pred_type
      
      if (pred_type=="predVar_s.lato") { ## re.form ignored so de facto NULL
        if (type=="(ranef|response)") {
          stop("meaningless argument type='(ranef|response)' for a fixed-effect model")
        } else if (verbtype) cat("Simulation from linear predictor variance | observed response:\n") 
        variances$cov <- (NROW(newdata)!=1L)
        control$fix_predVar <- NA
        point_pred_eta <- predict(object,newdata=newdata, type="link", control=control,
                                  variances=variances, verbose=verbose, ...) 
        predVar <- attr(point_pred_eta,"predVar")
        if (is.null(predVar)) stop("A 'variances' argument should be provided so that prediction variances are computed.") 
        rand_eta <- .mvrnorm(n=needed,mu=point_pred_eta[,1L], Sigma=predVar)
        if (needed>1L) rand_eta <- t(rand_eta) ## else mvrnorn value is a vector
        new_X_ZACblob <- attr(point_pred_eta,"new_X_ZACblob")
        cum_nobs <- new_X_ZACblob$cum_nobs # presence expected given control$simulate=TRUE; will be needed for .calc_phiW()
        mu <- .fv_linkinv(eta=rand_eta, family=object$family, families=object$families, cum_nobs=cum_nobs) 
      } else { # standard simulation withOUT ranefs
        control$fix_predVar <- FALSE
        mu <- predict(object,newdata=newdata,binding=NA,control=control, verbose=verbose)
        new_X_ZACblob <- attr(mu,"new_X_ZACblob")
        cum_nobs <- new_X_ZACblob$cum_nobs # presence expected given control$simulate=TRUE; will be needed for .calc_phiW()
        attr(mu,"new_X_ZACblob") <- NULL
        mu <- replicate(needed,mu,simplify=FALSE) # always a list at this stage
        is_mu_fix_btwn_sims <- TRUE
      }
      if (inherits(object,"fitmv") && is.null(newdata)) abyss <- 
          .provide_oriLinesInfo(object, newX_oldZACblob=new_X_ZACblob)
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
          control$fix_predVar <- NA
          point_pred_eta <- predict(object,newdata=newdata, type="link", control=control,
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
          } else rand_eta <- .mvrnorm(n=needed,mu=point_pred_eta[,1L], Sigma=predVar) # n=needed means we will get nsim distinct eta vectors
          if (needed>1L) rand_eta <- t(rand_eta) ## else mvrnorn value is a vector
          new_X_ZACblob <- attr(point_pred_eta,"new_X_ZACblob")
          cum_nobs <- new_X_ZACblob$cum_nobs # presence expected given control$simulate=TRUE
          mu <- .fv_linkinv(eta=rand_eta, family=object$family, families=object$families, 
                            cum_nobs=cum_nobs) ## ! freqs for binomial, counts for poisson: suitable for final code
        } else stop("This conditional simulation is not implemented for non-gaussian random-effects")
      } else if ( is.null(re.form)) { ## conditional on (all) predicted ranefs, type = "residual"
      # } else if (type=="residual") { ## conditional on (all) predicted ranefs,
        if (verbtype) cat("simulation of residuals, conditional on point predictions (hence on random effects):\n") 
        control$fix_predVar <- FALSE
        mu <- predict(object,newdata=newdata,binding=NA,control=control, verbose=verbose, ...)
        new_X_ZACblob <- attr(mu,"new_X_ZACblob")
        cum_nobs <- new_X_ZACblob$cum_nobs # presence expected given control$simulate=TRUE; will be needed for .calc_phiW()
        attr(mu,"new_X_ZACblob") <- NULL
        mu <- replicate(needed,mu,simplify=FALSE) #matrix(rep(mu,nsim),ncol=nsim)
        is_mu_fix_btwn_sims <- TRUE
      } else if ( inherits(re.form,"formula") || is.na(re.form) ) { ## explicit re.form; or {unconditional MIXED MODEL, type= "marginal" }
        if (verbtype) {
          if (inherits(re.form,"formula")) {
            cat("Simulation conditional on random effect(s) retained in 're.form':\n")
          } else cat("Unconditional simulation:\n") 
        }
        # mu_fixed <- predict(object, newdata=newdata, re.form=re.form,binding=NA,control=list(fix_predVar=FALSE), verbose=verbose) ## mu_T
        # eta_fixed <- .eta_linkfun(mu_fixed, object$family, object$families)
        control$fix_predVar <- FALSE
        # we will need a ZAL and below we need the newdata to construct it if they are not NULL:
        control$keep_ranef_covs_for_simulate <- ( ! is.null(newdata))
        # Includes the predicted value of ranefs conditioned upon:
        eta_fixed_cond <- predict(object, newdata=newdata, type="link", variances=variances,
                             re.form=re.form, #  
                             binding=NA,control=control, 
                             verbose=verbose, ...)
        # ... for that purpose, re.form is used to identify those ranefs and to provide the design matrices for them
        # These design matrices will a priori not be further used (or only in a trivial way,  times zero-valued ranefs),
        # but other elements of the new_X_ZACblob attribute will be used.
        new_X_ZACblob <- attr(eta_fixed_cond,"new_X_ZACblob")
        cum_nobs <- new_X_ZACblob$cum_nobs # presence expected given control$simulate=TRUE
        attr(eta_fixed_cond,"new_X_ZACblob") <- NULL
        if (is.null(newdata)) { ## we simulate with all ranefs (treated conditionally|ranef or marginally) hence no selection of matrix
          ZAL <- get_ZALMatrix(object, force_bind = ! (.is_spprec_fit(object)) )
          cum_n_u_h <- attr(object$lambda,"cum_n_u_h")
          vec_n_u_h <- diff(cum_n_u_h)
        } else { # new sampling design
          # if (inherits(object,"fitmv")) {
          #   newZAlist <- .calc_ZAlist_newdata_mv(object, new_X_ZACblob=new_X_ZACblob) # wherein new_X_ZACblob$newZAlist not clearly used 
          # } else 
            newZAlist <- new_X_ZACblob$newZAlist
          # contains matrices for all ranefs (cond or marg)
          # Likewise newV will contain elements for all ranefs (those conditioned-upon being 0)
          ZALlist <- .wrap_compute_ZALlist4simulate(new_X_ZACblob, newZAlist, object$strucList)       
          ##   
          # ZAL <- .ad_hoc_cbind(ZALlist, as_matrix=FALSE ) # inappropriate for IMRF andother ZAXlist stuff
          ZAL <- .compute_ZAL(XMatrix=NULL, ZAlist=ZALlist, as_matrix=FALSE, bind.=TRUE)  # may be a ZAXlist if ZALlist is "notBindable"
          #
          vec_n_u_h <- unlist(lapply(ZALlist,.ncol)) ## nb cols each design matrix = nb realizations each ranef
          cum_n_u_h <- cumsum(c(0,vec_n_u_h))
        }
        lcrandfamfam <- attr(object$rand.families,"lcrandfamfam") ## unlist(lapply(object$rand.families, function(rf) {tolower(rf$family)}))
        lcrandfamfam[which(conditioned_upon)] <- "conditional" 
        fittedLambda <- object$lambda.object$lambda_est
        newV <- vector("list", length(vec_n_u_h))
        for (rd in seq_along(newV)) {
          newV[[rd]] <- .simulate_ranef(object, rd=rd, vec_n_u_h=vec_n_u_h, cum_n_u_h=cum_n_u_h, 
                                        newdata=newdata, fittedLambda=fittedLambda, 
                                       nsim=needed, lcrandfamfam=lcrandfamfam)
        } ## one multi-rand.family simulation. ranefs "conditioned upon" should already be in the eta_fixed_cond value.
        newV <- do.call(rbind,newV) ## each column a simulation
        
        # Fix occasional effect of missing responses AND 
        # Provide info for patching NA's
        if (inherits(object,"pois4mlogit")) { 
          if (length(eta_fixed_cond) != cum_nobs[length(cum_nobs)]) {
            # Here ZAL and eta_fixed_cond already appear to match
            # cf simulate(byP3m)
            # But cum_nobs does not match them.
            # mv <- attr(eta_fixed_cond,"mv") # provided by .predict.pois4mlogit()
            cum_nobs <- c(0L, cumsum(attr(eta_fixed_cond,"nobs")))
            # eta_fixed_cond <- .unlist(mv) # no NA's
          }
        } else if (inherits(object,"fitmv")) { 
          # Here ZAL and eta_fixed_cond already may not match eta_fixed_cond
          # cf test simulate(P3):
          # nrow(ZAL) is 27 (old Z, cf 'newX_oldZACblob'), as in nrow 'old' X, 
          # but length(eta_fixed_cond) is 28 =  nrow(new_X_ZACblob$newX.pv) 
          # length(eta_fixed_cond) matches (new_X_ZACblob$)cum_nobs, hence
          # testing length(eta_fixed_cond) != cum_nobs[length(cum_nobs)] would be incorrect:
          if (nrow(ZAL) != length(eta_fixed_cond)) {
            if (is.null(newdata)) { 
              # Typical marginal-type simulate() ... as above, but quite different handling
              oriLines <-  .provide_oriLinesInfo(object, newX_oldZACblob=new_X_ZACblob)$oriLines
              eta_fixed_cond <- eta_fixed_cond[oriLines]
              cum_nobs <- attr(object$families,"cum_nobs") 
            } else stop("Unexpected case where nrow(ZAL) != length(eta_fixed_cond)")
          }
        }
        
        eta <-  matrix(rep(eta_fixed_cond,needed),ncol=needed) + as.matrix(ZAL %id*% newV) ## nobs rows, nsim col
        mu <- vector("list",needed)
        for (it in seq_len(needed)) mu[[it]] <- .fv_linkinv(eta[,it], object$family, object$families, cum_nobs = cum_nobs)
      } else stop("Unknown simulate 'type' value.")
    }
    ## ! mu := freqs for binomial, counts for poisson ## vector or matrix
    # phiW is always a matrix but mu cannot bc its elements may have attributes =>
    if ( ! inherits(mu,"list")) mu <- data.frame(mu) ## makes it always a list
    if (length(mu) != needed) stop("Programming error in simulate.HLfit() (ncol(mu) != needed).")
    #
    prior.weights <- .check_simulate_pw(prior.weights, mu, cum_nobs, famfams, 
                                        fit_pw=object$prior.weights, isNullnewData=is.null(newdata))
    phiW <- .get_phiW(object=object, newdata=newdata, 
                      newframes_info = new_X_ZACblob,  
                      dims=c(length(mu[[1]]), length(mu)), # (nrow= response length, ncol= # of replicates)
                      phi_type=phi_type, nsim=needed, 
                      prior.weights=prior.weights) # returned 'phiW' is always a matrix even for mv fits.
    family_par <- .get_family_parlist(object, newdata=newdata)

    # For some time, marginal simulation of newV's followed by mu[[it]] <- .fv_linkinv(eta[,it]...)
    # lost the "mv" attribute (and thus the included ZT info). This has been corrected.
    
    if (inherits(object,"pois4mlogit")) {
      if (is.null(newdata)) {
        linesInfo <- .provide_oriLinesInfo(object, newX_oldZACblob = new_X_ZACblob)
      } else linesInfo <- .provide_newLinesInfo(object, new_X_ZACblob = new_X_ZACblob, newdata)
      template_NA <- linesInfo$template_NAall
      # This re-inserts NAs where appropriate:
      block <- .r_resid_var_p4m(mu=mu, 
                                object=object,
                                sizes=sizes, 
                                newdata=newdata,
                                cum_nobs = cum_nobs, 
                                template=template_NA,
                                ...)
    } else {
      sizes <- .guess_new_BinomialDen(sizes=sizes, mu=mu, cum_nobs=cum_nobs, 
                                      isNullnewData=is.null(newdata), famfams=famfams) 
      # See comments on .r_resid_var_over_cols() for the format of mu
      block <- .r_resid_var_over_cols(mu, 
                                    phiW,
                                    family_par=family_par,
                                    sizes=sizes, family=object$family, families=object$families, is_mu_fix_btwn_sims=is_mu_fix_btwn_sims,
                                    nsim=needed, cum_nobs = cum_nobs, phi_type=phi_type, ...)
      
      # Re-insert NAs where appropriate:
      if (inherits(object,"fitmv")) {
        if (is.null(newdata)) {
          linesInfo <- .provide_oriLinesInfo(object, newX_oldZACblob = new_X_ZACblob)
        } else linesInfo <- .provide_newLinesInfo(object, new_X_ZACblob = new_X_ZACblob, newdata)
        template_NA <- linesInfo$template_NAall
        if ( ! is.null(template_NA)) {
          block <- apply(block, 2L, function(v) {
            template_NA[! is.na(template_NA)] <- v 
            template_NA
          })
        }
      }
    }
    
    if (is.null(resp_testfn)) {
      if (nsim==1L) block <- drop(block)
      attr(block,"nobs") <- diff(cum_nobs) 
      class(block) <- c("spaMM_simulations",class(block))
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
  class(resu) <- c("spaMM_simulations",class(resu))
  return(resu)    
} 

print.spaMM_simulations <- function (x, expanded=FALSE, ...) {
  # first print version without attribute
  if (is.null(dim(x))) {
    print(as.vector(x))
  } else print(structure(as.vector(x), dim=dim(x))) # as.matrix 
  cat("with attributes:")
  if (expanded) { # shows structure of attributes as in utils:::str.default
    .show_str_attrs(attributes(x))
  } else {
    cat(" ")
    std.attr <- c("names","dim","dimnames","class") ## attributes not to be shown
    a <- attributes(x)
    nam <- names(a)
    nam <- setdiff(nam,std.attr)
    cat(paste(nam,collapse=", "))  
    cat("\n")
  }
  invisible()
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
      ## it = 1 reduces to simulate(object[[1]])
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
