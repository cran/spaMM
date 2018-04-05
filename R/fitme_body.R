## next steps are COMPoisson optim, detection des modeles lambda, gestion des lambda multiples
#  et meilleurs default phi, lambda 
fitme_body <- function(processed,
                       init=list(),
#                       init.HLfit=list(),
                       fixed=list(), ## NULL not valid (should be handled in preprocess?)
                       lower=list(),upper=list(),
                       # control.dist=list(), ## info in processed
                       control=list(), ## optimizer, <optimizer controls>, precision
                       ... ## cf dotnames processing below 
) {
  dotlist <- list(...) ## forces evaluations, which makes programming easier...
  if (is.list(processed)) {
    proc1 <- processed[[1L]]
  } else proc1 <- processed
  verbose <-  proc1$verbose
  HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                names(formals(designL.from.Corr)),names(formals(make_scaled_dist))))  ## cf parallel code in HLCor.obj
  ## fill HLCor.args
  good_dotnames <- intersect(names(dotlist),HLnames) ## those specifically for the called fns as def'd by HLnames
  if (length(good_dotnames)) {
    HLCor.args <- dotlist[good_dotnames]
  } else HLCor.args <- list() 
  ## replace some HLCor.args members  
  if (  is.list(processed) ) {
    pnames <- names(processed[[1]])
  } else pnames <- names(processed)
  for (st in pnames) HLCor.args[st] <- NULL 
  # 'processed' may be modified below, then will be copied in HLCor.args (and then removed from this envir for safety)
  optim.scale <- control[["optim.scale"]] 
  if (is.null(optim.scale)) optim.scale="transformed" ## currently no public alternative
  sparse_precision <- processed$sparsePrecisionBOOL
  #
  corr_info <- proc1$corr_info 
  #
  # modify HLCor.args and <>bounds;   ## distMatrix or uniqueGeo potentially added to HLCor.args:
  user_init_optim <- init ## user_init_optim only for a check in new_locoptim, true initial value init.optim is modified below
  corr_types <- corr_info$corr_types
  fixed <- .post_process_fixed(fixed,corr_types=corr_types)
  init.optim <- .post_process_parlist(init,corr_types=corr_types)
  init.HLfit <- processed$init_HLfit #to be modified below ## dhglm uses fitme_body not (fitme-> .preprocess) => dhglm code modifies processed$init_HLfit
  init <- NaN ## make clear it's not to be used
  spatial_terms <- attr(proc1$ZAlist,'exp_spatial_terms')
  adjrho_ranges <- vector("list",length(corr_types))
  #
  ##### init.optim$phi/lambda will affect calc_inits -> calc_inits_dispPars.
  # outer estim seems useful when we can suppress all inner estim (thus the hatval calculations). 
  # ./. Therefore, we need to identify all cases where phi is fixed, 
  # ./. or can be outer optimized jointly with lambda:  
  phimodel <- proc1$models[['phi']]
  #if (phimodel=="phiGLM") {message("'fitme' not optimized for models with structured dispersion.")} ## FR->FR test this later"
  ranCoefs_blob <- proc1$ranCoefs_blob
  ## trying to guess all cases where optimization is useful. But FIXME: create all init and decide afterwardsS
  if ( (is_MixedM <- ( ! is.null(ranCoefs_blob) )) && (
    (var_ranCoefs <- (ranCoefs_blob$isRandomSlope & ! ranCoefs_blob$is_set)) ||
    phimodel == "phiScal" || # lambda + phi
    length(var_ranCoefs)>1L || # +s lambda
    length(corr_types[ ! is.na(corr_types)]) # lambda + corr pars
  )
  ) { ## prefer outer optimization if var_ranCoefs or not resid.model
    nrand <- length(processed$ZAlist)
    # (1) set (or not) outer optimization for phi: 
    phi.Fix <- proc1$phi.Fix
    ## All cases where phi is fixed, or can be outer optimized jointly with lambda:  
    
    ## Then determine a single phi value
    if (is.null(phi.Fix)) {
      if (phimodel == "phiScal") {
        if (NROW(processed$y)>200L) {
          ## default method is outer but can be reversed if init is NaN (and to test NaN we need to test NULL)
          not_inner_phi <- is.null(init.optim$phi) || ! is.nan(init.optim$phi) ## outer if NULL, NA or numeric
        } else {
          ## default is inner but can be reversed if numeric or NA (hence neither NULL nor NaN)
          not_inner_phi <- ! (is.null(init.optim$phi) || is.nan(init.optim$phi))
        }
      } else not_inner_phi <- FALSE
      if (not_inner_phi) {
        ## .get_init_phi takes more account of ranefs than .get_inits_by_glm(.)$phi_est, but it may return NULL.
        #  It seems (?) marginally better hence we try it first and check the result.
        if (is.null(init.optim$phi)) { 
          # if (is.call(processed$prior.weights)) {
          #   init.optim$phi <- .get_init_phi(processed,weights=NULL)
          # } else 
          init.optim$phi <- .get_init_phi(processed,weights=eval(processed$prior.weights)) ## _F I X M E_ test the eval()
          ## (fixme): fitme may fail obscurely if .get_init_phi(processed) fails silently
          if (is.null(init.optim$phi)) { ## .get_init_phi returned NULL if ./. 
            # ./. no replicate obs is available, or
            # ./. ULI too large ?? Grey area here (fixme)
            if (TRUE) {
              init.optim$phi <- .get_inits_by_glm(processed)$phi_est/(nrand+1L) ## at least one initial value should represent high guessed variance
              # if it's too low (as in min(.,2)) then fitme(Reaction ~ Days + AR1(1|Days) + (Days|Subject), data = sleepstudy) is poor
            } else {
              #init.optim$phi <- 2*sqrt(1.00001+.get_inits_by_glm(processed)$phi_est/(nrand+1L))-1 ## less coherent with general logic 
              init.optim$phi <- log(1.00001+.get_inits_by_glm(processed)$phi_est/(nrand+1L)) ## less coherent with general logic 
            }
          }  
          # init.optim$phi <- max(1e-4,init.optim$phi) ## Controlled by .calc_inits_dispPars()...
        }      
      }
    } else not_inner_phi <- TRUE ## not outer as well. Hence not_inner_phi != outer phi...
    # (2) set outer optimization for lambda and ranCoefs (handling incomplete ranFix$lambda vectors)
    if (is_MixedM) { 
      lFix <- proc1$lambda.Fix ## only for 'simple" ranefs with Xi_cols=1
      # Not super clear why I considered nranterms (user level ??) instead of nrand. FIXME.
      nranterms <- sum(var_ranCoefs | is.na(lFix)) ## var_ranCOefs has FALSE elements for non-ranCoefs (hence it is full-length)
      if (not_inner_phi) { ## Tests show it is very inefficient to use outer optim on lambda (at least) when phi must be inner optimized
        optim_lambda_with_NAs <- .reformat_init_lambda_with_NAs(init.optim$lambda, nrand=nrand, default=NA)
        which_NA_simplelambda <- which(is.na(lFix) & 
                                         (is.na(optim_lambda_with_NAs) & ! is.nan(optim_lambda_with_NAs)) & ## explicit NaN's will be inner-optimized
                                         ! ranCoefs_blob$isRandomSlope) ## exclude random slope whether set or not
        if (length(which_NA_simplelambda)) { ## but not NaN
          init_lambda <- .eval_init_lambda_guess(processed, stillNAs=which_NA_simplelambda, For="optim")
          optim_lambda_with_NAs[which_NA_simplelambda] <- init_lambda[which_NA_simplelambda]
          init.optim$lambda <- optim_lambda_with_NAs[ ! is.na(optim_lambda_with_NAs)] ## but NaN still there 
        }
      } ## else use inner optimization  for simple lambdas if inner_phi is necessary
      init.optim$lambda <- init.optim$lambda[ ! is.nan(init.optim$lambda)] ## removes users's explicit NaN, which effect is documented in help(fitme)
      #
      if (any(var_ranCoefs)) {
        guess_from_glm_lambda <- .get_inits_by_glm(processed)$lambda * (3L*nranterms)/((nranterms+1L)) # +1 for residual
        fam_corrected_guess <- .calc_fam_corrected_guess(guess=guess_from_glm_lambda, For="optim", processed=processed) ## divides by nrand...
        for (rt in which(var_ranCoefs)) {
          char_rt <- as.character(rt)
          if (is.null(init.optim$ranCoefs[[char_rt]])) {
            Xi_cols <- attr(proc1$ZAlist,'Xi_cols')
            Xi_ncol <- Xi_cols[rt]
            rc <- rep(0,Xi_ncol*(Xi_ncol+1L)/2L)
            rc[cumsum(seq(Xi_ncol))] <- fam_corrected_guess/(Xi_ncol)
            init.optim$ranCoefs[[char_rt]] <- rc ## see help(ran)
          }
        }
      }
    }
  }
  family <- proc1$family
  if (family$family=="COMPoisson") {
    checknu <- suppressWarnings(try(environment(family$aic)$nu,silent=TRUE))
    if (inherits(checknu,"try-error") && is.null(init.optim$COMP_nu)) init.optim$COMP_nu <- 1 ## FR->FR FIXME what if not incuded in the range ?
  } else if (family$family == "negbin") {
    checktheta <- suppressWarnings(try(environment(family$aic)$shape,silent=TRUE))
    if (inherits(checktheta,"try-error") && is.null(init.optim$NB_shape)) init.optim$NB_shape <- 1 ## FR->FR FIXME idem ...
  }
  user.lower <- .post_process_parlist(lower,corr_types)
  user.upper <- .post_process_parlist(upper,corr_types) ## keep user input 
  .check_conflict_init_fixed(fixed,init.optim, "given as element of both 'fixed' and 'init'. Check call.")
  .check_conflict_init_fixed(init.HLfit,init.optim, "given as element of both 'init.HLfit' and 'init'. Check call.") ## has quite poor effect on fits
  moreargs <- .calc_moreargs(processed=processed, # possibly a list of environments -> .calc_range_info -> scans then to compute a mean(nbUnique) 
                             corr_types=corr_types, fixed=fixed, init.optim=init.optim, control_dist=processed$control_dist, 
                             init.HLfit=init.HLfit, corr_info=corr_info, verbose=verbose, lower=lower, upper=upper)
  inits <- .calc_inits(init.optim=init.optim, init.HLfit=init.HLfit,
                       ranFix=fixed, corr_types=corr_types,
                       moreargs=moreargs,
                       user.lower=user.lower, user.upper=user.upper,
                       optim.scale=optim.scale, 
                       For="fitme"
  )
  #
  init <- inits$`init` ## list; keeps all init values, all in untransformed scale
  init.optim <- inits$`init.optim` ## list; subset of all estimands, as name implies, and in transformed scale
  init.HLfit <- inits$`init.HLfit` ## list; subset as name implies 
  #
  if ("lambda" %in% c(names(user.lower),names(user.lower)) 
      && is.null(init$lambda)) {
    stop("'lambda' in 'lower' or 'upper' has no effect if absent from 'init'.")
  }
  ################
  LUarglist <- list(canon.init=init,
                    init.optim=init.optim,
                    user.lower=user.lower,user.upper=user.upper,
                    corr_types=corr_types,
                    ranFix=fixed,
                    optim.scale=optim.scale, 
                    moreargs=moreargs) ## list needed as part of attr(,"optimInfo")
  LowUp <- do.call(".makeLowerUpper",LUarglist)
  ## LowUp: a list with elements lower and upper that inherits names from init.optim, must be optim.scale as init.optim is by construction
  lower <- LowUp$lower ## list ! which elements may have length >1 !
  upper <- LowUp$upper ## list !
  #
  processedHL1 <- proc1$HL[1] ## there's also HLmethod in processed<[[]]>$callargs
  needHLCor_specific_args <- (length(unlist(lower$corrPars)) || length(intersect(corr_types,c("Matern","adjacency","AR1","corrMatrix"))))
  if (needHLCor_specific_args) {
    HLcallfn.obj <- "HLCor.obj" 
    HLcallfn <- "HLCor"
    control.dist <- vector("list",length(moreargs))
    for (nam in names(moreargs)) control.dist[[nam]] <- moreargs[[nam]]$control.dist
    HLCor.args$control.dist <- control.dist # lapply(moreargs,getElement,name="control.dist")
    HLCor.args$ranPars <- fixed ## to be modified by objective function
  } else {
    HLcallfn.obj <- "HLfit.obj"
    HLcallfn <- "HLfit"
    HLCor.args$ranFix <- fixed  
  }
  HLCor.args$init.HLfit <- init.HLfit
  HLCor.args$processed <- processed ## for the <...>.obj and <...>_body functions  
  ## 
  anyHLCor_obj_args <- HLCor.args
  ## HLCor.obj uses a vector + skeleton
  initvec <- unlist(init.optim)
  if (needHLCor_specific_args) {
    anyHLCor_obj_args$skeleton <- structure(init.optim, moreargs=moreargs,
                                            type=relist(rep("fix",length(initvec)),init.optim),
                                            moreargs=moreargs)
  } else {
    anyHLCor_obj_args$skeleton <- structure(init.optim,
                                            type=relist(rep("fix",length(initvec)),init.optim))
  }
  .assignWrapper(anyHLCor_obj_args$processed,
                   paste("return_only <- \"",anyHLCor_obj_args$processed$objective,"APHLs\"",sep=""))
  if (length(initvec)) {
    
    if (identical(verbose["getCall"][[1L]],TRUE)) { ## toget an optim call with its initial value. Then HLcallfn is called and its call returned.
      ## confint -> get_HLCorcall needs an HLCor call with the following ranFix
      ranPars_in_refit <- structure(.modify_list(fixed,init.optim),
                                   type=.modify_list(attr(fixed,"type"), 
                                                     relist(rep("outer",length(initvec)),init.optim)) )
    } else {
      use_SEM <- (!is.null(processedHL1) && processedHL1=="SEM")
      if (use_SEM) {
        if (is.null(proc1$SEMargs$control_pmvnorm$maxpts)) {
          if (length(LowUp$lower)) {
            .assignWrapper(processed,"SEMargs$control_pmvnorm$maxpts <- quote(250L*nobs)") 
          } ## else default visible in SEMbetalambda
        }
        stop("reimplement iterateSEMSmooth() in fitme() later")
        ## and then we need to put back the code for logL from smoothing at the end of this function
        ## beware of SEM case with length(lower)=0
        refit_info <- control[["refit"]] ## may be a NULL/ NA/ boolean or a list of booleans 
        if (is.null(refit_info) || is.na(refit_info)) refit_info <- FALSE ## alternatives are TRUE or an explicit list, but no longer NULL
      } else {
        optPars <- .new_locoptim(init.optim, 
                                 LowUp, 
                                 control, objfn_locoptim=.objfn_locoptim, 
                                 HLcallfn.obj=HLcallfn.obj, anyHLCor_obj_args=anyHLCor_obj_args, 
                                 user_init_optim=user_init_optim)
        refit_info <- attr(optPars,"refit_info") ## 'derives' from control[["refit"]]
      }
      ranPars_in_refit <- structure(.modify_list(fixed,optPars),
                                   type=.modify_list(attr(fixed,"type"),
                                                     relist(rep("outer",length(unlist(optPars))),optPars)))
      if ( is.list(refit_info)) {refit_phi <- refit_info$phi} else refit_phi <- refit_info ## result may be NULL in list case
      init_refit <- list()
      if (identical(refit_phi,TRUE) && ! is.null(optPars$trPhi)) {
        init_refit$phi <- .dispInv(optPars$trPhi)
        ranPars_in_refit$trPhi <- NULL
        attr(ranPars_in_refit,"type")$trPhi <- NULL 
      }
      if ( is.list(refit_info)) {refit_lambda <- refit_info$lambda} else refit_lambda <- refit_info
      if (identical(refit_lambda,TRUE) && ! is.null(optPars$trLambda)) {
        lambdapos <- as.integer(names(optPars$trLambda))
        lambda <- rep(NA,max(lambdapos))
        lambda[lambdapos] <- .dispInv(optPars$trLambda)
        init_refit$lambda <- lambda
        ranPars_in_refit$trLambda <- NULL
        attr(ranPars_in_refit,"type")$trLambda <- NULL 
      }
      if (length(init_refit)) HLCor.args$init.HLfit <- .modify_list(HLCor.args$init.HLfit, init_refit)  
    } ## end if ...getCall... else
    #
    # refit_info is list if so provided by user, else typically boolean. An input NA should have been converted to something else (not documented).
    if (needHLCor_specific_args) {
      attr(ranPars_in_refit,"moreargs") <- moreargs 
      HLCor.args$ranPars <- ranPars_in_refit ## variable locally
    } else HLCor.args$ranFix <- ranPars_in_refit ## variable locally  
    
  }
  #
  # if processed is an envir, the following is not local to anyHLCor_obj_args$processed but change processed globally
  .assignWrapper(HLCor.args$processed,"return_only <- NULL") 
  .assignWrapper(HLCor.args$processed,"verbose['warn'] <- TRUE") ## important!
  hlcor <- do.call(HLcallfn,HLCor.args) ## recomputation post optimization, or only computation if length(initvec)=0
  if (is.call(hlcor)) {
    ## then do.call(HLcallfn,HLCor.args) has retuned the call, not the fit. 
    ## see def of get_HLCorcall() for further explanation
    return(hlcor) ## HLCorcall
  }
  # hlcor may have received attr(.,"info.uniqueGeo") from HLCor_body.
  if (length(initvec)) {
    attr(hlcor,"optimInfo") <- list(LUarglist=LUarglist, optim.pars=optPars, 
                                    objective=anyHLCor_obj_args$processed$objective) ## processed was erased for safety
  }
  ## substantial effect on object size! :
  lsv <- c("lsv",ls())
  rm(list=setdiff(lsv,"hlcor")) 
  ##
  return(hlcor)
}