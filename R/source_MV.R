## The file is in .Rbuildignore

.fitmv <- function(mv,...) {
  assign("spaMM_glm_conv_crit",list(max=-Inf) , envir=environment(spaMM_glm.fit))
  time1 <- Sys.time()
  oricall <- match.call(expand.dots=TRUE) ## mc including dotlist
  n_models <- length(mv) 
  processed_calls <- vector("list",n_models)
  for (mit in seq_along(processed_calls)) {
    coll <- oricall
    coll["mv"] <- NULL
    # I need to match the names of mv[[mit]] to those of a fitme call to make sure that they all named...
    for (st in setdiff(names(mv[[mit]]),"")) coll[[st]] <- mv[[mit]][[st]]
    coll$formula <- .preprocess_formula(coll$formula) 
    coll[[1L]] <- get(".check_args_fitme", asNamespace("spaMM")) 
    coll <- eval(coll,parent.frame()) # 
    coll[[1L]] <- get(".preprocess_fitme", asNamespace("spaMM")) 
    processed_calls[[mit]] <- eval(coll,parent.frame()) # returns modified call including an element 'processed'
  }
  { ## then I need to merge the processed_calls in some way
    merged <- .merge_processed(processed_calls) # assumes non-trivial 'mv'
  }
  # merged[[1L]] <- get("fitme_body", asNamespace("spaMM")) 
  # hlcor <- eval(merged,parent.frame()) 
  hlcor <- .fitmv_body(processed=merged)
  #
  for (mit in seq_along(processed_calls)) {
    if ("control.dist" %in% names(mv[[mit]])) {
      oricall[["mv"]][["control.dist"]] <- processed_calls[[mit]][["control_dist"]]
    } else oricall$"control.dist" <- processed_calls[[mit]][["control_dist"]] ## maybe # [[]] <- does not work if [["control_dist"]] orginally absent
    hlcor$call <- oricall ## this is a call to fitmv()
  }
  lsv <- c("lsv",ls())
  if ( ! .is.multi(family) && ! is.call(hlcor) ) {
    hlcor$how$fit_time <- .timerraw(time1)
    hlcor$how$fnname <- ".fitmv"
    hlcor$fit_time <- structure(hlcor$how$fit_time,
                                message="Please use how(<fit object>)[['fit_time']] to extract this information cleanly.")
  }
  rm(list=setdiff(lsv,"hlcor")) ## empties the whole local envir except the return value
  return(hlcor)
}

.merge_processed <- function(processed_calls) {
  nmodels <- length(processed_calls)
  unmerged <- vector("list",nmodels)
  for (mit in seq_len(nmodels)) {
    unmerged[[mit]] <- processed_calls[[mit]][["processed"]]
    #names_list[[mit]] <- names(unmerged[[mit]])
  }
  # allnames <- unique(unlist(names_list))
  # for (st in allnames) if (!identical(processed_calls[[1L]][["processed"]][[st]],processed_calls[[2L]][["processed"]][[st]])) print(st)
  vec_nobs <- integer(nmodels)
  merged <- list2env(list(ZAlist= unmerged[[1L]][["ZAlist"]],
                          AUGI0_ZX=unmerged[[1L]][["AUGI0_ZX"]],
                          REMLformula=unmerged[[1L]][["REMLformula"]],
                          X.Re=unmerged[[1L]][["X.Re"]],
                          verbose=unmerged[[1L]][["verbose"]], # needed for assignment to verbose['warn'] in fitmv_body()
                          family=unmerged[[1L]][["family"]], # FIXME quick patch
                          models=c(eta="etaHGLM",phi=unmerged[[1L]][["models"]][["phi"]]),
                          control.glm=unmerged[[1L]][["control.glm"]],
                          HL=unmerged[[1L]][["HL"]],
                          rand.families=unmerged[[1L]][["rand.families"]],
                          LevenbergM=unmerged[[1L]][["LevenbergM"]],
                          spaMM_tol=unmerged[[1L]][["spaMM_tol"]],
                          clik_fn=unmerged[[1L]][["clik_fn"]], ## FIXME... def'ly wrong when distinct families
                          residModel=unmerged[[1L]][["residModel"]], ## FIXME...
                          break_conv_logL=unmerged[[1L]][["break_conv_logL"]]
  ))
  prior.weights <- unmerged[[1L]][["prior.weights"]]
  off <- unmerged[[1L]][["off"]]
  y <- unmerged[[1L]][["y"]]
  vec_nobs[1L] <- length(y)  
  BinomialDen <- unmerged[[1L]][["BinomialDen"]]
  merged_X <- .unscale(unmerged[[1L]][["AUGI0_ZX"]]$X.pv)
  HLframes <- unmerged[[1L]][["HLframes"]]
  LMMbool <- unmerged[[1L]][["LMMbool"]]
  GLMMbool <- unmerged[[1L]][["GLMMbool"]]
  LLM_const_w <- unmerged[[1L]][["LLM_const_w"]]
  GLGLLM_const_w  <- unmerged[[1L]][["GLGLLM_const_w"]]
  iter_mean_dispVar <- unmerged[[1L]][["iter_mean_dispVar"]]
  canonicalLink <- unmerged[[1L]][["canonicalLink"]]
  max.iter <- unmerged[[1L]][["max.iter"]]
  for (proc_it in (seq_len(nmodels-1L)+1L)) {
    p2 <- unmerged[[proc_it]]
    # merged$predictor... I should try to get rid of its use in HLfit_body, at least... FIXME
    ### random effects stuff
    merged$ZAlist <- .merge_ZAlists(merged,p2)
    # rand.families, lcrandfamfam... FIXME
    ### response stuff:
    # family FIXME
    vec_nobs[proc_it] <- length(p2[["y"]])  
    y <- c(y,p2[["y"]])
    BinomialDen <- c(BinomialDen,p2[["BinomialDen"]])
    prior.weights <- c(prior.weights, p2[["prior.weights"]])
    off <- c(off, p2[["off"]])
    X2 <- .unscale(p2$AUGI0_ZX$X.pv)
    merged_X <- .merge_Xs(merged_X,X2)  
    HLframes[["Y"]]  <- rbind(HLframes[["Y"]], p2$HLframes[["Y"]]) # FIXME that for glm code... troubles expected
    LMMbool <- LMMbool && p2[["LMMbool"]]
    GLMMbool <- GLMMbool && p2[["GLMMbool"]]
    LLM_const_w <- LLM_const_w && p2[["LLM_const_w"]]
    GLGLLM_const_w <- GLGLLM_const_w && p2[["GLGLLM_const_w"]]
    canonicalLink <- canonicalLink && p2[["canonicalLink"]]
    iter_mean_dispVar <- max(iter_mean_dispVar, p2[["iter_mean_dispVar"]])
    max.iter <- max(max.iter, p2[["max.iter"]])
    # X.Re...
  }
  merged[["y"]] <- y
  merged[["BinomialDen"]] <- BinomialDen
  merged[["HLframes"]] <- structure(HLframes, vec_nobs=vec_nobs) 
  # FIXME next line probably note valid for quoted expressions (but .preprocess_pw() not usable here) 
  merged[["prior.weights"]] <- structure(prior.weights, unique=length(unique(prior.weights))==1L)
  merged[["off"]] <- off
  merged[["LMMbool"]] <- LMMbool
  merged[["GLMMbool"]] <- GLMMbool
  merged[["LLM_const_w"]] <- LLM_const_w
  merged[["GLGLLM_const_w"]] <- GLGLLM_const_w 
  merged[["canonicalLink"]] <- canonicalLink 
  merged[["iter_mean_dispVar"]] <- iter_mean_dispVar
  merged[["max.iter"]] <- max.iter
  # "method" should not be allowed in the 'mv' argument FIXME
  if (.spaMM.data$options$X_scaling) { ## use scaled X.pv by default v.2.4.83
    ##       standard REML    ||      ML 
    if ( is.null(merged$REMLformula) || ncol(merged$X.Re)==0L) merged_X <- .scale(merged_X) # scales and adds "scaled:scale" attribute
  }
  {
    # a true HLframes object would be needed for this to work: (in part to determine sparse_X)
    #X.pv <- .post_process_X(merged_X, merged$HL, merged$HLframes, control.HLfit=NULL) 
    #
    if (ncol(merged_X)) {
      Xattr <- attributes(merged_X)
      rankinfo <- .calc_rankinfo(merged_X, tol=spaMM.getOption("rankTolerance"))
      if (is.list(rankinfo)) merged_X <- .rankTrim(merged_X,rankinfo = rankinfo) # adds important attributes
      names_lostattrs <- setdiff(names(Xattr), c(names(attributes(merged_X)),"dim","dimnames"))
      attributes(merged_X)[names_lostattrs] <- Xattr[names_lostattrs] # as in .subcol_wAttr(). 
    }
  }
  # LMMbool..., GLMMbool...
  vec_n_u_h <- lapply(merged$ZAlist,ncol)
  merged[["cum_n_u_h"]] <- cumsum(c(0L, vec_n_u_h))
  nrd <- tail(merged[["cum_n_u_h"]],n=1)
  if ( ! .eval_as_mat_arg(merged)) {
    AUGI0_ZX <- list2env( list(I=.trDiagonal(n=nrd),
                               ZeroBlock= Matrix(0,nrow=nrd,ncol=ncol(merged_X)), X.pv=merged_X) )
  } else {
    AUGI0_ZX <- list2env( list(I=diag(nrow=nrd),ZeroBlock= matrix(0,nrow=nrd,ncol=ncol(merged_X)), 
                               X.pv=merged_X) )
  }
  #
  ##  control.HLfit$rankinfo cannot be used (not a big deal) unless 
  ## * it is procdied as an argument (it's not in processed)
  ## * something specific is done to merge different control.HLfit's... 
  merged$is_spprec <- unmerged[[1L]]$is_spprec # quick patch
  AUGI0_ZX <- .add_ZAfix_info(AUGI0_ZX, ZAlist=merged$ZAlist, sparse_precision=merged$is_spprec, 
                              as_mat=.eval_as_mat_arg(unmerged[[1L]])) # as_mat: quick patch
  #
  if (length(merged$ZAList) &&  ! merged$is_spprec ) { 
    if (inherits(AUGI0_ZX$ZeroBlock,"sparseMatrix")) {
      AUGI0_ZX$Zero_sparseX <- rbind2(AUGI0_ZX$ZeroBlock, as(AUGI0_ZX$X.pv,"CsparseMatrix"))
    } else AUGI0_ZX$Zero_sparseX <- rbind2(AUGI0_ZX$ZeroBlock, AUGI0_ZX$X.pv)
  }
  #
  merged$AUGI0_ZX <- AUGI0_ZX
  merged$vec_nobs <- vec_nobs
  merged$unmerged <- unmerged # A list whose each elemen is the processed arg of each processed call
  class(merged) <- c("mvarglist", class(merged))
  return(merged)
}

.merge_ZAlists <- function(p1,p2) {
  ZAlist1 <- p1$ZAlist
  ZAlist2 <- p2$ZAlist
  ranefs1 <- attr(ZAlist1,"exp_ranef_strings")
  ranefs2 <- attr(ZAlist2,"exp_ranef_strings")
  allranefs <- unique(c(ranefs1,ranefs2))
  ZAlist <- exp_ranef_terms <- namesTerms <- vector("list", length(allranefs))
  exp_ranef_types <- exp_ranef_strings <- type_attr <- character(length(allranefs))
  Xi_cols <- integer(length(allranefs))
  # exp_ranef_strings and exp_ranef_terms will have a global 'type' attibute
  for (ran_it in seq_along(allranefs)) {
    ranf <- allranefs[ran_it]
    in1 <- which(ranefs1==ranf)
    in2 <- which(ranefs2==ranf)
    if (length(in1) && length(in2)) { # true merge
      Zori <- ZAlist1
      ZA1 <- ZAlist1[[in1]]
      ZA2 <- ZAlist2[[in2]]
      # do I have levels info elsewhere ?
      levels1 <- colnames(ZA1)
      levels2 <- colnames(ZA2)
      if ( ! identical(levels1,levels2)) {
        alllevels <- unique(c(levels1,levels2))
        extracols <- setdiff(alllevels,levels1)
        nextcolnames <- c(levels1,extracols)
        ZA1 <- .adhoc_cbind_dgC_0(ZA1,length(extracols))
        colnames(ZA1) <- nextcolnames
        ZA1 <- ZA1[,alllevels]
        #
        extracols <- setdiff(alllevels,levels2)
        nextcolnames <- c(levels2,extracols)
        ZA2 <- .adhoc_cbind_dgC_0(ZA2,length(extracols))
        colnames(ZA2) <- nextcolnames
        ZA2 <- ZA2[,alllevels]
      }
      ZAlist[[ran_it]] <- rbind(ZA1,ZA2)
      ori <- in1
    } else  {
      if (length(in1)) {
        Zori <- ZAlist1
        ori <- in1
      } else {
        Zori <- ZAlist2
        ori <- in2
      }
      ZAlist[[ran_it]] <- Zori[[ori]]
    } 
    namesTerms[ran_it] <- attr(Zori,"namesTerms")[ori] # namesTerms is a *named* list... but pathetically that does not copy the name
    names(namesTerms)[ran_it] <- names(attr(Zori,"namesTerms")[ori])
    type_attr[ran_it] <- attr(attr(Zori,"exp_ranef_terms"),"type")[ori]
    exp_ranef_terms[[ran_it]] <- attr(Zori,"exp_ranef_terms")[[ori]]
    exp_ranef_types[ran_it] <- attr(Zori,"exp_ranef_types")[ori]
    exp_ranef_strings[ran_it] <- attr(Zori,"exp_ranef_strings")[ori]
    Xi_cols[ran_it] <- attr(Zori,"Xi_cols")[ori]
  }
  return(structure(ZAlist, 
                   namesTerms=namesTerms, 
                   exp_ranef_terms=structure(exp_ranef_terms,type=type_attr),
                   exp_ranef_types=exp_ranef_types,
                   exp_ranef_strings=structure(exp_ranef_strings,type=type_attr),
                   Xi_cols=Xi_cols
  ))
} 

.merge_Xs <- function(X1,X2) {
  if (attr(X1,"rankinfo")$rank < ncol(X1) || attr(X2,"rankinfo")$rank < ncol(X2)) {
    stop("case not yet handled")
    # for 'mv' case, I should inhibit .preprocess from .post_process_X()ing
  } 
  regr1 <- colnames(X1)
  regr2 <- colnames(X2)
  if ( ! identical(regr1,regr2)) {
    allregr <- unique(c(regr1,regr2))
    extracols <- setdiff(allregr,regr1)
    nextcolnames <- c(regr1,extracols)
    X1 <- cbind(X1,matrix(0, nrow=nrow(X1), ncol=length(extracols)))
    colnames(X1) <- nextcolnames
    X1 <- X1[,allregr, drop=FALSE]
    #
    extracols <- setdiff(allregr,regr2)
    nextcolnames <- c(regr2,extracols)
    X2 <- cbind(X2,matrix(0, nrow=nrow(X2), ncol=length(extracols)))
    colnames(X2) <- nextcolnames
    X2 <- X2[,allregr, drop=FALSE]
  }
  X <- rbind(X1,X2)
  return(X)
} 

.merge_optim_blobs <- function(optim_blobs) {
  optim_blob <- optim_blobs[[1L]] # FIXME quick patch 
  for (it in seq_along(optim_blobs)) {
    
  }
  return(optim_blob)
}

.fitmv_body <- function(processed,
                       init=list(),
                       #                       init.HLfit=list(),
                       fixed=list(), ## NULL not valid (should be handled in preprocess?)
                       lower=list(),upper=list(),
                       # control.dist=list(), ## info in processed
                       control=list(), ## optimizer, <optimizer controls>, precision
                       nb_cores=NULL,
                       ... ## cf dotnames processing below 
) {
  dotlist <- list(...) ## forces evaluations, which makes programming easier...
  if (is.list(processed)) {
    proc1 <- processed[[1L]]
  } else proc1 <- processed
  verbose <-  proc1$verbose
  HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                names(formals(mat_sqrt)),names(formals(make_scaled_dist))))  ## cf parallel code in HLCor.obj
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
  sparse_precision <- proc1$is_spprec
  #
  user_init_optim <- init ## user_init_optim only for a check in new_locoptim, true initial value init.optim is modified below
  
  unmerged <- processed$unmerged
  optim_blobs <- vector("list",length(unmerged))
  for (it in seq_along(unmerged)) optim_blobs[[it]] <-.calc_optim_args(proc1=unmerged[[it]], processed=unmerged[[it]],
                                                                       init=user_init_optim, fixed=fixed, lower=lower, upper=upper, 
                                                                       verbose=verbose, optim.scale=optim.scale, For="fitme") 
  optim_blob <- .merge_optim_blobs(optim_blobs)
  
  # modify HLCor.args and <>bounds;   ## distMatrix or uniqueGeo potentially added to HLCor.args:
  # init <- optim_blob$inits$`init` ## list; keeps all init values, all in untransformed scale
  init.optim <- optim_blob$inits$`init.optim` ## list; subset of all estimands, as name implies, and in transformed scale
  init.HLfit <- optim_blob$inits$`init.HLfit` ## list; subset as name implies 
  fixed <- optim_blob$fixed
  corr_types <- optim_blob$corr_types
  moreargs <- optim_blob$moreargs
  LUarglist <- optim_blob$LUarglist
  LowUp <- optim_blob$LowUp
  lower <- LowUp$lower ## list ! which elements may have length >1 !
  upper <- LowUp$upper ## list !
  #
  
  #
  if ( ! is.null(residproc1 <- proc1$residProcessed) && identical(spaMM.getOption("outer_optim_resid"),TRUE)) {
    ## Problem is that outer optim at the mean model level is useful if we can avoid computation of the leverages 
    ## But here anyway we need the leverages of the 'mean' model to define the resid model response.
    resid_optim_blob <- .calc_optim_args(proc1=residproc1, processed=proc1,
                                         init=proc1$residModel$init, fixed=proc1$residModel$fixed, ## all user input must be in proc1$residModel
                                         lower=proc1$residModel$lower, upper=proc1$residModel$upper, ## all user input must be in proc1$residModel
                                         verbose=c(SEM=FALSE), optim.scale=optim.scale, For="fitme") 
    resid_init.optim <- resid_optim_blob$inits$`init.optim` ## list; subset of all estimands, as name implies, and in transformed scale
    proc1$residModel$`init.HLfit` <- resid_optim_blob$inits$`init.HLfit` ## list; subset as name implies 
    proc1$residModel$fixed <- resid_optim_blob$fixed
    ##### residproc1$corr_types <- resid_optim_blob$corr_types
    # resid_user.lower <- resid_optim_blob$user.lower
    # resid_user.upper <- resid_optim_blob$user.upper
    init.optim <- c(init.optim,list(resid=resid_init.optim))
    lower <- c(lower,list(resid=resid_optim_blob$LowUp$lower))
    upper <- c(upper,list(resid=resid_optim_blob$LowUp$upper))
    LowUp <- list(lower=lower,upper=upper)
    moreargs <- c(moreargs,list(resid=resid_optim_blob$moreargs))
  }
  
  
  processedHL1 <- proc1$HL[1] ## there's also HLmethod in processed<[[]]>$callargs
  needHLCor_specific_args <- (length(unlist(lower$corrPars)) || 
                                length(intersect(corr_types,c("Matern","Cauchy","adjacency","AR1","corrMatrix", "IMRF"))))
  if (needHLCor_specific_args) {
    HLcallfn.obj <- "HLCor.obj" 
    HLcallfn <- "HLCor"
    control.dist <- vector("list",length(moreargs))
    for (nam in names(moreargs)) control.dist[[nam]] <- moreargs[[nam]]$control.dist 
    HLCor.args[["control.dist"]] <- control.dist ## always reconstructed locally, not in the fitme_body call
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
    anyHLCor_obj_args$skeleton <- structure(init.optim, moreargs=moreargs, ## moreargs is a list over ranefs 
                                            type=relist(rep("fix",length(initvec)),init.optim) )
  } else {
    anyHLCor_obj_args$skeleton <- structure(init.optim,
                                            type=relist(rep("fix",length(initvec)),init.optim))
  }
  .assignWrapper(anyHLCor_obj_args$processed,
                 paste0("return_only <- \"",proc1$objective,"APHLs\""))
  if (length(initvec)) {
    augZXy_phi_est <- NULL
    if (identical(verbose["getCall"][[1L]],TRUE)) { ## toget an optim call with its initial value. Then HLcallfn is called and its call returned.
      ## confint -> get_HLCorcall needs an HLCor call with the following ranFix
      ranPars_in_refit <- structure(.modify_list(fixed,init.optim),
                                    # I label this "F I X" as a TAG for this modif type attribute:
                                    type=.modify_list(relist(rep("fix",length(unlist(fixed))),fixed), #attr(fixed,"type"), 
                                                      relist(rep("outer",length(initvec)),init.optim)) )
    } else {
      use_SEM <- (!is.null(processedHL1) && processedHL1=="SEM")
      time2 <- Sys.time()
      if (use_SEM) {
        optimMethod <- "iterateSEMSmooth"
        if (is.null(proc1$SEMargs$control_pmvnorm$maxpts)) {
          .assignWrapper(processed,"SEMargs$control_pmvnorm$maxpts <- quote(250L*nobs)") 
        } ## else default visible in SEMbetalambda
        ## its names should match the colnames of the data in Krigobj = the  parameters of the likelihood surface. Current code maybe not general.
        loclist <- list(anyHLCor_obj_args=anyHLCor_obj_args,  ## contains $processed
                        LowUp=LowUp,init.corrHLfit=init.optim, ## F I X M E usage of user_init_optim probably not definitive
                        control.corrHLfit=control,
                        verbose=verbose[["iterateSEM"]],
                        nb_cores=nb_cores)
        optr <- .probitgemWrap("iterateSEMSmooth",arglist=loclist, pack="probitgem") # do.call("iterateSEMSmooth",loclist) 
        optPars <- relist(optr$par,init.optim)
        if (!is.null(optPars)) attr(optPars,"method") <-"optimthroughSmooth"
        refit_info <- control[["refit"]] ## may be a NULL/ NA/ boolean or a list of booleans 
        if ( is.null(refit_info) || is.na(refit_info)) refit_info <- FALSE ## alternatives are TRUE or an explicit list or NULL
      } else {
        optPars <- .new_locoptim(init.optim, ## try to use gradient? But neither minqa nor _LN_BOBYQA use gradients. optim() can
                                 LowUp, 
                                 control, objfn_locoptim=.objfn_locoptim, 
                                 HLcallfn.obj=HLcallfn.obj, anyHLCor_obj_args=anyHLCor_obj_args, 
                                 user_init_optim=user_init_optim,
                                 grad_locoptim=NULL, verbose=verbose[["TRACE"]])
        refit_info <- attr(optPars,"refit_info") ## 'derives' from control[["refit"]] with modification ## may be NULL but not NA
      }
      optim_time <- .timerraw(time2)
      ranPars_in_refit <- structure(.modify_list(fixed,optPars),
                                    type=.modify_list(relist(rep("fix",length(unlist(fixed))),fixed), #attr(fixed,"type"),
                                                      relist(rep("outer",length(unlist(optPars))),optPars)))
      ranPars_in_refit <- .expand_hyper(ranPars_in_refit, processed$hyper_info,moreargs=moreargs)
      augZXy_phi_est <- proc1$augZXy_env$phi_est ## may be NULL: if phi was not estimated by augZXy
      
      #any_nearly_singular_covmat <- FALSE
      if (! is.null(optPars$trRanCoefs)) {
        ranCoefs <- optPars$trRanCoefs # copy names...
        for (char_rd in names(optPars$trRanCoefs)) {
          ranCoefs[[char_rd]] <- .ranCoefsInv(optPars$trRanCoefs[[char_rd]], rC_transf=.spaMM.data$options$rC_transf)
          covmat <- .calc_cov_from_ranCoef(ranCoefs[[char_rd]])
          #if (kappa(covmat)>1e14 || min(eigen(covmat,only.values = TRUE)$values)<1e-6) any_nearly_singular_covmat <- TRUE # use of this removed 2019/12/16
        }
      }
      init_refit <- list()
      if ( is.null(refit_info) ) { ## not the case with SEM
        refit_info <- list(phi=FALSE,lambda=(! is.null(optPars$trRanCoefs)),ranCoefs=FALSE)
      } else if ( is.list(refit_info) ) { ## never the default
        refit_info <- .modify_list( list(phi=FALSE,lambda=(! is.null(optPars$trRanCoefs)),
                                         ranCoefs=FALSE), refit_info) 
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
        ranges <- attr(optPars$hyper,"ranges")
        for (rg_it in names(ranges)) {
          hyper_el <- ranPars_in_refit$hyper[[rg_it]]
          if ( ! is.null(hy_trL <- hyper_el$hy_trL)) {
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
            rancoef <- ranCoefs[[char_rd]]
            if ( ! is.null(augZXy_phi_est)) {
              Xi_ncol <- attr(rancoef,"Xi_ncol")
              lampos <- rev(length(rancoef) -cumsum(seq(Xi_ncol))+1L)  ## NOT cumsum(seq(Xi_cols))
              rancoef[lampos] <- rancoef[lampos] *augZXy_phi_est
              rancoef <- as.vector(rancoef) ## as.vector drops attributes
            } 
            if (identical(refit_info$ranCoefs,TRUE)) { 
              init_refit$ranCoefs[[char_rd]] <- rancoef
              ranPars_in_refit$trRanCoefs[char_rd] <- NULL
              attr(ranPars_in_refit,"type")$trRanCoefs[char_rd] <- NULL
            } else ranPars_in_refit$trRanCoefs[[char_rd]] <- .ranCoefsFn(rancoef, rC_transf=.spaMM.data$options$rC_transf_inner)  
          }
        }      
      }
      if (length(init_refit)) HLCor.args$init.HLfit <- .modify_list(HLCor.args$init.HLfit, init_refit)  
    } ## end if ...getCall... else
    #
    # refit_info is list if so provided by user, else typically boolean. An input NA should have been converted to something else (not documented).
    if (needHLCor_specific_args) {
      attr(ranPars_in_refit,"moreargs") <- moreargs 
      HLCor.args$ranPars <- ranPars_in_refit 
    } else HLCor.args$ranFix <- ranPars_in_refit 
  } else if (len_ranPars <- length(unlist(HLCor.args$ranPars))){ ## Set attribute
    HLCor.args$ranPars <- structure(HLCor.args$ranPars,
                                    type = relist(rep("fix", len_ranPars), HLCor.args$ranPars),
                                    moreargs=moreargs) ## moreargs needed if user handles fixed(<transformed params>) ('hyper' tests)
  }
  #
  # not local to anyHLCor_obj_args$processed: change processed globally
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
                                    objective=proc1$objective,
                                    augZXy_phi_est=augZXy_phi_est, ## gives info to interpret optim.pars in confint.HLfit()
                                    optim_time=optim_time) ## processed was erased for safety
    locoptr <- attr(optPars,"optr")
    if (attr(optPars,"method")=="nloptr") {
      if (locoptr$status<0L) hlcor$warnings$optimMessage <- paste0("nloptr() message: ",
                                                                   locoptr$message," (status=",locoptr$status,")")
    } else if ( attr(optPars,"method")=="optim" ) {
      if (locoptr$convergence) hlcor$warnings$optimMessage <- paste0("optim() message: ",locoptr$message,
                                                                     " (convergence=",locoptr$convergence,")")
    } else if ( attr(optPars,"method")== "optimthroughSmooth") {
      # provide logL estimate from the smoothing, to be used rather than the hlcor logL :
      logLapp <- optr$value
      attr(logLapp,"method") <- "  logL (smoothed)" 
      hlcor$APHLs$logLapp <- logLapp
    }
    if ( ! is.null(PQLdivinfo <- processed$envir$PQLdivinfo)) {
      hlcor$divinfo <- PQLdivinfo
      hlcor$warnings$divinfo <- "Numerical issue detected; see div_info(<fit object>) for more information."
      warning(hlcor$warnings$divinfo)
    }
  }
  ## substantial effect on object size! :
  lsv <- c("lsv",ls())
  rm(list=setdiff(lsv,"hlcor")) 
  ##
  return(hlcor)
}