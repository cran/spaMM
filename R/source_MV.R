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
  hlcor <- fitme_body(processed=merged)
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
  names_list <- HLframes <- vector("list",nmodels)
  for (mit in seq_len(nmodels)) names_list[[mit]] <- names(processed_calls[[mit]][["processed"]])
  allnames <- unique(unlist(names_list))
  # for (st in allnames) if (!identical(processed_calls[[1L]][["processed"]][[st]],processed_calls[[2L]][["processed"]][[st]])) print(st)
  vec_nobs <- integer(nmodels)
  merged <- processed_calls[[1L]][["processed"]]
  prior.weights <- merged[["prior.weights"]]
  vec_nobs[1L] <- length(merged[["y"]])  
  merged_X <- .unscale(merged$AUGI0_ZX$X.pv)
  HLframes[[1]]  <- merged$HLframes
  for (proc_it in (seq_len(nmodels-1L)+1L)) {
    # merged$predictor... I should try to get rid of its use in HLfit_body, at least... FIXME
    p2 <- processed_calls[[proc_it]][["processed"]]
    ### random effects stuff
    merged$ZAlist <- .merge_ZAlists(merged,p2)
    # rand.families, lcrandfamfam... FIXME
    ### response stuff:
    # family FIXME
    vec_nobs[proc_it] <- length(p2[["y"]])  
    merged[["y"]] <- c(merged[["y"]],p2[["y"]])
    merged_Y <- rbind
    merged[["BinomialDen"]] <- c(merged[["BinomialDen"]],p2[["BinomialDen"]])
    prior.weights <- c(prior.weights, p2[["prior.weights"]])
    X2 <- .unscale(p2$AUGI0_ZX$X.pv)
    merged_X <- .merge_Xs(merged_X,X2)  
    HLframes[[proc_it]]  <- p2$HLframes
    # X.Re...
  }
  merged$HLframes <- structure(HLframes, vec_nobs=vec_nobs)
  # FIXME next line probably note valid for quoted expressions (but .preprocess_pw() not usable here) 
  merged[["prior.weights"]] <- structure(prior.weights, unique=length(unique(prior.weights))==1L)
  # "method" should not be allowed in the 'mv' argument FIXME
  if (.spaMM.data$options$X_scaling) { ## use scaled X.pv by default v.2.4.83
    ##       standard REML    ||      ML 
    if ( is.null(merged$REMLformula) || ncol(merged$X.Re)==0L) merged_X <- .scale(merged_X) # scales and adds "scaled:scale" attribute
  }
  {
    # a true HLframes object would be needed for this to work: (in part to determine sparse_X)
    #X.pv <- .post_process_X(merged_X, merged$HL, merged$HLframes, control.HLfit=NULL) 
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
  AUGI0_ZX <- .add_ZAfix_info(AUGI0_ZX, ZAlist=merged$ZAlist, sparse_precision=merged$is_spprec, processed=merged)
  #
  if (length(merged$ZAList) &&  ! merged$is_spprec ) { 
    if (inherits(AUGI0_ZX$ZeroBlock,"sparseMatrix")) {
      AUGI0_ZX$Zero_sparseX <- rbind2(AUGI0_ZX$ZeroBlock, as(AUGI0_ZX$X.pv,"CsparseMatrix"))
    } else AUGI0_ZX$Zero_sparseX <- rbind2(AUGI0_ZX$ZeroBlock, AUGI0_ZX$X.pv)
  }
  #
  merged$AUGI0_ZX <- AUGI0_ZX
  merged$vec_nobs <- vec_nobs
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


