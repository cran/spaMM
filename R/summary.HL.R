## FR->FR accessors : http://glmm.wikidot.com/pkg-comparison

.get_methods_disp <- function(object) {
  iterativeEsts<-character(0)
  optimEsts<-character(0)
  ## What the HLfit call says:
  if ( object$models[["eta"]] != "etaGLM") { 
    if ( any(object$lambda.object$type=="outer")) optimEsts <- c(optimEsts,"lambda")
    if ( any(object$lambda.object$type=="inner")) iterativeEsts <- c(iterativeEsts,"lambda")
  }
  if ( ! (object$family$family %in% c("binomial","poisson","COMPoisson","negbin"))) {
    if ( ! is.null(object$phi.object$fixef) ) {
      iterativeEsts <- c(iterativeEsts,"phi")
    } else if ( identical(attr(object$phi.object$phi_outer,"type"),"var")) optimEsts <- c(optimEsts,"phi")
    # confint (fitme HLfitoide -> 3e cas)
  }
  optimNames <- names(attr(object,"optimInfo")$LUarglist$canon.init) ## will contain eg NB_shape
  corrPars <- get_ranPars(object,which="corrPars")
  optimEsts <- c(optimEsts,intersect(optimNames,c(names(corrPars),"COMP_nu","NB_shape")))
  iterativeEsts <- c(iterativeEsts,setdiff(optimNames,optimEsts))
  return(list(iterativeEsts=iterativeEsts,optimEsts=optimEsts))
}

.make_beta_table <- function(object, p_value="") {
  namesOri <- attr(object$X.pv,"namesOri")
  nc <- length(namesOri)
  betaOri_cov <- matrix(NA,ncol=nc,nrow=nc,dimnames=list(rownames=namesOri,colnames=namesOri))
  beta_cov <- .get_beta_cov_any_version(object) ## typically set by HLfit using get_from_MME(, which="beta_cov"), circumventing full beta_w_cov computation
  betaOri_cov[colnames(beta_cov),colnames(beta_cov)] <- beta_cov
  beta_se <- sqrt(diag(betaOri_cov))
  fixef_z <- object$fixef/beta_se
  beta_table <- cbind(Estimate=object$fixef,"Cond. SE"=beta_se,"t-value"=fixef_z)
  if (p_value=="Wald") {beta_table <- cbind(beta_table,"p-value"=1-pchisq(fixef_z^2,df=1))} ## FIXEM F_test: need to extract dfs
  rownames(beta_table) <- names(object$fixef)
  return(beta_table)
}


.legend_lambda <- function(object, in_table, type="") {
  ## (1) analyse distributions
  randfams <- object$rand.families
  rff <- sapply(randfams, getElement, name="family")
  # info about links:
  rfl <- sapply(rff, switch, 
    #        "Beta" = "log(lambda)", ## u ~ Beta(1/(2lambda),1/(2lambda))
    #        "inverse.Gamma" = "log(lambda)", ## u ~ I-Gamma(1+1/lambda,1/lambda)
    #        "log(lambda)" ## gaussian, or gamma with u ~ Gamma(lambda,1/lambda)
    "log(lambda)" ## currently constant except as modified below for adjd
  )
  check_adjd <- any(unlist(lapply(object$lambda.object$coefficients_lambdaS, names))=="adjd")
  if (check_adjd) {
    if (object$spaMM.version < "2.2.116") {
      whichadj <- attr(attr(object$ZAlist,"ranefs"),"type") %in% c("adjacency")  
    } else whichadj <- attr(object$ZAlist,"exp_ranef_types") %in% c("adjacency")  
    rfl[whichadj] <- "inverse[ lambda_i =var(V'u) ]"
    rff <- rff[!whichadj] # and remove such terms from the family info
  }
  #
  if (type=="family") {
    # print family info except for adjd
    urff <- unique(rff)
    abyss <- lapply(urff, switch,
                    "Beta" = cat("lambda = 4 var(u)/(1 - 4 var(u)) for u ~ Beta[1/(2*lambda),1/(2*lambda)]; \n"),
                    "inverse.Gamma" = cat("lambda = var(u)/(1 + var(u)) for u ~ inverse-Gamma(sh=1+1/lambda, rate=1/lambda); \n"),
                    "Gamma" = cat("lambda = var(u) for u ~ Gamma(sh=1/lambda, sc=1/lambda); \n"),
                    "gaussian" = cat("lambda = var(u) for u ~ Gaussian; \n")
    )
  } else if (type=="link") {
    # print link info
    urfl <- unique(rfl[in_table])
    cat(paste0("             --- Coefficients for ",paste(urfl,collapse=" or "),":"))
  }
  invisible(NULL)
}


.MLmess <- function(object,ranef=FALSE) {
  if (object$models[["eta"]]=="etaGLM") {
    return("by ML.")
  } else if (object$family$family=="gaussian" && all(attr(object$rand.families,"lcrandfamfam")=="gaussian")) { 
    return("by ML.") 
  } else {
    if (object$HL[1]=='SEM')  {
      return("by stochastic EM.")
    } else if (object$HL[1]==1L)  {
      return("by Laplace ML approximation (p_v).")
    } else if (object$HL[1]==0L)  {
      if (ranef) {
        return("by Laplace ML approximation (p_v).")
      } else return("by h-likelihood approximation.")
    } 
  }
}
## FR->FR il faudrait distinguer EQL approx of REML ?
.REMLmess <- function(object,return_message=TRUE) {
  ## 'object' has no 'processed' element but its processed$REMLformula was copied to 'REMLformula' element.
  ## ./. It is by default NULL if REML was used, but may be an explicit non-default formula
  ## ./. It is an explicit formula if ML was used
  if (return_message) {
    if (is.null(object$REMLformula)) { ## default REML case
      if (object$HL[1]=='SEM')  {
        resu <- ("by stochastic EM.")
      } else if (object$family$family !="gaussian" 
                 || (object$models[["eta"]]=="etaHGLM" && any(attr(object$rand.families,"lcrandfamfam")!="gaussian"))) { 
        resu <- ("by Laplace REML approximation (p_bv).") 
      } else {
        resu <- ("by REML.")
      }  
    } else {
      if (identical(attr(object$REMLformula,"isML"),TRUE)) { ## FALSE also if object created by spaMM <1.9.15 
        resu <- (.MLmess(object, ranef=TRUE))
      } else { ## if nontrivial REML formula was used...
        resu <- ("by non-standard REML")
        attr(resu,"fixeformFromREMLform") <- .stripRanefs(object$REMLformula)
      }
    }    
  } else return( is.null(object$REMLformula) && object$HL[1]!='SEM')
  return(resu)
}


summary.HLfitlist <- function(object, ...) {
  if (object[[1]]$spaMM.version<"1.11.57") {
    warning("This fit object was created with spaMM version<1.11.57, and is no longer supported.\n Please recompute it.")
  } ## warning added in v2.2.43 2017/10/28
  sapply(object,summary.HLfit) ## because summary(list object) displays nothing (contrary to print(list.object)) => rather call summary(each HLfit object)
  cat(" ======== Global likelihood values  ========\n")    
  zut <- attr(object,"APHLs")
  cat(paste0(names(zut),": ",signif(unlist(zut),6), collapse="; "),"\n")
  invisible(object)
}

.prettify_family <- function(family,linkstring="") {
  famfam <- family$family
  if ( ! is.null(withArgs <- attr(famfam,"withArgs"))) {
    withArgs <- eval(withArgs,envir=environment(family$aic))
    legend <- paste(withArgs, "( link =", family$link,")")
  } else legend <- paste(famfam, "( link =", family$link,")")
  if (identical(family$zero_truncated, TRUE)) legend <- paste("0-truncated", legend)
  return(legend)
}

.get_compact_cov_mats <- function(strucList,later) {
  isRandomSlope <- attr(strucList,"isRandomSlope")
  if (any(isRandomSlope)) {
    cov.mats <- vector("list",length(strucList))
    for (it in seq_along(strucList)) {
      if ( isRandomSlope[it] ) {
        cov.mats[[it]] <- attr(strucList[[it]],"latentL_blob")$compactcovmat
      }
      ## keep NULL slots for other elements as expected by summary.HLfit
    }
    return(cov.mats)
  } else return(NULL)
}

# conversion of columns of data frames:
.aschar_adhoc <- function(data) {
  for (colit in seq_len(ncol(data))) if (is.numeric(colval <- data[,colit])) {
    charvec <- as.character(signif(colval,4))
    charvec[is.na(charvec)] <- "" ## the reason why we create an ad hoc fn.
    data[,colit] <- charvec
  }
  return(data)
}


.display_raw_lambdas <- function(in_pointLambda, row_map, lambda.object) {
  # if (details$ranCoefs) {
  #   displaypos <- innerlambda_pos
  # } else displaypos <- setdiff(innerlambda_pos, random_slope_pos)
  # displayrows <- unlist(row_map[displaypos])
  #if ( ! details$ranCoefs) displaypos <- setdiff(displaypos, random_slope_pos)
  displaypos <- which(in_pointLambda)
  displayrows <- unlist(row_map[displaypos])
  nicertypes <- lambda.object$type[displaypos]
  nicertypes[ ! nicertypes=="fixed"] <- ""
  nicertypes[ nicertypes=="fixed"] <- "[fixed]"
  nicertypes <- rep(nicertypes, unlist(lapply(row_map[displaypos], length)))
  if ( ! is.null(displayrows)) {
    print_lambda <- unlist(lambda.object$lambda)
    cat(paste("  ",
              names(displayrows)," : ", 
              signif(print_lambda[displayrows],4),
              nicertypes,
              collapse="\n"),"\n")
  }
}

.lambda_table_fn <- function(namesTerms, object, lambda.object, linklam_coeff_list=lambda.object$coefficients_lambdaS) {
  nrand <- length(namesTerms)
  nrows <- unlist(lapply(namesTerms,length))
  cum_nrows <- cumsum(c(0,nrows))
  names(cum_nrows) <- NULL
  row_map <- lapply(nrows,seq)
  for (it in seq_len(length(row_map))) row_map[[it]] <- row_map[[it]]+cum_nrows[it]
  # first construct a table including NA's for some coeeficients (not "inner" estimated), then remove these rows
  repGroupNames <- rep(names(namesTerms),sapply(namesTerms,length))
  ## i.e. for namesTerms = list(Subject=c("(Intercept)", "Days")), repGroupNames[[1]] is c("Subject", "Subject")
  lambda_table <- data.frame(Group=repGroupNames,Term=unlist(namesTerms))
  in_table <- rep(FALSE,nrand)
  in_pointLambda <- rep(TRUE,nrand)
  cov.mats <- .get_compact_cov_mats(object$strucList,later=TRUE)
  if ( length(cov.mats)) { ## fixme ? rename cov.mats to refer to ranCoefs ?
    maxnrow <- cum_nrows[nrand+1] ## maxnrow should = nrow(lambda_table)
    #.varcorr <- function(nrows, maxnrow, cov.mats, in_table, in_pointLambda, cum_nrows) {
    summ_corr_cols <- data.frame(matrix(NA,ncol=max(nrows-1L),nrow=maxnrow))
    summ_variances <- data.frame(matrix(NA,ncol=1L,nrow=maxnrow))
    for (mt in seq_len(length(cov.mats))) { 
      m <- cov.mats[[mt]]
      if ( ! is.null(m)) {
        in_table[mt] <- TRUE
        in_pointLambda[mt] <- FALSE
        inrows <-  cum_nrows[mt]+(1:nrow(m))
        summ_variances[inrows,1] <- diag(m)
        m <- stats::cov2cor(m)
        for (it in (2:nrow(m))) {
          for (jt in (1:(it-1))) {
            summ_corr_cols[cum_nrows[mt]+it,jt] <- m[it,jt]
          }
        }
      }
    }
    colnames(summ_corr_cols) <- rep("Corr.",ncol(summ_corr_cols))
    colnames(summ_variances) <- "Var."
    # }
    random_slope_pos <- which( ! unlist(lapply(cov.mats,is.null))) ## fixme ? equivalentto isRandomSlope that might be available
    random_slope_rows <- unlist(row_map[ random_slope_pos ])
  } else random_slope_rows <- random_slope_pos <- integer(0)
  if ( ! is.null(linklam_coeff_list)) {
    in_table[which( ! unlist(lapply(linklam_coeff_list,is.null)))] <- TRUE
    lambda_table <- cbind(lambda_table, Estimate=unlist(linklam_coeff_list), 
                          "Cond.SE"=lambda.object$lambda_se)
    
    info_rows <- which(! is.na(lambda_table$Estimate)) ## must be evaluated before the next line sets more NAs
  } else info_rows <- NULL
  for (it in seq_len(length(namesTerms))) {
    if ("adjd" %in% namesTerms[[it]]) {
      in_pointLambda[it] <- FALSE
      in_table[it] <- TRUE 
    } 
  }
  if ( length(cov.mats)) lambda_table <- cbind(lambda_table, summ_variances, summ_corr_cols)
  lambda_table <- structure(lambda_table, 
                            class=c("lambda_table",class(lambda_table)), info_rows=info_rows,
                            random_slope_rows=random_slope_rows,random_slope_pos=random_slope_pos,
                            in_pointLambda=in_pointLambda,row_map=row_map, in_table=in_table)
  return(lambda_table)
}

.print_adj_ests <- function(object, namesTerms, linklam) {
  ncoeffs <- attr(object$ZAlist,"Xi_cols") ## RHS = 2 for random slope, else 1
  any_adjd <- FALSE
  for (it in seq_len(length(namesTerms))) if ("adjd" %in% namesTerms[[it]]) {
    ncoeffs[it] <- 2L
    any_adjd <- TRUE
  }
  if (any_adjd) {
    cat("           --- Variance parameters ('lambda'):\n") ## NOT the (Gamma-GLM) coefficients and SEs
  }
  cum_ncoeffs <- c(0,cumsum(ncoeffs))
  for (it in seq_len(length(namesTerms))) {
    if ("adjd" %in% namesTerms[[it]]) {
      namenames <- names(namesTerms[it])
      pos <- cum_ncoeffs[it]+1L
      cat(paste("Estimate of rho (",namenames,"CAR): ",
                signif( - linklam[pos+1L]/linklam[pos],4),"\n"))
      cat(paste("Estimate of lambda factor (",namenames,"CAR): ",
                with(object$lambda.object,signif(linkinvS[[rand_to_glm_map[it]]](linklam[pos]),4)),"\n"))
    } 
  }
  return(any_adjd)
}

.print_lambda_table <- function(object,lambda_table, details, namesTerms, linklam_coeff_list) {
  attribs <- attributes(lambda_table) 
  # strings for screen output:
  lambda_table <- .aschar_adhoc(lambda_table)
  if (length(attribs$random_slope_rows)) {
    keep <- which( ! colnames(lambda_table) %in% c("Estimate","Cond.SE")) ## may be null if only outer... + fixed
    cov_table <- lambda_table[attribs$random_slope_rows,keep] ## EXCLUDES the "Estimate","Cond.SE" cols
    cat("         --- Random-coefficients Cov matrices:\n")
    print(cov_table, digits = 4, row.names = FALSE)
  }
  #
  keep <- which( colnames(lambda_table) %in% c("Group","Term","Estimate","Cond.SE")) ## may be null if only outer... + fixed
  lambda_table <- lambda_table[,keep] 
  #
  # point lambda's (+adjd rho) estimates, NOT the (Gamma-GLM) coefficients and SEs
  linklam <- unlist(linklam_coeff_list)
  any_adjd <- .print_adj_ests(object, namesTerms, linklam)
  if (any(attribs$in_pointLambda)) {
    if ( ! any_adjd) cat("           --- Variance parameters ('lambda'):\n") 
    .legend_lambda(object, type="family") ## print info about the meaning of lambda according to the rand family (not rand link)
    .display_raw_lambdas(in_pointLambda=attribs$in_pointLambda, 
                         row_map=attribs$row_map, 
                         object$lambda.object) ## not the table with SEs / covariances
  }
  #
  # NOW the (Gamma-GLM) coefficients and SEs
  info_rows <- attribs$info_rows
  if ( ! details$ranCoefs) info_rows <- setdiff(info_rows, attribs$random_slope_rows) ## removes random-coef info except if detaisl are requested
  lambda_table <- lambda_table[ info_rows ,]
  in_table <- attribs$in_table
  if (nrow(lambda_table)) { 
    if (any(in_table) && ! is.null(linklam_coeff_list[in_table])) { ## (mixed cov.mat and) coefficients_lambdaS output
      .legend_lambda(object, in_table, type="link")
      cat("\n")
    # } else { ## only ranCoefs
    #   cat("Variances and correlation for random-coefficient terms:\n")
    } ## else only other terms with outer lambda (eg Matern term) 
    print(lambda_table,digits=4,row.names=FALSE)
  }
  invisible(lambda_table) ## not currently correct as it is truncated
}




`summary.HLfit` <- function(object, details=FALSE, max.print=100L, verbose=TRUE, ...) { 
  if ( ! verbose) {
    mc <- match.call(expand.dots = TRUE)
    mc$verbose <- TRUE
    capture.output({silent <- eval(mc)})
    return(silent)
  }
  oldopt <- options(max.print=max.print)
  if (is.null(names(details))) details <- structure(rep(details,2),names=c("ranCoefs","p_value")) ## handle FALSE or TRUE input
  details <- as.list(details)
  for (st in c("ranCoefs")) if (is.null(details[[st]])) details[st] <- FALSE
  for (st in c("p_value")) if (is.null(details[[st]])) details[st] <- "" ## a string such as "Wald"
  models <- object$models
  phi.object <- object$phi.object
  famfam <- object$family$family ## response !
  lcrandfamfam <- attr(object$rand.families,"lcrandfamfam") 
  randfamfamlinks <- unlist(lapply(object$rand.families, .prettify_family))
  randfamlinks <- unlist(lapply(object$rand.families, getElement, name="link"))
  summ <- list()
  cat("formula: ")
  form <- formula.HLfit(object)
  print(form,showEnv=FALSE)
  #
  #  HLchar <- paste(as.character(object$HL),collapse="")
  #  cat(paste0("[code: ",HLchar,"]"," method: "))
  messlist <- list()
  if (length(object$fixef)) messlist[["fixed"]] <- .MLmess(object)
  messlist[["ranef"]] <- .REMLmess(object)
  summ$formula <- object$formula
  summ$REMLformula <- object$REMLformula
  ## Distinguishing iterative algo within HLfit and numerical maximization outside HLfit 
  locblob <- .get_methods_disp(object)
  iterativeEsts <- locblob$iterativeEsts
  optimEsts <- locblob$optimEsts
  RE_Ests <- unique(c(iterativeEsts,optimEsts))
  len <- length(RE_Ests)
  lenIt <- length(iterativeEsts) ## p_v maxim, or p_bv(X.Re)
  lenOpt <- length(optimEsts) ## p_v maxim, or p_bv(X.Re)
  if (lenIt > 1) iterativeEsts <- paste(c(paste(iterativeEsts[-lenIt],collapse=", "),iterativeEsts[lenIt]),collapse=" and ")
  if (len > 1) RE_Ests <- paste(c(paste(RE_Ests[-len],collapse=", "),RE_Ests[len]),collapse=" and ")
  if (lenIt) {
    if (messlist[["ranef"]]=="by REML.") {## REMLmess has checked that this is a LMM
      cat("REML: Estimation of ")
      cat(RE_Ests);
      tab <- "      "
    } else if (messlist[["ranef"]]=="by ML.") {
      cat("ML: Estimation of ")      
      cat(RE_Ests);
      tab <- "    "
    } else {
      cat("Estimation of ")
      cat(RE_Ests);
      tab <-""
    }
    ## the aim of the tab is to align "Estimation of..." vertically
    cat(" ")
    cat(messlist[["ranef"]])
    if ( messlist[["ranef"]]=="by non-standard REML") {
      cat("\n");cat(tab);cat(" based on fixed-effects model: ")
      print(attr(messlist[["ranef"]],"fixeformFromREMLform"),showEnv=FALSE) 
    } else cat("\n") ## normal output for standard REML formula
  } else tab <- ""
  if (length(object$fixef)) {
    cat(tab)
    cat("Estimation of fixed effects ")
    cat(messlist[["fixed"]]);
    cat("\n")
  }  
  if (lenOpt > 1) optimEsts <- paste(c(paste(optimEsts[-lenOpt],collapse=", "),optimEsts[lenOpt]),collapse=" and ")
  if (lenOpt > 0) { 
    objective <- attr(object,"optimInfo")$objective  
    if(is.null(objective)) {
      #stop("attr(object,'optimInfo')$objective is missing: malformed object.")
      ## may happen when one refits an HLCorcall =>
      ## optimEsts ultimately deduced by its $processed$ranPars as in the original fit,
      ## but no optimInfo in the refit.
    } else {
      objString <- switch(objective,
                          p_bv= "'outer' REML, maximizing p_bv",
                          p_v= "'outer' ML, maximizing p_v",
                          cAIC= "'outer' minimization of cAIC",
                          paste("'outer' maximization of",objective)
      )
      outst <- paste0("Estimation of ",optimEsts," by ",objString,".\n")
      cat(outst) 
    }
  }
  cat("Family:", .prettify_family(object$family, linkstring = "link = ") , "\n") 
  summ$family <- object$family
  if (length(object$fixef)==0L) {
    cat("No fixed effect\n")
  } else {
    cat(" ------------ Fixed effects (beta) ------------\n")
    beta_table <- .make_beta_table(object, p_value=details$p_value)
    print(beta_table,4)
    summ$beta_table <- beta_table
  }
  if (models[["eta"]]=="etaHGLM") {
    cat(" --------------- Random effects ---------------\n") 
    urff <- unique(lcrandfamfam)
    urffl <- unique(randfamfamlinks)
    if (length(urffl)==1L) { 
      cat("Family:", .prettify_family(object$rand.families[[1]], linkstring = "link = ") , "\n") 
    } else {
      cat("Families(links):", paste(randfamfamlinks,collapse=", "), "\n")
    }
    corrPars <- get_ranPars(object,which="corrPars")
    cP <- unlist(corrPars)
    if ( ! is.null(cP) ) {
      moreargs <- attr(object,"optimInfo")$LUarglist$moreargs
      control_dists <- lapply(moreargs,getElement,name="control.dist")
      dist_methods <- lapply(control_dists,getElement,name="dist.method")
      if (any((udm <- unlist(dist_methods))!="Euclidean")) {
        ocd <- rep("Euclidean",length(dist_methods))
        names(ocd) <- names(dist_methods)
        ocd[names(udm)] <- udm
        cat("Distance(s):",paste(ocd, collapse=","),"\n")
      }
      cat("                   --- Correlation parameters:")
      corrFixNames <- names(unlist(corrPars[which(attr(corrPars,"type")=="fix")]))
      if (length(corrFixNames)>1) {
        cat(" [",paste0(corrFixNames,collapse=",")," were fixed]")
      } else if (length(corrFixNames)==1L) cat(" [",corrFixNames," was fixed]")
      cat("\n")
      print(cP)
    }
    lambda.object <- object$lambda.object
    namesTerms <- lambda.object$print_namesTerms ## list of vectors of variable length
    if (any(object$models[["lambda"]] == "lamHGLM")) { 
      stop("voir ici dans summary.HLfit")
    } else {
      linklam_coeff_list <- lambda.object$coefficients_lambdaS ## used beyond the next line
      summ$lambda_table <- lambda_table <- .lambda_table_fn(namesTerms, object, lambda.object,linklam_coeff_list)
      .print_lambda_table(object,lambda_table, details=details, namesTerms, linklam_coeff_list) 
    } 
    cat(paste0("# of obs: ",nrow(object$data),"; # of groups: ",
              paste0(names(namesTerms),", ",unlist(lapply(object$ZAlist,ncol)), collapse="; ")
              ), "\n")
  }
  ##
  if (object$family$family %in% c("gaussian","Gamma")) {
    if (object$family$family=="Gamma") {
      cat(" -- Residual variation ( var = phi * mu^2 )  --\n")
    } else cat(" ------------- Residual variance  -------------\n")    
    pw <- object$prior.weights
    if ( ! (identical(attr(pw,"unique"),TRUE) && pw[1]==1L)) cat(paste("Prior weights:",
                                                          paste(signif(pw[1:min(5,length(pw))],6),collapse=" "),
                                                          "...\n"))
    if ( ! is.null(phi_outer <- phi.object$phi_outer)) {
      if ( identical(attr(phi_outer,"type"),"fix") ) {
        if (length(phi_outer)==1L) {
          cat(paste("phi was fixed to",signif(phi_outer,6),"\n"))
        } else  cat(paste("phi was fixed.\n"))
      } else {
        if (length(phi_outer)==1L) {
          cat(paste("phi estimate was",signif(phi_outer,6),"\n"))
        } else  cat(paste("phi was estimated.\n"))
      }
      summ$phi_outer <- phi_outer
    } else {
      if (models[["phi"]]=="phiHGLM") {
        cat("Residual dispersion model includes random effects:\n  use summary(<fit object>$resid_fit) to display results.\n")       
      } else if ((loc_p_phi <- length(phi.object$fixef))) {
        glm_phi <- phi.object[["glm_phi"]]
        if (is.null(glm_phi)) glm_phi <- .get_glm_phi(object)
        phi_se <- summary(glm_phi,dispersion=1)$coefficients[(loc_p_phi+1L):(2L*loc_p_phi)]
        ## note dispersion set to 1 to match SmythHV's 'V_1' method, which for a log link has steps:
        #SmythHVsigd <- as.vector(sqrt(2)*phi_est);SmythHVG <- as.vector(phi_est); tmp <- SmythHVG / SmythHVsigd 
        ## tmp is here sqrt(2) !
        #if (length(tmp)>1) {SmythHVZstar <- diag(tmp) %*% X_disp} else SmythHVZstar <- tmp * X_disp
        #SmythHVcovmat <- solve(.ZtWZ(SmythHVZstar,(1-lev_phi))); phi_se <- sqrt(diag(SmythHVcovmat)) print(phi_se)
        #
        phi_table <- cbind(phi.object$fixef,phi_se)
        colnames(phi_table) <- c("Estimate", "Cond. SE")
        rownames(phi_table) <- namesX_disp <- names(phi.object$fixef)
        summ$phi_table <- phi_table
        phiform <- object$resid.predictor
        resid.family <- eval(object$resid.family)
        phiinfo <- resid.family$link 
        if (phiinfo=="identity") {phiinfo="phi "} else {phiinfo <- paste0(phiinfo,"(phi) ")}
        phiinfo <- paste("Coefficients for",phiinfo,paste0(phiform,collapse=" ")," :\n")
        cat(phiinfo)
        print(phi_table,4)
        dispOffset <- attr(object$resid.predictor,"off")
        if (!is.null(dispOffset)) dispOffset <- unique(dispOffset)
        if (length(namesX_disp)==1 && namesX_disp[1]=="(Intercept)" && length(dispOffset)<2L) {
          # constant phi: we can display it
          phi_est <- (phi.object$fixef)
          if (length(dispOffset)==1L) phi_est <- phi_est+dispOffset
          phi_est <- resid.family$linkinv(phi_est)
          if (object$family$family=="Gamma") {
            cat(paste("Estimate of phi: ",signif(phi_est,4),"\n"))
            ## la var c'est phi mu^2...
          } else cat(paste("Estimate of phi=residual var: ",signif(phi_est,4),"\n"))
        } ## else phi not constant; We don't try to display it
        wa <- glm_phi$warnmess
        if ( ! is.null(wa)) {
          if (wa=="glm.fit: algorithm did not converge") {
            cat("glm.fit for estimation of phi SE did not converge; this suggests\n")
            cat(" non-identifiability of some phi (and possibly also lambda) coefficients.\n")
          } else {
            cat("warning in glm.fit for estimation of phi SE: \n")
            cat(wa,"\n")
          }
        }
      } else {
        phiform <- object$resid.predictor
        if (length(phiform)==2) phiform <- as.formula(paste('"phi"',paste(phiform,collapse=" "))) ##FR->FR how does _dglm_ deal with this
        cat(paste("phi was fixed by an offset term: ",deparse(phiform) ,"\n")) ## quick fix 06/2016 
      }                                                 
    }
  } ## else binomial or poisson, no dispersion param
  if (object$HL[1]==0L) { 
    validnames <- intersect(names(object$APHLs),c("hlik","p_v","p_bv"))
  } else validnames <- intersect(names(object$APHLs),c("p_v","p_bv"))
  if (length(validnames)) { ## may be 0 in SEM...
    likelihoods <- unlist(object$APHLs[validnames]) # NULL if no validnames
    if ( models[["eta"]]=="etaHGLM"){
      APHLlegend <- c(hlik="       h-likelihood:",
                      p_v="p_v(h) (marginal L):",
                      p_bv="  p_beta,v(h) (ReL):")
    } else APHLlegend <- c(p_v="p(h)   (Likelihood):",
                           p_bv="  p_beta(h)   (ReL):")
    names(likelihoods) <- APHLlegend[validnames]
    if ( is.null(object$distinctX.Re)) {
      ## standard REML 
    } else {
      whichp_bv <- which(validnames=="p_bv")
      if (ncol(object$distinctX.Re)==0L) {
        likelihoods <- likelihoods[-whichp_bv] ## ML 
      } else names(likelihoods)[whichp_bv] <- "(Non-standard?)  ReL:"
    }
  } else likelihoods <- numeric(0)
  # messlist[["ranef"]]
  logLapp <- object$APHLs$logLapp
  if (!is.null(logLapp)) {
    locli <- list(logLapp[1]) ## [1] removes attribute
    names(locli)[1] <- attr(logLapp,"method")
    likelihoods <- c(likelihoods,locli)
  }
  cat(" ------------- Likelihood values  -------------\n")    
  astable <- as.matrix(likelihoods);colnames(astable)<-"logLik";
  print(astable)
  summ$likelihoods <- likelihoods
  if (length(object$warnings) ) { 
    silent <- sapply(object$warnings, cat, "\n") 
  }
  if (object$family$family=="Gamma" && 
      object$family$link=="inverse" &&
      identical(attr(object$muetablob$mu,"any_neg_eta"),TRUE) ## attribute added in version 2.4.100
  ) {
    cat("Non-positive predictions implied by final fitted coefficients in Gamma(inverse)-response model.\n")
  }
  options(oldopt)
  invisible(summ)
}

print.HLfit <- function(x,...) {
  summary(x,...)
  invisible(x)
}

print.HLfitlist <- function(x,...) {
  summary(x,...)
  invisible(x)
}


