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
  optimEsts <- c(optimEsts,intersect(optimNames,c(names(object$corrPars),"COMP_nu","NB_shape")))
  iterativeEsts <- c(iterativeEsts,setdiff(optimNames,optimEsts))
  return(list(iterativeEsts=iterativeEsts,optimEsts=optimEsts))
}


legend_lambda<- function(object) {
  ## (1) analyse distributions
  randfams <- object$rand.families
  rff <- sapply(seq_len(length(randfams)),function(rfit){tolower(randfams[[rfit]]$family)})
  rfl <- sapply(rff,function(st) {
    switch(st,
           "beta" = "log[ lambda = 4 var(u)/(1 - 4 var(u)) ]", ## u ~ Beta(1/(2lambda),1/(2lambda))
           "inverse.gamma" = "log[ lambda = var(u)/(1 + var(u)) ]", ## u ~ I-Gamma(1+1/lambda,1/lambda)
           "log[ lambda = var(u) ]" ## gaussian, or gamma with u ~ Gamma(lambda,1/lambda)
    )
  })
  check_adjd <- any(unlist(lapply(object$lambda.object$coefficients_lambdaS,function(v) any(names(v)=="adjd"))))
  if (check_adjd) {
    whichadj <- attr(attr(object$ZAlist,"ranefs"),"type") %in% c("adjacency","ar1")  
    rfl[whichadj] <- "inverse[ lambda_i =var(V'u) ]"
    rff <- rff[!whichadj]
  }
  urff <- unique(rff)
  urfl <- unique(rfl)
  if (length(urff)==1L && length(urfl)==1L) {
    cat(paste("Coefficients for ",urfl,":\n",sep=""))
  } else {
    cat(paste("Coefficients for ",paste(urfl,collapse=" or "),", with:\n",sep=""))
    abyss <- lapply(urff, function(st) {
      switch(st,
             "beta" = cat("lambda = 4 var(u)/(1 - 4 var(u)) for Beta distribution; \n"),
             "inverse.gamma" = cat("lambda = var(u)/(1 + var(u)) for inverse gamma distribution; \n"),
             "gamma" = cat("lambda = var(u) for gamma distribution; \n"),
             "gaussian" = cat("lambda = var(u) for Gaussian distribution; \n")
      )
    })     
  }  
  invisible(NULL)
}


.MLmess <-function(object,ranef=FALSE) {
  if (object$models[["eta"]]=="etaGLM") {
    return("by ML.")
  } else if (object$family$family=="gaussian" && all(attr(object$rand.families,"lcrandfamfam")=="gaussian")) { 
    return("by ML.") 
  } else {
    if (object$HL[1]=='SEM')  {
      return("by stochastic EM.")
    } else if (object$HL[1]==1L)  {
      return("by ML approximation (p_v).")
    } else if (object$HL[1]==0L)  {
      if (ranef) {
        return("by ML approximation (p_v).")
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
        resu <- ("by REML approximation (p_bv).") 
      } else {
        resu <- ("by REML.")
      }  
    } else {
      if (identical(attr(object$REMLformula,"isML"),TRUE)) { ## FALSE also if object created by spaMM <1.9.15 
        resu <- (.MLmess(object, ranef=TRUE))
      } else { ## if nontrivial REML formula was used...
        resu <- ("by non-standard REML")
        attr(resu,"fixeformFromREMLform") <- nobarsMM(object$REMLformula)
      }
    }    
  } else return( is.null(object$REMLformula) && object$HL[1]!='SEM')
  return(resu)
}


summary.HLfitlist <- function(object, ...) {
  sapply(object,summary.HLfit) ## because summary(list object) displays nothing (contrary to print(list.object)) => rather call summary(each HLfit object)
  cat(" ======== Global likelihood values  ========\n")    
  zut <- attr(object,"APHLs")
  cat(paste(names(zut),": ",signif(unlist(zut),6), sep="",collapse="; "),"\n")
  invisible(object)
}


`summary.HLfit` <- function(object, details=FALSE,...) {
  models <- object$models
  phi.object <- object$phi.object
  famfam <- object$family$family ## response !
  lcrandfamfam <- attr(object$rand.families,"lcrandfamfam") ## unlist(lapply(object$rand.families,function(rf) {tolower(rf$family)}))
  randfamfamlinks <- unlist(lapply(object$rand.families,function(rf) {paste(rf$family,"(",rf$link,")",sep="")}))
  randfamlinks <- unlist(lapply(object$rand.families,function(rf) {rf$link}))
  summ <- list()
  cat("formula: ")
  form <- attr(object$predictor,"oriFormula")
  if (is.null(form)) {
    form <- object$predictor ## valid case ?
    print(form)
  } else print(form,showEnv=FALSE)
  #
  #  HLchar <- paste(as.character(object$HL),collapse="")
  #  cat(paste("[code: ",HLchar,"]"," method: ",sep=""))
  messlist <- list()
  if (length(object$fixef)>0) messlist[["fixed"]] <- .MLmess(object)
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
  if (lenIt>0) {
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
  if (length(object$fixef)>0) {
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
                          paste("'outer' maximization of",objString)
      )
      outst <- paste("Estimation of ",optimEsts," by ",objString,".\n",sep="")
      cat(outst) 
    }
  }
  if ( ! is.null(withArgs <- attr(famfam,"withArgs"))) {
    withArgs <- eval(withArgs,envir=environment(object$family$aic))
    cat("Family:", withArgs, "( link =", object$family$link,")\n")
  } else cat("Family:", famfam, "( link =", object$family$link,")\n")
  #if (famfam=="COMPoisson") cat("COMPoisson's nu=",signif(environment(object$family$aic)$nu,4),"\n")
  #if (famfam=="negbin") cat("neg.binomial's shape=",signif(environment(object$family$aic)$shape,4),"\n")
  summ$family <- object$family
  if (length(object$fixef)==0L) {
    cat("No fixed effect\n")
  } else {
    cat(" ------------ Fixed effects (beta) ------------\n")
    namesOri <- attr(object$X.pv,"namesOri")
    nc <- length(namesOri)
    betaOri_cov <- matrix(NA,ncol=nc,nrow=nc,dimnames=list(rownames=namesOri,colnames=namesOri))
    beta_cov <- get_beta_cov_any_version(object)
    betaOri_cov[colnames(beta_cov),colnames(beta_cov)] <- beta_cov
    beta_se <- sqrt(diag(betaOri_cov))
    fixef_z <- object$fixef/beta_se
    beta_table <- cbind(object$fixef,beta_se,fixef_z)
    colnames(beta_table) <- c("Estimate", "Cond. SE", "t-value")
    rownames(beta_table) <- names(object$fixef)
    print(beta_table,4)
    summ$beta_table <- beta_table
  }
  if (models[["eta"]]=="etaHGLM") {
    cat(" --------------- Random effects ---------------\n") 
    urff <- unique(lcrandfamfam)
    urffl <- unique(randfamfamlinks)
    if (length(urffl)==1L) { 
      cat("Family:", urff , "( link =", object$rand.families[[1]]$link,")\n") 
    } else {
      cat("Families(links):", paste(randfamfamlinks,collapse=", "), "\n")
    }
    corrPars <- object$corrPars
    cP <- unlist(corrPars)
    if ( ! is.null(cP) ) {
      if (! is.null(ocd <- object$control.dist$dist.method)) cat(paste("Distance:",ocd,"\n"))
      cat("Correlation parameters:")
      corrFixNames <- names(unlist(corrPars[which(attr(corrPars,"type")=="fix")]))
      if (length(corrFixNames)>1) {
        cat(" [",paste(corrFixNames,collapse=",")," were fixed]",sep="")
      } else if (length(corrFixNames)==1L) cat(" [",corrFixNames," was fixed]",sep="")
      cat("\n")
      print(cP)
    }
    lambda.object <- object$lambda.object
    namesTerms <- lambda.object$namesTerms ## list of vectors of variable length
    print_lambda <- unlist(lambda.object$lambda)
    type <- lambda.object$type
    if (any(object$models[["lambda"]] == "lamHGLM")) { 
      stop("voir ici dans summary.HLfit")
    } else if ( ! is.null(linklam_coeff_list <- lambda.object$coefficients_lambdaS)) {
      # first construct a table including NA's for some coeeficients (not "inner" estimated), then remove these rows
      repGroupNames <- unlist(lapply(seq_len(length(namesTerms)),function(it) {
        names(namesTerms[[it]]) <- rep(names(namesTerms)[it],length(namesTerms[[it]]))
      })) ## makes group identifiers unique (names of coeffs are unchanged)
      lambda_table <- data.frame(Group=repGroupNames,Term=unlist(namesTerms),
                                 Estimate=unlist(linklam_coeff_list),
                                 "Cond.SE"=lambda.object$lambda_se)
      is_info <- ! is.na(lambda_table$Estimate) ## must be evaluated before the cov.mats block sets more NAs
      nrand <- length(namesTerms)
      nrows <- unlist(lapply(namesTerms,length))
      cum_nrows <- cumsum(c(0,nrows))
      names(cum_nrows) <- NULL
      row_map <- lapply(nrows,seq)
      for (it in seq_len(length(row_map))) row_map[[it]] <- row_map[[it]]+cum_nrows[it]
      cov.mats <- object$cov.mats
      if ( ! is.null(cov.mats)) {
        maxnrow <- cum_nrows[nrand+1] ## maxnrow should = nrow(lambda_table)
        corr_cols <- data.frame(matrix("",ncol=max(nrows-1),nrow=maxnrow),stringsAsFactors=FALSE)
        variances <- data.frame(matrix("",ncol=1,nrow=maxnrow),stringsAsFactors=FALSE)
        for (mt in seq_len(length(cov.mats))) { ## assumes cov.mats for all effects
          m <- cov.mats[[mt]]
          if ( ! is.null(m)) {
            inrows <-  cum_nrows[mt]+(1:nrow(m))
            variances[inrows,1] <- paste(signif(lambdas <- diag(m),4))
            covtocorr <- diag(1/sqrt(lambdas))
            m <- covtocorr %*% m %*% covtocorr
            for (it in (2:nrow(m))) {
              for (jt in (1:(it-1))) {
                corr_cols[cum_nrows[mt]+it,jt] <- paste(signif(m[it,jt],4))
              }
            }
          }
        }
        colnames(corr_cols) <- rep("Corr.",ncol(corr_cols))
        colnames(variances) <- "Var."
        lambda_table <- cbind(lambda_table,variances,corr_cols)
        random_slope_pos <- which( ! unlist(lapply(cov.mats,is.null)))
        random_slope_rows <- unlist(row_map[ random_slope_pos ])
        if ( ! details) lambda_table[(random_slope_rows),3:4] <- NA 
      } else random_slope_pos <- integer(0)
      lambda_table <- lambda_table[ is_info,]
      summ$lambda_table <- lambda_table
      legend_lambda(object)
      print(lambda_table,digits=4,row.names=FALSE)
      # wa <-attr(lambda.object,"warning")
      # if ( ! is.null(wa)) {
      #   if (wa=="glm.fit: algorithm did not converge") {
      #     cat("glm.fit for estimation of lambda SE did not converge; this suggests\n")
      #     cat(" non-identifiability of some lambda (and possibly also phi) coefficients.\n")
      #   } else {
      #     cat("warning in glm.fit for estimation of lambda SE: \n")
      #     cat(wa,"\n")
      #   }
      # } 
      innerlambda_pos <- which(type=="inner")
      linklam <- unlist(linklam_coeff_list)
      ncoeffs <- attr(object$ZAlist,"Xi_cols") ## RHS = 2 for random slope, else 1
      for (it in seq_len(length(namesTerms))) if ("adjd" %in% namesTerms[[it]]) ncoeffs[it] <- 2
      cum_ncoeffs <- c(0,cumsum(ncoeffs))
      for (it in seq_len(length(namesTerms))) {
        if ("adjd" %in% namesTerms[[it]]) {
          namenames <- names(namesTerms[it])
          pos <- cum_ncoeffs[it]+1L
          cat(paste("Estimate of rho (",namenames,"CAR): ",
                    signif( - linklam[pos+1L]/linklam[pos],4),"\n"))
          cat(paste("Estimate of lambda factor (",namenames,"CAR): ",
                    with(lambda.object,signif(linkinvS[[rand_to_glm_map[it]]](linklam[pos]),4)),"\n"))
          innerlambda_pos <- setdiff(innerlambda_pos,it) ## remove from  generic innerlambda_pos printing below
        } 
      }
      if (length(innerlambda_pos)>0L) {
        if (details) {
          displaypos <- innerlambda_pos
        } else displaypos <- setdiff(innerlambda_pos, random_slope_pos)
        displayrows <- unlist(row_map[displaypos])
        if ( ! is.null(displayrows)) {
          cat(paste("Estimate of lambda (",
                    names(displayrows),"): ", 
                    signif(print_lambda[displayrows],4)
                    ,collapse="\n"),"\n")
        }
      }
    } 
    outerlambda_pos <- which(type=="outer")
    if (length(outerlambda_pos)>0L) cat(paste("Outer estimate of lambda (",
                                              names(namesTerms)[outerlambda_pos],"): ", 
                                              signif(print_lambda[outerlambda_pos],4)
                                              ,collapse="\n"),"\n")
    fixedlambda_pos <- which(type=="fix")
    if (length(fixedlambda_pos)>0L) cat(paste("Fixed lambda value (",
                                              names(namesTerms)[fixedlambda_pos],"): ", 
                                              signif(print_lambda[fixedlambda_pos],4),
                                              collapse="\n"),"\n")
    cat(paste("# of obs: ",nrow(object$data),"; # of groups: ",
              paste(names(namesTerms),", ",unlist(lapply(object$ZAlist,ncol)),
                    collapse="; ",sep=""),
              sep=""),"\n")
  }
  ##
  if (object$family$family %in% c("gaussian","Gamma")) {
    if (object$family$family=="Gamma") {
      cat(" -- Residual variation ( var = phi * mu^2 )  --\n")
    } else cat(" ------------- Residual variance  -------------\n")    
    pw <- object$prior.weights
    if ( ! identical(attr(pw,"unique"),TRUE) && pw[1]!=1L) cat(paste("Prior weights:",
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
        cat("Random effects in residual dispersion model:\n  use summary(<fit object>$resid_fit) to display results.\n")       
      } else if ((loc_p_phi <- length(phi.object$fixef))>0L) {
        glm_phi <- phi.object[["glm_phi"]]
        if (is.null(glm_phi)) glm_phi <- .get_glm_phi(object)
        phi_se <- summary(glm_phi,dispersion=1)$coefficients[(loc_p_phi+1L):(2L*loc_p_phi)]
        ## note dispersion set to 1 to match SmythHV's 'V_1' method, which for a log link has steps:
        #SmythHVsigd <- as.vector(sqrt(2)*phi_est);SmythHVG <- as.vector(phi_est); tmp <- SmythHVG / SmythHVsigd 
        ## tmp is here sqrt(2) !
        #if (length(tmp)>1) {SmythHVZstar <- diag(tmp) %*% X_disp} else SmythHVZstar <- tmp * X_disp
        #SmythHVcovmat <- solve(ZtWZ(SmythHVZstar,(1-lev_phi))); phi_se <- sqrt(diag(SmythHVcovmat)) print(phi_se)
        #
        phi_table <- cbind(phi.object$fixef,phi_se)
        colnames(phi_table) <- c("Estimate", "Cond. SE")
        rownames(phi_table) <- namesX_disp <- names(phi.object$fixef)
        summ$phi_table <- phi_table
        phiform <- attr(object$resid.predictor,"oriFormula")
        resid.family <- eval(object$resid.family)
        phiinfo <- resid.family$link 
        if (phiinfo=="identity") {phiinfo="phi "} else {phiinfo <- paste(phiinfo,"(phi) ",sep="")}
        phiinfo <- paste("Coefficients for ",phiinfo,paste(phiform,collapse=" ")," :\n",sep="")
        cat(phiinfo)
        print(phi_table,4)
        dispOffset <- attr(object$resid.predictor,"offsetObj")$total
        if (!is.null(dispOffset)) dispOffset <- unique(dispOffset)
        if (length(namesX_disp)==1 && namesX_disp[1]=="(Intercept)" && length(dispOffset)<2) {
          phi_est <- (phi.object$fixef)
          if (length(dispOffset)==1L) phi_est <- phi_est+dispOffset
          phi_est <- resid.family$linkinv(phi_est)
          if (object$family$family=="Gamma") {
            cat(paste("Estimate of phi: ",signif(phi_est,4),"\n"))
            ## la var c'est phi mu^2...
          } else cat(paste("Estimate of phi=residual var: ",signif(phi_est,4),"\n"))
        } 
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
        phiform <- attr(object$resid.predictor,"oriFormula")
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
      } else names(likelihoods)[whichp_bv] <- "Non-standard 'ReL':"
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
  if (length(object$warnings)>0 ) { 
    silent<-sapply(length(object$warnings),function(i) {cat(object$warnings[[i]]);cat("\n")}) 
  }
  invisible(summ)
}

print.HLfit <-function(x,...) {
  summary(x,...)
  invisible(x)
}

print.HLfitlist <-function(x,...) {
  summary(x,...)
  invisible(x)
}




