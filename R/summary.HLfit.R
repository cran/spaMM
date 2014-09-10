summary.HLfit <-
function(object, ...) {
  MLmess <-function() {
    if (object$models[["eta"]]=="etaGLM") {
      cat("by ML.\n")
    } else if (famfam=="gaussian" && lcrandfamfam=="gaussian") { 
      cat("by ML.\n") 
    } else {
      if (object$HL[1]=='SEM')  {
        cat("by stochastic EM.\n")
      } else if (object$HL[1]==1)  {
        cat("by ML approximation (p_v).\n")
      } else if (object$HL[1]==0)  cat("by h-likelihood approximation.\n") 
    }
  }
  ## FR->FR il faudrait distinguer EQL approx of REML ?
  REMLmess <-function() {
    if (object$HL[1]=='SEM')  {
      cat("by stochastic EM.\n")
    } else if (famfam !="gaussian" 
               || (object$models[["eta"]]=="etaHGLM" && lcrandfamfam!="gaussian")) { 
      cat("by REML approximation (p_bv).\n") 
    } else {
      cat("by REML.\n")
    }
  }
  models <- object$models
  lambda.object <- object$lambda.object
  phi.object <- object$phi.object
  famfam <- object$family$family ## response !
  lcrandfamfam <- unlist(lapply(object$rand.families,function(rf) {tolower(rf$family)}))
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
  if (length(object$fixef)>0) {
    cat("Estimation of fixed effects ") 
    MLmess()
  }
  summ$formula <- object$formula
  summ$REMLformula <- object$REMLformula
  ## Distinguishing iterative algo within HLfit and numerical maximization outside HLfit 
  iterativeEsts<-character(0)
  optimEsts<-character(0)
  ## What the HLfit call says:
  if (any(object$models[["lambda"]] != "")) {
    if (is.null(lambda.object$lambda.fix)) {iterativeEsts <- c(iterativeEsts,"lambda")} else {optimEsts <- c(optimEsts,"lambda")}
  }
  ## c("binomial","poisson"): phi is 1, not NULL
  if (is.null(phi.object$phi.Fix)) {iterativeEsts <- c(iterativeEsts,"phi")} else {
    if (! (famfam %in% c("binomial","poisson"))) optimEsts <- c(optimEsts,"phi")
  }
  ## better conceived code for corrPars:
  corrPars <- object$corrPars
  iterativeEsts <- c(iterativeEsts,names(which(attr(corrPars,"type")=="var")))
  optimEsts <- c(optimEsts,names(which(attr(corrPars,"type")=="fix")))
  ## 
  ranFixNames <- attr(object,"ranFixNames") 
  if ( ! is.null(ranFixNames) ) { ## we take info from corrHLfit to know which were really fixed in the corrHLfit analysis
    optimEsts <- optimEsts[!optimEsts %in% ranFixNames] ## ie those not fixed in the corrHLfit call
  }
  len <- length(iterativeEsts)
  if (len > 1) iterativeEsts <- paste(c(paste(iterativeEsts[-len],collapse=", "),iterativeEsts[len]),collapse=" and ")
  if (len > 0) { 
    cat("Estimation of ")
    cat(iterativeEsts)
    cat(" ")
    if (is.null(object$REMLformula)) {
      REMLmess()
    } else {
      fixeformFromREMLform <- nobarsMM(object$REMLformula) 
      if (length(fixeformFromREMLform)<2 ## ie 0 if original formula was  <~(.|.)> or 1 if ... <lhs~(.|.)>
          || object$REMLformula[[3]]=="0" ## this is the whole RHS; for fixed effect models
         ) {
        MLmess()
      } else { ## if nontrivial REML formula was used...
        cat("by non-standard REML \n based on fixed effects model: ")
        print(fixeformFromREMLform)
      }
    }
  }
  len <- length(optimEsts)
  if (len > 1) optimEsts <- paste(c(paste(optimEsts[-len],collapse=", "),optimEsts[len]),collapse=" and ")
  if (len > 0) { 
    objective <- object$objective  
    if ( ! is.null(objective) ) { 
      objString <- switch(objective,
                          p_bv= "'outer' REML, maximizing p_bv",
                          p_v= "'outer' ML, maximizing p_v",
                          paste("'outer' maximization of",objString)
      )
      outst <- paste("Estimation of ",optimEsts," by ",objString,".\n",sep="")
      cat(outst) 
    } ## else no outer optimization
  }
  cat("Family:", famfam, "( link =", object$family$link,")\n")
  summ$family <- object$family
  if (length(object$fixef)==0) {
    cat("No fixed effect\n")
  } else {
    cat(" ------- Fixed effects (beta) -------\n")
    beta_se <- sqrt(diag(object$beta_cov))
    fixef_z <- object$fixef/beta_se
    beta_table <- cbind(object$fixef,beta_se,fixef_z)
    colnames(beta_table) <- c("Estimate", "Cond. SE", "t-value")
    rownames(beta_table) <- names(object$fixef)
    print(beta_table,4)
    summ$beta_table <- beta_table
  }
  if (models[["eta"]]=="etaHGLM") {
    cat(" ---------- Random effects ----------\n") 
    urff <- unique(lcrandfamfam)
    urffl <- unique(randfamfamlinks)
    if (length(urffl)==1) { 
      cat("Family:", urff , "( link =", object$rand.families[[1]]$link,")\n") 
    } else {
      cat("Families(links):", paste(randfamfamlinks,collapse=", "), "\n")
    }
    cP <- unlist(corrPars)
    ranFixNames <- attr(object,"ranFixNames")
    if ( ! is.null(cP) ) {
      cat("Correlation parameters:")
      corrFixNames <- intersect(ranFixNames,names(cP))
      if (length(corrFixNames)>2) {
        cat(" [(",paste(corrFixNames,collapse=",")," were fixed]",sep="")
      } else if (length(corrFixNames)==1) cat(" [",corrFixNames," was fixed]",sep="")
      cat("\n")
      print(cP)
    }
    if (any(object$models[["lambda"]] == "lamHGLM")) { 
      stop("voir ici dans summary.HL")
    } else if (is.null(lambda.object$lambda.fix)) {
      namesTerms <- lambda.object$namesTerms
      repGroupNames <- unlist(lapply(seq_len(length(namesTerms)),function(it) {
        names(namesTerms[[it]]) <- rep(names(namesTerms)[it],length(namesTerms[[it]]))
      })) ## makes group identifiers unique (names of coeffs are unchanged)
      coefficients <- unlist(lambda.object$namesTerms)
      lambda_table <- data.frame(Group=repGroupNames,Term=coefficients,
                                 Estimate=lambda.object$linkscale.lambda,"Cond.SE"=lambda.object$lambda_se)
      cov.mats <- object$cov.mats
      if ( ! is.null(cov.mats)) {
        nrand <- length(namesTerms)
        nrows <- unlist(lapply(namesTerms,length))
        cum_nrows <- cumsum(c(0,nrows))
        maxnrow <- cum_nrows[nrand+1] ## should be nrow(lambda_table)
        blob <- data.frame(matrix("",ncol=max(nrows-1),nrow=maxnrow),stringsAsFactors=FALSE)
        variances <- data.frame(matrix("",ncol=1,nrow=maxnrow),stringsAsFactors=FALSE)
        for (mt in length(cov.mats)) { ## assumes cov.mats for all effects
          m <- cov.mats[[mt]]
          variances[(cum_nrows[mt]+1):cum_nrows[mt+1],1] <- paste(signif(lambdas <- diag(m),4))
          covtocorr <- diag(1/sqrt(lambdas))
          m <- covtocorr %*% m %*% covtocorr
          for (it in (2:nrow(m))) {
            for (jt in (1:(it-1))) {
              blob[cum_nrows[mt]+it,jt] <- paste(signif(m[it,jt],4))
            }
          }
        }
        colnames(blob) <- rep("Corr.",ncol(blob))
        colnames(variances) <- "Var."
        lambda_table <- cbind(lambda_table,variances,blob)
      }
      summ$lambda_table <- lambda_table
      legend_lambda(urff)
      print(lambda_table,digits=4,row.names=FALSE)
      if (length(lambda.object$namesX_lamres)==1 && lambda.object$namesX_lamres[1]=="(Intercept)") {
        cat(paste("Estimate of lambda: ",signif(exp(lambda.object$linkscale.lambda),4),"\n"))
      } 
      wa <-attr(lambda.object,"warning")
      if ( ! is.null(wa)) {
        if (wa=="glm.fit: algorithm did not converge") {
          cat("glm.fit for estimation of lambda SE did not converge; this suggests\n")
          cat(" non-identifiability of some lambda (and possibly also phi) coefficients.\n")
        } else {
          cat("warning in glm.fit for estimation of lambda SE: \n")
          cat(wa,"\n")
        }
      }
      cat(paste("# of obs: ",nrow(object$data),"; # of groups: ",paste(repGroupNames,", ",unlist(lapply(object$ZAlist,ncol)),collapse="; ",sep=""),sep=""),"\n")
    } else {
      cat(paste("lambda was fixed to",paste(signif(lambda.object$lambda.fix,6),collapse=","),"\n"))
      summ$lambda.fix <- lambda.object$lambda.fix
    }        
  }
  if (object$family$family %in% c("gaussian","Gamma")) {
    cat(" -------- Residual variance  --------\n")    
    if ( ! is.null(phi.object$phi.Fix)) {
      if (length(phi.object$phi.Fix)==1) {
        cat(paste("phi was fixed to",signif(phi.object$phi.Fix,6),"\n"))
      } else  cat(paste("phi was fixed.\n"))
      summ$phi.Fix <- phi.object$phi.Fix
    } else {
      if (models[["phi"]]=="phiHGLM") {
        stop("From summary.HL: phiHGLM code not ready")
      } else {
        phi_table<-cbind(phi.object$beta_phi,phi.object$phi_se)
        colnames(phi_table) <- c("Estimate", "Cond. SE")
        rownames(phi_table) <- phi.object$namesX_disp
        summ$phi_table <- phi_table
        cat("phi formula: ")
        phiform <- attr(object$resid.predictor,"oriFormula")
        if (length(phiform)==2) phiform <- as.formula(paste('"phi"',paste(phiform,collapse=" "))) ##FR->FR how does _dglm_ deal with this
        print(phiform,showEnv=FALSE)
        phiinfo <- object$resid.family$link; if (phiinfo=="identity") phiinfo=""
        phiinfo <- paste("Coefficients for ",phiinfo,"[ phi= ",sep="")
        if (object$family$family=="Gamma") { ## response family to know if its a scale param; not phifam, which is always Gamma
          phiinfo <- paste(phiinfo,"scale param. ]\n",sep="")
        } else phiinfo <- paste(phiinfo,"residual var ]\n",sep="")
        cat(phiinfo)
        print(phi_table,4)
        dispoff <- attr(object$resid.predictor,"offset")
        if (!is.null(dispoff)) dispoff <- unique(dispoff)
        if (length(phi.object$namesX_disp)==1 && phi.object$namesX_disp[1]=="(Intercept)" && length(dispoff)<2) {
          phi_est <- (phi.object$beta_phi)
          if (length(dispoff)==1) phi_est <- phi_est+dispoff
          phi_est <- object$resid.family$linkinv(phi_est)
          if (object$family$family=="Gamma") {
            cat(paste("Estimate of phi: ",signif(phi_est,4)," (residual var = phi * mu^2)\n"))
            ## la var c'est phi mu^2...
          } else cat(paste("Estimate of phi=residual var: ",signif(phi_est,4),"\n"))
        } 
      }                                                 
    }
  } ## else binomial or poisson, no dispersion param
  if ( models[["eta"]]=="etaHGLM") {
    likelihoods <- c("p_v(h) (marginal L):"=object$APHLs$p_v,"  p_beta,v(h) (ReL):"=object$APHLs$p_bv)
  } else {
    likelihoods <- c("p(h)   (Likelihood):"=object$APHLs$p_v,"  p_beta(h)   (ReL):"=object$APHLs$p_bv)
  }
  logLapp <- object$APHLs$logLapp
  if (!is.null(logLapp)) {
    locli <- list(logLapp[1]) ## [1] removes attribute
    names(locli)[1] <- attr(logLapp,"method")
    likelihoods <- c(likelihoods,locli)
  }
  if (!is.null(object$APHLs$logLsmooth)) likelihoods <- c(likelihoods," Smoothed marginal L: "=object$APHLs$logLsmooth)
  if (!is.null(object$APHLs$cAIC)) likelihoods <- c(likelihoods,"               cAIC:"=object$APHLs$cAIC)
  cat(" -------- Likelihood values  --------\n")    
  astable <- as.matrix(likelihoods);colnames(astable)<-"logLik";
  print(astable)
  summ$likelihoods <- likelihoods
  if (length(object$warnings)>0 ) { 
    silent<-sapply(length(object$warnings),function(i) {cat(object$warnings[[i]]);cat("\n")}) 
  }

  invisible(summ)
}
