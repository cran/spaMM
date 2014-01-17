summary.HLfit <-
function(object, ...) {
  MLmess <-function() {
    if (famfam=="gaussian" && randfamfam=="gaussian") { 
      cat("by ML.\n") 
    } else {
      if (object$HL[1]==1)  {
        cat("by ML approximation (p_v).\n")
      } else if (object$HL[1]==0)  cat("by h-likelihood approximation.\n") 
    }
  }
  REMLmess <-function() {
    if (famfam=="gaussian" && randfamfam=="gaussian") { 
      cat("by REML.\n") 
    } else {
      cat("by REML approximation (p_bv).\n")
    }
  }
  models<-object$models
  lambda.object<-object$lambda.object
  phi.object<-object$phi.object
  famfam <- object$family$family
  randfamfam <- object$rand.family$family
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
  if (length(object$fixef)>=0) {
    cat("Estimation of fixed effects ") 
    MLmess()
  }
  summ$formula <- object$formula
  summ$REMLformula <- object$REMLformula
  ## Distinguishing iterative algo within HLfit and numerical maximization outside HLfit 
  iterativeEsts<-character(0)
  optimEsts<-character(0)
  ## What the HLfit call says:
  if (models[2] !="") {
    if (is.null(lambda.object$lambda.fix)) {iterativeEsts <- c(iterativeEsts,"lambda")} else {optimEsts <- c(optimEsts,"lambda")}
  }
  ## c("binomial","poisson"): phi is 1, not NULL
  if (is.null(phi.object$phi.Fix)) {iterativeEsts <- c(iterativeEsts,"phi")} else {if (! (famfam %in% c("binomial","poisson"))) optimEsts <- c(optimEsts,"phi")}
  ranFixNames <- object$ranFixNames 
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
                          p_bv= "'outer' REML",
                          p_v= "'outer' ML",
                          paste("'outer' maximization of",objString)
      )
      outst <- paste("Estimation of ",optimEsts," by ",objString,".\n",sep="")
      cat(outst) 
    } ## else no outer optimization
  }
  cat("Family:", famfam, "\n")
  cat("Link function:", object$family$link, "\n") ## this is print.family without some of the spaces...
  summ$family <- object$family
  if (length(object$fixef)==0) {
    cat("No fixed effect\n")
  } else {
    cat(" ------- Fixed effects (beta) -------\n")    
    fixef_z<-object$fixef/object$fixef_se
    beta_table<-cbind(object$fixef,object$fixef_se,fixef_z)
    colnames(beta_table) <- c("Estimate", "Std. Error", "t-value")
    rownames(beta_table) <- colnames(object$X)
    print(beta_table,4)
    summ$beta_table <- beta_table
  }
  if (models[1]=="etaHGLM") {
    cat(" ---------- Random effects ----------\n")    
    cat("Family:", randfamfam, "\n")
    cat("Link function:", object$rand.family$link, "\n") ## this is print.family wihtout some of the spaces...
    cP <- unlist(object$corrPars)
    if ( ! is.null(cP) ) {
      cat("Correlation parameters:\n")
      print(cP)
    }
    if (models[2]=="lamHGLM") { 
      stop("voir ici dans summary.HL")
    } else if (is.null(lambda.object$lambda.fix)) {
      lambda_table<-cbind(lambda.object$linkscale.lambda,lambda.object$lambda_se)
      colnames(lambda_table) <- c("Estimate", "Std. Error")
      rownames(lambda_table) <- lambda.object$namesRE ## that is for different '(1|namesRE)' terms
      summ$lambda_table <- lambda_table
      if (object$rand.family$family=="Beta") {
#        cat("Coefficients for log[ lambda ] for u ~ Beta(1/(2lambda),1/(2lambda))\n")
        cat("Coefficients for log[ lambda = 4 var(u)/(1 - 4 var(u)) ]:\n")
      } else if (object$rand.family$family=="inverse.Gamma") {
#        cat("Coefficients for log[ lambda ] for  u ~ I-Gamma(1+1/lambda,1/lambda)\n")
        cat("Coefficients for log[ lambda = var(u)/(1 + var(u)) ]:\n")
      } else if (object$rand.family$family=="Gamma") {
#        cat("Coefficients for log[ lambda = var(u) ] for  u ~ Gamma(lambda,1/lambda)\n")
        cat("Coefficients for log[ lambda = var(u) ]:\n")
      } else cat("Coefficients for log[ lambda = var(u) ]: \n") 
      print(lambda_table,4)
      if (length(lambda.object$namesX_lambda)==1 && lambda.object$namesX_lambda[1]=="(Intercept)") {
        cat(paste("Estimate of lambda: ",signif(exp(lambda.object$linkscale.lambda),4),"\n"))
      } 
    } else {
      cat(paste("lambda was fixed to",signif(lambda.object$lambda.fix,6),"\n"))
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
      if (models[3]=="phiHGLM") {
        stop("From summary.HL: phiHGLM code not ready")
      } else {
        phi_table<-cbind(phi.object$linkscale.phi,phi.object$phi_se)
        colnames(phi_table) <- c("Estimate", "Std. Error")
        rownames(phi_table) <- phi.object$namesX_disp
        summ$phi_table <- phi_table
        cat("phi formula: ")
        phiform <- attr(object$resid.predictor,"oriFormula")
        if (length(phiform)==2) phiform <- as.formula(paste('"phi"',paste(phiform,collapse=" ")))
        print(phiform,showEnv=FALSE)
        # cat("Link: "); cat(object$RespLink_disp) # not useful, the cat already says it is log
        if (object$family$family=="Gamma") {
          cat("Coefficients for log[ phi= scale param. ]\n")
        } else cat("Coefficients for log[ phi= residual var ]\n")
        print(phi_table,4)
        if (length(phi.object$namesX_disp)==1 && phi.object$namesX_disp[1]=="(Intercept)") {
          cat(paste("Estimate of phi=residual var: ",signif(exp(phi.object$linkscale.phi),4),"\n"))
          
        } 
      }                                                 
    }
  } ## else binomial or poisson, no dispersion param
  if ( models[1]=="etaHGLM") {
    likelihoods <- c("p_v(h) (marginal L):"=object$APHLs$p_v,"  p_beta,v(h) (ReL):"=object$APHLs$p_bv)
  } else {
    likelihoods <- c("p(h)   (Likelihood):"=object$APHLs$p_v,"  p_beta(h)   (ReL):"=object$APHLs$p_bv)
  }
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
