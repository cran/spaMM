summary.HLfit <-
function(object, ...) {
  models<-object$models
  lambda.object<-object$lambda.object
  phi.object<-object$phi.object
  famfam <- object$family$family
  summ <- list()
  cat("formula: ")
  form <- object$predictor$oriFormula
  if (is.null(form)) form <- object$predictor$formula
  if (is.null(form)) form <- object$predictor ## still normal for some HLfit calls
  print(form,showEnv=F)
  HLchar <- paste(as.character(object$HL),collapse="")
  HLmethod <- switch(HLchar,
      "001" = "PQL",
      "011" = "HL(0,1)",
      "111" = "HL(1,1)",
      paste("[code: ",HLchar,"]",sep="") ## ugly info
    )
  cat(HLmethod)
  cat(" approximation to the likelihood.\n") 
  summ$formula <- object$formula
  summ$REMLformula <- object$REMLformula
  iterativeEsts<-character(0)
  optimEsts<-character(0)
  if (is.null(lambda.object$lambda.fix)) {iterativeEsts <- c(iterativeEsts,"lambda")} else {optimEsts <- c(optimEsts,"lambda")}
  if (is.null(phi.object$phi.Fix)) {iterativeEsts <- c(iterativeEsts,"phi")} else {if (! (famfam %in% c("binomial","poisson"))) optimEsts <- c(optimEsts,"phi")}
  len <- length(iterativeEsts)
  if (len > 1) iterativeEsts <- paste(c(paste(iterativeEsts[-len],collapse=", "),iterativeEsts[len]),collapse=" and ")
  if (len > 0) { 
    cat("Iterative estimation of ")
    cat(iterativeEsts)
    if (is.null(object$REMLformula)) {
      cat(" by REML.\n")
    } else {
      fixeformFromREMLform <- nobarsMM(object$REMLformula) 
      if (length(fixeformFromREMLform)<2) { ## ie 0 if orginal formula was  <~(.|.)> or 1 if ... <lhs~(.|.)>
        cat(" by ML.\n")
      } else { ## if nontrivial REML formula was used...
        cat(" by non-standard REML \n based on fixed effects model: ")
        print(fixeformFromREMLform,showEnv=F)
      }
    }
  }
  len <- length(optimEsts)
  if (len > 1) optimEsts <- paste(c(paste(optimEsts[-len],collapse=", "),optimEsts[len]),collapse=" and ")
  if (len > 0) { 
      cat("Estimation of ")
      cat(optimEsts)
      cat(" by numerical maximization of marginal or restricted likelihood.\n") ## no info in an HLfit output... ## outer objective p_v or p_bv
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
  cat(" ---------- Random effects ----------\n")    
  cat("Family:", object$rand.family$family, "\n")
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
      cat("Coefficients for log[ lambda=var(ranef 'u') ]\n")
      print(lambda_table,4)
      if (length(lambda.object$namesX_lambda)==1 && lambda.object$namesX_lambda[1]=="(Intercept)") {
        cat(paste("Estimate of lambda=var(ranef 'u'): ",signif(exp(lambda.object$linkscale.lambda),4),"\n"))
      } 
  } else {
      cat(paste("lambda was fixed to",signif(lambda.object$lambda.fix,6),"\n"))
      summ$lambda.fix <- lambda.object$lambda.fix
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
        stop("voir LA dans summary.HL")
      } else {
         phi_table<-cbind(phi.object$linkscale.phi,phi.object$phi_se)
         colnames(phi_table) <- c("Estimate", "Std. Error")
         rownames(phi_table) <- phi.object$namesX_disp
         summ$phi_table <- phi_table
         cat("phi formula: ")
         print(object$resid.predictor$formula,showEnv=F)
         # cat("Link: "); cat(object$RespLink_disp) # not useful, the cat already says it is log
         cat("Coefficients for log[ phi= residual var ]\n")
         print(phi_table,4)
         if (length(phi.object$namesX_disp)==1 && phi.object$namesX_disp[1]=="(Intercept)") {
           print(paste("Estimate of phi=residual var: ",signif(exp(phi.object$linkscale.phi),4)),quote=F)
         } 
      }                                                 
    }
  } ## else binomial or poisson, no dispersion param
  if ( models[1]=="etaHGLM") {
    likelihoods <- c("p_v(h) (marginal L):"=object$APHLs$p_v,"  p_beta,v(h) (ReL):"=object$APHLs$p_bv)
  } else {
    likelihoods <- c("         Likelihood:"=object$APHLs$c.lik)
  }
  if (!is.null(object$APHLs$caic)) likelihoods <- c(likelihoods,"               cAIC:"=object$APHLs$caic)
  cat(" -------- Likelihood values  --------\n")    
  astable <- as.matrix(likelihoods);colnames(astable)<-"logLik";
  print(astable)
  summ$likelihoods <- likelihoods
  if (length(object$warnings)>0 ) { 
    silent<-sapply(length(object$warnings),function(i) {cat(object$warnings[[i]]);cat("\n")}) 
  }
  invisible(summ)
}
