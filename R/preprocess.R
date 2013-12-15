preprocess <-
function(control.HLfit,phi.Fix=NULL,HLmethod,predictor,resid.predictor,REMLformula,data,family,BinomialDen=BinomialDen,rand.family) {
    processed <- list()
    ## algorithms (control of defaults remains in the HLfit code)
    processed$vUpdating <- control.HLfit$vUpdating
    processed$LevenbergM <-control.HLfit$LevenbergM ## 
    processed$spam <-control.HLfit$spam ##  
    stop.on.error <-control.HLfit$stop.on.error ##  
    if (is.null(stop.on.error)) stop.on.error <- F
    processed$stop.on.error <- stop.on.error ##  
    AIC <-control.HLfit$AIC ##  
    if (is.null(AIC)) AIC <- F
    processed$AIC <- AIC ##  
    essai <-control.HLfit$essai ##  
    if (is.null(essai)) essai <- F
    processed$essai <- essai ##  
    ## numerical control parameters 
    conv.threshold <-control.HLfit$conv.threshold ## 
    if (is.null(conv.threshold)) conv.threshold <- 1e-05
    processed$conv.threshold <- conv.threshold  
    #
    HL2.threshold <- control.HLfit$HL2.threshold 
    if (is.null(HL2.threshold)) HL2.threshold <- 3
#    processed$HL2.threshold <- HL2.threshold  
    #
    iter.mean.dispFix <-control.HLfit$iter.mean.dispFix ## private control
    if (is.null(iter.mean.dispFix)) iter.mean.dispFix <-control.HLfit$max.iter.mean ## public control
    if (is.null(iter.mean.dispFix)) iter.mean.dispFix <- 40 ## control of inner loop when no disp param is estimated
    processed$iter.mean.dispFix <- iter.mean.dispFix  
    #
    iter.mean.dispVar <-control.HLfit$iter.mean.dispVar ## private control
    if (is.null(iter.mean.dispVar)) iter.mean.dispVar <-control.HLfit$max.iter.mean ## public control ## control of inner loop when some disp param is estimated
    if (is.null(iter.mean.dispVar)) iter.mean.dispVar <- 20 ## control of inner loop when some disp param is estimated
    processed$iter.mean.dispVar <- iter.mean.dispVar  
    #
    max.iter <-control.HLfit$max.iter  ## control of outer loop 
    if (is.null(max.iter)) max.iter <-200
    processed$max.iter <- max.iter  
    #
    ## only check, no modif of data -> not returned
    if ( ! ("data.frame" %in% class(data))) {
      mess <- pastefrom("'data' is not a data.frame.",prefix="(!) From ")
      stop(mess)
    }
    if (class(predictor)=="formula") predictor <- Predictor(predictor) 
    #
    if (class(resid.predictor)=="formula") resid.predictor <- Predictor(resid.predictor)
    if (length(resid.predictor$formula)==2) resid.predictor$formula <- as.formula(paste('"phi"',paste(resid.predictor$formula,collapse=" ")))
    processed$resid.predictor <- resid.predictor  
    #
   # comme l'EQL est dans un monde quasi gaussien, il se ramene aux LMM et utilise les leverage standard,
   # Pour les GLMM non LMM, Dans la mesure ou HL(0,.) utilise h et non q+, ce n'est pas l'EQL pour les fixed params
   # Dans la mesure ou on utilise les leverages standard pour les dispersion param, c'est de l'EQL
    ## glm() syntax:  
    ## With binomial data the response can be either a vector or a matrix with two columns.
    ## If the response is a vector, it is treated as a binary factor with the first level representing "success" and all others representing "failure". 
    ## In this case R generates a vector of ones to represent the binomial denominators.
    ## Alternatively, the response can be a matrix where the first column shows the number of "successes" and the second column shows the number of "failures". 
    ## In this case R adds the two columns together to produce the correct binomial denominator. 
    if (family$family=="binomial") {
      if ( ! is.null(predictor$BinDenForm)) { ## the cbind syntax was used in the formula 
        BinomialDen <- eval(parse(text=predictor$BinDenForm),envir=data)
        ## la suite ducode suppose que pas cbind => ie essentially obsolete syntax 
      } else if (missing(BinomialDen) || is.null(BinomialDen)) { ## then this should be a binary response
        checkResp <- eval(parse(text=as.character(predictor$formula[2])),envir=data) ## 'positives'
        if (any(checkResp>1)) {
          mess <- pastefrom("binomial, non-binary response. Please use the",prefix="(!) From ")
          message(mess)
          message("    standard glm() syntax with _cbind_: 'cbind(<successes>, <failures>) ~ <etc>'")
          stop()
        } else BinomialDen <- rep(1,nrow(data)) ## response appears to be binary...
      }
      no.info <- (BinomialDen == 0)
      if (any(no.info)) {
        mess <- pastefrom("please remove missing data (i.e. for which binomial sample size is 0).",prefix="(!) From ")
        stop(mess)
      }
      ## It's not really possible to remove data at this stage as this may not match the dimension of the distance matrices
      ## moreover one cannot simply remove ros of a matrix "root"...
      ## it _could_ be useful to be able to hand BinomilaDen=0 by the general code but...
    } else {
      BinomialDen <- rep(1,length(nrow(data)))
    }
    processed$BinomialDen <- BinomialDen  
    ## conversion from user-friendly format to standard 'XX(...)' format
    ## first index is for (0) h / (1) p_v(h) / (2) p^s_v(h)
    ## second index is for no (0) or yes (1) D log Det / d log lambda correction (2) further p^s_bv correction
    ## third index is for use of (0) QL deviance residuals (this affects the leverages) or not (1)
    ## NohL07 table 1 has interesting terminology and further tables show even more cases
    if (HLmethod=="ML") {
       if (family$family=="binomial" && mean(BinomialDen)<HL2.threshold) {
         HLmethod <- "ML(1,1,1)" ## back to routine (1,1)
       } else {HLmethod <- "ML(1,1,1)"} ## that makes a log det correction on the log det of p_v, not p_bv...
       ## the 'ML' part is used only locally to determine the REMLformula ## this could be done !here! in the code
    } else if (HLmethod=="REML") {
      if (family$family=="binomial" && mean(BinomialDen)<HL2.threshold) {
         HLmethod <- "HL(1,1,1)" ## back to routine (1,1)
      } else {HLmethod <- "HL(1,1,1)"} 
    } else if (HLmethod=="REPQL") { ## (0,0 : no D log Det / d log lambda correction NohL07 table 7; ,1): does not use quasi deviance residuals
      HLmethod <- "HL(0,0,1)" ## currently no distinction between PQL and REPQL here
    } else if (HLmethod=="PQL/L") { ## was HL(0,0,1); again no D log Det / d log lambda correction
      HLmethod <- "ML(0,0,1)" 
    } else if (HLmethod=="PQL") { 
      message("Old-style 'PQL' is ambiguous, and is interpreted as 'PQL/L'; better use 'REPQL' or 'PQL/L' ")
      HLmethod <- "ML(0,0,1)" 
    } else if (HLmethod=="EQL") { ## tentative interpretation EQ has leverage corr hence is an 'HL' method
      HLmethod <- "HL(0,0,0)" ## (0,...): gaussianise everything, hence no a(1) correction ## there is no HL(1) a(1) correction in GLM.MME
    }
    HL <- eval(parse(text=paste("c",substr(HLmethod,3,100),sep=""))) ## extracts the (...) part into a vector
    if (length(HL)==2) HL <- c(HL,1)
    processed$HL <- HL  
    if (substr(HLmethod,0,2)=="ML") {
      if ( ! is.null(REMLformula)) {
        message("Confusing combination of arguments: 'HLmethod=ML(...)' with non-null 'REMLformula'.")
        stop("  Make sure what you mean and simplify the arguments.")
      }
      bars <- findbarsMM(predictor$formula)  ## extract random effects
      lhs <- paste(predictor$formula)[[2]] ## extract response, either cbind or not
      ## build formula with only random effects
      REMLformula <- as.formula(paste(lhs,"~",paste(paste("(",as.character(bars),")"),collapse="+")))
    }
    processed$REMLformula <- REMLformula  
    #
    processed$loglfn.fix <- selectLoglfn(family$family)
    MeanFrames <- HLframes(formula=predictor$formula,data=data) ## design matrix X, Y, fixef names
    ## in particular, $mf contains values of expl.vars and levels (indices) of ranefs for all observations. Indices are/should be ordered as uniqueGeo in ULI()
    X.pv <- MeanFrames$X  
    colnames(X.pv) <- names(MeanFrames$fixef)
#   namesX <- names(MeanFrames$fixef)
    processed$X.pv <- X.pv
    y <- MeanFrames$Y
    if ( family$family %in% c("binomial","poisson")) {
      ## the response variable should always be Counts
      if (max(abs(y-as.integer(y)))>1e-05) {
        mess <- pastefrom("response variable should be integral values.",prefix="(!) From ")
        stop(mess)
      }
    }
    if ( family$family == "binomial" && ncol(y)==2) y <- y[,1,drop=F] ## that is, we have the cbind syntax up to this fn
    ## code derived from the glm() function in the safeBinaryRegression package
#    if(family$family == "binomial" && length(unique(y)) == 2 && require(lpSolveAPI)) {
    if(family$family == "binomial" && length(unique(y)) == 2) {
      separation <- separator(X.pv, as.numeric(y), purpose = "test")$separation
      if(separation) {
        message("Separation exists among the sample points.\n\tThis model cannot be fit by maximum likelihood.")
        message("The following terms are causing separation among the sample points:")
        separation <- separator(X.pv, as.numeric(y), purpose = "find")$beta
        separating.terms <- dimnames(X.pv)[[2]][abs(separation) > 1e-09]
        if(length(separating.terms)) message(paste(separating.terms, collapse = ", "))
        stop()
      }
    }
    processed$y <- y
    #
    if ( ! is.null(REMLformula) ) { ## differences affects only REML estimation of dispersion params, ie which p_bv is computed
       REMLFrames <- HLframes(formula=REMLformula,data=data) ## design matrix X, Y, fixef names
       X.Re <- REMLFrames$X  
    } else {
       X.Re <- X.pv
    }
    processed$X.Re <- X.Re
    #
    canonicalLink <- FALSE
    if (family$family=="gaussian" && family$link=="identity") {
      canonicalLink <- TRUE
    } else if (family$family=="poisson" && family$link=="log") {
      canonicalLink <- TRUE
    } else if (family$family=="binomial" && family$link=="logit") {
      canonicalLink <- TRUE
    } else if (family$family=="Gamma" && family$link=="inverse") {
      canonicalLink <- TRUE
    }
    processed$canonicalLink <- canonicalLink
    #
    models <- c(eta="",lambda="",phi="")
    #
    list_ranefs<-findbarsMM(predictor$formula) ## a list with the random effect terms in the mean formula
    if (!is.null(list_ranefs)) { ## Design matriCES for v, typically a match between the levels of v and the observ. Ie Z, not ZL 
      models[["eta"]] <- "etaHGLM" 
      FL <- spMMFactorList(predictor$formula, MeanFrames$mf, 0L, 0L) ## this uses the spatial information in the formula, ieven if an explicit distm was used elsewhere
      namesRE <- FL$namesRE
      Zlist <- FL$Design ## : Zlist is a list of design matrices 
    } else { ## if is.null(random.mean): no random effect in linear predictor
      models[["eta"]] <- "etaGLM"
      namesRE <- NULL
      Zlist <- NULL
    } ## end if (!is.null(random.mean)) ... else ... endif
    processed$namesRE <- namesRE ## not sure what is the best way to pass this info
    processed$Zlist <- Zlist
    #
    nrand <- length(list_ranefs)
    if (nrand>0) {   
      nrand_lambda <- 0
      if (!is.null(predictor[["LinPredRandVariance"]])) { ## linear predictor for variance of random effects (lambda) (lamGLM or lamHGLM) 
        formulaLambda<-predictor[["LinPredRandVariance"]]
        if (length(formulaLambda)==2) formulaLambda <- as.formula(paste('"lambda"',paste(formulaLambda,collapse=" ")))
        lambda_ranefs <- findbarsMM(formulaLambda)
        if (!is.null(lambda_ranefs)) {  ## lamHGLM
          nrand_lambda <- length(lambda_ranefs)
          models[["lambda"]] <- "lamHGLM"
        } else models[["lambda"]] <- "lamGLM"  
        fr_lambda <- HLframes(formula=formulaLambda,data=data) ## but the "data" should probably be distinct data here, with nrow=number of reals of ranefs 
        X_lambda <- fr_lambda$X
        colnames(X_lambda) <- names(fr_lambda$fixef)
        mess <- pastefrom("LIKELY missing code to handle linear predictor for lambda.")
        stop(mess)
      } else { ## the day we have true 'data' for lambda we could avoid the following ad hoc code
        formulaLambda <- "lambda" ~ 1 
        models[["lambda"]] <- "lamScal"
        q <- rep(0, nrand)
        for (i in 1:nrand) q[i] <- ncol(Zlist[[i]]) ## (=nrow) of each design matrix = nb realizations each ranef
        qcum <- cumsum(c(0, q)) ## if two ranef,  with q=(3,3), this is 0,3,6. qcum[nrand+1] is then 6, the total # of realizations
        X_lambda <- matrix(0,qcum[nrand+1L],nrand)
        colnames(X_lambda) <- rep("(Intercept)",nrand)
        for (i in seq(nrand)) {
          if (i==1) {
            X_lambda[1:q[i],i]<-1 ## design X_lambda is a column of 1
          } else { ## the design X_lambda has two columns, each containing a block of 0's and one of 1's
            ## so that the two lambda estimates are obtained through one fit
            X_lambda[(qcum[i]+1L):qcum[i+1L],i]<-1L
          }
        }
      }
    }
    processed$X_lambda <- X_lambda
    #
    if ( is.null(phi.Fix)) {
         formulaDisp <- resid.predictor$formula
         fr_disp <- HLframes(formula=formulaDisp,data=data) 
         X_disp <- fr_disp$X
         namesX_disp <- names(fr_disp$fixef)
         colnames(X_disp) <- namesX_disp
         random_dispersion<-findbarsMM(formulaDisp) ## random effect in mean predictor of dispersion phi
         if (!is.null(random_dispersion)) {
           FL_disp <- spMMFactorList(formulaDisp, fr_disp$mf, 0L, 0L)
           namesRE_disp <- FL_disp$namesRE
           Z_disp <- FL_disp$Design
           models[["phi"]] <- "phiHGLM"
           mess <- pastefrom("LIKELY missing code to handle random effect for linear predictor for phi.")
           stop(mess)
         } else {
           if (length(namesX_disp)==1 && namesX_disp[1]=="(Intercept)") {
              models[["phi"]] <- "phiScal"
           } else { 
             models[["phi"]] <- "phiGLM"
           }
         } 
      processed$X_disp <- X_disp
    } 
    #
    Offset<-predictor$offset
    if (is.null(Offset)) {
       off <- rep(0,nrow(data)) ## model.frame.default(formula = locform, offset = off, drop.unused.levels = TRUE) expects a vector...
    } else {off <- Offset}
    processed$off <- off
    #
    processed$models <- models
    if (class(rand.family)=="family") { ## then check what is feasible
          RandDist<-tolower(rand.family$family) ## tolower once and for all
          oklink <- F
          if (RandDist=="gaussian" && rand.family$link=="identity") oklink <- T          
          if (RandDist=="gamma" && rand.family$link=="log") oklink <- T          
          if (RandDist=="inverse.gamma" && rand.family$link=="log") oklink <- T
          if (RandDist=="inverse.gamma" && rand.family$link=="-1/mu") oklink <- T
          if (RandDist=="beta" && rand.family$link=="logit") oklink <- T
          if ( ! oklink) stop("(!) Link not handled for rand.family")
    } else { ## rand.family is a string ## for private use only...
          RandDist<-tolower(rand.family) ## tolower once and for all
          rand.family <- switch(RandDist,
            gaussian = gaussian(),
            gamma = Gamma(link=log), ## NOT the default link
            beta = Beta(), ## as defined in hglm package
            "inverse.gamma" = inverse.Gamma(),
            stop("rand.family argument not valid")
          )
    }
    processed$RandDist <- RandDist
    class(processed) <- c("arglist","list")
    return(processed)
}
