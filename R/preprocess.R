checkRandLink <- function(rand.family) {
  if (class(rand.family)=="family") { ## then check what is feasible
    lcrandfamfam <- tolower(rand.family$family) ## tolower once and for all
    oklink <- F
    ## cases where g(u)=th(u)
    if (lcrandfamfam=="gaussian" && rand.family$link=="identity") oklink <- T          
    if (lcrandfamfam=="gamma" && rand.family$link=="log") oklink <- T          
    if (lcrandfamfam=="inverse.gamma" && rand.family$link=="-1/mu") oklink <- T
    if (lcrandfamfam=="beta" && rand.family$link=="logit") oklink <- T
    ## cases where g(u)!=th(u)
    if (lcrandfamfam=="inverse.gamma" && rand.family$link=="log") oklink <- T 
    if ( ! oklink) {
      allowed <- switch(lcrandfamfam,
                        gaussian= "is 'identity'",
                        gamma= "is 'log'",
                        beta= "is 'logit'",
                        "inverse.gamma" = "are '-1/mu' and 'log'"
                        )
      mess <- paste("(!) rand.family/link combination not handled;\nallowed link(s) for rand.family '",rand.family$family,"' ",allowed,sep="")
      stop(mess)
    }
  } else { ## rand.family is a string ## preprocess -> checkRandLinkS -> here : OK
    lcrandfamfam<-tolower(rand.family) ## tolower once and for all
    rand.family <- switch(lcrandfamfam,
                          gaussian = gaussian(),
                          gamma = Gamma(link=log), ## NOT the default link
                          beta = Beta(), 
                          "inverse.gamma" = inverse.Gamma(),
                          stop("rand.family argument not valid")
    )
  }
  lcrandfamfam
}

checkRespFam <- function(family) {
  ## four lines from glm()
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  family
}

checkRandLinkS <- function(rand.families) {
  unlist(lapply(rand.families,checkRandLink)) ## a vector of lcrandfamfam := tolower(rand.family$family)
}

getProcessed <- function(object,element,idx=1) {
  if ( ! is.null(attr(object,"multiple"))) object <- object[[idx]]  
  eval(parse(text=paste("object$",element,sep="")))
}

setProcessed <- function(object,element,value=1) {
  if ( ! is.null(attr(object,"multiple"))) {
    for (nam in names(object)) {
      eval(parse(text=paste("object[[",nam,"]]$",element," <- ",value,sep="")))
    }
  } else eval(parse(text=paste("object$",element," <- ",value,sep="")))
  return(object)
}

################# translate HLfit.args from user input to directly manipulated variables
## if called without non(default=NULL) phi.Fix (ie outside HLFit) => processes formulaDisp
## if called with nonNULL phi.Fix, then no need to process formulaDisp
preprocess <- function(control.HLfit,phi.Fix=NULL,HLmethod,predictor,resid.predictor,REMLformula,data,family,
                       BinomialDen,rand.families) {
  callargs <- match.call() ## to make easy to change these arguments in a later fit
  #####################################################################
  if (inherits(data,"list")) {
    locargs <- as.list(callargs)
    famfam <- family$family
    processed <- lapply(data,function(dd) {
      locargs$data <- dd
      if ( ! is.null(famfam) && famfam=="multi") locargs$family <- family$binfamily  
      eval(as.call(locargs))
    })
    attr(processed,"multiple") <- TRUE ## but NULL otherwise hence test it as not null  
    return(processed)
  }
  #####################################################################
  processed <- list()
  stop.on.error <-control.HLfit$stop.on.error ##  
  if (is.null(stop.on.error)) stop.on.error <- FALSE
  processed$stop.on.error <- stop.on.error ##  
  AIC <-control.HLfit$AIC ##  
  if (is.null(AIC)) AIC <- FALSE
  processed$AIC <- AIC ##  
  resid.family <-control.HLfit$resid.family ##  
  if (is.null(resid.family)) {
    resid.family <- Gamma(log)
  } else if (resid.family$family!= "Gamma") stop("resid.family must be Gamma.")
  processed$resid.family <- resid.family ##  
  essai <-control.HLfit$essai ##  
  if (is.null(essai)) essai <- FALSE
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
  if (is.null(iter.mean.dispFix)) iter.mean.dispFix <- control.HLfit$max.iter.mean ## public control
  if (is.null(iter.mean.dispFix)) iter.mean.dispFix <- 200 ## control of inner loop when no disp param is estimated ## was 40, 06/2014
  processed$iter.mean.dispFix <- iter.mean.dispFix  
  #
  iter.mean.dispVar <-control.HLfit$iter.mean.dispVar ## private control
  if (is.null(iter.mean.dispVar)) iter.mean.dispVar <- control.HLfit$max.iter.mean ## public control ## control of inner loop when some disp param is estimated
  if (is.null(iter.mean.dispVar)) iter.mean.dispVar <- 50 ## control of inner loop when some disp param is estimated  ## was 20, 06/2014
  processed$iter.mean.dispVar <- iter.mean.dispVar  
  #
  max.iter <- control.HLfit$max.iter  ## control of outer loop 
  if (is.null(max.iter)) max.iter <-200
  processed$max.iter <- max.iter  
  #
  ## only check, no modif of data -> not returned
  if ( ! (inherits(data,"data.frame") || inherits(data,"environment"))) {
    mess <- pastefrom("'data' is not a data.frame not an environment.",prefix="(!) From ")
    stop(mess)
  }
  if (! "predictor" %in% class(predictor)) predictor <- Predictor(predictor) 
  if (! "predictor" %in% class(resid.predictor)) resid.predictor <- Predictor(resid.predictor)
  MeanFrames <- HLframes(formula=predictor,data=data) ## design matrix X, Y, fixef names
  y <- MeanFrames$Y
  if (var(y)==0) {
    stop("(!) var(response) = 0.")
  }  
  X.pv <- MeanFrames$X
  off <- model.offset(MeanFrames$mf) ## look for offset from (ori)Formula 
  if ( is.null(off) ) { ## no offset (ori)Formula term. Check attribute (Predictor ensures offset is not both in formula and attr)
    off <- attr(predictor,"offsetObj")$vector 
  } else {
    keepOriForm <- attr(predictor,"oriFormula")
    predictor <- noOffset(predictor)
    attr(predictor,"oriFormula") <- keepOriForm
  }
  if (!is.null(off)) {
    off <- pmax(log(.Machine$double.xmin),off) ## handling log(0) ## but if input off were NULL, output off would be is numeric(0) where it should remain NULL
  }
  roff <- model.offset(model.frame(resid.predictor,data=data)) ## look for offset from (ori)Formula 
  if ( is.null(roff) ) { ## ## no offset (ori)Formula term. Check attribute (Predictor ensures offset is not both in formula and attr)
    roff <- attr(resid.predictor,"offsetObj")$vector 
  } else {
    keepOriForm <- attr(resid.predictor,"oriFormula")
    resid.predictor <- noOffset(resid.predictor)
    attr(resid.predictor,"oriFormula") <- keepOriForm
  }
  colnames(X.pv) <- names(MeanFrames$fixef)
  processed$X.pv <- X.pv ## further modified in HLfit...
  ## in particular, $mf contains values of expl.vars and levels (indices) of ranefs for all observations. Indices are/should be ordered as uniqueGeo in ULI()
  nobs <- NROW(MeanFrames$Y)
  # MODIFICATION of formula and offset (but not oriFormula)     
  if (is.null(off)) { ## model.frame.default(formula = locform, offset = off,...) expects a vector....
    attr(predictor,"offsetObj") <- list(vector=rep(0,nobs),nonZeroInfo=FALSE)
  } else {
    attr(predictor,"offsetObj") <- list(vector=off,nonZeroInfo=TRUE) ## subtly, will be of length 1 if   original offset was a constant...
  }
  processed$predictor <- predictor  
  #
  if (is.null(roff)) { ## model.frame.default(formula = locform, offset = off,...) expects a vector....
    attr(resid.predictor,"offsetObj") <- list(vector=rep(0,nobs),nonZeroInfo=FALSE)
  } else {
    attr(resid.predictor,"offsetObj") <- list(vector=roff,nonZeroInfo=TRUE) ## subtly, will be of length 1 if   original offset was a constant...
  }
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
    BinDenForm <- attr(predictor,"BinDenForm")
    if ( ! is.null(BinDenForm)) { ## the cbind syntax was used in the formula 
      BinomialDen <- eval(parse(text=BinDenForm),envir=data)
      ## la suite ducode suppose que pas cbind => ie essentially obsolete syntax 
    } else if (missing(BinomialDen) || is.null(BinomialDen)) { ## then this should be a binary response
      checkResp <- eval(parse(text=as.character(predictor[2])),envir=data) ## 'positives'
      if (any(checkResp>1)) {
        mess <- pastefrom("binomial, non-binary response. Please use the",prefix="(!) From ")
        message(mess)
        message("    standard glm() syntax with _cbind_: 'cbind(<successes>, <failures>) ~ <etc>'")
        stop()
      } else BinomialDen <- rep(1,nobs) ## response appears to be binary...
    }
    no.info <- (BinomialDen == 0)
    if (any(no.info)) {
      mess <- pastefrom("please remove missing data (i.e. for which binomial sample size is 0).",prefix="(!) From ")
      stop(mess)
    }
    ## It's not really possible to remove data at this stage as this may not match the dimension of the distance matrices
    ## moreover one cannot simply remove rows of a matrix "root"...
    ## it _could_ be useful to be able to hand BinomilaDen=0 by the general code but...
  } else {
    BinomialDen <- rep(1,nobs)
  }
  processed$BinomialDen <- BinomialDen  
  ## conversion from user-friendly format to standard 'XX(...)' format
  ## first index is for (0) h / (1) p_v(h) / (2) p^s_v(h) ie whether h lik or marginal lik is used for fixed effect estimation
  ## Other components determine three options wrt to leverages, some for ReML correction, some for notEQL.
  ## ML() vs HL() determines whether the hatval leverages are computed ie whether some basic ReML correction in applied
  ## second index is for further ReML correction: no (0) or yes (1) D log Det / d log lambda correction (2) further p^s_bv correction
  ## third index is for use of (0) EQL deviance residuals (this affects the leverages) or not (1) (this is not an ReML correction... but impatcs only dispersion estimation..)
  ## thus overall we have <ReML/not>( <h/l> , <more ReML/not> , <not/EQL> )
  ## NohL07 table 1 has interesting terminology and further tables show even more cases
  terms_ranefs <- parseBars(predictor) ## a vector of char strings
  nbars <- length(terms_ranefs) ## not nbars != nrand because nested effects will be further expanded
  models <- list(eta="",lambda="",phi="")
  if (nbars>0) {
    processed$lambdaFamily <- Gamma(link="log")
    models[["eta"]] <- "etaHGLM" 
    FL <- spMMFactorList(predictor, MeanFrames$mf, 0L, drop=TRUE) ## this uses the spatial information in the formula, even if an explicit distMatrix was used elsewhere
    ZAlist <- FL$Design ## : is a list of design matrices (temporarily only Z)
    attr(ZAlist,"ranefs") <- terms_ranefs
    attr(ZAlist,"Groupings") <- FL$Groupings
    attr(ZAlist,"namesTerms") <- FL$namesTerms
    AMatrix <- attr(predictor,"AMatrix")
    if (!is.null(AMatrix)) {
      ## logic is Z[nresp,nUniqueRespLoc].A[nUniqueRespLoc,nHiddenv].L[nHiddenv,nHiddenv]
      for (iMat in seq(length(ZAlist))) {
        ZAlist[[iMat]] <- ZAlist[[iMat]] %*% AMatrix[[iMat]]
      }
    }
  } else {
    models[["eta"]] <- "etaGLM" 
    ZAlist <- NULL
  }
  processed$ZAlist <- ZAlist
  nrand <- length(ZAlist)
  #
  if (inherits(rand.families,"family")) rand.families <- list(rand.families) ## I should pass rand.families to preprocess
  if (nrand != 1 && length(rand.families)==1) rand.families <- rep(rand.families,nrand) 
  lcrandfamfam <- checkRandLinkS(rand.families)  
  if (HLmethod=="ML") {
    HLmethod <- "ML(1,1,1)" ## here there could be a special hack for (family$family=="binomial" && mean(BinomialDen)<HL2.threshold) 
  } else if (HLmethod=="SEM") {
    HLmethod <- "ML('SEM',NA,NA)"  
    SEMseed <-control.HLfit$SEMseed ##  
    if (is.null(SEMseed)) SEMseed <- 123 ## OK pour SEM *unique* mais remplace par NULL dans locoptimthroughSmooth
    processed$SEMseed <- SEMseed ##
    nMCint <-control.HLfit$nMCint ##  
    #      if (is.null(nMCint)) nMCint <- 10000
    processed$nMCint <- nMCint ##
    nSEMiter <-control.HLfit$nSEMiter ##  
    #      if (is.null(nSEMiter)) nSEMiter <- 100
    processed$nSEMiter <- nSEMiter ##
    ngibbs <-control.HLfit$ngibbs ##  
    #      if (is.null(ngibbs)) ngibbs <- 20
    processed$ngibbs <- ngibbs ##
    processed$SAEMsample <- control.HLfit$SAEMsample ## stays NULL if NULL
  } else if (HLmethod=="REML") {
    HLmethod <- "HL(1,1,1)" ## here there could be a special hack for (family$family=="binomial" && mean(BinomialDen)<HL2.threshold) 
  } else if (HLmethod=="REPQL" || HLmethod=="PQL") {
    if (any(lcrandfamfam!="gaussian"))  stop("PQL is not defined for HGLMs in general. Do you mean 'EQL-'?") 
    HLmethod <- "HL(0,0,1)" ## (=REPQL, equivalent to HL(0,0,0) ie EQL- for GLMMs )
  } else if (HLmethod=="PQL/L") { ## again no D log Det / d log lambda correction
    if (any(lcrandfamfam!="gaussian"))  stop("PQL is not defined for HGLMs in general. Do you mean 'EQL-'?") 
    HLmethod <- "ML(0,0,1)" ## (equivalent to ML(0,0,0) for GLMMs)
  } else if (HLmethod=="EQL-") { ## version LeeNP06 p.212 incomplete 1st order ## probably hglm package
    ## thus overall we have <ReML->HL >( <h->0> , <not more ReML->0> , <EQL -> 0> )
    HLmethod <- "HL(0,0,0)" ## 
  } else if (HLmethod=="EQL+") { ## version LeeN01 complete 1st order
    ## thus overall we have <ReML->HL >( <h->0> , <more ReML->1> , <EQL -> 0> )
    HLmethod <- "HL(0,1,0)" ## (0,...): gaussianise everything, hence no a(1) correction ## there is no HL(1) a(1) correction in GLM.MME
  }
  HL <- eval(parse(text=paste("c",substr(HLmethod,3,100),sep=""))) ## extracts the (...) part into a vector
  if (length(HL)==2) HL <- c(HL,1)
  ## HLmethod is not a member of the return object
  processed$HL <- HL  
  if (substr(HLmethod,0,2)=="ML" && HL[1]!="SEM") { ## FR->FR c'est bizarre d'exclure le SEM là... p_bv est il vraiment utilisé ?
    if ( ! is.null(REMLformula)) {
      message("Confusing combination of arguments: 'HLmethod=ML(...)' with non-null 'REMLformula'.")
      stop("  Make sure what you mean and simplify the arguments.")
    }
    lhs <- paste(predictor)[[2]] ## extract response, either cbind or not
    if ( ! is.null(terms_ranefs)) {
      ## build formula with only random effects
      REMLformula <- as.formula(paste(lhs,"~",paste(terms_ranefs,collapse="+")))
    } else REMLformula <- as.formula(paste(lhs,"~ 0")) 
  }
  processed$REMLformula <- REMLformula  
  #
  processed$loglfn.fix <- selectLoglfn(family$family)
  if ( family$family %in% c("binomial","poisson")) {
    ## the response variable should always be Counts
    if (max(abs(y-as.integer(y)))>1e-05) {
      mess <- pastefrom("response variable should be integral values.",prefix="(!) From ")
      stop(mess)
    }
    if (DEPARSE(resid.predictor) != "~1") {
      warning("resid.formula is ignored in binomial- or Poisson-response models")
    }
  }
  if ( family$family == "binomial" && ncol(y)==2) y <- y[,1,drop=F] ## that is, we have the cbind syntax up to this fn
  ## code derived from the glm() function in the safeBinaryRegression package
  #    if(family$family == "binomial" && length(unique(y)) == 2 && require(lpSolveAPI)) {
  if(family$family == "binomial" && length(unique(y)) == 2 && ncol(X.pv)>0) {
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
  if (nrand>0 && family$family=="gaussian" && all(lcrandfamfam=="gaussian") && canonicalLink) {
    LMMbool <- TRUE
    ## identifiability checks cf modular.R -> checkNlevels() in lmer:
    ## FR->FR one could also add an isNested check as in https://github.com/lme4/lme4/blob/master/R/utilities.R, but presumably not here
    LMatrix <- attr(predictor,"LMatrix")
    for (iMat in seq(length(ZAlist))) {
      nc <- ncol(ZAlist[[iMat]])
      if (nc < 2 && is.null(phi.Fix)) {
        mess <- paste("Only ",nc," level for random effect ",terms_ranefs[iMat],
                      ";\n   this model cannot be fit unless phi is fixed.",sep="")
        warning(mess)
      }
      if (is.null(LMatrix) && substr(terms_ranefs[iMat], 1, 1)=="(" && nc == nobs && is.null(phi.Fix)) {
        mess <- paste("Number of levels = number of observations \n   for random effect ",
                      terms_ranefs[iMat],
                      ";\n   this model cannot be fit unless phi is fixed\n   or a correlation matrix is given.",sep="")
        stop(mess)
      }
    }
  } else LMMbool <- FALSE
  processed$LMMbool <- LMMbool
  ## algorithms (control of defaults remains in the HLfit code)
  betaFirst <-control.HLfit$betaFirst ##  
  if (is.null(betaFirst)) {
    betaFirst <- FALSE
  } else if (betaFirst && HL[1]=="SEM") {
    message("betaFirst && HLmethod= SEM: betaFirst turned to FALSE")
    betaFirst <- FALSE
  }
  processed$betaFirst <- betaFirst ##
  LevenbergM <- control.HLfit$LevenbergM ## 
  if (is.null(LevenbergM)) { 
    if (HL[1]=="SEM") {
      LevenbergM <- FALSE 
    } else if (betaFirst) {
      LevenbergM <- FALSE 
    } else {
      if (LMMbool) {  
        LevenbergM <- FALSE ## because no reweighting when beta_eta changes => no IRWLS necess   
      } else LevenbergM <- TRUE
    }
  }
  processed$LevenbergM <- LevenbergM
  #
  if (nrand>0) {   
    nrand_lambda <- 0
    models[["lambda"]] <- FL$termsModels
    ################################################################################
    # for a random slope term X.v in X(beta+ v), the X went into the general ZAL matrix 
    # (construction of ZAlist by spMMFactorList), and
    # we are still estimating the lambda's using a X_lamres with 0/1's only
    # unless there is a non-trivial model for the lambdas
    ################################################################################
    if (all(models[["lambda"]]=="lamScal")) { 
      Xi_cols <- unlist(lapply(FL$namesTerms,length))
      if (any(Xi_cols>1 & !lcrandfamfam=="gaussian")) {
        stop("(!) random slope models with correlated non-gaussian random effects are not fitted.")
      }
      cum_Xi_cols <- cumsum(c(0, Xi_cols)) ## if two ranef,  with n_u_h=(3,3), this is 0,3,6. cum_h_u_h[nrand+1] is then 6, the total # of realizations
      n_u_h <- rep(0, sum(Xi_cols))
      for (i in 1:nrand) n_u_h[(cum_Xi_cols[i]+1L):cum_Xi_cols[i+1L]] <-  nlevels(FL$Subject[[i]]) 
      cum_h_u_h <- cumsum(c(0, n_u_h)) ## if two "Intercept" ranefs,  with n_u_h=(3,3), this is 0,3,6. cum_h_u_h[nrand+1] is then 6, the total # of realizations
        ## if (1+X|...) +(1|...),  with n_u_h=(3,4), this is 0,3,6,10. cum_h_u_h[sum(Xi_cols)+1] is then 10, the total # of realizations
      X_lamres <- matrix(0,cum_h_u_h[sum(Xi_cols)+1L],sum(Xi_cols))
      colnames(X_lamres) <- unlist(FL$namesTerms)
      for (i in seq(nrand)) {
        for (j in (cum_Xi_cols[i]+1L):cum_Xi_cols[i+1L]) {
          X_lamres[(cum_h_u_h[j]+1L):cum_h_u_h[j+1L],j] <- 1L ## this maps the deviance residuals to the lambda's to be estimated from them
        }
      } 
      attr(X_lamres,"Xi_cols") <- Xi_cols ## will provide info for covariance computations
      attr(X_lamres,"cum_Xi_cols") <- cum_Xi_cols ## will provide info for covariance computations
    } else {  ## linear predictor for variance of random effects (lambda) (lamGLM or lamHGLM) 
      if (any(models[["lambda"]]=="lamHGLM")) { ##need distinct calls... to fit each lambda model  
        if (length(formulaLambda)==2) formulaLambda <- as.formula(paste('"lambda"',paste(formulaLambda,collapse=" ")))
        lambda_ranefs <- findbarsMM(formulaLambda)
        if (!is.null(lambda_ranefs)) {  ## lamHGLM
          nrand_lambda <- length(lambda_ranefs)
          models[["lambda"]] <- list("lamHGLM")
        } else models[["lambda"]] <- list("lamGLM")  
        colnames(X_lambda) <- names(fr_lambda$fixef)
      } else { ## can use a single design matrix for all random effects, which is convenient.
        mess <- pastefrom("LIKELY missing code to handle linear predictor for lambda.")
        stop(mess)
        # la suite c'est dexu residus de code a assembler: il faut une liste de HLframes du type
        fr_lambda <- HLframes(formula=formulaLambda,data=data) ## but the "data" should probably be distinct data here, with nrow=number of reals of ranefs 
        # (pobablement calculee en amont pour determiner lamScal aussi...) ensuite extraire les design matrices
        
      }
    } 
  } else X_lamres <- NULL
  processed$X_lamres <- X_lamres
  #
  if ( is.null(phi.Fix)) {
    formulaDisp <- resid.predictor
    fr_disp <- HLframes(formula=formulaDisp,data=data) 
    X_disp <- fr_disp$X
    ## if formula= ~1 and data is an environment, there is no info about nobs, => X_disp has zero rows, which is a problem later 
    if(nrow(X_disp)==0) X_disp=matrix(1,nrow=nobs)
    namesX_disp <- names(fr_disp$fixef)
    colnames(X_disp) <- namesX_disp
    random_dispersion<-findbarsMM(formulaDisp) ## random effect in mean predictor of dispersion phi
    if (!is.null(random_dispersion)) {
      FL_disp <- spMMFactorList(formulaDisp, fr_disp$mf, 0L, drop=TRUE)
      namesRE_disp <- FL_disp$Groupings
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
  processed$models <- models
  attr(rand.families,"lcrandfamfam") <- lcrandfamfam
  processed$rand.families <- rand.families
  processed$callargs <- callargs
  class(processed) <- c("arglist","list")
  return(processed)
}

eval.update.call <- function(mc,...) {
  mc <- as.list(mc)
  dotlist <- list(...)
  # mc[names(dotlist)] <- dotlist
  eval(as.call(mc))  
}
