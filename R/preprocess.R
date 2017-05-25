checkRandLink <- function(rand.family) {
  lcrandfamfam <- tolower(rand.family$family) ## tolower once and for all
  oklink <- F
  ## cases where g(u)=th(u)
  if (lcrandfamfam=="gaussian" && rand.family$link=="identity") oklink <- T          
  if (lcrandfamfam=="gamma" && rand.family$link=="log") oklink <- T          
  if (lcrandfamfam=="inverse.gamma" && rand.family$link=="-1/mu") oklink <- T
  if (lcrandfamfam=="beta" && rand.family$link=="logit") oklink <- T
  ## cases where g(u)!=th(u)
  if (lcrandfamfam=="inverse.gamma" && rand.family$link=="log") oklink <- T 
  if (lcrandfamfam=="gamma" && rand.family$link=="identity") oklink <- T ## gamma(identity)
  if ( ! oklink) {
    allowed <- switch(lcrandfamfam,
                      gaussian= "is 'identity'",
                      gamma= "is 'log' ('identity' may be *tried*)", ## gamma(identity) not yet working
                      #  and the above code protext against Gamma(inverse); it could perhaps be allowed for trying ? 
                      beta= "is 'logit'",
                      "inverse.gamma" = "are '-1/mu' and 'log'"
    )
    mess <- paste("(!) rand.family/link combination not handled;\nallowed link(s) for rand.family '",rand.family$family,"' ",allowed,sep="")
    stop(mess)
  }
  lcrandfamfam
}


checkRandLinkS <- function(rand.families) {
  rand.families <- lapply(rand.families, function(rf) {
    if (is.character(rf)) {
      rf <-switch(tolower(rf),
                  gaussian = gaussian(),
                  gamma = Gamma(link=log), ## NOT the default link
                  beta = Beta(), 
                  "inverse.gamma" = inverse.Gamma(),
                  stop("rand.family argument not valid"))
    }
    return(rf) 
  })
  lcrandfamfam <- unlist(lapply(rand.families,checkRandLink)) ## a _vector_ of lcrandfamfam := tolower(rand.family$family)
  return(list(rand.families=rand.families,lcrandfamfam=lcrandfamfam))
}

.checkRespFam <- function(family) {
  envname <- environmentName(environment(family))
  if ( ! envname %in% c("stats","spaMM","")) { ## " typically occurs when family has laready been checked...
    message(paste("family from namespace environment",envname,"possibly not correctly handled"))
    if (envname=="mgcv" && identical(paste(formals(family)$theta)[1],"stop")) { ## a call via family=family prevents testing the name
      mess <- "spaMM::negbin is masked by mgcv::negbin. Use explicitly 'family=spaMM::negbin', or assign 'negbin <- spaMM::negbin'"
      stop(mess)
    }
  }
  ## four lines from glm(), which should have no effect on processed$family, spec. those with a param.
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  # if (is.name(family)) family <- eval(family)
  if (family$family=="Gamma") {
    family <- spaMM_Gamma(link=family$link) 
  }
  return(family) ## input negbin(...) or COMPoisson(...) are returned as is => nu/shape unaffected
}

# to convert back a family object to a call, EXCEPT for "multi". 
# get_param should be TRUE only in cases we are sure that the param has a value
.as_call_family <- function(family,get_param=FALSE) {
  # checkRespFam should have been run. 
  # Then we have either a call or an evaluated family object
  if (identical(family[[1L]],"multi")) {
    return(family) ## FR->FR tempo fix
  } else {
    as_call <- environment(family$aic)$mc ## shouldn't be called on multi...
    if (is.null(as_call)) { ## then we have a stats::  family 
      as_call <- call(family$family, link=family$link) ## need a substitute() ? 
    } else if (get_param) { ## param values are constant or dynamically assigned to processed$family
      famfam <- family$family 
      if (famfam=="negbin") {
        as_call$shape <- environment(family$aic)$shape
      } else if (famfam=="COMPoisson") {
        as_call$nu <- environment(family$aic)$nu
      }
    }
    return(as_call)
    ## avoids display of family def in fit results
  }
}

get_ZAlist <- function(predictor, model_frame) { ## model_frame is an $mf element
  terms_ranefs <- parseBars(predictor) ## a vector of char strings with attribute(s), 
  ##    otherwise similar to bars <- .spMMexpandSlash(findbarsMM(formula[[length(formula)]])): FR->FR maybe simplification of code possible here?
  if ( length(terms_ranefs) > 0L ) {
    FL <- .spMMFactorList(predictor, model_frame, 0L, drop=TRUE) ## this uses the spatial information in the formula, even if an explicit distMatrix was used elsewhere
    ZAlist <- FL$Design ## : is a list of design matrices (temporarily only Z)
    attr(ZAlist,"ranefs") <- terms_ranefs
    ## FR->FR ca serait utile d'appeler l'attribut 'terms_ranefs' mais predict() etc deviendrait 
    ##   back incompat avec anciens fits. 
    ## Toutefois rend visible le pb que les attr 'ranefs ' de LMatrix et de ZAlist paraissent calcules
    ##   de deux facons differentes (parseBars() vs unlist(lapply(spatial.terms,.DEPARSE)) ) 
    attr(ZAlist,"namesTerms") <- FL$namesTerms ## list of predictor vars for each element of ZAList
    AMatrix <- attr(predictor,"AMatrix")
    if (!is.null(AMatrix)) {
      ## logic is Z[nresp,nUniqueRespLoc].A[nUniqueRespLoc,nHiddenv].L[nHiddenv,nHiddenv]
      for (iMat in seq(length(ZAlist))) {
        ZAlist[[iMat]] <- ZAlist[[iMat]] %*% AMatrix[[iMat]]
      }
    }
  } else ZAlist <- structure(list(),anyRandomSlope=FALSE)
  return(ZAlist)
}

getProcessed <- function(object,element,from=NULL) {
  if ( ! is.null(attr(object,"multiple"))) {
    if (is.null(from)) from <- seq_len(length(object))
    if (length(from)>1L) {
      resu <- lapply(from,function(id) {eval(parse(text=paste("object[[",id,"]]$",element,sep="")))})
      names(resu) <- names(object[from])
      return(resu) ## a list, e.g a data list
    } else {
      object <- object[[from]]
    } 
  } 
  return(eval(parse(text=paste("object$",element,sep=""))))
}


setProcessed <- function(object,element,value=1) {
  if ( ! is.null(attr(object,"multiple"))) {
    #    for (it in seq_len(length((object))) {eval(parse(text=paste("object[[",it,"]]$",element," <- ",value,sep="")))} 
    ## fails bc nam loses its enclosing "s :
    #for (nam in names(object)) {eval(parse(text=paste("object[[",nam,"]]$",element," <- ",value,sep="")))} 
    for (nam in names(object)) {eval(parse(text=paste("object[[\"",nam,"\"]]$",element," <- ",value,sep="")))} 
  } else eval(parse(text=paste("object$",element," <- ",value,sep="")))
  return(object)
}

generateInitPhi <- function(formula,data,family,weights=NULL) {
  ## look for replicates to estimate phi
  form <- subbarsMM(asStandardFormula(formula)) ## removes all random effect decorum (but retains its variables)
  # lhs
  mf <- model.frame(form,data=data)
  Y <- model.response(mf) ## evaluated rhs (e.g. log(y) rather than y...)
  # rhs
  rhs_terms <- stats::delete.response(terms(form))
  mf <- model.frame(rhs_terms, data)
  # builds a model which indexes responses with identical predictor [ULI(mf)] 
  # and retains data that are replicates for each level of this index [uli]
  # but (g)lm complains about a model where uli has a single level [though this is meaningful]
  uli <- ULI(mf) 
  # selection of data for replicates
  mf <- data.frame(y=Y,uli=as.factor(uli)) ## what's needed for both sides
  tuli <- table(uli)
  phiinfo <- which(tuli>1L) ## which levels of uli are informative
  ############################################
  # -> may generate a large # of levels -> of parameters to estimate -> glm takes time ! ## FR->FR may take time => quick patch    
  if (length(phiinfo)>30) return(0.1)  ## RETURN !
  ############################################
  whichRows <- uli %in% phiinfo
  mf <- mf[whichRows,] ## keep lines with replicates of all predictor variables
  if (!is.null(weights)) weights <- weights[whichRows]
  if (length(phiinfo)>1L) {
    locform <- y ~ uli 
  } else locform <- y ~ 1 ## only one level of uli had replicates
  # estimation of residual var
  if (NROW(mf)>1L) { ## if replicate observations available
    # formula must retain any operation on the lhs
    locglm <- spaMM_glm(formula=locform,data=mf,family=family,weights=weights)  
    phi_est <- as.numeric(deviance(locglm)/locglm$df.residual)
  } else phi_est <- NULL
  return(phi_est)
}

## function for init.optim hence for fitme.
.get_init_phi <- function(processed,weights=NULL) {
  if (is.null(processed$envir$init_phi)) { 
    processed$envir$init_phi <- generateInitPhi(formula=processed$predictor,data=processed$data,
                                family=processed$family,weights=weights) 
  } 
  return(processed$envir$init_phi)
}

.get_init_lambda <- function(processed,reset=FALSE,stillNAs,init_lambda_by_glm=NULL) {
  if (is.null(processed$envir$init_lambda)) { ## which occurs either if reset is TRUE or if $get_init_lambda has not yet been called
    nrand <-  length(processed$ZAlist)
    preproFix <- processed$lambda.Fix
    family <- processed$family
    rand.families <- processed$rand.families
    lcrandfamfam <- attr(rand.families,"lcrandfamfam")
    init_lambda <- rep(NA,nrand)
    
    ### (1) Generate a single value for all ranefs ## FR->FR this could  be simplified to generate missing values only
    init.lambda <- init_lambda_by_glm
    ## les tests et examples sont mieux sans la $variance, mais ce n'est pas attendu theor pour large N
    ## with variance: cf CoullA00 p. 78 top for motivation, et verif par simul:
    ## the excess rel var (given by their corr_rr' for r=r') is of the order of (their rho_rr=)lambda
    ## hence excess var (overdisp) is lambda Npq hence lambda ~ overdisp/Npq ~ overdisp/(resglm$prior.weights*family$variance(fv)) ?
    ## pas convainquants en temps:
    #init.lambda <- max(0.01,sum(((resid(resglm,type="pearson")/resglm$prior.weights)^2)/resglm$family$variance(fv)-1)/resglm$df.residual )
    #init.lambda <- max(0.01,sum((resid(resglm,type="pearson")/resglm$prior.weights)^2-resglm$family$variance(fv))/resglm$df.residual )
    ## pas loin du temps de réference ?:
    #init.lambda <- sum(resid(resglm,type="pearson")^2/(resglm$prior.weights*family$variance(fv)))/resglm$df.residual
    if (family$link=="log") { ## test of family, not rand.family...
      init.lambda <- log(1.00001+init.lambda) ## max(0.0001,log(init.lambda))
    } else init.lambda <- init.lambda/5 ## assume that most of the variance is residual
    init.lambda <- init.lambda/nrand        
    ###
    #
    ### (2) generate missing values from this unique value, each adapted to a family and link
    # allows for different rand.family
    valuesforNAs <- unlist(lapply(stillNAs, function(it) {
      if(lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="identity" && init.lambda==1) {
        adhoc <- 0.9999
      } else if(lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="log") {
        objfn <- function(lambda) {psigamma(1/lambda,1)-init.lambda}
        adhoc <- uniroot(objfn,interval=c(1e-8,1e8))$root
      } else if(lcrandfamfam[it]=="beta" && rand.families[[it]]$link=="logit") {
        #ad hoc approximation which should be quite sufficient; otherwise hypergeometric fns.
        objfn <- function(lambda) {8* lambda^2+3.2898*lambda/(1+lambda)-init.lambda}
        adhoc <- uniroot(objfn,interval=c(2.5e-6,1e8))$root
      } else if(lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="log") {
        ## this tries to controle var(v)<-init.lambda of v=log(u) by defining the right lambda for var(u)=lambda/(1-lambda)
        ## log (X~inverse.G) = - log (~Gamma), ie
        ##  log (X~inverse.G(shape=1+1/init.lambda,scale=1/init.lambda))~ - log (rgamma(shape=1+1/init.lambda,scale=1/init.lambda)
        ## where + log (~Gamma) has known mean (=...) and variance=trigamma(shape)
        ## (pi^2)/6=1.644934... is upper bound for trigamma(x->1) ie for lambda ->infty.
        ## however, there is something wrong with setting very large initial lambda values. Hence
        # if (init.lambda > 1.64491 ) { ## trigamma(1+1/100000) 
        #   adhoc <- 100000 
        if (init.lambda > 0.3949341 ) { ## trigamma(3) [for lambda=0.5]
          adhoc <- 0.5 
        } else { ## loer init.lambda -> lower lambda
          objfn <- function(lambda) {psigamma(1+1/lambda,1)-init.lambda}
          adhoc <- uniroot(objfn,interval=c(1e-8,1e8))$root
        }
        ## but the mean of v  is function of lambda and I should rather try to control the second moment of v
      } else if(lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="-1/mu") {
        adhoc <- (sqrt(1+4*init.lambda)-1)/2 # simple exact solution
      } else adhoc <- init.lambda
      adhoc
    }))
    init_lambda[stillNAs] <- valuesforNAs
    processed$envir$init_lambda <- init_lambda
  }
  return(processed$envir$init_lambda)
}


# function called within HLfit for missing inits
.get_inits_by_glm <- function(processed,
                              family=processed$family,
                              reset=FALSE) {
  if (reset || is.null(processed$envir$inits_by_glm)) { 
    y <- processed$y ## requested by the formula
    if (family$family=="binomial" && NCOL(y)==1L) { 
      ##  && ncol(y)==1: attempt to implement the cbind() for y itself syntax throughout. But fails later on 'y - mu'...
      BinomialDen <- processed$BinomialDen ## requested by the formula
      begform <-"cbind(y,BinomialDen-y)~"  
    } else {begform <-"y~"}
    ###################################################if (pforpv==0) {endform <-"0"} else 
    X.pv <- processed$`X.pv` ## possibly requested by the formula
    pforpv <- ncol(X.pv)
    if(pforpv>0) {
      endform <-"X.pv-1" ## pas besoin de rajouter une constante vue qu'elle est deja dans X
    } else {
      if (family$family %in% c("binomial","poisson")) {
        endform <- "1" ## no meaningful glm without fixed effect in this case !
      } else {endform <- "0"}
    }
    locform <- as.formula(paste(begform, endform))
    off <- attr(processed$predictor,"offsetObj")$total
    prior.weights <- processed$prior.weights   
    if (is.call(prior.weights)) {
      resglm <- spaMM_glm(locform,family=family,offset=off,
                          weights=NULL,control=processed[["control.glm"]])
    } else resglm <- spaMM_glm(locform,family=family,offset=off,
                               weights=eval(prior.weights),control=processed[["control.glm"]])
    resu <- list(glm=resglm,phi_est=as.numeric(deviance(resglm)/resglm$df.residual))
    if (family$family=="binomial" && max(resglm$prior.weights)==1L) { ## binary response
      resu$lambda <- 1
    } else {
      fv <- fitted(resglm)
      resu$lambda <- sum(resid(resglm)^2/(resglm$prior.weights*family$variance(fv)))/resglm$df.residual
    }
    if (pforpv>0) {
      ## Two potential problems (1) NA's pour param non estimables (cas normal); 
      ## (2) "glm.fit: fitted probabilities numerically 0 or 1 occurred" which implies separation or large offset
      if (max(abs(c(coefficients(resglm))),na.rm=TRUE)>1e10) { ## na.rm v1.2 
        message("(!) Apparent divergence of estimates in a *GLM* analysis of the data.")
        message("    Check your data for extreme values, separation or bad offset values.")
        stop("    I exit.") 
      } 
      beta_eta <- c(coefficients(resglm)) ## this may include NA's. Testcase: HLfit(Strength ~ Material*Preheating+Method,data=weld)
      if (all(names(beta_eta)=="X.pv")) { ## si la formula etait y ~X.pv-1
        names(beta_eta) <- colnames(resglm$model$X.pv)
      } else names(beta_eta) <- unlist(lapply(names(beta_eta),substring,first=5)) ## removes "X.pv" without guessing any order or length
      resu$beta_eta <- beta_eta
    } 
    #parent.env(environment(processed$get_inits_by_glm)) <- environment(stats::glm)
    processed$envir$inits_by_glm <- resu
  } 
  return(processed$envir$inits_by_glm)
} 


preprocess <- function(control.HLfit, ranFix=NULL, HLmethod, 
                       predictor, resid.model,
                       REMLformula, data, family,
                       BinomialDen, rand.families, etaFix, prior.weights,
                       objective=NULL,
                       control.glm, adjMatrix=NULL
                       ) {
  callargs <- match.call() 
  #
  if ( ! is.list(resid.model)) resid.model <- list(formula=resid.model,family=control.HLfit$resid.family)
  resid.model$fixed <- as.list(resid.model$fixed) ## converts NULL to list() as exp'd for 'fixed' in fitme_body()
  ################ handling list of data #######################
  if (inherits(data,"list")) {
    locargs <- as.list(callargs)
    famfam <- family$family
    processed <- lapply(data,function(dd) {
      locargs$data <- dd
      if ( ! is.null(famfam) && famfam=="multi") locargs$family <- family$binfamily  
      eval(as.call(locargs)) ## call("preprocess",...) on each data set
    })
    attr(processed,"multiple") <- TRUE ## but NULL otherwise hence test it as not null  
    return(processed)
  }
  ###############################################################
  
  # remove rows with NA's in required variables:
  validdata <- getValidData(formula=predictor, resid.formula=resid.model$formula,
                            data=data, callargs=callargs["prior.weights"]) ## OK for all prior weights
  if ( inherits(data,"data.frame")) {
    data <- data[rownames(validdata),,drop=FALSE]  
  } else if  (inherits(data,"environment")) {
    data <- validdata
  } else {
    mess <- pastefrom("'data' is not a data.frame not an environment.",prefix="(!) From ")
    stop(mess)
  }
  # add easily testable family name
  famfam <- family$family
  processed <- list(data=data,family=family,
                    envir=list2env(list(),parent=environment(HLfit)))
  #
  stop.on.error <- control.HLfit$stop.on.error ##  
  if (is.null(stop.on.error)) stop.on.error <- FALSE
  processed$stop.on.error <- stop.on.error ##  
  break_conv_logL <- control.HLfit$break_conv_logL ##whether to stop if logL (p_v) appears to have converged  
  if (is.null(break_conv_logL)) break_conv_logL <- FALSE
  processed$break_conv_logL <- break_conv_logL ##  
  ## numerical control parameters 
  conv.threshold <-control.HLfit$conv.threshold ## 
  if (is.null(conv.threshold)) conv.threshold <- 1e-05 ## changing to 1e-04 (08/2016) reduces tests/testthat time by nearly 30s
  processed$conv.threshold <- conv.threshold  
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
  resid.family <- resid.model$family
  if (is.null(resid.family)) {
    resid.family <- spaMM_Gamma(log)
  } else {
    if (resid.family$family!= "Gamma") stop("resid.family must be Gamma.")## spaMM_Gamma also valid by this test
    resid.family <- spaMM_Gamma(resid.family$link) ## we will need the returned link, not the promise 
  } ## and "quoted" suitable for saving in HLfit return object
  attr(resid.family,"quoted") <- substitute(spaMM_Gamma(link),
                                            list(link=resid.family$link))
  resid.model$family <- resid.family ## resid.predictor will also be returned
  #
  if (! inherits(predictor,"predictor")) predictor <- Predictor(predictor) 
  MeanFrames <- HLframes(formula=predictor,data=data) ## design matrix X, Y...
  #
  y <- MeanFrames$Y
  if ( family$family == "binomial" && NCOL(y)==2L) y <- y[,1L,drop=FALSE] ## that is, we have the cbind syntax up to this fn
  ## binomial denom determined later, not using y => FR->FR not clear why we need two column format up to this point.
  processed$y <- y
  nobs <- NROW(y)
  #
  X.pv <- MeanFrames$X
  ###   OFFSET
  off <- model.offset(MeanFrames$mf) ## look for offset from (ori)Formula 
  if ( is.null(off) ) { ## no offset (ori)Formula term. Check attribute (Predictor ensures offset is not both in formula and attr)
    offsetObj <- attr(predictor,"offsetObj")
    off <- offsetObj$offsetArg 
  } else {
    offsetObj <- list()
    keepOriForm <- attr(predictor,"oriFormula")
    predictor <- noOffset(predictor)
    attr(predictor,"oriFormula") <- keepOriForm
  }
  if (!is.null(off)) {
    off <- pmax(log(.Machine$double.xmin),off) ## handling log(0) ## but if input off were NULL, output off would be is numeric(0) where it should remain NULL
  }
  colnames(X.pv) <- colnames(MeanFrames$X)
  if (ncol(X.pv)>0L) { 
    checkNAs <- coefficients(lm(processed$y ~ X.pv-1))   
    if (anyNA(checkNAs)) {   
      validbeta <- which(!is.na(checkNAs))
      X.pv <- structure(X.pv[,validbeta,drop=FALSE],namesOri=colnames(X.pv)) 
      # etaFix$beta |         variables 
      #             |  valid vars | invalid vars
      #     (1)           (2)           (3)
      # (2): colnames(<HLfit>$beta_cov) = colnames (<HLfit>$X.pv)
      # (1+2+3): namesOri, names(<HLfit>$fixef)
    } else attr(X.pv,"namesOri") <- colnames(X.pv)  
  } 
  # next extracts column for fixed oefficients
  X.Re <- X.pv ## may be updated from etaFix$beta or REMLformula, but this line makes the opposite possible
  ## reimplementation of etaFix$beta (2015/03)
  if ( length(betaFix <- etaFix$beta)>0 ) {
    namesbetafix <- names(betaFix)
    if (is.null(namesbetafix)) {
      message("The elements of etaFix$beta should be named and the names should match the column names of the design matrix.")
    }
    if (length(setdiff(namesbetafix,colnames(X.pv)))==0L) { ## if no incorrect name
      offFromEtaFix <- X.pv[,namesbetafix,drop=FALSE] %*% betaFix
      namesbetavar <- setdiff(colnames(X.pv),namesbetafix)
      X.pv <- structure(X.pv[,namesbetavar,drop=FALSE],namesOri=attr(X.pv,"namesOri"))
      offsetObj$betaFix <- betaFix
      if (is.null(off)) {
        off <- offFromEtaFix
      } else off <- off + offFromEtaFix
      ## TRUE by default:
      if ( is.null( keepInREML <- attr(betaFix,"keepInREML") ) ||  ( ! keepInREML) ) X.Re <- X.pv ## can be overwritten according to REMLformula
    } else {
      stop("The names of elements of etaFix$beta should all match column names of the design matrix.")
    }
  } # X.Re may again be modified later
  if (is.null(off)) { ## model.frame.default(formula = locform, offset = off,...) expects a vector....
    offsetObj$total <- rep(0,nobs)
    offsetObj$nonZeroInfo <- FALSE
  } else {
    offsetObj$total <- off
    offsetObj$nonZeroInfo <- TRUE
  }
  #
  attr(predictor,"offsetObj") <- offsetObj ## has $total, $nonZeroInfo, possibly $offsetArg (argument of Predictor()) and $betaFix. As nothing from the formula itself 
  processed$predictor <- predictor  
  processed$X.pv <- X.pv 
  pforpv <- ncol(X.pv)
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
    if (var(y)==0 && var(BinomialDen)==0) { warning("var(response) = 0.") }  
  } else {
    BinomialDen <- rep(1,nobs)
    if ( var(y)==0 ) { warning("var(response) = 0.") }  
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
  ZAlist <- get_ZAlist(predictor=predictor, model_frame=MeanFrames$mf)
  processed$ZAlist <- ZAlist ## attributes will be added later !
  nrand <- length(ZAlist)
  if (inherits(rand.families,"family")) rand.families <- list(rand.families) 
  if (nrand != 1L && length(rand.families)==1L) rand.families <- rep(rand.families,nrand) 
  rfblob <- checkRandLinkS(rand.families)  
  rand.families <- rfblob$rand.families
  lcrandfamfam <- rfblob$lcrandfamfam
  if (HLmethod=="ML") {
    HLmethod <- "ML(1,1,1)"  
  } else if (HLmethod=="SEM") {
    #if ( ! requireNamespace("probitgem",quietly = TRUE)) {## will pass CHECK if CRAN knows probitgem
    # if ( ! ("probitgem" %in% .packages()) ) { ## passed CRAN checks; appropriate sif probitgem explicitly loaded
    if ( ! ("probitgem" %in% .packages(all.available = TRUE)) ) { ## 
        stop("Package 'probitgem' not available for fitting by SEM.")
    } else eval(as.call(c(quote(require),list(package="probitgem", quietly = TRUE)))) ## passes CRAN checks (but temporary)
    #      eval(as.call(c(quote(requireNamespace),list(package="probitgem", quietly = TRUE)))) ## is not sufficient
    HLmethod <- "ML('SEM',NA,NA)" 
    if ( ! (family$family=="binomial" && family$link=="probit")) {
      stop("SEM is applicable only to binomial(probit) models.")
    }
    if (sum(BinomialDen) != nobs) {
      stop("(!) SEM procedure: the data do not seem binary; non-binary data are not handled.")
      
      # repsucces <- rep(seq(nrow(currentSample)),currentSample$succes)
      # repechec <- rep(seq(nrow(currentSample)),currentSample$echec)
      # currentSample <- currentSample[c(repsucces,repechec),]
      # currentSample$echec <- rep(c(0,1),c(length(repsucces),length(repechec)))
      # currentSample$succes <- 1L-currentSample$echec
      
    }
    SEMseed <- control.HLfit$SEMseed ##  
    # if (is.null(SEMseed)) SEMseed <- 123 ## OK pour SEM *unique* mais remplace par NULL dans optimthroughSmooth
    SEMargs <- list(SEMseed=SEMseed)
    SEMargs$nMCint <- control.HLfit$nMCint ##  as is => SEM procedure must handle null value
    SEMargs$get_SEM_info <- control.HLfit$get_SEM_info ##  as is => SEM procedure must handle null value
    SEMargs$control_pmvnorm$maxpts <- control.HLfit$pmvnorm_maxpts ##  as is => SEM procedure must handle null value
    #SEMlogL <- control.HLfit$SEMlogL
    #if (is.null(SEMlogL)) SEMlogL <- "pMVN" ## FR->FR 1.8.25 !
    #SEMargs$SEMlogL <- SEMlogL
    SEMargs$SEMlogL <- control.HLfit$SEMlogL
    nSEMiter <- control.HLfit$nSEMiter ##  
    if ( (! is.null(nSEMiter)) && nSEMiter < 10) {
      stop(" 'nSEMiter' should be >9")
    } else SEMargs$nSEMiter <- nSEMiter
    SEMargs$ngibbs <-control.HLfit$ngibbs ##  
    SEMargs$SEMsample <- control.HLfit$SEMsample ## stays NULL if NULL
    SEMargs$whichy1 <- (y==1) 
    processed$SEMargs <- SEMargs
  } else if (HLmethod=="REML") {
    HLmethod <- "HL(1,1,1)"  
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
  if (substr(HLmethod,0,2)=="ML") { # && HL[1]!="SEM") { ## FR->FR c'est bizarre d'exclure le SEM là... p_bv est il vraiment utilisé ?
    if ( ! is.null(REMLformula)) {
      message("Confusing combination of arguments: 'HLmethod=ML(...)' with non-null 'REMLformula'.")
      stop("  Make sure what you mean and simplify the arguments.")
    }
    lhs <- paste(predictor)[[2]] ## extract response, either cbind or not
    if ( nrand > 0L ) {
      ## build formula with only random effects  ##FR->FR  why have I done that ??
      REMLformula <- as.formula(paste(lhs,"~",paste(attr(ZAlist,"ranefs"),collapse="+")))
    } else REMLformula <- as.formula(paste(lhs,"~ 0")) 
    attr(REMLformula,"isML") <- TRUE
  } ## else do nothing: keeps input REMLformula, which may be NULL or a non-trivial formula
  # REMLformula <- .stripFormula(REMLformula)
  processed$REMLformula <- REMLformula  
  processed$loglfn.fix <- selectLoglfn(family)
  ## code derived from the glm() function in the safeBinaryRegression package
  if(family$family == "binomial" && length(unique(y)) == 2L && ncol(X.pv)>0L) {
    isSeparated <- .is_separated(X.pv, as.numeric(y))
    if(isSeparated) {
      warning(paste("Separation or quasi-separation exists among the sample points:",
                    "\n\tsome estimates of fixed-effect coefficients could be infinite,",
                    "\n\tcausing numerical issues in various functions."))
    }
    # separation <- separator(X.pv, as.numeric(y), purpose = "test")$separation
    # if(separation) {
    #   message("Separation exists among the sample points.\n\tThis model cannot be fit by maximum likelihood.")
    #   message("The following terms are causing separation among the sample points:")
    #   separation <- separator(X.pv, as.numeric(y), purpose = "find")$beta
    #   separating.terms <- dimnames(X.pv)[[2]][abs(separation) > 1e-09]
    #   if(length(separating.terms)) message(paste(separating.terms, collapse = ", "))
    #   stop()
    # }
  }
  #
  if ( ! is.null(REMLformula) ) { ## differences affects only REML estimation of dispersion params, ie which p_bv is computed
    REMLFrames <- HLframes(formula=REMLformula,data=data) ## design matrix X, Y...
    X.Re <- REMLFrames$X
    # wAugX will have lost its colnames...
    unrestricting_cols <- which(colnames(X.pv) %in% setdiff(colnames(X.pv),colnames(X.Re))) ## not in X.Re
    extra_vars <- setdiff(colnames(X.Re),colnames(X.pv)) ## example("update") tests this.
    distinct.X.ReML <- c(length(unrestricting_cols)>0L, length(extra_vars)>0L) ## TWO booleans
    if (any(distinct.X.ReML)){
      attr(X.Re,"distinct.X.ReML") <- distinct.X.ReML 
      if (attr(X.Re,"distinct.X.ReML")[1L]) attr(X.Re,"unrestricting_cols") <- unrestricting_cols
      if (attr(X.Re,"distinct.X.ReML")[2L])attr(X.Re,"extra_vars") <- extra_vars ## example("update") tests this.
      processed$X.Re <- X.Re
    } ## else processed$X.Re  is NULL
  } ## else keep previously computed X.Re = X.pv
  # if standard ML: there is an REMLformula ~ 0 (or with ranefs ?); local X.Re and processed$X.Re is 0-col matrix
  # if standard REML: REMLformula is NULL: $X.Re is X.pv, processed$X.Re is NULL
  # non standard REML: other REMLformula: $X.Re and processed$X.Re identical, and may take essentially any value
  if (is.null(objective)) {
    if (ncol(X.Re)>0) { ## standard or non-standard REML
      objective <- "p_bv"  ## info for fitme_body and corrHLfit_body, while HLfit instead may use return_only="p_bvAPHLs"
    } else objective <- "p_v"
  }
  processed$objective <- objective
  #
  models <- list(eta="",lambda="",phi="")
  if ( nrand > 0L ) {
    models[["eta"]] <- "etaHGLM" 
    vec_n_u_h <- unlist(lapply(ZAlist,ncol)) ## nb cols each design matrix = nb realizations each ranef
    processed$cum_n_u_h <- cum_n_u_h <- cumsum(c(0, vec_n_u_h)) ## if two ranef,  with q=(3,3), this is 0,3,6 ;
    nrd <- cum_n_u_h[nrand+1L]
    if (is.null(ZAfix <- attr(predictor,"ZALMatrix"))) { 
      ## then we could look further whether an LMatrix is available
      ##   but that does not look interesting, bc presumably LMatrix will be variable 
      ##   and in all case reconstruction will be done in HLfit_body
      processed$QRmethod <- .choose_QRmethod(ZAlist, predictor)
      # alternative (commented) code removed from version 1.11.35 
      ZAfix <- .post_process_ZALlist(ZAlist,as_matrix=.eval_as_mat_arg(processed))  
      # QRmethod is dense for adjacency and then ZAfix is dense here. 
      # ZAfix may be refined as sparse by post-processing in fitting function(s) if sparse_precision
      # but still it seems better to keep QRmethod dense in that case.
    } 
    if (inherits(ZAfix,"Matrix")) {
      AUGI0_ZX <- list(I=suppressWarnings(as(Diagonal(n=nrd),"CsparseMatrix")), ## avoids repeated calls to as() through rbind2...
                       ZeroBlock=Matrix(0,nrow=nrd,ncol=pforpv),
                       ZAfix=ZAfix, ## either ZA, or ZAL if the latter is fixed through attr(predictor,"ZALMatrix")
                       X.pv=X.pv)
    } else {
      AUGI0_ZX <- list(I=diag(nrow=nrd),ZeroBlock=matrix(0,nrow=nrd,ncol=pforpv),
                       ZAfix=ZAfix, ## either ZA, or ZAL if the latter is fixed through attr(predictor,"ZALMatrix")
                       X.pv=X.pv)
    }
    AUGI0_ZX$adjMatrix <- adjMatrix ## may be NULL
    AUGI0_ZX$envir <- list2env(list(), parent=environment(preprocess))
    processed$AUGI0_ZX <- AUGI0_ZX
    #
    # processed info for u_h inference
    ## builds box constraints either NULL or non-trivial, of length n_u_h
    lower.v_h <- rep(-Inf,cum_n_u_h[nrand+1L])
    boxConstraintsBool <- FALSE
    gaussian_u_ranges <- integer(0)  
    for (it in seq_len(nrand)) {
      if (lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="identity") { ## gamma(identity)
        u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
        lower.v_h[u.range] <- 1e-6 ## 1e-12 is disastrous
        boxConstraintsBool <- TRUE
      }
      if (lcrandfamfam[it]=="gaussian") gaussian_u_ranges <- c(gaussian_u_ranges,(cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])) 
    }
    if ( ! boxConstraintsBool ) lower.v_h <- NULL 
    boxConstraintsBool <- FALSE
    upper.v_h <- rep(Inf,cum_n_u_h[nrand+1L])
    for (it in seq_len(nrand)) {
      if (lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="-1/mu") { ## v ~ -Gamma
        u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
        upper.v_h[u.range] <- -1e-06
        boxConstraintsBool <- TRUE
      } else if (lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="log") {
        ## Gamma(log)$linkinv is pmax(exp(eta), .Machine$double.eps), ensuring that all gamma deviates are >= .Machine$double.eps
        ## we ensure that log(u_h) has symmetric bounds on log scale (redefine Gamma()$linkfun ?)
        upper.v_h <- Gamma(log)$linkfun(1/.Machine$double.eps)
      } 
    }
    if ( ! boxConstraintsBool ) upper.v_h <- NULL
    processed$u_h_info <- list(lower.v_h=lower.v_h,upper.v_h=upper.v_h) 
    
  } else {
    models[["eta"]] <- "etaGLM" 
  }
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
  } else if (family$family=="COMPoisson" && family$link=="loglambda") {
    canonicalLink <- TRUE
  } ## no implemented canonical link case for negbin
  processed$canonicalLink <- canonicalLink  
  #
  GLMMbool <- (nrand>0 && all(lcrandfamfam=="gaussian") ) ## only allowed gaussian rand.family is gaussian(identity) 
  if (GLMMbool) {
    LMMbool <- (family$family=="gaussian" && canonicalLink) 
  } else LMMbool <- FALSE
  processed$LMMbool <- LMMbool
  processed$GLMMbool <- GLMMbool
  processed$coef12needed <- ! ((family$family=="gaussian" && family$link=="identity")
                              || (family$family=="Gamma" && family$link=="log")  ) ## two ad hoc cases
  #
  if (is.null(subs_p_weights <- substitute(prior.weights))) {
    prior.weights <- structure(rep(1L,nobs),unique=TRUE) ## <- 1L prevented by glm -> model.frame(... prior.weights)
  } else if ( ! (inherits(subs_p_weights,"call") && subs_p_weights[[1L]] == "quote") )  {
    prior.weights <- as.vector(stats::model.weights(validdata)) ## as.vector as in say lm() protects against array1d
    if ( ! is.numeric(prior.weights)) 
      stop("'weights' must be a numeric vector")
    if (any(prior.weights < 0)) 
      stop("negative weights not allowed")
    attr(prior.weights,"unique") <- (length(unique(prior.weights))==1L) 
    #attr(prior.weights,"only1") <- all(upw==1L)
  } else attr(prior.weights,"unique") <- FALSE ## when 'prior.weights' is a quoted expression   
  processed$prior.weights <- prior.weights
  #
  ## algorithms (control of defaults remains in the HLfit code)
  betaFirst <- control.HLfit$betaFirst ##  
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
    } else if (LMMbool) {  
      LevenbergM <- FALSE ## because no reweighting when beta_eta changes => no IRWLS necess   
    } else if (HLmethod=="ML(0,0,1)") { ## ie, PQL/L in standardized syntax
      LevenbergM <- TRUE
    } else LevenbergM <- .spaMM.data$options$LevenbergM
  }
  processed$LevenbergM <- LevenbergM
  #
  if (nrand>0) {   
    lFix <- getPar(ranFix,"lambda")
    if (is.null(lFix)) lFix <- rep(NA,nrand)
    if (length(lFix)!=nrand) stop("length of 'fixed lambda' vector does not match number of random effects")
    names(lFix) <- as.character(seq(nrand))
    processed$lambda.Fix <- lFix
    nrand_lambda <- 0
    models[["lambda"]] <- rep("lamScal",nrand) ## even for adjacnency, random slope...
    ################################################################################
    # for a random slope term, ie v= v_1+x v_2 , the x went into the general ZAL matrix 
    # (construction of ZAlist by .spMMFactorList), and
    # we are still estimating the lambda's using a X_lamres with 0/1's only
    # unless there is a non-trivial model for the lambdas
    ################################################################################
    if (all(models[["lambda"]]=="lamScal")) { ## all mixed models handled in 06/2015 (random slope, adjacency...) hence currently always TRUE
      Xi_cols <- unlist(lapply(attr(ZAlist,"namesTerms"),length))
      if (any(Xi_cols>1 & !lcrandfamfam=="gaussian")) {
        stop("(!) random slope models with correlated non-gaussian random effects are not fitted.")
      }
      cum_Xi_cols <- cumsum(c(0, Xi_cols)) ## if two ranef,  with n_u_h=(3,3), this is 0,3,6. cum_h_u_h[nrand+1] is then 6, the total # of realizations
      n_u_h <- rep(0, sum(Xi_cols))
      #for (i in 1:nrand) n_u_h[(cum_Xi_cols[i]+1L):cum_Xi_cols[i+1L]] <- ncol(FL$Design[[i]]) ##  nlevels(FL$Subject[[i]])
      # if 18 responses in a random slope model ncol(FL$Design[[i]]) is 36 while nlevels(FL$Subject[[i]]) was 18
      for (i in 1:nrand) n_u_h[cum_Xi_cols[i]+seq(Xi_cols[i])] <- ncol(ZAlist[[i]])/Xi_cols[i]
      # h_u_h not n_u_h ...
      cum_h_u_h <- cumsum(c(0, n_u_h)) ## if two "Intercept" ranefs,  with n_u_h=(3,3), this is 0,3,6. cum_h_u_h[nrand+1] is then 6, the total # of realizations
        ## if (1+X|...) +(1|...),  with n_u_h=(3,4), this is 0,3,6,10. cum_h_u_h[sum(Xi_cols)+1] is then 10, the total # of realizations
      X_lamres <- matrix(0,cum_h_u_h[sum(Xi_cols)+1L],sum(Xi_cols))
      colnames(X_lamres) <- unlist(attr(ZAlist,"namesTerms"))
      for (i in seq(nrand)) {
        for (j in (cum_Xi_cols[i]+1L):cum_Xi_cols[i+1L]) {
          X_lamres[(cum_h_u_h[j]+1L):cum_h_u_h[j+1L],j] <- 1L ## this maps the deviance residuals to the lambda's to be estimated from them. None of the random-slope columns is a constant full column because each dev res is used for estimating only one lambda. Nevertheless, the first col will be called "(Intercept)", and this makes a valid output.
        }
      } 
    } else {  ## linear predictor for variance of random effects (lambda) (lamGLM or lamHGLM) 
      if (any(models[["lambda"]]=="lamHGLM")) { ##need distinct calls... to fit each lambda model  
        if (length(formulaLambda)==2) formulaLambda <- as.formula(paste('"lambda"',paste(formulaLambda,collapse=" ")))
        lambda_ranefs <- findbarsMM(formulaLambda)
        if (!is.null(lambda_ranefs)) {  ## lamHGLM
          nrand_lambda <- length(lambda_ranefs)
          models[["lambda"]] <- list("lamHGLM")
        } else models[["lambda"]] <- list("lamGLM")  
        colnames(X_lambda) <- colnames(fr_lambda$X) ## but code not effective, fr_lambda not computed
      } else { ## can use a single design matrix for all random effects, which is convenient.
        mess <- pastefrom("LIKELY missing code to handle linear predictor for lambda.")
        stop(mess)
        # la suite c'est dexu residus de code a assembler: il faut une liste de HLframes du type
        fr_lambda <- HLframes(formula=formulaLambda,data=data) ## but the "data" should probably be distinct data here, with nrow=number of reals of ranefs 
        # (pobablement calculee en amont pour determiner lamScal aussi...) ensuite extraire les design matrices
        #X_lamres ? Xi_cols ?
      }
    } 
    processed$X_lamres <- X_lamres ## for glm for lambda, and SEMbetalambda
    attr(processed$ZAlist,"Xi_cols") <- Xi_cols ## vector of ncol of X design for LHSs of (|) 
    attr(processed$ZAlist,"anyRandomSlope") <- any(Xi_cols>1L) ## used to handle random slope...
  } 
  #
  phi.Fix <- getPar(ranFix,"phi")
  if (is.null(phi.Fix)) {
    if (family$family %in% c("poisson","binomial","COMPoisson","negbin")) phi.Fix <- 1 
  } else if (any(phi.Fix==0)) stop("phi cannot be fixed to 0.")
  processed$phi.Fix <- phi.Fix
  #
  resid.predictor <- resid.model$formula
  if ( is.null(phi.Fix)) {
    phi_ranefs <- findbarsMM(resid.predictor)
    if ( ! is.null(phi_ranefs)) {
      if (is.null(resid.model$rand.family)) resid.model$rand.family <- gaussian() # avoids rand.families being NULL in call below.
      resid.predictor <- as.formula(paste(".phi",.DEPARSE(resid.predictor)))
      data$.phi <- seq(nobs) ##  I need a dummy variabel response else var(y)==0 test will catch it.
      preprocess_arglist <- list(control.HLfit=control.HLfit, ## constrained
                         ranFix=resid.model$fixed, ## not constrained, but should rather use 'resid.model$fixed' (but anyway preprocess checks only  phi and lambda)
                         HLmethod=HLmethod, ## constrained
                         predictor=resid.predictor, ## obvious
                         resid.model=~1, # no resid.resid.model... 
                         REMLformula=NULL, # constrained
                         data=data, # obvious (?) 
                         family=resid.family, # obvious
                         BinomialDen=NULL, # obviously no binomial response
                         rand.families=resid.model$rand.family, # (NULL not handled by preprocess); 
                         #   outer preprocess calls *receive* adefautl value from formals(HLfit)
                         etaFix=resid.model$etaFix, ## not constrained, but should rather use 'resid.model$fixed'
                         prior.weights=NULL, ## currently defined  dynamically using lev_phi...
                         control.glm=control.glm ## constrained
      )
      residProcessed <- do.call(preprocess,preprocess_arglist)
      # preprocess here plays the role of fitme as wrapper bringing the following info to fitme_body:
      processed$residProcessed <- residProcessed
      models[["phi"]] <- "phiHGLM" 
      p_phi <- NA
    } else {
      if (! inherits(resid.predictor,"predictor")) resid.predictor <- Predictor(resid.predictor)
      residFrames <- HLframes(formula=resid.predictor,data=data)
      dispOffset <- model.offset(residFrames$mf) ## look for offset from (ori)Formula 
      if ( is.null(dispOffset) ) { ## ## no offset (ori)Formula term. Check attribute (Predictor ensures offset is not both in formula and attr)
        dispOffset <- attr(resid.predictor,"offsetObj")$total 
      } else {
        keepOriForm <- attr(resid.predictor,"oriFormula")
        resid.predictor <- noOffset(resid.predictor)
        attr(resid.predictor,"oriFormula") <- keepOriForm
      }
      if (is.null(dispOffset)) { ## model.frame.default(formula = locform, offset = off,...) expects a vector....
        attr(resid.predictor,"offsetObj") <- list(total=rep(0,nobs),nonZeroInfo=FALSE)
      } else {
        attr(resid.predictor,"offsetObj") <- list(total=dispOffset,nonZeroInfo=TRUE) ## subtly, will be of length 1 if   original offset was a constant...
      }
      ## if formula= ~1 and data is an environment, there is no info about nobs, => fr_disp$X has zero rows, which is a problem later 
      p_phi <- NCOL(residFrames$X)
      namesX_disp <- colnames(residFrames$X)
      if (p_phi==1 && namesX_disp[1]=="(Intercept)"
          && is.null(dispOffset) ## added 06/2016 (bc phiScal does not handle offset in a phi formula) 
      ) {
        models[["phi"]] <- "phiScal"
      } else { 
        models[["phi"]] <- "phiGLM"
      }
      resid.model$predictor <- resid.predictor  ## absent  if phiHGLM has been detected
    } 
  } else p_phi <- 0
  processed$p_phi <- p_phi # no X_disp is saved in processed
  if ( family$family %in% c("binomial","poisson","COMPoisson","negbin")) {
    ## the response variable should always be Counts
    if (max(abs(y-as.integer(y)))>1e-05) {
      mess <- pastefrom("response variable should be integral values.",prefix="(!) From ")
      stop(mess)
    }
    if ( .DEPARSE(resid.predictor) != "~1") {
      warning(paste("resid.model is ignored in ",family$family,"-response models",sep=""))
    }
  }
  #
  if (LMMbool) {
    ## identifiability checks cf modular.R -> checkNlevels() in lmer:
    vec_n_u_h <- diff(processed$cum_n_u_h)
    if (any(vec_n_u_h<2L) && is.null(phi.Fix)) {
      problems <- which(vec_n_u_h<2L) 
      for (iMat in problems) {
        mess <- paste("Only ",vec_n_u_h[iMat]," level for random effect ",
                      attr(ZAlist,"ranefs")[iMat],
                      ";\n   this model cannot be fitted unless phi is fixed.",sep="")
        warning(mess)
      }
    }
    if (any(vec_n_u_h==nobs) && models[["phi"]] %in% c("phiScal","phiGLM")) { ## tests (mean)ranefs if intecept in resid formula
      resid.mf <- residFrames$mf 
      if (attr(attr(resid.mf, "terms"),"intercept")!=0L) { ## there is an intercept in the resid.formula
        # ideally for random-coefficients models we should compare the design columns... 
        ## FR->FR cf isNested check as in https://github.com/lme4/lme4/blob/master/R/utilities.R, 
        problems <- which(vec_n_u_h==nobs) 
        for (iMat in problems) {
          term_ranef <- attr(ZAlist,"ranefs")[iMat]
          if (# is.null(LMatrix) && ## does not seem useful
            substr(term_ranef, 1, 1)=="(" ## excludes spatial ranefs 
          ) {
            mess <- paste("Number of levels = number of observations",
                          "\n   for random effect ", term_ranef,
                          ";\n   this model cannot be fitted unless phi is fixed",
                          "\n   or a correlation matrix is given.",sep="")
            stop(mess)
          }          
        }
      }
    }
  } 
  #
  processed$models <- models
  attr(rand.families,"lcrandfamfam") <- lcrandfamfam
  attr(rand.families,"unique.psi_M") <- sapply(lcrandfamfam, function(v) {
    switch(v, 
           gaussian = 0,
           gamma = 1, 
           beta = 1/2, 
           "inverse.gamma" = 1
    )
  })
  processed$residModel <- resid.model 
  processed$rand.families <- rand.families
  processed$fixef_terms <- MeanFrames$fixef_terms ## added 2015/12/09 for predict
  processed$fixef_levels <- MeanFrames$fixef_levels ## added 2015/12/09 for predict
  processed[["control.glm"]] <- do.call("glm.control", control.glm) ## added 04/2016 (LHS <- RHS list)
  processed$port_env <- new.env(parent=emptyenv()) ## added 09/2016
  class(processed) <- c("arglist","list")
  return(processed)
}

.eval.update.call <- function(mc,...) { # not currently used
  mc <- as.list(mc)
  dotlist <- list(...)
  mc[names(dotlist)] <- dotlist ## a un moment j'ai mis cette ligne en commentaire, ce qui rend la fonction ineffective !
  eval(as.call(mc))  
}