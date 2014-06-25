HLCor <-
function(formula,ranPars, ## all dispersion and correlation params ideally provided through ranPars
                  data,
                  distMatrix,uniqueGeo=NULL,adjMatrix,corrMatrix,
                  verbose=c(warn=TRUE,trace=FALSE,summary=FALSE),
                  ...) { 
  mc <- match.call() ## potentially used by getCall(object) in update.HL...
  if (is.na(verbose["trace"])) verbose["trace"] <- FALSE
  if (is.na(verbose["warn"])) verbose["warn"] <- TRUE
  if (is.na(verbose["summary"])) {
    verbose["HLCorSummary"] <- FALSE
  } else verbose["HLCorSummary"] <- verbose["summary"]
  verbose["summary"] <- FALSE ## this is for HLfit
  ## either these two lines, or a family argument and   HL.info$family <- family
  dotlist <- list(...)
  family <- dotlist$family
  famfam <- family$family
  if ( ! is.null(famfam) && famfam=="multi") {
    if ( ! inherits(data,"list")) {
      if(family$binfamily$family=="binomial") {
        familyargs <- family
        familyargs$family <- NULL
        familyargs$binfamily <- NULL
        data <- do.call(binomialize,c(list(data=data),familyargs)) ## if data not already binomialized
      }
    }
  }
  ################# data LIST
  if ( inherits(data,"list")) {
    processed <- dotlist$processed
    if ( ( ! is.null(processed) ) && ( ! inherits(processed,"list") ) ) {
      stop("(!) 'processed' is not NULL but not a list, while data is a list.")
    }
    fitlist <- lapply(seq_len(length(data)),function(it){
      locmc <- mc
      if (family$family=="multi") locmc$family <- family$binfamily
      locmc$data <- data[[it]]
      locmc$distMatrix <- mc$distMatrix[[it]]
      locmc$uniqueGeo <- mc$uniqueGeo[[it]]
      if (inherits(processed,"list")) locmc$processed <- processed[[it]]
      eval(locmc)
    }) ## a pure list of HLCor objects
    liks <- sapply(fitlist,function(v) {unlist(v$APHLs)})
    liks<- apply(liks,1,sum)
    attr(fitlist,"APHLs") <- as.list(liks)
    class(fitlist) <- c("HLfitlist",class(fitlist)) 
    return(fitlist) ## list of HLfit object + one attribute
  }
  ################# 
  if ( ! is.null(dotlist$ranFix)) { ##debugging code will become obsolete some day
    stop("!From HLCor: ranFix found in '...'. Make sure to use ranPars only")
  }
  if (!is.null(HLmethod <- dotlist$HLmethod)) {
    if (HLmethod=="SEM") dotlist$`try.chol` <- FALSE
  }
  processed <- dotlist$processed ## FR->FR suggests we should add it as argument of HLCor...
  if ( ! is.null(processed)) {
    predictor <- processed$predictor 
  } else { ## preprocess has not been called hence validRows probably not called before preprocess...
    validrows <- validRows(formula=formula,resid.formula=dotlist$resid.formula,data=data) ## will remove rows with NA's in required variables
    data <- data[validrows,,drop=FALSE] ## ## before Predictor is called and an LMatrix is added, etc.
    predictor <- formula
    if (! "predictor" %in% class(predictor)) predictor <- Predictor(formula) 
  }
  spatial.terms <- findSpatial(predictor)
  spatial.model <- spatial.terms[[1]] 
  if ( ! is.null(spatial.model)) {
    corr.model <- as.character(spatial.model[[1]]) 
  } else {
    if ( ! missing(corrMatrix)) {
      mess <- pastefrom("corrMatrix argument despite no corrMatrix term in formula:",prefix="(!) From ")
      message(mess)
      stop("This syntax is obsolete; add a corrMatrix(...) term in the formula.")
    } ## ELSE more generic message: 
    stop("Call to 'HLCor' without a spatial term in the formula is suspect.")
  }
  trueCorrpars <- list()
  if (corr.model %in% c("Matern","corMatern")) {
    if (!is.null(ranPars$trNu)) { ## either we have nu,rho or trNu,trRho 
        ranPars$nu <- nuInv(ranPars$trNu,ranPars$trRho) ## before trRho is removed...
        ranPars$trNu <- NULL
        attr(ranPars,"type")$nu <- attr(ranPars,"type")$trNu
        attr(ranPars,"type")$trNu <- NULL
    } 
    nu <- ranPars$nu
    if (is.null(nu)) {
      mess <- pastefrom("nu missing from ranPars (or correlation model mis-identified).",prefix="(!) From ")
      stop(mess)
    }
    trueCorrpars$nu <- nu 
  } 
  if (corr.model=="AR1") {
    ARphi <- ranPars$ARphi
    if (is.null(ARphi)) {
      mess <- pastefrom("ARphi missing from ranPars.",prefix="(!) From ")
      stop(mess)
    }
    trueCorrpars$ARphi <- ARphi    
  } else if (corr.model != "corrMatrix") { ## all models with a 'rho' parameter
    if (!is.null(ranPars$trRho)) {
      ranPars$rho <- rhoInv(ranPars$trRho)
      ranPars$trRho <- NULL
      attr(ranPars,"type")$rho <- attr(ranPars,"type")$trRho
      attr(ranPars,"type")$trRho <- NULL
    } ## else there may simply be rho rather than trRho (including for adjacency model through optim procedure !)
    rho <- ranPars$rho
    if (is.null(rho)) {
      mess <- pastefrom("rho missing from ranPars.",prefix="(!) From ")
      stop(mess)
    }
    trueCorrpars$rho <- rho
  }
  coordinates <- NULL
  ##
  test.in <- FALSE
    if (corr.model %in% c("BreslowC93","adjacency")) { 
      ## no nugget in the adjacency model...
      if ( missing(adjMatrix) ) stop("missing 'adjMatrix' for adjacency model")
      m<-solve(diag(rep(1,nrow(adjMatrix)))-rho*(adjMatrix))
    }  else if (corr.model=="AR1") {
      coordinates <- extract.check.coords(spatial.model=spatial.model,datanames=names(data))
      uniqueGeo <- unique(data[,coordinates,drop=F]) ## keeps the names of first instances of the coordinates in data
      #browser()
      txt <- paste(spatial.model[[2]][[3]]) ## the RHS of the ( . | . ) 
      if (length(grep("%in%",txt))>0) {
        scaled.dist <- as.blockDiag.bar(spatial.model[[2]],formula,data=uniqueGeo)
        test.in <- TRUE
      } else scaled.dist <- proxy::dist(uniqueGeo)
      m <- trueCorrpars$ARphi^scaled.dist
    } else  if (corr.model %in% c("Matern","corMatern")) {
      txt <- paste(spatial.model[[2]][[3]]) ## the RHS of the ( . | . ) 
      if (length(grep("%in%",txt))>0) {
        stop("(!) Matern( . | <coord> %in% <grp>) is not yet handled.")
        test.in <- TRUE ## should be useful when this case will be handled
      } 
      ## in a typical call from corrHLfit the following test should be FALSE because uniqueGeo and maybe distMatrix should have been precomputed
      if ((length(rho)>1 || missing(distMatrix)) && is.null(uniqueGeo)) { ## all cases where we need uniqueGeo
        coordinates <- extract.check.coords(spatial.model=spatial.model,datanames=names(data))
        uniqueGeo <- unique(data[,coordinates,drop=F]) ## keeps the names of first instances of the coordinates in data
      } 
      ## then compute scaled distances from unscaled info, for HLfit call
      msd.arglist <- list(rho = rho)
      if (length(rho)>1) {
        msd.arglist <- c(msd.arglist,list(uniqueGeo=uniqueGeo))
        if ( ! is.null(dotlist$`rho.mapping`)) msd.arglist <- c(msd.arglist,list(`rho.mapping`=dotlist$`rho.mapping`))
      } else {
        if ( missing(distMatrix)) { 
          distMatrix <- proxy::dist(uniqueGeo)
        }
        msd.arglist <- c(msd.arglist,list(distMatrix=distMatrix))
      }
      m <- do.call("make.scaled.dist",msd.arglist)
      ## at this point is a single location, m should be dist(0) and make.scaled.dist was modified to that effect
      if ( nrow(m)>1 ) { ## >1 locations
        Nugget <- ranPars$Nugget
        if (! is.null(Nugget)) trueCorrpars$Nugget <- Nugget ## else trueCorrpars keeps its default nu value
        nunu <- trueCorrpars; nunu$rho <- NULL ## because the Matern.corr input will be an already scaled distance 'm'
        m <- do.call(Matern.corr,args=c(nunu,list(d=m)))        
      } 
    } else if (corr.model== "corrMatrix") {
      if (missing(corrMatrix)) {
        mess <- pastefrom("missing(corrMatrix) argument despite corrMatrix term in formula.",prefix="(!) From ")
        stop(mess)
      } ## ELSE:
      m <- corrMatrix
      ranPars <- NULL
    } 
    if (verbose["trace"] && length(trueCorrpars)>0) { 
      print(unlist(trueCorrpars))
    }
    # print(paste("Correlation params in HLCor",paste(c(rho,nu),collapse=" ")))
    argsfordesignL <- dotlist[intersect(names(dotlist),names(formals(designL.from.Corr)))] 
    if ("dist" %in% class(m)) {
      m <- as.matrix(m)
      diag(m) <- 1L ## always a correlation matrix
    }
    if (is.null(attr(predictor,"LMatrix"))) { ## test FR 11/2013
      Lunique <- try(do.call(designL.from.Corr,c(list(m=m),argsfordesignL)))
      if (class(Lunique)=="try-error") { 
        print("correlation parameters were:") ## makes sense if designL.from.Corr already issued some warning
        print(unlist(trueCorrpars))    
        stop()
      }
      attr(Lunique,"ranefs") <- unlist(lapply(spatial.terms,DEPARSE))
      attr(predictor,"LMatrix") <- Lunique
      attr(predictor,"%in%") <- test.in
    }
###
  dotlist$verbose <- verbose[intersect(names(verbose),c("warn","trace","summary"))] ## all printing in HLfit is suppressed by default
  HLFormals <- names(formals(HLfit))
  HL.info <- dotlist[intersect(names(dotlist),HLFormals)]
  HL.info$data <- data
  if (! is.null(processed)) {
    processed$predictor <- predictor
    HL.info$processed <- processed
  } else HL.info$formula <- predictor
## this should become obsolete. see handling of ranPars below:
  if (!is.null(dotlist$LamFix)) {
    stop("argument LamFix of HLCor is obsolete")
  }
  if (!is.null(dotlist$PhiFix)) {
    stop("argument PhiFix of HLCor is obsolete")
  }  
##
  if (!is.null(ranPars$trPhi)) {
    ranPars$phi <- dispInv(ranPars$trPhi)
    ranPars$trPhi <- NULL
    attr(ranPars,"type")$phi <- attr(ranPars,"type")$trPhi
    attr(ranPars,"type")$trPhi <- NULL
  } else if (!is.null(ranPars$logphi)) { ## debug code
    ## HL.info$ranFix$phi <- exp(ranPars$logphi)
    stop("logphi in HLCor...")
  } #####################  else HL.info$ranFix$phi <- ranPars$phi ## y st deja !?
  if (!is.null(ranPars$trLambda)) {## 
    ranPars$lambda <- dispInv(ranPars$trLambda)
    ranPars$trLambda <- NULL
    attr(ranPars,"type")$lambda <- attr(ranPars,"type")$trLambda
    attr(ranPars,"type")$trLambda <- NULL
  } else if (!is.null(ranPars$loglambda)) { ## debug code
    stop("loglambda in HLCor...")
  } ##################### else HL.info$ranFix$lambda <- ranPars$lambda
  ## convert ranPars to ranFix + init.HLfit
  ## allows log and not log:
  varNames <- names(which(attr(ranPars,"type")=="var"))
  HL.info$init.HLfit[varNames] <- ranPars[varNames]
  fixNames <- setdiff(names(ranPars),varNames) ## names(which(attr(ranPars,"type")=="fix")) pas correct si pas d'attribute (appel direct HLCor) -> on suppose que c'est fix
  HL.info$ranFix[fixNames] <- ranPars[fixNames]
  hlfit <- do.call(HLfit,HL.info) ## with a _list_ of arguments -> do.call ## perhaps should create a list of unevaluated arguments ???? 
  if ( ! is.null(hlfit$error)) {
    errfile <- generateFileName("HLfitCall")
    errfile <- paste(errfile,".RData",sep="")
    save(HL.info,file=errfile)
    mess <- pastefrom("'do.call(HLfit,HL.info)' failed:",prefix="(!) From ")
    message(mess)
    message(hlfit$error)
    message("'HL.info' is saved in the ",errfile," file",sep="")
    stop("I exit.")
  } ## ELSE:
  attr(hlfit,"info.uniqueGeo") <- uniqueGeo ## more spatial info is to be found in hlfit$predictor (Lunique = corrmat^1/2) and hlfit$ZALMatrix
  # attr(hlfit,"call") <- mc
  if (verbose["HLCorSummary"]) { ## useful in final call from corrHLfit
    summary(hlfit) ## input corr pars have been printed at the beginning...   
  }
  return(hlfit) ## 
}
