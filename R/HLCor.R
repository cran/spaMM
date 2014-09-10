HLCor <-
function(formula,
                  ranPars, ## all dispersion and correlation params ideally provided through ranPars
                  data,
                  distMatrix,uniqueGeo=NULL,adjMatrix,corrMatrix,
                  verbose=c(warn=TRUE,trace=FALSE,summary=FALSE),
                  ...) { 
  mc <- match.call() ## potentially used by getCallHL(object) in update.HL...
  if (is.na(verbose["trace"])) verbose["trace"] <- FALSE
  if (is.na(verbose["warn"])) verbose["warn"] <- TRUE
  if (is.na(verbose["summary"])) {
    verbose["HLCorSummary"] <- FALSE
  } else verbose["HLCorSummary"] <- verbose["summary"]
  verbose["summary"] <- FALSE ## this is for HLfit
  ## either these two lines, or a family argument and   HL.info$family <- family
  dotlist <- list(...)
  family <- dotlist$family
  family <- checkRespFam(family)
  dotlist$family <- family
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
  } else { ## preprocess has not been called hence validData probably not called before preprocess...
    validdata <- validData(formula=formula,resid.formula=dotlist$resid.formula,data=data) ## will remove rows with NA's in required variables
    if (!inherits(data,"environment")) {
      data <- data[rownames(validdata),,drop=FALSE] ##     before Predictor is called and an LMatrix is added, etc. 
    } else data <- validdata
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
  ## convert back ranPars to canonical scale:
  if (corr.model== "corrMatrix") {
    ranPars <- NULL
  } else {
    rpblob <- toCanonical(ranPars=ranPars,corr.model=corr.model) 
    ranPars <- rpblob$ranPars
    trueCorrpars <- rpblob$trueCorrpars
    rho <- ranPars$rho
  }
  coordinates <- NULL
  ##
  test.in <- FALSE
    if (corr.model %in% c("BreslowC93","adjacency")) { 
      ## no nugget in the adjacency model...
      if ( missing(adjMatrix) ) stop("missing 'adjMatrix' for adjacency model")
      m <- solve(diag(rep(1,nrow(adjMatrix)))-rho*(adjMatrix))
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
        norho <- trueCorrpars; norho$rho <- NULL ## because the Matern.corr input will be an already scaled distance 'm'
        m <- do.call(Matern.corr,args=c(norho,list(d=m)))        
      } 
    } else if (corr.model== "corrMatrix") {
      if (missing(corrMatrix)) {
        mess <- pastefrom("missing(corrMatrix) argument despite corrMatrix term in formula.",prefix="(!) From ")
        stop(mess)
      } ## ELSE:
      m <- corrMatrix
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
  ## convert ranPars to ranFix + init.HLfit
  ## allows log and not log:
  varNames <- names(which(attr(ranPars,"type")=="var"))
  HL.info$init.HLfit[varNames] <- ranPars[varNames]
  fixNames <- setdiff(names(ranPars),varNames) ## names(which(attr(ranPars,"type")=="fix")) pas correct si pas d'attribute (appel direct HLCor) -> on suppose que c'est fix
  HL.info$ranFix[fixNames] <- ranPars[fixNames]
  hlfit <- do.call("HLfit",HL.info) ## with a _list_ of arguments -> do.call ## perhaps should create a list of unevaluated arguments ???? 
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
  attr(hlfit,"HLCorcall") <- mc
  if (verbose["HLCorSummary"]) { ## useful in final call from corrHLfit
    summary(hlfit) ## input corr pars have been printed at the beginning...   
  }
  return(hlfit) ## 
}
