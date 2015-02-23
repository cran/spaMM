`extract.check.coords` <- function(spatial.model,datanames) {
  if ( ! is.null(spatial.model)) {
    bars <- spatial.model[[2]] 
    coordinates <- DEPARSE(bars[[3]]) ## "x + y"
    coordinates <-  strsplit(coordinates," ")[[1]]
    coordinates <- setdiff(coordinates,c("+","%in%",":","/"))
  } else {
    stop("Call to 'HLCor' without a spatial term in the formula is suspect.")
    ## very old code handling old syntax with (1|pos) and default values of the coordinates argument
    coordinates <- c("x","y") ## back compat
  }
  coordcols <- which(datanames %in% coordinates)
  if ( length(coordcols) != length(coordinates) ) {
    stop("variables 'coordinates' not all found in the 'data'")
  }
  return(coordinates) ## should be ordered as bars[[3]] (important for predict)
}

## better for development to avoir name conflicts with OKsmooth :toCanonical and :canonize
canonizeRanPars <- function(ranPars,corr.model,checkComplete=TRUE) {
  trueCorrpars <- list()
  if (corr.model %in% c("Matern","corMatern")) {
    if (!is.null(ranPars$trNu)) { ## either we have nu,rho or trNu,trRho 
      ranPars$nu <- nuInv(ranPars$trNu,ranPars$trRho) ## before trRho is removed...
      ranPars$trNu <- NULL
      attr(ranPars,"type")$nu <- attr(ranPars,"type")$trNu
      attr(ranPars,"type")$trNu <- NULL
    } 
    nu <- ranPars$nu
    if (is.null(nu) && checkComplete) {
      mess <- pastefrom("nu missing from ranPars (or correlation model mis-identified).",prefix="(!) From ")
      stop(mess)
    }
    trueCorrpars$nu <- nu 
  } 
  if (corr.model=="AR1") {
    ARphi <- ranPars$ARphi
    if (is.null(ARphi) && checkComplete) {
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
    if (is.null(rho) && checkComplete) {
      mess <- pastefrom("rho missing from ranPars.",prefix="(!) From ")
      stop(mess)
    }
    trueCorrpars$rho <- rho
  }
  Nugget <- ranPars$Nugget
  if (! is.null(Nugget)) trueCorrpars$Nugget <- Nugget 
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
  return(list(trueCorrpars=trueCorrpars,ranPars=ranPars))
}

HLCor <- function(formula,
                  ranPars, ## all dispersion and correlation params ideally provided through ranPars
                  data,
                  distMatrix,uniqueGeo=NULL,adjMatrix,corrMatrix,
                  verbose=c(warn=TRUE,trace=FALSE,summary=FALSE),control.dist=list(),
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
    rpblob <- canonizeRanPars(ranPars=ranPars,corr.model=corr.model) 
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
      decomp <- attr(adjMatrix,"symSVD")
      if (is.null(decomp)) {
        m <- solve(diag(rep(1,nrow(adjMatrix)))-rho*(adjMatrix))
      } else {
        m <- ZWZt(decomp$u,1/(1-rho*decomp$d))
      }
    }  else if (corr.model=="AR1") {
      coordinates <- extract.check.coords(spatial.model=spatial.model,datanames=names(data))
      uniqueGeo <- unique(data[,coordinates,drop=F]) ## keeps the names of first instances of the coordinates in data
      txt <- paste(spatial.model[[2]][[3]]) ## the RHS of the ( . | . ) 
      if (length(grep("%in%",txt))>0) {
        stop("HLCor code should be allowed again to handle blockDiag objects")
        #scaled.dist <- as.blockDiag.bar(spatial.model[[2]],formula,data=uniqueGeo)
        #test.in <- TRUE
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
      if ( ! is.null(dist.method <- control.dist$`dist.method`)) {
        msd.arglist$`dist.method` <- dist.method
      }
      
      if (length(rho)>1L) {
        msd.arglist <- c(msd.arglist,list(uniqueGeo=uniqueGeo))
        if ( ! is.null(rho.mapping <- control.dist$`rho.mapping`)) {
          msd.arglist$`rho.mapping` <- rho.mapping
        }
      } else {
        if ( missing(distMatrix)) { 
          dist.arglist <- list(x=uniqueGeo)
          if(!is.null(dist.method <- control.dist$dist.method)) dist.arglist$method <- dist.method
          distMatrix <- do.call(proxy::dist,dist.arglist)
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
  fixNames <- setdiff(names(ranPars),varNames) 
  if (!is.null(fixNames)) { ## could be NULL for corrMatrix case
    ranFix <- ranPars[fixNames] ## 11/2014 as there is no other source for ranFix
    typelist <- list() 
    typelist[fixNames] <- "fix" 
    if (!is.null(rPtype <- attr(ranPars,"type"))) { ## it may not exist, or elements may be "fix" or "outer"
      typelist[names(rPtype)] <- rPtype
    }
    attr(ranFix,"type") <- typelist 
    HL.info$ranFix <- ranFix
  }
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
  hlfit$control.dist <- control.dist
  attr(hlfit,"info.uniqueGeo") <- uniqueGeo ## more spatial info is to be found in hlfit$predictor (Lunique = corrmat^1/2) and hlfit$ZALMatrix
  if (corr.model %in% c("Matern","corMatern")) attr(hlfit,"msd.arglist") <- msd.arglist ## more organized, easier to reuse. 
  ## FR->FR but info.uniqueGeo more general (eg AR1) -> a revoir
  attr(hlfit,"HLCorcall") <- mc
  if (verbose["HLCorSummary"]) { ## useful in final call from corrHLfit
    summary(hlfit) ## input corr pars have been printed at the beginning...   
  }
  return(hlfit) ## 
}


## wrapper for HLCor, suitable input and output for optimization
`HLCor.obj` <- function(ranefParsVec,skeleton,HLCor.obj.value="p_bv",trace=NULL,family=gaussian(),...) { ## name of first arg MUST differ from names in dotlist...
  mc <- match.call(expand.dots=TRUE) ## (1) expand.dots added 11/04/2014 for the mutlinomial... eval 
  dotlist <- list(...)
  if ( inherits(mc$data,"list")) {
    ## then processed should already be a list
    family <- mc$family
    fitlist <- lapply(seq_len(length(mc$data)),function(it){
      locmc <- mc
      locmc[[1L]] <- as.name("HLCor.obj") ## replaces "f" !
      locmc$ranefParsVec <- ranefParsVec ## replaces "arg" !
      if (family$family=="multi") locmc$family <- family$binfamily
      locmc$distMatrix <- mc$distMatrix[[it]]
      locmc$uniqueGeo <- mc$uniqueGeo[[it]]
      locmc$data <- mc$data[[it]]
      locmc$processed <- mc$processed[[it]]
      eval(locmc) ## this will execute all the code below starting from dotlist <- list(...) 
    })
    resu <- sum(unlist(fitlist))
    if (is.character(trace)) {
      verif <- paste("#global:",ranefParsVec,resu) 
      write(verif,file=trace,append=T) ## the file is unlink'ed in corrHLfit()  
    }
    return(resu)
  }
  HLCor.formals <- names(formals(HLCor))
  HLfit.formals <- names(formals(HLfit))
  designL.formals <- names(formals(designL.from.Corr))
  makescaled.formals <- names(formals(make.scaled.dist))
  HLnames <- (c(HLCor.formals,HLfit.formals,designL.formals,makescaled.formals))  ## cf parallel code in corrHLfit
  HLCor.args <- dotlist[intersect(names(dotlist),HLnames)]
  forGiven <- relist(ranefParsVec,skeleton) ## given values of the optimized variables
  HLCor.args$ranPars[names(forGiven)] <- forGiven ## do not wipe out other fixed, non optimized variables
  HLCor.args$family <- family 
  if (is.character(trace)) {
    if(.spaMM.data$options$TRACE.UNLINK) unlink("HLCor.args.*.RData")
    zut <- paste(ranefParsVec,collapse="")  
    save(HLCor.args,file=paste("HLCor.args.",zut,".RData",sep="")) ## for replicating the problem
  }
  hlfit <- do.call("HLCor",HLCor.args)
  aphls <- hlfit$APHLs
  resu <- aphls[[HLCor.obj.value]]
  readable <- unlist(canonizeRanPars(ranPars=forGiven,corr.model=dotlist$`corr.model`,checkComplete=FALSE)$ranPars) ## FR->FR use of dotlist...
  verif <- c(unlist(aphls),hlfit$lambda,hlfit$phi,readable,ranefParsVec) ## hlfit$phi may be NULL
  if (is.character(trace)) {
    write(verif,file=trace,ncolumns=length(verif),append=T) ## the file is unlink'ed in corrHLfit()  
  }
  return(resu) #
}


