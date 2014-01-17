HLCor <-
function(formula,ranPars, ## all dispersion and correlation params ideally provided through ranPars
                  data,
                  distMatrix,uniqueGeo,adjMatrix,corrMatrix,
                  verbose=rep(FALSE,3),...) { 
  mc <- match.call()
  lv <- length(verbose)
  if (lv<3) verbose<- c(verbose,c(F,F,F)[(lv+1):3]) ## verbose[1] -> printing from HLCor, [2:3] -> argument for HLfit 
  dotlist <- list(...)
  if ( ! is.null(dotlist$ranFix)) { ##debugging code will become obsolete some day
    stop("!From HLCor: ranFix found in dotlist. Make sure to use ranPars only")
  }
  processed <- dotlist$processed ## FR->FR suggests we should add it as argument of HLCor...
  if (! is.null(processed)) {
    predictor <- processed$predictor 
  } else {
    predictor <- formula
    if (! "predictor" %in% class(predictor)) predictor <- Predictor(formula) 
  }
  spatial.model <- findSpatial(predictor)[[1]] 
  if ( ! is.null(spatial.model)) {
    corr.model <- as.character(spatial.model[[1]]) 
  } else {
    stop("Call to 'HLCor' without a spatial term in the formula is suspect.")
    ## very early versions instead assumed:
    corr.model <- "Matern" 
  }
  trueCorrpars <- list()
  if (corr.model %in% c("Matern","corMatern")) {
    if (!is.null(ranPars$trNu)) { ## either we have nu,rho or trNu,trRho 
        ranPars$nu <- nuInv(ranPars$trNu,ranPars$trRho) ## before trRho is removed...
        ranPars$trNu <- NULL
    } 
    nu <- ranPars$nu
    if (is.null(nu)) {
      mess <- pastefrom("nu missing from ranPars (or correlation model misidentifie).",prefix="(!) From ")
      stop(mess)
    }
    trueCorrpars$nu <- nu 
  } 
  if (!is.null(ranPars$trRho)) {
    ranPars$rho <- rhoInv(ranPars$trRho)
    ranPars$trRho <- NULL
  } ## else there may simply be rho rather than trRho (including for adjacency model through optim procedure !)
  rho <- ranPars$rho
  if (is.null(rho)) {
    mess <- pastefrom("rho missing from ranPars.",prefix="(!) From ")
    stop(mess)
  }
  trueCorrpars$rho <- rho
  coordinates <- NULL
  if (missing(corrMatrix)) {
    if (corr.model %in% c("Matern","corMatern")) {
    ## no nugget in the adjacency model cf a few lines below
    Nugget <- ranPars$Nugget
    if (! is.null(Nugget)) trueCorrpars$Nugget <- Nugget 
    msd.arglist <- list(rho=rho)
    datanames <- names(data)
    ## in a typical call from corrHLfit the following test should be FALSE
      if ( (length(rho)>1 || missing(distMatrix)) && missing(uniqueGeo)) { ## all cases where we need uniqueGeo
        if (corr.model=="adjacency") {
          if ( missing(adjMatrix) ) stop("missing 'adjMatrix' for adjacency model")
        } else {
          if ( ! is.null(spatial.model)) {
            bars <- spatial.model[[2]] 
            coordinates <- deparse(bars[[3]]) ## "x + y"
            coordinates <-  strsplit(coordinates," ")[[1]]
            coordinates <- coordinates[coordinates != "+"]
          } else {
            stop("Call to 'HLCor' without a spatial term in the formula is suspect.")
            ## very old code handling old syntax with (1|pos) and default values of the coordinates argument
            coordinates <- c("x","y") ## back compat
          }
          coordcols <- which(datanames %in% coordinates)
          if ( length(coordcols) != length(coordinates) ) {
            stop("variables 'coordinates' not all found in the 'data'")
          }
        } ## 
        #### then uniqueGeo
        uniqueGeo <- unique(data[,coordinates,drop=F]) ## keeps the names of first instances of the coordinates in data
      } 
      ## then compute scaled distances from unscaled info, for HLfit call
      if (length(rho)>1) {
        msd.arglist <- c(msd.arglist,list(uniqueGeo=uniqueGeo))
      } else {
        if ( missing(distMatrix)) { ## we need distMatrix here, even if length(rho)>1
          ##        distMatrix <- as.matrix(dist(uniqueGeo))
          distMatrix <- dist(uniqueGeo)
        }
        msd.arglist <- c(msd.arglist,list(distMatrix=distMatrix))
      }
      scaled.dist <- do.call("make.scaled.dist",msd.arglist)
      nunu <- trueCorrpars;nunu$rho <- NULL ## because the Matern.corr input will be an already scaled distance
      m <- do.call(Matern.corr,args=c(nunu,list(d=scaled.dist)))  
    } else if (corr.model %in% c("BreslowC93","adjacency")) { 
      m<-solve(diag(rep(1,nrow(adjMatrix)))-rho*(adjMatrix))
    }  
    if (verbose[1]) { 
      print(unlist(trueCorrpars))
    }
    # print(paste("Correlation params in HLCor",paste(c(rho,nu),collapse=" ")))
    argsfordesignL <- dotlist[intersect(names(dotlist),names(formals(designL.from.Corr)))] 
    if (class(m)=="dist") {
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
      attr(predictor,"LMatrix") <- Lunique
    }
  } else { ## we have a corrMatrix
    Lunique <- designL.from.Corr(corrMatrix)
  }  
  if (length(verbose)>0) {
    dotlist$verbose <- verbose[-1] ## all printing in HLfit is suppressed by default
  }
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
    dotlist$ranFix$lambda <- dotlist$LamFix
    dotlist$LamFix <- NULL
  }
  if (!is.null(dotlist$PhiFix)) {
    stop("argument PhiFix of HLCor is obsolete")
    dotlist$ranFix$phi <- dotlist$PhiFix
    dotlist$PhiFix <- NULL
  }  
##
  ## allows log and not log:
  HL.info$ranFix <- ranPars  ## dispersion and correlation params, will permit rho,nu,... in the HLfit results in the same way that $phi... 
                                   ## anticipates the day where corrpars will be treated as disp pars
                                   ## and simplifies the corrMM.LRT code
  if (!is.null(ranPars$trPhi)) {
    HL.info$ranFix$phi <- dispInv(ranPars$trPhi)
    ranPars$trPhi <- NULL
  } else if (!is.null(ranPars$logphi)) { ## debug code
    ## HL.info$ranFix$phi <- exp(ranPars$logphi)
    stop("logphi in HLCor...")
  } #####################  else HL.info$ranFix$phi <- ranPars$phi ## y st deja !?
  if (!is.null(ranPars$trLambda)) {## 
    HL.info$ranFix$lambda <- dispInv(ranPars$trLambda)
    ranPars$trLambda <- NULL
  } else if (!is.null(ranPars$loglambda)) { ## debug code
    ## HL.info$ranFix$lambda <- exp(ranPars$loglambda)
    stop("loglambda in HLCor...")
  } ##################### else HL.info$ranFix$lambda <- ranPars$lambda
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
  ## we probably do not need to keep hlfit$call (cf resu$call below)
  hlfit$call <- NULL
  if (verbose[1]) { ## useful in final call from corrHLfit
    summary(hlfit) ## input corr pars have been printed at the beginning...   
  }
  resu <-list(hlfit=hlfit)
  HLCor.info <-list(coordinates=coordinates) ## more spatial info is to be found in hlfit$predictor (Lunique = corrmat^1/2) and hlfit$ZALMatrix
  resu$HLCor.info <- HLCor.info
  resu$call <- mc
  class(resu) <- c("HLCor",class(hlfit)) 
  return(resu) ## 
}
