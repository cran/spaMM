`%id*%` <- function(A,b) {
  if (is.identity(A)) return(b)
  return(A %*% b)
}


`%id*id%` <- function(A,B) {
  if (is.identity(A)) return(B)
  if (is.identity(B)) return(A)
  return(A %*% B)
}

calc_beta_w_cov <- function(object) {
  beta_w_cov <- attr(object$beta_cov,"beta_v_cov")
  invL <- calc_invL(object)
  if ( ! is.null(invL)) {
    pforpv <- ncol(object$beta_cov)
    v.range <- pforpv+seq(ncol(invL))
    beta_w_cov[v.range,] <- t(invL) %*% beta_w_cov[v.range,]
    beta_w_cov[,v.range] <- beta_w_cov[,v.range] %*% invL
  }
  return(beta_w_cov)
}

calcPredVar <- function(Coldnew,X.pv,newZAC=NULL,newZA,beta_w_cov,logdispObject=NULL,covMatrix=FALSE,blockSize=100L) {
  nrX <-  nrow(X.pv)
  if ( ( ! covMatrix ) && nrX > blockSize) {
    slices <- unique(c(seq(0L,nrX,blockSize),nrX))
    sliceVar <- function(it) {
      slice <- (slices[it]+1L):slices[it+1L]
      if (is.list(newZA)) {
        nrand <- length(newZA) ## or of any other of the lists of matrices
        newZAClist <- lapply(seq_len(nrand), function(it) {
          terme <- newZA[[it]][slice,,drop=FALSE] %id*id% t(Coldnew[[it]])[] 
          as.matrix(terme) ## loses names but they are not useful here 
        })
        newZAC <- do.call(cbind,newZAClist)
      } else newZAC <- newZA[slice,,drop=FALSE] %id*id% t(Coldnew)[]
      calcPredVar(X.pv=X.pv[slice,,drop=FALSE],newZAC=newZAC,
                  beta_w_cov=beta_w_cov,logdispObject=NULL,covMatrix=covMatrix,blockSize=blockSize)
    }
    unlist(sapply(seq_len(length(slices)-1L),sliceVar))
  }
  if (is.null(newZAC)) {
    if (is.list(newZA)) {
      nrand <- length(newZA) ## or of any other of the lists of matrices
      newZAClist <- lapply(seq_len(nrand), function(it) {
        terme <- newZA[[it]] %id*id% t(Coldnew[[it]])[] ## %id*id% further keeps col names of C if tZA==I
        as.matrix(terme) ## loses names but they are not useful here 
      })
      newZAC <- do.call(cbind,newZAClist)
    } else newZAC <- newZA %id*id% t(Coldnew)[]
  }
  newAugX <- CBIND(X.pv,newZAC) 
  if (covMatrix) {
    predVar <- newAugX %*% beta_w_cov %*% t(newAugX)
  } else {
    premul <- newAugX %*% beta_w_cov
    predVar <- rowSums(premul * newAugX)
  }
  if (! is.null(logdispObject)) {
    newZACw <- newZAC %*% logdispObject$dwdlogdisp ## typically nnew * 2 hence small 
    if (covMatrix) {
      disp_effect_on_newZACw <- newZACw %*% logdispObject$logdisp_cov %*% t(newZACw)  
    } else {
      premul <- newZACw %*% logdispObject$logdisp_cov
      disp_effect_on_newZACw <- rowSums(premul * newZACw)
    }
    predVar <- predVar + disp_effect_on_newZACw
  }
  return(predVar) ## may be a Matrix
}




calcResidVar <- function(object,newdata=NULL) {
  phi.object <- object$phi.object
  if ( ! is.null(phi.object$phi.Fix)) {
    if (length(phi.object$phi.Fix)==1L) {
      residVar <- rep(phi.object$phi.Fix,nrow(newdata))
    } else if (is.null(newdata)) {
      residVar <- rep(phi.object$phi.Fix)          
    } else {
      stop("Unable to compute 'residVar' given 'newdata' argument and given vector of phi values in fit object.")
    }
  } else {
    ## creates pseudo HLfit object
    phipred <- object$resid.predictor
    if (!inherits(phipred,"predictor")) phipred <- Predictor(phipred) ## but the offset is not always set...
    pseudoObject <- list(fixef=object$phi.object$beta_phi, predictor=phipred,family=object$resid.family)
    residVar <- predict.HLfit(pseudoObject, newdata=newdata,
                              residVar=FALSE) ## avoid recursive call        
    attributes(residVar) <- NULL
  }
  residVar
}  

calcNewCorrs <- function(object,locdata,#predVar,
                         spatial.model) {
  olduniqueGeo <- attr(object,"info.uniqueGeo")
  geonames <- colnames(olduniqueGeo)
  newuniqueGeo <- calcUniqueGeo(data=locdata[,geonames,drop=FALSE]) ## It is essential that it keeps the same order as spMMfactorlist -> ULI -> unique. 
  ### rho only used to compute scaled distances
  rho <- getPar(object$ranFix,"rho")
  if( !is.null(rho.mapping <- object$control.dist$rho.mapping)) rho <- rho[rho.mapping] # has been missing
  ## rows from newuniqueGeo, cols from olduniqueGeo:
  msd.arglist <- list(uniqueGeo=newuniqueGeo,uniqueGeo2=olduniqueGeo,rho=rho,return_matrix=TRUE)
  if ( ! is.null(dist.method <- object$control.dist$dist.method)) msd.arglist$dist.method <- dist.method
  uuCnewold <- do.call(make_scaled_dist,msd.arglist) ## ultimately allows products with Matrix ## '*cross*dist' has few methods, not even as.matrix
#  if (predVar)  {
#    msd.arglist$uniqueGeo2 <- NULL
#    Cnewnew <- do.call(make_scaled_dist,msd.arglist) ## dist matrix and then call to correl fn
#    Cnewnew <- as.matrix(Cnewnew) 
#  }
  if (spatial.model[[1]]=="AR1") {
    args <-object$ranFix[which(names(object$ranFix) %in% c("ARphi"))]
    uuCnewold <- args$ARphi^uuCnewold  
#    if (predVar) Cnewnew <- args$ARphi^Cnewnew  
  } else {
    args <-object$ranFix[which(names(object$ranFix) %in% c("nu","Nugget"))] ## so that rho=1 in MaternCorr
    uuCnewold <- do.call(MaternCorr,args=c(args,list(d=uuCnewold)))  
#    if (predVar) Cnewnew <- do.call(MaternCorr,args=c(args,list(d=Cnewnew)))  
  }
  return(uuCnewold)
}

makenewname <- function(base,varnames) { ## FR->FR post CRAN 1.4.1
  varnames <- varnames[which(substring(varnames,1,nchar(base))==base)] 
  allremainders <- substring(varnames,nchar(base)+1) 
  allnumericremainders <- as.numeric(allremainders[which( ! is.na(as.numeric(allremainders )))  ]) ## as.numeric("...")
  ## FR->FR 2015/03/04
  if (length(allremainders) == 0L && length(allnumericremainders) == 0L) { ## if base = allremainders => length(allnumericremainders) == 0 not sufficient
    fvname <- base
  } else {
    num <- max(allnumericremainders)+1
    fvname <-paste ( base , num , sep="") 
  }
  fvname
}

## (1) for surface prediction: (developed in InferentialSimulation/InferentialSimulation.R)
## (2) But also for generation of fixed effects in simulation of nested-effects models
predict.HLfit <- function(object,newdata=newX,newX = NULL,re.form = NULL,
                          variances=list(fixef=FALSE,linPred=FALSE,dispVar=FALSE,resid=FALSE,sum=FALSE,cov=FALSE),
                          predVar=variances$linPred,residVar=variances$resid,
                          binding = FALSE, 
                          ...) { ## but not new Y
  if ( ! is.null(variances$ranef)) {
    message("'variances$ranef' is obsolete: I interpret this as 'variances$linPred'.")
    variances$linPred <- variances$ranef
    variances$ranef <- NULL
  }
  if (is.null(variances$sum)) variances$sum <- FALSE ## transient value; so that the following lines do not leave the lhs as NULL
  if (is.null(variances$linPred)) variances$linPred <- variances$sum
  if (is.null(variances$resid)) variances$resid <- variances$sum
  if (is.null(variances$fixef)) variances$fixef <- FALSE ## 12/2014 post 1.4.4
  if (is.null(variances$cov)) variances$cov <- FALSE
  if (is.null(variances$dispVar)) variances$dispVar <- FALSE
  if ( ! is.null(residVar) && residVar) variances$resid <- TRUE
  if (predVar == "Cov") {
    variances$linPred <- TRUE; variances$cov <- TRUE
  } else if ( ! is.null(predVar) && predVar) variances$linPred <- TRUE
  variances$sum <- (variances$linPred && variances$resid)
  ##
  locform <- attr(object$predictor,"oriFormula")
  ## possible change of random effect terms
  if (noReForm(re.form)) {
    locform <- nobarsMM(locform)  ## ie. (if re.form implies that there is no random effect)
  } else if (inherits(re.form,"formula")) locform <- re.form
  # checking variables in the data
  if (length(locform)==3L) locform <- locform[-2] ## removes RHS for checking  vars on RHS etc
  allvars <- all.vars(locform) 
  if (is.vector(newdata)) { ## ## less well controlled case, but useful for maximization
    binding <- binding ## :-) binding must be evaluated before newdata is modified
    newdata <- data.frame(matrix(newdata,nrow=1))
    if (length(allvars)==ncol(newdata)) {
      names(newdata) <- allvars
    } else {
      stop(paste("(!) newdata has incorrect length. It should match the following variables:\n",paste(allvars,collapse=" ")))
    }
  } 
  ## it is important that newdata remains NULL if it was so initially because it is tested below. Hence copy in locdata
  if (is.null(newdata)) {
    locdata <- object$data[,allvars,drop=FALSE]
  } else {
    if( is.matrix(newdata) ) newdata <- as.data.frame(newdata)  
    # so that matrix 'newdata' arguments can be used as in some other predict methods.
    locdata <- newdata ## FR->FR [,allvars,drop=FALSE] ?
    if (any(is.na(locdata))) {
      stop("NA's in required variables from 'newdata'. Prediction not possible.")
    }
  }
  ## preparation for fixed effects
  allFrames <- HLframes(formula=locform,data=locdata)
  newX.pv <- allFrames$X
  if ( ! is.null(betaFix <- attr(object$predictor,"offsetObj")$betaFix)) { ## suppress betaFix cols so that this is consistent with <object>$X.pv 
    newX.pv <- newX.pv[,colnames(object$`X.pv`),drop=FALSE]
  }
  needNewEta <- ( ( ! is.null(newdata) ) || variances$linPred || ! is.null(re.form))
  if (needNewEta) etaFix <- newetaFix(object,allFrames)  ## new fixed effects (or [part of] new linear predictor if re.form)      
  ## preparation for random effects
  spatial.terms <- findSpatial(locform) ## list of formula terms
  spatial.model <- spatial.terms[[1]] ## one formula term, e.g Matern(1|...)
  if ( ! is.null(spatial.model)) { 
    if (! is.null(newdata) && as.character(spatial.model[[1]]) %in% c("adjacency","ar1")) {
      stop("Prediction in newdata not implemented or not possible in the 'adjacency' model")
    } ## FR->FR would be possible for new non-spatial predictor values in the original locations... il faudrait un test sur les elements de la distance matrix
  } 
  oldLMatrix <- attr(object$predictor,"LMatrix") ## may be NULL
  if ( ! is.null(spatial.model) && needNewEta ) { ## all cases where we need distance matrices
    uuCnewold <- calcNewCorrs(object=object,locdata=locdata,#predVar=variances$linPred,
                              spatial.model=spatial.model)
    attr(uuCnewold,"ranefs") <- attr(oldLMatrix,"ranefs")  
  } 
  ## matching ranef terms of re.form
  if (noReForm(re.form)) {
    nrand <- 0
  } else {
    if (inherits(re.form,"formula")) {
      ranefs <- parseBars(locform)
      nrand <- length(ranefs)
      newinold <- sapply(ranefs, function(v) {which(v==attr(object$ZAlist,"ranefs"))})
      subZAlist <- object$ZAlist[newinold] ## and reordered
      sublambda <- object$lambda[newinold]
    } else {
      ranefs <- attr(object$ZAlist,"ranefs")
      nrand <- length(ranefs)
      newinold <- seq(nrand)
      subZAlist <- object$ZAlist
      sublambda <- object$lambda
    }    
    spatialOne <- which(ranefs == spatial.model) ## strictly a spatial one, not other correlated ones      
  }
  ## (1) computes fv (2) compute predVar
  ##### fv
  if (noReForm(re.form)) {
    fv <- object$family$linkinv(etaFix) 
  } else if ( is.null(newdata) && ! inherits(re.form,"formula")) {
    fv <- object$fv ## same length including replicates
    newZAlist <- subZAlist ## useful if predVar
  } else { ## 
    if ( nrand==0L ) {
      eta <- etaFix
      newZAlist <- NULL
    } else {
      FL <- spMMFactorList(locform, allFrames$mf, 0L, drop=TRUE) 
      newZAlist <- FL$Design ## must be ordered as parseBars result...
      attr(newZAlist,"ranefs") <- ranefs ## utile pour computeZAXlist to match the ranefs of LMatrix
      if ( is.null(spatial.model)) {
        newZAClist <- computeZAXlist(XMatrix=NULL,ZAlist=newZAlist) 
      } else {
        newZAClist <- computeZAXlist(XMatrix=uuCnewold,ZAlist=newZAlist) ## ZAL's for ZA's and L's (typically some ZA's are unaffected)
      }  
      #### precomputation of coeffs
      ## on the gaussian scale, L.v_ori ~ lam C (lam C + phi I)^{-1}y 
      ## new random autocorr term ~ lam c (lam C + phi I)^{-1}y = c C^{-1} L_ori.v_ori = c [t(L_ori)]^{-1} v_ori
      ## [t(L_ori)]^{-1} v_ori can be computed once for all predictions => 'w_h_coeffs'
      w_h_coeffs <- object$get_w_h_coeffs() ## should work for nonspatial models 
      for (it in seq_len(nrand)) { 
        if ( inherits(newZAClist[[it]],"Matrix")) { ## this is not tested in compute ZALlist because the latter does not runs other all list elements, only those with a matching L
          cnames <- colnames(newZAClist[[it]]) ## used a few lines belwo to match levels of ranef
          rnames <- rownames(newZAClist[[it]]) ## used a few lines belwo to match levels of ranef
          newZAClist[[it]] <- as.matrix(newZAClist[[it]]) ## pfff  loses names
          colnames(newZAClist[[it]]) <- cnames
          rownames(newZAClist[[it]]) <- rnames
        }
      }
      old_cum_n_u_h <- attr(object$lambda,"cum_n_u_h")
      lcrandfamfam <- attr(object$`rand.families`,"lcrandfamfam")[newinold]
      augm_w_h_coeffs <- lapply(seq_len(nrand),function(it) {
        oldu.range <- (old_cum_n_u_h[newinold[it]]+1L):(old_cum_n_u_h[newinold[it]+1L])
        if (it %in%  spatialOne) {     # %in% handles zero-length spatialOne...
          return(w_h_coeffs[oldu.range])          
        } else {
          oldlevels <- colnames(subZAlist[[it]])
          newlevels <- colnames(newZAClist[[it]])
          interlevels <- intersect(oldlevels,newlevels)
          oldv <- w_h_coeffs[oldu.range]
          names(oldv) <- oldlevels
          psi_M <- switch(lcrandfamfam[it], 
                           gaussian = 0,
                           gamma = 1, 
                           beta = 1/2, 
                           "inverse.gamma" = 1
          )
          vpsi_M <- object$rand.families[[newinold[it]]]$linkfun(psi_M) 
          ## since vpsi_M can be non-zero, the expectation of the response can be modified in a re.form model compared to the original
          newv <- rep(vpsi_M,length(newlevels))
          names(newv) <- newlevels
          newv[interlevels] <- oldv[interlevels] 
          return(newv)
        }
      })
      if (nrand>1L) {
        ZACw <- sapply(seq_len(nrand),function(it) {newZAClist[[it]][] %*% augm_w_h_coeffs[[it]]}) ## matrix
        ZACw <- rowSums(ZACw)
      } else ZACw <- newZAClist[[1]][] %*% augm_w_h_coeffs[[1]]
      eta <- etaFix + ZACw ## (length(eta)) col vector from coeffs = length(eta) row vector...
    }
    # done with eta
    fv <- object$family$linkinv(eta) ## ! freqs for binomial, counts for poisson
  }
  resu <- fv ## suitable for objective function of optim() etc ## matrix ! maybe more suitable than data frame as objective function
  if ( ! is.logical(binding) ) { ## expecting a string
    binding <- makenewname(base=binding,varnames=colnames(locdata)) ## FR->FR 09/11/2014 = posterior to CRAN 1.4.1 
    resu <- data.frame(resu)
    colnames(resu) <- binding
    resu <- cbind(locdata,resu) ## becomes a data frame !
    attr(resu,"fittedName") <- binding
  } else { ## expecting binding= FALSE
    if (ncol(locdata)>0)  attr(resu,"frame") <- locdata 
  }
  ##### (2) predVar
  if(variances$linPred) {
    beta_w_cov <- object$get_beta_w_cov()
    ## list for Cnewnew, which enters in  newZA %*% Cnewnew %*% tnewZA, hence should not represent newZA itself 
    #     if (nrand>0L) newnewClist <- lapply(seq_len(nrand),function(it) {
    #       if ( it %in% spatialOne) { ## avoids running the next algo which is slow on large matrices
    #         Cnn <- Cnewnew ## already computed for point prediction
    #       } else { 
    #         Cnn <- Diagonal(ncol(newZAlist[[it]])) ## diag(ncol(newZAlist[[it]]))
    #       }
    #       return(Cnn)
    #     }) ## => completely ignores effects removed in re.form 
    ## list for Coldnew, which enters in ZA %id*id% Coldnew %id*id% tnewZA 
    if (nrand>0L) oldnewClist <- lapply(seq_len(nrand),function(it) {
      if (it %in% spatialOne) { ## avoids running the next algo which is slow on large matrices
        oldnewC <- t(uuCnewold)
      } else {
        oldlevels <- colnames(subZAlist[[it]])
        newlevels <- colnames(newZAlist[[it]])
        if (identical(oldlevels,newlevels)) {
          oldnewC <- Diagonal(length(oldlevels)) ## replaces old identityMatrix
        } else {
          oldornew <- unique(c(oldlevels,newlevels))
          oldnewC <- diag(length(oldornew))
          colnames(oldnewC) <- rownames(oldnewC) <- oldornew
          oldnewC <- oldnewC[oldlevels,newlevels]
        }
      }
      return(oldnewC)
    })
    if (nrand>0L) {
      loclist <- list(X.pv=newX.pv,beta_w_cov=beta_w_cov,covMatrix=variances$cov)
      if (nrand==1L) {
        loclist$Coldnew=oldnewClist[[1]]
        loclist$newZA=newZAlist[[1]]
        ## but the code for nrand >1 should work for nrand==1:
      } else {
        loclist$Coldnew=oldnewClist
        loclist$newZA=newZAlist
      }
      if (variances$dispVar) loclist$logdispObject <- object$logdispObject
      if (variances$cov) {
        sumVar <- as.matrix(do.call(calcPredVar,loclist)) ## matrix, not Matrix (assumed below)
        rownames(sumVar) <- colnames(sumVar) <- rownames(locdata)
      } else {
        sumVar <- do.call(calcPredVar,loclist) 
        names(sumVar) <- rownames(locdata)
      }
    } else sumVar <- matrix(0,nrow=nrow(locdata),ncol=nrow(locdata))
    attr(resu,"predVar") <- sumVar ## vector or matrix
  } else sumVar <- 0  
  if (variances$resid) {
    if (object$family$family %in% c("poisson","binomial")) {
      attr(resu,"residVar") <- object$family$variance(fv)
    } else attr(resu,"residVar") <- calcResidVar(object,newdata=locdata) 
    if (inherits(sumVar,"matrix")) {
      diag(sumVar) <- diag(sumVar) + attr(resu,"residVar")
    } else sumVar <- sumVar + attr(resu,"residVar")
  }
  if ( variances$fixef && ! is.null(object$beta_cov)) {
    fixefcov <- newX.pv %*% object$beta_cov %*% t(newX.pv)
    if (variances$cov) {
      attr(resu,"fixefVar") <- fixefcov 
    } else attr(resu,"fixefVar") <- diag(fixefcov)
    ## sumVar <- sumVar + attr(resu,"fixefVar") ## there is already such a term in predVar
  }
  if (variances$sum) attr(resu,"sumVar") <- sumVar
  class(resu) <- c(class(resu),"predictions") ## for print method
  return(resu)
}

print.predictions <- function (x, expanded=FALSE, ...) 
{
  asvec <- as.vector(x)
  rnames <- rownames(x)
  toolong <- nchar(rnames)>9
  rnames[toolong] <- paste(substr(rnames[toolong],0,8),".",sep="")
  names(asvec) <- rnames
  print(asvec)
  cat("*stored as* 1-col matrix with attributes:")
  std.attr <- c("names","dim","dimnames") ## attributes not to be shown
  a <- attributes(x)
  nam <- names(a)
  if (expanded) { # shows structure of attributes as in utils:::str.default
    cat("\n")
    nest.lev <- 0
    indent.str <- paste(rep.int(" ", max(0, nest.lev + 1)), collapse = "..")
    strO <- getOption("str")
    strO <- modifyList(strOptions(), strO) ## seems to provide a format.fun function
    `%w/o%` <- function(x, y) x[is.na(match(x, y))]
    nfS <- names(fStr <- formals())
    ## this scans the substructure of each attribute
    strSub <- function(obj, ...) {
      nf <- nfS %w/o% c("object", "give.length", "comp.str", 
                        "no.list", names(match.call())[-(1:2)], "...")
      aList <- as.list(fStr)[nf]
      aList[] <- lapply(nf, function(n) eval(as.name(n)))
      strObj <- function(...) str(obj, ...)
      do.call(strObj, c(aList, list(...)), quote = TRUE)
    }
    for (i in seq_along(a)) if (all(nam[i] != std.attr)) {
      cat(indent.str, paste0("- attr(*, \"", nam[i], "\")="), 
          sep = "")
      strSub(a[[i]], give.length = TRUE, indent.str = paste(indent.str, 
                                                            ".."), nest.lev = nest.lev + 1)
    }
  } else {
    cat(" ")
    nam <- setdiff(nam,std.attr)
    cat(paste(nam,collapse=", "))  
    cat("\n")
  }
  invisible()
}