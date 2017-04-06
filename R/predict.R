`%id*%` <- function(A,b) {
  if (is.identity(A)) return(b)
  return(A %*% b)
}

`%*id%` <- function(a,B) {
  if (is.identity(B)) return(a)
  return(a %*% B)
}


`%id*id%` <- function(A,B) {
  if (is.identity(A)) return(B)
  if (is.identity(B)) return(A)
  return(A %*% B)
}

calcPredVar <- function(ColdnewList,X.pv,newZAC=NULL,newZAlist,beta_w_cov,
                        CnewnewList=NULL,
                        invColdoldList=NULL, ## may be null is no newdata.
                        lambda=NULL, ## may be null if no newdata.
                        logdispObject=NULL, ## should remain NULL is disp not requested
                        covMatrix=FALSE,blockSize=100L) {
  nrX <-  nrow(X.pv)
  if ( ( ! covMatrix ) && nrX > blockSize) {
    slices <- unique(c(seq(0L,nrX,blockSize),nrX))
    sliceVar <- function(it) {
      slice <- (slices[it]+1L):slices[it+1L]
      ## here  the problem is that newZA should map the new response levels 
      ## to the 'new' levels of random effects  
      nrand <- length(newZAlist) ## or of any other of the lists of matrices
      requiredLevelsList <- lapply(seq_len(nrand), function(rd) {
        which(colSums(newZAlist[[rd]][slice,,drop=FALSE])>0L) })
      locnewZA <- lapply(seq_len(nrand), function(rd) { 
        newZAlist[[rd]][slice,requiredLevelsList[[rd]],drop=FALSE] })
      # predVar in observed points uses C rather than L hence we need to compute ZA.C in all cases ('newZAC")
      #  andthen we need newZA and ColdnewList
      locColdnewList <- lapply(seq_len(nrand), function(rd) { 
        locC_on <- ColdnewList[[rd]][,requiredLevelsList[[rd]],drop=FALSE]
        attr(locC_on,"isEachNewInOld") <- attr(ColdnewList[[rd]],"isEachNewInOld")[requiredLevelsList[[rd]]] ## for non-spatial effects; (qualifies sub cols of sub Cnewold)
        return(locC_on)
      })
      locnewZAClist <- lapply(seq_len(nrand), function(rd) { 
        terme <- locnewZA[[rd]] %id*id% t(locColdnewList[[rd]])[] 
        as.matrix(terme) ## loses names but they are not useful here 
      })
      if (nrand>1L) {newZAC <- do.call(cbind,locnewZAClist)} else {newZAC <- locnewZAClist[[1]]}
      ## the only one that corresponds to real newdata in predict()
      if (!is.null(CnewnewList)) {
        locCnewnewList <- lapply(seq_len(nrand), function(rd) { 
          CnewnewList[[rd]][requiredLevelsList[[rd]],requiredLevelsList[[rd]],drop=FALSE] }) 
      } else locCnewnewList <- NULL
      calcPredVar(ColdnewList=locColdnewList, ## needed only if newdata in predict() call
                  X.pv=X.pv[slice,,drop=FALSE],## problem is that this creates the apparence of new data end more calculations
                  newZAC=newZAC,
                  newZAlist=locnewZA, ## either newZA or newZAC needed even if no newdata in predict() call
                  beta_w_cov=beta_w_cov,CnewnewList=locCnewnewList,
                  invColdoldList=invColdoldList, ## needed only if newdata in predict() call
                  lambda=lambda,
                  logdispObject=logdispObject,
                  covMatrix=covMatrix,blockSize=blockSize)
    }
    unlist(sapply(seq_len(length(slices)-1L),sliceVar))
  }
  ############################################################
  if (is.null(newZAC)) {
    nrand <- length(newZAlist) ## or of any other of the lists of matrices
    newZAClist <- lapply(seq_len(nrand), function(it) {
      terme <- newZAlist[[it]] %id*id% t(ColdnewList[[it]])[] ## %id*id% further keeps col names of C if tZA==I
      as.matrix(terme) ## loses names but they are not useful here 
    })
    if (nrand>1L) {newZAC <- do.call(cbind,newZAClist)} else {newZAC <- newZAClist[[1]]}
  }
  newAugX <- cbind2(X.pv,newZAC) ## mais en fait pas un AugX since it uses C (in C.w) rather than L (in L.v)
  ## First component of predVar
  # variance of expectation of Xbeta+Zb due to var of (hat(beta),hat(v)) using E[b] as function of hat(v)
  calcZWZt_mat_or_diag <- function(Z,W,returnMat) { ## fixefVar or fixefVar + a bit of ranefVar
    if (returnMat) {
      return(Z[] %id*id% W[] %id*id% t(Z)[])
    } else { ## returns only the diagonal
      if (is.vector(W)) {
        return(rowSums(Matrix_times_Dvec(Z[],W[]) * Z[]))
      } else {
        rownames(W) <- colnames(W) <- NULL ## inhibits a useless warning from Matrix:: dimNamesCheck 
        premul <- Z[] %id*id% W[]
        return(rowSums(suppressMessages(premul * Z[]))) ## suppress message("method with signature...") [found by debug(message)]
      }
    }
  }
  predVar <- calcZWZt_mat_or_diag(newAugX,beta_w_cov,covMatrix)
  ## Second component of predVar:
  # Evar: expect over distrib of (hat(beta),hat(v)) of [variance of Xbeta+Zb given (hat(beta),hat(v))]
  if ( ! is.null(invColdoldList)) { ## must be equivalent to the presence of newdata
    nrand <- length(newZAlist) ## or of any other of the lists of matrices
    Evarlist <- lapply(seq_len(nrand), function(it) {
      isEachNewInOld <- attr(ColdnewList[[it]],"isEachNewInOld")
      if ( ! is.null(isEachNewInOld)) { ## non spatial effect: a vector of booleans indicating whether new is in old (qualifies cols of Cnewold)
        Evar <- Diagonal(x=lambda[it]*as.numeric(! isEachNewInOld))
      } else { ## spatial effect
        if (is.null(CnewnewList[[it]])) { ## we compute only the var vector
          Cno_InvCoo_Con <- colSums(ColdnewList[[it]] * ( (invColdoldList[[it]])[] %id*id% (ColdnewList[[it]])[] ) )
          Evar <- lambda[it] * (1 - Cno_InvCoo_Con)
        } else { ## full cov matrix
          Cno_InvCoo_Con <- t(ColdnewList[[it]])[] %id*id% (invColdoldList[[it]])[] %id*id% (ColdnewList[[it]])[]
          Evar <- lambda[it] * (CnewnewList[[it]] - Cno_InvCoo_Con)
        }
      } 
      terme <- calcZWZt_mat_or_diag(newZAlist[[it]],Evar,covMatrix)
      if (covMatrix) {
        return(as.matrix(terme))  ## loses names but they are not useful here
      } else return(terme)
    })
    if (nrand>1L) {Evar <- Reduce("+",Evarlist)} else {Evar <- Evarlist[[1]]}
    predVar <- predVar + Evar
  }
  # If components for uncertainty in dispersion params were requested,
  #   logdispObject is not NULL
  # If some components ere computable, $$dwdlogdisp should not be NULL
  # Former approach (changed 08/2016) was to test logdispObject and then 
  #   for any 'problem'. But there may some 'problem' and still a valid logdispObject
  # => In new version, dwdlogdisp should be either NULL or a conforming matrix;
  #  'problems" should not be tested.
  if ( ! is.null(logdispObject$dwdlogdisp) ) {
    newZACw <- newZAC %*% logdispObject$dwdlogdisp ## typically (nnew * n_u_h) %*% (n_u_h * 2) = nnew * 2 hence small 
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
  if (is.null(phi_outer <- phi.object$phi_outer)) { ## valid whether newdata are NULL or not:
    glm_phi <- phi.object[["glm_phi"]]
    if (is.null(glm_phi)) glm_phi <- .get_glm_phi(object)
    residVar <- predict(glm_phi, newdata=newdata, type="response")
  } else { ## phi, but not glm_phi
    if (length(phi_outer)==1L) {
      if (is.null(newdata)) {
        residVar <- rep(phi_outer,nrow(object$X.pv)) ## assumes (length(phi_outer)==1L)           
      } else residVar <- rep(phi_outer,nrow(newdata))
    } else stop("Unable to compute 'residVar' given length(phi_outer)!=1L.") ## and no glm_phi
    ## FR->FR we could improve this if we get a glm_phi when phi was estimed by outer iterations
  }
  residVar
}  

calcNewCorrs <- function(object,locdata,which,
                         spatial.model) {
  resu <- list(uuCnewold=NULL,uuCnewnew=NULL)
  if (any(unlist(which))) {
    # (1) code  common to different ranef models 
    olduniqueGeo <- attr(object,"info.uniqueGeo")
    geonames <- colnames(olduniqueGeo)
    newuniqueGeo <- calcUniqueGeo(data=locdata[,geonames,drop=FALSE]) ## It is essential that it keeps the same order as .spMMfactorlist -> ULI -> unique. 
    ### rho only used to compute scaled distances
    rho <- getPar(object$ranFix,"rho")
    if ( ! is.null(rho_mapping <- attr(object,"dist_info")$rho.mapping)
        && length(rho)>1L ) rho <- fullrho(rho=rho,coordinates=geonames,rho_mapping=rho_mapping)
    ## rows from newuniqueGeo, cols from olduniqueGeo:
    msd.arglist <- list(uniqueGeo=newuniqueGeo,uniqueGeo2=olduniqueGeo,
                        rho=rho,return_matrix=TRUE)
    if ( ! is.null(dist.method <- object$control.dist$dist.method)) msd.arglist$dist.method <- dist.method
    if (which$no) resu$uuCnewold <- do.call(make_scaled_dist,msd.arglist) ## ultimately allows products with Matrix ## '*cross*dist' has few methods, not even as.matrix
    if (which$nn)  {
      msd.arglist$uniqueGeo2 <- NULL
      if (nrow(msd.arglist$uniqueGeo)==1L) {
        resu$uuCnewnew <- matrix(0)
      } else resu$uuCnewnew <- do.call(make_scaled_dist,msd.arglist) 
    }
    ## distance matrix and then call to correl fn:
    # (2) code specific to each ranef model
    if ( ! is.null(spatial.model)) {
      resu <- .calc_corr_from_dist(resu, object, spatial.model, which)
    }
  }
  return(resu)
}


makenewname <- function(base,varnames) { ## post CRAN 1.4.1
  varnames <- varnames[which(substring(varnames,1,nchar(base))==base)] 
  allremainders <- substring(varnames,nchar(base)+1) 
  allnumericremainders <- as.numeric(allremainders[which( ! is.na(as.numeric(allremainders )))  ]) ## as.numeric("...")
  ## 2015/03/04
  if (length(allremainders) == 0L && length(allnumericremainders) == 0L) { ## if base = allremainders => length(allnumericremainders) == 0 not sufficient
    fvname <- base
  } else {
    if (length(allnumericremainders)>0L) {
      num <- max(allnumericremainders)+1L
    } else num <- 1L
    fvname <-paste ( base , num , sep="") 
  }
  fvname
}

.process_variances <- function(variances) {
  # $respVar implies components
  if (is.null(variances$predVar)) variances$predVar <- variances$respVar ## may still be NULL
  if (is.null(variances$residVar)) variances$residVar <- variances$respVar ## may still be NULL
  if (is.null(variances$respVar)) variances$respVar <- FALSE 
  # $predVar implies components
  if (is.null(variances$linPred)) variances$linPred <- variances$predVar ## may still be NULL
  if (is.null(variances$disp)) variances$disp <- variances$predVar ## may still be NULL
  if (is.null(variances$predVar)) variances$predVar <- FALSE 
  # Do not let any component empty
  if (is.null(variances$fixefVar)) variances$fixefVar <- FALSE 
  if (is.null(variances$linPred)) variances$linPred <- FALSE 
  if (is.null(variances$disp)) variances$disp <- FALSE ## uncertaintly on dispersion parameters
  if (is.null(variances$residVar)) variances$residVar <- FALSE ## uncertaintly on dispersion parameters
  if (is.null(variances$cov)) variances$cov <- FALSE
  ##
  return(variances)
}

## (1) for surface prediction: (developed in InferentialSimulation/InferentialSimulation.R)
## (2) But also for generation of fixed effects in simulation of nested-effects models
predict.HLfit <- function(object, newdata = newX, newX = NULL, re.form = NULL,
                          variances=list(),
                          binding = FALSE, 
                          intervals = NULL,
                          level = 0.95,
                          ...) { ## but not new Y
  if ( ! is.null(variances$ranef)) {
    message("'variances$ranef' is obsolete: I interpret this as 'variances$linPred'.")
    variances$linPred <- variances$ranef
    variances$ranef <- NULL
  }
  if ( ! is.null(variances$sum)) {
    message("'variances$sum' is obsolete: I interpret this as 'variances$respVar'.")
    variances$respVar <- variances$sum
    variances$sum <- NULL
  }
  if (is.null(object$envir)) object$envir <- list2env(list(), ## back-compatibility fix for old objects
                                                     parent=environment(HLfit_body))
  ## the final components returned as attributes have names ...Var, other terms should be named differently
  checkIntervals <- (substr(x=intervals, nchar(intervals)-2, nchar(intervals))=="Var")
  if (any(!checkIntervals)) warning("Element(s)",intervals[!checkIntervals],"are suspect, not ending in 'Var'.")
  # possible elements in return value: fixefVar, predVar, residVar, respVar
  variances[intervals] <- TRUE 
  variances <- .process_variances(variances)
  new_X_ZACblob <- .calc_new_X_ZAC(object=object, newdata=newdata, re.form = re.form,
                 variances=variances, needNewNew=variances$cov)
  locdata <- new_X_ZACblob$locdata
  newX.pv <- new_X_ZACblob$newX.pv
  newZAlist <- new_X_ZACblob$newZAlist
  nrand <- length(newZAlist) ## may be reduced if non trivial re.form
  newZAClist <- new_X_ZACblob$newZAClist
  if (nrand>1L) {
    newZAC <- do.call(cbind,newZAClist)
  } else {newZAC <- newZAClist[[1]]}
  uuCnewold <- new_X_ZACblob$uuCnewold
  uuCnewnew <- new_X_ZACblob$uuCnewnew
  oldnewClist <- new_X_ZACblob$oldnewClist  # cf calcNewCorrs
  subZAlist <- new_X_ZACblob$subZAlist
  spatialOne <- new_X_ZACblob$spatialOne
  spatial.model <- new_X_ZACblob$spatial.model
  newinold <- new_X_ZACblob$newinold
  etaFix <- new_X_ZACblob$etaFix
    #
  ## (1) computes fv (2) compute predVar
  ##### fv
  if (.noRanef(re.form)) {
    fv <- object$family$linkinv(etaFix) 
  } else if ( is.null(newdata) && ! inherits(re.form,"formula")) {
    fv <- object$fv ## same length including replicates
    newZAlist <- subZAlist ## useful if predVar
  } else { ## 
    if ( nrand==0L ) {
      eta <- etaFix
      newZAlist <- NULL
    } else {
      ## for random-slope model the original eta's can be recomputed as 
      # object$X.pv %*% fixef(object) + 
      #    object$ZAlist[[1]] %*% object$strucList[[1]] %*% object$v_h
      #
      #### precomputation of coeffs
      ## on the gaussian scale, L.v_ori ~ lam C (lam C + phi I)^{-1}y 
      ## new random autocorr term ~ lam c (lam C + phi I)^{-1}y = c C^{-1} L_ori.v_ori = c [t(L_ori)]^{-1} v_ori
      ## [t(L_ori)]^{-1} v_ori can be computed once for all predictions => 'w_h_coeffs'
      if (is.null(w_h_coeffs <- object$envir$w_h_coeffs)) { 
        object$envir$w_h_coeffs <- w_h_coeffs <- .calc_invL_coeffs(object,object$v_h)
      }
      old_cum_n_u_h <- attr(object$lambda,"cum_n_u_h")
      lcrandfamfam <- attr(object$`rand.families`,"lcrandfamfam")[newinold]
      match_old_new_levels <- function(it) {
        oldu.range <- (old_cum_n_u_h[newinold[it]]+1L):(old_cum_n_u_h[newinold[it]+1L])
        if (it %in%  spatialOne) {     # %in% handles zero-length spatialOne...
          return(w_h_coeffs[oldu.range])          
        } else {
          oldlevels <- colnames(subZAlist[[it]])
          newlevels <- colnames(newZAClist[[it]])
          interlevels <- intersect(oldlevels,newlevels)
          oldpos <- which(oldlevels %in% interlevels) ## positions: handle replicates for random-coef
          newpos <- which(newlevels %in% interlevels)
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
          newv <- rep(vpsi_M,length(newlevels)) ## fills new levels with psi_M
          names(newv) <- newlevels
          newv[newpos] <- oldv[oldpos] 
          return(newv)
        }
      }
      augm_w_h_coeffs <- lapply(seq_len(nrand),match_old_new_levels)
      if (nrand>1L) {
        ZACw <- lapply(seq_len(nrand),function(it) {
          drop(newZAClist[[it]][] %*% augm_w_h_coeffs[[it]]) ## sapply preforms a cbind if everything is matrix (not Matrix)
        }) ## not sapply which binds 1*1 matrices into a vector of length nrand  
        ZACw <- do.call(cbind,ZACw)
        ZACw <- rowSums(ZACw)
      } else ZACw <- drop(newZAClist[[1]][] %*% augm_w_h_coeffs[[1]])
      eta <- etaFix + ZACw ## (length(eta)) col vector from coeffs = length(eta) row vector...
    }
    # done with eta
    fv <- object$family$linkinv(eta) ## ! freqs for binomial, counts for poisson
  }
  resu <- as.matrix(fv) ## suitable for objective function of optim() etc ## matrix ! maybe more suitable than data frame as objective function
  if ( ! is.logical(binding) ) { ## expecting a string
    binding <- makenewname(base=binding,varnames=colnames(locdata)) ## 09/11/2014 = posterior to CRAN 1.4.1 
    resu <- data.frame(resu)
    colnames(resu) <- binding
    resu <- cbind(locdata,resu) ## becomes a data frame !
    attr(resu,"fittedName") <- binding
  } else { ## expecting binding= FALSE
    if (ncol(locdata)>0)  attr(resu,"frame") <- locdata 
  }
  ##### (2) predVar
  if(variances$predVar && identical(attr(object$ZAlist, "anyRandomSlope"),TRUE)) {
    warning("This prediction variance computation is not yet implemented for random-coefficient models")
  }
  
  if(variances$linPred) {
    beta_cov <- get_beta_cov_any_version(object)
    beta_w_cov <- .get_beta_w_cov(object)
    if (inherits(re.form,"formula")) {
      # identifies an selects columns for the [retained ranefs, which are given by newinold 
      subrange <- unlist(lapply(newinold,function(it) {(old_cum_n_u_h[it]+1L):(old_cum_n_u_h[it+1L])}))
      Xncol <- ncol(beta_cov)
      subrange <- c(seq_len(Xncol),subrange + Xncol)
      beta_w_cov <- beta_w_cov[subrange,subrange]
    }
    if ( ! is.null(newdata)) {
      invColdoldList <- .get_invColdoldList(object)
      ## list for Cnewnew, which enters in  newZA %*% Cnewnew %*% tnewZA, hence should not represent newZA itself 
      if (nrand>0L) newnewClist <- lapply(seq_len(nrand),function(it) {
        if ( it %in% spatialOne) { ## avoids running the next algo which is slow on large matrices
          Cnn <- uuCnewnew ## already computed for point prediction
        } else { 
          Cnn <- Diagonal(ncol(newZAlist[[it]])) ## diag(ncol(newZAlist[[it]]))
        }
        return(Cnn)
      }) ## => completely ignores effects removed in re.form 
    } else {
      invColdoldList <- NULL
      newnewClist <- NULL
    }
    if (nrand>0L) {
      loclist <- list(X.pv=newX.pv,beta_w_cov=beta_w_cov,covMatrix=variances$cov,lambda=object$lambda)
      if (!is.null(uuCnewnew)) {
        uuCnewnewList <- lapply(seq_len(nrand),function(it) {
          if (it %in% spatialOne) { ## avoids running the next algo which is slow on large matrices
            return(uuCnewnew)
          } else return(NULL)
        }) ## should make sure it has the same structure as the other matrix lists
        loclist$CnewnewList <- uuCnewnewList 
      } ## else no loclist$CnewnewList (tested in calcPredVar) 
      #
      if ( ! is.null(newdata)) loclist$invColdoldList <- invColdoldList
      loclist$ColdnewList <- oldnewClist
      loclist$newZAlist <- newZAlist
      if (variances$disp) loclist$logdispObject <- .get_logdispObject(object)
      if (variances$cov) {
        respVar <- as.matrix(do.call(calcPredVar,loclist)) ## matrix, not Matrix (assumed below)
        rownames(respVar) <- colnames(respVar) <- rownames(locdata)
      } else {
        respVar <- do.call(calcPredVar,loclist) 
        names(respVar) <- rownames(locdata)
      }
    } else {
      if (variances$cov) {
        respVar <- matrix(0,nrow=nrow(locdata),ncol=nrow(locdata))
      } else respVar <- rep(0,nrow(locdata))
    }
  } else if (any(unlist(variances))) {
    respVar <- rep(0,nrow(locdata))
  } else respVar <- NULL 
  beta_cov <- get_beta_cov_any_version(object)
  if (! is.null(beta_cov)) {
    if ( variances$fixefVar || (nrand==0L && variances$linPred) ) {
      fixefcov <- newX.pv %*% beta_cov %*% t(newX.pv)
      if (variances$cov) {
        attr(resu,"fixefVar") <- fixefcov 
      } else attr(resu,"fixefVar") <- diag(fixefcov)
      if (nrand==0L) { ## otherwise there is already such a term in predVar
        respVar <- respVar + attr(resu,"fixefVar") 
      }
    }
  }
  attr(resu,"predVar") <- respVar ## vector or matrix
  if (variances$residVar) {
    pw <- object$prior.weights
    if ( ! (attr(pw,"unique") && pw[1L]==1L)) {
      if (! identical(spaMM.getOption("prior_weights_residvar_warned"),TRUE)) {
        warning("Prior weights are not taken in account in residVar computation.")
        spaMM.options(prior_weights_residvar_warned=TRUE)
      }
    }
    if (object$family$family %in% c("poisson","binomial","COMPoisson","negbin")) {
      attr(resu,"residVar") <- object$family$variance(fv)
    } else attr(resu,"residVar") <- calcResidVar(object,newdata=locdata) 
    if (inherits(respVar,"matrix")) {
      nc <- ncol(respVar)
      diagPos <- seq.int(1L,nc^2,nc+1L)
      respVar[diagPos] <- respVar[diagPos] + attr(resu,"residVar")
    } else respVar <- respVar + attr(resu,"residVar")
  }
  if (variances$respVar) attr(resu,"respVar") <- respVar
  if ( is.matrix(resu) && NCOL(resu)==1L) {
    class(resu) <- c("predictions",class(resu))
  } ## for print.predictions method which expects a 1-col matrix
  # intervals
  checkVar <- setdiff(intervals,names(attributes(resu)))
  if (length(checkVar)>0L) {
    warning(paste("Variances",paste(checkVar,collapse=", "),
                  "not available for interval computation.\n Check arguments."))
    intervals <- intersect(intervals,names(attributes(resu)))
  } 
  if(length(intervals)>0L) {
    intervalresu <- NULL
    for (st in intervals) {
      varcomp <- attr(resu,st)
      if (is.null(varcomp)) warning(paste("Prediction variance component",st,"requested but not available: check input."))
      if (is.matrix(varcomp)) varcomp <- diag(varcomp)
      eta <- object$family$linkfun(resu[,1L])
      pv <- 1-(1-level)/2
      ## special case for simple LM
      if (length(object$rand.families)==0L &&
          object$family$family=="gaussian" &&
          deparse(object$resid.predictor)=="~1" 
          ) { 
        resdf <- length(object$y) - ncol(object$X.pv) ## don't use fixef here, that contains bot NAs and etaFix$beta! 
        ## FR->FR (use a type attribute for fixef ?)
        sd <- qt(pv,df=resdf)*sqrt(varcomp)
      } else sd <- qnorm(pv)*sqrt(varcomp)
      interval <- cbind(object$family$linkinv(eta-sd),object$family$linkinv(eta+sd))
      colnames(interval) <- paste(st,c(signif(1-pv,4),signif(pv,4)),sep="_")
      intervalresu <- cbind(intervalresu,interval)
    }
    attr(resu,"intervals") <- intervalresu
  }
  return(resu)
}

print.vcov.HLfit <-function(x, expanded=FALSE, ...) {
  a <- attributes(x)
  attr(x,"beta_v_cov") <- NULL  
  print.default(x)
  cat("with additional attribute(s):")
  std.attr <- c("names","dim","dimnames","class") ## attributes not to be shown
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
  invisible() ## do not return x since it has lost a useful attribute
}


`[.predictions` <- function (x, i, j, 
                             drop = TRUE ## by default, this function will return scalar/vector  
                             ) {
  class(x) <- "matrix" ## removes "predictions" => set back later
  #   if (is.data.frame(x)) {
  #     resu <- x[i,j]
  #   } else 
  resu <- x[i,j,drop=drop]
  if ( ! drop) {
    fixefVar <- attr(x, "fixefVar")
    if ( ! is.null(fixefVar)) {
      if (is.null(dim(fixefVar))) {
        fixefVar <- fixefVar[x]
      } else fixefVar <- fixefVar[x,x,drop=FALSE]
    }
    predVar <- attr(x, "predVar")
    if ( ! is.null(predVar)) {
      if (is.null(dim(predVar))) {
        predVar <- predVar[x]
      } else predVar <- predVar[x,x,drop=FALSE]
    }
    frame <- attr(x, "frame")
    if ( ! is.null(frame)) frame <- frame[x,] ## dataframe => nodrop
    residVar <- attr(x, "residVar")
    if ( ! is.null(frame)) residVar <- residVar[x,drop=FALSE]
    respVar <- attr(x, "respVar")
    if ( ! is.null(respVar)) {
      if (is.null(dim(respVar))) {
        respVar <- respVar[x]
      } else respVar <- respVar[x,x,drop=FALSE]
    }
    class(resu) <- c("predictions","matrix")
    structure(resu,fixefVar=fixefVar,predVar=predVar,residVar=residVar,frame=frame,fittedName=attr(x, "fittedName"))
  } else return(resu)
} # Use unlist() to remove attributes from the return value

print.predictions <- function (x, expanded=FALSE, ...) {
  asvec <- as.vector(x) ## important to remove names and keep them separately
  rnames <- rownames(x)
  toolong <- nchar(rnames)>9
  rnames[toolong] <- paste(substr(rnames[toolong],0,8),".",sep="")
  names(asvec) <- rnames
  cat("Point predictions:\n")
  print(asvec)
  cat("*stored as* 1-col matrix with attributes:")
  std.attr <- c("names","dim","dimnames","class") ## attributes not to be shown
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