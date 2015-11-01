iterateSEMSmooth <- function(processed, anyHLCor_obj_args, LowUp, init.corrHLfit, preprocess.formal.args, 
                            control.corrHLfit,verbose=interactive()) {
  pargrid <- sampleGridFromLowUp(LowUp,n=init.corrHLfit$nSmoothed) ## n may be NULL
  ## bits of codes needed whether PQL is run or not
  prevPredVars <- predVar <- 0
  lower <- LowUp$lower ## list ! which elements may have length >1 !
  upper <- LowUp$upper ## list !
  if (interactive()) {
    predi <- getProcessed(processed,"predictor")
    Xpv <- getProcessed(processed,"X.pv")
  }  
  #
  if ( ! getProcessed(processed,"SEMargs$SEMlogL") %in% c("GHK","pmvnorm")) { ## ie if var of logL estimates is large, so that their smoothing is difficult,
    ##                                                  use PQL to find a good starting region
    PQLarglist <- list(pargrid=pargrid,anyHLCor.args=anyHLCor_obj_args) ## copies anyHLCor_obj_args$`HLCor.obj.value` = objective
    locargs <- preprocess.formal.args
    locargs$HLmethod <- "PQL/L"
    PQLarglist$anyHLCor.args$processed <- do.call("preprocess",locargs)
    ###### 
    PQLoptr <- do.call("optimthroughSmooth",PQLarglist)       ############## CALL (with screen outputs)
    Krigobj <- PQLoptr$Krigobj
    predVar <- as.numeric(attr(predict(Krigobj,newdata=PQLoptr$par,variances=list(linPred=TRUE,dispVar=TRUE)),"predVar"))
    ## new sampling **************for SEM**************** guided by the PQL results
    blocksize <- 30 ## FR->FR comparer à mes settings pour Infusion...
    ## expand = 1 uses the fact that PQL is informative even if the smoothing must be redone.
    nextpoints <- sampleNextPoints(n=blocksize,Xpredy=PQLoptr$predictions,minPtNbr=3,expand=1, 
                                   D.resp=sqrt(predVar)/2) ## random sample
    info <- attr(nextpoints,"info") ## only used if diagnostic plot but removed by the rbind
    nearbypts <- sampleNearby(nextpoints,n=6,stepsizes=(unlist(upper)-unlist(lower))/100)      
    nextpoints <- rbind(nextpoints,nearbypts)
    ## diagnostic plot for previous and next computation
    if (interactive() && length(lower)==2L) {
      zut <- signif(unlist(Krigobj$corrPars$rho),4)
      titlemain <- bquote(paste(.(DEPARSE(predi))))
      if (nchar(eval(titlemain))>57) {
        titlemain <- bquote(paste(.(DEPARSE(nobarsNooffset(predi))),"+..."))
      }
      if (nchar(eval(titlemain))>57) {
        titlemain <- bquote(paste(.(substr(aschar,0,42)),"+... [length(",beta,")=",.(ncol(Xpv)),"]"))
      }
      if (length(zut)>1) {
        titlesub <- bquote(paste("PQL initialization: ",rho[f(rho)],"=",.(zut[1]),", ",rho[f(nu)],"=",.(zut[2]),"; predVar=",.(signif(predVar,4))))
      } else titlesub <- bquote(paste("PQL initialization: ",rho,"=",.(zut),"; predVar=",.(signif(predVar,4))))
      SEMdiagnosticPlot2D(Krigobj, smoothingOK=TRUE, titlemain=titlemain, titlesub=titlesub, nextpoints=nextpoints, info=info, optrPar=PQLoptr$par)
    }
    pargrid <- rbind(pargrid,nextpoints)
  }  ## end optional PQL block
  ## now the SEM computations
  allsmooths <- list(initSmooth=control.corrHLfit$initSmooth, ## NULL by default
                     resid.family=control.corrHLfit[["smooth.resid.family"]], ## NULL by default -> default controlled by optimthroughSmooth 
                     resid.formula=control.corrHLfit[["smooth.resid.formula"]] ## idem
  )  
  control.smooth <- allsmooths ## distinction between what goes in allsmooths and others is important ! nrepl will vary
  control.smooth$nrepl <- 20 ## number of points for which replicate estimates of likelihood are computed (modified later)
  processedHL1 <- getProcessed(processed,"HL[1]") ## there's also HLmethod in processed<[[]]>$callargs
  # if(processedHL1 =="SEM") 
  anyHLCor_obj_args$`HLCor.obj.value` <- "logLapp"
  arglist <- list(pargrid=pargrid,anyHLCor.args=anyHLCor_obj_args,control.smooth=control.smooth)
  it <- 1
  optr <- do.call("optimthroughSmooth",arglist)  ## first SEM iteration, "it=1"   ############## CALL (with screen outputs)
    ##
  smoothingOK <- FALSE
  EIfac <- control.corrHLfit$EIfac
  if (is.null(EIfac)) EIfac <- 1.96
  dit <- control.corrHLfit$dit ## NULL by default
  if (is.null(dit)) dit <- 0 ## default: controls test predVar < prevPredVars[it-dit] for smoothingOK or not
  continue <- TRUE
  nPredictors <- length(lower)
  while ( continue ) { ## note that some SEM results have already been analyzed previous to the loop
    control.smooth <- allsmooths ## reinitialize optimthroughSmooth arguments with constant ones 
    prevPredVars <- c(prevPredVars,predVar) ## has length it+1 ## predVar from previous Krigobj, not current Krigobj !
    Krigobj <- optr$Krigobj
    predVar <- as.numeric(attr(predict(Krigobj,newdata=optr$par,variances=list(linPred=TRUE,dispVar=TRUE)),"predVar"))
    if ( interactive() ) {
      if (it==1L) print(paste("residual variance formula:",
                              deparse(attr(Krigobj$resid.predictor,"oriFormula")), ## FR->FR maybe define safe extractor for predictor objects ?
                              ", with",Krigobj$resid.family$link,"link")) 
      cat(it," ");cat(paste(signif(optr$value,4),
                                                 "+/-",signif(sqrt(predVar),4),
                                                 "; n_points=",nrow(Krigobj$data),
                                                 "; smooth.lambda=",signif(Krigobj$lambda,4),sep=""))
    } 
    prevPtls <- optr$forSmooth
    smoothRho <- unlist(Krigobj$corrPars$rho)
    ## tests whether some correlation structure has been detected and adjust smoothing controls accordingly:
    # ... the best way it to perform some LRT on the smoothing parameters...
    smoothtest <- as.list(attr(Krigobj,"HLCorcall"))
    smoothrho <- smoothtest$ranPars$rho
    ## there is trRho or rho whether smoothing was performed or not ## FR->FR how to ensure info is in only one place ???  
    ## if both trRho and rho, trRho is used if HLCor call
    if (is.null(smoothrho)) {
      testedvalue <- rhoFn(rhoInv(smoothtest$ranPars$trRho)*2)
      if (is.nan(testedvalue)) { ## *2 exceeds max value => no real smoothing
        smoothtest <- FALSE
      } else {
        smoothtest <- eval(as.call(smoothtest))
        Krigobj$APHLs$p_bv> (smoothtest$APHLs$p_bv+1.92) ## test of information about rho_smooth 
      }
    } else {
      smoothtest$ranPars$rho <- smoothrho*2
      smoothtest <- eval(as.call(smoothtest))
      Krigobj$APHLs$p_bv> (smoothtest$APHLs$p_bv+1.92) ## test of information about rho_smooth
    } 
    tests <- c( it > dit, 
                predVar < prevPredVars[it-dit], ## ie for iter it-1-dit; for default dit=0, penultimate value 
                smoothtest
                )
    nextpoints <- sampleNextPoints(n=6,Xpredy=optr$predictions,expand=1,D.resp=sqrt(predVar)/2) ## always these 22/08/2014     
    if ( all(tests) ){ 
      smoothingOK <- TRUE 
      info <- attr(nextpoints,"info") ## only used if diagnostic plot but removed by the rbind
      nextpoints <- rbind(nextpoints,optr$par,optr$par) ## inferred maximum added in nextpoints ## 22/08/2014 
      control.smooth$nrepl <- 0 ## the above enforces duplicate of optr$par
      control.smooth$ranFix <- Krigobj$corrPars
    } else {
      smoothingOK <- FALSE
      if ( verbose) cat(" ",paste(c(paste("iter <=",dit),"high predvar","low LRT for scale parameters")[!tests], sep=" & "))
      nextpoints <- rbind(nextpoints,optr$par) ## inferred maximum added in nextpoints (not trypointsas EI might not retain it) ## 04/2015 ...
      trysize <- 180 ## a set from which blocksize points will be chosen for likelihood computation
      blocksize <- 9
      ## get rid of some possibly aberrant points that prevent good smoothing 
      prevPtls <- prevPtls[order(prevPtls$logLobj)[-c(1:2)],] ## FR->FR but aberrant points may not be the lowest... 
      trypoints <- sampleNextPoints(n=trysize,Xpredy=optr$predictions,expand=Inf,D.resp=sqrt(predVar)/2) ## random sample
      #info <- attr(trypoints,"info") ## might not be useful in this case (?)
      ## ... because otherwise algo may be fooled by incorrect maximum with high predVar, but missed by the EI procedure    
      ## => keeps inferrinf this incorrect maximum with high predVar over iterations.  
      ###### selection of points by improvement function with measurement error BinghamRW14 p. 121
      obspred <- predict(Krigobj,predVar=TRUE)
      obsSE <- attr(obspred,"predVar")
      obsSE[obsSE<0] <- 0
      obsSE <- sqrt(obsSE)
      Qmax <- max(obspred[,1]+EIfac * obsSE) ## best improvement function for already computed points 
      # 
      trypred <- predict(Krigobj,trypoints,predVar=TRUE)
      trySE <- attr(trypred,"predVar")
      trySE[trySE<0] <- 0
      trySE <- sqrt(trySE)
      tryQ <- trypred[,1] + EIfac*trySE ## improvement function for candidate points
      #
      expectedImprovement <- trySE*dnorm((Qmax-tryQ)/trySE)+(tryQ-Qmax)*pnorm((tryQ-Qmax)/trySE) ## 7.5 p. 121
      trypoints <- cbind(trypoints,EI=expectedImprovement)
      trypoints <- trypoints[order(trypoints[,"EI"],decreasing=TRUE)[seq_len(blocksize)],,drop=FALSE]
      trypoints <- trypoints[which(trypoints[,"EI"]>0),names(lower),drop=FALSE] ## maybe no point...
      nextpoints <- rbind(nextpoints,trypoints)
      ## need close pairs to estimate better the smoothing parameters
      ulower <- unlist(lower)
      uupper <- unlist(upper)
      epsilon <- (uupper-ulower)/1000
      ulower <- ulower+epsilon
      uupper <- uupper-epsilon ##useful for pmin, pmax 
      nearbypts <- sampleNearby(nextpoints,n=min(nrow(nextpoints),6),stepsizes=(uupper-ulower)/(100*smoothRho))     
      ## FR->FR problem: nearbypts may extrapolate... particularly for small smoothRho. We correct:
      for (ii in seq_len(length(ulower))) {
        nearbypts[,ii] <- pmax(nearbypts[,ii],ulower[ii])
        nearbypts[,ii] <- pmin(nearbypts[,ii],uupper[ii])
      }
      control.smooth$ranFix <- Krigobj$corrPars["nu"] ## passing original nu, fixed for this smoothing 
      control.smooth$nrepl <- ceiling(20/it - 0.0001)
      nextpoints <- rbind(nextpoints,nearbypts)
    }
    ## and a bit of extrapolation
    #       if (it>1) {
    #         cS <- connectedSets(info$simplicesTable)
    #         outerpoints <- lapply(cS, function(v){
    #           v <- intersect(v,info$innerVertexIndices) ## only the really good points in the set
    #           pts <- info$vertices[v,,drop=FALSE]
    #           if (nrow(pts)>length(lower)+1) { ## more vertices than a simplex => can be redundant
    #             return(pts[unique(as.vector(convhulln(info$vertices[v,],"Pp"))),])
    #           } else return(pts) ## extrapolhull will handle special cases
    #         })
    #         extrap <- lapply(outerpoints,extrapolhull)
    #         extrap <- do.call(rbind,extrap)
    #         nextpoints <- rbind(nextpoints,extrap)
    #       }
    ##
    if (interactive() ) {
      titlemain <- bquote(paste(.(DEPARSE(predi)),", iter=",.(it)))
      if (nchar(eval(titlemain))>50) {
        titlemain <- bquote(paste(.(DEPARSE(nobarsNooffset(predi))),"+..., iter=",.(it)))
      }
      if (nchar(eval(titlemain))>50) {
        titlemain <- bquote(paste(.(substr(aschar,0,35)),"+... [length(",beta,")=",.(ncol(Xpv)),"], iter=",.(it)))
      }
      if (nPredictors==2) {
        zut <- signif(smoothRho,4)
        if (length(zut)>1) {
          titlesub <- bquote(paste(rho[f(rho)],"=",.(zut[1]),", ",rho[f(nu)],"=",.(zut[2]),
                                   "; max=",.(signif(optr$value,4)),"; predVar=",.(signif(predVar,4))))
        } else titlesub <- bquote(paste(rho[smooth],"=",.(zut),"; max=",.(signif(optr$value,4)),"; predVar=",.(signif(predVar,4))))
        SEMdiagnosticPlot2D(Krigobj, smoothingOK=smoothingOK, titlemain=titlemain, titlesub=titlesub, nextpoints=nextpoints, info=info, optrPar=optr$par)
        SEMdiagnosticPlot(Krigobj,"Raw profiles",optr) ## as the title says
        ## it woudl be nice to have contour lines on a SEMdiagnosticPlot2D -> spaMMplot2D but this requires a grid of values+ smoothing as in spaMM.filled.contour   
      } else { SEMdiagnosticPlot(Krigobj, titlemain=titlemain, optr)      }
    }
    #browser()
    arglist <- list(pargrid=nextpoints,control.smooth=control.smooth,anyHLCor.args=anyHLCor_obj_args,prevPtls=prevPtls)
    it <- it+1
    optr <- do.call("optimthroughSmooth",arglist) ## it>1    ############## CALL (with screen outputs)
    ## terminates if either of these two considtions are reached *...* :
    if (predVar < 0.02) continue <- FALSE
    #if (nrow(Krigobj$data) > 1000) continue <- FALSE ## FR->FR can only be tempo
    ## ... UNLESS one of these conditions are true
    #if (it < 10) continue <- TRUE
    if (predVar > prevPredVars[it]) continue <- TRUE
  } ## end 'while' loop
  Krigobj <- optr$Krigobj
  predVar <- as.numeric(attr(predict(Krigobj,newdata=optr$par,variances=list(linPred=TRUE,dispVar=TRUE)),"predVar"))
  if (interactive() && length(lower)==2) {
    zut <- signif(unlist(Krigobj$corrPars$rho),4)
    if (length(zut)>1) {
      titlemain <- bquote(paste(rho[f(rho)],"=",.(zut[1]),", ",rho[f(nu)],"=",.(zut[2]),
                                "; max=",.(signif(optr$value,4)),"; predVar=",.(signif(predVar,4))))
    } else titlemain <- bquote(paste(rho[smooth],"=",.(zut),"; max=",.(signif(optr$value,4)),"; predVar=",.(signif(predVar,4))))
    if (nPredictors==2) {
      SEMdiagnosticPlot2D(Krigobj, smoothingOK=FALSE, titlemain=titlemain, titlesub="", nextpoints=NULL, info=NULL, optrPar=optr$par)
    } else { SEMdiagnosticPlot(Krigobj, titlemain=titlemain, optr)      }
  }
  attr(optr$value,"predVar") <- predVar
  return(optr)
}