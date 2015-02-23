sampleNearby <- function(focalPts,n=NULL,stepsizes) {
  d <- ncol(focalPts)
  nr <- nrow(focalPts)
  if (n<nr) {
    subfocal <- focalPts[sample(seq_len(nr),n),,drop=FALSE]
  } else subfocal <- focalPts
  randpts <- apply(subfocal,1,function(v) {
    rdxy <- runif(d)* sample(c(-1,1),d,replace=TRUE)* stepsizes/2
    v + rdxy 
  })
  if (d==1L) {
    randpts <- as.matrix(randpts) 
  } else randpts <- t(randpts)
  colnames(randpts) <- colnames(focalPts) ## (in 1D at least) rbind checks that the names match...
  randpts
}

## uses gridSteps or n depending on sampling
## default sampling depends on length(lower)
sampleGridFromLowUp <- function(LowUp,n=NULL,gridSteps=NULL,sampling=NULL) {
  ## grid from lower, upper
  lower <- LowUp$lower
  upper <- LowUp$upper
  d <- length(lower)
  if (is.null(sampling)) {
    if (d<4) {sampling="grid"} else {sampling="rvolTriangulation"}
  }
  byvar <- t(rbind(unlist(lower),unlist(upper))) 
  byvar <- 0.999 * byvar + 0.001 *rowMeans(byvar)
  grillelist <- list()
  if (is.null(gridSteps)) {
    gridSteps <- c(20,6,4,3)[min(d,4)] ## => 20  36  64  81 243 729 points 
  }
  for(name in rownames(byvar)) {grillelist[[name]] <- seq(byvar[name,1],byvar[name,2],length.out=gridSteps)}
  pargrid <- expand.grid(grillelist)
  if (sampling=="rvolTriangulation") { ## while randomly sample the simplices defined by the regular grid
    vT <- volTriangulation(pargrid) ## note no convhulln call -> many internal simplices
    if (is.null(n)) n <- gridSteps^min(d,6) ## default maxi 729 for high d
    pargrid <- rbind(pargrid,rvolTriangulation(n,vT)) ## regular + random
  } else { ## more structured sampling: will sample the grid 
    attr(pargrid,"regularGrid") <- seq_len(nrow(pargrid))
    ## redefines a grid
    insides <- lapply(grillelist, function(v) {(v[2]-v[1])/2+v[-c(gridSteps)]}) ## 19 5 3 2 2 2...
    stepsizes <- unlist(lapply(insides, function(v) {v[2]-v[1]})) 
    insides <- expand.grid(insides) ## 19 25 27 16 32 64 ...
    randgrid <- sampleNearby(insides,n=nrow(insides),stepsizes=stepsizes) ## not really nearby given the large stepsizes
    pargrid <- rbind(pargrid,randgrid) ## regular + random
  }
  ##
  pargrid
}

#pargrid <- gridFromLowUp(LowUp)

# doc provisoire:
# pargrid : a matrix of fittedPars
# anyHLCor.args = arguments for evaluation of likelihood by HLCor
# prevPtls : previous points (fittedPars, respName, lambda)
# control.smooth : smoothing parameters (typically $rho )and/or number of duplicates ($nrepl)
#
# This function estimates likelihod in the input points + a second estimate in nrepl of these points 
# The computes a smoothed likelihood surface (default $nu=4) by ordinary kriging
# Then finds the maximum of this smoothed likelihood surface
#
# retruns a list with elements
# $par: cf optimize
# $value: cf optimize
# $predictions: predicted responses values at smoothed point values; attributes include an estimate of the variance of [likelihood estimation by HLCor]
# $Krigobj: ... 
# $forSmooth: the input data of the smoothing computation, i.e. the estimates of likelihood by HLCor
#
#

locoptimthroughSmooth <- function(pargrid,anyHLCor.args,prevPtls=NULL,control.smooth=list()) {
  anyHLCor.args$processed <- setProcessed(anyHLCor.args$processed,"SEMseed",value="NULL") ## SEMseed bien pour controler individuellement un SEM mais pas une série 
  ranges <- apply(pargrid,2,range)
  LowUp <- apply(ranges,1,as.list)
  names(LowUp) <- c("lower","upper")
  lower <- LowUp$lower
  upper <- LowUp$upper
  ## 
  processedHL1 <- getProcessed(anyHLCor.args$processed,"HL[1]") ## there's also HLmethod in processed<[[]]>$callargs
  logLobj <- anyHLCor.args$`HLCor.obj.value`
  grid.obj <- apply(pargrid,1,function(v) {
    anyHLCor.args$ranPars[names(lower)] <- relist(v,lower)
    hlcor <- do.call("HLCor",anyHLCor.args) ## this reconstructs a list of the form of initvec, then adds it to other anyHLCor$ranPars information
    c("logLobj"=hlcor$APHLs[[logLobj]],lambda=hlcor$lambda,
      seInt =attr(hlcor$APHLs[[logLobj]],"seInt")) ## seInt attr may be NULL then no seInt element
  })
  ## apply does not name things as oen would wish, and moreover has/had inconsistent behaviouris some R-devel version
  rownames(grid.obj) <- c("logLobj","lambda","seInt")[seq_len(nrow(grid.obj))]
  nrepl <- control.smooth$nrepl
  processedHL1 <- getProcessed(anyHLCor.args$processed,"HL[1]") ## there's also HLmethod in processed<[[]]>$callargs
  if (processedHL1 == "SEM") {
    if (is.null(nrepl)) {
      nrepl <- Inf ## default is to duplicate all simuls !
    }
    if (nrepl == 0 && (length(control.smooth$ranFix) != length(lower))) {
      message("(!) From locoptimthroughSmooth: suspect nrepl==0 for stochastic simulation with correlation parameters to be estimated.")
    }
  } else {
    if (is.null(nrepl)) {
      nrepl <- 0
    }
    if (nrepl>0) message("(!) From locoptimthroughSmooth: suspect nrepl>0 for deterministic simulation.")
  }  
  if (nrepl>0) {
    subpargrid <- pargrid[sample(nrow(pargrid),min(nrepl,nrow(pargrid))),,drop=FALSE]
    subgrid.obj <- apply(subpargrid,1,function(v) {
      anyHLCor.args$ranPars[names(lower)] <- relist(v,lower)
      hlcor <- do.call("HLCor",anyHLCor.args) ## this reconstructs a list of the form of initvec, then adds it to other anyHLCor$ranPars information
      c("logLobj"=hlcor$APHLs[[logLobj]],lambda=hlcor$lambda,
        seInt =attr(hlcor$APHLs[[logLobj]],"seInt")) ## seInt attr may be NULL then no seInt element
    })
    rownames(subgrid.obj) <- c("logLobj","lambda","seInt")[seq_len(nrow(grid.obj))]
    forSmooth <- cbind(rbind(pargrid,subpargrid),rbind(t(grid.obj),t(subgrid.obj)))      
  } else {
    forSmooth <- cbind(pargrid,t(grid.obj))      
  }
  forSmooth <- data.frame(forSmooth)
  if ( ! is.null(prevPtls)) forSmooth <- rbind(prevPtls,forSmooth) ## older points first!
  RHS <- paste(names(lower),collapse="+")
  form <- as.formula(paste("logLobj ~ 1+Matern(1|",RHS,")"))
  ## ****** predict a (profile) likelihood surface for the correlation parameters ******
  ## phi here represent varSEMandMCint the variance of MC integration and that of the SEM algo
  ##  phi = a [SEM]+ b lambda [MCint]: 
  ## i e log(phi) = log(a+ b lambda), not a +b log(lambda) => code pas ideal
  ## fortunately lambda may be roughly indep of the corrpars
  ## Further the IRLS in glm.fit is fussy... 
  ## et le fit des dev resids peut etre surprenant, cf plot(var explicative,log(resp)) dans dispGammaGLM -> glm(resp ~...)
  #Krigobj <- corrHLfit(form,data=forSmooth,resid.formula= ~log(lambda),init.corrHLfit=list(rho=rep(1,length(lower))),ranFix=list(nu=4))
  ranFix <- list(nu=4)
  processedHL1 <- getProcessed(anyHLCor.args$processed,"HL[1]") ## there's also HLmethod in processed<[[]]>$callargs
  if (processedHL1 != "SEM") ranFix$phi <- 1e-06 ## interpolation
  ranFix[names(control.smooth$ranFix)] <- control.smooth$ranFix
  initSmooth <- control.smooth$initSmooth ## corr pars du smoothing
  if (is.null(initSmooth)) initSmooth <- rep(1,length(lower))
  init.corrHLfit <- list(rho=initSmooth) ## important as it gives the length of rho to corrHLfit
  init.corrHLfit[names(ranFix)] <- NULL
  if (is.null(forSmooth$seInt)) {
    resid.formula <- ~1
  } else {
    resid.formula <- control.smooth[["resid.formula"]]
    if(is.null(resid.formula)) resid.formula <- ~1+offset(seInt^2) ## default... ######   ~log(seInt)+I(log(seInt)^2)
  }  
  resid.family <- control.smooth[["resid.family"]]
  if(is.null(resid.family)) resid.family <- Gamma(identity) ## default
  Krigobj <- corrHLfit(form,data=forSmooth,resid.formula= resid.formula ,init.corrHLfit=init.corrHLfit,ranFix=ranFix,
                       control.HLfit=list(resid.family=resid.family))
  ### quick check of variance of logL estimation:
  if (FALSE) {
    essai <- forSmooth[do.call(order,forSmooth[,seq_len(length(lower))]),,drop=FALSE]
    diffs <- apply(essai,2,diff)
    ## next line selects diff lines for identical coordinates and takes from them the diff value for logL 
    logLdiffs <- diffs[apply(diffs[,seq_len(length(lower))]==rep(0,length(lower)),1,all),length(lower)+1]
    var(logLdiffs)/2 ## since diff= sum(eps_1+eps_2)
  }
  ###
  ## ****** optimize in the predicted likelihood surface ******
  predictions <- predict(Krigobj,binding="fitted")
  predictions <- predictions[order(predictions[,1],decreasing=TRUE),]
  initvec <- predictions[1,names(lower)]
  ## redefines lower, upper, for maximization
  ranges <- apply(forSmooth[,names(lower),drop=FALSE],2,range)
  LowUp <- apply(ranges,1,as.list)
  names(LowUp) <- c("lower","upper")
  lower <- LowUp$lower
  upper <- LowUp$upper
  ##
  if (length(initvec)==1L) {
    optr <- optimize(function(v) {predict(Krigobj,v)[,1]},maximum=TRUE,lower=unlist(lower),upper=unlist(upper)) 
    optr$par <- optr$maximum
    optr$value <- optr$objective
    optr$maximum <- NULL
    optr$objective <- NULL
  } else { 
    optr <- optim(initvec,function(v) {predict(Krigobj,v)[,1]},method="L-BFGS-B",
                  control=list(fnscale=-1),lower=unlist(lower),upper=unlist(upper))
  }
  names(optr$par) <- names(lower)
  attr(predictions,"fittedPars") <- names(lower) 
  ####### #attr(predictions,"respName") <- "fitted" ## predict() has provided attribute "fittedName" (a bit confusing)
  attr(predictions,"MSy") <- Krigobj$phi ## FR->FR ecrire un extracteur pour phi... 
  optr$predictions <- predictions 
  #  ranges <- apply(uppersurf,2,range)
  #  nextLowUp <- apply(ranges,1,as.list)
  #  names(nextLowUp) <- c("lower","upper")
  ######
  #  optr$nextLowUp <- nextLowUp
  optr$Krigobj <- Krigobj ## 
  optr$forSmooth <- forSmooth
  return(optr)
} ## end def locoptimthroughSmooth
