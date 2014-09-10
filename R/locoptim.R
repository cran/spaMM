locoptim <-
function(init.optim,LowUp,anyObjfnCall.args,trace=list(file=NULL,append=T),Optimizer="L-BFGS-B",optimizers.args,corners=TRUE,objfn=HLCor.obj,maximize=FALSE) {
  processedHL1 <- getProcessed(anyObjfnCall.args$processed,"HL[1]") ## there's also HLmethod in processed<[[]]>$callargs
  initvec <- unlist(init.optim)
  if ( ! is.null(anyObjfnCall.args$ranefParsVec) ) stop("! is.null(anyObjfnCall.args$ranefParsVec)") ## FR->FR tempo debugging
  ## anyObjfnCall.args$ranefParsVec <- NULL ## removes this for optimization! otherwise fatal for nlminb! ## but should become obsolete ?
  optimizerObjfnCall.args <- notoptimObjfnCall.args <- anyObjfnCall.args
  notoptimObjfnCall.args$ranefParsVec <- initvec ## logscale, inherits names from init.optim
  if (length(initvec)>1) {
    ## old version before 2014/09/03
    #parscale <- rep(1,length(unlist(LowUp$lower))) ## unlist because some list elements may have length >1 (eg if rho has several components)
    ## new version 2014/09/03
    parscale <- (unlist(LowUp$upper)-unlist(LowUp$lower)) ## unlist because some list elements may have length >1 (eg if rho has several     ###### basic unrefined fit for initvec...
    if (is.character(trace$file)) write("## Call for initial values",file=trace$file,append=T)   
    if (corners) {      
      init.obj <- do.call(objfn,notoptimObjfnCall.args)
      ###### look in the corners 
      ## L-BFGS-B tends to find the local max closest to the initial point. Search for good initial point:
      byvar <- t(rbind(unlist(LowUp$lower),unlist(LowUp$upper))) ## HLCor expects trLambda...
      ### only the optim target variables are retained from the REML fit, although the lambda valuefrom this fit might be the most useful.  
      ## But if corrHLfit was called with an ranFix, we do want this (as part of ranPars) in the HLCor.obj calls. 
      byvar <- 0.999 * byvar + 0.001 *rowMeans(byvar)
      grillelist <- list()
      gridSteps <- floor(35^(1.05/length(initvec))) ## 6 for 2 pars,  3 for 3 pars, then 2  2  1  1 
      gridSteps <- max(2,gridSteps)
      for(name in rownames(byvar)) {grillelist[[name]] <- seq(byvar[name,1],byvar[name,2],length.out=gridSteps)}
      corners <- expand.grid(grillelist)
      ## ranefParsVec corresponds to ranPars but as vector, not as list
      ## uses HLCor.obj because of the return value...
      if (is.character(trace$file)) write("## Calls for 'corners'",file=trace$file,append=T)  
      if (FALSE) { ## version for debug can use info attribute
        scorners <- split(corners,rownames(corners))
        ll <- lapply(scorners,function(v) {
          notoptimObjfnCall.args$ranefParsVec <- v
          ## FR->FR serait mieux de pouvoir supprimer que certain warnings
          bla <- suppressWarnings(do.call(objfn,notoptimObjfnCall.args)) 
          c(bla,attr(bla,"info"))
        })
        corners.objWithAttr <- t(sapply(ll,function(v) {c(v,info=attr(v,"info"))}))
        corners.obj <- corners.obj[,1]
      } else {
        corners.obj <- apply(corners,1,function(v) {
          notoptimObjfnCall.args$ranefParsVec <- v
          ## FR->FR serait mieux de pouvoir supprimer que certain warnings
          suppressWarnings(do.call(objfn,notoptimObjfnCall.args)) 
        })
      }
      #      browser()
      if (maximize) {
        if (max(corners.obj) > init.obj) {
          initvec <- corners[which.max(corners.obj),] 
        }
      } else {
        if (min(corners.obj) < init.obj) {
          initvec <- corners[which.min(corners.obj),] 
        }
      }
    }
    ## Search for good initial point done. 
    if (is.character(trace$file)) write("## Optimization call",file=trace$file,append=T)   
    if (Optimizer=="nlminb") {
      ## nlminb code
      nlminbArgs <- optimizerObjfnCall.args 
      nlminbArgs$skeleton <- init.optim
      if (maximize) {
        nlminbArgs$objective <- function(...) {-objfn(...)}
      } else nlminbArgs$objective <- function(...) {objfn(...)}
      nlminbArgs$start <- initvec    
      nlminbArgs$lower <- unlist(LowUp$lower)
      nlminbArgs$upper <- unlist(LowUp$upper)
      nlminbArgs$scale <- 1/(nlminbArgs$upper - nlminbArgs$lower)
      nlminbArgs$control <- optimizers.args$nlminb$control
      optr <- do.call(nlminb,nlminbArgs) ## optimize HLCor.obj.value (=objective=p_bv by default, except for SEM)
    } else {
      ndepsFac <- optimizers.args$optim$ndepsFac
      ## if (is.null(ndepsFac)) ndepsFac <- 2000 ## was 10000 up to 060213, then 2000 up to 04092014; then parscale was modified
      ## length(unlist(LowUp$upper)) != length(LowUp$upper) because list elements may have length > 1
      ndeps <- rep(1e-03,length(unlist(LowUp$upper))) # this is actually optim's default ## previously was ( unlist(LowUp$upper) -  unlist(LowUp$lower))/ndepsFac until parscale was modified
      control <- list(parscale=parscale,factr=1e9,ndeps=ndeps) ## default values... factr was the stricter 1e8 up to 23/01/13
      if (maximize) control$fnscale <- -1     
      control[names(optimizers.args$optim$control)] <- optimizers.args$optim$control ## ...which may be overwritten 
      anyBaseOptim.args <- c(optimizerObjfnCall.args,list(par=initvec,control=control,method=Optimizer))
      optr <- do.call(optim,c(anyBaseOptim.args,list(fn=objfn))) ## optimize HLCor.obj.value (=p_bv by default, except for SEM)
    }
    optPars <- relist(optr$par,init.optim)
    ## HLCor.args$ranPars[names(optPars)] <- optPars 
  } else if (length(initvec)==1) { ## one-dimensional optimization
    byvar <- c(unlist(LowUp$lower),unlist(LowUp$upper)) ## HLCor expects trLambda...
    byvar <- 0.999 * byvar + 0.001 *mean(byvar)
    corners <- seq(byvar[1],byvar[2],length.out=11)
    init.obj <- do.call(objfn,notoptimObjfnCall.args)
    if (is.character(trace$file)) write("## Calls for 'corners'",file=trace$file,append=T)   
    corners.obj <- sapply(corners,function(v) {
      notoptimObjfnCall.args$ranefParsVec <- v
      do.call(objfn,notoptimObjfnCall.args) ## this reconstructs a list of the form of initvec, then adds it to other anyHLCor$ranPars information
    })
    if (maximize) {
      if (max(corners.obj) > init.obj) {
        initvec <- corners[which.max(corners.obj)] 
      }
    } else {
      if (min(corners.obj) < init.obj) {
        initvec <- corners[which.min(corners.obj)] 
      }
    }
    ## optimize uses only an interval, no initial value (!!)
    intervlower <- max(unlist(LowUp$lower),corners[corners<initvec])
    intervupper <- min(unlist(LowUp$upper),corners[corners>initvec])
    anyOptimize.args <- optimizerObjfnCall.args
    if (maximize) anyOptimize.args$maximum <- TRUE
    anyOptimize.args$interval <-c(intervlower,intervupper)
    anyOptimize.args$tol <- optimizers.args$optimize$tol
    if (is.character(trace$file)) write("## Optimization call",file=trace$file,append=T)   
    optr <- do.call(optimize,c(anyOptimize.args,list(f=objfn)))  ## optimize p_bv
    optPars <- relist(optr$maximum,init.optim)
    ## HLCor.args$ranPars[names(optPars)] <- optPars 
  } else { ## nothing to optimize
    ## HLCor.args$ranPars should already contain all it needs
    optPars <- NULL
  }
  return(optPars)
}
