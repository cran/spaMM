locoptim <-
function(init.optim,LowUp,anyOptim.args,trace,optim.method,optim.args) {
    initvec <- unlist(init.optim)
    if (length(initvec)>1) {
      parscale <- rep(1,length(unlist(LowUp$lower))) ## unlist because some list elements may have length >1 (eg if rho has several components)
      ## optim uses par (vector) AND lower, upper
      anyOptim.args$lower <- unlist(LowUp$lower)
      anyOptim.args$upper <- unlist(LowUp$upper)
      ###### basic unrefined fit for initvec...
      anyHLCor.args <- anyOptim.args
      anyHLCor.args$ranefParsVec <- initvec ## logscale, inherits names from init.optim
      if (is.character(trace$file)) write("## Call for initial values",file=trace$file,append=T)   
      init.obj <- do.call(HLCor.obj,anyHLCor.args)
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
      corners.obj <- apply(corners,1,function(v) {
        anyHLCor.args$ranefParsVec <- v
        do.call(HLCor.obj,anyHLCor.args) ## this reconstructs a list of the form of initvec, then adds it to other anyHLCor$ranPars information
      })
      if (max(corners.obj) > init.obj) {
        initvec <- corners[which.max(corners.obj),] 
      }
      ## Search for good initial point done. 
      if (is.character(trace$file)) write("## Optimization call",file=trace$file,append=T)   
      if (optim.method=="nlminb") {
        ## nlminb code
        nlminbArgs <- anyHLCor.args ## not tested since long...
        nlminbArgs$skeleton <- init.optim
        nlminbArgs$objective <- function(...) {-HLCor.obj(...)}
        nlminbArgs$start <- initvec    
        nlminbArgs$lower <- unlist(LowUp$lower)
        nlminbArgs$upper <- unlist(LowUp$upper)
        nlminbArgs$scale <- 1/(anyOptim.args$upper - anyOptim.args$lower)
        optr <- do.call(nlminb,nlminbArgs) ## optimize 'HLCor.obj.value <- objective' (p_bv by default...?)
      } else {
        ndepsFac <- optim.args$ndepsFac
        if (is.null(ndepsFac)) ndepsFac <- 2000 ## was 10000 up to 060213
        ndeps <- (anyOptim.args$upper - anyOptim.args$lower)/ndepsFac
        control <- list(fnscale=-1,parscale=parscale,factr=1e9,ndeps=ndeps) ## default values... factr was the stricter 1e8 up to 23/01/13
        control[names(optim.args$control)] <- optim.args$control ## ...which may be overwritten 
        anyOptim.args <- c(anyOptim.args,list(par=initvec,control=control,method=optim.method))
        optr <- do.call(optim,anyOptim.args) ## optimize 'HLCor.obj.value <- objective' (p_bv by default...?)
      }
      optPars <- relist(optr$par,init.optim)
      ## HLCor.args$ranPars[names(optPars)] <- optPars 
    } else if (length(initvec)==1) { ## one-dimensional optimization
        byvar <- c(unlist(LowUp$lower),unlist(LowUp$upper)) ## HLCor expects trLambda...
        byvar <- 0.999 * byvar + 0.001 *mean(byvar)
        corners <- seq(byvar[1],byvar[2],length.out=11)
        anyHLCor.args <- anyOptim.args
        anyHLCor.args$ranefParsVec <- initvec ## logscale, inherits names from init.optim
        init.obj <- do.call(HLCor.obj,anyHLCor.args)
        if (is.character(trace$file)) write("## Calls for 'corners'",file=trace$file,append=T)   
        corners.obj <- sapply(corners,function(v) {
          anyHLCor.args$ranefParsVec <- v
          do.call(HLCor.obj,anyHLCor.args) ## this reconstructs a list of the form of initvec, then adds it to other anyHLCor$ranPars information
        })
        if (max(corners.obj) > init.obj) {
          initvec <- corners[which.max(corners.obj)] 
        }
        ## optimize uses only an interval, no initial value (!!)
        intervlower <- max(unlist(LowUp$lower),corners[corners<initvec])
        intervupper <- min(unlist(LowUp$upper),corners[corners>initvec])
        anyOptim.args$interval <-c(intervlower,intervupper)
        anyOptim.args <- c(anyOptim.args,list(maximum=T))
        if (is.character(trace$file)) write("## Optimization call",file=trace$file,append=T)   
        optr <- do.call(optimize,anyOptim.args)  ## optimize p_bv
        optPars <- relist(optr$maximum,init.optim)
        ## HLCor.args$ranPars[names(optPars)] <- optPars 
    } else { ## nothing to optimize
      ## HLCor.args$ranPars should already contain all it needs
        optPars <- NULL
    }
    return(optPars)
}
