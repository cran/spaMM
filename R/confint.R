confint.HLfit <- function(object,parm,level=0.95,verbose=TRUE,...) {
  dlogL <- qchisq(level,df=1)/2
  znorm <- qnorm((1+level)/2)
  if (is.character(parm)) {
    whichcol <- which(names(object$fixef)==parm)
    if (length(whichcol)==0L) stop("Parameter not in the model")
    attr(parm,"col") <- whichcol 
  } else {
    parmcol <- parm
    if (parm > length(object$fixef)) stop("'parm' not compatible with # of fixed-effects coefficients")
    parm <- names(object$fixef)[parmcol]
    attr(parm,"col") <- parmcol
  }
  llc <- as.list(getCallHL(object))
  if (llc[[1]]=="corrHLfit") {
    lc <- as.list(attr(object,"HLCorcall"))
  } else lc <- llc
  lc$control.HLfit$intervalInfo$fitp_v <- object$APHLs$p_v
  lc$control.HLfit$intervalInfo$targetp_v <- object$APHLs$p_v-dlogL
  lc$control.HLfit$intervalInfo$MLparm <- object$fixef[parm]
  lc$control.HLfit$intervalInfo$parm <- parm
  lc$control.HLfit$LevenbergM <- FALSE ## tempo mais simple
  beta_se <- sqrt(diag(object$beta_cov))
  lc$control.HLfit$intervalInfo$init <- (object$fixef-znorm* beta_se)[parm]
  if (llc[[1]]=="corrHLfit") {
    olc <- lc
    ## good starting values are important... important to use canonizeRanPars as in HLCor
    trTemplate <- attr(object,"optimInfo")$`optim.pars`
    objfn <- function(ranefParsVec) { 
      ## FR->FR bc currently locoptim expects a fn with first arg ranefParsVec
      ## ca serait mieux de pas avoir de contrainte la dessus et de pvr nommer l'arg trParsVec
      ## bc HLCor call uses transformed scale for ranPars
      olc$ranPars[names(trTemplate)] <- ranefParsVec 
      locfit <- eval(as.call(olc)) ## HLCor call with optimized corrpars vars and fixed ones kept
      resu <- locfit$fixef[parm]
      attr(resu,"info") <- locfit$APHLs$p_v 
      ## attribute lost by optim but otherwise useful for debugging 
      return(resu)
    }
    canonTemplate <- canonizeRanPars(ranPars=trTemplate,
                                 corr.model=lc$`corr.model`,
                                 checkComplete=FALSE
    )$ranPars
    LUarglist <- list(canon.init=canonTemplate,
                      lower=attr(object,"optimInfo")$`init.optim`,
                      upper=attr(object,"optimInfo")$`init.optim`,
                      user.lower=attr(object,"optimInfo")$`user.lower`,
                      user.upper=attr(object,"optimInfo")$`user.upper`, 
                      corr.model=lc$`corr.model`,nbUnique=attr(lc$distMatrix,"Size"),
                      ranFix=llc$ranFix,
                      optim.scale="transformed") ## FR->FR transformed is a guess
    LowUp <- do.call("makeLowerUpper",LUarglist)
    loclist <- list(init.optim=trTemplate,LowUp=LowUp,objfn=objfn,anyObjfnCall.args=list(),corners=TRUE,optimizers.args=list())
    optr <- do.call(locoptim,loclist) ## minimizes the confint parameter
    ## recover HLCor fit for optimized params
    olc$ranPars[names(trTemplate)] <- optr[names(trTemplate)] ## should be optr, but order ?
    lowerfit <- eval(as.call(olc))
    ##
  } else lowerfit <- eval(as.call(lc))
  ## upperfit:
  lc$control.HLfit$intervalInfo$init <- (object$fixef+znorm* beta_se)[parm]
  if (llc[[1]]=="corrHLfit") {
    olc <- lc
    loclist$maximize <- TRUE ## 
    optr <- do.call(locoptim,loclist) ## maximizes the confint parameter
    ## recover HLCor fit for optimized params
    olc$ranPars[names(trTemplate)] <- optr[names(trTemplate)]
    upperfit <- eval(as.call(olc))
    ##
  } else upperfit <- eval(as.call(lc))
  interval <- c(lowerfit$fixef[parm],upperfit$fixef[parm])
  names(interval) <- paste(c("lower","upper"),parm)
  if (verbose) print(interval)
  invisible(list(lowerfit=lowerfit,upperfit=upperfit,interval=interval))
}
