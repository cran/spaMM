update.HL <-
function(object,...) {
  ## process object call
  gcall <- object$call
  if (is.null(gcall)) {
    mess <- pastefrom("No 'call' information in 'object'.",prefix="(!) From ")
    stop(mess)
  } else HL.args <- as.list(gcall)[-1]
  ## process dotlist; first the arguments ith sub-arguments...
  dotlist <- list(...)
  ranFix <- dotlist$ranFix
  if (!is.null(ranFix)) {
    HL.args$ranFix[names(ranFix)] <- ranFix
  }
  dotlist$ranFix <- NULL
  v_h <- dotlist$etaFix$v_h
  if (!is.null(v_h)) HL.args$etaFix$v_h <- v_h
  beta <- dotlist$etaFix$beta
  ## if (!is.null(beta)) HL.args$etaFix$beta[names(beta)] <- beta ## this does not work. Either beta and v_h are fully known => computation APHLs, or:
  #
  if (!is.null(beta)) {
    if (!is.null(dotlist$predictor) ) {
      mess <- pastefrom("Conflicting 'beta' and 'predictor' arguments.",prefix="(!) From ")
      stop(mess)
    }   
    predictor <- HL.args$predictor   
    if (class(predictor)=="call") predictor <- Predictor(formula=eval(predictor)) ## the class of elements in a match.call is not what one could expect...
    if (class(predictor)=="formula") predictor <- Predictor(formula=predictor)
    form <- predictor$formula
    if (is.null(form) || class(form)!="formula") {
      mess <- pastefrom("No formula information in 'predictor' argument.",prefix="(!) From ")
      stop(mess)
    } ## ELSE
    MeanFrames <- HLframes(formula=form,data=HL.args$data) ## design matrix X, Y, fixef names
    X.pv <- MeanFrames$X
    bars <- findbarsMM(form)  ## extract random effects
    lhs <- paste(form[[2]]) ## extract response
    ## build formula with only random effects
    formulaRE <- as.formula(paste(lhs,"~",paste(paste("(",as.character(bars),")"),collapse="+")))
    offset <- predictor$offset
    if (!is.null(offset)) {
      offset <- offset + X.pv %*% beta
    } else offset <- X.pv %*% beta
    predictor$formula <- formulaRE
    predictor$offset <- offset
    HL.args$predictor <- predictor
    ##
  }
  #
  dotlist$etaFix <- NULL
  ## then all other arguments
  HL.args[names(dotlist)] <- dotlist
  if (class(object)[1]=="HLCor") {
    update <- do.call(HLCor,HL.args)  
  } else update <- do.call(HLfit,HL.args)
  update
}
