dispGammaGLM <-
function(dev.res,lev,X,offset=NULL,family=GammaForDispGammaGLM(link=log),method="glm") {
  ## do not 'filter' the dev.res and lev (in any way different from the lamScal one) here otherwise the lambda estimate may be inconsistent with the v_h
  ## except that fatally, a 0L response value occurs in normal use; eg Poisson, y ~ 1+ ranef(lambda->very small), mean(y) happens to be equal to some y value -> u_h =0L
  resp <- dev.res/((1-lev)) 
  resp[resp==0L] <- 1e-100
  resp[resp>1e150] <- 1e150 ## v1.2 fix for glm -> glm.fit -> .Call(C_Cdqrls, x[good, , drop = FALSE]*w problem
  weight <- (1-lev)/2 
  etastart <- rep(family$linkfun(mean(resp)),nrow(resp))   ## glm needs a bit help...
  if (is.null(offset)) offset <- rep.int(0, nrow(resp))
  etastart <- etastart - offset
  if (method=="glm"){
    if (interactive()) {
      resglm <- glm(resp~X-1,family=family,weights = weight,etastart=etastart,,offset=offset)
    } else {
      resglm_wrap <- tryCatch.W.E(glm(resp~X-1,family=family,weights = weight,etastart=etastart,offset=offset))    
      resglm <- resglm_wrap$value
      attr(resglm,"warning") <- resglm_wrap$warning$message ## may be NULL
    }
    # if (verbose["trace"]) {print(paste("phi coefficients=",paste(signif(resglm$coefficients,4),collapse=", ")),quote=F)}
    Qr <- resglm$qr  
    beta_disp <- resglm$coefficients[Qr$pivot[1L:ncol(X)]] ## As in summary.glm.
  } else {
    ## this will remain similar to glm aslong as HLM -> provide.resglm -> glm.fit 
    stop("need to put back HLM.R into the sources")
    #resglm <- HLM(resp~X-1,family=family,prior.weights = weight,offset=offset) 
    beta_disp <- resglm$fixef
    summ <- NULL
  }
  return(list(beta_disp=beta_disp,next_disp_est=fitted(resglm),resglm=resglm))
}
