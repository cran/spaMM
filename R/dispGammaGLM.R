dispGammaGLM <-
function(dev.res,lev,X,offset=NULL) {
  ## do not 'filter' the dev.res and lev (in any way different from the lamScal one) here otherwise the lambda estimate may be inconsistent with the v_h
  ## except that fatally, a 0L response value occurs in normal use; eg Poisson, y ~ 1+ ranef(lambda->very small), mean(y) happens to be equal to some y value -> u_h =0L
  resp <- dev.res/(1-lev)
  resp[resp==0L] <- 1e-100
  weight <- (1-lev)/2
  etastart <- rep(log(mean(resp)),nrow(resp))   ## glm needs a bit help...
  if (is.null(offset)) offset <- rep.int(0, nrow(resp))
  etastart <- etastart - offset
  resglm<-glm(resp~X-1,family=GammaForDispGammaGLM(link=log),weights = weight,etastart=etastart)
  return(resglm)
}
