old.nextPoints <-
function(n=1,optr,replace=TRUE) { ## random sampling of volume defined from previous fit
  uP <- upperPoints(optr$predictions) ## indices
  uP <- optr$predictions[uP,attr(optr$predictions,"fittedPars")]
  uP <- rbind(uP,optr$par) ## not sure this is useful for volumetric sampling
  erV <- elim.redundant.V(uP)  
  vT <- volTriangulation(erV)
  rvolTriangulation(n,vT,replace=replace)
}
