corFactor.corMatern <-
function(object, ...)
{
  corD <- Dim(object)
  val <- .C("matern_factList",
	    as.double(as.vector(object)),
	    as.integer(attr(object, "nugget")),
            as.integer(attr(object, "nuScaled")),
	    as.double(unlist(getCovariate(object))),
	    as.integer(unlist(corD)),
	    as.double(attr(object, "minD")),
	    factor = double(corD[["sumLenSq"]]),
	    logDet = double(1))[c("factor", "logDet")]
  lD <- val[["logDet"]]
  val <- val[["factor"]]
  attr(val, "logDet") <- lD
  val
}
