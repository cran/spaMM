`coef<-.corMatern` <-
function(object, ..., value)
{
  if (length(value) != length(object)) {
    stop("cannot change the length of the parameter after initialization")
  }
#  cat("tutu",object[1],"->",value[1],object[2],"->",value[2],"\n")
  object[] <- value
  corD <- attr(object, "Dim")
  ## updating the factor list and logDet
  aux <- .C("matern_factList",
	    as.double(as.vector(object)),
	    as.integer(attr(object, "nugget")),
            as.integer(attr(object, "nuScaled")),
	    as.double(unlist(getCovariate(object))),
	    as.integer(unlist(corD)),
	    as.double(attr(object, "minD")),
	    factor = double(corD[["sumLenSq"]]),
	    logDet = double(1))[c("factor", "logDet")]
  attr(object, "factor") <- aux[["factor"]]
  attr(object, "logDet") <- -aux[["logDet"]]
  object
}
