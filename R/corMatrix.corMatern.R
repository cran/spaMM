corMatrix.corMatern <-
function(object, covariate = getCovariate(object), corr = TRUE, ...)
{

  if (data.class(covariate) == "list") { ## groups are defined
    if (is.null(names(covariate))) {
      names(covariate) <- 1:length(covariate)
    }
    corD <- Dim(object, rep(names(covariate),
			    unlist(lapply(covariate,
                                          function(el) round((1 + sqrt(1 + 8 * length(el)))/2)))))
  } else { ## no groups
    corD <- Dim(object, rep(1, round((1 + sqrt(1 + 8* length(covariate)))/2)))
  }

  if (corr) {
    val <- .C("matern_matList",
	      as.double(as.vector(object)),
	      as.integer(attr(object, "nugget")),
              as.integer(attr(object, "nuScaled")),
	      as.double(unlist(covariate)), #distances ?
	      as.integer(unlist(corD)),
	      as.double(attr(object, "minD")),
	      mat = double(corD[["sumLenSq"]]))[["mat"]]
    lD <- NULL
  } else {
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
  }
  if (corD[["M"]] > 1) {
    val <- split(val, rep(1:corD[["M"]], (corD[["len"]])^2))
    val <- lapply(val, function(el) {
      nel <- round(sqrt(length(el)))
      array(el, c(nel, nel))
    })
    names(val) <- names(corD[["len"]])
    val <- as.list(val)
  } else {
    val <- array(val, c(corD[["N"]], corD[["N"]]))
  }
  attr(val, "logDet") <- lD
  val
}
