Initialize.corMatern <-
function(object, data, ...)
{
    if (!is.null(attr(object, "minD"))) { #already initialized
      return(object)
    }
    
    form <- formula(object)
    ## obtaining the groups information, if any
    if (!is.null(getGroupsFormula(form))) {
      attr(object, "groups") <- getGroups(object, form, data = data)
      attr(object, "Dim") <- Dim(object, attr(object, "groups"))
    } else {
      # no groups
      attr(object, "Dim") <- Dim(object, as.factor(rep(1, nrow(data))))
    }
    ## obtaining the covariate(s)
    ## is this where the distance matrix is actually computed ?
    attr(object, "covariate") <- getCovariate(object, data = data)

    nug  <- attr(object, "nugget")
    #nusc <- attr(object, "nuScaled")
    
    val <- as.vector(object)            # how many parameters?
    if (length(val) > 0) {		# is initialized
    if (val[1] <= 0) {                  # test for values of range
      stop("'range' must be > 0 in \"corMatern\" initial value")
    }
    if (val[2] <= 0) {                  # test for values of nu
      stop("'nu' must be > 0 in \"corMatern\" initial value")
    }    
    if (nug) {				# with nugget effect
      if (length(val) == 2) {		# assuming nugget effect not given
	    val <- c(val, 0.1)		# setting it to 0.1
      } else {
	if (length(val) != 3) {
	  stop("initial value for \"corMatern\" parameters of wrong dimension")
	}
      }
      if ((val[3] <= 0) || (val[3] >= 1)) {
	stop("initial value of nugget ratio must be between 0 and 1")
      }
    } else {				# only range and nu
      if (length(val) != 2) {
	stop("initial value for \"corMatern\" parameters of wrong dimension")
      }
    }
  } else { ## val of length 0
    val <- min(unlist(attr(object, "covariate"))) * 0.9
    val <- c(val,0.5) ## FR 18/03/13
    if (nug) val <- c(val, 0.1)
  }
  val[1] <- log(val[1]) ## range
  val[2] <- log(val[2]) ## nu    
  if (nug) val[3] <- log(val[3]/(1 - val[3]))
  oldAttr <- attributes(object)
  object <- val
  attributes(object) <- oldAttr
  attr(object, "minD") <- min(unlist(attr(object, "covariate")))
  attr(object, "factor") <- corFactor(object)
  attr(object, "logDet") <- -attr(attr(object, "factor"), "logDet")

  object
}
