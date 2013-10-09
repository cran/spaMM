coef.corMatern <-
function(object, unconstrained = TRUE, ...)
{
  ##cat("tata",object[1],object[2],"\n")
  if (attr(object, "fixed") && unconstrained) {
    return(numeric(0))
  }
  val <- as.vector(object)
  if (length(val) == 0) {               # uninitialized
    return(val)
  }
  if (!unconstrained) {
    val <- exp(val)
    if (attr(object, "nugget")) val[3] <- val[3]/(1+val[3])
  }
  if (attr(object, "nugget")) names(val) <- c("range", "nu", "nugget")
  else names(val) <- c("range","nu")
  val
}
