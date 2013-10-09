recalc.corMatern <-
function(object, conLin, ...)
{
  ##print(sys.function(1))
  val <-
    .C("matern_recalc",
       Xy = as.double(conLin[["Xy"]]),
       as.integer(unlist(Dim(object))),
       as.integer(ncol(conLin[["Xy"]])),
       as.double(as.vector(object)),
       as.double(unlist(getCovariate(object))),
       as.double(attr(object, "minD")),
       as.integer(attr(object, "nugget")),
       as.integer(attr(object, "nuScaled")),
       logLik = double(1))[c("Xy", "logLik")]
  conLin[["Xy"]][] <- val[["Xy"]]
  conLin[["logLik"]] <- conLin[["logLik"]] + val[["logLik"]]
  conLin
}
