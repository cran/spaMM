corMatern <-
function(value = c(1,0.5), form = ~1, nugget=FALSE, nuScaled=FALSE, metric = c("euclidean", "maximum", "manhattan"), fixed = FALSE)
{
    attr(value, "formula") <- form
    attr(value, "nugget") <- nugget
    attr(value, "nuScaled") <- nuScaled
    attr(value, "metric") <- match.arg(metric)
    attr(value, "fixed") <- fixed
    class(value) <- c("corMatern", "corStruct")
    ##class(value) <- c("corMatern", "corSpatial", "corStruct")

    value
}
