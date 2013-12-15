Dim.corMatern <-
function(object, groups, ...)
{
  if (missing(groups)) return(attr(object, "Dim"))
  ugrp <- unique(groups)
  groups <- factor(groups, levels = ugrp)
  len <- table(groups)
  val <- list(N = length(groups),
              M = length(len),
              maxLen = max(len),
              sumLenSq = sum(len^2),
              len = len,
              start = NA,
              spClass=0)
  val[["start"]] <- c(0, 4, cumsum(val[["len"]] * (val[["len"]] - 1)/2)[-val[["M"]]])
  val
}
