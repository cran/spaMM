chol.blockDiag <-
function(x,...) {
  bd <- lapply(1:length(x),function(v){chol(x[[v]],...)})
  attributes(bd) <- attributes(x)
  bd
}
