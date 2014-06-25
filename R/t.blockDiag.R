t.blockDiag <-
function(x) {
  bd <- lapply(1:length(x),function(v){t(x[[v]])})
  attributes(bd) <- attributes(x)
  bd
}
