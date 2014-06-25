as.blockDiag.matrix <-
function(mat) {
  partition <- findblocks(mat)
  partition <- cumsum(c(0,partition))
  as.blockDiag.partition(partition,mat)
}
