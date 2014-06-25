as.blockDiag.partition <-
function(partition,mat,fill=0L) { ## 
  bd <- lapply(seq_len(length(partition)-1),function(v) {
    sequ <- (partition[v]+1):partition[v+1] 
    mat[sequ,sequ,drop=FALSE]
  })
  attr(bd,"rowNames") <- rownames(mat)  
  attr(bd,"cumsum") <- partition 
  attr(bd,"fill") <- fill
  class(bd) <- c("blockDiag",class(bd))
  bd  
}
