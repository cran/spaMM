`%*%.blockDiag` <-
function(mat,bdobject) {
  identifyCols <- c(0,cumsum(unlist(lapply(bdobject,ncol))))
  bd <- lapply(1:length(bdobject),function(v){mat[,seq(identifyCols[[v]]+1,identifyCols[[v+1]])] %*% bdobject[[v]]})
  attributes(bd) <- attributes(bdobject)
  bd
}
