RcppChol.blockDiag <-
function(bdobject) {
  ## chol for each block
  blob <- lapply(1:length(bdobject),function(v){
    glop <- RcppChol(bdobject[[v]])
    if (glop$Status==TRUE) {colnames(glop$L) <- rownames(glop$L) <- rownames(bdobject[[v]])}
    glop
  }) 
  ## collates results
  statuses <- lapply(blob,function(v) {v$Status})
  Status <- all(unlist(statuses)) 
  if (Status==TRUE) {
    bd <- lapply(blob,function(v) {v$L})   ## $L !
    attributes(bd) <- attributes(bdobject)
    ## return(list(L=bd,Status=Status)) ## FR->FR 03/2014 on veut ça et ensuite une ZAL en blocks pour pouvoir travailler sur des blocks de wAugZALI
    ## mais compute.ZAL n'est pas up to this. 
    L <- as.matrix(bd)
    return(list(L=L,Status=Status))
  } else return(list(Status=Status))
}
