fitted.HLfitlist <-
function(object,...) {
  allrownames <- unique(unlist(lapply(object,function(hl){rownames(hl$data)})))
  fv <- lapply(object,function(hl){
    fv <- hl$fv; ## FR->FR ideally I should use fitted.HLfit(hl,...) here
    rownames(fv) <- rownames(hl$data); fv
  })
  resu <- matrix(0,nrow=length(allrownames),ncol=length(object))
  rownames(resu) <- allrownames
  for (it in seq_len(length(fv))) {
    resu[rownames(fv[[it]]),it] <- fv[[it]]
  }
  for (it in seq_len(ncol(resu)-1)) {
    for (col in (it+1):ncol(resu)) resu[,col] <- resu[,col] * (1-resu[,it])
  }
  colnames(resu) <- names(object)
  resu <-cbind(resu,1-rowSums(resu))
  resu
}
