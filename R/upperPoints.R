upperPoints <-
function(Xpredy,minPtNbr=1#sqrt(nrow(Xpredy)+1)
                        ,D.resp) {
  respName <- attr(Xpredy,"respName")
  orderResp <- order(Xpredy[,respName], decreasing = TRUE)
  orderedXy <- Xpredy[orderResp,]
  Dy <- orderedXy[1,respName]-orderedXy[,respName]
  ptNbr <- max(c(which(Dy<D.resp)),minPtNbr) 
  orderResp[seq_len(ptNbr)] ## now returns indices
}
