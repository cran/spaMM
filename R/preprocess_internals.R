# even though the Z's were sparse postmultplication by LMatrix leads some of the ZAL's to dgeMatrix (dense)
choose_QRmethod <- function(ZAlist, predictor, trySparse=TRUE) {
  if ( is.null(QRmethod <- .spaMM.data$options$QRmethod) ) { ## user setting. The code should NOT write into it. 
    if (trySparse) { ## 
      nrand <- length(ZAlist)
      if (nrand>0L && any(attr(attr(ZAlist,"ranefs"),"type") %in% c("Matern","corrMatrix") ) ) {
        if (sum(unlist(lapply(ZAlist,ncol)))>200L) { ## ad hoc
          QRmethod <- "lmwithQ_denseZAL"
        } else QRmethod <- "Matrix::qr"
      } else if (nrand==1L && is.identity(ZAlist[[1]])) {
        QRmethod <- "Matrix::qr" ## special case for poisson or binomial with saturated ranef
      } else if (nrand==1L && inherits(ZAlist[[1]],"sparseMatrix") && length(ZAlist[[1]]@x)/prod(dim((ZAlist[[1]])))<1/133) {
        QRmethod <- "Matrix::qr" ## recurrent case
      } else if (all(unlist(lapply(ZAlist,inherits,what="sparseMatrix")))) { ## OK pour Female/Male si Pdiag util sparse
        totsize <- prod(colSums(do.call(rbind,lapply(ZAlist,dim))))
        nonzeros <- sum(unlist(lapply(ZAlist,function(spm) {length(spm@x)})))          
        if (nonzeros/totsize<1/133) { ## of sparse class and effectively sparse 
          ## without the test, double le temps de calcul de exemples:
          QRmethod <- "Matrix::qr"
        } else {
          QRmethod <- "lmwithQ_denseZAL" ## of sparse class but actually not so sparse
          ## denseZAL still seems ~ 4 % faster on the tests
        }
      } else QRmethod <- "lmwithQ_denseZAL" 
      ## mais dans ce cas (08/2016) Matrix::qr ne change pas visiblement les temps
      ## even with CAR (etc), this is currently dense...
    } else QRmethod <- "lmwithQ_denseZAL"  
  }
  return(QRmethod)
}
