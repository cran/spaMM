.nonzeros <- function(spm) {
  if (inherits(spm,"ddiMatrix") && spm@diag=="U") {
    return(ncol(spm))
  } else if (inherits(spm,"sparseMatrix")) {
    nz <- length(spm@x)
    if ( methods::.hasSlot(spm, "diag")) nz <- nz+ncol(spm)
    return(nz) 
  } else return(sum(spm !=0)) ## AMatrix reaches here
}

# even though the Z's were sparse postmultplication by LMatrix leads some of the ZAL's to dgeMatrix (dense)
.choose_QRmethod <- function(ZAlist, predictor, trySparse=TRUE) {
  if ( is.null(QRmethod <- .spaMM.data$options$QRmethod) ) { ## user setting. The code should NOT write into it. 
    nrand <- length(ZAlist)
    if (trySparse && nrand>0L) {
      # adjacency speed to be tested on 2nd example from test-spaMM.R
      if (any(attr(ZAlist,"exp_ranef_types") %in% c("adjacency") ) ) {
        QRmethod <- "dense" ## la corr matrix est dense !
      } else {
        totdim <- colSums(do.call(rbind,lapply(ZAlist,dim)))
        if (any(attr(ZAlist,"exp_ranef_types") %in% c("Matern") ) ) {
          if (totdim[1L]>4*totdim[2L]) {
            QRmethod <- "sparse"
          } else QRmethod <- "dense" 
        } else if (nrand==1L && any(attr(ZAlist,"exp_ranef_types") %in% c("corrMatrix") ) ) {
          ## quik patch for this case, should be rethought
          if (totdim[2L]>200L) { ## ad hoc: we should use the type of corrMatrix or else scan its contents
            QRmethod <- "dense"
          } else QRmethod <- "sparse"
        } else if (nrand==1L && .is_identity(ZAlist[[1]])) { ## test pertinent slmt pour non-spatial models !
          QRmethod <- "sparse" ## special case for poisson or binomial with saturated ranef
        } else {
          totsize <- prod(totdim)
          nonzeros <- sum(unlist(lapply(ZAlist, .nonzeros)))          
          if (nonzeros/totsize < .spaMM.data$options$sparsity_threshold) { 
            QRmethod <- "sparse"
          } else {
            QRmethod <- "dense" ## ZAlist actually not so sparse
          }
        }
      } 
    } else QRmethod <- "dense" 
  }
  return(QRmethod)
}
