## to get the dim of the matrix of an object def as     corrMatrix <- structure(list(matrix=matrix_),class=c("list","precision")) 
dim.precision <- function(x) return(dim(x[["matrix"]]))

`[.precision` <- function(x,i,j, 
                          drop = TRUE ## by default, this function will return scalar/vector  
) {
  x[["matrix"]] <- x[["matrix"]][i,j,drop=drop]
  return(x)
}

## Inelegant but safe function to construct a sub-precision matrix. But best is not to have to use it.
.subset_prec <- function(prec_mat, corr_mat=NULL, rowcol_ids, tol=1e-10, verbose=TRUE) {
  rowcol_ids <- as.character(rowcol_ids) ## make sure it is not integer
  if (is.null(corr_mat)) {
    if (verbose) message("'corr_mat' has to be computed using solve(): this may be slow/and or exceed memory limits.\n See help('subset_prec') for Details.")
    corr_mat <- solve(prec_mat) ## corr_mat is better precomputed if many subsets are to be taken}
  }
  corr_mat <- drop0(corr_mat, tol=tol)
  sub_rel_mat <- corr_mat[rowcol_ids,rowcol_ids]
  prec_mat <- solve(sub_rel_mat) ## Matrix::solve; more elaborate methods might be implemented
  prec_mat <- drop0(prec_mat, tol=tol)
  colnames(prec_mat) <- rownames(prec_mat) <- rowcol_ids 
  return(prec_mat)
}
