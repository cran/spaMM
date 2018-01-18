## to get the dim of the matrix of an object def as     corrMatrix <- structure(list(matrix=matrix_),class=c("list","precision")) 
dim.precision <- function(x) return(dim(x[["matrix"]]))

`[.precision` <- function(x,i,j, 
                          drop = TRUE ## by default, this function will return scalar/vector  
) {
  x[["matrix"]] <- x[["matrix"]][i,j,drop=drop]
  return(x)
}


