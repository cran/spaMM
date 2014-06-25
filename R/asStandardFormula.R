asStandardFormula <-
function(formula) {
  aschar <- DEPARSE(formula)
  aschar <- gsub("adjacency(","(",aschar,fixed=T)
  aschar <- gsub("Matern(","(",aschar,fixed=T)
  aschar <- gsub("corMatern(","(",aschar,fixed=T)
  aschar <- gsub("AR1(","(",aschar,fixed=T)
  aschar <- gsub("corrMatrix(","(",aschar,fixed=T)
  as.formula(aschar)
}
