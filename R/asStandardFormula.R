asStandardFormula <-
function(formula) {
  aschar <- deparse(formula)
  aschar <- gsub("adjacency(","(",aschar,fixed=T)
  aschar <- gsub("Matern(","(",aschar,fixed=T)
  aschar <- gsub("corMatern(","(",aschar,fixed=T)
  as.formula(aschar)
}
