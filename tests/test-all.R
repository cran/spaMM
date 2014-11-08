if(require("testthat", quietly = TRUE)) {
  pkg   <- "spaMM"
  require(pkg, character.only=TRUE, quietly=TRUE)
  ## test_package(pkg) ## for an installed package
  test_check(pkg) ## for R CMD check
  print(warnings()) # TODO? catch most of these by expect_warning(..)
} else {
  cat( "package 'testthat' not available, cannot run unit tests\n" )
}
