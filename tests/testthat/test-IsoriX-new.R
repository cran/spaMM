if (Sys.getenv("_LOCAL_TESTS_")=="TRUE" && spaMM.options("example_maxtime")>40 && requireNamespace("IsoriX",quietly=TRUE)) { 
  cat(crayon::yellow("\ntest of IsoriX compatibility:"))
  testthat::test_local(path=paste0(projpath(), "/../IsoriXplus/IsoriX/IsoriX/tests/testthat/")) # test_NA_handling notably...
}
