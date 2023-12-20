if (Sys.getenv("_LOCAL_TESTS_")=="TRUE" && spaMM.getOption("example_maxtime")>40 && 
    requireNamespace("IsoriX",quietly=TRUE)) { 
  cat(crayon::yellow("\nRuning the IsoriX/tests/testthat/ tests (with IsoriX option spaMM_debug = TRUE):"))
  IsoriX::options_IsoriX(spaMM_debug = TRUE) # to show the warnings from spaMM.
  testthat::test_local(path=paste0(spaMM::projpath(), "/../IsoriXplus/whatGitPulls/IsoriX/tests/testthat/")) # test_NA_handling notably...
} else if (Sys.getenv("_LOCAL_TESTS_")!="TRUE") message("IsoriX tests not run: check whether _LOCAL_TESTS_=TRUE in .Renviron") # usethis::edit_r_environ() ...
# older tests:
if (FALSE) {
  source(file=paste0(spaMM::projpath(), "/package/tests_private/test-IsoriX.R"))
  # itself sourcing files in /nestedFiles
}
