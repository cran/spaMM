if (file.exists((privtest <- paste0(projpath(),"/package/tests_private/test_outer.R")))) {
  source(privtest)
}