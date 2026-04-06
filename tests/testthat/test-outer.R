if (file.exists((privtest <- paste0(spaMM::projpath(),"/package/tests_private/test_outer.R")))) {
  source(privtest)
}