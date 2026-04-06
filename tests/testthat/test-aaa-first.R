# test this before anything else to catch what would be a very confusing bug.

foo <- list2env(list(a=666))
delayedAssign("b", {333}, assign.env = foo)
# ls(foo)
..is_evaluated <- get(".is_evaluated", asNamespace("spaMM"), inherits=FALSE) 
testthat::expect_true(..is_evaluated("a",foo))
testthat::expect_true( ! ..is_evaluated("b",foo))
