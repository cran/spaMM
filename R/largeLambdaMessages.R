largeLambdaMessages <-
function() {
  message("A too high lambda may indicate a very poorly fitting fixed-effect part of the model.")
  message("To control the maximum lambda, use e.g. 'spaMM.options(maxLambda=1e06)'.")
  message("It may also indicate convergence issues, possibly improved by altering the initial value through the 'init.HLfit' argument.")
}
