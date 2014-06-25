legend_lambda <-
function(urff) {
  if (length(urff)==1) {
    if (urff=="beta") {
      #        cat("Coefficients for log[ lambda ] for u ~ Beta(1/(2lambda),1/(2lambda))\n")
      cat("Coefficients for log[ lambda = 4 var(u)/(1 - 4 var(u)) ]:\n")
    } else if (urff=="inverse.gamma") {
      #        cat("Coefficients for log[ lambda ] for  u ~ I-Gamma(1+1/lambda,1/lambda)\n")
      cat("Coefficients for log[ lambda = var(u)/(1 + var(u)) ]:\n")
    } else if (urff=="gamma") {
      #        cat("Coefficients for log[ lambda = var(u) ] for  u ~ Gamma(lambda,1/lambda)\n")
      cat("Coefficients for log[ lambda = var(u) ]:\n")
    } else cat("Coefficients for log[ lambda = var(u) ]: \n") 
  } else {
    cat("Coefficients for log[ lambda ], with:\n")
    lapply(urff, function(st) {
      switch(st,
             "beta" = cat("lambda = 4 var(u)/(1 - 4 var(u)) for Beta distribution; \n"),
             "inverse.gamma" = cat("lambda = var(u)/(1 + var(u)) for inverse gamma distribution; \n"),
             "gamma" = cat("lambda = var(u) for gamma distribution; \n"),
             "gaussian" = cat("lambda = var(u) for Gaussian distribution; \n")
      )
    })     
  }  
}
