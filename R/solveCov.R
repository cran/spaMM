solveCov <-
function(Sa,Sb,Sab,tra,trb,trab,trba,lambda_a,lambda_b,cov) {
#print(c(Sa,Sb,Sab,tra,trb,trab,trba,lambda_a,lambda_b,cov))  
  for (it in seq_len(1)) {
    next_lambda_a <- (lambda_b^2 *Sa - 2*cov*lambda_b*Sab+cov^2 * (Sb+lambda_b*tra)+cov^3 *trba) /(lambda_b^2*tra+lambda_b*cov*trba)
    next_lambda_b <- (lambda_a^2 *Sb - 2*cov*lambda_a*Sab+cov^2 * (Sa+lambda_a*trb)+cov^3 *trab) /(lambda_a^2*trb+lambda_a*cov*trab)
    ## doute sur le sign(): ambiguite du calcul analytique
    next_cov <- (2*Sab+next_lambda_a*trba+next_lambda_b*trab)/(tra+trb) ## approx lin for small cov
    lambda_a <- next_lambda_a
    lambda_b <- next_lambda_b
print(c(next_cov,lambda_a,lambda_b))    ## jamais marche
    cov <- next_cov 
  }
  resu <- c(lambda_a=lambda_a,lambda_b=lambda_b,cov=cov)
  return(resu)
}
