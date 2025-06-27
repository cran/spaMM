data(wafers) 
{
  glmfit <- glm(I(y/1000)~X1, family=gaussian(), data=wafers)
  pw_a <- deviance(glmfit) #           3...                       (a)
  pw_b <- sum(residuals(glmfit)^2) #   3...                       (b) 
  
  # Same model, with different parametrization of residual variance 
  glmfit2 <- glm(I(y/1000)~X1, family=gaussian(), data=wafers, weights=rep(2,198))
  pw_c <- deviance(glmfit2) #          6...                       (c)  
  pw_d <- sum(residuals(glmfit2)^2) #  6...                       (d)
  
  # Same comparison but for HLfit objects:
  spfit <- fitme(I(y/1000)~X1, family=gaussian(), data=wafers)
  pw_e <- deviance(spfit) #            3...                       (e)
  pw_f <- sum(residuals(spfit)^2) #    3...                       (f)  
  pw_g. <- sum(dev_resids(spfit)) #    3...                         
  pw_r_a <- residVar(spfit)[1] #      0.015...
  pw_r_b <- get_residVar(spfit)[1] #  0.015...
  
  spfit2 <- fitme(I(y/1000)~X1, family=gaussian(), data=wafers, prior.weights=rep(2,198))
  pw_g <- deviance(spfit2) #           6...                       (g) ~ (c,d) # post v4.2.0
  pw_h <- sum(residuals(spfit2)^2) #   6...                       (h) ~ (c,d) 
  pw_i <- sum(dev_resids(spfit2)) #    3...                         
  pw_r_c <- head(residVar(spfit2)) #     0.015...
  pw_r_d <- head(get_residVar(spfit2)) # 0.030...
  
  # Unscaled residuals should not depend on arbitrarily fixed residual variance:
  spfit3 <- fitme(I(y/1000)~X1, family=gaussian(), data=wafers, fixed=list(phi=2),
                  prior.weights=rep(2,198))
  pw_j <- deviance(spfit3) #           6...                       (j) ~ (g)
  pw_k <- sum(residuals(spfit3)^2) #   6...                       (k) ~ (h)
  pw_l <- sum(dev_resids(spfit3)) #    3...                         
  residVar(spfit3)[1] #     1
  get_residVar(spfit3)[1] # 2
  
  crit <- diff(range(c(pw_a,pw_b,pw_e,pw_f,pw_g.,pw_i,pw_l,
               c(pw_c,pw_d,pw_g,pw_h,pw_j,pw_k)/2)))
  testthat::test_that("Correct handling of pw in computations of residuals",
                      testthat::expect_true(crit<1e-14))
  crit <- diff(range(c(pw_r_a,pw_r_b,pw_r_c,pw_r_d/2)))
  testthat::test_that("Correct handling of pw in computations of residVar/get_residVar",
                      testthat::expect_true(crit<1e-14))
  
}
