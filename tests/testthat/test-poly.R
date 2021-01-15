cat(crayon::yellow("\ntest poly():\n"))

set.seed(123)
d <- data.frame(x = 1:10, y = rnorm(10), z=rnorm(10))

m1 <- fitme(y ~ x + I(x^2), data = d)
m2 <- fitme(y ~ poly(x, 2, raw = TRUE), data = d) # default is raw=FALSE
m3 <- fitme(y ~ poly(x, 2), data = d)
testthat::expect_true(diff(range(c(-11.1162027435,logLik(m1),logLik(m2),logLik(m3))))<1e-7) 

p1 <- predict(m1, newdata = data.frame(x = 2)) 
p2 <- predict(m2, newdata = data.frame(x = 2)) 
p3 <- predict(m3, newdata = data.frame(x = 2))  
testthat::expect_true(diff(range(c(0.1018852,p1,p2,p3)))<1e-7) 

m4 <- fitme(y ~ poly(cbind(x,z), 2, raw = TRUE), data = d)
m5 <- fitme(y ~ poly(cbind(x,z), 2), data = d)
testthat::expect_true(diff(range(c(-7.60735012544,logLik(m4),logLik(m5))))<1e-7) 


## with a random effect
set.seed(123)
data_test <- data.frame(x = rnorm(100), id = as.factor(1:100))
d <- data_test[rep(1:nrow(data_test), each = 10), ]
d$y <- rnorm(nrow(d))
modtp <- fitme(y ~ poly(x, 2) + (1|id), data = d)
new_data <- data.frame(x = 2, id = d$id[1])
abyss <- simulate(modtp, newdata = new_data[rep(1L,100),], type="residual")

## 
if (spaMM.getOption("example_maxtime")>6) {
  set.seed(123)
  data_test <- data.frame(x = factor(runif(100)>0.5), id = as.factor(1:100))
  d <- data_test[rep(1:nrow(data_test), each = 10), ]
  ranef <- structure(rnorm(length(unique(d$id)),sd = 0.666),names=levels(d$id))
  d$y <- rnorm(nrow(d))+ranef[d$id]
  modtp <- fitme(y ~ x + (1|id), data = d)
  new_data <- expand.grid(x = "FALSE", id = d$id[1])
  var(drop(simulate(modtp, newdata = new_data, nsim = 100000, type = "residual"))) ## OK ~phi=1.02722
  var(drop(simulate(modtp, newdata = new_data, nsim = 100000, type = "marginal"))) ## OK ~phi=1.02722+lambda= 0.3777
  mean(simulate(modtp, newdata = new_data, nsim = 100000, type = "marginal")) ## OK ~Intercept =0.06874 as x is "FALSE"
  mean(simulate(modtp, newdata = new_data, nsim = 100000, type = "residual")) ## OK ~Intercept =0.06874+ 0.178734=ranef(modtp)[[1L]][paste(d$id[1])]
}

# poly in the ranefs
if(requireNamespace("lme4", quietly = TRUE)) {
  data("sleepstudy",package = "lme4")
  # merfit <- lme4::lmer(Reaction ~ 1 + (poly(Days,2)|Subject), data = sleepstudy)
  m1 <- fitme(Reaction ~ 1 + (poly(Days,2)|Subject), data = sleepstudy, method="REML")
  m2 <- fitme(Reaction ~ 1 + (poly(Days,2, raw=TRUE)|Subject), data = sleepstudy, method="REML")
  testthat::expect_true(diff(range(c(-877.895899884,logLik(m1),logLik(m2))))<1e-9) 
  p1 <- predict(m1)
  p1n <- predict(m1, newdata=m1$data)
  testthat::expect_true(diff(range(p1-p1n))<1e-12) 
  p2 <- predict(m2)
  p2n <- predict(m2, newdata=m1$data)
  testthat::expect_true(diff(range(p2-p2n))<1e-12) 
  testthat::expect_true(diff(range(p1-p2))<3e-4) # not precise, but the variance estimates are large
}
