cat(crayon::yellow("\ntest-nested-geostat:"))

data("blackcap")
grouped <- cbind(blackcap,grp=c(rep(1,7),rep(2,7))) 

(fit1 <- fitme(migStatus ~ 1 + Matern(1|longitude+latitude %in% grp),data=grouped,
               fixed=list(nu=4,rho=0.4,phi=0.05)))
(crit <- abs(diff(c(-13.02094, logLik(fit1)))))
testthat::test_that("Check whether any bug is being introduced in nested-geostat code (fit)",testthat::expect_true(crit<1e-05)) 

p1 <- predict(fit1,newdata=fit1$data)[,1]
p2 <- predict(fit1)[,1]

(fit2 <- fitme(migStatus ~ 1 + Matern(1|longitude+latitude %in% grp),data=grouped,
               fixed=list(nu=4,rho=c(0.4,0.4),phi=0.05)))
p3 <- predict(fit2,newdata=fit2$data)[,1]
p4 <- predict(fit2)[,1]

(crit <- abs(diff(range(p2-p1,p3-p1,p4-p1))))
testthat::test_that("Check whether any bug is being introduced in nested-geostat code (predict)",
                    testthat::expect_true(crit<1e-13)) 

# see devel/nested_Matern/ for a more memory-challenging example.

