cat(crayon::yellow("\ntest distMatrix:\n"))

data("blackcap")
MLdistMat <- as.matrix(proxy::dist(blackcap[,c("latitude","longitude")]))
fit_d1 <- fitme(migStatus ~ means+ Matern(1|latitude+longitude),data=blackcap,
      distMatrix=MLdistMat, method="ML", fixed =list(nu=0.6285603))
fit_d2 <- fitme(migStatus ~ means+ Matern(1|latitude+longitude),data=blackcap,
                distMatrix=2*MLdistMat, method="ML", fixed =list(nu=0.6285603))
testthat::expect_true(diff(c(get_ranPars(fit_d1, which="corrPars")[[1]]$rho,
                             2*get_ranPars(fit_d2, which="corrPars")[[1]]$rho))<1e-8)
