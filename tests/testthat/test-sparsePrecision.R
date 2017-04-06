cat("\ntest of sparse Precision method:")

data("scotlip")
fit1 <- fitme(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
               fixed=list(lambda=0.1), 
               adjMatrix=Nmatrix,family=poisson(),data=scotlip)
spaMM.options(wDEVEL=TRUE)
fit2 <- fitme(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
              fixed=list(lambda=0.1), 
              adjMatrix=Nmatrix,family=poisson(),data=scotlip)
spaMM.options(wDEVEL=FALSE)
expect_equal(logLik(fit1),logLik(fit2),tolerance=1e-05)