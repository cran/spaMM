cat(crayon::yellow("\ntest singular design matrices:\n"))

# Singular matrix from ?Matrix::qr :
singX <- cbind(int = 1,
           b1=rep(1:0, each=3), b2=rep(0:1, each=3),
           c1=rep(c(1,0,0), 2), c2=rep(c(0,1,0), 2), c3=rep(c(0,0,1),2))
rownames(singX) <- paste0("r", seq_len(nrow(singX)))
donn <- as.data.frame(singX)
set.seed(123)
donn$y <- runif(6)
fv1 <- fitted(lm(y~int+ b1+b2+c1+c2+c3,data=donn))
spaMM.options(rankMethod="qr")
fv2 <- fitme(y~int+ b1+b2+c1+c2+c3,data=donn)$fv
spaMM.options(rankMethod=".rankinfo")
fv3 <- fitme(y~int+ b1+b2+c1+c2+c3,data=donn)$fv ## fixef is not unique, but fitted values must be equivalent
spaMM.options(rankMethod="qr")
testthat::expect_true( max(apply(cbind(fv1,fv2,fv3),1L,var))<1e-16)

donn$dummy<- c(0,0,0,1,1,1)
chk <- try(nlme::lme(y~int+ b1+b2+c1+c2+c3, data = donn, random = ~ 1 | dummy), silent=TRUE)
testthat::test_that("check claim, from the Description of spaMM::rankinfo, that nlme::lme does not handle singular X",
                    testthat::expect_true(inherits(chk,"try-error"))
)
# chk <- try(lme4::lmer(y~int+ b1+b2+c1+c2+c3 + (1 | dummy), data = donn), silent=TRUE), warnings, but this fits contrary to lme.

