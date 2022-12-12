cat(crayon::yellow("\ntest singular design matrices:\n"))

# Singular matrix from ?Matrix::qr :
singX <- cbind(int = 1,
           b1=rep(1:0, each=3), b2=rep(0:1, each=3),
           c1=rep(c(1,0,0), 2), c2=rep(c(0,1,0), 2), c3=rep(c(0,0,1),2))
rownames(singX) <- paste0("r", seq_len(nrow(singX)))
donn <- as.data.frame(singX)
set.seed(123)
donn$y <- runif(6)
fv1 <- fitted(singlm <- lm(y~int+ b1+b2+c1+c2+c3,data=donn)) # (glmmTMB returns a full vector with numeric values instead of NA's)... vcov will be full of NaN's
spaMM.options(rankMethod="qr")
fv2 <- (singfit <- fitme(y~int+ b1+b2+c1+c2+c3,data=donn))$fv
spaMM.options(rankMethod=".rankinfo")
fv3 <- fitme(y~int+ b1+b2+c1+c2+c3,data=donn)$fv ## fixef is not unique, but fitted values must be equivalent
spaMM.options(rankMethod="qr")
testthat::expect_true( max(apply(cbind(fv1,fv2,fv3),1L,var))<1e-16)
crit <- diff(range(anova(singlm) -anova(singfit), na.rm=TRUE))
testthat::test_that("whether anova.lm and spaMM:::.anova.lm give equivalent results for singular X",
                    testthat::expect_true(crit<1e10))

# the way singular matrices are handled differ for LM vs GLMs: for the latter there are rows with 0 df in the Table
fv1 <- fitted(singglm <- glm(y~int+ b1+b2+c1+c2+c3,data=donn, family=Gamma(log)))
spaMM.options(rankMethod="qr")
fv2 <- (singfit <- fitme(y~int+ b1+b2+c1+c2+c3,data=donn, family=Gamma(log)))$fv
spaMM.options(rankMethod=".rankinfo")
fv3 <- fitme(y~int+ b1+b2+c1+c2+c3,data=donn, family=Gamma(log))$fv ## fixef is not unique, but fitted values must be equivalent
spaMM.options(rankMethod="qr")
testthat::expect_true( max(apply(cbind(fv1,fv2,fv3),1L,var))<1e-9)
crit <- diff(range(anova(singglm) -anova(singfit), na.rm=TRUE)) #
testthat::test_that("whether anova.lm and spaMM:::.anova.lm give equivalent results for singular X",
                    testthat::expect_true(crit<1e10))

donn$dummy<- c(0,0,0,1,1,1) # completely équivalent to b1 or b2!
chk <- try(nlme::lme(y~int+ b1+b2+c1+c2+c3, data = donn, random = ~ 1 | dummy), silent=TRUE)
testthat::test_that("whether nlme::lme does not handle singular X, as stated in the Description of spaMM::rankinfo",
                    testthat::expect_true(inherits(chk,"try-error"))
)

if (spaMM.getOption("example_maxtime")>0.7) { # a check accumulating possible causes of problems
  { # check on minipal pathological case
    ## b1 is redundant with | dummy ...
    trivml <- fitme(y~b1 + (1 | dummy), data = donn, #verbose=c(TRACE=interactive()), 
                    method="ML") # There is info about lambda, which is estimated as zero
    numInfo(trivml) # the nonzero gradient at the boundary is detected. lambda is removed from target params, so that num derivation does not generate negative lambdas
    trivreml <- fitme(y~b1 + (1 | dummy), data = donn, #verbose=c(TRACE=interactive()), 
                      method="REML") # There is NO info kept in Re.L about lambda
    numInfo(trivreml,which=c("lambda","phi","beta")) # the haphazard 'estimate' of lambda is far from 0, so is included, with vanishing 2nd derivative. 
                      # The phi result is consistent with the next one
                      # The beta result is different from the next one as expected from the different lambda.
    trivreml <- fitme(y~b1 + (1 | dummy), data = donn, control=list(refit=TRUE), method="REML") 
    numInfo(trivreml, which=c("lambda","phi","beta")) # values at boundary are excluded from computation to avoid negative lambda
  }
  
  ## Each of b1 and b2 are redundant with | dummy ...
  # Here lambda appears unidentifiable by *RE*ML: RE.L is 'completely' flat wrt it (marginal L is not). It would be nice if a message pointed that (___F I X M E___):
  # There is one only for singfitT bc refit=TRUE -> leverages are computed -> detection of unit leverages.
  (singfitF <- fitme(y~int+ b1+b2+c1+c2+c3 + (1 | dummy), data = donn, method="REML", 
                    # verbose=c(TRACE=TRUE), # quite useful to check what's going on in numInfo computation too
                    init=list(lambda=0.01), 
                    control=list(refit=FALSE))) ## Final value of lambda dependent on initial value
  suppressWarnings(anova(singfitF, type="2")) # lambda not low enough for removal, but numInfo singular since no info on lambda => warning (suppressed here)
  ranef(singfitF) #  =0, which is consistent with the refit below. The latter de facto replaces an undefined lambda estimate ((v=0)/(1-lev=0)) by v/(tiny value) = 0
  (singfitT <- fitme(y~int+ b1+b2+c1+c2+c3 + (1 | dummy), data = donn, method="REML", 
                    # verbose=c(TRACE=TRUE), # quite useful to check what's going on
                    init=list(lambda=0.01), 
                    control=list(refit=TRUE))) # Note effect on lambda
  suppressWarnings(anova(singfitT, type="2")) # Distinct warning for removal of lambda at boundary
  suppressWarnings(chk <- lme4::lmer(y~int+ b1+b2+c1+c2+c3 + (1 | dummy), data = donn)) # warnings, but this fits contrary to lme.
  suppressWarnings(anova(lmerTest::as_lmerModLmerTest(chk), type="II")) 
  # Trying to get the same result as lme4... not completely successful, but there are good reasons for that: 
  (singfit <- fitme(y~int+ b1+b2+c1+c2+c3 + (1 | dummy), data = donn, method="REML", 
                    init=list(lambda=0.0164^2), 
                    lower=list(lambda=0.0163^2),
                    upper=list(lambda=0.0165^2),
                    control=list(refit=FALSE))) ## Final value of lambda = initial value
  suppressWarnings(anova(singfit, type="2")) # Again, lambda not low enough for removal, but numInfo singular since no info on lambda => warning
  suppressWarnings(anova(singfit, type="2",transf=FALSE)) # Note the difference in DenDf hence in p-value
  # Since lambda is not quite low, there is no check of the num derivs, despite all the identifiability issues. The difference is on the b1 coeff, redundant with the ranef...
}

if (spaMM.getOption("example_maxtime")>0.7) { # Test of correct handling of singular fits in anova lead to an odd F test| cX2 test comparison on non-singular fits.
  # Here for devel purposes a single collinearity is retained ('int' is equiv to the Intercept). 
  #
  ## Also b1 was redundant with | dummy and was removed here (this does not prevent lambda -> 0 unless we modify y: see commented-out line. 
  donn$yy <- donn$y
  # donn$yy[donn$dumm==1] <- donn$y[donn$dumm==1]+1 # if we want to avoid the zero lambda estimate
  ##
  okchk <- lme4::lmer(yy~int+c1+c2 + (1 | dummy), data = donn, REML=FALSE) 
  (okfit <- fitme(yy~c1+c2 + (1 | dummy), data = donn, method="ML")) 
  (singfit <- fitme(yy~int+ c1+c2 + (1 | dummy), data = donn, method="ML")) 
  # library("lmerTest")
  a1 <- anova(lmerTest::as_lmerModLmerTest(okchk), type="1") # consistent (despite possibly different treatment of lambda) with:
  suppressWarnings(a2 <- anova(okfit, type="1")) 
  suppressWarnings(a3 <- anova(singfit, type="1")) 
  crit <- max(abs(c(a1[1:2,6]-a2[1:2,6],a1[1:2,6]-a3[1:2,6])))
  testthat::test_that("whether LMM type 1 anova for singular X is consistent with trimmed X",
                      testthat::expect_true(crit<1e8))
  #
  anova(okfit, type="1", method="t.Chisq")
  # Remarkably, the X2 stats are identical to the previous F statz. 
  # The X2 stat comes from .anova_fallback -> Wald's method, hard to mess with.
  # The F stat is the square of the t-stat computed by lmerTest:::anova.lmerModlmerTest() -> single_anova() -> contest1D() -> [locally defined]mk_ttable()
  
  # same compar for type 2, quite different from type 1 ...
  a1 <- anova(lmerTest::as_lmerModLmerTest(okchk), type="2")
  suppressWarnings(a2 <- anova(okfit, type="2")) # 
  suppressWarnings(a3 <- anova(singfit, type="2"))
  crit <- max(abs(c(a1[1:2,6]-a2[1:2,6],a1[1:2,6]-a3[1:2,6])))
  testthat::test_that("whether LMM type 2 anova for singular X is consistent with trimmed X",
                      testthat::expect_true(crit<1e8))
  #
  anova(okfit, type="2", method="t.Chisq")
  
}
