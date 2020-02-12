cat(crayon::yellow("\ntest old examples and new tests:\n"))
# spaMM

data("scotlip") ## loads 'scotlip' data frame, but also 'Nmatrix'

# Bug before 2.3.70 (NB_shape requested before optimization)
(hl <- fitme(I(1+cases)~I(prop.ag/10)+offset(log(expec))+adjacency(1|gridcode),
             control.HLfit = list(LevenbergM=FALSE), # otherwise "!LM" occurs and costs 7.5s; it's in test-LevM.R
            family=negbin(), adjMatrix=Nmatrix, data=scotlip)) ## explicit spaMM::negbin() may be needed.

(hl1 <- corrHLfit(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
          data=scotlip,family=poisson(),
          adjMatrix=Nmatrix)) ## 1D optimization -> optimize
testthat::expect_equal(hl1$APHLs$p_v,-161.5140,tolerance=1e-4)

(hl2 <- HLCor(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
            data=scotlip,family=poisson(), adjMatrix=Nmatrix))
testthat::expect_equal(hl2$APHLs$p_v,-161.5141,tolerance=1e-4)
testthat::expect_true(diff(range(AIC(hl2,verbose=FALSE)-AIC(hl1,verbose=FALSE)))<1e-2)

data("salamander")
hl <- HLfit(cbind(Mate,1-Mate)~1+(1|Female)+(1|Male),family=binomial(),
      rand.family=list(gaussian(),Beta(logit)),data=salamander,HLmethod="ML",control.HLfit = list(LevenbergM=FALSE))

testthat::expect_equal(hl$APHLs$p_v,-238.715,tolerance=1e-3)

## Nested effects
# lmer syntax allowing several degrees of nesting
hl <- HLfit(cbind(Mate,1-Mate)~1+(1|Female/Male),
      family=binomial(),rand.family=Beta(logit),data=salamander,HLmethod="ML",control.HLfit = list(LevenbergM=FALSE))

testthat::expect_equal(hl$APHLs$p_v,-243.6668,tolerance=1e-4)

# A syntax described in ?formula ## removed from the example()
hl <- HLfit(cbind(Mate,1-Mate)~1+(1|Female)+(1|Male %in% Female),
            ranFix=list(lambda=c('Female' = 0.127517,'Male %in% Female' = 4.64595e-07)),
            family=binomial(),rand.family=Beta(logit),data=salamander,HLmethod="ML",control.HLfit = list(LevenbergM=FALSE))

testthat::expect_equal(hl$APHLs$p_v,-243.6668,tolerance=1e-4)

### check NULL auglinmodblob
d <- data.frame(y = 1:10)
summary(fitme(y ~ 0, data = d))


# test of scaling of ranCoef predictor
(ranSlope1 <- fitme(I(1+cases)~I(prop.ag/10)+adjacency(0+expec|gridcode),
                    family=poisson(), adjMatrix=Nmatrix, data=scotlip))
scotlip$verif <- scotlip$expec/2
(ranSlope2 <- fitme(I(1+cases)~I(prop.ag/10)+adjacency(0+verif|gridcode),
                    family=poisson(), adjMatrix=Nmatrix, data=scotlip)) ## explicit spaMM::negbin() may be needed.
testthat::expect_true(abs(ranSlope2$lambda-4*ranSlope1$lambda)<1e-5)

