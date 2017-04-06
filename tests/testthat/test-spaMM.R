cat("\ntest spaMM.Rd old donttest examples:\n")
# spaMM

data(scotlip) ## loads 'scotlip' data frame, but also 'Nmatrix'
hl <- corrHLfit(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
          data=scotlip,family=poisson(),
          adjMatrix=Nmatrix,control.corrHLfit = list(maxcorners=0))
expect_equal(hl$APHLs$p_v,-161.5140,tolerance=1e-4)
AIC(hl)

## Adding a Gamma random effect to fit a negative-binomial response:
# (slow ~~25s with corners) 
# (and the best fit is not Gamma ranef, which is sadly missed with the default bounds: TODO here)
hl <- corrHLfit(cases~I(prop.ag/10) +(1|gridcode)+adjacency(1|gridcode)
          +offset(log(expec)),
          data=scotlip,family=poisson(),rand.family=list(Gamma(log),gaussian()),
          adjMatrix=Nmatrix,lower=list(rho=0),upper=list(rho=0.1745),control.corrHLfit = list(maxcorners=0))

expect_equal(hl$APHLs$p_v,-161.5141,tolerance=1e-4)

data(salamander)
hl <- HLfit(cbind(Mate,1-Mate)~1+(1|Female)+(1|Male),family=binomial(),
      rand.family=list(gaussian(),Beta(logit)),data=salamander,HLmethod="ML")

expect_equal(hl$APHLs$p_v,-238.715,tolerance=1e-3)

## Nested effects
# lmer syntax allowing several degrees of nesting
hl <- HLfit(cbind(Mate,1-Mate)~1+(1|Female/Male),
      family=binomial(),rand.family=Beta(logit),data=salamander,HLmethod="ML")

expect_equal(hl$APHLs$p_v,-243.6668,tolerance=1e-4)

# A syntax described in ?formula ## removed from the example()
hl <- HLfit(cbind(Mate,1-Mate)~1+(1|Female)+(1|Male %in% Female),
            ranFix=list(lambda=c('Female' = 0.127517,'Male %in% Female' = 4.64595e-07)),
            family=binomial(),rand.family=Beta(logit),data=salamander,HLmethod="ML")

expect_equal(hl$APHLs$p_v,-243.6668,tolerance=1e-4)
