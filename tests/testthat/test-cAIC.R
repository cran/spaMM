cat(crayon::yellow("test cAIC:\n"))

if (spaMM.getOption("example_maxtime")>3) {
  if(requireNamespace("lme4", quietly = TRUE)) { 
    data("sleepstudy",package = "lme4")
    # lb <- lme4::lmer(Reaction ~ Days + (1 | Subject), data=sleepstudy)
    # cAIC4::cAIC(lb)
    # Conditional log-likelihood:  -864.53
    # Degrees of freedom:    19.03
    # Conditional Akaike information criterion:  1767.12
    # 2*864.53+2*19.03 # that's it, but one df seems to be missing for the phi estimation
    
    #cAIC(lb, method="conditionalBootstrap")
    #b <- fitme(Reaction ~ Days + (1 | Subject), data=sleepstudy, method="REML")
    b <- fitme(Reaction ~ Days + (1 | Subject), data=sleepstudy, method="REML", fixed=list(lambda=1378)) # conditional AIC: 1766.8439
    get_any_IC(b) # cAIC is  -2*forAIC$clik + 2*(pd+p_phi) while eff df is length(object$y) - pd ; pd = 17.89231
    get_any_IC(b,nsim=100)
  }
  
  #### Non-canonical link
  # z <- fitme(y ~ 1+(1|batch), family=Gamma(log), data=wafers)
  z <- fitme(y ~ 1+(1|batch), family=Gamma(log), data=wafers, fixed=list(lambda=0.01212))
  get_any_IC(z, nsim=100L)
}
