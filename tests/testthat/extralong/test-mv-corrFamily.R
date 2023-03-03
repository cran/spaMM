cat(crayon::yellow(" -> test-mv-corrFamily:")) # not part of the testthat.R tests (neither test-composite-extra.R)

library(spaMM)
options(error=recover)

{
  data("blackcap")
  MLdistMat2 <- as.matrix(proxy::dist(blackcap[,c("latitude","longitude")]))
  MLcorMat2 <- MaternCorr(proxy::dist(blackcap[,c("latitude","longitude")]),
                          nu=0.6285603,rho=0.0544659)
  cap_mv <- blackcap
  cap_mv$name <- as.factor(rownames(blackcap))                
  cap_mv$grp <- 1L+(blackcap$migStatus>1)   
  set.seed(123)
  cap_mv$status2 <- blackcap$migStatus+ rnorm(14,sd=0.001)
}

{
  cat(crayon::yellow("MaternIMRFa; ")) 
  { # create IMRF model
    ## Creating the mesh 
    oldMDCopt <- options(Matrix.warnDeprecatedCoerce = 0) # # INLA issue
    mesh <- INLA::inla.mesh.2d(loc = blackcap[, c("longitude", "latitude")], 
                               cutoff=30,
                               max.edge = c(3, 20)) 
    mesh$n ## 40
    matern <- INLA::inla.spde2.matern(mesh)
    options(oldMDCopt)
  }
  { # first corrFamily mv fit
    (zut_cF <- fitmv(submodels=list(mod1=list(migStatus ~ 1),
                                    mod2=list(status2 ~ 1+ corrFamily(1|longitude+latitude))), # verbose=c(TRACE=TRUE), 
                     fixed=list(phi=c(0.02,0.02)), covStruct=list(corrFamily=MaternIMRFa(mesh=mesh, fixed=c(alpha=2))),
                     data=cap_mv))
    (zut_cF_lambdaFix <- fitmv(submodels=list(mod1=list(migStatus ~ 1),
                                              mod2=list(status2 ~ 1+ corrFamily(1|longitude+latitude))), # verbose=c(TRACE=TRUE), 
                               fixed=list(phi=c(0.02,0.02), lambda=30), covStruct=list(corrFamily=MaternIMRFa(mesh=mesh, fixed=c(alpha=2))),
                               data=cap_mv))
    (zut_cF_kappa_outer_fix <- fitmv(submodels=list(mod1=list(migStatus ~ 1),
                                                    mod2=list(status2 ~ 1+ corrFamily(1|longitude+latitude))), # verbose=c(TRACE=TRUE), 
                                     fixed=list(phi=c(0.02,0.02), corrPars=list("1"=c(kappa=0.5))), covStruct=list(corrFamily=MaternIMRFa(mesh=mesh)),
                                     data=cap_mv)) 
    
    (zut_cF_alpha_est <- fitmv(submodels=list(mod1=list(migStatus ~ 1),
                                              mod2=list(status2 ~ 1+ corrFamily(1|longitude+latitude))), # verbose=c(TRACE=TRUE), 
                               fixed=list(phi=c(0.02,0.02)), covStruct=list(corrFamily=MaternIMRFa(mesh=mesh)),
                               data=cap_mv)) 
  }
}