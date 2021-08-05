cat(crayon::yellow("test-predVar-Matern-corrMatrix"))
if (spaMM.getOption("example_maxtime")>0.7) { ##  not based on real timing
  
  # checks predVar w/o permutation, + Matern vs corrMatrix, + spprec T/F  ... + perm_Q
  
  data("blackcap")
  MLcorMat <- MaternCorr(proxy::dist(blackcap[,c("latitude","longitude")]),
                         nu=4,rho=0.4)
  blackcap$name <- as.factor(rownames(blackcap))
  bidon <- list(T=TRUE)
  old_perm_Q <- spaMM.getOption("perm_Q")
  for (nam in c("T","N")) { # so it ends on default, NULL
    spaMM.options(perm_Q=bidon[[nam]])
    ## (1) Single variable present in the data 
    (f1 <- HLCor(migStatus ~ means+ corrMatrix(1|name),data=blackcap,
                 corrMatrix=MLcorMat,method="ML", 
                 control.HLfit=list(sparse_precision=TRUE)))
    f2 <- corrHLfit(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap,
                    ranFix=list(corrPars=list("1"=list(nu=4,rho=0.4))),method="ML")
    f3 <- HLCor(migStatus ~ means+ corrMatrix(1|name),data=blackcap,
                corrMatrix=MLcorMat,HLmethod="ML", 
                control.HLfit=list(sparse_precision=FALSE))
    f4 <- corrHLfit(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap,
                    ranFix=list(corrPars=list("1"=list(nu=4,rho=0.4))),method="ML")
    f5 <- HLCor(migStatus ~ means+ corrMatrix(1|longitude+latitude),data=blackcap,
                corrMatrix=MLcorMat,method="ML") # Check that order of data is respected in the Zmatrix for this "unsafe" input.
    # imput as precision matrix
    f6 <- fitme(migStatus ~ means + corrMatrix(1|name), data=blackcap,
                covStruct=list(precision=as_precision(MLcorMat)))
    # Manual version of the same:
    as_mat <- proxy::as.matrix(MLcorMat, diag=1) 
    prec_mat <- solve(as_mat) ## precision factor matrix
    f7 <- fitme(migStatus ~ means + corrMatrix(1|name), data=blackcap,
                covStruct=list(precision=prec_mat))
    # the last two fitme spprec fits ssslightly differ from previous ones, even from f1,f2 spprec fits => its a fitme vs others difference. 
    testthat::expect_true(diff(range(c(logLik(f1),logLik(f2),logLik(f3),logLik(f4),logLik(f5),logLik(f6),logLik(f7))))<1e-8) # 
    testthat::expect_true(diff(range(c(predict(f1, newdata=f1$data[1,]), predict(f2, newdata=f1$data[1,]), predict(f3, newdata=f1$data[1,]),
                                       predict(f4, newdata=f1$data[1,]), predict(f5, newdata=f1$data[1,]), predict(f6, newdata=f1$data[1,]), 
                                       predict(f7, newdata=f1$data[1,]))))<1e-5)  
    #  I once had a problem with perm_Q=TRUE in test-predVar-Matern-corrMatrix -> predict(f2,.) 
    # so make sure to check this when changing relevant code. Seems OK now.
    # New ZA was Z_1x1 * A_nxn[identified row] -> 1xn and cov_newLv_oldv_list was 1xn
    # so the product failed. The new ZA should instead be 1x1.
    (c4a <- get_predVar(f4,variances=list(cov=TRUE))[3:2,3:2])
    (c4b <- get_predVar(f4,newdata=f4$data[c(3,2),],variances=list(cov=TRUE)))
    (c1 <- get_predVar(f1,newdata=f4$data[c(3,2),],variances=list(cov=TRUE)))
    (c2 <- get_predVar(f2,newdata=f4$data[c(3,2),],variances=list(cov=TRUE)))
    (c3 <- get_predVar(f3,newdata=f4$data[c(3,2),],variances=list(cov=TRUE)))
    (c5 <- get_predVar(f5,newdata=f4$data[c(3,2),],variances=list(cov=TRUE)))
    (c6 <- get_predVar(f6,newdata=f4$data[c(3,2),],variances=list(cov=TRUE)))
    (c7 <- get_predVar(f7,newdata=f4$data[c(3,2),],variances=list(cov=TRUE)))
    testthat::expect_true(diff(range(c(c1[1],c2[1],c3[1],c4a[1],c4b[1],c5[1],c6[1],c7[1])))<1e-5)  
    
  }
  spaMM.options(perm_Q=old_perm_Q)
}



data("blackcap") 
## Here we manually reconstruct the correlation matrix 
##  of the ML fit produced by corrHLfit:
MLcorMat <- MaternCorr(proxy::dist(blackcap[,c("latitude","longitude")]),
                       nu=0.6285603,rho=0.0544659)
blackcap$name <- as.factor(rownames(blackcap))                
## (1) Single variable present in the data 
HLCor(migStatus ~ means+ corrMatrix(1|name),data=blackcap,
      corrMatrix=MLcorMat,method="ML")
## (2) Same, permuted: still gives correct result
perm <- sample(14)
# Permuted matrix (with permuted names)
pmat <- as.matrix(MLcorMat)[perm,perm]
HLCor(migStatus ~ means+ corrMatrix(1|name),data=blackcap,
      corrMatrix=as.dist(pmat),method="ML")
## (3) Other grouping terms (note the messages):
HLCor(migStatus ~ means+ corrMatrix(1|longitude+latitude),data=blackcap,
      corrMatrix=MLcorMat,method="ML")
