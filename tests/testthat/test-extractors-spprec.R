cat(crayon::yellow("test-extractors-spprec.R (notably ranCoefs)"))
{
  data("blackcap")
  ## make sure phi estimates are high.
  fake <- blackcap
  fake$grp <- rep(1:2,7)
  fake$migStatus <- fake$migStatus +(fake$grp==2)
  
  fake$ID <- gl(7,2)
  fake$grp <- factor(fake$grp)
  (dd <- fitme(migStatus ~ 1 +  (0+grp|ID),data=fake, control.HLfit=list(sparse_precision=FALSE), fixed=list(phi=0.1)))
  (ss <- fitme(migStatus ~ 1 +  (0+grp|ID),data=fake, control.HLfit=list(sparse_precision=T), fixed=list(phi=0.1)))
  ranef(ss)
  ranef(dd)
  
  # check coherence of structList, ZAL matrix and predict : for spprec and correl algo:
  p1 <- predict(dd)
  p2 <- predict(ss)
  p3 <- fixef(dd)[[1]] +(dd$ZAlist[[1]] %*% dd$strucList[[1]] %*% (t(ranef(dd, type="uncorrelated")[[1]]))[1:14])[,1] # (don't forget the t())
  p4 <- fixef(dd)[[1]]  +(get_ZALMatrix(dd) %*% t(ranef(dd, type="uncorrelated")[[1]])[1:14])[,1] 
  p5 <- fixef(dd)[[1]]  +(ss$ZAlist[[1]] %*% ss$strucList[[1]] %*% (t(ranef(ss, type="uncorrelated")[[1]]))[1:14])[,1] 
  p6 <- fixef(dd)[[1]]  +(get_ZALMatrix(ss) %*% t(ranef(ss, type="uncorrelated")[[1]])[1:14])[,1] 
  (crit <- diff(range(p1-p2,p1-p3,p1-p4,p1-p5,p1-p6)))
  testthat::test_that(paste0("ranef corrMatrix(mv()...): criterion was ",signif(crit,4)," >1e-8"),
                      testthat::expect_true(crit<1e-8) )
  
  # Reasons for setting all latent d's =1 in the fitted object:
  # predict use a cov matrix deduced from compactcovmat (not modified by .post_process_v_h_LMatrices), and the result of .calc_invL_coeffs(),
  # which uses the $v_h and he latent $d (equally unaffected by .post_process_v_h_LMatrices()
  # => when $d !=1 the attributes of the strucList were unaffected and used, but the main matrices in strucList were not commensurate with the $v_h
  # => get_ZALMatrix as far as being equivalent to ss$ZAlist[[1]] %*% ss$strucList[[1]] was not commensurate with the $v_h. Illustration:
  # => code using the ZAL from get_ZALmatrix had to use d's !1 from the latent_d_list to provide a true ZAL of the $v_h
  # => some code pre v3.8.34 may have be incorrect (and some code using .get_LMatrix too)
  
}
