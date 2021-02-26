.spaMM.data <- new.env(parent = emptyenv())
.spaMM.data$options <- list(
  F_I_X_M_E=FALSE,
  store_data_as_mf=FALSE, # Surely many issues.
  Rcpp_crossprod=TRUE, # integer with usual bool interp., and >1: .crossprod() prints types when .Rcpp_crossprod() not called; >2: always prints types;
  TRY_R_new=TRUE, # 
  update_CHM=TRUE, # measurable benefits only if Cholesky(., perm=TRUE)
  perm_G=TRUE, 
  perm_Q=NULL,  
  use_ZA_L=TRUE, # NULL may act as TRUE when augZxy_cond=TRUE
  bind_ZAL=TRUE, # set it to FALSE to use ZAXlist beyond spprec 
  sparsity_threshold=0.05,
  spprec_threshold=50, # ohio small by correlation algo, large by spprec: threshold is n>=140 has crit 'near' 62 (varying betw replicates). 
  separation_max=10,
  sep_solver="glpk",
  spprec_method="def_AUGI0_ZX_sparsePrecision", 
  matrix_method="def_sXaug_EigenDense_QRP_Chol_scaled", 
  Matrix_method= "def_sXaug_Matrix_QRP_CHM_scaled", 
  EigenDense_QRP_method=".lmwithQR", # .lmwithQR seems fast cf bootstrap
  use_spprec_QR=FALSE, # TRUE visibly slows several of the long tests (incl fitar1) 
  #Matrix_method="def_sXaug_Matrix_cholP_scaled", 
  #Matrix_method= "def_sXaug_Matrix_QRP_scaled", 
  ## possible values: matches to def_sXaug_
  #
  spaMM_glm_conv_silent=FALSE,
  LevenbergM=NULL, 
  LevM_HL11_method=list(b_step="v_b", # [1]= "v_b" or "b", or "v_in_b",
                        rescue_thr_null=c(rescue=TRUE,strictv=0L,V_IN_B=2L,re_V_IN_B=Inf), ## affects LevM.negbin test
                        rescue_thr_AoN=c(rescue=TRUE,strictv=0L,V_IN_B=2L,re_V_IN_B=Inf) # binomial All-or-None case
                        ), 
  # cf also spaMM_tol
  use_G_dG=TRUE, # meaningful only for spprec
  spprec_LevM_D="1", # form of the perturbation of Md2hdv2 in .calc_G_dG() (alternatives are "colSums" or "rowSums")
  #
  USEEIGEN=TRUE, # Whether to use the Eigen C++ library for some matrix computations. The source code should be consulted for further information. 
  X_scaling=TRUE,
  maxLambda=1e10,
  regul_lev_lambda=1e-8,
  ############## augZXy stuff (see also ranCoefs settings)
  allow_augZXy=NULL, ## interpreted as TRUE if phiScal (=>not phiFix) before further conditions are applied, and FALSE otherwise 
  # allow_augZXy=2L forces augZXy usage with non-constant prior weights, if other conditions for its usage are satisfied.
  augZXy_solver=c("chol","EigenQR"), # "chol", "QR" (currently = "EigenQR"), "EigenQR" (dense or sparse), or "qr" (=base::qr)
  augZXy_fitfn=".HLfit_body_augZXy", # safe version, no specific singularity, but no refinement beyond augmentation by y
  check_alt_augZXy=FALSE, ## private, effective only if alternative augZXy fitfn is set to TRUE
  ##############
  optimizer1D="optimize", 
  optimizer=".safe_opt", ## "nloptr", ##  "bobyqa", ## "L-BFGS-B",
  optim_boot=".safe_opt",
  #
  fpot_tol= - Inf, #newer criterion for IRLS NOT LevM # cf implicit ftol_abs for coherence
  optimize_tol=.Machine$double.eps^0.25, ## default tol value for optimize
  bobyqa_margin=1e-14, # horrors may happen if bobyqa's init is precisely at a boundary
  bobyqa_rhofn= function(lower,upper) min(155.2475, 0.1*min(upper-lower)), # min() to avoid infiniy here; 155.2475 is diff(spaMM:::.dispFn(c(1e-6,1e6)))/10
  # For a long time this was effectively min(0.95.2475, 0.1*min(upper-lower)), where 0.95 is mysterious value from minqa doc.                                    
  #                                   0.1 seems the right balance between time (increase with lower value) and effective minimization (cf 21 selected examples from test-nloptr) 
  bobyqa=list(), 
  nloptr=list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=5e-6, print_level=0), # nloptr options only control the termination criteria, not the step sizes. Nothing like rhobeg
  maxeval=quote(10^(3+(log(length(initvec))-log(5))/log(4))), # nloptr; *modified for bobyqa (which recommends > 10 * npar^2)
  maxeval_corr=1, # devel: for easy control of maxeval in .safe_opt()
  xtol_abs_factors=c(rcLam=5e-7,rcCor=5e-6,others=5e-11,abs=1e-7), # nloptr! # laxer rcLam=5e-5 strongly affect tests spherical transfo (+minor effect in test-poly)
  xtol_abs=quote(.xtol_abs_fn(LowUp,  rC_transf = rC_transf)), # nloptr; zero's for max precision (?)
  ############## ranCoefs settings: (see also xtol_abs_factors)
  optim_inner=".safe_opt",
  recheck_at_bound=FALSE, # control of .safe_opt()
  rC_unbounded=FALSE, # unbounded parametrization as in PinheiroB96
  rC_transf="chol", 
  rC_transf_inner="chol", 
  rC_transf_fac=1, # should be made dependent on link __F I X M E__ # >1 to reduce effect in canon sp of steps in transf sp
  tol_ranCoefs_inner=c(inner=TRUE,
                       cholmax=Inf, # in "chol" case # so it's always Inf... 
                       lo_lam=1e-6,up_lam=1e6,corr=1e-12,tol=1e-5 # in "sph" case ## corr must be < 1 ! 
                       ,regul=1e-09
                       ), # controls .makeCovEst1() bounds (and .calc_latentL(.)$lambda_est)
  tol_ranCoefs_outer=c(inner=FALSE,lo_lam=1e-4,up_lam=1e5,corr=1e-4, # controls fitme() bounds in "sph" case, not in default rc_transf="chol"; corr must be < 1 ! 
                       regul=1e-08), ## regul used also for chol. 
  max_bound_ranCoefs=1, # adjustment in transformed coordinates; 1 means no adjustment# previously 0.99999, affects HLfit(Reaction ~ Days + (Days|Subject), data = sleepstudy)
  regul_ranCoefs=c(10*.Machine$double.eps), ## used to avoid zero eigenvalue after correction in .smooth_regul()
                   #,covregulcor=0L), # OL inhibits correction by .CovRegulCor(); otherwise covregulcor give the correction in that fn) )
  # .calc_latentL() control: (triangular facto. is ALWAYS used for spprec)
  use_tri_for_augZXy=FALSE, # *!*spprec; fitme -> .HLfit_body_augZXy[_W]; Seems marginally faster with no drawback.
  use_tri_for_makeCovEst=FALSE, # *!*spprec; TRUE WAS required for acceptable result in HLfit3 rC_transf_inner="sph" test! Affects numerical precision of calc_latentL() in .makeCovEst1() [ultimately using sXaug, not augZXy method].
  replace_design_u=TRUE, # the strucList[[<ranCoef>]] is not $design_u if $d were not all 1
  #•
  invL_threshold=1e6, ## for devel code .HLfit_body_augZXy_invL(); compare to prod(sqrt(lambda)) ## F I X_invL but test "set.seed(666)" fails for invL_threshold>100
  ###############
  example_maxtime=0.7,
  COMP_maxn=1e4,
  bin_mu_tol=.Machine$double.eps, # was 1e12 for a long time
  sanitize_eta=c(gauslog=20,otherlog=20), # otherlog value affects difficult Gamma(log) fits (cf Leucadendron_hard.R)
  Gamma_min_y = 1e-10, ## for warnings in .preprocess(), and automatic correction in simulate() -> .r_resid_var(); .calc_dispGammaGLM() has indep, and much less strict, correction
  QRmethod=NULL, ## For user-provided values. The code does not and should not change this.
  #
  # Most useful tests of LevM controls: initial updates in tests_private/optim_LevM.R, and the tough newy in Leucadendron_hard.R
  #
  spaMM_tol=list(Xtol_rel=1e-5, # .calc_etaGLMblob (beta coeffs); and HLfit_body (dispersion params)
                 Xtol_abs=1e-6, # HLfit_body (dispersion params)
                 # not LevMar HL11 (but still LevMar PQL)
                 d_relV_b_tol=1e-05, # for solve_IRLS; critical termination test for outer IRLS loop; strong impact on timings of test_all.
                 #                     Lax d_relV_b_tol not good for large data as small d(v,b) has large p_v impact;
                 # LevMar HL11: (cf optim_LevM.R for optimization of this)
                 v_pot_tol_rescue=1e-06, # critical Q/P trade-off (e.g. test-nloptr 362)
                 v_pot_tol_noresc=1e-04, # 
                 b_pot_tol=1e-05, # wrap_do_damped.; sensitive and 1e-04 was previously bad but seems mostly OK now (v3.0.35) though still suspect in a COMPoisson difficult case
                 d_relV_b_tol_LM=1e-05, # critical termination test for outer !HL11! IRLS loop; 1e-04 => too good lik in some lines of LevM.Frailty test.
                 d_relV_tol=1e-05,# for v_h IRLS; critical "not_moving" test; optim_LevM -> 1e-04 is slower and less accurate!
                 loose_fac=5, # need to test whether is useful or not
                 loose_resc=0.1,
                 dampings_env_v=list(v=list("b"=1e-7,"v"=1e-7,"v_b"=1e-7,"V_IN_B"=1e-7,"v_in_b"=1e-7, 
                                            "strict_v|b"=1e-7, "b_&_v_in_b"=1e-7, "b_from_v_b"=1e-7))
                 # residues from devel, NO IMPACT:
                 #,Ftol_v_in_b=1e-6, # NO IMPACT on current fits : only for diagnosis of not_moving in IRLS_v_h fns
                 #Ftol_LM=1e-5 # some not_moving diagnostics 
  ), 
  #
  CMP_asympto_cond=quote((#nu<1 && ## affects test-COMPoisson
    pow_lam_nu > 10/nu) || 1+pow_lam_nu+6*sqrt(pow_lam_nu/(nu)) > .spaMM.data$options$COMP_maxn),
  rankMethod="qr", ## private
  rankTolerance=quote(max(1e-7,.Machine$double.eps*10*ncol(X.pv))), ## private, used  by preprocess
  qrTolerance=1e-10, ## private, used by select qr() calls for predVar computation
  # , sparse_X=NULL## private
  uGeo_levels_type="mf", # same type to be used by .calc_AMatrix_IMRF() and .calc_Zmatrix() for IMRFs. Explicit names, useful for debugging. ALternative is ".ULI" (which is faster?)
  INLA_A=TRUE,
  #
  stylefns=list(v_in_loop=crayon::green, 
                v_in_last=crayon::green$underline, # final output of v_h .do_damped_WLS_v_in_b
                rescue=crayon::red,
                strictv=crayon::blue,
                vloop=crayon::cyan,
                v_out_last=crayon::cyan$underline, # final output of v_h .do_damped_WLS_outer; also also bracketing each .solve_v_h_IRLS loop for v_h( tentative beta(damping) ) 
                # colors tell what the numbers are for: grad of objective for v, versus grad of objective for beta (or joint beta,v) 
                betaloop=crayon::yellow, # also bracketing the damped_WLS loop for new beta when  which_LevMar_step=="b_&_v_in_b"
                betalast=crayon::yellow$underline,
                # 
                # not for LevM:
                hardwarn=crayon::bold),
  H_scale_regul=1e-4,
  ## only to avoid warnings when using spaMM.options()
  nb_cores=NULL,
  barstyle=quote(if(interactive()) {3L} else {0L}),
  #        add control=list(fix_predVar=NA) in predict() calls in the following calls? Probably not worth the mess.
  fix_predVar=list("NA"="MSL|bboptim|isoscape|isofit|calibfit|optimthroughSmooth|spaMM_rhullByEI|sampleByResp",
                   "TRUE"=NULL,"FALSE"=NULL), 
  thr_backsolve=0L # for devel testing of .backsolve(); 0L means that .Rcpp_backsolve may be called irrespective of matrix dimension
)

