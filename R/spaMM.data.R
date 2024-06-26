.spaMM.data <- new.env(parent = emptyenv())
.spaMM.data$options <- list(
  F_I_X_M_E=FALSE,
  obsInfo=TRUE, # !! Don't forget to inactivate .obsInfo_warn() when the default is TRUE !!
  store_data_as_mf=FALSE, # Surely many issues.
  Rcpp_crossprod=TRUE, # integer with usual bool interp., and >1: .crossprod() prints types when .Rcpp_crossprod() not called; >2: always prints types;
  update_CHM=TRUE, # measurable benefits only if Cholesky(., perm=TRUE)
  perm_G=TRUE, 
  perm_Q=NULL,  
  use_ZA_L=TRUE, # NULL may act as TRUE when augZxy_cond=TRUE
  bind_ZAL=TRUE, # set it to FALSE to use ZAXlist beyond spprec 
  sparsity_threshold=0.05,
  spprec_threshold=50, # ohio small by correlation algo, large by spprec: threshold is n>=140 has crit 'near' 62 (varying betw replicates). 
                       # 2023/08: such high threshold still relevant (see Nmatrix example in test-spaMM)
  separation_max=10,
  sep_solver="glpk",
  HLfit_body="HLfit_body",
  spprec_method="def_AUGI0_ZX_spprec", 
  matrix_method="def_sXaug_EigenDense_QRP_Chol_scaled", # handling negative weights  
  Matrix_method= "def_sXaug_Matrix_QRP_CHM_scaled", 
  #Hobs_Matrix_method= "def_sXaug_Matrix_QRP_CHM_scaled", # may have a patch for handling negative weights, but not exactly.  
  Hobs_Matrix_method= "def_sXaug_Matrix_CHM_H_scaled", # handling negative weights  
  force_LLF_CHM_QRP=FALSE, # if set to TRUE, LevM no longer uses the exact Hessian with signs (CHM_H methods) when this Hessian is SPD.
  force_LLM_nosigns_CHM_H=FALSE, # set it to TRUE to test CHM_H methods when there are no $signs (quite slow)
  LLgeneric=TRUE,
  EigenDense_QRP_method=".lmwithQR", # .lmwithQR seems fast cf bootstrap
  use_spprec_QR=FALSE, # TRUE visibly slows several of the long tests (incl fitar1) 
  presolve_cond=quote(ncol(BLOB$R_scaled)>20000L), # may be slightly faster and more memory efficient for huge data (nested-Matern 40000L subset) 
  #Matrix_method="sXaug_Matrix_CHM_H_scaled", 
  #Matrix_method= "def_sXaug_Matrix_QRP_scaled", 
  ## possible values: matches to def_sXaug_
  #
  LevenbergM=NULL, 
  LevM_HL11_method=list(b_step="v_b", # [1]= "v_b" or "b", or "v_in_b",
                        rescue_thr_null=c(rescue=TRUE,strictv=0L,V_IN_B=2L,re_V_IN_B=Inf), ## affects LevM.negbin test
                        rescue_thr_AoN=c(rescue=TRUE,strictv=0L,V_IN_B=2L,re_V_IN_B=Inf) # binomial All-or-None case
                        ), 
  # cf also spaMM_tol
  use_G_dG=TRUE, # meaningful only for spprec
  spprec_LevM_D="1", # form of the perturbation of Md2hdv2 in .calc_G_dG() (alternatives are "colSums" or "rowSums")
  #
  Utri_chol_method="RcppEigen", # 
  USEEIGEN=TRUE, # Whether to use the Eigen C++ library for some matrix computations. The source code should be consulted for further information. 
  X_scaling=TRUE,
  minLambda=1e-8,
  maxLambda=1e8,
  regul_lev_lambda=1e-16, #has now seemed OK since at least v3.12.13 (July 2022)
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
  fpot_tol= - Inf, #newer criterion for IRLS NOT LevM; coherent with the following (implicit) ftol_abs:
  # ftol_abs= - Inf, # distinct from nloptr potential control; could be made explicit whereas is implicit NULL as-is. 
  optimize_tol=.Machine$double.eps^0.25, ## default tol value for optimize
  bobyqa_margin=1e-14, # horrors may happen if bobyqa's init is precisely at a boundary
  bobyqa_rhofn= function(lower,upper) min(155.2475, 0.1*min(upper-lower)), # min() to avoid infiniy here; 155.2475 is diff(spaMM:::.dispFn(c(1e-6,1e6)))/10
  # For a long time this was effectively min(0.95.2475, 0.1*min(upper-lower)), where 0.95 is mysterious value from minqa doc.                                    
  #                                   0.1 seems the right balance between time (increase with lower value) and effective minimization (cf 21 selected examples from test-nloptr) 
  bobyqa=list(), 
  nlminb=list(), 
  # default value for nloptr() 'opts':
  nloptr=list(algorithm="NLOPT_LN_BOBYQA",
              xtol_rel=4e-6, # cf comment on compMatfit
              print_level=0), # nloptr options only control the termination criteria, not the step sizes. Nothing like rhobeg
  ## further control of nloptr 'opts' (but not suitable input for 'opts'):
  xtol_abs_factors=c(abs=1e-8, # That's the general one when next ones are not used. # cf comments in .xtol_abs_fn()
                     rcLam=5e-7,rcCor=5e-6,others=5e-11), # ___F I X M E____ all only when there are ranCoefs...
  # laxer rcLam=5e-5 strongly affect tests spherical transfo (+minor effect in test-poly)
  xtol_abs=quote(.xtol_abs_fn(LowUp,  rC_transf = rC_transf)), # nloptr; zero's for max precision (?)
  maxeval=quote(as.integer(10^(3+(log(length(initvec))-log(5))/log(4)))), # nloptr; *modified for bobyqa (which recommends > 10 * npar^2)
  maxeval_corr=1, # devel: for easy control of maxeval in .safe_opt()
  ## 
  ############## ranCoefs settings: (see also xtol_abs_factors)
  optim_inner=".safe_opt",
  recheck_at_bound=FALSE, # control of .safe_opt()
  reuse_bobyqa=FALSE, # whether to reuse bobyqa when it reached its maxit at a boundary
  rC_unbounded=FALSE, # unbounded parametrization as in PinheiroB96
  rC_transf="chol", 
  rC_transf_inner="chol", 
  rC_transf_fac=1, # should be made dependent on link _F I X M E__ # >1 to reduce effect in canon sp of steps in transf sp
  tol_ranCoefs_inner=c(inner=TRUE,
                       cholmax=Inf, # in "chol" case # so it's always Inf... 
                       lo_lam=1e-6,up_lam=1e6,corr=1e-12,tol=1e-5 # in "sph" case ## corr must be < 1 ! 
                       ,regul=1e-09
                       ), # controls .makeCovEst1() bounds (and .calc_latentL(.)$lambda_est)
  tol_ranCoefs_outer=c(inner=FALSE,lo_lam=1e-4,up_lam=1e5,corr=1e-4, # controls fitme() bounds in "sph" case, not in default rC_transf="chol"; corr must be < 1 ! 
                       regul=1e-08), ## regul used also for chol. 
  max_bound_ranCoefs=1, # adjustment in transformed coordinates; 1 means no adjustment# previously 0.99999, affects HLfit(Reaction ~ Days + (Days|Subject), data = sleepstudy)
  regul_ranCoefs=c(10*.Machine$double.eps), ## used to avoid zero eigenvalue after correction in .smooth_regul()
                   #,covregulcor=0L), # OL inhibits correction by .CovRegulCor(); otherwise covregulcor give the correction in that fn) )
  # .calc_latentL() control: (triangular facto. is ALWAYS used for spprec)
  use_tri_for_augZXy=FALSE, # *NOT-spprec*; fitme -> .HLfit_body_augZXy[_W]; Seems marginally faster with no drawback.
  use_tri_for_makeCovEst=FALSE, # *NOT-spprec*; TRUE WAS required for acceptable result in HLfit3 rC_transf_inner="sph" test! Affects numerical precision of calc_latentL() in .makeCovEst1() [ultimately using sXaug, not augZXy method].
  #•
  invL_threshold=1e6, ## for devel code .HLfit_body_augZXy_invL(); compare to prod(sqrt(lambda)) ## F I X_invL but test "set.seed(666)" fails for invL_threshold>100
  ## Family-speific numerical controls
  # Gaunt et al:
  # Based on the numerical results
  # of Section 4, we consider that a safe rule of thumb for obtaining accurate approximations
  # using the asymptotic approximation (1.4), ******with the first three terms,**** is for both [lambda] >= 1.5
  # and [pow_lam_nu] >= 1.5 to hold (the absolute error was always less than 0.5% in our tests),
  #
  # but this is only for the limited set of  values in the Table.
  CMP_asympto_cond=function(pow_lam_nu, nu, lambda){pow_lam_nu > 3 && lambda>5},
  # CMP_asympto_cond=quote((#nu<1 && ## affects test-COMPoisson
  #   pow_lam_nu > 10/nu) || 1+pow_lam_nu+6*sqrt(pow_lam_nu/nu) > .spaMM.data$options$COMP_maxn),
  COMP_maxn=1e4,
  sanitize_eta=c(gauslog=20,COMPlog=16,otherlog=20), # otherlog value affects difficult Gamma(log) fits (cf Leucadendron_hard.R)
  Gamma_min_y = 1e-10, ## for warnings in .preprocess(), and automatic correction in simulate() -> .r_resid_var(); .calc_dispGammaGLM() has indep, and much less strict, correction
  beta_min_y = 1e-8, ##  for warnings in .preprocess(), and automatic correction in simulate() -> .r_resid_var();  ___F I X M E____ value is quick ad hoc fix
  ###############
  example_maxtime=0.7,
  bin_mu_tol=.Machine$double.eps, # was 1e12 for a long time
  QRmethod=NULL, ## For user-provided values. The code does not and should not change this. Cf control.HLfit$algebra too
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
                                            "strict_v|b"=1e-7, "b_&_v_in_b"=1e-7, "b_from_v_b"=1e-7)),
                 # residues from devel, NO IMPACT:
                 #,Ftol_v_in_b=1e-6, # NO IMPACT on current fits : only for diagnosis of not_moving in IRLS_v_h fns
                 #Ftol_LM=1e-5 # some not_moving diagnostics 
                 logL_tol=5e-5 # heuristic but effective (precision ~2*logL_tol )
  ), 
  rankMethod="qr", ## private
  rankTolerance=quote(max(1e-7,.Machine$double.eps*10*ncol(X.pv))), ## private, used  by preprocess
  qrTolerance=1e-10, ## private, used by select qr() calls for predVar computation
  # , sparse_X=NULL## private
  uGeo_levels_type="data_order", # same type to be used by .calc_AMatrix_IMRF() and .calc_Zmatrix() for IMRFs. Explicit names, useful for debugging.
  #
  stylefns=list(v_in_loop=crayon::green, 
                v_in_last=crayon::green$underline, # final output of v_h .do_damped_WLS_v_in_b
                rescue=crayon::red,
                strictv=crayon::blue,
                vloop=crayon::cyan,
                v_out_last=crayon::cyan$underline, # seen a lot for v steps... # old comment : final output of v_h .do_damped_WLS_outer; also also bracketing each .solve_v_h_IRLS loop for v_h( tentative beta(damping) ) 
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
  wrap_parallel="dopar",
  #        add control=list(fix_predVar=NA) in predict() calls in the following calls? Probably not worth the mess.
  fix_predVar=list("NA"="MSL|bboptim|isoscape|isofit|calibfit|optimthroughSmooth|spaMM_rhullByEI|sampleByResp",
                   "TRUE"=NULL,"FALSE"=NULL), 
  tr_beta=FALSE, # whether to optim on transformed scale in speculative outer-optimization of beta
  #thr_backsolve=0L # for devel testing of .backsolve(); 0L means that .Rcpp_backsolve may be called irrespective of matrix dimension
  xLM_conv_silent=FALSE,
  xLM_conv_crit=list(max=-Inf),
  use_terms_info_attr=FALSE, # controls data vs model frame updating post fit.
  NbThreads=1L,
  diagnose_conv=2000L,
  #
  # devl
  .betaFn=function(v) {sign(v)*log1p(abs(v))},
  .betaInv=function(v) {sign(v)*(exp(abs(v))-1)},
  n_names2expr=FALSE # I fail to reproduce the problem that motivated the devel of ..n_names2expr()
)

.spaMM.data$keywords <- new.env(parent = emptyenv())
# What is a "special ranef"? Originally, presumably those without a corrFamily constructor. 
# But more importantly now: $special_ranefs controls corr_info$corr_types and therefore the ad-hoc code for specific corr_type's:
# .assign_corr_types_families() will call special functions such as IMRF() to evaluate corr_info$corr_families, even in fitmv -> .preprocess() case
#  While for non-special families, the corr_info$corr_families will be stubs in this case.
# .canonizeRanPars() assumes that special types have a $canonize() function, for other types it handles transformations using its own code...

.spaMM.data$keywords$special_ranefs <- c("adjacency", "Matern", "Cauchy", "AR1", "corrMatrix", "IMRF", "corrFamily") 
.spaMM.data$keywords$all_cF <- .spaMM.data$keywords$built_in_cF <- c("ARp", "ARMA", "diallel", "ranGCA", "MaternIMRFa", "antisym") 
.spaMM.data$keywords$all_ranefs <- .spaMM.data$keywords$built_in_ranefs <- 
  unlist(c(.spaMM.data$keywords$special_ranefs,.spaMM.data$keywords$built_in_cF),recursive = FALSE, use.names = FALSE) 
.spaMM.data$keywords$all_keywords <- .spaMM.data$keywords$built_in_keywords <- c(.spaMM.data$keywords$built_in_ranefs,
                                                                                 "multIMRF", "mv")
.spaMM.data$keywords$user_defined <- character(0L)
# For comparing names and strings:
# as.name if fastest (~8e-7 s versus as.vector(, "character") ~ 3e-6 and paste() ~ 4e-6)
# so to conpare single string to single name, convert the string.
# But comparing to a sequence of names requires a 'for' loop bc even unlist( < keywords as names > ) is a list. Hence no %in% test possible
#  Then the for loop costs ~1e-5s compared to the as.vector() ~ 3e-6 => We avoid using sequences of names

register_cF <- function(corrFamilies=NULL, reset=FALSE) {
  keywords <- .spaMM.data$keywords
  if (reset) {
    keywords$user_defined <- character(0L)
    keywords$all_cF <- keywords$built_in_cF    
    keywords$all_ranefs <- keywords$built_in_ranefs
    keywords$all_keywords <- keywords$built_in_keywords
  }
  if ( ! is.null(corrFamilies) ) {
    user_defined <- unique(c(keywords$user_defined, corrFamilies))
    if (length(shared <- intersect(keywords$built_in_keywords, user_defined))) {
      stop(paste0("Keyword(s) '", paste(shared,collapse="', '"),"' are already defined in spaMM and cannot be redeclared."))
    } else{ 
      keywords$user_defined <- user_defined
      keywords$all_cF <- unique(c(keywords$all_cF, user_defined))    
      keywords$all_ranefs <- unique(c(keywords$all_ranefs, user_defined))
      keywords$all_keywords <- unique(c(keywords$all_keywords, user_defined))
    }
  }
}

unregister_cF <- function(corrFamilies) {
  keywords <- .spaMM.data$keywords
  corrFamilies <- setdiff(corrFamilies, keywords$built_in_keywords)
  keywords$user_defined <- setdiff(keywords$user_defined, corrFamilies)
  keywords$all_ranefs <- setdiff(keywords$all_ranefs, corrFamilies)
  keywords$all_cF <- setdiff(keywords$all_cF, corrFamilies)
  keywords$all_keywords <- setdiff(keywords$all_keywords, corrFamilies)
}
 
.spaMM.data$class_cache <- new.env(parent = emptyenv())
