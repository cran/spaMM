.spaMM.data <- new.env(parent = emptyenv())
.spaMM.data$options <- list(
  Rcpp_crossprod=TRUE, # integer with usual bool interp., and >1: .crossprod() prints types when .Rcpp_crossprod() not called; >2: always prints types;
  rC_transf="chol", 
  rC_transf_inner="chol", 
  TRY_R_new=TRUE, # 
  TRY_update=TRUE, # measurable benefits only if Cholesky(., perm=TRUE) hence Q_CHMfactor do not benefit from it.... (fixme)
  perm_G=TRUE, 
  TRY_ZAX=NULL, # default may be TRUE for augZxy_cond=TRUE
  sparsity_threshold=0.05,
  separation_max=1000,
  spprec_method="def_AUGI0_ZX_sparsePrecision", 
  matrix_method="def_sXaug_EigenDense_QRP_Chol_scaled", 
  Matrix_method= "def_sXaug_Matrix_QRP_CHM_scaled", 
  EigenDense_QRP_method=".lmwithQR", # .lmwithQR seems fast cf bootstrap
  use_spprec_QR=FALSE, # TRUE visibly slows several of the long tests (incl fitar1) 
  #Matrix_method="def_sXaug_Matrix_cholP_scaled", 
  #Matrix_method= "def_sXaug_Matrix_QRP_scaled", 
  ## possible values: matches to def_sXaug_
  #
  LevenbergM=NULL, 
  use_G_dG=TRUE, # good
  spprec_LevM_D="1", # form of the perturbation of Md2hdv2 in .calc_G_dG() (alternatives are "colSums" or "rowSums")
  #
  mat_sqrt_fn="mat_sqrt",
  USEEIGEN=TRUE, 
  lev_by_sparse_Q=20000L, # switch to sparse QR in .leveragesWrap()
  X_scaling=TRUE,
  maxLambda=1e10,
  regul_lev_lambda=1e-8,
  ############## augZXy stuff (see also ranCoefs settings)
  allow_augZXy=NULL, ## interpreted as TRUE if phiScal (=>not phiFix) before further conditions are applied, and FALSE otherwise 
  # allow_augZXy=2L forces augZXy usage with non-constant prior weights, if other conditions for its usage are satisfied.
  augZXy_solver=c("chol","EigenQR"), # "chol", "QR" (currently = "EigenQR"), "EigenQR" (dense or sparse), or "qr" (=base::qr)
  augZXy_fitfn=".HLfit_body_augZXy", # safe version, no specific singularity, but no refinement beyond augmentation by y
  check_alt_augZXy=FALSE, ## private, effective only if alternative augZXy fitfn is set to TRUE
  ############## ranCoefs settings: (see also xtol_abs_factors)
  covEstmethod=".makeCovEst1",
  rC_unbounded=FALSE, # unbounded parametrization as in PinheiroB96
  tol_ranCoefs_inner=c(inner=TRUE,
                       cholmax=Inf, # in "chol" case # so it's always Inf... 
                       lo_lam=1e-6,up_lam=1e6,corr=1e-12,tol=1e-5 # in "sph" case ## corr must be < 1 ! 
                       ), # controls .makeCovEst1() bounds (and .calc_latentL(.)$lambda_est)
  tol_rel_ranCoefs=c(inner=FALSE,lo_lam=1e-4,up_lam=1e5,corr=1e-4), ## corr must be < 1 ! # controls fitme() bounds in "sph" case
  max_corr_ranCoefs=0.99999, ## actually not the max corr but an adjustement in transformed coordinates; affects HLfit(Reaction ~ Days + (Days|Subject), data = sleepstudy)
  regul_ranCoefs=c(0,10*.Machine$double.eps),
  # Lower condnum for inner estimation (bc of hatval computation)? 
  # svd and qr *allow* lower condnum than eigen and Cholesky, but they also allow good fits with low condnum
  condnum_for_latentL=1e30, ## fitme(Clim ~ year + (yearStd|id), resid.model = ...) is an important test here
  condnum_for_latentL_spprec=1e11, ## fitme
  condnum_for_latentL_inner=1e12, ## for "sph" ! Code now avoids dependence on this for "chol" (where it needed 1e14)
  #trychol_for_latentL_inner=FALSE, ## (effective only if ! spprec)
  # condnum_for_latentL_lu=1e11,
  invL_threshold=1e6, ## compare to prod(sqrt(lambda)) ## F I X_invL but test "set.seed(666)" fails for invL_threshold>100
  # .calc_latentL() control:
  use_tri_for_augZXy=FALSE, # fitme; Seems marginally faster with no drawback.
  use_tri_for_makeCovEst=TRUE, # (F I X M E) HLfit; check HLfit3, HLfit6 and sleepstudy test! Affects numerical precision of calc_latentL() in .makeCovEst1() [ultimately using sXaug, not augZXy method].
  ###############
  example_maxtime=0.7,
  COMP_maxn=1e4,
  #ff_threshold=Inf, ## removable with given blackbox> 1.1.25 now on CRAN
  #wRegularization=FALSE,
  QRmethod=NULL, ## For user-provided values. The code does not and should not change this.
  spaMM_tol=list(Xtol_rel=1e-5, Xtol_abs=1e-6, Ftol_LM=1e-1), # Ftol_LM for LevMar code
  optimizer1D="optimize", 
  optimizer="nloptr", ## "nloptr" vs "bobyqa" or "L-BFGS-B"
  #
  Gamma_min_y = 1e-10,
  optimize_tol=.Machine$double.eps^0.25, ## default tol value for optimize
  bobyqa=list(),
  nloptr=list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=5e-6, print_level=0),
  maxeval=quote(10^(3+(log(length(initvec))-log(5))/log(4))), # nloptr; *modified for bobyqa (which recommends > 10 * npar^2)
  xtol_abs_factors=c(rcLam=5e-7,rcCor=5e-6,others=5e-11,abs=1e-7), # nloptr!
  xtol_abs=quote(.xtol_abs_fn(LowUp,  rC_transf = rC_transf)), # nloptr; zero's for max precision (?)
  #
  CMP_asympto_cond=quote((#nu<1 && ## affects test-COMPoisson
    pow_lam_nu > 10/nu) || 1+pow_lam_nu+6*sqrt(pow_lam_nu/(nu)) > .spaMM.data$options$COMP_maxn),
  rankMethod="qr", ## private
  rankTolerance=quote(max(1e-7,.Machine$double.eps*10*ncol(X.pv))), ## private, used  by preprocess
  qrTolerance=1e-10, ## private, used by select qr() calls for predVar computation
  # , sparse_X=NULL## private
  Zcolsbyrows=FALSE, ## private; diferent values in fit and post fit would currently generate wrong results
  levels_type="mf" # explicit names, useful for debugging. ALternative is ".ULI" (which is faster?)
)

