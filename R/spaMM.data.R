.spaMM.data <- new.env(parent = emptyenv())
.spaMM.data$options <- list(MESSAGES.FULL.STACK=TRUE,
                            sparsity_threshold=0.05,
                            separation_max=1000,
                            matrix_method="def_sXaug_EigenDense_QRP_Chol_scaled", 
                            #Matrix_method="def_sXaug_Matrix_cholP_scaled", 
                            #Matrix_method= "def_sXaug_Matrix_QRP_scaled", 
                            Matrix_method= "def_sXaug_Matrix_QRP_CHM_scaled", 
                            ## possible values: matches to def_sXaug_
                            LevenbergM=NULL, 
                            USEEIGEN=TRUE,
                            maxLambda=1e10,
                            example_maxtime=0.7,
                            covEstmethod=".makeCovEst1",
                            COMP_maxn=1e4,
                            ff_threshold=Inf, ## ! this affects tryn in blackbox::rhullByEI ! F I X M E remove when new blackbox on CRAN
                            wRegularization=FALSE,
                            # sparse_precision=FALSE,
                            wDEVEL2=FALSE, ## m'a servi a tester le remplacement de ranFix par une liste un peu diff ? laissé en chantier
                            QRmethod=NULL, ## For user-provided values. The code does not and should not change this.
                            spaMM_glm_conv_crit=list(max=-Inf),
                            spaMM_tol=list(Xtol_rel=1e-5, Xtol_abs=1e-6, Ftol_LM=1e-1), 
                            optimizer1D="optimize", 
                            optimizer="nloptr", ## "nloptr" vs "bobyqa" or "L-BFGS-B"
                            optimize_tol=.Machine$double.eps^0.25, ## default tol value for optimize
                            bobyqa=list(),
                            nloptr=list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=5e-6, xtol_abs= 1e-7,
                                        maxeval=1000,print_level=0),
                            CMP_asympto_cond=quote((#nu<1 && ## affects test-COMPoisson
                              pow_lam_nu > 10/nu) || 1+pow_lam_nu+6*sqrt(pow_lam_nu/(nu)) > .spaMM.data$options$COMP_maxn),
                            rankMethod="qr", ## private
                            rankTolerance=quote(max(1e-7,.Machine$double.eps*10*ncol(X.pv))) ## private
                            # , sparse_X=NULL## private
)

