.spaMM.data <- new.env(parent = emptyenv())
.spaMM.data$options <- list(MESSAGES.FULL.STACK=TRUE,
                            sparsity_threshold=0.05,
                            matrix_method="def_sXaug_EigenDense_QR_scaled", 
                            Matrix_method="def_sXaug_Matrix_QRP_scaled", # _EigenSparse_QR_ is slower !
                            ## possible values: matches to def_sXaug_
                            LevenbergM=FALSE, 
                            USEEIGEN=TRUE,
                            maxLambda=1e10,
                            example_maxtime=0.8,
                            covEstmethod="makeCovEst1",
                            COMP_maxn=1e4,
                            ff_threshold=1e07, ## ! this affects tryn in OKsmooth::rhullByEI !
                            wRegularization=FALSE,
                            wDEVEL=FALSE,
                            wDEVEL2=FALSE,
                            QRmethod=NULL, ## For user-provided values. The code does not and should not change this.
                            spaMM_glm_conv_crit=list(max=-Inf)
)

