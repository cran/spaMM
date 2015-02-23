## .spaMM.data could be declared as a list() but would not be kosher (unLockBinding requested to modify options)
.spaMM.data <- new.env(parent = emptyenv())
.spaMM.data$Constants <- list(Version = NA)
.spaMM.data$options <- list(RHOMAX=100000,
                            NUMAX=50,
                            TRACE.UNLINK=FALSE,
                            MESSAGES.FULL.STACK=TRUE,
                            LevenbergM=TRUE, ## not much used
                            INIT.HLFITNAME=NA,
                            USEEIGEN=TRUE,
                            maxLambda=1e10,
                            example_maxtime=1,
                            covEstmethod="makeCovEst1"                            
                            ## default QRmethod is NULL
                            ##QRmethod="lmwithQ_denseZAL" ## meaningful: "Matrix::qr" "lmwithQ_sparseZAL" "lmwithQ_denseZAL" "lmwithSparseQ"
                            )

