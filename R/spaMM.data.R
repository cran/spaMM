## .spaMM.data could be declared as a list() but would not be kosher (unLockBinding requested to modify options)
.spaMM.data <- new.env(parent = emptyenv())
#.spaMM.data$Constants <- list(Version = NA)
.spaMM.data$options <- list(TRACE.UNLINK=FALSE,
                            MESSAGES.FULL.STACK=TRUE,
                            LevenbergM=TRUE, ## not much used
                            INIT.HLFITNAME=NA,
                            USEEIGEN=TRUE,
                            maxLambda=1e10,
                            example_maxtime=1,
                            covEstmethod="makeCovEst1",
                            COMP_maxn=1e4,
                            ff_threshold=1e07, ## ! this affects tryn in OKsmooth::rhullByEI !
                            wRegularization=FALSE,
                            wDEVEL=FALSE,
                            QRmethod=NULL ## For user-provided values. The code does not andshould not change this.
)

