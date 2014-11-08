## .spaMM.data could be declared as a list() but would not be kosher (unLockBinding requested to modify options)
.spaMM.data <- new.env(parent = emptyenv())
.spaMM.data$Constants <- list(Version = NA)
.spaMM.data$options <- list(RHOMAX=100000,NUMAX=50,TRACE.UNLINK=FALSE,MESSAGES.FULL.STACK=TRUE,
                 INIT.HLFITNAME=NA,USEEIGEN=TRUE,maxLambda=1e10,USE_SPARSE_QR=FALSE)

