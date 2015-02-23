spaMM.options <- function(...) {
  if (nargs() == 0) return(.spaMM.data$options)
  current <- .spaMM.data$options
  temp <- list(...)
  if (length(temp) == 1 && is.null(names(temp))) {
    arg <- temp[[1]]
    switch(mode(arg),
           list = temp <- arg,
           character = return(.spaMM.data$options[arg]),  ## return here for eg ... = "NUMAX"
           stop("invalid argument: ", sQuote(arg)))
  }
  if (length(temp) == 0) return(current)
  n <- names(temp)
  if (is.null(n)) stop("options must be given by name")
  current[n] <- temp
  .spaMM.data$options <- current
  ## : al palce dans la version 'list' il y avait
  #   if (sys.parent() == 0) {
  #     env <- asNamespace("spaMM")
  #   } else env <- parent.frame()
  #   assign(".SpaMM", current, envir = env) 
  invisible(current)
}

spaMM.getOption <- function (x) {spaMM.options(x)[[1]]}


#`.SpaMM` <- list(RHOMAX=100000,NUMAX=50,TRACE.UNLINK=FALSE,MESSAGES.FULL.STACK=TRUE,
#                 INIT.HLFITNAME=NA,USEEIGEN=TRUE,maxLambda=1e10)
## large rho is not a problem
## large nu is a problem the more so when rho is small (=> 'small scaled distance gaussian')
# lme did not manage to go beyond nu=17.48 in case Vn phiFix...

".onAttach" <- function (lib, pkg) {
  version <- utils::packageDescription("spaMM")$Version
  packageStartupMessage("spaMM (version ", version, 
                          ## not sure this will always work and makes sense only for devel version :
                          # ", packaged ", utils::packageDescription("spaMM")$Packaged,
                        ") is loaded.", 
    "\nType 'help(spaMM)' for a short introduction,\nand news(packages='spaMM') for news.")
  #unlockBinding(".SpaMM", asNamespace("spaMM")) ## required when a .SpaMM list was used instead of an envir
}


".onLoad" <- function (lib, pkg) {
  .spaMM.data$Constants$Version <- utils::packageDescription("spaMM")$Version
}  

largeLambdaMessages <- function() {
  message("A too high lambda may indicate a very poorly fitting fixed-effect part of the model.")
  message("To control the maximum lambda, use e.g. 'spaMM.options(maxLambda=1e06)'.")
  message("It may also indicate convergence issues, possibly improved by altering the initial value through the 'init.HLfit' argument.")
}

