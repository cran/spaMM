.onAttach <-
function (lib, pkg) {
#   packageStartupMessage( spam.version$version.string," is loaded.",
#       "\nType 'help( Spam)' or 'demo( spam)' for a short introduction ",
#       "\nand overview of this package.",
#       "\nHelp for individual functions is also obtained by ",
#       "adding the\nsuffix '.spam' to the function name, e.g. 'help( chol.spam)'.")
  unlockBinding(".SpaMM", asNamespace("spaMM")) ## required for spaMM.options to be able to change values
}
.SpaMM <-
structure(list(RHOMAX = 1e+05, NUMAX = 50, TRACE.UNLINK = FALSE, 
    MESSAGES.FULL.STACK = TRUE, INIT.HLFITNAME = NA), .Names = c("RHOMAX", 
"NUMAX", "TRACE.UNLINK", "MESSAGES.FULL.STACK", "INIT.HLFITNAME"
))
