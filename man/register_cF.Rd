\name{register_cF}
\alias{register_cF}
\title{
Declare corrFamily constructor for use in formula
}
\description{
This function registers the name of a new corrFamily constructor so that it can be used as the keyword of a random effect in a formula (as in \code{y ~ 1 + ARp()})
}
\usage{
register_cF(corrFamilies = NULL, reset = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{corrFamilies}{
NULL, or character vector of names of corrFamily constructors. 
}
  \item{reset}{
Boolean. Set it to \code{TRUE} in order to reset the list of registered constructors to the \pkg{spaMM} built-in default, before registering the ones specified by \code{corrFamilies}.
}
}
\value{
No value; operates through side-effects on internal variables. 
}
\examples{
ts <- data.frame(lh=lh,time=seq(48)) ## using 'lh' data from 'stats' package

myARp <- ARp                   #  defines 'new' corrFamily from built-in one

# Now, this would not yet work:

# fitme(lh ~ 1 + myARp(1|time), data=ts, method="REML")

# but thisworks if we first register "myARp"

register_cF("myARp")           #  registers it

fitme(lh ~ 1 + myARp(1|time), data=ts, method="REML")
#
# same as 
#
fitme(lh ~ 1 + corrFamily(1|time), data=ts, method="REML",
                  covStruct=list(corrFamily=myARp()))
# 
# showing it's possible not to register myARp, 
# although this has limitations (see Details in help("corrFamily")).

## Specifying arguments of the corrFamily constructor:

fitme(lh ~ 1 + myARp(1|time, p=3), data=ts, method="REML")
#                  
# same as
#
fitme(lh ~ 1 + corrFamily(1|time), data=ts, method="REML",
                  covStruct=list(corrFamily=ARp(p=3)))
                  
register_cF(reset=TRUE) # Tidy things before leaving.                  
                  
}
