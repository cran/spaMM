\name{stripHLfit}
\alias{stripHLfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reduce the size of fitted objects 
}
\description{
Large matrices and other memory-expensive objects may be stored in a fit object. This function removes them in order to reduce the size of the object, particularly when stored on disk. In principle, the removed objects can be regenerated automatically when needed (e.g., for a predict()). 
}
\usage{
stripHLfit(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
The result of a fit (an object of class \code{HLfit}).
}
  \item{\dots}{
  Further arguments, not currently used.
}
}
\value{
The input fit objects with some elements removed.
}
\note{The effect may change without notice between versions as the efficiency of the operation is highly sensitive to implementation details.}
\examples{
\dontrun{
## rather unconvincing example : quantitative effect is small.

# measure size of saved object:
saveSize <- function (object,...) {
    tf <- tempfile(fileext = ".RData")
    on.exit(unlink(tf))
    save(object, file = tf,...)
    file.size(tf)
  }
data("Loaloa")  
lfit <- fitme(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
                   +Matern(1|longitude+latitude), method="HL(0,1)",
              data=Loaloa, family=binomial(), fixed=list(nu=0.5,rho=1,lambda=0.5))
saveSize(lfit)                 
pfit <- predict(lfit,newdata=Loaloa,variances=list(cov=TRUE)) # increases size!
saveSize(lfit)
lfit <- stripHLfit(lfit)
saveSize(lfit)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
