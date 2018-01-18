\name{salamander}
\alias{salamander}
\docType{data}
\title{
Salamander mating data
}
\description{
Data from a salamander mating experiment discussed by McCullagh and Nelder (1989, Ch. 14). Twenty males and twenty females from two populations (Rough Butt and Whiteside) were each paired with 6 individuals from their own or from the other population. The experiments were later published by Arnold et al. (1996).
}
\usage{data("salamander")}
\format{
  The data frame includes 360 observations on the following variables:
  \describe{
    \item{Female}{Index of the female;}
    \item{Male}{Index of the male;}
    \item{Mate}{Whether the pair successfully mated or not;}
    \item{TypeF}{Population of origin of female;}
    \item{TypeM}{Population of origin of male;}
    \item{Cross}{Interaction term between \code{TypeF} and \code{TypeM};}
    \item{Season}{A factor with levels \code{Summer} and \code{Fall};}
    \item{Experiment}{Index of experiment}
  }
}
\source{
The data frame was borrowed from the \code{HGLMMM} package (Molas and Lesaffre, 2011), version 0.1.2.
}
\references{
Arnold, S.J., Verrell, P.A., and Tilley S.G. (1996) The evolution of asymmetry in sexual isolation: a model and a test case. Evolution 50, 1024-1033.

McCullagh, P. and Nelder, J.A. (1989). Generalized Linear Models, 2nd edition. London: Chapman & Hall.

Molas, M., Lesaffre, E. (2011) Hierarchical Generalized Linear Models: The R Package HGLMMM. Journal of Statistical Software 39, 1-20.
}
\examples{
data("salamander")

\dontrun{ 
% in test-salamander.R
if (spaMM.getOption("example_maxtime")>0.7) {
  HLfit(cbind(Mate,1-Mate)~TypeF+TypeM+TypeF*TypeM+(1|Female)+(1|Male),
    family=binomial(),data=salamander,HLmethod="ML",control.HLfit=list(LevenbergM=FALSE))
}
}
if (spaMM.getOption("example_maxtime")>0.7) {
 fitme(cbind(Mate,1-Mate)~TypeF+TypeM+TypeF*TypeM+(1|Female)+(1|Male),
    family=binomial(),data=salamander,control.HLfit=list(LevenbergM=FALSE))
}
}
\keyword{datasets}
