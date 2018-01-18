\name{welding}
\alias{welding}
\docType{data}
\title{
 Welding data set
}
\description{
The data give the results of an unreplicated experiment for factors affecting welding
quality conducted by the National Railway Corporation of Japan (Taguchi and Wu, 1980, cited in Smyth et al., 2001). It is a toy example 
for heterocedastic models and is also suitable for illustrating fit of overparameterized models.}
\usage{data("welding")}
\format{
  The data frame includes 16 observations on 10 variables:
  \describe{
    \item{Strengh}{response variable;}
    \item{...}{nine two-level factors.}
  }
}
\source{
The data were downloaded from http://www.statsci.org/data/general/welding.txt on 2014/08/19 and are consistent with those shown
in table 5 of Bergman and Hynén (1997).
}
\references{
Bergman B, Hynén A (1997) Dispersion effects from
unreplicated designs in the \eqn{2^{k-p}} series.
Technometrics, 39, 191–98.

Smyth GK, Huele AF, Verbyla AP (2001). Exact and approximate REML for heteroscedastic regression. Statistical Modelling 1, 161-175. 

Taguchi G, Wu Y (1980) Introduction to off-line
quality control. Nagoya, Japan: Central Japan
Quality Control Association.
}
\examples{
data("welding")
## toy example from Smyth et al.
HLfit(Strength ~ Drying + Material,resid.model = ~ Material+Preheating ,data=welding)
## toy example of overparameterized model
HLfit(Strength ~ Rods+Thickness*Angle+(1|Rods),resid.model = ~ Rods+Thickness*Angle ,data=welding)
}
\keyword{datasets}
