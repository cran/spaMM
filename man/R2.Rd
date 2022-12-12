\name{pseudoR2}
\alias{pseudoR2}
\alias{LR2R2}
\title{
Pseudo R-squared
}
\description{
Generalization of R-squared based on likelihood ratios, called pseudo-R2 below, and variously attributed to Cragg & Uhler (1970), Cox & Snell (1989), Magee (1990) and some other authors (see comments in the References section). The null model used in the definition of R2 can be modified by the user.
}
\usage{
pseudoR2(fitobject, nullform = . ~ 1, R2fun = LR2R2, rescale=FALSE, verbose=TRUE)
}
\arguments{
  \item{fitobject}{
  The fitted model object, obtained as the return value of a \pkg{spaMM} fitting function.
}
  \item{nullform}{
  Main-response formula for the null model. The default value (including only an intercept) represents the traditional choice in R2 computation for linear models. Alternative formulas (including, e.g., random effects) can be specified using either the \code{\link{update.formula}} syntax (e.g., with a \code{'.'} on the right hand side; note that \pkg{spaMM}'s updating conventions differ from those implemented by \code{stats::update.formula}, see \code{\link{update.HLfit}}), or a full formula (which may be a safer syntax).
}
  \item{R2fun}{
  The backend function computing R2 given the fitted and null model. The default implements the pseudo-R2. For linear models, it reduces to the canonical R2 and the value adjusted as in \code{summary.lm} is also returned.
}
\item{rescale}{Boolean or formula, controlling whether and how to rescale R2 so that its maximum possible value is 1 (often considered for discrete-response models). If a formula, it should specify the model with maximal R2. If TRUE, rescaling is performed in a way meaningful only for binary logistic regression (see  Examples for how this is implemented).}
  \item{verbose}{
  Boolean; whether to display various informations about the procedure (most notably, to warn about some potential problem in applying the default procedure to \code{fitobject}).
}
}
\details{
None of the R2-like computations I am aware of helps in addressing, for the general class of models handled by \pkg{spaMM}, a well-defined inference task (comparable to, say, formally testing goodness of fit, or measuring accuracy of prediction of new data as considered for AIC). This problem has been well-known (e.g., Magee, 1990), and the canonical R2 itself for linear models is not devoid of weaknesses (e.g., \url{https://stats.stackexchange.com/questions/13314/is-r2-useful-or-dangerous}). Ultimately the main reason for computing R2 may be to deal with requests by misguided reviewers. Given these problems, strong statements about the properties that R2 should have are difficult to follow (and this includes the claim that it should always have maximum value 1), and no attempt has been made here to implement the wide diversity of R2-like descriptors discussed in the literature. The \code{LR2R2} backend function implements the pseudo-R2, chosen on the basis that this is the simplest general method that makes at least as much sense as any other computation I have seen, and implementation of rescaling by maximal R2 is minimal (the examples explain some of its details). Adaptation of the R2 definition for mixed-effect models can be performed by including some random effect(s) in the null model, using the  \code{nullform} argument.
}
\value{
As returned by the function specified by argument \code{R2fun}. The default function returns a numeric vector of length 2 for linear models and a single value otherwise.
}
\references{
Cox, D.R., Snell, E.J. (1989). The analysis of binary data (2nd ed.). Chapman and Hall.\cr
Often cited in this context, but they barely mention the topic, in an exercise p. 208-209.

Pseudo-R2 is known to go back at least to\cr 
Cragg, J. G., & Uhler, R. S. (1970). The demand for automobiles. The Canadian Journal of Economics, 3(3), 386. \doi{10.2307/133656}\cr
where they already discussed its rescaling by a maximum value, in the context of binary regression.

Magee, L. (1990) R2 Measures based on Wald and likelihood ratio joint significance tests. The American Statistician, 44, 250-253. \doi{10.1080/00031305.1990.10475731}\cr
also often cited for the pseudo-R2, this paper reformulates some related descriptors and concisely reviews earlier literature. 

Nagelkerke, N.J.D. (1991) A note on a general definition of the coefficient of determination. Biometrika, Vol. 78, No. 3. (Sep., 1991), pp. 691-692. \doi{10.1093/biomet/78.3.691}\cr
details the properties of pseudo-R2 (including the way it \dQuote{partition} variation). Argues emphatically for its rescaling, for which it is often cited.
}
\examples{

#### Pseudo-R2 *is* R2 for linear models:
#
# lmfit <- lm(sr ~ pop15+pop75+dpi+ddpi , data = LifeCycleSavings)
# summary(lmfit) # Multiple R-squared = 0.3385, adjusted = 0.2797
#
spfit <- fitme(sr ~ pop15+pop75+dpi+ddpi , data = LifeCycleSavings)
pseudoR2(spfit)  # consistent with summary(lmfit)

#### Toy example of pseudo-R2 for binary data
#
set.seed(123)
toydf <- data.frame(x=seq(50), y=sample(0:1,50,TRUE))
#
##   Binary logistic regression:
#
binlog <- fitme(y~x, data=toydf, family=binomial())
(blR2 <- pseudoR2(binlog)) # quite low, despite the model being correct
#
##   Rescaling by 'maximum possible' R2 for binary logistic regression:
#
pseudoR2(binlog, rescale=TRUE)
#
#    which is acheived by silently computing the maximum possible R2 value 
#    by the following brutal but effective way:
#
perfbinlog <- fitme(y~I(y), data=toydf, family=binomial())
(maxblR2 <- pseudoR2(perfbinlog)) # = 0.7397...
#
# (this 'maximum possible' value would be modified if the null model were modified).
#
blR2/maxblR2       # again, rescaled value 
#
##   Same by more general syntax:
#
pseudoR2(binlog, rescale=y~I(y)) 
}
