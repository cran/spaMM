\name{plot_effects}
\alias{plot_effects}
\alias{pdep_effects}
\title{
Partial-dependence effects and plots
}
\description{

  
\code{pdep_effects} evaluates the effect of a given fixed-effect variable, as (by default, the average of) predicted values on the response scale, over the empirical distribution of all other fixed-effect variables in the data, and of inferred random effects. This can be seen as the result of an experiment where specific treatments (given values of the focal variable) are applied over all conditions defined by the other fixed effects and by the inferred random effects. Thus, apparent dependencies induced by associations between predictor variables are avoided (see Friedman, 2001, from which the name \dQuote{partial dependence plot} is taken; or Hastie et al., 2009, Section 10.13.2). This also avoids biases of possible alternative ways of plotting effects. In particular, such biases occur if the response link is not identity, and if averaging is performed on the linear-predictor scale or when other variables are set to some conventional value other than its average. 

\code{pdep_effects} also compute intervals of the type defined by its \code{intervals} argument (by default, prediction intervals). By default, it returns a data frame of average values of point predictiosn and interval bounds for each value of the focal variable, but it can also return lists of all predictions.
  
\code{plot_effects} calls \code{pdep_effects} and produces a simple plot (using only base graphic functions) of its results, including prediction bands representing the two average one-sided widths of intervals.  If added to the plot, the raw data may appear to depart from the partial-dependence predictions, since the data are a priori affected by the associations between variables which the predictions free themselves from.  
}
\usage{
pdep_effects(object, focal_var, newdata = object$data, length.out = 20, 
             levels = NULL, intervals = "predVar", indiv = FALSE, ...)
plot_effects(object, focal_var, newdata = object$data, effects = NULL,
             xlab = focal_var, ylab = NULL, rgb.args = col2rgb("blue"), 
             add = FALSE,  ylim=NULL, ...)
}
\arguments{
  \item{object}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{focal_var}{
   Character string: the name of the predictor variable whose effect is to be represented
}
  \item{newdata}{
   If non-NULL, a data frame passed to \code{predict.HLfit}, whose documentation should be consulted for further details.
}
  \item{effects}{
    If non-NULL, a data frame to substitute to the one produced by default by \code{pdep_effects}.
  }
  \item{xlab}{
   If non-NULL, a character string: X-axis label for the plot.
}
  \item{ylab}{
   If non-NULL, a character string: Y-axis label for the plot.
}
\item{ylim}{
   The \code{plot}'s \code{ylim} argument. Default is based on the (0.025,0.975) quantiles of the response.
}
  \item{rgb.args}{
   Color control arguments, in the format produced by \code{\link[grDevices]{col2rgb}}.
}
  \item{add}{
    Boolean: whether to add graphic elements of a previous plot produced by \code{plot_effects}
}
  \item{length.out}{
    Numeric: for a numeric predictor variable, the number of values at which predictions are evaluated.
}
  \item{levels}{
   If non-NULL, a character vector: for a factor  predictor variable, the levels for which which predictions are evaluated.
 }
  \item{intervals}{
   Passed to \code{predict.HLfit}, whose documentation should be consulted for further details.
}
  \item{indiv}{
   Boolean: whether to return all predictions given the values of other predictors in the \code{newdata}, or only their means.
}
  \item{\ldots}{
   Further arguments passed by \code{plot_effects} to \code{pdep_effects}, or by \code{pdep_effects} to \code{predict.HLfit}.
}

}
\value{
For \code{pdep_effects}, a nested list, or a data frame storing values of the \code{focal_var}, average point predictions \code{pointp} and bounds \code{low} and \code{up} of intervals, depending on the \code{indiv} argument. When \code{indiv} is \code{TRUE}, each sublist contains vectors for \code{pointp}, \code{low} and \code{up}.  

For \code{plot_effects}, the same value, returned invisibly.
}
\references{
J.H. Friedman (2001). Greedy Function Approximation: A Gradient Boosting
Machine. Annals of Statistics 29(5):1189-1232.

J. Friedman, T. Hastie and R. Tibshirani (2009) The Elements of Statistical Learning, 2nd ed. Springer.
}
\examples{
data("scotlip")
hlcor <- HLCor(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
           adjMatrix=Nmatrix,family=poisson(),data=scotlip) 
plot_effects(hlcor,focal_var="prop.ag",ylim=c(0,max(scotlip$cases)))  
points(cases~prop.ag, data=scotlip, col="blue",pch=20)
}
