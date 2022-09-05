\name{plot_effects}
\alias{plot_effects}
\alias{pdep_effects}
\title{
Partial-dependence effects and plots
}
\description{
The following functions evaluate or plot \emph{partial-dependence} effects. 
  
\code{pdep_effects} evaluates the effect of a given fixed-effect variable, as (by default, the average of) predicted values on the response scale, over the empirical distribution of all other fixed-effect variables in the data, and of inferred random effects. This can be seen as the result of an experiment where specific treatments (given values of the focal variable) are applied over all conditions defined by the other fixed effects and by the inferred random effects. Thus, apparent dependencies induced by associations between predictor variables are avoided (see Friedman, 2001, from which the name \dQuote{partial dependence plot} is taken; or Hastie et al., 2009, Section 10.13.2). This also avoids biases of possible alternative ways of plotting effects. In particular, such biases occur if the response link is not identity, and if averaging is performed on the linear-predictor scale or when other variables are set to some conventional value other than its average. 

\code{pdep_effects} also compute intervals of the type defined by its \code{intervals} argument (by default, prediction intervals). By default, it returns a data frame of average values of point predictions and interval bounds for each value of the focal variable, but it can also return lists of all predictions.
  
A plot function is available for numeric or factor predictors: \code{plot_effects} calls \code{pdep_effects} and produces a simple plot (using only base graphic functions) of its results, including prediction bands representing the two average one-sided widths of intervals. The last section of the Examples shows how to obtain more elaborate plots including the same information using \pkg{ggplot2}. 

If added to the plot, the raw data may appear to depart from the partial-dependence predictions, since the data are a priori affected by the associations between variables which the predictions free themselves from. An adapted plot of fit residuals may be then be more useful, and the Examples also show how it can be performed.   
}
\usage{
pdep_effects(object, focal_var, newdata = object$data, length.out = 20, 
             focal_values=NULL, levels = NULL, intervals = "predVar", indiv = FALSE, 
             ...)
plot_effects(object, focal_var, newdata = object$data, focal_values=NULL, 
             effects = NULL, xlab = focal_var, ylab = NULL, 
             rgb.args = col2rgb("blue"), add = FALSE,  ylim=NULL, ...)
}
\arguments{
  \item{object}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{focal_var}{
   Character string: the name of the predictor variable whose effect is to be represented. The variable must be numeric for \code{plot_effects} but not necessarily so for \code{pdep_effects}. 
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
    Integer: for a numeric predictor variable, this controls the number of values at which predictions are evaluated. By default, predictions are made at regular intervals over the range of the predictor variable. If \code{length.out=0}, predictions are made for the actual values of the focal predictor in the data. The default behaviour is also overriden by using \code{focal_values}, in which case predictions are evaluated at the given \code{focal_values} (as if \code{length.out=0}), unless a non-zero \code{length.out} is also specified. In the latter case, predictions are evaluated at regular intervals over the range of \code{focal_values}.
}
  \item{focal_values, levels}{ \code{focal_values} may be used to specify the values of the focal variable at which predictions are evaluated. For factor variables, \code{levels} is an older implementation of this control, and is now redundant.
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

# Impose specific values of a numeric predictor using 'focal_values':
plot_effects(hlcor, focal_var="prop.ag", focal_values=1:5)

### Adding 'partial residuals' [residuals relative to predict(<fit object>),
###  but plotted relative to pdep_effects() predictions]:

# One first needs predictions for actual values of the predictor variable,
# provided by pdep_effects(.,length.out=0L):
#
pdep_points <- pdep_effects(hlcor,focal_var="prop.ag",length.out=0L)

# Rename for easy prediction for each observation, and add the residuals 
# of the actual fit, using the default residuals() i.e. deviance ones: 
#
rownames(pdep_points) <- pdep_points$focal_var
pdep_res <- pdep_points[paste(hlcor$data$prop.ag),"pointp"] + 
              residuals(hlcor)

points(x = hlcor$data$prop.ag, y = pdep_res, col = "red", pch = 20)

\dontrun{ 

## Plotting pdep-effects for different categories, using ggplot.
library(ggplot2)

data("Gryphon")
tmp <- na.omit(Gryphon_df)
spfit <- spaMM::fitme(TARSUS ~ BWT*sex, data = tmp)

tmp$sex <- "1"
pdep_1 <- pdep_effects(spfit,"BWT", newdata=tmp)
tmp$sex <- "2"
pdep_2 <- pdep_effects(spfit,"BWT", newdata=tmp)
pdep_1$sex <- "1" ; pdep_2$sex <- "2"  
pdep <- rbind(pdep_1,pdep_2)

ggplot(pdep,aes(y = pointp , x = focal_var ,col = sex, fill=sex)) + geom_point() +
  geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.3) + xlab("BWT") +
  ylab("TARSUS")


}
}
