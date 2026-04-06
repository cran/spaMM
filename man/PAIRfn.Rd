\name{dyad}
\alias{dyad}
\alias{mm}
\alias{mmfn}
\alias{PAIRfn}
\title{
Correlated random effects in dyadic interactions
}
\description{
\code{mmfn} is a function to be used for specifying random effects in asymmetric multi-membership interactions. This feature is experimental and has only been tested for some dyadic interactions. 
%A pair of effects is not the same thing as a factor with two levels, but base R procedures constructing model matrices from factors are used internally, which might have some visible side effects. 

The simplest case is that of dyadic interactions between one focal individual and its partner in each dyad. The models declared in this way include correlation parameter(s) for the correlation of effects expressed by a given individual in different possible roles (such as focal vs. partner). 

At the same time, correlation of effects expressed among individuals in the same role can also be specified. For example, the random effects may represent genetic influences on focal and partner effects. The correlation between effects expressed in the focal role may depend on the relatedness between individuals, and likewise for effects expressed in the partner role. This can be specified by a random-effect term 
\code{corrMatrix( <LHS> | <RHS> )} where:

\code{ * }the \code{LHS} expression includes an expression of the form 
\code{mm(<id1>, <id2>)}, \code{<id1>} and \code{<id2>} being two \code{factor}s in the \code{data}; and the \code{RHS} expression is of the form
\code{mmfn(<id1>, <id2>)}. A dummy variable can be used instead of the \code{mm(<id1>, <id2>)} expression (see Examples), the factors to be used being determined from the RHS expression;\cr

\code{ * }The \code{LHS} expression is interpreted as a random-coefficient specification, meaning that distinct random effects \eqn{u_{\code{<id1>}(i)}},\eqn{v_{\code{<id2>}(i)}},... affect the \eqn{i}th response, and that correlations between the effects are fitted by default, as for other random-coefficient terms. The  elements \eqn{u_{\code{<id1>}(i)}},\eqn{v_{\code{<id2>}(j)}},... are correlated when \eqn{\code{<id1>}(i)} and \eqn{\code{<id2>}(j)} are identical, rather than when \eqn{i=j}. This is illustrated in the Examples, where the two factors represents roles of individuals as focal individuals and as their mothers, so  distinct but correlated random effects \eqn{u_{k}},\eqn{v_{k}}, affecting different levels of the response (\eqn{i\neq j}), are assigned to a same individual \eqn{k} in these two roles;

\code{ * }a \code{corrMatrix} specifies the correlations of random effects among individuals in either role (correlations as focal, or correlations as partner). The whole \code{corrMatrix(  <LHS> | <RHS> )} expression thus describes a \code{\link{composite-ranef}} controlled by two correlation models, one for within-role correlations, the other for among-role correlations.

If the \code{LHS} is a more complex expression than only the \code{mm(...)} term (as in the Example with interaction with sex), different random effects vectors are assigned to each term implied by the LHS (as for random-coefficient terms specified by other syntaxes), and the \code{corrMatrix} specifies the correlations between the effects within each vector; 

\code{ * }the \code{RHS} of the form \code{mmfn(...)} specifies the factors containing the levels to be matched to rows and columns of the \code{corrMatrix}, in the same way as other forms of RHS serve to identify rows and columns of the correlation matrix of correlated random effects.  

\code{PAIRfn} is an alias for \code{mmfn} (its name suggests that it handles only dyadic interactions, but the \code{mmfn} aims to be more general). 

The syntax using a dummy variable assumes that this variable is present in the input data, as the numeric constant \code{1} (see Examples). By default, this dummy variable is named \code{PAIR}, though this can be changed by the \code{VAR} argument of \code{mmfn} or \code{PAIRfn}. 
}
\usage{
mmfn(..., VAR = "PAIR", only.vars=TRUE)
PAIRfn(..., VAR = "PAIR", only.vars=TRUE)

## formula term:
# corrMatrix(  < LHS using mm(<.>, <.>) >| PAIRfn(<.>, <.>) ) 
# corrMatrix(  < LHS using dummy variable >| PAIRfn(<.>, <.>) ) 
}
\arguments{
  \item{\ldots}{
   \code{factor}s present in the \code{data} provided for the fit, identifying the individuals in the dyadic or multi-membership interaction. 
}
\item{VAR}{character: name of dummy variable possibly used in LHS of random-effect term.  
}
  \item{only.vars}{For programming purposes, not documented.}
}
%\details{
%In principle, one can define alternatives du \code{PAIRfn} that implement other random effect structures.
%Such a function must have a \code{only.vars} set to \code{TRUE} by default. 
%For \code{only.vars=TRUE}, it must return a matrix of the variables ... 
%}
\seealso{\code{\link{diallel}} for (anti-)symmetric dyadic interactions.}
\examples{
if (spaMM.getOption("example_maxtime")>1) {

data("Gryphon")
fitme(BWT ~ 1+ corrMatrix(mm(ID,mother) | PAIRfn(ID,mother)), 
      data=Gryphon_df, corrMatrix=Gryphon_A)

# : same as using a PAIR dummy variable as follows:
#
# Gryphon_df$PAIR <- 1 # Do not use values other than 1.
# fitme(BWT ~ 1+ corrMatrix(PAIR | PAIRfn(ID,mother)), 
#       data=Gryphon_df, corrMatrix=Gryphon_A)

if (spaMM.getOption("example_maxtime")>15) {
  # Distinct random effects for each role, and each sex of the focal: 
  #
  fitme(BWT ~ 1+ corrMatrix(mm(ID,mother)*sex | PAIRfn(ID,mother)), 
        data=Gryphon_df, corrMatrix=Gryphon_A)
  #
  # This is controlled as other rancom-coefficient terms. E.g.,
  #   fixed=list(ranCoefs=list("1"=c(NA,0,0,0,NA,0,0,NA,0,NA)))
  # could be used to fit only the variances.

}
}

}
