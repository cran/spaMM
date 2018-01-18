\name{multinomial}
\alias{multinomial}
\alias{binomialize}
\alias{multi}
\alias{fitted.HLfitlist}
\alias{logLik.HLfitlist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Analyzing multinomial data}
\description{
  These functions facilitate the conversion and analysis of multinomial data as as series of nested binomial data.
  The main function is \code{multi}, to be used in the \code{family} argument of the fitting  functions.
  It calls \code{binomialize}, which can be called directly to check how the data are converted to nested binomial data.
  The \code{fitted.HLfitlist} method of the \code{fitted} generic function returns a matrix of fitted multinomial probabilities.
  The \code{logLik.HLfitlist} method of the \code{logLik} generic function returns a log-likelihood for the joint fits.
}
\usage{
multi(binResponse=c("npos","nneg"),binfamily=binomial(),input="types",...)
binomialize(data,responses,sortedTypes=NULL,binResponse=c("npos","nneg"),
             depth=Inf,input="types")
\method{fitted}{HLfitlist}(object,...)
\method{logLik}{HLfitlist}(object,which,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
   The data frame to be analyzed.
  }
  \item{object}{
    A list of binomial fits returned by a multinomial analysis
  } 
  \item{responses}{
     column names of the data, such that \code{<data>[,<responses>]} contain the multinomial response data, as levels of factor variables.   
}
  \item{sortedTypes}{
     Names of multinomial types, i.e. levels of the multinomial response factors. Their order determines  which types are taken first to define the nested binomial samples. By default, the most common types are considered first.   
}
  \item{binResponse}{
  The names to be given to the number of \dQuote{success} and \dQuote{failures} in the binomial response. 
}
  \item{depth}{
  The maximum number of nested binomial responses to be generated from the multinomial data.
  }
  \item{binfamily}{
  The family applied to each binomial response.
  }
  \item{input}{
  If \code{input="types"}, then the \code{responses} columns must contain factor levels of the binomial response.
  If \code{input="counts"}, then the \code{responses} columns must contain counts of different factor levels, and the column names are the types.
  }
  \item{which}{Which element of the \code{APHLs} list to return. The default depends on the fitting method.In particular, if it was REML or one of its variants, the function returns the log restricted likelihood (exact or approximated).}
  \item{\dots}{
    Other arguments passed from or to other functions.
  } 
}
\details{
  A multinomial response, say counts 17, 13, 25, 8, 3, 1 for types \code{type1} to \code{type6} can be represented as a series of nested binomials
  e.g. \code{type1} against others (17 vs 50) then among these 50 others, \code{type2} versus others (13 vs 37), etc.
  The \code{binomialize} function generates such a representation. By default the representation considers types in decreasing order of the number of positives, i.e. first \code{type3} against others (25 vs 42), then \code{type1} against others within these 42, etc. It stops if it has reached \code{depth} nested binomial responses. This can be modified by the \code{sortedTypes} argument, e.g. \code{sortedTypes=c("type6","type4","type2")}. 
\code{binomialize} returns a list of data frames which can be directly provided as a \code{data} argument for the fitting functions, with binomial response.

Alternatively, one can provide the multinomial response data frame, which will be internally converted to nested binomial data if the \code{family} argument is a call to \code{multinomial} (see examples).    

For mixed models, the multinomial data can be fitted to a model with the same correlation parameters, and either the same or different variances of random effects, for all binomial responses. Which analysis is performed depends on the \code{init.corrHLfit} argument (see \code{\link{corrHLfit}} and the Examples).   
}
\value{
\code{binomialize} returns a list of data frames appropriate for analysis as binomial response. Each data frame contains the original one plus
Two columns named according to \code{binResponse}. 
\code{multi} returns a list.
}
\examples{
## An example considering pseudo-data at one diploid locus for 50 individuals 
set.seed(123)
genecopy1 <- sample(4,size=50,prob=c(1/2,1/4,1/8,1/8),replace=TRUE)
genecopy2 <- sample(4,size=50,prob=c(1/2,1/4,1/8,1/8),replace=TRUE)
alleles <- c("122","124","126","128")
genotypes <- data.frame(type1=alleles[genecopy1],type2=alleles[genecopy2])
## Columns "type1","type2" each contains an allele type => input is "types" (the default)
datalist <- binomialize(genotypes,responses=c("type1","type2"))

## two equivalent fits:
f1 <- HLfit(cbind(npos,nneg)~1,data=datalist, family=binomial())
f2 <- HLfit(cbind(npos,nneg)~1,data=genotypes, family=multi(responses=c("type1","type2")))
fitted(f2)

## distinct fits for spatial data
\dontrun{
genoInSpace <- data.frame(type1=alleles[genecopy1],type2=alleles[genecopy2],x=runif(50),y=runif(50))
## Fitting distinct variances of random effects for each binomial response
corrHLfit(cbind(npos,nneg)~1+Matern(1|x+y),data=genoInSpace, 
          family=multi(responses=c("type1","type2")),
          ranFix=list(rho=1,nu=0.5))
## Fitting the same variance for all binomial responses           
corrHLfit(cbind(npos,nneg)~1+Matern(1|x+y),data=genoInSpace, 
          family=multi(responses=c("type1","type2")),
          ranFix=list(rho=1,nu=0.5),init.corrHLfit=list(lambda=1))
}
}
\keyword{ family }
\keyword{ model }
\keyword{ manip }
