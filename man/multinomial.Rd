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
  The main interface is the \code{multi} \dQuote{family}, to be used in the \code{family} argument of the fitting functions.
  Fits using it call \code{binomialize}, which can be called directly to check how the data are converted to nested binomial data, and to use these data directly. 
  The \code{fitted.HLfitlist} method of the \code{fitted} generic function returns a matrix of fitted multinomial probabilities.
  The \code{logLik.HLfitlist} method of the \code{logLik} generic function returns a log-likelihood for the joint fits.
}
\usage{
multi(binResponse=c("npos","nneg"),binfamily=binomial(),input="types",...)
binomialize(data,responses,sortedTypes=NULL,binResponse=c("npos","nneg"),
             depth=Inf,input="types")
\method{fitted}{HLfitlist}(object, version=2L, ...)
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
  \item{version}{
    Integer, for \code{fitted.HLfitlist} (i.e. for multinomial fits using \code{\link{multi}}); \code{1} will provide the result of past versions up to 3.5.0 (See Value). 
  }
  \item{\dots}{
    Other arguments passed from or to other functions.
  } 
}
\details{
  A multinomial response, say counts 17, 13, 25, 8, 3, 1 for types \code{type1} to \code{type6}, can be represented as a series of nested binomials
  e.g. \code{type1} against others (17 vs 50) then among these 50 others, \code{type2} versus others (13 vs 37), etc.
  The \code{binomialize} function generates such a representation. By default the representation considers types in decreasing order of the number of positives, i.e. first \code{type3} against others (25 vs 42), then \code{type1} against others within these 42, etc. It stops if it has reached \code{depth} nested binomial responses. This can be modified by the \code{sortedTypes} argument, e.g. \code{sortedTypes=c("type6","type4","type2")}. 
\code{binomialize} returns a list of data frames which can be directly provided as a \code{data} argument for the fitting functions, with binomial response.

Alternatively, one can provide the multinomial response data frame, which will be internally converted to nested binomial data if the \code{family} argument is a call to \code{multinomial} (see Examples).    

For mixed models, the multinomial data can be fitted to a model with the same correlation parameters, and either the same or different variances of random effects, for all binomial responses. Which analysis is performed depends on whether the variances are fitted by \dQuote{outer optimization} or by \code{HLfit}'s \dQuote{inner iterative} algorithm, as controlled by the \code{init} or \code{init.corrHLfit} arguments (see Examples). These initial values therefore affect the definition of the model being fitted.      
\code{corrHLfit} will fit different variances by default. Adding an \code{init.corrHLfit} will force estimation of a single variance across models. \code{fitme}'s default optimization strategy is more complex, and has changed and still change over versions. This creates a \bold{back-compatibility issue} where the model to be fitted may change over versions of spaMM. To avoid that, it is strongly advised to use an explicit initial value when fitting a \code{multi} model by \code{fitme}.        

}
\value{
\code{binomialize} returns a list of data frames appropriate for analysis as binomial response. Each data frame contains the original one plus
two columns named according to \code{binResponse}. 

The main fitting functions, when called on a model with \code{family=multi(.)}, return an object of  class \code{HLfitlist}, which is a list with attributes. The list elements are fits of the nested binomial models (objects of class \code{HLfit}). The attributes provide additional information about the overall multinomial model, such as global log-likelihood values and other information properly extracted by the \code{how()} function.  

\code{multi} is a function that returns a list, but users may never need to manipulate this output.

\code{fitted.HLfitlist} returns a matrix. The current default \code{version=2L} provides meaningful fitted values (predicted multinomial frequencies for each response type) even for data rows where the nested binomial fit for a type had no response information remaining. By contrast, the first version provided a matrix with \code{0}s for these row*fit combinations, except for the last column; in many cases this may be confusing.  
}
\examples{
## Adding colour to the famous 'iris' dataset:
iriscol <- iris
set.seed(123) # Simulate colours, then fit colour frequencies:
iriscol$col <- sample(c("yellow", "purple", "blue"),replace = TRUE, 
                      size = nrow(iriscol), prob=c(0.5,0.3,0.2))
colfit <- fitme(cbind(npos,nneg) ~ 1+(1|Species), family=multi(responses="col"), 
                data=iriscol, init=list(lambda=NA)) # note warning if no 'init'...
head(fitted(colfit))

# To only generate the binomial datasets:
binomialize(iriscol,responses="col")

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

if (spaMM.getOption("example_maxtime")>1.7) {

##### Control of lambda estimation over different binomial submodels

genoInSpace <- data.frame(type1=alleles[genecopy1],type2=alleles[genecopy2],
                          x=runif(50),y=runif(50))
method <- "PQL" # for faster exampple

## Fitting distinct variances for all binomial responses:           

multifit <- corrHLfit(cbind(npos,nneg)~1+Matern(1|x+y),data=genoInSpace, 
                      family=multi(responses=c("type1","type2")),
                      ranFix=list(rho=1,nu=0.5), method=method)
length(unique(unlist(lapply(multifit, get_ranPars, which="lambda")))) # 3   

multifit <- fitme(cbind(npos,nneg)~1+Matern(1|x+y),data=genoInSpace, 
                  family=multi(responses=c("type1","type2")),
                  init=list(lambda=NaN), # forcing 'inner' estimation for fitme 
                  fixed=list(rho=1,nu=0.5), method=method)
length(unique(unlist(lapply(multifit, get_ranPars, which="lambda")))) # 3          

## Fitting the same variance for all binomial responses:           

multifit <- fitme(cbind(npos,nneg)~1+Matern(1|x+y),data=genoInSpace, 
                  family=multi(responses=c("type1","type2")),
                  init=list(lambda=NA), # forcing 'outer' estimation for fitme 
                  fixed=list(rho=1,nu=0.5), method=method)
length(unique(unlist(lapply(multifit, get_ranPars, which="lambda")))) # 1          

multifit <- 
  corrHLfit(cbind(npos,nneg)~1+Matern(1|x+y),data=genoInSpace, 
            family=multi(responses=c("type1","type2")),
            init.corrHLfit=list(lambda=1), # forcing 'outer' estimation for corrHLfit 
            ranFix=list(rho=1,nu=0.5), method=method)
length(unique(unlist(lapply(multifit, get_ranPars, which="lambda")))) # 1          
}
}
\keyword{ family }
\keyword{ model }
\keyword{ manip }
