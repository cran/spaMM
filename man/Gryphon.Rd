\name{Gryphon}
\alias{Gryphon}
\alias{Gryphon_A}
\alias{Gryphon_df}
\alias{Gryphon_pedigree}
\docType{data}
\title{
Gryphon data
}
\description{
Loading these data loads three objects describing a mythical 'Gryphon' population used by Wilson et al. to illustrate mixed-effect modelling in quantitative genetics. These objects are a data frame \code{Gryphon_df} containing the model variables, a genetic relatedness matrix \code{Gryphon_A}, and another data frame \code{Gryphon_pedigree} containing pedigree information (which can be used by some packages to reconstruct the relatedness matrix).
}
\usage{data("Gryphon")}
\format{
\code{Gryphon_df} is 
\preformatted{
'data.frame':	1084 obs. of  6 variables:
 $ ID    : int  1029 1299 ...:               individual identifier 
 $ sex   : Factor w/ 2 levels "1","2":       sex, indeed
 $ year  : Factor w/ 34 levels "968","970", ...: birth year
 $ mother: Factor w/ 429 levels "1","2",..:  individual's mother identifier 
 $ BWT   : num  10.77 9.3  ...:              birth weight 
 $ TARSUS: num  24.8 22.5 12 ...:            tarsus length
}
\code{Gryphon_A} is a genetic relatedness matrix, in sparse matrix format, for 1309 individuals.

\code{Gryphon_pedigree} is 
\preformatted{
'data.frame':	1309 obs. of  3 variables:
 $ ID  : int  1306 1304 ...: individual identifier 
 $ Dam : int  NA NA ...:     individual's mother    
 $ Sire: int  NA NA ...:     individual's father
}
}
\references{
Wilson AJ, et al. (2010) An ecologist's guide to the animal model. Journal of Animal Ecology 79(1): 13-26. 
\doi{10.1111/j.1365-2656.2009.01639.x}
}
\examples{
# See help("sparse_precision")
}
\keyword{datasets}
