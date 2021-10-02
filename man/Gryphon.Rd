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

#### Bivariate-response model used as example in Wilson et al. (2010):
# joint modelling of birth weight (BWT) and tarsus length (TARSUS).

# The relatedness matrix is specified as a 'corrMatrix'. The random 
# effect 'corrMatrix(0+mv(1,2)|ID)' then represents genetic effects 
# correlated over traits and individuals (see help("composite-ranef")).
# The ...(0+...) syntax avoids contrasts being used in the design 
# matrix of the random effects, as it would not does make much sense 
# to represent TARSUS as a contrast to BWT. 

# The relatedness matrix will be specified through its inverse,
# using as_precision(), so that spaMM does not have to find out and 
# inform the user that using the inverse is better (as is typically 
# the case for relatedness matrices). But using as_precision() is 
# not required. See help("algebra") for Details.

# The second random effect '(0+mv(1,2)|ID)' represents correlated 
# environmental effects. Since measurements are not repeated within 
# individuals, this effect also absorbs all residual variation. The 
# residual variances 'phi' must then be fixed to some negligible values 
# in order to avoid non-identifiability.

if (spaMM.getOption("example_maxtime")>7) { 
  data("Gryphon")
  gry_prec <- as_precision(Gryphon_A)
  gry_GE <- fitmv(
    submodels=list(BWT ~ 1 + corrMatrix(0+mv(1,2)|ID)+(0+mv(1,2)|ID), 
                  TARSUS ~ 1 + corrMatrix(0+mv(1,2)|ID)+(0+mv(1,2)|ID)), 
    fixed=list(phi=c(1e-6,1e-6)), 
    corrMatrix = gry_prec, 
    data = Gryphon_df, method = "REML")
    
  # Estimates are practically identical to those reported for package 
  # 'asreml' (https://www.vsni.co.uk/software/asreml-r) 
  # according to Supplementary File 3 of Wilson et al., p.7:

  lambda_table <- summary(gry_GE, digits=5,verbose=FALSE)$lambda_table 
  by_spaMM <- na.omit(unlist(lambda_table[,c("Var.","Corr.")]))[1:6]
  by_asreml <- c(3.368449,12.346304,3.849875,17.646017,0.381463,0.401968)
  by_spaMM/by_asreml-1  # relative differences ~ O(1e-4)

}

}
\keyword{datasets}
