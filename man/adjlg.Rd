\name{adjlg}
\alias{adjlg}
\alias{adjlgMat}
\docType{data}
\title{
Simulated data set for testing sparse-precision code
}
\description{
This is used in \code{tests/test-adjacency-long.R}
}
\usage{data("adjlg")}
\format{
  Includes an adjacency matrix \code{adjlgMat}. and a data frame \code{adjlg} with 5474 observations on the following 8 variables.
  \describe{
    \item{\code{ID}}{a factor with levels \code{1} to \code{1000}}
    \item{\code{months}}{a numeric vector}
    \item{\code{GENDER}}{a character vector}
    \item{\code{AGE}}{a numeric vector}
    \item{\code{X1}}{a numeric vector}
    \item{\code{X2}}{a numeric vector}
    \item{\code{month}}{a numeric vector}
    \item{\code{BUY}}{a numeric vector}
  }
}
\source{
The simulation code show in Example was suggested by Jeroen van den Ochtend. 
}
\examples{
data(adjlg)
## See further usage in tests/test-adjacency-long.R
\dontrun{
# as produced by:
  library(data.table) ## Included data produced using version 1.10.4.3
  library(igraph) ## Included data produced using version 1.2.1
  
  rsample <- function(N=100, ## size of implied adjacency matrix
                      month_max=10,seed) {
    if (is.integer(seed)) set.seed(seed)
    dt <- data.table(ID=factor(1:N))
    dt$months <- sample(1:month_max,N,replace=T) ## # of liens for each level of ID
    dt$GENDER <- sample(c("MALE","FEMALE"),N,replace=TRUE)
    dt$AGE <- sample(18:99,N,replace=T)
    dt$X1 <- sample(1000:9900,N,replace=T)
    dt$X2 <-  runif(N)
    
    dt <- dt[, c(.SD, month=data.table(seq(from=1, to=months, by = 1))), by = ID] 
    dt[,BUY := 0]
    dt[month.V1==months,BUY := sample(c(0,1),1),by=ID]
    setnames(dt,"month.V1","month")
    
    #### create adjacency matrix
    Network <- data.table(OUT=sample(dt$ID,N*month_max*4/10))
    Network$IN <- sample(dt$ID,N*month_max*4/10)
    Network <- Network[IN != OUT]
    Network <- unique(Network)
    g <- graph.data.frame(Network,directed=F)
    g <- add_vertices(g,sum(!unique(dt$ID) \%in\% V(g)),name=unique(dt[!dt$ID \%in\% V(g),list(ID)]))
    Network <- as_adjacency_matrix(g,sparse = TRUE,type="both")
    return(list(data=dt,adjMatrix=Network))
  }

  set.seed(123)
  adjlg_sam <- rsample(N=1000,seed=NULL) 
  adjlg <- as.data.frame(adjlg_sam$data)
  adjlgMat <- adjlg_sam$adjMatrix
}
}
\keyword{datasets}
