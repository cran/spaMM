cat("\ntest of sparse Precision method:")

data("scotlip")
fit1 <- fitme(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
               fixed=list(lambda=0.1), 
               adjMatrix=Nmatrix,family=poisson(),data=scotlip)
expectedMethod <- "sXaug_EigenDense_QRP_Chol_scaled" ## bc data too small to switch to sparse
if (interactive()) {
  if (! (expectedMethod %in% fit1$MME_method)) {
    message(paste('! ("',expectedMethod,'" %in% fit1$MME_method): was a non-default option selected?'))
  }
} else testthat::expect_true(expectedMethod %in% fit1$MME_method) 
oldop <- spaMM.options(sparse_precision=TRUE)
fit2 <- fitme(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
              fixed=list(lambda=0.1), 
              adjMatrix=Nmatrix,family=poisson(),data=scotlip)
testthat::expect_true("AUGI0_ZX_sparsePrecision" %in% fit2$MME_method)
spaMM.options(oldop)
testthat::expect_equal(logLik(fit1),logLik(fit2),tolerance=1e-05)

if (spaMM.getOption("example_maxtime")>535) { ## approx 116 by IRLS and 417 by LevM (v2.1.107) 
  ## fixme This is sparse_precision LevM: improve efficiency? 
  ## example suggested by Jeroen van den Ochtend jeroenvdochtend@gmail.com Jeroen.vandenochtend@business.uzh.ch
  library(data.table)
  library(igraph)
  
  rsample <- function(N=100, ## size of implied adjacency matrix
                      month_max=10,seed) {
    if (is.integer(seed)) set.seed(seed)
    dt <- data.table(ID=factor(1:N))
    dt$months <- sample(1:month_max,N,replace=T) ## nombres de lignes qui vont être créées pour chaque ligne originelle de dt
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
    g <- add_vertices(g,sum(!unique(dt$ID) %in% V(g)),name=unique(dt[!dt$ID %in% V(g),list(ID)]))
    Network <- as_adjacency_matrix(g,sparse = TRUE,type="both")
    return(list(data=dt,adjMatrix=Network))
  }

  ## it's no use to try sparse_precision=FALSE bc bc the augmented sigma matrix is huge
  ###################### spaMM.options(sparse_precision=FALSE) 
  ## make sure that the wrong valu is not set:
  oldop <- spaMM.options(sparse_precision=TRUE) ## but this should be fitme()'s default given the data. 
  
  set.seed(123)
  blob <- rsample(N=1000,seed=NULL)
  system.time({
    IRLS.Frailty <- fitme(BUY ~ factor(month) + AGE + GENDER + X1*X2 + adjacency(1|ID),
                          data=blob$data,family = binomial(link = cloglog),method = "ML",
                          control.HLfit=list(LevenbergM=FALSE), ## inhibits default for binary data 
                          verbose=c(TRACE=TRUE), # to trace convergence 
                          adjMatrix=blob$adjMatrix
    )
  }) ##  116.56
  system.time({
    LevM.Frailty <- fitme(BUY ~ factor(month) + AGE + GENDER + X1*X2 + adjacency(1|ID),
                          data=blob$data,family = binomial(link = cloglog),method = "ML",
                          verbose=c(TRACE=TRUE), # to trace convergence 
                          adjMatrix=blob$adjMatrix
    )
  }) ## 393.45 (v2.1.108)
  spaMM.options(oldop)
  untrace(HLfit_body)
  untrace(def_AUGI0_ZX_sparsePrecision)
}