cat("\ntest of sparse Precision method:")

data("scotlip")
fit1 <- fitme(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
               fixed=list(lambda=0.1), 
               adjMatrix=Nmatrix,family=poisson(),data=scotlip)
spaMM.options(sparse_precision=TRUE)
fit2 <- fitme(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
              fixed=list(lambda=0.1), 
              adjMatrix=Nmatrix,family=poisson(),data=scotlip)
spaMM.options(sparse_precision=FALSE)
expect_equal(logLik(fit1),logLik(fit2),tolerance=1e-05)

if (spaMM.getOption("example_maxtime")>323) {
  ## example suggested by Jeroen van den Ochtend jeroenvdochtend@gmail.com Jeroen.vandenochtend@business.uzh.ch
  library(data.table)
  library(igraph)
  
  rsample <- function(N=100,month_max=10,seed) {
    if (is.integer(seed)) set.seed(seed)
    dt <- data.table(ID=factor(1:N))
    dt$months <- sample(1:month_max,N,replace=T) ## nombres de lignes qui vont être créées pour chaque ligne originelle de dt
    dt$GENDER <- sample(c("MALE","FEMALE"),N,replace=T)
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
  
  # to trace convergence 
  trace(HLfit_body,print=FALSE,tracer=quote(try(cat(paste(names(ranFix),"=",signif(unlist(ranFix),6))," "))),
        exit=quote(cat(paste(" ",signif(res$APHLs$p_v,6),"\n"))))
  trace(def_AUGI0_ZX_sparsePrecision,print=FALSE,tracer=quote(cat(".")))
  
  spaMM.options(sparse_precision=NULL) 
  spaMM.options(sparse_precision=FALSE) 
  ## IMPORTANT :
  spaMM.options(sparse_precision=TRUE) ## but this should be fitme()'s default given the data. 
  
  set.seed(123)
  blob <- rsample(N=1000,seed=NULL)
  system.time({
    MLfit.Frailty <- fitme(BUY ~ factor(month) + AGE + GENDER + X1*X2 + adjacency(1|ID),
                           #fixed=list(rho=0.01,lambda=0.1),
                           data=blob$data,family = binomial(link = cloglog),method = "ML",
                           control.HLfit=list(LevenbergM=TRUE),
                           adjMatrix=blob$adjMatrix
    )
  }) ## 323s
  untrace(HLfit_body)
  untrace(def_AUGI0_ZX_sparsePrecision)
}