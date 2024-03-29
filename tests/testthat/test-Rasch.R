cat(crayon::yellow("Multilevel Rasch Model example by Doran et al. (2007):"))

## example from:
## Harold Doran, Douglas Bates, Paul Bliese, Maritza Dowling
## 	Estimating the Multilevel Rasch Model: With the lme4 Package
## j stat software 2007 DOI:  	10.18637/jss.v020.i02

if( requireNamespace("multilevel", quietly = TRUE)) {
  if (spaMM.getOption("example_maxtime")>7) {
    data("lq2002", package = "multilevel")
    wrk <- lq2002
    for (i in 3:16) wrk[[i]] <- ordered(wrk[[i]])
    for (i in 17:21) wrk[[i]] <- ordered(5 - wrk[[i]])
    lql <- reshape(wrk, varying = list(names(lq2002)[3:21]), v.names = "fivelev",
                   idvar = "subj", timevar = "item", drop = names(lq2002)[c(2,22:27)],
                   direction = "long")
    lql$itype <- with(lql, factor(ifelse(item < 12, "Leadership",
                                         ifelse(item < 15, "Task Sig.", "Hostility"))))
    for (i in c(1,2,4,5)) lql[[i]] <- factor(lql[[i]])
    lql$dichot <- factor(ifelse(lql$fivelev < 4, 0, 1))
    attr(lql,"reshapeLong") <- NULL
    str(lql)
    summary(lql)
    
    ## spaMM expects integer response
    lql$resp <- as.integer(lql$dichot)-1
    # ~9.77s  see tests_other_pack/test-Rasch-lme4.R for other timings
    hl1 <- fitme(cbind(resp,1-resp) ~ 0+itype+(1|subj)+(1|COMPID)+(1|item), 
                 #verbose=c(TRACE=interactive()),
                 data=lql, family=binomial()) 
    testthat::expect_equal(logLik(hl1),c(p_v=-20354.8828))
  } else cat(" requires example_maxtime>7")
} else cat(" requires 'multilevel' package")
cat("\n")
  

