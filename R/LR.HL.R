LR.HL <-
function(HL1,HL2,bootRepl=0) { ## compare two HM objects
   df1 <- ncol(HL1$X)
   df2 <- ncol(HL2$X)
   if (df1==df2){ 
      stop("code needed for test of dispersion params")
   }
   ## else test of fixed effects based on p_v
   if (df1 > df2) {
      fullm <- HL1
      nullm <- HL2
      df <- df1-df2
   } else if (df1 < df2) {
      fullm <- HL2
      nullm <- HL1
      df <- df2-df1
   } 
   ## checking the comparability of REML fits
   if ( ! is.null(fullm$X.Re) ) {
     df.f.Re <-ncol(fullm$X.Re)
   } else df.f.Re <-ncol(fullm$X)
   if ( ! is.null(nullm$X.Re) ) {
     df.n.Re <-ncol(nullm$X.Re)
   } else df.n.Re <-ncol(nullm$X)
   if ( df.f.Re !=  df.n.Re ) {
     warning("LRT comparing REML fits with different designs is highly suspect")
   }
   raw.LR <- 2*(fullm$APHLs$p_v-nullm$APHLs$p_v)
   raw.pvalue <- 1-pchisq(raw.LR,df=df)
   resu <- list(raw.LR=raw.LR,raw.pvalue=raw.pvalue)
   if (bootRepl==0) {
      print("Recommended bootstrap correction:",quote=F)
      print("rerun this analysis with bootRepl=200 (or more)",quote=F)
   } else {
      ## extract model info from both models  -- but we need the calls such as corrHLfit
      computeBootRepl <- function() {
         ## draw sample
         ## analyze under both models
         ## compute and return LR 
      }
      LRtable <- replicate(bootRepl,computeBootRepl)
      LR <- raw.LR * df/mean(LRtable) ## Bartlett-corrected LR
      resu$bootRepl <- bootRepl
      resu$LR <- LR
      resu$pvalue <- 1-pchisq(LR,df=df)
   }
   return(resu)
}
