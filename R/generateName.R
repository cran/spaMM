# maybe a way to reuse base::tempfile  (-> .Internal) ?

# two functions that diverged from a common idea... probably a merged one should remain.
.makenewname <- function(base,varnames) { ## post CRAN 1.4.1
  varnames <- varnames[which(substring(varnames,1,nchar(base))==base)] 
  allremainders <- substring(varnames,nchar(base)+1) 
  allnumericremainders <- as.numeric(allremainders[which( ! is.na(as.numeric(allremainders )))  ]) ## as.numeric("...")
  ## 2015/03/04
  if (length(allremainders) == 0L && length(allnumericremainders) == 0L) { ## if base = allremainders => length(allnumericremainders) == 0 not sufficient
    fvname <- base
  } else {
    if (length(allnumericremainders)) {
      num <- max(allnumericremainders)+1L
    } else num <- 1L
    fvname <-paste0( base , num) 
  }
  fvname
}

## generateName() was defined here for some time but was not used and was poorly conceived. See Infusion version... 

## generateFileName() was defined here until  1.10.1 (09/2016) but is moved to blackbox
##base::tempfile generates more random names
#generateFileName <- function(...) stop("generateFileName() has been removed from spaMM:\n update the package calling it (presumably blackbox).")
