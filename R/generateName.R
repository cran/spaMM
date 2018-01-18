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
    fvname <-paste ( base , num , sep="") 
  }
  fvname
}

#makenewname <- function(...) .makenewname(...) ## removed since new (>08/2017) Infusion version on CRAN

## does not seem to be used
# generateName <- function(base="tmp",
#                             oldnames, ## typically for creating name of new variable
#                             pos=".GlobalEnv" ## typically not used: for creating name of object for debugging
#                             ) {
#   if ( ! missing(oldnames)) { ## because presumably the test should be TRUE if preexisting was provided yet is NULL
#     allmatches <- pmatch(x=base,oldnames)
#   } else {
#     pattern <- paste(base,"*",sep="")
#     allmatches <- ls(pattern=pattern,pos=pos)
#   }
#   allremainders <- substring(allmatches,nchar(base)+1) 
#   allremainders <- as.numeric(allremainders[which( ! is.na(as.numeric(allremainders )))  ]) ## as.numeric("...")
#   if (length(allremainders) == 0L) {
#     validname <- base 
#   } else {
#     num <- max(allremainders)+1
#     validname <-paste ( base , num , sep="") 
#   }
#   return(validname)
# }



## generateFileName() was defined here until  1.10.1 (09/2016) but is moved to blackbox
##base::tempfile generates more random names
#generateFileName <- function(...) stop("generateFileName() has been removed from spaMM:\n update the package calling it (presumably blackbox).")
