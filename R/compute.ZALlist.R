compute.ZALlist <-
function(LMatrix=NULL,ZAlist,Groupings) {## we need to check for user's confusion when we multiply Z by LMatrix
  ## ZAL is nobs * (# levels ranef) and ZA too
  ## LMatrix is (# levels ranef) * (# levels ranef)
  ## the levels of the ranef must match eachother in the two matrices
  ## the only way to check this is to have the levels as rownames and colnames and to check these
  if (is.null(ZAlist)) return(list())
  ## ELSE
  ZAL <- ZAlist
  if ( ! is.null(LMatrix) && length(ZAlist)>0 ) {
    if (inherits(LMatrix,"blockDiag")) {
      stop("compute.ZALlist code should be revised to handle blockDiag objects")
    } ## ELSE
    if ( ! is.list(LMatrix)) LMatrix <- list(LMatrix)
    LMlen <- length(LMatrix)
    for (ii in seq_len(LMlen)) {
      lmatrix <- LMatrix[[ii]]
      ## find ZAlist elements affected by LMatrix element
      affecteds <- which(attr(ZAlist,"ranefs") %in% attr(lmatrix,"ranefs"))
      for (it in affecteds) {
        ZA <- ZAlist[[it]]
        if (attr(ZA,"identityMatrix")) {
          ZAL[[it]] <- lmatrix          
        } else {
          locnc <- ncol(ZA)
          locnr <- nrow(lmatrix)
          if ( locnc %% locnr !=0) {
            mess <- paste("The number of levels of the grouping variable in random term (...|",Groupings[it],")",sep="")
            mess <- paste(mess,"\n  is not the dimension of the correlation matrix.") ## by distMatrix checking in corrHLfit or no.info check somewhere...
            stop(paste(mess," I exit."))
          }         
          nblocks <- locnc %/% locnr 
          if (nblocks>1) {
            locZA <- ZA
            for (bt in 1:nblocks) 
              locZA[,locnr*(bt-1)+(1:locnr)] <- locZA[,locnr*(bt-1)+(1:locnr)] %*% lmatrix
            ZAL[[it]] <- locZA
          } else {
            ## rownames(LMatrix) are the names of first occurrences of unique geographic locations, 
            ## or (if user provided distMatrix) whatever was in this distMatrix. But with a distMatrix, it is likely that ZA was = I and we don't reach this code
            ### it's difficult to make checks on names at this step
            ## LMatrix inherits its names from thos of uniqueGeo. These are the names of first occurrences of unique geographic locations, or (if user provided distMatrix) whatever was in this distMatrix
            ## ZAlist inherits anything from the spMMFactorList call which input does not include info about rownames of data
            #             if ( ! all(attr(ZA,"colnames")==rownames(lmatrix))) {
            #               stop("The colnames of the design matrix Z in eta=...+Zv should be the rownames of the design matrix L  in v=Lu")
            #             }
            ZAL[[it]] <- ZA %*% lmatrix
          }
        }
        attr(ZAL[[it]],"userLfixed") <- attr(lmatrix,"userLfixed") ## TRUE or NULL
      }
    }
  }
  attr(ZAL,"userLfixeds") <- unlist(lapply(ZAL,function(mat) { 
    att <- attr(mat,"userLfixed") ## TRUE or NULL   
    if (is.null(att)) att <- FALSE
    att
  })) ## vector of TRUE or FALSE
  return(ZAL)
}
