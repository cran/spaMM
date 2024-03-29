## pseudo family syntax: converts to a list with a $family member
## explicit argument <=> default values visible in doc
## ... <=> no default value eg for responses
multi <- function(binResponse=c("npos","nneg"),binfamily=binomial(),input="types",...) {
  attr(binfamily,"multi") <- TRUE
  return(c(list(family="multi",binResponse=binResponse,binfamily=binfamily,input=input),list(...)))
}

binomialize <- function(data,responses,sortedTypes=NULL,binResponse=c("npos","nneg"),depth=Inf,input="types") {
  if (input=="types") {
    counts <- table(unlist(data[,responses]))
  } else {
    counts <- colSums(data[,responses])
  }
  if(is.null(sortedTypes)) {
    sortedCounts <- sort(counts,decreasing=TRUE) ## most frequent types first  
    sortedTypes <- names(sortedCounts)
  }
  fullsortedTypes <- sortedTypes
  sure_depth <- min(c(length(sortedTypes)-1,depth))
  seq_sure_depth <- seq_len(sure_depth)
  if (input=="types") { ## contents are types
    resu <- lapply(seq_sure_depth, function(v){
      npos <- apply(data[,responses,drop=FALSE]==sortedTypes[1],1,sum) ## 1.2.1 21/07
      nneg <- apply(data[,responses,drop=FALSE],c(1,2), function(v) {v  %in% sortedTypes[-1]})
      nneg <- apply(nneg,1,sum)    
      locdata <- cbind(npos,nneg,data) 
      names(locdata)[1:2] <- binResponse
      sortedTypes <<- sortedTypes[-1]
      locdata[npos+nneg>0,] ## discards "missing data"
    })
  } else { ## colnames are types
    resu <- lapply(seq_len(min(c(length(sortedTypes)-1,depth))), function(v){
      npos <- data[,sortedTypes[1],drop=FALSE]
      nneg <- data[,sortedTypes[-1],drop=FALSE]
      nneg <- apply(nneg,1,sum)    
      locdata <- cbind(npos,nneg,data) 
      names(locdata)[1:2] <- binResponse
      sortedTypes <<- sortedTypes[-1]
      locdata[npos+nneg>0,]
    })    
  }
  ## resu is then a list with unnamed elements. Regular code will not use names but these may be useful for users.
  typeIds <- fullsortedTypes[seq_sure_depth] ## the names of the elements of the list; helpful in case of problem
  names(resu) <- typeIds
  ## then also provide attribute for each element of the list
  for (ll in seq_len(length(resu))) {
    attr(resu[[ll]],"identifier") <- structure(typeIds[ll], names=paste(responses,collapse="&"))
  }
  othername <- setdiff(make.unique(c(names(resu),"others")),names(resu)) ## simply "others" unless it is already in names(resu)
  attr(resu,"sortedTypes") <- c(names(resu),othername) 
  attr(resu,"responses") <- responses
  resu ## a list of data.frames
}

#essai <- binomialize(popBasse,c("Z1_1","Z1_2"))

## function to construct a table of fitted values = multinomial frequencies 
fitted.HLfitlist <- function(object, version=2L, ...) {
  if (version==1L) {
    allrownames <- unique(unlist(lapply(object, function(hl){rownames(hl$data)})))
    fv <- lapply(object, function(hl){
      fv <- fitted(hl,...)
      names(fv) <- rownames(hl$data); fv
    })
    resu <- matrix(0,nrow=length(allrownames),ncol=length(object))
    rownames(resu) <- allrownames
    for (it in seq_len(length(fv))) {
      resu[names(fv[[it]]),it] <- fv[[it]] # Using names(fv[[it]]) leaves 0's where the BinomialDen[[it]] was zero since these cases are excluded from the data
    }
    colnames(resu) <- names(object) ## not convincing... no specific names
  } else {
    resu <- do.call(cbind,lapply(object,predict,newdata=object[[1L]]$data)) 
    colnams <- character(length(object))
    for (it in seq_along(object)) colnams[it] <- attr(object[[it]]$data,"identifier")[[1L]]
    colnames(resu) <- colnams
  }
  for (it in seq_len(ncol(resu)-1L)) {
    for (col in (it+1L):ncol(resu)) resu[,col] <- resu[,col] * (1-resu[,it])
  }
  resu <- cbind(resu,"other(s)"=1-rowSums(resu))
  resu
}

logLik.HLfitlist <- function(object, which = NULL, ...) {
  item <- .getHLfit(object[[1]])
  if (is.null(which)) {
    mess <- .REMLmess(item)
    which <- switch(mess, 
                    "by stochastic EM."= "logLapp",
                    "by Laplace ML approximation (p_v)."= "p_v",
                    "by h-likelihood approximation."= "p_v",
                    "by ML."= "p_v",
                    "by Laplace REML approximation (p_bv)."= "p_bv",
                    "by REML."= "p_bv",
                    "by non-standard REML"= "p_bv",
                    stop(paste0("No default '",which,"' value for '",mess,"' estimation method."))
    ) 
  }
  resu  <- attr(object,"APHLs")[[which]]
  names(resu) <- which
  return(resu)
}  
