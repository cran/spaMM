## S3 class blockDiag

##### constructors from al ist of diagonal blocks, from a cumsum partition, from a matrix, from a formula:

as.blockDiag.list <- function(matlist,fill=0L) { ## FR->FR but no check that the matrix is square 
  bd <- matlist
  attr(bd,"rowNames") <- unlist(lapply(seq_len(length(bd)), 
                                       function(v) rownames(bd[[v]], do.NULL = FALSE,prefix=paste(v,".",sep=""))
                                      )
                                ) ## do.NULL=FALSE 10/2014 as non.null names requested by e.g. as.matrix.blockDiag()
  attr(bd,"cumsum") <- cumsum(c(0,unlist(lapply(bd,ncol)))) 
  attr(bd,"fill") <- fill
  class(bd) <- c("blockDiag",class(bd))
  bd  
}

## constructor
as.blockDiag.partition <- function(partition,mat,fill=0L) { ## 
  bd <- lapply(seq_len(length(partition)-1),function(v) {
    sequ <- (partition[v]+1):partition[v+1] 
    mat[sequ,sequ,drop=FALSE]
  })
  attr(bd,"rowNames") <- rownames(mat)  
  attr(bd,"cumsum") <- partition 
  attr(bd,"fill") <- fill
  class(bd) <- c("blockDiag",class(bd))
  bd  
}

## this finds diagonal blocks in a square matrix, progressing from 1,1 to n,n
findblocks <- function(mat,symm=FALSE) {
  if (!symm) {
    mat <- abs(mat)
    mat <- mat+t(mat)
  }
  nc <- ncol(mat)
  if (nc>1) {
    for (it in (2:nc)) {
      chk <- max(mat[1:(it-1),it:nc,drop=FALSE])
      if(chk==0) break
    }  
    if (it==nc) {
      return(nc)
    } else {
      return(c(it-1,findblocks(mat[it:nc,it:nc,drop=FALSE],symm=TRUE)))
    }
  } else {
    return(1)
  }
}

#zut <-diag(1:6);zut[1,3] <-1;zut[5,6] <- 1;  findblocks(zut)
# cf usage in the following def:

## constructor
as.blockDiag.matrix <- function(mat) {
  partition <- findblocks(mat)
  partition <- cumsum(c(0,partition))
  as.blockDiag.partition(partition,mat)
}

#as.blockDiag.matrix(zut)

#### cf remarks on as.blockDiag.bar; as.blockDiag.formula is not even used currently
## constructor
as.blockDiag.formula <-function(formula,data,fill=0) { 
  ############# here I need to get the RHS of the fixed effects expression
  ## frame.form <- subbarsMM(formula) ## this comes from lme4 and converts (...|...) terms to some "+" form
  reExpr <- findbarsMM(formula)
  #  reExpr <- findbarsMM(y~z+AR1(1 | ran/x))
  #  reExpr <- findbarsMM(y~z+AR1(1 | x %in% ran ))
  if (length(reExpr)>1) {
    stop("spaMM does not yet handle multiple random effets. Check the next version.")
  }
  reStruct <- sapply(reExpr,function(v) {v[[3]][[1]]}) 
  ######################################
  for (reidx in seq_len(length(reStruct))) { ## over random effect terms
    relationAschar <- as.character(reStruct[[reidx]])
    if (relationAschar=="/") {
      form <- as.formula(paste("~",paste(reExpr[[reidx]])[[3]])) ## formula of the RHS of the | expr
      mf <- model.frame(formula=formula,data=data) 
      stop("more code needed for nested groups ./. in as.blockDiag") 
      grps <- getGroups(mf, form) ## nlme:::getGroups.data.frame voir le resultat et essayer de l'obtenir sans getGroups (voir %in% ci dessous) FR->FR
      ## et rajouter du code pour obtenir le blockdiagonal object ? inversement remplacer %in% par : et tjrs util getGroups...
    } else if (relationAschar=="%in%") {
      st <- paste(reExpr[[reidx]])[[3]]
      st <- gsub("%in%","/",st)
      fakeform <- as.formula(paste("~",st)) ## so that getGroups can be used
      ## FR->FR mais ce code doit pouvoir etre nettoye pusque getGroups est maintenant evite
      mf <- model.frame(formula=fakeform,data=data) 
      #      grps <- getGroups(mf, fakeform) ## this returns factors and the distance between factors won't be appropriate    
      varlist <-   as.list(unlist(strsplit(paste(reExpr[[reidx]])[[3]],"%in%",fixed=TRUE))) 
      formlist <- lapply(varlist,function(st){as.formula(paste("~",st))})
      names(formlist) <- unlist(lapply(formlist, function(el) DEPARSE(el[[length(el)]]))) ## important !
      grps <- getgroups(formlist,data) ## spaMM:::getgroups ...
      groupingvar <- as.character(fakeform[[2]][[3]])
      groupedvar <- as.character(fakeform[[2]][[2]])
      ## reuse a bit of code for multiple columns = >1 nesting but whole function is not up to this
      for (i in 2:ncol(mf)) {
        grps[, i] <- paste(as.character(grps[, i - 
                                               1]), as.character(grps[, i]), sep = "/")
        NULL ## sic
      } 
      rowNames <-  grps[,groupingvar]
      ##
      rowNames <- rownames(data) ## much more logical 
      rownames(mf) <- rowNames
      ## 
      bd <- lapply(unique(mf[,groupingvar]),function(v) {
        goodrows <- mf[,groupingvar]==v
        xx <- mf[goodrows,groupedvar,drop=FALSE]    
        value <- as.matrix(proxy::dist(xx))
        rownames(value) <- colnames(value) <- rownames(xx)
        value
      })
      attr(bd,"rowNames") <- rowNames  
    }
  }
  attr(bd,"cumsum") <- cumsum(c(0,unlist(lapply(bd,ncol)))) 
  attr(bd,"fill") <- fill
  class(bd) <- c("blockDiag",class(bd))
  bd
}

#bd <- as.blockDiag.formula( formula= yt1~(1|x %in% ran),data=zut )

#### fn def only for handling complex spatial structures eq Matern(1| a %in% b)
## bar is a (.|.) expression
as.blockDiag.bar <-function(bar,formula,data,fill=0) { 
  ## FR->FR extensive revision expected 
  ############# here I need to get the RHS of the fixed effects expression
  #  reExpr <- findbarsMM(y~z+AR1(1 | ran/x))
  #  reExpr <- findbarsMM(y~z+AR1(1 | x %in% ran ))
  reStruct <- bar[[3]][[1]] 
  ######################################
    relationAschar <- as.character(reStruct)
    if (relationAschar=="/") {
      form <- as.formula(paste("~",paste(bar)[[3]])) ## formula of the RHS of the | expr
      mf <- model.frame(formula=formula,data=data) 
      stop("more code needed for nested groups ./. in as.blockDiag") 
      grps <- getGroups(mf, form) ## nlme:::getGroups.data.frame voir le resultat et essayer de l'obtenir sans getGroups (voir %in% ci dessous)
      ## et rajouter du code pour obtenir le blockdiagonal object ?
    } else if (relationAschar=="%in%") {
      st <- paste(bar)[[3]]
      st <- gsub("%in%","/",st) ## FR->FR 06/2014: but try ':' since a %in% b <=> a:b cf code in spMMFactorList
      fakeform <- as.formula(paste("~",st)) ## so that getGroups can be used
      ## FR->FR mais ce code doit pouvoir etre nettoye pusque getGroups est maintenant evite
      mf <- model.frame(formula=fakeform,data=data) 
      #      grps <- getGroups(mf, fakeform) ## this returns factors and the distance between factors won't be appropriate    
      varlist <-   as.list(unlist(strsplit(paste(bar)[[3]],"%in%",fixed=TRUE))) 
      formlist <- lapply(varlist,function(st){as.formula(paste("~",st))})
      names(formlist) <- unlist(lapply(formlist, function(el) DEPARSE(el[[length(el)]]))) ## important !
      grps <- getgroups(formlist,data) ## spaMM:::getgroups ...
      groupingvar <- as.character(fakeform[[2]][[3]])
      groupedvar <- as.character(fakeform[[2]][[2]])
      ## reuse a bit of code for multiple columns = >1 nesting but whole function is not up to this
      for (i in 2:ncol(mf)) {
        grps[, i] <- paste(as.character(grps[, i - 
                                               1]), as.character(grps[, i]), sep = "/")
        NULL ## sic
      } 
      rowNames <-  grps[,groupingvar]
      ##
      rowNames <- rownames(data) ## much more logical 
      rownames(mf) <- rowNames
      ## 
      bd <- lapply(unique(mf[,groupingvar]),function(v) {
        goodrows <- mf[,groupingvar]==v
        xx <- mf[goodrows,groupedvar,drop=FALSE]    
        value <- as.matrix(proxy::dist(xx))
        rownames(value) <- colnames(value) <- rownames(xx)
        value
      })
      attr(bd,"rowNames") <- rowNames          
    }
  attr(bd,"cumsum") <- cumsum(c(0,unlist(lapply(bd,ncol)))) 
  attr(bd,"fill") <- fill
  class(bd) <- c("blockDiag",class(bd))
  bd
}


## as.matrix method
### cf nlme:::pdMatrix.pdBlocked
as.matrix.blockDiag <- function(x,fill=attr(x,"fill"),...) {
  rowNames <- attr(x,"rowNames")
  cumsum <- attr(x,"cumsum")
  Ncol <- length(rowNames)
  value <- array(fill, c(Ncol, Ncol), list(rowNames,rowNames)) ## full matrix of 0's with rows and col names
  for (i in seq_along(x)) {
    aux <- as.matrix(x[[i]],...) ## recursif...
    i.range <- (cumsum[i]+1L):(cumsum[i+1])
    value[i.range,i.range] <- as.vector(aux) ## copie du block dans le block indexé par ses col et row names
  }
  rownames(value) <- colnames(value) <- rowNames
  value
}

#as.matrix(bd)
#as.matrix(bd,NA)

## `^` method (x^blockDiag)
`^.blockDiag` <- function(scal,bdobject) {
  bd <- lapply(1:length(bdobject),function(v){scal^bdobject[[v]]})
  attributes(bd) <- attributes(bdobject)
  bd
}

## this may work but was not tested in a useful way :
## `%*%` method (matrix %*% blockDiag)
`%*%.blockDiag` <- function(mat,bdobject) {
  identifyCols <- c(0,cumsum(unlist(lapply(bdobject,ncol))))
  bd <- lapply(1:length(bdobject),function(v){mat[,seq(identifyCols[[v]]+1,identifyCols[[v+1]])] %*% bdobject[[v]]})
  attributes(bd) <- attributes(bdobject)
  bd
}

#bdpow <- 0.4^bd
#as.matrix(bdpow,0) ## corr matrix; note 0 fill argument in the correlation matrix

## chol method (chol is generic !)
chol.blockDiag <- function(x,...) {
  bd <- lapply(1:length(x),function(v){chol(x[[v]],...)})
  attributes(bd) <- attributes(x)
  bd
}

## RcppChol : uneasy to define as a method for a generic fn 'RcppChol' since RcppChol def is constructed by RcppEigen
`RcppChol.blockDiag` <- function(bdobject) {
  ## chol for each block
  blob <- lapply(1:length(bdobject),function(v){
    glop <- RcppChol(bdobject[[v]])
    if (glop$Status==TRUE) {colnames(glop$L) <- rownames(glop$L) <- rownames(bdobject[[v]])}
    glop
  }) 
  ## collates results
  statuses <- lapply(blob,function(v) {v$Status})
  Status <- all(unlist(statuses)) 
  if (Status==TRUE) {
    bd <- lapply(blob,function(v) {v$L})   ## $L !
    attributes(bd) <- attributes(bdobject)
    ## return(list(L=bd,Status=Status)) ## FR->FR 03/2014 on veut ça et ensuite une ZAL en blocks pour pouvoir travailler sur des blocks de wAugZALI
    ## mais compute.ZAL n'était pas up to this. A revoir
    L <- as.matrix(bd)
    return(list(L=L,Status=Status))
  } else return(list(Status=Status))
}

#chol(bdpow)
#as.matrix(chol(bdpow))-chol(as.matrix(bdpow)) ## 0

## transpose method
t.blockDiag <- function(x) {
  bd <- lapply(1:length(x),function(v){t(x[[v]])})
  attributes(bd) <- attributes(x)
  bd
}
