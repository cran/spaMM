as.blockDiag.bar <-
function(bar,formula,data,fill=0) { 
  ## FR->FR extension revision expected 
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
      st <- gsub("%in%","/",st) ## FR->FR 06/2014: but try ':' since a %in% b <=> b:a
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
