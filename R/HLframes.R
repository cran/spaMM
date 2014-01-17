HLframes <-
function (formula, data, vnms = character(0)) {
    ## m gives either the position of the matched term in the matched call 'mc', or 0
    formula <- asStandardFormula(formula) ## strips out the spatial information
    if (is.character(formula[[2]])) { ## implies that the var designated by a string (phi, for example) should not be in the data frame 
       respname <- formula[[2]]
       if (is.null(data[[respname]])) {
          validname <- respname
       } else {
         datanames <- names(data)
         allphis <- datanames[which(substr(datanames,0,nchar(respname))==respname)] ## "phi..."
         allphis <- substring(allphis,nchar(respname)+1) ## "..."
         allphis <- as.numeric(allphis[which( ! is.na(as.numeric(allphis)))  ]) ## as.numeric("...")
         if (length(allphis) == 0) {
            num <- 0
         } else num <- max(allphis)+1
         validname <-paste ( respname , num , sep="") 
       }
       data[validname] <- 1 ## adds a column $phi of 1 
       formula[[2]] <- as.name(validname) ## now the formula is standard
    }
    mf <- call("model.frame",data=data) ## it adds the formula argument below....
    frame.form <- subbarsMM(formula) ## this comes from lme4 and converts (...|...) terms to some "+" form
    if (length(vnms) > 0) 
        frame.form[[3]] <- substitute(foo + bar, list(foo = parse(text = paste(vnms, 
            collapse = " + "))[[1]], bar = frame.form[[3]]))
    norand.form <- nobarsMM(formula) ## nobars removes the (...|...) terms...
    environment(norand.form) <- environment(frame.form) <- environment(formula)
    mf$formula <- frame.form
    mf$drop.unused.levels <- TRUE
    fe <- mf ## copy before further modif of mf
    mf <- eval(mf)
    Y <- model.response(mf, "any")
    if (! is.null(Y)) {
      ## if binomial Y (may be) a numeric vector and length(dim(Y)) = length(NULL) = 0 
      ## if poisson Y (may be) an integer(!) vector and length(dim(Y)) = length(NULL) = 0
      ## The following looks like an exception
      # if (length(dim(Y)) == 1) { ## 1-D array ? conversion to numeric with problem of keeping the name. 
      ## But as.matrix does this anyway 
      #  nm <- rownames(Y)
      #  dim(Y) <- NULL
      #  if (!is.null(nm)) names(Y) <- nm
      #}
      Y <- as.matrix(Y) ## also useful in binomial case because preprocess tests ncol(Y) later. Revise...
    }
    fixef.form <- nobarsNooffset(formula) ## nobars removes the (...|...) terms...
    if (inherits(fixef.form, "name")) { ## either a var name of a math expression to evaluate, but not cbind...
       X <- matrix(, NROW(Y), 0)
    } else {
       fe$formula <- fixef.form
       fe <- eval(fe)
       mt <- attr(fe, "terms")
       if (mt[[length(mt)]]==0) { 
    	     X <- matrix(nrow=nrow(mf),ncol=0) ## model without fixed effects, not even an Intercept 
       } else if (!is.empty.model(mt)) {
          X <- model.matrix(mt, mf, contrasts)
       } else {
         mess <- pastefrom("no variables identified. Check model formula.",prefix="(!) From ")
         stop(mess)
       }
    }
    storage.mode(X) <- "double" ## not clear what for...
    fixef <- numeric(ncol(X))
    names(fixef) <- colnames(X)
    dimnames(X) <- NULL
    ## creates problem with non-gaussian ranef... because the format of <rand.family>$<member fn>(<matrix>) may be vector or matrix depending on <rand.family>
    ## on the other hand X.pv %*% <vector> is a matrix, so either we try to convert everything to vector (hum) or we just care where it matters...
    list(Y = Y, X = X, wts = NULL, off = NULL, 
        mf = mf, fixef = fixef)
}
