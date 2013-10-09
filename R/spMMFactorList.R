spMMFactorList <-
function (formula, mf, rmInt, drop) {
    ## record dimensions and algorithm settings
    ## create factor list for the random effects
    bars <- spMMexpandSlash(findbarsMM(formula[[3]])) ## lme4::: refs removed
    if (!length(bars)) stop("No random effects terms specified in formula")
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]]))) 
    fl <- lapply(bars, function(x) {
        ## le bidule suivant evalue le bout de formule x[[3]] et en fait un facteur. Useless for spatial effects like longitude + latitude
        ## but fac may be any vector returned by the evaluation of x[[3]] in the envir mf
        txt <- deparse(x[[3]])
        if (length(grep("\\+",txt))>0) { ## coordinates is a vector with a single string; grep is 1 if  + was found in this single string and numeric(0) otherwise
          aslocator <- parse(text=paste("ULI(",gsub("\\+", "\\,", txt),")"))
          ff <- as.factor(eval(expr=aslocator,envir=mf))
        } else if (length(grep("c(\\w*)",txt))>0) { ## c(...,...) was used
          aslocator <-  parse(text=gsub("c\\(", "ULI(", txt)) 
          ff <- as.factor(eval(expr=aslocator,envir=mf))
        } else { ## standard ( | ) rhs (if a singl variables, it does not matter whether it is spatial or not )
          ff <- eval( expr = substitute(as.factor(fac)[, drop = TRUE], list(fac = x[[3]])), envir = mf)  
          ## or ff <- as.factor(mf[,deparse(x[[3]])])  ## deparse(x[[3]]) should be the rhs of (|) cleanly converted to a string by terms(formula,data) in HLframes
        }
        im <- as(ff, "sparseMatrix")
        ##(lme4) Could well be that we should rather check earlier .. :
        if (!isTRUE(validObject(im, test = TRUE))) {
            stop("invalid conditioning factor in random effect: ", format(x[[3]]))
        }
        if (is.name(x[[2]])) { ## additional bit from HGLMfactorList
            tempexp <- paste("~", as.character(x[[2]]), "-1")
            tempexp <- as.formula(tempexp)[[2]]
        } else tempexp <- x[[2]]
        mm <- model.matrix(eval(substitute(~expr, list(expr = tempexp))), mf)
        if (rmInt) {
            if (is.na(icol <- match("(Intercept)", colnames(mm)))) {
                ## break  ##FR "break may be used in wrong context: no loop is visible"
                ## but the break was in lmer code ! 
                ##FR introduced the 'else'  
            } else {
              if (ncol(mm) < 2) 
                stop("lhs of a random-effects term cannot be an intercept only")
              mm <- mm[, -icol, drop = FALSE]
            }
        }
        ##FR at this point the code diverges from lmerFactorList  
        ans <- list(f = ff, A = do.call(rBind, lapply(seq_len(ncol(mm)), 
            function(j) im)), Zt = do.call(rBind, lapply(seq_len(ncol(mm)), 
            function(j) {
                im@x <- mm[, j]
                im
            })), ST = matrix(0, ncol(mm), ncol(mm), dimnames = list(colnames(mm), 
            colnames(mm))))
        if (drop) {
            ans$A@x <- rep(0, length(ans$A@x))
            ans$Zt <- drop0(ans$Zt)
        }
        ans
    })
    Design <- list(0)
    Subject <- list(0)
    for (i in 1:length(fl)) {
        Subject[[i]] <- as.factor(fl[[i]]$f)
        tempmat <- fl[[i]]$Zt
        tempmat <- as.matrix(Matrix::t(tempmat))
        Design[[i]] <- tempmat
    }
    list(Design = Design, Subject = Subject, namesRE = names(bars))
}
