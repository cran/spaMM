spMMFactorList <-
function (formula, mf, rmInt, drop) {
    ## record dimensions and algorithm settings
    ## create factor list for the random effects
    bars <- spMMexpandSlash(findbarsMM(formula[[3]])) ## lme4::: refs removed
    if (!length(bars)) stop("No random effects terms specified in formula")
    names(bars) <- unlist(lapply(bars, function(x) DEPARSE(x[[3]]))) 
    fl <- lapply(bars, function(x) {
        ## le bidule suivant evalue le bout de formule x[[3]] et en fait un facteur. Useless for spatial effects like longitude + latitude
        ## but fac may be any vector returned by the evaluation of x[[3]] in the envir mf
        txt <- DEPARSE(x[[3]])
        ## evaluating expression(ULI(...)) serves as a way to identify unique combinations of its arguments
        if (length(grep("\\+",txt))>0) { ## coordinates is a vector with a single string; grep is 1 if  + was found in this single string and numeric(0) otherwise
          aslocator <- parse(text=paste("ULI(",gsub("\\+", "\\,", txt),")"))
          ff <- as.factor(eval(expr=aslocator,envir=mf))
        } else if (length(grep("%in%",txt))>0) { ## coordinates is a vector with a single string; grep is 1 if  + was found in this single string and numeric(0) otherwise
          aslocator <- parse(text=paste("ULI(",gsub("%in%", "\\,", txt),")"))
          ff <- as.factor(eval(expr=aslocator,envir=mf))
          # FR->FR  but Design below will provide the Z matrix in ZAL and ideally it would be better to construct a blockDiag object... ?  
        } else if (length(grep("c(\\w*)",txt))>0) { ## c(...,...) was used
          aslocator <-  parse(text=gsub("c\\(", "ULI(", txt)) ## slow pcq ULI() est slow
          ff <- as.factor(eval(expr=aslocator,envir=mf))
        } else { ## standard ( | ) rhs (if a singl variables, it does not matter whether it is spatial or not )
          ff <- eval( expr = substitute(as.factor(fac)[, drop = TRUE], list(fac = x[[3]])), envir = mf)  
          ## or ff <- as.factor(mf[,DEPARSE(x[[3]])])  ## DEPARSE(x[[3]]) should be the rhs of (|) cleanly converted to a string by terms(formula,data) in HLframes
        }
        im <- as(ff, "sparseMatrix") ##FR->FR slow step; but creates S4 objects with slots as assumed in following code
        ##(lme4) Could well be that we should rather check earlier .. :
        if (!isTRUE(validObject(im, test = TRUE))) {
            stop("invalid conditioning factor in random effect: ", format(x[[3]]))
        }
        tempexp <- x[[2]]
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
        ans <- list(f = ff, 
                    A = do.call(rBind, lapply(seq_len(ncol(mm)),function(j) im)),
                    Zt = do.call(rBind, lapply(seq_len(ncol(mm)), 
                                               function(j) {
                                                 im@x <- mm[, j]
                                                 im
                                               })), ## obs <-> levels of ranef
                    ST = matrix(0, ncol(mm), ncol(mm), dimnames = list(colnames(mm), 
                                                                       colnames(mm))),
                    lambda_X=mm ## design matrix for predictor of lambda
                   )
        if (drop) {
            ans$A@x <- rep(0, length(ans$A@x))
            ans$Zt <- drop0(ans$Zt)
        }
        ans
    })
    termsModels <- c()
    Design <- list(0)
    Subject <- list(0)
    namesTerms <- list(0)
    lambda_Xs <- list(0)
    GrpNames <- names(bars)
    for (i in 1:length(fl)) {
      termsModels[i] <- "lamScal" ## FR->FR tempo fix because it's not here that this should be determined
      Subject[[i]] <- as.factor(fl[[i]]$f) # levels of grouping var for all obs
      Design[[i]] <- as.matrix(Matrix::t(fl[[i]]$Zt)) ## nobs *(nr*npar) matrix => used to compute *ZAL* not a model matrix from a formula for lambda 
      nt <- colnames(fl[[i]]$ST) ## length=npar
      namesTerms[[i]] <- nt ## eg intercept of slope... possibly several variables
      names(namesTerms)[i] <- GrpNames[i] ## eg intercept of slope... possibly several variables
    }
    ## to each (.|.) corresponds a Grouping (not necess distinct) and an element of namesTerms each identifying one or more terms
    list(Design = Design, Subject = Subject, Groupings = GrpNames,namesTerms=namesTerms,termsModels=termsModels
         #,lambda_Xs=lambda_Xs
         )
}
