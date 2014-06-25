spaMMLRT <-
function(null.formula=NULL,formula,
                     null.disp=list(),REMLformula=NULL,
                     method="corrHLfit",boot.repl=0,
                     ## currently trace always false; this is not an argument t be forwarded as is to corrHLfit! 
                     trace=FALSE, ## T means lead to calls of corrHLfit(... trace=list(<file name>,<over/append>))
                     verbose=c(trace=FALSE,warn=NA,summary=FALSE),  
                     ...) {
  if (is.na(verbose["trace"])) verbose["trace"] <- FALSE
  if (is.na(verbose["warn"])) verbose["warn"] <- FALSE ## will be unconditionally ignored by the final fit in corrHLfit  
  if (is.na(verbose["summary"])) verbose["summary"] <- FALSE ## this is for HLCor
  dotlist <-list(...)
  ## birth pangs :
  if ("predictor" %in% names(dotlist)) {
    stop("'spaMMLRT' called with 'predictor' argument which should be 'formula'" )
  }
  if ("null.predictor" %in% names(dotlist)) {
    stop("'spaMMLRT' called with 'null.predictor' argument which should be 'null.formula'" )
  }
  ## here we makes sure that *predictor variables* are available for all data to be used under both models
  data <- dotlist$data
  if ( inherits(data,"list")) {
    data <- lapply(data,function(dt) {
      null.validrows <- validRows(formula=null.formula[-2],resid.formula=dotlist$resid.formula,data=dt) ## will remove rows with NA's in required variables
      full.validrows <- validRows(formula=formula[-2],resid.formula=dotlist$resid.formula,data=dt) ## will remove rows with NA's in required variables
      dt[intersect(null.validrows,full.validrows),,drop=FALSE]     
    })
  } else {
    null.validrows <- validRows(formula=null.formula[-2],resid.formula=dotlist$resid.formula,data=data) ## will remove rows with NA's in required variables
    full.validrows <- validRows(formula=formula[-2],resid.formula=dotlist$resid.formula,data=data) ## will remove rows with NA's in required variables
    data <- data[intersect(null.validrows,full.validrows),,drop=FALSE]     
  }  
  dotlist$data <- data
  predictor <- formula   
  if (! "predictor" %in% class(formula)) predictor <- Predictor(formula)
  null.predictor <- null.formula   
  if (! "predictor" %in% class(null.formula)) null.predictor <- Predictor(null.formula)
  form <- predictor
  if (!is.null(dotlist$LamFix)) {
    dotlist$ranFix$lambda <- dotlist$LamFix
    dotlist$LamFix <- NULL
  }
  if (!is.null(dotlist$PhiFix)) {
    dotlist$ranFix$phi <- dotlist$PhiFix
    dotlist$PhiFix <- NULL
  }  
  dotlist$formula <- predictor 
  dotlist$verbose <- verbose 
  if (method=="corrHLfit") {
    if (is.null(dotlist$objective)) { ## implements default objective function
      if ( ! is.null(null.predictor)) {dotlist$objective <- "p_v"} else {dotlist$objective <- "p_bv"}
    }
  } else dotlist$objective <- NULL
  ## do not change dotlist afterwards !
  fullm.list <- dotlist
  nullm.list <- dotlist
  fullm.list$formula <- predictor
  #### "limited use, for initializing bootstrap replicates:"
  which.iterative.fit <-character(0)
  which.optim.fit <-character(0)
  if ( ! is.null(dotlist$init.corrHLfit$lambda) ) {
    which.optim.fit <- c(which.optim.fit,"lambda")
  } else if ( is.null (fullm.list$ranFix$lambda)) which.iterative.fit <- c(which.iterative.fit,"lambda")
  if (is.null(fullm.list$family) ## (=gaussian)
      || fullm.list$family$family %in% c("gaussian","Gamma")) {
    if ( ! is.null(dotlist$init.corrHLfit$phi) ) { ## if user estimates it by ML...
      which.optim.fit <- c(which.optim.fit,"phi")
    } else which.iterative.fit <- c(which.iterative.fit,"phi")
  }
  if (method=="corrHLfit") {
    if ( ! "rho" %in% names(dotlist$ranFix)) which.optim.fit <- c(which.optim.fit,"rho")
    if ( ! "nu" %in% names(dotlist$ranFix)) which.optim.fit <- c(which.optim.fit,"nu")
    if ( ! "ARphi" %in% names(dotlist$ranFix)) which.optim.fit <- c(which.optim.fit,"ARphi") # mpf
    if ( ! "Nugget" %in% names(dotlist$ranFix)) which.optim.fit <- c(which.optim.fit,"Nugget") ## by default not operative given later code and init.optim$Nugget is NULL
  }
  if ( ! is.null(null.predictor)) { ## ie if test effet fixe
    testFix <- T
    if (dotlist$HLmethod =="SEM") {
      test.obj <- "estlogL"
    } else test.obj <- "p_v"
    ## check fullm.list$REMLformula, which will be copied into nullm in all cases of fixed LRTs
    if (dotlist$HLmethod %in% c("ML","PQL/L","SEM") || substr(dotlist$HLmethod,0,2) == "ML") {
      fullm.list$REMLformula <- NULL
      nullm.list$REMLformula <- NULL
    } else { ## an REML variant
      if (is.null(REMLformula)) { ## default
        fullm.list$REMLformula <- NULL ## HLfit will reconstruct proper REML from this
        nullm.list$REMLformula <- fullm.list$formula ## 
      } else { ## allows alternative choice, but still the same *both* the full and the null fit
        fullm.list$REMLformula <- REMLformula
        nullm.list$REMLformula <- REMLformula
      }
      namesinit <- names(fullm.list$init.corrHLfit)
      namesinit <- setdiff(namesinit,c("rho","nu","Nugget","ARphi"))
      len <- length(namesinit)
      if ( len>0) {
        if (len > 1) namesinit <- paste(c(paste(namesinit[-len],collapse=", "),namesinit[len]),collapse=" and ")
        message("Argument 'init.corrHLfit' is used in such a way that")
        message(paste("  ",namesinit," will be estimated by maximization of p_v.",sep=""))
        message("  'REMLformula' will be inoperative if all dispersion")
        message("  and correlation parameters are estimated in this way.")
      }
    }
    ############## deja converti if (! "predictor" %in% class(null.formula)) null.predictor <- Predictor(null.predictor)
    nullm.list$formula <- null.predictor
  } else if ( length(null.disp)>0 ) { ## test disp/corr param
    testFix <- F
    test.obj <- "p_bv"
    #
    mess <- pastefrom("changes to objective fn required to renullfit or refull fit computation.")
    stop(mess)
    #
    namescheck <- names(fullm.list$lower)
    namescheck <- namescheck[ ! namescheck %in% names(null.disp)]  
    nullm.list$lower <- nullm.list$lower[namescheck]
    nullm.list$upper <- nullm.list$upper[namescheck]
    nullm.list$ranFix <- c(nullm.list$ranFix,null.disp) ## adds fixed values to any preexisting one
  } else testFix <- NA
  
#  trace.info <- NULL
  fullfit <- do.call(method,fullm.list)
  nullfit <- do.call(method,nullm.list)
  if (testFix) {df <- length(fullfit$fixef)-length(nullfit$fixef)} else {df <- length(null.disp)}
  if (df<0) {
    tmp <- fullfit
    fullfit <- nullfit
    nullfit <- tmp
  }
  if (inherits(fullfit,"HLfitlist")) {
    fullL <- attr(fullfit,"APHLs")[[test.obj]]
  } else fullL <- fullfit$APHLs[[test.obj]]
  if (inherits(nullfit,"HLfitlist")) {
    nullL <- attr(nullfit,"APHLs")[[test.obj]]
  } else nullL <- nullfit$APHLs[[test.obj]]
  LRTori <- 2*(fullL-nullL)
  ## BOOTSTRAP
  if ( ! is.na(testFix)) {
    if (boot.repl>0) {
      bootlist <- dotlist ## copies ranFix
      bootlist <- c(bootlist,list(null.formula=null.predictor,null.disp=null.disp,REMLformula=REMLformula,method=method)) ## unchanged user REMLformula forwarded
      bootlist$verbose <- c(trace=FALSE,summary=FALSE)
      bootlist$trace <- FALSE 
      bootlist$boot.repl <- 0 ## avoids recursive call of bootstrap
      all.estim.ranvars <- c(which.optim.fit) ## FR->FR to be modified if optimFits are reintroduced, see corresponding code in corrMM.LRT 
      for (st in all.estim.ranvars) {
        if (dotlist$HLmethod=="corrHLfit") {
          if ( st %in% names(nullfit$corrPars)) {
            bootlist$init.corrHLfit[st] <- nullfit$corrPars[st]
          } else if ( st %in% names(nullfit)) {
            bootlist$init.corrHLfit[st] <- nullfit[st] ## handled in dotlist by corrHLfit, cf notes 090113
          } ## it's also possible that st is nowhere (10/2013: currently for Nugget)
        } else {
          if ( st %in% names(nullfit)) {
            bootlist$init.HLfit[st] <- nullfit[st] ## 
          } 
        }
      }        
      bootreps<-matrix(,nrow=boot.repl,ncol=2) 
      colnames(bootreps) <- paste(c("full.","null."),test.obj,sep="")
      cat("bootstrap replicates: ")
      simbData <- nullfit$data
      if (tolower(nullfit$family$family)=="binomial") {
        form <- attr(nullfit$predictor,"oriFormula") ## this must exists...  
        if (is.null(form)) {
          mess <- pastefrom("a 'predictor' object must have an 'oriFormula' member.",prefix="(!) From ")
          stop(mess)
        }
      }
      ## the data contain any original variable not further used; e.g original random effect values in the simulation tests  
      thisFnName <- as.character(sys.call()[[1]]) ## prevents a bug when we change "this" function name
      for (ii in 1:boot.repl) {
        locitError <- 0
        repeat { ## for each ii!
          newy <- simulate(nullfit) ## cannot simulate all samples in one block since some may not be analyzable  
          if (tolower(nullfit$family$family)=="binomial") {
            ## c'est bouseux: soit j'ai (pos, neg) et le remplacement est possible
            ##    soit j'ai (pos,ntot -pos) et le 2e remplacment n'est pas poss (et pas necess)
            ##    aussi (ntot - pos, pos) ...
            ## would be simple if always ntot-pos, but how to control this ? 
            ## simbData[[as.character(form[[2]][[2]])]] <- newy
            ## simbData[[as.character(form[[2]][[3]])]] <- nullfit$weights - newy    
            exprL <- as.character(form[[2]][[2]]) 
            exprR <- as.character(form[[2]][[3]]) 
            if (length(exprL)==1) simbData[[exprL]] <- newy 
            if (length(exprR)==1) simbData[[exprR]] <- nullfit$weights - newy                    
            ## if (length(exprR)! =1) exprRdoes not correspond to a column in the data;frmae so there is no column to replace                     
          } else {simbData[[as.character(nullfit$predictor[[2]])]] <- newy}
          bootlist$data <- simbData
          bootrepl <- try(do.call(thisFnName,bootlist)) ###################### CALL ##################
          if (class(bootrepl)[1] != "try-error") { ## eg separation in binomial models... alternatively, test it here (require full and null X.pv... )
            if (inherits(bootrepl$fullfit,"HLfitlist")) {
              fullL <- attr(bootrepl$fullfit,"APHLs")[[test.obj]]
            } else fullL <- bootrepl$fullfit$APHLs[[test.obj]]
            if (inherits(bootrepl$nullfit,"HLfitlist")) {
              fullL <- attr(bootrepl$nullfit,"APHLs")[[test.obj]]
            } else fullL <- bootrepl$nullfit$APHLs[[test.obj]]
            bootreps[ii,] <- c(bootrepl$fullfit$APHLs[[test.obj]],bootrepl$nullfit$APHLs[[test.obj]])
            break ## replicate performed, breaks the repeat
          } else { ## there was one error
            locitError <- locitError + 1
            if (locitError>10) { ## to avoid an infinite loop
              stop("Analysis of bootstrap samples fails repeatedly. Maybe no statistical information in them ?")
            } ## otherwise repeat!
          }
        } 
        cat(ii);cat(" ")
        if ((ii %% 50)==0) cat("\n")
      } ## end main bootstrap loop
      cat("\n") ##  
    } ## end bootstrap 
  } else { ## nothing operativ yet
    bootreps<-matrix(,nrow=boot.repl,ncol=length(unlist(fullfit$APHLs))) 
    colnames(bootreps) <- names(unlist(fullfit$APHLs))
    ## more needed here ?
  }
  ## prepare output
  if ( ! is.na(testFix)) {
    if (testFix) {df <- length(fullfit$fixef)-length(nullfit$fixef)} else {df <- length(null.disp)}
    ## pvalue <- 1-pchisq(LRT,df=df)
    resu <- list(fullfit=fullfit,nullfit=nullfit)
    LRTinfo <- list(df=df,LRTori = LRTori)
    if (boot.repl>0) {
      meanbootLRT <- 2*mean(bootreps[,1]-bootreps[,2])  
      LRTinfo$meanbootLRT <- meanbootLRT
      LRTinfo$bootreps <- bootreps
    }
  } else {
    resu <- list(fullfit=fullfit)
    LRTinfo <- list()
  }
#  LRTinfo$trace.info <- trace.info 
  ##  resu$LRTinfo <- LRTinfo ## pas compatible avec hglmjob.R...
  resu <- c(resu,LRTinfo) ## loses the sublist structure, which wouldnot be compatible with hglmjob.R...  
  class(resu) <- c("fixedLRT",class(resu)) 
  return(resu)
}
