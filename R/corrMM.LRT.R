corrMM.LRT <-
function(null.predictor=NULL,predictor,
                       null.disp=list(),REMLformula=NULL,
                       method=corrHLfit,boot.repl=0,
                       which.iterative=c(),
                       trace=F, ## T means lead to calls of corrHLfit(... trace=list(<file name>,<over/append>))
                       control=list(), ## profiles=Inf,REMLfits=T,optimFits,restarts,maxit...
                       control.boot=list(), ## REMLfits=F,optimFits=T,profiles=0
                       verbose=c(F,F),  
                       ...) {
  ##################################################                     
  ## methods handled:
  ## RE<...> (ie. 'REML', 'RE(0,1)'...) => REML fits using the iterative algorithm, both with design matrix from full models; a bit slow...
  ## Other methods: no REML correction but the same REML fits are still possible to provide starting values for other fits
  ## if an init.corrHLfit argument is provided, then there is no REML estimation of the given parameter, even if the given method is RE<...>
  ## Original PQL would be RE(0,0)
  ## monPQL was ML(0,1)
  ## my HL(...) were ML(...) because of init.corrHLfit, but would otherwise be RE(...)
  ## FR->FR still problems ML(...) + REML formula
  ##################################################                     
  profiles <- control$profiles
  if (is.null(profiles)) profiles <- Inf
  restarts <- control$restarts
  if (is.null(restarts)) restarts <- TRUE
  maxit <- control$maxit
  if (is.null(maxit)) maxit <- 1
  # added 10/2013
  bootcontrol <- list()
  bootcontrol$REMLfits <- control.boot$REMLfits
  if (is.null(bootcontrol$REMLfits)) bootcontrol$REMLfits <- FALSE
  bootcontrol$optimFits <- control.boot$optimFits 
  if (is.null(bootcontrol$optimFits)) bootcontrol$optimFits <- TRUE
  bootcontrol$profiles <- control.boot$profiles
  if (is.null(bootcontrol$profiles)) bootcontrol$profiles <- 0
  #
  bootFix <- control$bootFix
  if (is.null(bootFix)) bootFix <- c()
  if (class(predictor)=="formula") predictor <- Predictor(predictor)
  form <- predictor$formula
  dotlist <-list(...)
  dotlist$minimalOutput <- NULL ## FR->FR no longer used 
  if (!is.null(dotlist$LamFix)) {
    dotlist$ranFix$lambda <- dotlist$LamFix
    dotlist$LamFix <- NULL
  }
  if (!is.null(dotlist$PhiFix)) {
    dotlist$ranFix$phi <- dotlist$PhiFix
    dotlist$PhiFix <- NULL
  }  
  dotlist$predictor <- predictor ## this is here so that the user sees 'predictor' as an argument
  dotlist$verbose <- verbose 
  if (is.null(dotlist$objective)) { ## implements default objective function
    if ( ! is.null(null.predictor)) {dotlist$objective <- "p_v"} else {dotlist$objective <- "p_bv"}
  }
  ## do not change dotlist afterwards !
  fullm.list <- dotlist
  ####
  which.iterative.fit <-character(0)
  which.optim.fit <-character(0)
  if ( ! is.null(dotlist$init.corrHLfit$lambda) ) { ## if user estimates it by ML...
    which.optim.fit <- c(which.optim.fit,"lambda")
  } else if ( is.null (fullm.list$ranFix$lambda)) which.iterative.fit <- c(which.iterative.fit,"lambda")
  if (is.null(fullm.list$family) ## (=gaussian)
      || fullm.list$family$family %in% c("gaussian","Gamma")) {
    if ( ! is.null(dotlist$init.corrHLfit$phi) ) { ## if user estimates it by ML...
      which.optim.fit <- c(which.optim.fit,"phi")
    } else which.iterative.fit <- c(which.iterative.fit,"phi")
  }
  if ( ! "rho" %in% names(dotlist$ranFix)) which.optim.fit <- c(which.optim.fit,"rho")
  if ( ! "nu" %in% names(dotlist$ranFix)) which.optim.fit <- c(which.optim.fit,"nu")
  if ( ! "Nugget" %in% names(dotlist$ranFix)) which.optim.fit <- c(which.optim.fit,"Nugget") ## currently not operative
  ######## Determines iterativeFits, optimFits and the corresponding REML formulas
  ## *by default* the iterative fits are more REML than the optim fits; true REML is possible only through the iterative fits.
  #### optimREMLformula is a (false, typically ML) REML formula for nullfit, fullfit, not for REML fits
  if (dotlist$HLmethod %in% c("ML","PQL","PQL/L")) { ## corrMM.LRT assumes a non-standard form of PQL
    iterativeFits <- control$REMLfits ## FR->FR 'REMLfits' should become obsolete syntax
    optimFits <- control$optimFits ## 
    if (is.null(optimFits)) optimFits <- TRUE 
    if (optimFits) optimREMLformula <- NULL  ## HLfit (translate...) will reject any non null input REMLformula and reconstruct it internally  
  } else if (substr(dotlist$HLmethod,0,2) == "RE") { ## explicit request for REML; only iterativeFits should be run and optimREMLformula shoul not be used 
    iterativeFits <- TRUE
    optimFits <- FALSE
  } else {  ## HL<...>: in corrMM.LRT, this is interpreted as ML fits; although this is also non-standard   
    iterativeFits <- control$REMLfits ## FR->FR 'REMLfits' should become obsolete syntax
    bars <- findbarsMM(form)  ## extract random effects
    lhs <- paste(form[[2]]) ## extract response
    ## build formula with only random effects
    optimFits <- control$optimFits ## 
    if (is.null(optimFits)) optimFits <- TRUE 
    if (optimFits) optimREMLformula <- as.formula(paste(lhs,"~",paste(paste("(",as.character(bars),")"),collapse="+"))) ## ML in iterative fit
  }
  ####
  if (is.null(iterativeFits)) { ## ie default for non RE<...>
    if (length(which.iterative)==0) {
      iterativeFits <- TRUE ## we can construct a non trivial iterative fit that usefully can be used to initiate the optim fits
    } else {
      iterativeFits <- FALSE ## then the optim fits are iterative wrt to lambda in *G*LMMs; then we usually don't need the REML fits to find a good starting lambda, 
                             ## although we can still explicitely request them 
    }
  } 
  ####
  #### REMLformula for iterativeFits is a false (by default full-model formula) REML formula for both iterative fits; that makes them *somewhat* more suitable for LR tests
  if (iterativeFits) {
    if (is.null(REMLformula)) { ## default
      ## this one is the formula for the true REML full fit and the hacked (nullm.list$REMLformula=full...) REML null fit that serve to initiate the later ML fits
      fullm.list$REMLformula <- fullm.list$predictor$formula ## REML correction by default under full model (this is also the HLfit default)
    } else { ## allows alternative choice, but still the same *both* the full and the null fit
      fullm.list$REMLformula <- REMLformula
      namesinit <- names(fullm.list$init.corrHLfit)
      if ( ! is.null(namesinit)) {
        message("Argument 'init.corrHLfit' is used in such a way that")
        message(paste("  ",namesinit," will be estimated by maximization of p_v.",sep=""))
        message("  'REMLformula' will be inoperative if all dispersion")
        message("  and correlation parameters are estimated in this way.")
      }
    }  
  }
  if (verbose[1]) {
    if (iterativeFits) {
      cat("corrMM.LRT will perform iterative fits of dispersion parameters with REML formula ")
      print(fullm.list$REMLformula,showEnv=F)
    } else print("corrMM.LRT will not perform iterative fits of dispersion parameters.")  
    if (optimFits) {
      cat("corrMM.LRT will perform generic optimization fits of dispersion parameters with REML formula ")
      print(optimREMLformula,showEnv=F)
    } else print("corrMM.LRT will not perform generic optimization fits of dispersion parameters.")
  }
  #########
  ## definitions for updating parameters
  ## local fn:
  `update.ranef.pars` <- function(from.fit,to.arglist,which.pars=c("rho","nu")) {
     ## default 'which.pars' prevents the dispersion params from being initialized by this fn
     ## using lambda leads to some poor results
     if ("lambda" %in% which.pars) {
        if ("lambda" %in% which.optim.fit ) {
          to.arglist$init.corrHLfit$lambda <- from.fit$lambda
        } else if ("lambda" %in% which.iterative.fit ) {
          to.arglist$init.HLfit$lambda <- from.fit$lambda
        } 
     } ## ELSE keeps any preexisting to.arglist$init.corrHLfit$lambda
     if ("phi" %in% which.pars) {
        if ("phi" %in% which.optim.fit ) {
          to.arglist$init.corrHLfit$phi <- from.fit$phi
        } else if ("phi" %in% which.iterative.fit ) {
          to.arglist$init.HLfit$phi <- from.fit$phi
        }
     }
     if ("rho" %in% which.pars) {
        if ("rho" %in% which.optim.fit )  { ## always by ML whether user provided initial values or not
          to.arglist$init.corrHLfit$rho <- from.fit$corrPars$rho   
        }    
     }  
     if ("nu" %in% which.pars) {
        if ("nu" %in% which.optim.fit) { ## always by ML whether user provided initial values or not
          to.arglist$init.corrHLfit$nu <- from.fit$corrPars$nu   
        }
     }      
     if ("v_h" %in% which.pars) { ## it is a very bad idea to update with v_h by default. Cf Ln... (binary)
        to.arglist$init.HLfit$v_h <- from.fit$v_h ## added 04/12/12
     }
    return(to.arglist)
  }
  ## another lcoal fn. ! different default which.pars
  init.fixed.ranefpars <- function(from.fit,to.arglist,which.pars=c("rho","nu","lambda","phi")) {
     ## default 'which.pars' prevents the dispersion params from being initialized by this fn
     ## using lambda leads to some poor results
     if ("lambda" %in% which.pars) {
        if ("lambda" %in% which.optim.fit ) {
          to.arglist$ranFix$lambda <- from.fit$lambda
        } else if ("lambda" %in% which.iterative.fit ) {
          to.arglist$ranFix$lambda <- from.fit$lambda
        } 
        to.arglist$init.corrHLfit$lambda <- NULL
     } ## ELSE keeps any preexisting to.arglist$init.corrHLfit$lambda
     if ("phi" %in% which.pars) {
        if ("phi" %in% which.optim.fit ) {
          to.arglist$ranFix$phi <- from.fit$phi
        } else if ("phi" %in% which.iterative.fit ) {
          to.arglist$ranFix$phi <- from.fit$phi
        }
        to.arglist$init.corrHLfit$phi <- NULL
     }
     if ("rho" %in% which.pars) {
        if ("rho" %in% which.optim.fit )  { ## always by ML whether user provided initial values or not
          to.arglist$ranFix$rho <- from.fit$corrPars$rho   
          to.arglist$init.corrHLfit$rho <- NULL
        }    
     }  
     if ("nu" %in% which.pars) {
        if ("nu" %in% which.optim.fit) { ## always by ML whether user provided initial values or not
          to.arglist$ranFix$nu <- from.fit$corrPars$nu   
          to.arglist$init.corrHLfit$nu <- NULL
        }
     }      
    return(to.arglist)
  }  
  ##
  trace.info <- NULL
  nullm.list <- dotlist
  ####
  #### ALWAYS p_v for testing fixed effects, p_bv for random effects 
  if ( ! is.null(null.predictor)) { ## ie if test effet fixe
    testFix <- T
    test.obj <- "p_v"
    if (class(null.predictor) == "formula") null.predictor <- Predictor(null.predictor)
    nullm.list$predictor <- null.predictor
    nullm.list$REMLformula <- fullm.list$REMLformula ## REML correction still under full model
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
  ####
  #### run iterative fit
  ## FIRST NULL FIT; iterative then optim
  if (iterativeFits) { ## user requested not iterative (default), but in this case the algo still tries something for iterative lambda 
    nullREML.list <- nullm.list
    if (nullREML.list$HLmethod=="ML") nullREML.list$HLmethod <- "REML" ## but do not over write more explicit HL(.,.) statements 
    if (nullREML.list$HLmethod %in% c("PQL","PQL/L")) nullREML.list$HLmethod <- "REPQL" ## but do not over write more explicit HL(.,.) statements 
    nullREML.list$init.corrHLfit$lambda <- NULL 
    nullREML.list$init.corrHLfit$phi <- NULL ## 21/02/2013... changes something only for LMMs, which are not an issue
    nullREML.list$control.HLfit$conv.threshold <- 1e-04 
    if (trace) nullREML.list$trace <- list(file="trace.lamREMLnullfit.txt",append=F)
    lamREMLnullfit <- do.call(method,nullREML.list)$hlfit ## must return something that simulate() can manipulate
    trace.info <- data.frame(iter=0,step="lamREMLnullfit",obj=lamREMLnullfit$APHLs[[test.obj]])
    ## for ML 
    nullm.list <- update.ranef.pars(from.fit=lamREMLnullfit,to.arglist=nullm.list) 
    nullm.list <- update.ranef.pars(from.fit=lamREMLnullfit,to.arglist=nullm.list,which.pars=c("lambda")) 
  }  else lamREMLnullfit <- NULL 
  ####
  #### if not RE<...>, run optim fit
  if (optimFits) {
    if ("phi" %in% which.iterative) nullm.list$init.corrHLfit$phi <- NULL
    if ("lambda" %in% which.iterative) nullm.list$init.corrHLfit$lambda <- NULL
    nullm.list$REMLformula <- optimREMLformula ## ie ML fit... always the case anyway but needed if iterative...
    if (trace) nullm.list$trace <- list(file="trace.nullfit.txt",append=F)
    nullfit <- do.call(method,nullm.list)$hlfit ## must return something that simulate() can manipulate
    if ( ! is.null(lamREMLnullfit)) { ## if an iterative fit was performed, we check whether it provided a better fit
        ## since lamREMLnullfit use full model formula and nullfit uses ML, one should not compare their p_v
        ## instead we fix the ranPars [since HLCor, not corrHLfit, call] and refit by ML
        MLforREMLranPars <- nullREML.list
        MLforREMLranPars$ranFix <- NULL
        MLforREMLranPars$REMLformula <- optimREMLformula
        MLforREMLranPars$ranPars<- with(lamREMLnullfit,c(corrPars,list(phi=phi,lambda=lambda)))
        MLforREMLranPars$coordinates <- dotlist$coordinates
        gnrf <- do.call(HLCor,MLforREMLranPars)$hlfit$APHLs[[test.obj]]
        if (gnrf > nullfit$APHLs[[test.obj]]) { ## annoying as it means ML has failed to maximize p_v, but can occur
          trace.info <- rbind(trace.info,data.frame(iter=0,step="(!) nullfit (-)",obj=nullfit$APHLs[[test.obj]]))
          nullfit <-lamREMLnullfit ## step back
          ## this suggests a divergence of the inner loop; cf notes for 29/12/12 for case Hn
          ## Using which.iterative may be a fix
        } else trace.info <- rbind(trace.info,data.frame(iter=0,step="nullfit (+)",obj=nullfit$APHLs[[test.obj]]))
    } else trace.info <- rbind(trace.info,data.frame(iter=0,step="nullfit",obj=nullfit$APHLs[[test.obj]]))
  } else {
    nullfit <- lamREMLnullfit ### 
  }
## THEN FIRST FULL FIT, iterative then optim
  ## iterative
  fullm.list <- update.ranef.pars(from.fit=nullfit,to.arglist=fullm.list) 
  if (iterativeFits ) { ## in which case still tries somethng for iterative lambda 
    REML.list <- fullm.list
    if (REML.list$HLmethod=="ML") REML.list$HLmethod <- "REML" ## but do not over write more explicit HL(.,.) statements 
    if (REML.list$HLmethod %in% c("PQL","PQL/L")) REML.list$HLmethod <- "REPQL" ## but do not over write more explicit HL(.,.) statements 
    REML.list$init.corrHLfit$lambda <- NULL 
    REML.list$init.corrHLfit$phi <- NULL ## 21/02/2013
    REML.list$control.HLfit$conv.threshold <- 1e-04 
    if (trace) REML.list$trace <- list(file="trace.lamREMLfullfit.txt",append=F)
    lamREMLfullfit <- do.call(method,REML.list)$hlfit ## this performs an REML fit for lambda, provide the initial value for the ML fit
    trace.info <- rbind(trace.info,data.frame(iter=0,step="lamREMLfullfit",obj=lamREMLfullfit$APHLs[[test.obj]]))
    ## for ML
    fullm.list <- update.ranef.pars(from.fit=lamREMLfullfit,to.arglist=fullm.list)  
    fullm.list <- update.ranef.pars(from.fit=lamREMLfullfit,to.arglist=fullm.list,which.pars=c("lambda")) 
  } else lamREMLfullfit <- NULL 
  if (optimFits) {
    if ("phi" %in% which.iterative) fullm.list$init.corrHLfit$phi <- NULL
    if ("lambda" %in% which.iterative) fullm.list$init.corrHLfit$lambda <- NULL
    fullm.list$REMLformula <- optimREMLformula ## ie ML fit... always the case anyway but needed if iterative...
    if (trace) fullm.list$trace <- list(file="trace.fullfit.txt",append=F)
    fullfit <- do.call(method,fullm.list)$hlfit ## must return something that simulate() can manipulate
    if ( ! is.null(lamREMLfullfit)) { ## if an iterative fit was performed, we check whether it provided a better fit
        MLforREMLranPars <- REML.list
        MLforREMLranPars$ranFix <- NULL
        MLforREMLranPars$REMLformula <- optimREMLformula
        MLforREMLranPars$ranPars<- with(lamREMLfullfit,c(corrPars,list(phi=phi,lambda=lambda)))
        MLforREMLranPars$coordinates <- dotlist$coordinates
        gnrf <- do.call(HLCor,MLforREMLranPars)$hlfit$APHLs[[test.obj]]
        if (gnrf > fullfit$APHLs[[test.obj]]) { ## annoying as it means ML has failed to maximize p_v, but can occur
          trace.info <- rbind(trace.info,data.frame(iter=0,step="(!) fullfit (-)",obj=fullfit$APHLs[[test.obj]]))
          fullfit <-lamREMLfullfit
        } else trace.info <- rbind(trace.info,data.frame(iter=0,step="fullfit (+)",obj=fullfit$APHLs[[test.obj]]))
    } else trace.info <- rbind(trace.info,data.frame(iter=0,step="fullfit",obj=fullfit$APHLs[[test.obj]]))
  } else {
    fullfit <- lamREMLfullfit
  }
  ## ITERATIONS
  nullnamesX <- colnames(nullfit$X)
  fullnamesX <- colnames(fullfit$X)
  nullVarsPos <- which( fullnamesX %in% nullnamesX) 
  profileVars <- which( ! fullnamesX %in% nullnamesX) ## les parametres testes
  profileX <- fullfit$X[,profileVars,drop=F]
  locit <- 1
  if (dotlist$HLmethod %in% c("REML","REPQL")) {
    conv <- 0 ## will skip the loop
    LRTori <- 2*(fullfit$APHLs[[test.obj]]-nullfit$APHLs[[test.obj]])
  } else  conv <- 1e08
  newnull <- T; newfull <- T; 
  newforprof <- T
  LRTprof <- NA
  while ( optimFits && conv>0.001 && (newnull || newfull) && locit <= maxit ) { ## optimFits required since the later fits use optimREMLformula; but the premise could be reconsidered
    ####### this loop aims to ensure that maximal nullfit and fullfit are found
    ## currently implemented only in cases where some random effect parameters are fit
    ## (test below; no alternative yet)
    #######
    conv <-0 
    prof.needed <- F ## must be reset at the beginning of each loop
    newnull <- F
    newfull <- F
    if (length(c(which.iterative.fit,which.optim.fit))>0) { ## some random effect parameters are fit 
      ####### while fullfit or nullfit has been changed in the last iter
      # if (fullfit seems OK) then use it to refine nullfit: {
      #   renullfit with update.ranef.pars
      #   renullfit with further update lambda, and iter.mean.dispFix <- 10
      #   retain the best (so that update lambda, and iter.mean.dispFix <- 10 can contaminate later renullfit and the profile) 
      # }
      # if (fullfit lower than refined nullfit) then use nullfit to refine fullfit: {
      #   refullfit with update.ranef.pars and init.HLfit$fixef[nullVarsPos] <- as.numeric(nullfit$fixef)
      #   renullfit with further update lambda, and iter.mean.dispFix <- 10
      #   retain the best (so that update lambda, and iter.mean.dispFix <- 10 can contaminate later renullfit) 
      # }
      #######      
if (restarts) { 
      if ( fullfit$APHLs[[test.obj]] > nullfit$APHLs[[test.obj]]) { 
        ## ideally always the case. Then we look whether we can improve nullfit using fullfit as starting point
        ## (else, we will improve fullfit first)
        ## REFIT NULL FIT
        renullm.list <- nullm.list
        if (trace) renullm.list$trace <- list(file="trace.renullfit.txt",append=F)
        renullm.list <- update.ranef.pars(from.fit=fullfit,to.arglist=renullm.list)
        if ("phi" %in% which.iterative) renullm.list$init.corrHLfit$phi <- NULL
        if ("lambda" %in% which.iterative) renullm.list$init.corrHLfit$lambda <- NULL
        renullfit <- do.call(method,renullm.list)$hlfit 
        if (renullfit$APHLs[[test.obj]] > nullfit$APHLs[[test.obj]]) {
          conv <- conv + renullfit$APHLs[[test.obj]] - nullfit$APHLs[[test.obj]]
          nullfit <-renullfit
          nullm.list <- renullm.list ## updated for immediate use and use by profilize
          newnull <- T
          trace.info <- rbind(trace.info,data.frame(iter=locit,step="renullfit.1 (+)",obj=nullfit$APHLs[[test.obj]]))
        } 
        if (trace) renullm.list$trace <- list(file="trace.renullfit.txt",append=T)
        renullm.list <- update.ranef.pars(from.fit=fullfit,to.arglist=renullm.list)
        renullm.list <- update.ranef.pars(from.fit=fullfit,to.arglist=renullm.list,which.pars=c("lambda"))
        ############################ renullm.list$control.HLfit$iter.mean.dispFix <- 10
        if ("phi" %in% which.iterative) renullm.list$init.corrHLfit$phi <- NULL
        if ("lambda" %in% which.iterative) renullm.list$init.corrHLfit$lambda <- NULL
        renullfit <- do.call(method,renullm.list)$hlfit 
        if (renullfit$APHLs[[test.obj]] > nullfit$APHLs[[test.obj]]) {
          conv <- conv + renullfit$APHLs[[test.obj]] - nullfit$APHLs[[test.obj]]
          nullfit <-renullfit
          nullm.list <- renullm.list ## updated for use by profilize
          newnull <- T
          trace.info <- rbind(trace.info,data.frame(iter=locit,step="renullfit.2 (+)",obj=nullfit$APHLs[[test.obj]]))
        } 
      } 
      ## two heuristic tries for large improvements, else then a sure small improvement
      if ( fullfit$APHLs[[test.obj]] < nullfit$APHLs[[test.obj]]) { ## ideally never the case
        ## REFIT FULL FIT
        ## then we try other starting values
        ## the starting lambda value seems to be quite important
        refullm.list <- fullm.list
        if (trace) refullm.list$trace <- list(file="trace.refullfit.txt",append=F)
        refullm.list <- update.ranef.pars(from.fit=nullfit,to.arglist=refullm.list)
        refullm.list$init.HLfit$fixef <- rep(0,length(fullnamesX)) ## fixef must be modified here
#        refullm.list$init.HLfit$fixef[nullVarsPos] <- as.numeric(nullfit$fixef) 
        refullm.list$init.HLfit$fixef[nullVarsPos] <- as.numeric(mean(nullfit$eta)) 
        if ("phi" %in% which.iterative) refullm.list$init.corrHLfit$phi <- NULL
        if ("lambda" %in% which.iterative) refullm.list$init.corrHLfit$lambda <- NULL
        refullfit <- do.call(method,refullm.list)$hlfit
        if ( refullfit$APHLs[[test.obj]] > fullfit$APHLs[[test.obj]] ) { ## 
          conv <- conv + refullfit$APHLs[[test.obj]] - fullfit$APHLs[[test.obj]] ## only D(FULL)
          fullfit <- refullfit
          fullm.list <- refullm.list ## updated for immediate use 
          newfull <- T
          trace.info <- rbind(trace.info,data.frame(iter=locit,step="refullfit.1 (+)",obj=fullfit$APHLs[[test.obj]]))
        }
      }
      if ( fullfit$APHLs[[test.obj]] < nullfit$APHLs[[test.obj]]) { ## if previous heuristic attempt not good enough
        ####### try something else
        if (trace) refullm.list$trace <- list(file="trace.refullfit.txt",append=T)
        refullm.list <- update.ranef.pars(from.fit=nullfit,to.arglist=refullm.list,which.pars=c("lambda"))
        ## ie, essentially, refullm.list$init.corrHLfit$lambda <- nullfit$lambda, which works *sometimes*          
        ################## refullm.list$control.HLfit$iter.mean.dispFix <- 10
        if ("phi" %in% which.iterative) refullm.list$init.corrHLfit$phi <- NULL
        if ("lambda" %in% which.iterative) refullm.list$init.corrHLfit$lambda <- NULL
        refullfit <- do.call(method,refullm.list)$hlfit
        if ( refullfit$APHLs[[test.obj]] > fullfit$APHLs[[test.obj]] ) { ## 
          conv <- conv + refullfit$APHLs[[test.obj]] - fullfit$APHLs[[test.obj]] ## only D(FULL)
          fullfit <- refullfit
          fullm.list <- refullm.list ## maybe not useful
          newfull <- T
          trace.info <- rbind(trace.info,data.frame(iter=locit,step="refullfit.2 (+)",obj=fullfit$APHLs[[test.obj]]))
        } 
      }
      if ( fullfit$APHLs[[test.obj]] < nullfit$APHLs[[test.obj]] ) { ## if previous heuristic attemptS not good enough
        ##FR->FR formally a call to corrHLfit, but all optim parameters are fixed. This is is only an HLCor call, and it generates no valid output in the trace file
        fullm.con.list <- fullm.list
        fullm.con.list <- init.fixed.ranefpars(from.fit=nullfit,to.arglist=fullm.con.list)
        fullm.con.list <- update.ranef.pars(from.fit=nullfit,to.arglist=fullm.con.list,which.pars=c("v_h"))
        fullm.con.list$init.HLfit$fixef <- rep(0,length(fullnamesX))
#        fullm.con.list$init.HLfit$fixef[nullVarsPos] <- as.numeric(nullfit$fixef) 
        fullm.con.list$init.HLfit$fixef[nullVarsPos] <- as.numeric(mean(nullfit$eta)) 
        if ("phi" %in% which.iterative) refullm.list$init.corrHLfit$phi <- NULL
        if ("lambda" %in% which.iterative) refullm.list$init.corrHLfit$lambda <- NULL
        if (trace) fullm.con.list$trace <- list(file="trace.refullfit.txt",append=T)
        refullfit <- do.call(method,fullm.con.list)$hlfit
        if ( refullfit$APHLs[[test.obj]] > fullfit$APHLs[[test.obj]] ) { ## but not necess better than nullfit
          conv <- conv + refullfit$APHLs[[test.obj]] - fullfit$APHLs[[test.obj]] ## only D(FULL)
          fullfit <- refullfit
          fullm.list <- update.ranef.pars(from.fit=refullfit,to.arglist=fullm.list)
          fullm.list <- update.ranef.pars(from.fit=refullfit,to.arglist=fullm.list,which.pars=c("lambda"))
          fullm.list$init.HLfit$fixef <- as.numeric(refullfit$fixef) 
          newfull <- T
          trace.info <- rbind(trace.info,data.frame(iter=locit,step="refullfit.3 (+)",obj=fullfit$APHLs[[test.obj]]))
        } 
        if ( fullfit$APHLs[[test.obj]] < nullfit$APHLs[[test.obj]]-0.001) { ## if still BAD
          ## serious problem with iterative algo in the simplest case
          ## profiling is our last chance
          prof.needed <- T ## => fullfit$APHLs[[test.obj]] < nullfit$APHLs[[test.obj]]-0.001 
          trace.info <- rbind(trace.info,data.frame(iter=locit,step="prof.needed",obj=fullfit$APHLs[[test.obj]]))
        }
      }
}     
      if (locit==1) LRTori <- 2*(fullfit$APHLs[[test.obj]]-nullfit$APHLs[[test.obj]]) ## first iteration only       
      if (newfull || newnull) {
        newforprof <- T ## T at initialization and after any 'unused' new full/null
      } else if (locit>1) trace.info <- rbind(trace.info,data.frame(iter=locit,step="no.new(f/n)ull",obj=fullfit$APHLs[[test.obj]]))
      if ( ! newforprof) {
        if (prof.needed) {
          print(trace.info)
          mess <- pastefrom("prof.needed && ! (newnull || newfull).")
          message(mess)
          ## conv must be null (else one the 'new...' is true) hence we will exit the loop...
        }
      } else if ( prof.needed || profiles>0 ) { ## given newforprof
        if ( ! prof.needed ) profiles <- profiles - 1 ## 'pays' one profile
        # then we need to maximize p_v(tau(beta)) over different beta values 
        # and therefore we need to be able to fit for fixed beta values 
        ## il faut (1) identifier les fixef sur lesquel le profile est fait
        prof.list <- nullm.list ## uses last null fit here
        if (trace) {
          try(unlink("trace.profile.txt")) ## the file is written in by HLCor()                   
          prof.list$trace <- list(file="trace.profile.txt",append=T)
        }
        profilize <- function(offset) {
          prof.list$predictor$offset <- profileX %*% offset ## semble marcher pour scalaire comme pour vecteur
          prof.fit <- do.call(method,prof.list)$hlfit ## 
          return(prof.fit$APHLs[[test.obj]])
        }
        if (length(profileVars)>1) {
          ## optim syntax
          ## uses last full fit here:
          lowup <- fullfit$fixef[profileVars] %*% t(c(-1,5)/2) ## should be of dim [,2]
          lowup <- apply(lowup,1,sort) ## should be of dim [2,]
          stop("optim call missing 'here' in corrMM.LRT")
          profmax <- optim(par=fullfit$fixef[profileVars],
                              fn=profilize,lower=lowup[1,],upper=lowup[2,],
                              control=list(fnscale=-1),method="L-BFGS-B")
          mess <- pastefrom("code needed here for optim output.")
          stop(mess)
        } else {
          ## one dim optimize syntax
          if (fullfit$fixef[profileVars] != 0) { ## if fullfit sufficiently different from null fit
            interval <- sort(fullfit$fixef[profileVars]*c(-1,5)/2)
          } else interval <- 0.001*c(-1,1) ##FR->FR very quick ad hoc patch
          profmax <- optimize(profilize,interval=interval,maximum=T)
          prof.list$predictor$offset <- profileX %*% profmax$maximum 
          proffit <- do.call(method,prof.list)$hlfit ## to recover the ranef parameters etc
          if ( proffit$APHLs[[test.obj]] > fullfit$APHLs[[test.obj]] ) { ## if profile max should improve fullfit
            trace.info <- rbind(trace.info,data.frame(iter=locit,step="proffit (+)",obj=proffit$APHLs[[test.obj]]))
            conv <- conv + proffit$APHLs[[test.obj]] - fullfit$APHLs[[test.obj]]
            refullm.list <- update.ranef.pars(from.fit=proffit,to.arglist=fullm.list)
            refullm.list <- update.ranef.pars(from.fit=proffit,to.arglist=refullm.list,which.pars=c("lambda","phi"))
            ## apparently no a good idea to provide init.HLfit$fixef without $v_h
#            refullm.list$init.HLfit$fixef[nullVarsPos] <- as.numeric(proffit$fixef) 
            refullm.list$init.HLfit$fixef[nullVarsPos] <- as.numeric(mean(proffit$eta)) 
            refullm.list$init.HLfit$fixef[profileVars] <- as.numeric(profmax$maximum) 
            refullm.list$init.HLfit$v_h <- proffit$v_h 
            if ("phi" %in% which.iterative) refullm.list$init.corrHLfit$phi <- NULL
            if ("lambda" %in% which.iterative) refullm.list$init.corrHLfit$lambda <- NULL
            refullfit <- do.call(method,refullm.list)$hlfit 
            if ( refullfit$APHLs[[test.obj]] > fullfit$APHLs[[test.obj]] ) { ## if confirmed improvement of fullfit 
              conv <- conv + refullfit$APHLs[[test.obj]] - fullfit$APHLs[[test.obj]]
              fullfit <- refullfit
              fullm.list <- refullm.list ## maybe not useful
              newfull <- T
              trace.info <- rbind(trace.info,data.frame(iter=locit,step="prof.refullfit (+)",obj=fullfit$APHLs[[test.obj]]))
            } else if ( refullfit$APHLs[[test.obj]] < proffit$APHLs[[test.obj]]-0.001 ) { ## if failure to confirm improvement of fullfit
              trace.info <- rbind(trace.info,data.frame(iter=locit,step="prof.refullfit (-)",obj=refullfit$APHLs[[test.obj]]))
              mess <- pastefrom("refullfit not as good as profile fit.")
              message(mess)
            } else { ## refull < full < proffit but all within 0.001 logL units 
              trace.info <- rbind(trace.info,data.frame(iter=locit,step="prof.refullfit (=)",obj=refullfit$APHLs[[test.obj]]))
            } 
          } else if ( proffit$APHLs[[test.obj]] < fullfit$APHLs[[test.obj]]-0.001 ) { ## somewhat problematic
            trace.info <- rbind(trace.info,data.frame(iter=locit,step="proffit (-)",obj=proffit$APHLs[[test.obj]]))
          } else {
            trace.info <- rbind(trace.info,data.frame(iter=locit,step="proffit (=)",obj=proffit$APHLs[[test.obj]]))
          }
        }
        newforprof <- F ## last newfull/newnull have been 'used'
        LRTprof <- 2*(fullfit$APHLs[[test.obj]]-nullfit$APHLs[[test.obj]])
      }
    } ## dispersion/correlation params were fit
    locit <- locit + 1
  } ## end loop
  ## BOOTSTRAP
  if ( ! is.na(testFix)) {
    if (boot.repl>0) {
      bootlist <- dotlist ## copies ranFix
      bootlist$control <- bootcontrol ## (a list)
      bootlist <- c(bootlist,list(null.predictor=null.predictor,null.disp=null.disp,REMLformula=REMLformula,method=method)) ## unchanged user REMLformula forwarded
      bootlist$verbose <- c(FALSE,FALSE)
      bootlist$trace <- FALSE 
      bootlist$boot.repl <- 0 ## avoids recursive call of bootstrap
      ## Give known simulated valus as starting values
      all.optim.vars <- c(which.optim.fit,which.iterative.fit)
      for (st in all.optim.vars) {
        if ( st %in% names(nullfit$corrPars)) {
          bootlist$init.corrHLfit[st] <- nullfit$corrPars[st]
        } else if ( st %in% names(nullfit)) {
          bootlist$init.corrHLfit[st] <- nullfit[st] ## handled in dotlist by corrHLfit, cf notes 090113
        } ## it's also possible that st is nowhere (10/2013: currently for Nugget)
      }
      bootreps<-matrix(,nrow=boot.repl,ncol=2) 
      if (testFix) {
        colnames(bootreps) <- c("full.p_v","null.p_v")
      } else colnames(bootreps) <- c("full.p_bv","null.p_bv")
      cat("bootstrap replicates: ")
      simbData <- nullfit$data
      if (tolower(nullfit$family$family)=="binomial") {
        form <- nullfit$predictor$oriFormula ## this must exists...  
        if (is.null(form)) {
          mess <- pastefrom("a processed binomial model formula must have an 'oriFormula' member.",prefix="(!) From ")
          stop(mess)
        }
      }
      ## the data contain any original variable not further used; e.g original random effect values in the simulation tests  
      for (ii in 1:boot.repl) {
        locitError <- 0
        repeat { ## for each ii!
          newy <- simulate(nullfit) ## cannot simulate all samples in one block since some may not be analyzable  
          if (tolower(nullfit$family$family)=="binomial") {
            simbData[[as.character(form[[2]][[2]])]] <- newy
            simbData[[as.character(form[[2]][[3]])]] <- nullfit$weights - newy    
          } else {simbData[[as.character(nullfit$predictor$formula[[2]])]] <- newy}
          bootlist$data <- simbData
          bootrepl <- try(do.call(corrMM.LRT,bootlist))
          if (class(bootrepl)[1] != "try-error") { ## eg separation in binomial models... alternatively, test it here (require full and null X.pv... )
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
    LRTinfo <- list(df=df,LRTprof = LRTprof,LRTori = LRTori)
    if (boot.repl>0) {
      meanbootLRT <- 2*mean(bootreps[,1]-bootreps[,2]) 
      LRTinfo$meanbootLRT <- meanbootLRT
      LRTinfo$bootreps <- bootreps
    }
  } else {
    resu <- list(fullfit=fullfit)
    LRTinfo <- list()
  }
  LRTinfo$trace.info <- trace.info 
##  resu$LRTinfo <- LRTinfo ## pas compatible avec hglmjob.R...
  resu <- c(resu,LRTinfo) ## loses the sublist structure, which wouldnot be compatible with hglmjob.R...  
  class(resu) <- c("fixedLRT",class(resu)) 
  return(resu)
}
