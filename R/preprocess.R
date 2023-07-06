.get_inits_preprocess_args <- local({
  inits_list <- list() # list of lists of default arguments for .preprocess calls(); expedites repeated calls of this fn.
  function(For) {
    if (is.null(inits_list[[For]])) {
      FP <- formals(.preprocess)
      nFP <- names(FP)
      inits <- FP # with .preprocess defaults (e.g., 'init')
      FHF <- formals(HLfit)
      nFHF <- names(FHF)
      if (For=="HLfit") {
        # push defaults from HLfit
        validnames <- intersect(nFP, nFHF)
        inits[validnames] <- FHF[validnames] 
      } else {
        # push defaults from HLfit or else HLCor
        FHC <- formals(HLCor)
        nFHC <- names(FHC)
        nHCspecific <- setdiff(nFHC, nFHF)
        validnames <- intersect(nFP, unique(c(nFHF,nFHC)))
        inits[validnames] <- c(FHF,FHC[nHCspecific])[validnames] 
        if (For=="is_separated") inits$family <- binomial() 
      } 
      inits$rand.families <- FHF$rand.family ## because preprocess expects $rand.families 
      inits$For <- For 
      inits_list[[For]] <<- inits
    } 
    return(inits_list[[For]])
  }
})

.checkRandLink <- function(rand.family) {
  lcrandfamfam <- tolower(rand.family$family) ## tolower once and for all
  oklink <- FALSE
  ## cases where g(u)=th(u)
  if (lcrandfamfam=="gaussian" && rand.family$link=="identity") oklink <- TRUE          
  if (lcrandfamfam=="gamma" && rand.family$link=="log") oklink <- TRUE          
  if (lcrandfamfam=="inverse.gamma" && rand.family$link=="-1/mu") oklink <- TRUE
  if (lcrandfamfam=="beta" && rand.family$link=="logit") oklink <- TRUE
  ## cases where g(u)!=th(u)
  if (lcrandfamfam=="inverse.gamma" && rand.family$link=="log") oklink <- TRUE 
  if (lcrandfamfam=="gamma" && rand.family$link=="identity") oklink <- TRUE ## gamma(identity)
  if ( ! oklink) {
    allowed <- switch(lcrandfamfam,
                      gaussian= "is 'identity'",
                      gamma= "is 'log' (explicit 'identity' may be *tried*)", ## gamma(identity) not yet working
                      #  and the above code protext against Gamma(inverse); it could perhaps be allowed for trying ? 
                      beta= "is 'logit'",
                      "inverse.gamma" = "are '-1/mu' and 'log'"
    )
    mess <- paste0("(!) rand.family/link combination not handled;\nallowed link(s) for rand.family '",
                   rand.family$family,"' ",allowed)
    stop(mess)
  }
  lcrandfamfam
}


.checkRandLinkS <- function(rand.families) {
  for (it in seq_len(length(rand.families))) {
    rf <- rand.families[[it]]
    if (is.character(rf)) {
      rand.families[[it]] <- switch(tolower(rf),
                  gaussian = gaussian(),
                  gamma = spaMM_Gamma(link=log), 
                  beta = Beta(), 
                  "inverse.gamma" = inverse.Gamma(),
                  stop("rand.family argument not valid"))
    }
  }
  lcrandfamfam <- unlist(lapply(rand.families, .checkRandLink)) ## a _vector_ of lcrandfamfam := tolower(rand.family$family)
  unique.psi_M <- sapply(lcrandfamfam, switch,
                         gaussian = 0,
                         gamma = 1, 
                         beta = 1/2, 
                         "inverse.gamma" = 1
  )
  names(unique.psi_M) <- NULL
  return(structure(rand.families,lcrandfamfam=lcrandfamfam,unique.psi_M=unique.psi_M))
}

.checkRespFam <- function(family, spaMM.=TRUE) {
  family <- tryCatch(family,error=function(e) e)
  if (inherits(family, "simpleError")) { # presumably 'mgcv::negbin()' ['mgcv::negbin' handled below]
    if ( family$message == "'theta' must be specified") {
      mess <- "spaMM::negbin is masked by mgcv::negbin. Unload mgcv, or use 'family=spaMM::negbin'."
      stop(mess)
      #family <- spaMM::negbin()
    }
  } 
  ## lines derived from glm().
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    envname <- environmentName(environment(family))
    if ( ! envname %in% c("stats","spaMM","")) { ## " typically occurs when family has already been checked...
      if (envname=="mgcv") {
        mess <- "spaMM::negbin is masked by mgcv::negbin. Unload  mgcv, or use 'family=spaMM::negbin'."
        stop(mess)
      } else message(paste("family from namespace environment",envname,"possibly not correctly handled"))
    }
    family <- family()
  }  
  if (is.language(family)) { # to bypass spaMM_Gamma by using family=quote(stats::Gamma(log))
    spaMM. <- (! length(grep("stats::",paste(family)[1L]))) ## spaMM. becomes FALSE if explicit quote(stats::...)
    family <- eval(family) 
  } 
  if (family$family=="Gamma" && spaMM.) {
    family <- spaMM_Gamma(link=family$link) 
  }
  return(family) ## input negbin(...) or COMPoisson(...) are evaluated without error and returned as is => nu/shape unaffected
}

.eval_v_h_bounds <- function(cum_n_u_h, rand.families) {
  ## builds box constraints either NULL or non-trivial, of length n_u_h
  n_u_h <- cum_n_u_h[length(cum_n_u_h)]
  lower.v_h <- rep(-Inf,n_u_h)
  boxConstraintsBool <- FALSE
  for (it in seq_along(rand.families)) {
    if (rand.families[[it]]$family=="Gamma" && rand.families[[it]]$link=="identity") { ## gamma(identity)
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      lower.v_h[u.range] <- 1e-6 ## 1e-12 is disastrous
      boxConstraintsBool <- TRUE
    }
  }
  if ( ! boxConstraintsBool ) lower.v_h <- NULL 
  boxConstraintsBool <- FALSE
  upper.v_h <- rep(Inf,n_u_h)
  for (it in seq_along(rand.families)) {
    if (rand.families[[it]]$family=="inverse.Gamma" && rand.families[[it]]$link=="-1/mu") { ## v ~ -Gamma
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      upper.v_h[u.range] <- -1e-06
      boxConstraintsBool <- TRUE
    } else if (rand.families[[it]]$family=="Gamma" && rand.families[[it]]$link=="log") {
      ## Gamma(log)$linkinv is pmax(exp(eta), .Machine$double.eps), ensuring that all gamma deviates are >= .Machine$double.eps
      ## we ensure that log(u_h) has symmetric bounds on log scale (redefine Gamma()$linkfun ?)
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      upper.v_h[u.range] <- Gamma(log)$linkfun(1/.Machine$double.eps)
      boxConstraintsBool <- TRUE
    } 
  }
  if ( ! boxConstraintsBool ) upper.v_h <- NULL
  return(list(lower.v_h=lower.v_h,upper.v_h=upper.v_h))
}

.def_u_h_v_h_from_v_h <- function(processed) {
  rand.families <- processed$rand.families
  lowup <- .eval_v_h_bounds(processed$cum_n_u_h, rand.families)
  cum_n_u_h <- processed$cum_n_u_h
  u_list <- vector("list", length(rand.families)) # to avoid reepated initializations within .u_h_v_h_from_v_h()
  .u_h_v_h_from_v_h <- function(v_h, lower.v_h=lowup$lower.v_h, upper.v_h=lowup$upper.v_h) {
    if(!is.null(lower.v_h)) {v_h[v_h<lower.v_h] <- lower.v_h}
    if(!is.null(upper.v_h)) {v_h[v_h>upper.v_h] <- upper.v_h}
    nrand <- length(rand.families)
    for (it in seq_len(nrand)) {
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      u_list[[it]] <- rand.families[[it]]$linkinv(v_h[u.range])
      if (any(is.infinite(u_list[[it]]))) {
        warning("infinite random values ('u_h') were constrained to finite range.") 
        u_list[[it]] <- pmin(.Machine$double.xmax, pmax(-.Machine$double.xmax,u_list[[it]]) )
      }
    }
    u_h <- .unlist(u_list)
    ## if there were box constr, v_h may have been modified, we put it in return value
    if ( ! (is.null(lower.v_h) && is.null(upper.v_h))) attr(u_h,"v_h") <- v_h
    return(u_h)
  }
  .u_h_v_h_from_v_h
}

.def_updateW_ranefS <- function(processed) {
  rand.families <- processed$rand.families
  nrand <- length(rand.families)
  w.ranef_list <- dlogWran_dv_h_list <- dvdu_list <- vector("list", nrand)
  cum_n_u_h <- processed$cum_n_u_h
  .updateW_ranefS <- function(u_h, v_h, lambda) {
    for (it in seq_len(nrand)) {
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      blob <- .updateWranef(rand.family=rand.families[[it]],lambda[u.range],u_h[u.range],v_h[u.range])
      w.ranef_list[[it]] <- blob$w.ranef
      dlogWran_dv_h_list[[it]] <- blob$dlogWran_dv_h
      dvdu_list[[it]] <- blob$dvdu
    }
    resu <- list(w.ranef=.unlist(w.ranef_list),dlogWran_dv_h=.unlist(dlogWran_dv_h_list),dvdu=.unlist(dvdu_list))
    ## the test is invalid for ranCoefs:
    # if (nrand==1L && rand.families[[1L]]$family=="gaussian") resu$unique_w.ranef <- w.ranef[[1L]] # used in sparsePrecision code
    #if (length(unique_w.ranef <- unique(w.ranef))==1L) resu$unique_w.ranef <- unique_w.ranef # used in sparsePrecision code
    resu
  }
  .updateW_ranefS
}


# to convert back a family object to a call, EXCEPT for "multi". 
# get_param should be TRUE only in cases we are sure that the param has a value
.as_call_family <- function(family,get_param=FALSE) {
  # checkRespFam should have been run. 
  # Then we have either a call or an evaluated family object
  if (identical(family[[1L]],"multi")) {
    return(family) ## FR->FR tempo fix
  } else {
    as_call <- environment(family$aic)$mc ## shouldn't be called on multi...
    if (is.null(as_call)) { ## then we have a stats::  family 
      as_call <- call(family$family, link=family$link) ## need a substitute() ? 
    } else if (get_param) { ## param values are constant or dynamically assigned to processed$family
      famfam <- family$family 
      if (famfam %in% c("negbin1","negbin2")) {
        as_call$shape <- environment(family$aic)$shape
      } else if (famfam %in% c("beta_resp","betabin")) {
        as_call$prec <- environment(family$aic)$prec
      } else if (famfam=="COMPoisson") {
        as_call$nu <- environment(family$aic)$nu
      }
    }
    return(as_call)
    ## avoids display of family def in fit results
  }
}

.calc_ZAlist <- function(Zlist, AMatrices) { 
  if ( length(Zlist) > 0L ) {
    if (length(AMatrices)) {
      Znames <- names(Zlist) ## explicit Znames necessary for prediction with re.form: Otherwise this should be OK:
      if (is.null(names(Zlist))) names(Zlist) <- Znames <- paste(seq_len(length(Zlist)))
      ## logic is Z[nresp,nUniqueRespLoc].A[nUniqueRespLoc,nHiddenv].L[nHiddenv,nHiddenv]
      for (char_rd in Znames) {
        Amatrix <- AMatrices[[char_rd]]
        if ( ! is.null(Amatrix)) {
          is_incid <- attr(Zlist[[char_rd]],"is_incid")
          if (inherits(Amatrix,"pMatrix")) {
            Amatrix <- as(as(Amatrix, "nMatrix"), "TsparseMatrix") # => ngTMatrix 
          } else if ( ! is.null(is_incid)) {
            if (is_incid) is_incid <- attr(Amatrix,"is_incid") # .spaMM_spde.make.A() provides this attr. Otherwise, may be NULL, in which case ./.
            # ./. a later correct message may occur ("'is_incid' attribute missing, which suggests inefficient code in .calc_new_X_ZAC().)
          } 
          ZAnames <- colnames(Zlist[[char_rd]])
          if ( ! setequal(rownames(Amatrix),ZAnames)) {
            mess <- paste0("Any 'A' matrix must have row names that match the levels of the random effects\n (",
                           paste0(ZAnames[1L:min(5L,length(ZAnames))], collapse=" "),if(length(ZAnames)>5L){"...)."} else{")."})
            stop(mess)
          }
          Zlist[[char_rd]] <- Zlist[[char_rd]] %*% Amatrix[ZAnames,] # for the IMRF model, Z must be identity
          rownames(Zlist[[char_rd]]) <- NULL
          attr(Zlist[[char_rd]],"is_incid") <- is_incid
        }
      }
      attr(Zlist,"AMatrices") <- AMatrices 
    } else attr(Zlist,"AMatrices") <- list() # to allow any later fn such as .assign_geoinfo.... to safely set named elements to it. 
  } 
  return(Zlist)
}

.assign_ZAfix <- function(processed) { # _F I X M E_ extend to all cases of ZAfix <- .ad_hoc_cbind(...) ?
  ZAfix <- .ad_hoc_cbind(processed$ZAlist, as_matrix=FALSE)   
  if (processed$is_spprec) {
    if ( ! inherits(ZAfix,"sparseMatrix")) {
      if (.spaMM.data$options$Matrix_old) { # this block appears to evade the long tests
        ZAfix <- as(ZAfix,"dgCMatrix") # .Rcpp_as_dgCMatrix(ZAfix) ## 
      } else ZAfix <- as(as(ZAfix,"generalMatrix"),"CsparseMatrix") # .Rcpp_as_dgCMatrix(ZAfix) ## 
    }
    processed$AUGI0_ZX$is_unitary_ZAfix <- FALSE
  }
  processed$AUGI0_ZX$ZAfix <- ZAfix
}

.evalWrapper <- function(object,
                         element, ## using try() has a visible impact on speed. Hence, the element must exist in the envir.
                         from=NULL) { ## exported for programming purposes
  if (  is.list(object)) {
    if (is.null(from)) from <- seq_len(length(object))
    if (length(from)>1L) {
      resu <- vector("list",length(from))
      for (id in seq_len(length(from))) {
        resu[[id]] <- eval(parse(text=element),envir=object[[id]])
      }
      names(resu) <- names(object[from])
      return(resu) ## a list, e.g a data list
    } else {
      object <- object[[from]]
    } 
  } 
  resu <- eval(parse(text=element),envir=object) ## try(eval(parse(text=element),envir=object),silent = TRUE)
  ## if (inherits(resu,"try-error")) resu <- NULL
  return(resu)
}

# .assignWrapper <- function(object,assignment) {
#   if (  is.list(object)) {
#     #    for (it in seq_len(length((object))) {eval(parse(text=paste0("object[[",it,"]]$",element," <- ",value)))} 
#     ## fails bc nam loses its enclosing "s :
#     #for (nam in names(object)) {eval(parse(text=paste0("object[[",nam,"]]$",element," <- ",value)))} 
#     for (nam in names(object)) {eval(parse(text=paste0("object[[\"",nam,"\"]]$",assignment)))} 
#   } else eval(parse(text=paste0("object$",assignment)))
#   ## no need to return the modified environment
# }
.assignWrapper <- function(object,assignment) { ## not using assign(), to allow eg. verbose['warn'] <- ... instead of simple variable name    # older version ".setProcessed" (<2.2.119) did not use de envir argument; this was tricky
  if (  is.list(object)) {
    for (it in seq_along(object)) eval(parse(text=assignment),envir=object[[it]]) 
  } else eval(parse(text=assignment),envir=object)
  ## no need to return the modified environment
}

.calc_fam_corrected_guess <- function(guess, For, processed, link_=NULL, trunc_=NULL, nrand=NULL) {
  if (is.null(link_)) link_ <- processed$family$link
  link_[link_=="loglambda"] <- "log" # For mv, link_ is actually several links
  if (is.null(nrand)) nrand <- length(processed$ZAlist)
  if ( any(link_ != "identity")) {
    if (For=="optim" || ## to avoid high initial values of lambda with spuriously high logLik by Laplace approx.
        processed$HL[1L]=="SEM") { 
      if (any(link_=="log")) { ## test of family, not rand.family... 
        fam_corrected_guess <- log(1.00001+guess/nrand) 
      } else {
        if (processed$bin_all_or_none) {
          maxinit <- 0.1 ##  a low init value is better even if final lambda estimate is high.
        } else maxinit <- 0.2 
        fam_corrected_guess <- min(guess, maxinit/nrand)
      }
    } else { ## iterative: allow larger values, but within some limits
      if (any(link_=="log")) { ## test of family, not rand.family...
        fam_corrected_guess <- log(1.00001+guess/nrand) 
      } else {
        maxinit <- 2 ## even for processed$bin_all_or_none, test_all does not support lower value; but this is all based on ~ nothing
        fam_corrected_guess <- min(guess, maxinit/nrand)
      }
    }
  } else { ## identity link
    if (For=="optim") {
      maxinit <- 2
    } else {
      maxinit <- Inf
    }
    fam_corrected_guess <- min(guess, maxinit/nrand) 
  }
  # if (any(trunc_)) fam_corrected_guess <- fam_corrected_guess/2 # quick patch but has drawbacks
  return(fam_corrected_guess)
}

.preprocess_valuesforNAs <-  function(it, lcrandfamfam, rand.families, init.lambda, rd) {
  if(lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="identity" && init.lambda==1) {
    adhoc <- 0.9999
  } else if(lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="log") {
    objfn <- function(lambda) {psigamma(1/lambda,1)-init.lambda}
    adhoc <- uniroot(objfn,interval=c(1e-8,1e8))$root
  } else if(lcrandfamfam[it]=="beta" && rand.families[[it]]$link=="logit") {
    #ad hoc approximation which should be quite sufficient; otherwise hypergeometric fns.
    objfn <- function(lambda) {8* lambda^2+3.2898*lambda/(1+lambda)-init.lambda}
    adhoc <- uniroot(objfn,interval=c(2.5e-6,1e8))$root
  } else if(lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="log") {
    ## this tries to controle var(v)<-init.lambda of v=log(u) by defining the right lambda for var(u)=lambda/(1-lambda)
    ## log (X~inverse.G) = - log (~Gamma), ie
    ##  log (X~inverse.G(shape=1+1/init.lambda,scale=1/init.lambda))~ - log (rgamma(shape=1+1/init.lambda,scale=1/init.lambda)
    ## where + log (~Gamma) has known mean (=...) and variance=trigamma(shape)
    ## (pi^2)/6=1.644934... is upper bound for trigamma(x->1) ie for lambda ->infty.
    ## however, there is something wrong with setting very large initial lambda values. Hence
    # if (init.lambda > 1.64491 ) { ## trigamma(1+1/100000) 
    #   adhoc <- 100000 
    if (init.lambda > 0.3949341 ) { ## trigamma(3) [for lambda=0.5]
      adhoc <- 0.5 
    } else { ## loer init.lambda -> lower lambda
      objfn <- function(lambda) {psigamma(1+1/lambda,1)-init.lambda}
      adhoc <- uniroot(objfn,interval=c(1e-8,1e8))$root
    }
    ## but the mean of v  is function of lambda and I should rather try to control the second moment of v
  } else if(lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="-1/mu") {
    adhoc <- (sqrt(1+4*init.lambda)-1)/2 # simple exact solution
  } else adhoc <- init.lambda
  adhoc
}

.preprocess_HLmethod <- function(HLmethod, family, lcrandfamfam) {
  ## conversion from user-friendly format to standard 'XX(...)' format
  ## first index is for (0) h / (1) p_v(h) / (2) p^s_v(h) ie whether h lik or marginal lik is used for fixed effect estimation
  ## Other components determine three options wrt to leverages, some for ReML correction, some for notEQL.
  ## ML() vs HL() determines whether the hatval leverages are computed ie whether some basic ReML correction in applied
  ## second index is for further ReML correction: no (0) or yes (1) D log Det / d log lambda correction (2) further p^s_bv correction
  ## third index is for use of (0) EQL deviance residuals (this affects the leverages) or not (1) (this is not an ReML correction... but impatcs only dispersion estimation..)
  ## thus overall we have <ReML/not>( <h/l> , <more ReML/not> , <not/EQL> )
  ## NohL07 table 1 has interesting terminology and further tables show even more cases
  # comme l'EQL est dans un monde quasi gaussien, il se ramene aux LMM et utilise les leverage standard,
  # Pour les GLMM non LMM, Dans la mesure ou HL(0,.) utilise h et non q+, ce n'est pas l'EQL pour les fixed params
  # Dans la mesure ou on utilise les leverages standard pour les dispersion param, c'est de l'EQL
  if (HLmethod=="ML") {
    HLmethod <- "ML(1,1,1)"  
  } else if (HLmethod=="SEM") {
    #if ( ! requireNamespace("probitgem",quietly = TRUE)) {## will pass CHECK if CRAN knows probitgem
    # if ( ! ("probitgem" %in% .packages()) ) { ## passed CRAN checks; appropriate sif probitgem explicitly loaded
    if ( ! ("probitgem" %in% .packages(all.available = TRUE)) ) { ## 
      stop("Package 'probitgem' not available for fitting by SEM.")
    } else do.call("require", list(package="probitgem", quietly = TRUE)) ## passes CRAN checks and transparent...
    #      eval(as.call(c(quote(require),list(package="probitgem", quietly = TRUE)))) ## passes CRAN checks (but temporary)
    #      eval(as.call(c(quote(requireNamespace),list(package="probitgem", quietly = TRUE)))) ## is not sufficient
    HLmethod <- "ML('SEM',NA,NA)" 
    if ( ! (family$family=="binomial" && family$link=="probit")) {
      stop("SEM is applicable only to binomial(probit) models.")
    }
  } else if (HLmethod=="REML") {
    HLmethod <- "HL(1,1,1)"  
  } else if (HLmethod=="REPQL" || HLmethod=="PQL") {
    if (any(lcrandfamfam!="gaussian"))  stop("PQL is not defined for HGLMs in general. Do you mean 'EQL-'?") 
    HLmethod <- "HL(0,0,1)" ## (=REPQL, equivalent to HL(0,0,0) ie EQL- for GLMMs )
  } else if (HLmethod=="PQL/L") { ## again no D log Det / d log lambda correction
    if (any(lcrandfamfam!="gaussian"))  stop("PQL is not defined for HGLMs in general. Do you mean 'EQL-'?") 
    HLmethod <- "ML(0,0,1)" ## (equivalent to ML(0,0,0) for GLMMs)
  } else if (HLmethod=="EQL-") { ## version LeeNP06 p.212 incomplete 1st order ## probably hglm package
    ## thus overall we have <ReML->HL >( <h->0> , <not more ReML->0> , <EQL -> 0> )
    HLmethod <- "HL(0,0,0)" ## 
  } else if (HLmethod=="EQL+") { ## version LeeN01 complete 1st order
    ## thus overall we have <ReML->HL >( <h->0> , <more ReML->1> , <EQL -> 0> )
    HLmethod <- "HL(0,1,0)" ## (0,...): gaussianise everything, hence no a(1) correction ## there is no HL(1) a(1) correction in GLM.MME
  }
  return(HLmethod)
}


.preprocess_SEMargs <- function(BinomialDen, nobs, control.HLfit, y) {
  if (sum(BinomialDen) != nobs) {
    stop("(!) SEM procedure: the data do not seem binary; non-binary data are not handled.")
    
    # repsucces <- rep(seq(nrow(currentSample)),currentSample$succes)
    # repechec <- rep(seq(nrow(currentSample)),currentSample$echec)
    # currentSample <- currentSample[c(repsucces,repechec),]
    # currentSample$echec <- rep(c(0,1),c(length(repsucces),length(repechec)))
    # currentSample$succes <- 1L-currentSample$echec
    
  }
  SEMseed <- control.HLfit$SEMseed ##  
  # if (is.null(SEMseed)) SEMseed <- 123 ## OK pour SEM *unique* mais remplace par NULL dans optimthroughSmooth
  SEMargs <- list(SEMseed=SEMseed)
  SEMargs$nMCint <- control.HLfit$nMCint ##  as is => SEM procedure must handle null value
  SEMargs$get_SEM_info <- control.HLfit$get_SEM_info ##  as is => SEM procedure must handle null value
  SEMargs$control_pmvnorm$maxpts <- control.HLfit$pmvnorm_maxpts ##  as is => SEM procedure must handle null value
  #SEMlogL <- control.HLfit$SEMlogL
  #if (is.null(SEMlogL)) SEMlogL <- "pMVN" ## FR->FR 1.8.25 !
  #SEMargs$SEMlogL <- SEMlogL
  SEMargs$SEMlogL <- control.HLfit$SEMlogL
  nSEMiter <- control.HLfit$nSEMiter ##  
  if ( (! is.null(nSEMiter)) && nSEMiter < 10) {
    stop(" 'nSEMiter' should be >9")
  } else SEMargs$nSEMiter <- nSEMiter ## default given by probitgem, typic 200 
  SEMargs$ngibbs <- control.HLfit$ngibbs ## default given by probitgem, typic 20 
  SEMargs$SEMsample <- control.HLfit$SEMsample ## stays NULL if NULL
  SEMargs$whichy1 <- (y==1) 
  return(SEMargs)
}

.is_link_canonical <- function(family) {
  if (is.null(canonicalLink <- family$flags$canonicalLink)) {
    canonicalLink <- FALSE
    if (family$family=="gaussian" && family$link=="identity") {
      canonicalLink <- TRUE
    } else if (family$family=="poisson" && family$link=="log") {
      canonicalLink <- TRUE
    } else if (family$family=="binomial" && family$link=="logit") {
      canonicalLink <- TRUE
    } else if (family$family=="Gamma" && family$link=="inverse") {
      canonicalLink <- TRUE
    } else if (family$family=="COMPoisson" && family$link=="loglambda") {
      canonicalLink <- TRUE
    } ## no implemented canonical link case for negbin
  }
  return(canonicalLink)
}

.is_LM <- function(family) {
  if (is.null(LMbool <- family$flags$LMbool)) {
    LMbool <- (family$family=="gaussian" && family$link=="identity")
  }
  return(LMbool)
}

.setattr_G_LMMbool <- function(models, processed, family=processed$family, lcrandfamfam=processed$lcrandfamfam, 
                               prior.weights=processed$prior.weights) {
  GLMMbool <- (length(lcrandfamfam) && all(lcrandfamfam=="gaussian") ) ## only allowed gaussian rand.family is gaussian(identity) 
  const_pw <- ( ! inherits(prior.weights,"call"))
  unit_GLMweights <- (
                     (family$family=="gaussian" && family$flags$canonicalLink ) ||
                     (family$family=="Gamma" && family$link=="log") 
                  )
  unit_Hobs_weights <- (
    (family$family=="gaussian" && family$flags$canonicalLink ) ||
      (family$family=="Gamma" && family$link=="log" && ! processed$how$obsInfo) 
  )
  const_Hobs_wresid <- unit_Hobs_weights && const_pw # constant non-unit GLM weights do not occur in actual families otherwise that case might need to be distinguished (cd muetafn).
  if (GLMMbool) { # const w.ranef
    GLGLLM_const_w <- const_Hobs_wresid # const w.ranef && [const w.resid: including Gamma(log) GLMM...]
    # GLGLLM_const_w controls whether the weights and augmented matrix need to be updated over iterations of IRLS.
    # it is not that the weights are constant across augmented 'levels' (hence const_pw is not determined by unique(pw))
    # Likewise it does not mean that phi is not reestimated between IRLSs. phi valeus are always constant within IRLS.
    LMMbool <- (family$family=="gaussian" && family$flags$canonicalLink ) 
    LLM_const_w <- (LMMbool && const_Hobs_wresid) 
  } else LMMbool <- LLM_const_w <- GLGLLM_const_w <- FALSE
  return(structure(models, LMMbool=LMMbool, GLMMbool=GLMMbool, LLM_const_w=LLM_const_w,  
                   GLGLLM_const_w=GLGLLM_const_w, const_Hobs_wresid=const_Hobs_wresid, 
                   unit_GLMweights=unit_GLMweights, # for .muetafn()
                   unit_Hobs_weights=unit_Hobs_weights) # for .vecdisneeded()
         )
}

.calc_rankinfo <- function(X.pv, verbose=FALSE, tol) {
  nc <- ncol(X.pv)
  nobs <- nrow(X.pv)
  tol <- eval(tol)
  if (inherits(X.pv,"sparseMatrix")) {
    if (nc>nobs) { ## qr does not handle "wide " matrices
      zerofill <- sparseMatrix(dims = c(nc-nobs,nc), i={}, j={})
      if (verbose) print("qr(rbind2(<sparseMatrix>,zerofill))")
      rankinfo <- qr(x=rbind2(X.pv,zerofill),tol=tol)
    } else {
      if (verbose) print("qr(<sparseMatrix>)")
      rankinfo <- qr(x=X.pv,tol=tol)
    }
  } else  {
    if ((nc^2)*(max(nc,nobs)-nc/3)>1e10) {
      message("The design matrix for fixed effects is quite large. Checking its rank will take time. 
              See help('rankinfo') for possible ways to circumvent this.")
    }
    ## .rankinfo request much memory, but is clearly faster for matrices of size ~ 5000*2900
    ## HOWEVER it does not always yield *p_bv* results consistent with qr !!!!!
    ## +
    ## FIXME this may evolve when  ColPivHouseholderQR implements "blocking strategies".
    ## Waiting for such improvements, spaMM.getOption("rankMethod") is "qr" by default.
    rankmethod <- spaMM.getOption("rankMethod")
    if (is.null(rankmethod)) {
      .rankinfo_usable <- (prod(dim(X.pv))<1e8)
      rankmethod <- c("qr",".rankinfo")[.rankinfo_usable+1L] 
    }
    if (verbose)  print(paste0(rankmethod,"(<matrix>)"))
    rankinfo <- do.call(rankmethod,list(x=X.pv,tol=tol))
  }
  if (inherits(rankinfo,"qr")) {
    rankinfo <- list(rank=rankinfo$rank, 
                     whichcols=sort(rankinfo$pivot[1:rankinfo$rank]), ## cf lm.fit or base::qr.coef. sort() is cosmetic, \neq sort.list() !
                     method="qr") 
  } else if (inherits(rankinfo,"sparseQR")) { ## result of calling qr() when inherits(X.pv,"sparseMatrix")
    checkd <- (abs(diag(rankinfo@R))>tol) ## 1:rankinfo$rank in base::qr; but here the zero-filled rows are not the last of R
    whichcols <- (rankinfo@q+1L)[checkd] ## $pivot[ <checkd> ] in base::qr
    rankinfo <- list(rank=sum(checkd), whichcols=sort((rankinfo@q+1L)[checkd]), method="sparseQR")
  } else if (rankinfo$method==".rankinfo") { ## .rankinfo:
    checkd <- seq_len(rankinfo$rank) ## not quite clear  
    piv <- (rankinfo$pivot+1L)[checkd] ## +1L bc .rankinfo()$pivot indexes columns from 0
    rankinfo <- list(rank=rankinfo$rank, whichcols= sort(piv), # "not sure about correct code here", #  
                     method=".rankinfo")
  } else stop("Unknown rankinfo$method.")
  return(rankinfo)
}

.provide_CHMfactor <- function(corrMatrix, CHMtemplate) {
  if (inherits(corrMatrix,"dist")) { corrMatrix <- proxy::as.matrix(corrMatrix, diag=1) }
  if ( ! inherits(corrMatrix,"sparseMatrix")) { corrMatrix <- as(corrMatrix, "sparseMatrix") }
  if ( ! inherits(corrMatrix,"dsCMatrix")) { corrMatrix <- forceSymmetric(corrMatrix) } # essential for _correct_ .updateCHMfactor() call
  if (is.null(CHMtemplate)) { 
     Matrix::Cholesky(corrMatrix, perm= TRUE, LDL=FALSE)
  } else{ 
    Matrix::.updateCHMfactor(CHMtemplate, parent= corrMatrix, mult=0)
  }   
}
# =>
# Lunique can be deduced (as seen in geo_info code)
# 'solve(corrMatrix)'would be computed as 
# Linv <- Matrix::solve(corr_CHMfactor, as(corr_CHMfactor,"pMatrix"), system="L")
# precision <- .crossprod(Linv) # ~ solve(corrMatrix)


.get_corr_prec_from_covStruct <- function(covStruct,it, required) { 
  # this comes after .preprocess_covStruct() so covStruct must already be a list
  # Look whether [it] is a sublist with name element named "corrMatrix", with corr_prec_blob being the element, not the sublist!
  # If not, then look whether [it] is a sublist with name element named "precision" ... ditto
  if (is.null(corrMatrix <- covStruct[it][["corrMatrix"]])) {
    if (is.null(prec_blob <- covStruct[it][["precision"]])) {
      if (required) stop("missing covariance structure for corrMatrix model") ## no info in any form
    } else {
      ## else either the user provided covStruct$precision, .preprocess_covStruct()'ed (=>  reformatting no longer necessary ?)
      ## or spaMM produced precision from corrMat (=> already a list inheriting from "precision" (??))
      # => _F I X M E_ is next line necess ?
      if ( ! inherits(prec_blob,"precision"))  prec_blob <- structure(list(matrix=prec_blob),
                                                                      class=c("list","precision")) 
      return(prec_blob) # a list inheriting from "precision"
    }
  } else return(corrMatrix) # a matrix
}

.get_adjMatrix_from_covStruct <- function(covStruct,it, required=TRUE) { ## compatible with spaMM3.0 extended syntax
  if (length(covStruct)>1L) {
    adjMatrix <- covStruct[it][["adjMatrix"]]
  } else {
    adjMatrix <- covStruct[["adjMatrix"]]
  }
  if (required && is.null(adjMatrix)) stop("missing 'adjMatrix' for adjacency model") ## or SAR_WWt model...
  return(adjMatrix)
}

.sym_checked <- function(mMatrix, mname="matrix") { ## should return a dsCMatrix; user input may be a matrix
  if (ncol(mMatrix)!=nrow(mMatrix)) {stop(paste0("Some ",mname," that should be symmetric is found to be not even square."))}
  if ( ! identical(colnames(mMatrix), rownames(mMatrix))) {
    warning("Forcing colnames(mMatrix) <- rownames(mMatrix) before calling isSymmetric().")
    colnames(mMatrix) <- rownames(mMatrix)
  }
  if (isSymmetric(mMatrix)) {
    return(forceSymmetric(drop0(mMatrix))) # dsCMatrix! spaMM:::.sym_checked(diag(3)) is dsCMatrix
  } else {
    stop(paste0("Some ",mname," appears not to be symmetric."))
  }
}

.assign_cov_matrices__from_covStruct <- function(corr_info, covStruct=NULL, corrMatrix=NULL, adjMatrix=NULL,
                                                 required) {
  if ( ! is.null(covStruct)) covStruct <- .preprocess_covStruct(covStruct)
  corr_info$AMatrices <- attr(covStruct,"AMatrices")
  corr_types <- corr_info$corr_types
  namedlist <- structure(vector("list",length(corr_types)), names=seq_along(corr_types))
  if (is.null(corr_info$adjMatrices)) corr_info$adjMatrices <- namedlist
  if (is.null(corr_info$corrMatrices)) corr_info$corrMatrices <- namedlist
  # if (is.null(corr_info$corr_families)) corr_info$corr_families <- namedlist
  for (it in seq_along(corr_types)) {
    corr_type <- corr_types[[it]]
    if ( ! is.na(corr_type)) {
      if (corr_type=="adjacency" || corr_type=="SAR_WWt") {
        if ( is.null(adjMatrix) ) adjMatrix <- .get_adjMatrix_from_covStruct(covStruct,it, required=required)
        if (required) {
          nc <- ncol(adjMatrix)
          dsCdiag <- .symDiagonal(nc, x = rep.int(1,nc), uplo = "U",   kind="d")
          corr_info$adjMatrices[[it]] <- structure(.sym_checked(adjMatrix,"adjMatrix"), # dsCMatrix
                                                   dsCdiag=dsCdiag)
        }
      } else if (corr_type=="corrMatrix") {
        if (is.null(corrMatrix)) corrMatrix <- .get_corr_prec_from_covStruct(covStruct,it, required=required) 
        if ( (is.matrix(corrMatrix) || inherits(corrMatrix,"Matrix")) && 
             # : need to exclude "dist" and "precision" objects
             # (I defined a dim.precision() method) so dim() would not exclude precision
             .calc_denseness(sparseCorr <- drop0(corrMatrix), relative=TRUE) < 0.15) corrMatrix <- sparseCorr
        if ( ! is.null(corrMatrix)) corr_info$corrMatrices[[it]] <- corrMatrix 
        if (required) .check_corrMatrix(corr_info$corrMatrices[[it]], element=1) 
      } # else if (corr_type=="corrFamily")  # handled later.
      # IMRF AMatrices are assigned later from Zlist info, not from covStruct info
    }
  }
}

.calc_terms_heuristic_denseness <- function(terms_info, fixef_off_terms=terms_info$fixef_off_terms,
                                            fixef_levels=terms_info$fixef_levels) {
  vars_terms_table <- attr(fixef_off_terms,"factors") ## ""factors"" is a misleading name as table includes quantitative predictors
  if (length(vars_terms_table)) {
    n_levels <- sapply(fixef_levels,length)
    terms_heuristic_denseness <- rep(1,ncol(vars_terms_table)) 
    names(terms_heuristic_denseness) <- colnames(vars_terms_table)
    for (term in colnames(vars_terms_table)) { ## the cols of vars_terms_table should match the terms in the order used by attr(X.pv,"assign")...
      vars_in_term <- names(which(vars_terms_table[ ,term]>0))
      factors_in_terms <- intersect(names(fixef_levels),vars_in_term)
      if (any(vars_terms_table[ ,term]>0)) terms_heuristic_denseness[term] <- 1/prod(n_levels[factors_in_terms])
    }
  } else terms_heuristic_denseness <- NA
  return(terms_heuristic_denseness)
}

# Not often TRUE; test-cloglog is a case.
.determine_sparse_X <- function(terms_info, X.pv=terms_info$X) {
  sparse_X <- spaMM.getOption("sparse_X") 
  ## forcing sparse_X may (1) be slow for small problems 
  ## (2) entails the use of Matrix::Cholesky, which is less accurate => small bu visible effect on predVar in singular 'twolambda' case
  if (is.null(sparse_X)) {
    asgn <- attr(X.pv,"assign") ## "for each column in the matrix ... the term in the formula which gave rise to the column"
    if (any(asgn>0L) &&  length(terms_info$fixef_levels)) { # tests asgn bc terms_info$fixef_levels may be spuriously non-empty for terms that are added and removed within the formula
      col_heuristic_denseness <- rep(1,ncol(X.pv))
      terms_heuristic_denseness <- .calc_terms_heuristic_denseness(terms_info)
      for (it in seq_along(asgn)) if (asgn[it]>0L) col_heuristic_denseness[it] <- terms_heuristic_denseness[asgn[it]]
      sparse_X <- (mean(col_heuristic_denseness)<0.11) ## FIXME not enough tests of threshold; could use data.test in test-predVar which has mean density=0.19
    } else sparse_X <- FALSE
  }
  return(sparse_X)
}


.process_corr_info_spprec <- function(corr_info, For, sparse_precision) {
  corr_types <- corr_info$corr_types
  for (it in seq_along(corr_types)) {
    corr_type <- corr_types[[it]]
    if ( ! is.na(corr_type)) {
      if (corr_type=="corrMatrix" && sparse_precision) {
        if (! inherits(corr_info$corrMatrices[[it]],"precision")) {
          message("Conversion of corrMatrix to precision matrix")
          corr_info$corrMatrices[[it]] <- as_precision(corr_info$corrMatrices[[it]]) ## test by test-Matern-spprec
        }
      } else if (corr_type=="adjacency" || corr_type=="SAR_WWt") {
        if (For %in% c("fitme","corrHLfit")) {
          ## This block cannot be in .assign_cov_matrices__from_covStruct bc it requires sparse_precision value
          # We will need the $d to determine rhorange.
          # given the decomp is computed, we try to make it available for other computations.
          # But others computs should not assume it is available.
          # HLCor_body can also compute "symSVD" attribute and add in each processed environment
          decomp <- .provide_AR_factorization_info(corr_info$adjMatrices[[it]], sparse_precision, corr_type)
          if (corr_type=="SAR_WWt") attr(corr_info$adjMatrices[[it]],"UDU.") <- decomp
          if (corr_type=="adjacency") {
            attr(corr_info$adjMatrices[[it]],"symSVD") <- decomp
          } 
        }
      }
    }
  }
}


.preprocess_control.dist <- function(control.dist,corr_types) {
  # handles either a list with possible elements "dist.method","rho.mapping", or a nested list with elements for the different ranefs
  if ( ! is.null(control.dist)) {
    if (length(control.dist)) { # list with elements (or NULL_s_)
      matchnames <- intersect (names(control.dist), c("dist.method","rho.mapping"))
      if (length(matchnames)) {
        control_dist <- vector("list",length(corr_types))
        for (rd in seq_along(corr_types)) {
          if (corr_types[[rd]] %in% c("Matern","Cauchy","AR1", "IMRF")) { ## not (yet) clearly useful for AR1 ?
            control_dist[[as.character(rd)]] <- control.dist[matchnames]
          }
        }
        return(control_dist) # return nested list with elements for the different ranefs built from the 'matchnames'
      } else return(control.dist) ## No 'matchnames': control.dist should then be a nested list with elements for the different ranefs 
    } else return(vector("list", length(corr_types))) # Return "empty nested list" built from empty list() (could return NULL ?)
  } else return(NULL)
}


.preprocess_augZXy <- function(processed, init, ranFix) { 
  augZXy_cond_inner <- augZXy_cond <- spaMM.getOption("allow_augZXy")
  if (is.null(augZXy_cond)) {
    augZXy_cond <- (processed$models[["phi"]]=="phiScal") ## allows prior weights, but there must be a single phi scaling factor
    # => default is FALSE if phi is fixed but can be overcome by allow_augZXy=TRUE for devel purpose at least
    augZXy_cond_inner <- TRUE ## for .makeCovEst1() ## I have not thought about ths case with [phi fixed but overcome by allow_augZXy=TRUE]
  }
  # Conditions common to outer and inner optim
  if (augZXy_cond_inner) augZXy_cond_inner <- (is.null(processed$intervalInfo)) 
  if (augZXy_cond_inner) augZXy_cond_inner <-   attr(processed$models,"LMMbool")
  if (augZXy_cond_inner) augZXy_cond_inner <- ( is.null(processed$X.Re) || ! ncol(processed$X.Re)) ## exclude non-standard REML (avoiding NCOL(NULL)=1)
  if (augZXy_cond_inner) augZXy_cond_inner <- is.null(processed$X_off_fn) # not outer beta estimation
  augZXy_cond <- (augZXy_cond_inner && augZXy_cond) 
  # Conditions specific to outer optim
  if (augZXy_cond) augZXy_cond <- (processed$For=="fitme")
  if (augZXy_cond) augZXy_cond <- (length(processed$init_HLfit)==0L) ## excludes (among others) internal rho estimation
  # Fixed lambda must be excluded (bc then augZXy would fit a model with another lambda that the declared one)
  # And then we should not forget to exclude (partially)-fixed lambda in ranCoefs, and fixed hy_lam in multIMRF.
  # Exclude fixed lambda:  
  if (augZXy_cond) augZXy_cond <- all(is.na(processed$lambda.Fix))
  # then we also wish (as discussed in .calc_optim_args()) to exclude *forced* inner estimation of lambde, that is, NaN in init$lambda.
  if (augZXy_cond) augZXy_cond <- ! any(is.nan(init$lambda)) 
  # Exclude (partially)-fixed lambda in ranCoefs: 
  if (augZXy_cond) augZXy_cond <- (! any(processed$ranCoefs_blob$is_set)) ## from is_set at preprocessing stage, this excludes cases with ranCoefs set by user
  # There is no equivalent for corrFamily bc (even if daigonal elements are partially fixed) there is a single lambda which is fixed or not but not 'partially'.
  if (augZXy_cond && length(partiallyfixed <- ranFix$ranCoefs)) { #
    # we reach here if NO ranCoefs_blob$is_set => if there is a partially fixed ranCoef, which must be in fixed$ranCoefs
    # We then look to exclude cases where lambda's are fixed
    Xi_cols <- attr(processed$ZAlist,'Xi_cols')
    for (char_rd in names(partiallyfixed)) if (augZXy_cond) {
      Xi_ncol <- Xi_cols[as.integer(char_rd)]
      lampos <- rev(Xi_ncol*(Xi_ncol+1L)/2L -cumsum(seq(Xi_ncol))+1L)  ## NOT cumsum(seq(Xi_cols))
      augZXy_cond <- all(is.na(partiallyfixed[[char_rd]][lampos]))
    }
  }
  # Exclude fixed hy_lam in multIMRF:
  if (augZXy_cond) augZXy_cond <- ! any(grepl("hy_lam", names(unlist(ranFix$hyper)), fixed=TRUE))
  #
  #
  if (augZXy_cond == 1L) augZXy_cond <- ( ! is.call(processed$prior.weights) && attr(processed$prior.weights,"unique")) ## cf APHLs_by_augZXy code
  # which means that allow_augZXy=2L allows augZXy usage with non-constant prior weights 
  if (augZXy_cond)  {
    ######### exclude case where there is an init lambda / hy_lam / ranCoefs[<lambda positions>] 
    # bc it does not have the expected meaning when augZXY is used and it not clear how to correct that.
    # Poor init ranCoefs are quite a problem (refitting fitme6 with get_inits_from_fit(fitme6)$init["ranCoefs"]) would be quite poor by augZXy)
    #
    # For effect of phi controls, the problem is quite different (we could imagine using the augZXy code for given phi within HLfit)
    # But we decide... (see see .calc_optim_args() for handling of phi):
    augZXy_cond <- is.null(init$phi) 
    #
    if (augZXy_cond) augZXy_cond <- is.null(init$lambda) 
    if (augZXy_cond)  {
      if ( ! is.null(ini_rC <- init$ranCoefs)) { # We then look to exclude cases where lambda's have inits
        Xi_cols <- attr(processed$ZAlist,'Xi_cols')
        for (char_rd in names(ini_rC)) if (augZXy_cond) {
          Xi_ncol <- Xi_cols[as.integer(char_rd)]
          lampos <- rev(Xi_ncol*(Xi_ncol+1L)/2L -cumsum(seq(Xi_ncol))+1L)  ## NOT cumsum(seq(Xi_cols))
          augZXy_cond <- all(is.na(ini_rC[[char_rd]][lampos]))
        }
      }  
    } 
    #
    if (augZXy_cond) augZXy_cond <- ! any(grepl("hy_lam", names(unlist(init$hyper)), fixed=TRUE))
    #  Conditional on verbose bc otherwise various bootstraps will repeat this message:
    if ( ! augZXy_cond && processed$verbose["TRACE"]) message("Providing initial values for some variance parameter(s)\n   prevents use of an efficient algorithm to fit LMMs.")
    # : that must come after all other conditions otherwise the message is issued despite other reasons for not using augZXy.
  }
  processed$augZXy_cond <- structure(augZXy_cond, inner=augZXy_cond_inner)
  if (augZXy_cond) {
    processed$augZXy_env <- new.env(parent=emptyenv()) 
    processed$augZXy_env$objective <- -Inf 
  }
}

.reformat_resid_model <- function(resid.model
                                  #  check_old_syntax for back compatibility,
                                  # to check whether control.HLfit$resid.family was used, when resid.model is only a formula
                                  # This use of control.HLfit is no longer documented
                                  ,check_old_syntax=NULL) { 
  fixed <- as.list(resid.model$fixed) ## converts NULL to list() as exp'd for 'fixed' in fitme_body()
  if ( ! is.null(resid.model)) { 
    if ( ! is.null(form <- resid.model$formula)) {
      resid.model$formula <- .preprocess_formula(form)
    } else { ## including resid.model=~1
      resid.model <- list(formula=.preprocess_formula(resid.model), ## otherwise it retains the local environment of the fn in which it is match.call()ed!
                          family=check_old_syntax) # may be NULL, in which case it is later processed as spaMM_Gamma(log)
    }
    # control of phi matters for phiHGLM but the phi model has not yet been determined
    if (is.null(fixed[["phi"]])) {
      fixed[["phi"]] <- c(default=1)
      #        message("'phi' of residual dispersion model set to 1 by default") ## inappropriate when resid.model=~1
    } else if (is.na(fixed[["phi"]])) fixed[["phi"]] <- NULL ## to force estimation of this phi; 
    resid.model$fixed <- fixed
    if (is.null(resid.model$resid.model)) resid.model$resid.model <- .preprocess_formula(~1)
  } 
  resid.family <- resid.model$family
  if (is.null(resid.family)) {
    resid.family <- spaMM_Gamma(log)
  } else {
    if (resid.family$family!= "Gamma") stop("resid.family must be Gamma.")## spaMM_Gamma also valid by this test
    resid.family <- spaMM_Gamma(resid.family$link) ## we will need the returned link, not the promise 
  } ## and "quoted" suitable for saving in HLfit return object
  attr(resid.family,"quoted") <- substitute(spaMM_Gamma(link),
                                            list(link=resid.family$link))
  resid.model$family <- resid.family ## resid.predictor will also be returned
  if (is.null(resid.model$rand.family)) resid.model$rand.family <- gaussian() # apparently to avoid it to be NULL in .preprocess_phi_model() call
  return(resid.model)
}

.rankTrim <- function(X.pv, rankinfo, verbose=FALSE) {  # tests in /test-rank.R
  if (verbose)  str(rankinfo)
  if (rankinfo$rank < ncol(X.pv)) {   
    X.pv <- structure(X.pv[,rankinfo$whichcols,drop=FALSE], namesOri=colnames(X.pv),
                      assign=attr(X.pv, "assign")[rankinfo$whichcols], # that's the trimmed value that what lmerTest:::term2colX expect
                      assignOri=attr(X.pv, "assign") # this one only for .anova.glm() ...
                      ) 
    # etaFix$beta |         variables 
    #             |  valid vars | invalid vars
    #     (1)           (2)           (3)
    # (2): colnames(<HLfit>$envir$beta_cov_info$beta_cov) = colnames (<HLfit>$X.pv)
    # (1+2+3): namesOri, names(<HLfit>$fixef)
  } else {
    attr(X.pv,"namesOri") <- colnames(X.pv)
  }  
  attr(X.pv,"rankinfo") <- rankinfo
  return(X.pv)
}

.post_process_X <- function(X.pv, HL, 
                            rankinfo=NULL, # NULL default => compute it locally and trim matrix
                                           # <pre-existing list> => use it to trim matrix
                                           # FALSE => don't trim matrix, in .merge_processed() for mv fits
                            sparse_X) {
  Xattr <- attributes(X.pv)
  if ( ncol(X.pv)) {
    if (sparse_X) { ## sparse_X is useful for rankinfo bc Matrix::qr can be much faster
      if (.spaMM.data$options$Matrix_old) { # ugly... but such versions do not handle as(, "dMatrix"))
        X.pv <- as(X.pv,"dgCMatrix") # .Rcpp_as_dgCMatrix(X.pv) # 
      } else X.pv <- as(as(X.pv,"generalMatrix"),"CsparseMatrix") # .Rcpp_as_dgCMatrix(X.pv) #
    }
  } 
  if (ncol(X.pv)) {
    if (is.null(rankinfo)) rankinfo <- .calc_rankinfo(X.pv, tol=spaMM.getOption("rankTolerance")) # always return a list
    if (is.list(rankinfo)) {
      X.pv <- .rankTrim(X.pv,rankinfo = rankinfo)
    } else attr(X.pv,"namesOri") <- colnames(X.pv) # input 'rankinfo' was FALSE.
    # the unmodified "assign" attribute refers to the cols of the untrimmed matrix.
  }
  names_lostattrs <- setdiff(names(Xattr), c(names(attributes(X.pv)),"dim","dimnames"))
  attributes(X.pv)[names_lostattrs] <- Xattr[names_lostattrs] # as in .subcol_wAttr(). 
  return(X.pv)
}

.preprocess_formula <- function(formula, control.HLfit=NULL, ...) {
  if (inherits(formula,"predictor")) { 
    return(formula) ## happens eg in confint
    # stop("Do not call '.preprocess_formula' on a predictor object.")
  }
  env <- control.HLfit$formula_env 
  ## By partial matching, if control.HLfit was absent from the parent call,
  ## 'control' is matched to the present control.HLfit. Then control$formula_env will be effective
  ## (interesting feature here as control can hold all the info from everything from the control.[...] args.
  if (is.null(env)) env <- new.env() # emptyenv() does not work as base functions must be accessible...
  formlen <- length(formula)
  formula[[formlen]] <- .expandDoubleVerts(formula[[formlen]])
  rhs <- .expand_multIMRFs(formula[[formlen]]) 
  hyper_info <- attr(rhs,"hyper_info")
  if ( !is.null(hyper_info)) {
    hyper_form <- formula 
    environment(hyper_form) <- env
    hyper_info$formula <- hyper_form
    attr(rhs,"hyper_info") <- NULL
  }
  formula[[formlen]] <- rhs
  environment(formula) <- env # zut <- y~x; rezut <- spaMM:::.preprocess_formula(zut); ls(environment(zut)) confirms that the original 'formula's environment is unaffected
  predictor <- structure(formula,hyper_info=hyper_info) 
  class(predictor) <- c("predictor",class(formula))
  return(predictor)
}

.calc_vec_normIMRF <- function(exp_ranef_terms, corr_info, 
                               corr_types=corr_info$corr_types, corr_families=corr_info$corr_families) {
  nrand <- length(exp_ranef_terms)
  vec_normIMRF <- rep(FALSE, nrand)
  for (rd in seq_along(exp_ranef_terms)) { 
    corr_type <- corr_types[[rd]]
    if (! is.na(corr_type)) {
      if (corr_type== "IMRF" ) {
        useNorm <- attr(attr(exp_ranef_terms[[rd]],"type"),"pars")$no 
        if (is.null(useNorm)) {
          stop("is.null(useNorm)")  # DEBUGGING
        } else vec_normIMRF[rd] <- useNorm
      } else if (corr_type== "corrFamily" ) {
        useNorm <- corr_families[[1]]$normIMRF
        if ( ! is.null(useNorm)) {
          vec_normIMRF[rd] <- useNorm
        } # else leave it FALSE 
      } 
    }
  }
  return(vec_normIMRF)
}

.add_ZAfix_info <- function(AUGI0_ZX, ZAlist, sparse_precision, as_mat) {
  if ( ! any(AUGI0_ZX$vec_normIMRF)) {
    ZAfix <- .ad_hoc_cbind(ZAlist, as_matrix=FALSE)  
    if (sparse_precision) {
      if ( ! inherits(ZAfix,"sparseMatrix"))  {
        if (.spaMM.data$options$Matrix_old) { # this block appears to evade the long tests
          ZAfix <- as(ZAfix,"dgCMatrix") # .Rcpp_as_dgCMatrix(ZAfix) ## 
        } else ZAfix <- as(as(ZAfix,"generalMatrix"),"CsparseMatrix") # .Rcpp_as_dgCMatrix(ZAfix) ## 
      }
      rsZA <- rowSums(ZAfix) ## test that there a '1' per row and '0's otherwise:  
      AUGI0_ZX$is_unitary_ZAfix <- (all(unique(rsZA)==1L) && all(rowSums(ZAfix^2)==rsZA)) ## $ rather than attribute to S4 ZAfix
    } else {
      if (as_mat) ZAfix <- as.matrix(ZAfix)  
    }
    AUGI0_ZX$ZAfix <- ZAfix # Used in code for ZAL in HLfit_body(); and extensively to fit  by spprec. 
                            # Special care is required when A is later modified (as by .calc_normalized_ZAlist())
  }
  return(AUGI0_ZX)
}

.assign_X.Re_objective <- local({
  obj_warned <- FALSE
  function(processed, 
           XReinput, # used only if  ! is.null(REMLformula)  && identical(attr(REMLformula,"isML"),TRUE)
           REMLformula, data, 
           X.pv, # used only if  ! is.null(REMLformula) !! may be different from XReinput: seek cases where 'keepInREML' is TRUE 
           objective) {
    if ( ! is.null(REMLformula) ) { # ML or non-standard REML... 
      if (identical(attr(REMLformula,"isML"),TRUE)) {
        processed$X.Re <- matrix(ncol=0, nrow=nrow(XReinput))
      } else {# non-standard REML... 
        if (length(REMLformula)==3L) REMLformula <- REMLformula[-2] # no need to handle any complicated LHS, which moreover is pb for reusing a previous terms_info
        REMLFrames <- .get_terms_info(formula=REMLformula, data=data, famfam="") ## design matrix X, offset...
        # where we don't need famfam since we don't use REMLFrames$Y
        X.Re <- REMLFrames$X
        if (ncol(X.Re)) { 
          # wAugX will have lost its colnames...
          unrestricting_cols <- which(colnames(X.pv) %in% setdiff(colnames(X.pv),colnames(X.Re))) ## not in X.Re
          extra_vars <- setdiff(colnames(X.Re),colnames(X.pv)) ## example("update") tests this. # should occur too with etaFix + REMLformula
          distinct.X.ReML <- c(length(unrestricting_cols), length(extra_vars)) ## TWO integers often used as booleans 
          attr(X.Re,"distinct.X.ReML") <- distinct.X.ReML 
          if (any(distinct.X.ReML)){ # *subcase* of non-standard REML (other subcase when formula= explicit REMLformula... with ranFix...)
            if (attr(X.Re,"distinct.X.ReML")[1L]) attr(X.Re,"unrestricting_cols") <- unrestricting_cols # cols of X.pv not in X.Re
            if (attr(X.Re,"distinct.X.ReML")[2L]) attr(X.Re,"extra_vars") <- extra_vars # cols of X.Re not in X.pv ## example("update") tests this.
          } 
        } else message("Possibly inefficient code: REMLformula not recognized as representing an REML specification.") 
        processed$X.Re <- X.Re 
      }
    } ## else this should be standard REML; let effectively processed$X.Re <- NULL 
    # if standard ML: there is an REMLformula ~ 0 (or with ranefs ?); local X.Re and processed$X.Re is 0-col matrix
    # if standard REML: REMLformula is NULL: $X.Re is X.pv, processed$X.Re is NULL
    # non standard REML: other REMLformula: $X.Re and processed$X.Re identical, and may take essentially any value
    # So a testing pattern is if ( is.null(X.Re)) {<REML standard>} else if ( ncol(X.Re)==0L) {<ML standard>} else {<non-standard REML>}
    # NCOL(.$X.Re) identifies all REML cases
    if (is.null(objective)) {
      if (NCOL(processed$X.Re)) { ## standard or non-standard REML
        processed$objective <- "p_bv"  ## info for fitme_body and corrHLfit_body, HLfit may return_only="p_bvAPHLs" but use $objective in logL_tol convergence test
      } else processed$objective <- "p_v"
    } else {
      if ( ! obj_warned) {
        warning("Non-NULL 'objective' is deprecated except for development purposes.", immediate. = TRUE)
        obj_warned <<- TRUE
      }
      processed$objective <- objective
    }
  }
})

.assign_corr_types_families <- function(covStruct, # may be NULL on input
                                          corr_info, exp_ranef_types, exp_barlist) {
  
  corr_families <- vector('list',length(exp_ranef_types)) # full list and no char_rd so far
  special_ranefs <- .spaMM.data$keywords$special_ranefs # c("adjacency","Matern","Cauchy","AR1","corrMatrix", "IMRF", "corrFamily")
  # But note special case for IMRF... except that I commented it out...
  corr_types <- special_ranefs[match(exp_ranef_types, special_ranefs)] ## full length with (so far) NA's when no match
  is_cF_internally <- exp_ranef_types=="corrFamily" # unregistered ones; registered ones set to TRUE below 
  
  for (rd in seq_along(exp_ranef_types)) {
    ranef_type <- exp_ranef_types[rd]
    if (ranef_type=="(.|.)") {
      # do nothing
    # } else if (ranef_type=="IMRF") { # This one was an experimental test of the gridIMRF corrFamily constructor
                                       # It works when the IMRF is really gridIMRF, not MaternIMRFa; and it fails to distinguish between them. 
    #   term <- exp_barlist[[rd]]
    #   cov_term <- term[-2L]
    #   cov_term[[1L]] <- as.name("gridIMRF")
    #   covStruct[rd] <- list(corrFamily=cov_term) # see two comments below
    #   is_cF_internally[rd] <- TRUE 
    } else if ( ! is.na(corr_types[[rd]])) { # those not implemented as corrFamily: Matern, AR1... but not ARp
      corr_families[[rd]] <- do.call(corr_types[[rd]],list()) # corr_types' elements are function names!
    } else { # ranef_type is NA
      if (ranef_type %in% .spaMM.data$keywords$all_cF) { # registered corrFamily
        corr_types[[rd]] <- "corrFamily" # replaces NA
        term <- exp_barlist[[rd]]
        covStruct[rd] <- list(corrFamily=term[-2L]) # name is lost when covStruct was NULL... does it occur given condition on ranef_type? 
        is_cF_internally[rd] <- TRUE 
        corr_families[[rd]] <- do.call(corr_types[[rd]],list()) # well this calls corrFamily() which is a stub...
        # For these ones, corr_families[[rd]] will be set later by .preprocess_corrFamily() (though possibly delayed in fitmv case)
      } else stop("Unknown or unregsitered correlation model.")
    }
  }
  corr_info$is_cF_internally <- is_cF_internally
  corr_info$corr_types <- corr_types
  corr_info$corr_families <- corr_families
  
  covStruct
}


.def_off_fn <- function(X_off, ori_off) { # outer beta
  force(X_off)
  force(ori_off)
  if (is.null(ori_off)) ori_off <- 0
  function(beta) {
    off <- drop(ori_off + X_off %*% beta)
    attr(off,"beta") <- beta
    off
  }
}

.preprocess <- function(control.HLfit, ranFix=NULL, HLmethod, 
                       predictor, resid.model,
                       REMLformula, data, family,
                       BinomialDen, rand.families, etaFix, prior.weights,
                       objective=NULL, 
                       control.glm, # under user control! => impacts inits_by_xLM
                       adjMatrix=NULL, verbose=NULL, For,
                       init.HLfit=list(), # user's init.HLfit here
                       corrMatrix=NULL,covStruct=NULL,
                       distMatrix=NULL, 
                       control.dist=NULL,
                       init=NULL, # for .preprocess_augZXy() ... and outer-beta
                       X2X,
                       ADFun=NULL # default for private argument; any non-default value input by a user is copied in 'processed';
                                  # Possible ADFun values are then defined by the .wrap_MakeADFun() code.
                       ) {
  callargs <- match.call() 
  #
  ################ handling list of data #######################
  if (inherits(data,"list")) {
    locargs <- as.list(callargs)
    famfam <- family$family
    processed <- lapply(data, function(dd) {
      locargs$data <- dd
      if ( ! is.null(famfam) && famfam=="multi") locargs$family <- family$binfamily  
      eval(as.call(locargs)) ## call("preprocess",...) on each data set
    })
    return(processed) ## a list of environments
  }
  ############### initiate 'processed' envir ##############
  if (For_fitmv <- (For=="fitmv")) For <- "fitme" # i.e. treat as "fitme" except where identified by For_fitmv or callargs$For
  resid.model <- .reformat_resid_model(resid.model,check_old_syntax=control.HLfit$resid.family) ## calls .preprocess_formula() ## the list(...) is used even for poisson, binomial...
  
  ##### some family processing
  canonicalLink <- .is_link_canonical(family)
  if (not_spaMM_fam <- is.null(family$flags)) { # standard stats:: family or perhaps from 3rd package
    obs <- ( # According to the documentation:
      canonicalLink || 
      family$family=="binomial" ||## implemented through ad hoc bit of code outside the family object so perhaps that could have been done for Poisson)
        (family$family %in% c("Gamma","gaussian") && family$link=="log") #" Same remark...
    )
    family$flags <- list(exp=TRUE,obs=obs, LLgeneric=FALSE) # represent capacity for these exp or obs methods, resp., 
                                                            # whether through code in the family object or not
  }
  family$flags$canonicalLink <- canonicalLink
  family$flags$LMbool <- (canonicalLink && family$family=="gaussian")
  #
  # nrand now needed early...
  exp_barlist <- .process_bars(predictor,as_character=FALSE) ## but default expand =TRUE; also -> .parseBars() -> .process_IMRF_bar() parses RHS info
  exp_ranef_strings <- .process_bars(barlist=exp_barlist,expand=FALSE, as_character=TRUE) ## no need to expand again 
  nrand <- length(exp_ranef_strings)
  # ...for...
  obsAlgo_needed <- .need_obsAlgo(HLmethod, family, canonicalLink, nrand=nrand) 
  # : is .spaMM option set to FALSE, then only ad hoc obsInfo algos are available, in principle as follows:
  # binomial(), negbin() and poisson() (all links, and including zero-truncated variants), Gamma(log), and gaussian(log).
  # At least that was so in version 3.13.0. For poisson(), one may need to use *P*oisson(., LLgeneric=FALSE).
  # Maintaining this feature is not a priority.
  if (not_spaMM_fam && obsAlgo_needed && .spaMM.data$options$LLgeneric) family <- .statsfam2LLF(family)
  #####
  
  processed <- list2env(list(family=family, clik_fn=.get_clik_fn(family), For=For,
                             envir=list2env(list(), parent=environment(HLfit)),
                             port_env=new.env(parent=emptyenv()),
                             verbose=.reformat_verbose(verbose,For=For),
                             control.glm=do.call("glm.control", control.glm),
                             how=list(obsInfo=obsAlgo_needed),
                             intervalInfo = control.HLfit$intervalInfo, # used by .preprocess_augZXy and beyond
                             ADFun=ADFun
                             ))
  #
  if (is.null(main_terms_info <- attr(data,"updated_terms_info"))) { # standard case for primary fit
    if ( inherits(data,"data.frame")) {
      if (length(class(data))>1L) {
        data <- as.data.frame(data) # The pb addressed by this line is that the test data Orthodont is a 'groupedData' object and that 
        # .calc_AR1_sparse_Q_ranges() -> by() -> `[` creates a factor with levels "MO1" and values 1, 2... when the variable was char
        # where `[` on a generic data.frame (be it that data or the model frame created from it) returns a char variable.
        # So (although it inherits(.,"data.frame")) we use as.data.frame() to remove the groupedData class... (grumbles)
        if (length(class(data))>1L) {
          message(paste0("The data appear to inherit from more classes than 'data.frame'. Unexpected results may happen..."))
        } else paste0("The data inherited from more classes than 'data.frame' and have been converted to 'data.frame' only.")
      }
      #
      loccall <- callargs
      loccall[["formula"]] <- predictor
      loccall[["resid.formula"]] <- resid.model$formula
      # contains the data, potentially the prior.weights, plus irrelevant stuff
      # A more standard idiom might be 
      # m <- match(c(<some relevant names>), names(callargs), 0L)
      # loccall <- callargs[c(1L, m)]  
      loccall[[1L]] <- get(".GetValidData_info", asNamespace("spaMM"), inherits=FALSE)  ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
      validData_info <- eval(loccall,parent.frame()) # 
      
      data <- data[validData_info$rownames,,drop=FALSE] 
    } else {
      stop("'data' is not a data.frame.")
    }
    main_terms_info <- .get_terms_info(formula=predictor,data=data, famfam=family$family, weights=validData_info$weights) ## design matrix X, Y... 
  } else { # results from update response -> .update_main_terms_info() provide response-update model frame for mean response.
    # The data (although not any response variable) may still be needed for other purposes (=> the resid model) 
    # (what if the response of the main model were a predictor for the residual dispersion model?)
    attr(data,"updated_terms_info") <- NULL # immediately clean the 'data' attribute to avoid mix-ups.
    main_terms_info$Y <- .get_Y(full_frame=main_terms_info$mf, famfam=family$family)
    fixef_terms <- main_terms_info$fixef_terms
    if (fixef_terms[[length(fixef_terms)]]==0L) { ## check that the fixef are only an explicit '0' (odd that it compares to 0L, but it does)
      main_terms_info$X <- matrix(nrow=nrow(main_terms_info$mf),ncol=0L) ## model without fixed effects, not even an Intercept 
    } else {
      main_terms_info$X <- model.matrix(fixef_terms, main_terms_info$mf, contrasts.arg = NULL) ## always valid, but slower
    } 
  }
  #
  nobs <- NROW(main_terms_info$X) ## not using Y which may be NULL
  if (nobs==0L) stop("No line in the data have all the variables required to fit the model.")
  #
  ##### processed$data for post-fit, but also to evaluate glm_phi if not previously done... 
  if (.spaMM.data$options$store_data_as_mf) { ## *FALSE* : TRUE may not be far from working for univariate, but breaks predict(., newdata=<object>$data)
    # mv might be an issue (it's better to store raw data rather than several mf)
    # poly (incl. in LHS ranef term) might be a pb (mf stores monomials)
    processed$data <- structure(main_terms_info$mf, 
                                # rawvarnames=colnames(data), # not yet used
                                # fixef_terms=main_terms_info$fixef_terms, #in the main_terms_info
                                # fixef_levels=main_terms_info$fixef_levels, #in the main_terms_info
                                fixefvarnames=rownames(attr(main_terms_info$fixef_off_terms,"factors")), 
                                fixefpredvars=attr(main_terms_info$fixef_off_terms,"predvars"))
  } else processed$data <- structure(data,  # this will be used by preprocess_phi_model() so we cannot put main_terms_info$mf here. 
                                     # rawvarnames=main_terms_info(data), # not yet used
                                     # fixef_terms=main_terms_info$fixef_terms, #in the main_terms_info
                                     # fixef_levels=main_terms_info$fixef_levels, #in the main_terms_info
                                     fixefvarnames=rownames(attr(main_terms_info$fixef_off_terms,"factors")), 
                                     fixefpredvars=attr(main_terms_info$fixef_off_terms,"predvars"))
  # 
  processed$prior.weights <- .preprocess_pw(subs_p_weights=substitute(prior.weights), nobs, model_frame=main_terms_info$mf)
  ##### Storing main_terms_info EXCEPT $mf:
  ## Things needed in the update_response case, where processed$main_terms_info will serve as template for new main_terms_info,
  # include those needed for preprocessing the updated call, and those  needed in .update_main_terms_info():
  # Y for family$initialize() (binomial spec.)
  # fixef_off_terms for .calc_terms_heuristic_denseness() (we need the full attr(fixef_off_terms,"factors"))
  # fixef_terms was used above to rebuild the main_terms_info$X
  # fixef_levels for .determine_sparse_X() and .calc_terms_heuristic_denseness() ; [but also in postfit: .calc_newFrames_fixed() -> model.frame(., xlev=...)]
  processed$main_terms_info <- main_terms_info[c("Y","fixef_off_terms","fixef_terms","fixef_levels","specials_levels")] 
  #
  ### class(processed$main_terms_info) <- "HLframes" # Convenience: when we update its $mf, the class is kept
  #####
  processed$BinomialDen <- BinomialDen <- .calc_Binomial_Den(main_terms_info$Y, family, nobs)
  processed$y <- y <- main_terms_info$Y[,1L,drop=FALSE] # seems marginally faster as matrix; vector may 'fit' (cf comment in .merge_processed()) but formal extractor response.HLfit() expects a 1-col matrix. 
  processed$bin_all_or_none <- .check_y(family, y, BinomialDen)
  #### Various control parameters
  processed$spaMM_tol <-   .preprocess_spaMM_tol(processed$bin_all_or_none, control.HLfit)
  #
  processed$max.iter <- c(control.HLfit$max.iter,200L)[1]## control of outer loop 
  ## control of maxit.mean in inner loop i.e. in .solve_IRLS_as_...(), dependent on y (caseStudies/Jeanne/Leucadendron_hard.R was revelatory)
  processed$iter_mean_dispFix <- .calc_iter_mean_dispFix(control.HLfit, family, y)
  # ... when some disp param is inner estimated: (HLfit_body() tells cases apart) 
  processed$iter_mean_dispVar <- .calc_iter_mean_dispVar(control.HLfit, family, y)
  #
  if (nrand) {
    if (family$link %in% c("log","loglambda"))  {
      maxLambda <- log(sqrt(.spaMM.data$options$maxLambda))^2
    } else maxLambda <- .spaMM.data$options$maxLambda
    processed$maxLambda <- rep(maxLambda, nrand) # rep() useful only for unified syntax with mv
    #
    ## Initialize $corr_info (ASAP to assign_cov_matrices ASAP):
    processed$corr_info <- corr_info <- new.env() ## do not set parent=emptyenv() else with(corr_info,...) will not find trivial fns such as `[`
    exp_ranef_types <- attr(exp_ranef_strings,"type") ## expanded
    covStruct <- .assign_corr_types_families(covStruct=covStruct, corr_info=corr_info, exp_ranef_types=exp_ranef_types, 
                                             exp_barlist=exp_barlist) # provides 'corr_types' and 'is_cF_internally', soon necessary, and 'corr_families'
    processed$control_dist <- .preprocess_control.dist(control.dist, corr_info$corr_types)
    ## Assigns $corr_info$corrMatrices, $adjMatrices, $AMatrices using $corr_info$corr_type: (much) BEFORE determining sparse precision:
    .assign_cov_matrices__from_covStruct(corr_info, covStruct=covStruct, corrMatrix=corrMatrix, adjMatrix=adjMatrix,
                                         required= ! For_fitmv)
    exp_ranef_terms <- .process_bars(predictor[[length(predictor)]], 
                                     barlist=exp_barlist, 
                                     expand=TRUE, which. = "exp_ranef_terms")
    ## 
    # # Possible special treatment for IMRF term. We instead do that in .process_bars if passing the formula env recursively
    # # But if we were to avoid that. we need to reproduce all operations otherwise performed by .process_bars(): 
    # for (rd in seq_len(nrand)) {
    #   # RHS must still use the exp_barlist (terms with their type keyword as first element) instead of exp_ranef_terms
    #   if (attr(exp_ranef_terms,"type")[rd]=="IMRF") {
    #     full_term_with_attr <- .process_IMRF_bar(exp_barlist[[rd]], env=environment(predictor)) # 'predictor' has received the (evaluated) control.HLfitformula_env
    #     exp_ranef_terms[[rd]] <- .lhs_rhs_bars(list(full_term_with_attr))[[1]]
    #       }
    # }
    #####
    if (inherits(rand.families,"family")) rand.families <- list(rand.families) 
    if (nrand != 1L && length(rand.families)==1L) rand.families <- rep(rand.families,nrand) 
    names(rand.families) <- seq(nrand)
    # (fixme): this will make two copies of lcrandfamfam in processed
    rand.families <- .checkRandLinkS(rand.families)  
    processed$lcrandfamfam <- attr(rand.families,"lcrandfamfam") ## else remains NULL
    #
    #
    # (1) preprocess to make sure of the $levels_type:
    if ( ! For_fitmv) for (rd in which(corr_info$is_cF_internally)) { # (allows NA in $corr_types)
      # For the hard coded Matern(), AR1() etc. $corr_families[[it]] is already a list of functions $calc_moreargs, $canonize...
      # for corrFamily() by the next line it will be the corrFamily descriptor as an *environment* with $f, $tpar, $type, $template, $Af... $calc_moreargs, $canonize...
      corr_info$corr_families[[rd]] <- .preprocess_corrFamily(corrfamily=eval(covStruct[[rd]])) 
    } 
    # (2) Create Z matrices using corrFamilies' levels_type:
    Zlist <- .calc_Zlist(exp_ranef_terms=exp_ranef_terms, data=data, rmInt=0L, drop=TRUE, 
                         corr_info=corr_info,
                         lcrandfamfam=processed$lcrandfamfam) 
    # (3) initialize things in the corrFamilies using the Z matrices:
    if ( ! For_fitmv) for (rd in which(corr_info$is_cF_internally)) {
      .initialize_corrFamily(corr_info$corr_families[[rd]], Zmatrix=Zlist[[rd]])
    } 
    #
    #
    # rand.fam.fam info needed for Zlist and Zlist needed for next line => prevents merging it in .preprocess_rand_families()
    for (rd in seq_len(nrand)) rand.families[[rd]]$prior_lam_fac <- attr(Zlist[[rd]],"prior_lam_fac")
    processed$rand.families <- rand.families
    #####
    ## Assigns $corr_info$AMatrices for IMRFs:
    # __F I X M E__ there's an inefficiency: blob$uniqueGeo etc is computed twice : 
    # once for the new Z and once for the new A.
    # processed$geo_info does not appear to be used for IMRF nor corrFamily. 
    # If we were to use geo_envir <- .get_geo_info(processed, which_ranef=rd, which=c("uniqueGeo")) :
    # $geo_info is created by .init_assign_geoinfo() which is called much later than
    # both .preprocess_corrFamily() and .assign_AMatrices_corrFamily().
    # Should we create geo_info earlier than in .init_assign_geoinfo()? Currently this would not be compatible 
    # with the .preprocess_mv() workflow, for reasons detailed in the source of .init_assign_geoinfo().
    # Moreover, to take advantage of an early-defined geo_info, the formals of .assign_AMatrices_IMRF() and 
    # .calc_AMatrix_IMRF() would need to be modified. This seems a bit too messy currently. 
    
    .assign_AMatrices_IMRF(corr_info, Zlist, exp_barlist=exp_barlist, processed$data, control_dist=processed$control_dist) 
    ## Same for corrFamilies that have some:
    .assign_AMatrices_corrFamily(corr_info, Zlist, exp_barlist=exp_barlist, processed$data, control_dist=processed$control_dist)
    ZAlist <- .calc_ZAlist(Zlist=Zlist, AMatrices=corr_info$AMatrices)
    attr(ZAlist,"exp_ranef_strings") <- exp_ranef_strings ## expanded 
    attr(ZAlist,"exp_ranef_types") <- exp_ranef_types ## expanded

  } else processed$rand.families <- list() ## ## corrects the default [gaussian()] when nrand=0
  #
  HLmethod1 <- .preprocess_HLmethod(HLmethod[[1L]], family, processed$lcrandfamfam) ## not a member of the 'processed' object
  .preprocess_HL_REMLformula(HLmethod1, processed, BinomialDen, nobs, control.HLfit, y, REMLformula) # $HL, $REMLformula
  ####
  ####
  ## determine true # cols before determining sparse_precision (!= sparse_X)
  ## X.pv post-processing; later extract column for fixed coefficients, but we will need the offset for that
  if  (! is.null(X2X)) {
    if (is.null(colnames(X2X))) stop("'X2X' *must* have column names")
    
    merged_X <- .post_process_X(X.pv=main_terms_info$X %*% X2X,  # loss of attributes, "assign" notably
                                HL=processed$HL, rankinfo=control.HLfit$rankinfo,
                            sparse_X=.determine_sparse_X(main_terms_info) )
  } else {
    X.pv <- .post_process_X(X.pv=main_terms_info$X, HL=processed$HL, rankinfo=control.HLfit$rankinfo,
                            sparse_X=.determine_sparse_X(main_terms_info) )
  }
  # .post_process_X() calls .rankTrim() which adds attributes "namesOri" and "rankinfo"
  #
  ##################
  .check_init.HLfit(init.HLfit)
  if (nrand) {
    # Standardize init.HLfit in the one case where it may have corrPars (simplified version of .reformat_corrPars)
    processed$init_HLfit <-  .preprocess_init.HLfit(init.HLfit, corr_info) 
    ## : this, as residProcessed[["init_HLfit"]], can be used and modified by dhglm-specific code;
    ##  and  processed$init_HLfit is used by .deternine_augZXy() and ..calc_optim_args()
    if (For_fitmv) {
      sparse_precision <- FALSE
    } else {
      sparse_precision <- .wrap_determine_spprec(control.HLfit, ZAlist, processed, X.pv) # -> .-> uses init_HLfit
      ## F I X M E: ZAL diagnosis uses elements that may be modified afterwards and before .choose_QRmethod() reuses $G_diagnosis
      ## post-processing of corr_info depending on sparse_precision
      .process_corr_info_spprec(corr_info=corr_info, For=processed$For,sparse_precision=sparse_precision)
    }
    processed$is_spprec <- sparse_precision
    #
  } else {
    processed$is_spprec <- FALSE ## for .do_TRACE()
    processed$init_HLfit <- init.HLfit
  }
  if (processed$is_spprec) {
    processed$solve_IRLS_fn <- .solve_IRLS_as_spprec
  } else  processed$solve_IRLS_fn <- .solve_IRLS_as_ZX
  ###   
  # assigns $X.Re, $off, $objective but there is NO $X.pv
  X.pv <- .preprocess_X_XRe_off(main_terms_info, predictor, processed, X.pv, etaFix, data, objective, nobs) 
  #
  if (For=="is_separated") {
    return(ncol(X.pv) && is_separated(X.pv, as.numeric(y),verbose=FALSE)) 
  } else if (family$family  %in% c("binomial","betabin") && processed$bin_all_or_none) {
    abyss <- (ncol(X.pv) && is_separated(X.pv, as.numeric(y)))
  }
  ## Now X.pv has its final ncol := pforpv
  #
  if (.spaMM.data$options$X_scaling) { ## use scaled X.pv by default v.2.4.83
    ##       standard REML    ||      ML 
    if ( is.null(processed$REMLformula) || ncol(processed$X.Re)==0L) X.pv <- .scale(X.pv) # scales and adds "scaled:scale" attribute
  }

  if ( ! is.null(init$beta)) { # outer beta
    betanames <- names(init$beta)
    if (length(intersect(colnames(X.pv),betanames))!=length(init$beta)) stop("init$beta must have names matching those of the design matrix")
    # {
    #   warning("init$beta must have names matching those of the design matrix.\n init$beta names will be automatically set.",
    #           immediate. = TRUE)
    #   betanames <- names(init$beta) <- colnames(X.pv)
    # }
    # It's useless to rename init$beta here bc init is not preprocessed (___F I X M E___?): the init used later is the one in the call of the parent fn
    X_off <-.subcol_wAttr(X.pv, j=betanames, drop=FALSE)
    X.pv <- .subcol_wAttr(X.pv, j=setdiff(colnames(X.pv),betanames), drop=FALSE)
    processed$X_off_fn <- .def_off_fn(X_off, ori_off=processed$off)
  }
  models <- list(eta="",lambda="",phi="", rdispar="")
  if ( ! For_fitmv) {
    thread_nbr <- control.HLfit$NbThreads
    if (is.null(thread_nbr)) thread_nbr <- .spaMM.data$options$NbThreads # should be 1L by default
  }
  if ( nrand) {
    models[["eta"]] <- "etaHGLM" 
    ZAlist <- .init_assign_geoinfo(processed=processed, ZAlist=ZAlist, For=callargs$For,
                               exp_barlist=exp_barlist, nrand=nrand, distMatrix=distMatrix)
    vec_normIMRF <- .calc_vec_normIMRF(exp_ranef_terms=attr(ZAlist, "exp_ranef_terms"), corr_info=corr_info)   
    if (any(vec_normIMRF)) attr(ZAlist,"Zlist") <- Zlist
    processed$ZAlist <- ZAlist 
    # This, together with two commented lines in .is_identity(), is not clearly useful:
    #for (rd in seq_along(ZAlist)) attr(ZAlist[[rd]], "is_identity") <- .is_identity(ZAlist[[rd]], matrixcheck=TRUE)
    vec_n_u_h <- .eval_check_vec_n_u_h(ZAlist, nobs, processed)
    processed$psi_M <- rep(attr(rand.families,"unique.psi_M"),vec_n_u_h) # in position suitable to make it available for fit and to drop it from fit object
    processed$cum_n_u_h <- cumsum(c(0L, vec_n_u_h)) ## if two ranef,  with q=(3,3), this is 0,3,6 ;
    processed$reserve <- .preprocess_arglists(processed)
    processed$hyper_info <- .preprocess_hyper(processed=processed) # uses$ZAlist and $predictor
    if (For_fitmv) {
      processed$AUGI0_ZX <- list(X.pv=X.pv) # not most elegant, but avoids distinguishing subcases when no ranef too.
    } else {
      processed$QRmethod <- .choose_QRmethod(processed$ZAlist, corr_info=corr_info,
                                             is_spprec=processed$is_spprec, processed=processed, control.HLfit=control.HLfit)
      algebra <- .set_augX_methods(processed) # sets processed$corr_method In spprec case, both this and processed$spprec_method may be needed so must be distinct.
      .check_time_G_diagnosis(.provide_G_diagnosis, processed, algebra)  
      nrd <- processed$cum_n_u_h[nrand+1L]
      if (nrd==1L) warning("Found a single random effect with a *single level*. Check formula?", immediate.=TRUE)
          processed$AUGI0_ZX <- .init_AUGI0_ZX(X.pv, vec_normIMRF, processed$ZAlist, nrand, n_u_h=nrd, sparse_precision, 
                                           as_mat=.eval_as_mat_arg(processed))
      if (algebra=="decorr" && nrd>900L && thread_nbr==1L) message('Using paralellisation might be useful. See help("setNBThreads")')
      # The larger subsamples (nrd=1000 for the 3rd) in useR2021_spatial_timings.R may be used to test the effect IF method uses obsInfo (otherwise there is in particular no .tcrossprodCpp)
    }
  } else {
    models[["eta"]] <- "etaGLM"
    processed$AUGI0_ZX <- list(X.pv=X.pv)
  }
  #
  ## models: ## progressively add info to 'processed' and to 'models' before storing final 'models' in 'processed.'
  if (is.null(ranFix)) ranFix <- list()
  processed$lambda.Fix <- .reformat_lambda(.getPar(ranFix,"lambda"), nrand, namesTerms=attr(ZAlist,"namesTerms"), full_lambda = TRUE) # should always have nrand elements
  models <- .preprocess_lam_rC_models(processed, models, ranFix, nrand, ZAlist=processed$ZAlist)
  #processed$fixed <- ranFix # always a list() # would be used as proc1$fixed... ## : difficult to use since fitme_body(.,fixed=...) can be called.
  processed$phi.Fix <- .check_phi_Fix(phi.Fix=.getPar(ranFix,"phi"), family)
  control.HLfit$intervalInfo <- NULL # because control.HLfit will be passed to .preprocess_resid() where we don't want this info.
  models <- .preprocess_phi_model(processed, models, resid.model, 
                                  control.HLfit, # elements controlling may fit are by default passed to resid fit !!!!
                                  HLmethod=HLmethod, # that's the user-level input not the processed value => can keep "exp"
                                  data, processed$control.glm, family) # assigns things to 'processed' (incl $models, $residModel)
  #
  processed$models <- .setattr_G_LMMbool(models, processed=processed) # after modif by .preprocess_phi_model()
  if (obsAlgo_needed) {
    processed$etaxLM_fn <- .calc_etaLLMblob
  } else processed$etaxLM_fn <- .calc_etaGLMblob
  
  if (attr(processed[["models"]],"LMMbool") && ! For_fitmv) .check_identifiability_LMM(processed, nobs=nobs) # depending on mdoel booleans.
  processed$vecdisneeded <- .vecdisneeded(pforpv=ncol(X.pv), family, processed) # after setting $models
  #
  ##         ##
  #### Ultimate controls of algorithms
  processed$break_conv_logL <- c(control.HLfit$break_conv_logL, FALSE)[1]   ##whether to stop if logL (p_v) appears to have converged  
  # confint on a GLM assumes processed$LevenbergM exists
  processed$LevenbergM <- .preprocess_LevM(control.HLfit$LevenbergM, processed, nrand=nrand) # uses $models, $HL & optionally $bin_all_or_non & $cum_n_u_h
  #####
  .preprocess_augZXy(processed, init=init, ranFix=ranFix)  # may be modified later.
  # the problem of changing augZXy_cond later would be that .do_TRACE() -> may then be tracing the wrong function
  # So we minimize such later changes. But see handling of phi user controls in .calc_optim_args()
  #####
  if ( ! For_fitmv) {
    if (processed$augZXy_cond) {
      processed$HLfit_body_fn <- ".HLfit_body_augZXy"
      processed$HLfit_body_fn2 <- .spaMM.data$options$HLfit_body
      .do_TRACE(processed)
      processed$HLfit_body_fn <- get(processed$HLfit_body_fn, asNamespace("spaMM"), inherits=FALSE) 
      processed$HLfit_body_fn2 <- get(processed$HLfit_body_fn2, asNamespace("spaMM"), inherits=FALSE) 
    } else {
      processed$HLfit_body_fn <- .spaMM.data$options$HLfit_body
      .do_TRACE(processed)
      processed$HLfit_body_fn <- processed$HLfit_body_fn2 <- get(processed$HLfit_body_fn, asNamespace("spaMM"), inherits=FALSE) 
    }
    delayedAssign("HLCor_body", get("HLCor_body", asNamespace("spaMM"), inherits=FALSE), assign.env = processed) 
    delayedAssign("HLCor", get("HLCor", asNamespace("spaMM"), inherits=FALSE), assign.env = processed) 
    processed$HLfit <- get("HLfit", asNamespace("spaMM"), inherits=FALSE) 
    
    .check_nb_cores(nb_cores = thread_nbr)
    .setNbThreads(thr=thread_nbr)
  }
  processed$fitenv <- list2env(list(prevmsglength=0L))
  class(processed) <- c("arglist",class(processed))
  return(processed)
}

# MEMO : idiom for easy delayedAssign() calling functions with arguments not function in its eval.env or assign.env:
# (1) Define a function which as the assign.env in its own definition environmant, e.g. 
# .assign_get_inits_by_glm <- function(envir) {
#   envir$inits_by_glm <- NULL
#   envir$get_inits_by_glm <- function(processed,
#                                X.pv=processed$AUGI0_ZX$X.pv, family=processed$family, reset=FALSE) {
#     #if (is.null(envir$inits_by_glm) || eval(reset)) envir$inits_by_glm <<- .calc_inits_by_xLM(processed=processed)
#     
#     envir$inits_by_glm
#   }
# }
# (2) where this is to be called and as function of argument found there, un a delayedAssign calling this function, e.g.
# delayedAssign("inits_by_glm", {
#   processed$envir$get_inits_by_glm$get_inits_by_glm(processed, family=family,
#                       reset=quote(family$family %in% c("COMPoisson","negbin2")) )
#   }

.eval.update.call <- function(mc,...) { # not currently used
  mc <- as.list(mc)
  dotlist <- list(...)
  mc[names(dotlist)] <- dotlist ## a un moment j'ai mis cette ligne en commentaire, ce qui rend la fonction ineffective !
  eval(as.call(mc))  
}

.options.processed <- function(processed,...) {  # based on spaMM.options()
  if (nargs() == 0) return(processed)
  temp <- list(...)
  if (length(temp) == 1 && is.null(names(temp))) {
    arg <- temp[[1]]
    switch(mode(arg),
           list = temp <- arg,
           character = return(mget(arg, envir=processed, ifnotfound = list(NULL))),  ## return here for eg ... = "NUMAX"
           stop("invalid argument: ", sQuote(arg)))
  }
  if (length(temp) == 0) return(processed)
  argnames <- names(temp)
  if (is.null(argnames)) stop("options must be given by name")
  old <- mget(argnames, envir=processed, ifnotfound = list(NULL))
  names(old) <- argnames ## bc names are not valid for previously absent elements
  list2env(temp[argnames], envir=processed) # for (st in argnames) assign(st, temp[[st]], envir=processed)
  invisible(old)
}


