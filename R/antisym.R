antisym <- function() {
  
  oldZlevels <- NULL
  
  initialize <- function(Zmatrix, ...) {
    oldZlevels <<- colnames(Zmatrix)
  }
  
  Af <- function(newdata, 
                 term,
                 fit.=FALSE,
                 ...) {
    if (fit.) {
      Zlevels <- oldZlevels # provided by initialize()
    } else {
      rhs <- term[[2L]][[3L]]
      txt <- .DEPARSE(rhs) ## should be the rhs of (|) cleanly converted to a string by terms(formula,data) in .get_terms_info()
      RHS_info <- .as_factor(txt=txt,mf=newdata, type=parent.env(environment())$levels_type) # 'levels_type' as provided there by .preprocess_corrFamily
      # it is the same level tpe that .cacl_Zmatrix uses for corrFamily types
      Zlevels <- levels(RHS_info$factor)
    }
    Z2Ablob <- .dyad_Z2A(Zlevels, ord_pairs=TRUE) # orde_pairs different from ranGCA
    Z2A <- Z2Ablob$Z2A
    Afactor <- Z2Ablob$Afactor
    Alevels <- levels(Z2Ablob$Afactor) # order determined in .dyadZ2A. See comments there before trying to change.
    if (length(Alevels)==1L) { # that must be a single 'homozygote'
      Amatrix <- .trivial_incidMat * 0  # different from ranGCA
    } else {
      # __F I X M E__ A formula like interface would be nice, such as  
      # rhs <- term[[2L]][[3L]]
      # txt <- .DEPARSE(rhs) ## should be the rhs of (|) cleanly converted to a string by terms(formula,data) in .get_terms_info()
      # sparse.model.matrix(as.formula(paste("~","I(",txt,")")), newdata)
      # but the ids are not interpreted as factors ***even when*** id1 and id2 have been converted to factors in the newdata
      # and the .dyad_Z2A() processing is necessar to match the levels of the two factors => newdata cannot be directly used.
      Amatrix <- sparse.model.matrix(~ID1-1,Z2A) - sparse.model.matrix(~ID2-1,Z2A) # different from ranGCA :difference* of the two individual random effects
    }
    colnames(Amatrix) <- Alevels
    rownames(Amatrix) <- Zlevels # ... matches the levels of Z with reciprocal ordered pairs i:j and j:i...
    # Now using the "global" or the local Afactor depending on 'fit.':
    attr(Amatrix,"is_incid") <- FALSE
    
    Amatrix
  }
  
  Cf <- function(tpar) NULL 
  
  make_new_corr_lists <- function(newLv_env, which_mats, newZAlist, new_rd, ...) { 
    # Cf code for (.|.) ranefs: (OK even for fix_info)
    newLv_env$cov_newLv_oldv_list[new_rd] <- list(NULL) # .calc_sub_diagmat_cov_newLv_oldv() will be called as for (.|.)
    if (which_mats$nn[new_rd]) {
      newLv_env$cov_newLv_newLv_list[new_rd] <- list(NULL) 
    } else { 
      newLv_env$diag_cov_newLv_newLv_list[[new_rd]] <- rep(1,ncol(newZAlist[[new_rd]])) # allowing subsetting
    }
  }
  
  list(Af=Af, tpar=numeric(0L), Cf=Cf, initialize=initialize, 
       make_new_corr_lists=make_new_corr_lists,
       levels_type="data_order",
       sparsePrec=TRUE, possiblyDenseCorr=FALSE,
       tag="antisym")
}

if (FALSE) {
  register_cF("antisym")
  set.seed(123)
  nind <- 10
  x <- runif(nind^2) 
  id12 <- expand.grid(id1=seq(nind),id2=seq(nind))
  id1 <- id12$id1
  id2 <- id12$id2
  u <-  rnorm(nind,mean = 0, sd=0.5)
  
  ## additive individual effects:
  y <-  0.1 + 1*x + u[id1] -  u[id2] + rnorm(nind^2,sd=0.2)
  
  ## Same with non-additive individual effects:
  dist.u <- abs(u[id1] -  u[id2])
  z <- 0.1 + 1*x + dist.u + rnorm(nind^2,sd=0.2)
  
  dyaddf <- data.frame(x=x, y=y, z=z, id1=id1,id2=id2)
  # : note that this contains two rows per dyad, which avoids identifiability issues.
  
  # Enforce that interactions are between distinct individuals (not essential for the fit):
  dyaddf <- dyaddf[- seq.int(1L,nind^2,nind+1L),] 
  
  # Fits:
  #
  (addfit <- fitme(y ~x +antisym(1|id1+id2), data=dyaddf, control.HLfit=list(algebra="spcorr")))
  (addfit <- fitme(y ~x +antisym(1|id1-id2), data=dyaddf, control.HLfit=list(algebra="spprec"))) 
}
