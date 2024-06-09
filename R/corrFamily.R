

.stub <- function(errmess) {
  force(errmess)
  fn <- function(...) {
    stop(errmess)
  }
  class(fn) <- c("stub", class(fn))
  fn
}

.define_corrFamily_from_map <- function(corrfamily) {
  # private first attempt. 
  # build $template$parnames, $diagpars, $poslist from user input $map, $fixed ; a corr CHM template may by added by other functions  
  map <- corrfamily$map
  if ( ! inherits(map,"dgCMatrix")) map <- as(map,"dgCMatrix") # for sp|de corr algo early conversion to dsC leads to 
  #  call of calc_Lunique_for_correl_algos -> mat_sqrt_dsCMatrix on unregularized matrix. Rethink later (_F I X M E__)
  corrfamily$map <- map
  fixed <- corrfamily$fixed
  if ( ! is.null(fixed)) {
    if (identical(fixed,"unit diag")) fixed <- Matrix::.symDiagonal(ncol(map))
    if (! inherits(fixed,"dgCMatrix")) fixed <- as(fixed,"dgCMatrix")
    corrfamily$fixed <- fixed # maybe no necessary
    augmap <- fixed
    augmap@x <- rep(NA_real_,length(augmap@x))
    augmap <- augmap+map # NA's is fixed positions
  } else augmap <- map
  ## 
  levels_ <- unique(map@x) # those from the map 
  poslist <- vector("list", length(levels_))
  for (lev in seq_along(levels_)) poslist[[lev]] <- which(augmap@x==levels_[lev])
  #
  parnames <- corrfamily$parnames
  if (is.null(parnames)) parnames <- paste0("p", levels_)
  names(poslist) <- corrfamily$parnames <- parnames # parnames needed in the envir of corrFamily$canonize, at least (or restrict that function to the $Cf case)
  corrfamily$poslist <- poslist
  #
  diaglevels <- unique(diag(map)) 
  corrfamily$diagpars <- corrfamily$parnames[which(levels_ %in% diaglevels)] # for .calc_inits_corrFamily()
  #
  template <- map
  template@x <- rep(NA_real_,length(template@x))
  if ( ! is.null(fixed)) template <- template+fixed # NA's in map position
  dim_template <- dim(template)
  if (dim_template[1L]!=dim_template[2L])  stop("corrMatrix's template is not square")
  corrfamily$template <- template
  corrfamily$map <- corrfamily$fixed <- NULL
}

.preprocess_corrFamily <- function(corrfamily) { 
  ## Standardize the structure of the object:
  if ( ! inherits(corrfamily,"list")) corrfamily <- list(map=corrfamily)
  corrfamily <- list2env(corrfamily) # so that .corrfamily2Matrix() can write into it...
  
  if (is.null(corrfamily$"Cf")) {
    corrfamily <- .define_corrFamily_from_map(corrfamily) # private first attempt. 'Real' code is in the alternative
  } else {
    if (is.null(corrfamily$tpar)) stop("'tpar' is required in the corrFamily descriptor (i.e., the 'covStruct' element for the random effect).")
    # so as to distinguish between forgotten tpar and no parameter case (=>numeric(0))
    #
    parnames <- corrfamily$parnames
    if (length(corrfamily$tpar)) { # allow zero-length tpar
      if (is.null(parnames)) parnames <- names(corrfamily$tpar)
      if (is.null(parnames)) {
        parnames <- paste0("p", seq_along(corrfamily$tpar)) # that means that predefined names must be assigned to tpar by the constructor in all cases
      }
    }
    names(corrfamily$tpar) <- parnames
    
    fixed <- corrfamily$fixed
    np <- length(parnames)

    if ( ! is.null(fixed)) {
      fixnames <- names(fixed)
      if (is.null(fixnames)) {
        if (length(fixed) != np) {
          stop("Ambiguous 'fixed' vector for corrFamily parameters: provide parameter names or a complete vector.")
        } else {
          names(fixed) <- parnames
          corrfamily$fixed <- fixed
        }
      } else if (length(setdiff(fixnames, parnames))) {
        stop("Invalid 'fixed' vector for corrFamily parameters: their names do not match the 'tpar' names.")
      } 
      
      tpar <- corrfamily$tpar
      which_varpars <- which(is.na(match(parnames, names(fixed))))
      
      fulltpar <- tpar
      fulltpar[names(fixed)] <- fixed
      
      oldCf <- corrfamily$"Cf"
      
      corrfamily$"Cf" <- function(parvec) {
        fullparvec <- fulltpar
        fullparvec[which_varpars] <- parvec
        oldCf(fullparvec)
      }
      
      corrfamily$tpar <- tpar[which_varpars] # redefined to match the reducd init, lower, upper...
      parnames <- parnames[which_varpars] # idem
      
      # now the environment of Cf is redefined, it contains a local copy of tpar, which_varpars, fulltpar, and oldCf
      # the environment of oldCf is the envir of the original $f so oldCf will still be aware of whatever was in the local envir of the corrFamily constructor
    }  
    #
    corrfamily$parnames <- parnames 
    
  }
  
  # Provide default functions with uniform API
  if (is.null(corrfamily$canonize)) corrfamily$canonize <- function(corrPars_rd, checkComplete, ...) {
    list(corrPars_rd=corrPars_rd[parnames]) # ensure ordering of merged fixed and variable parameters; 'parnames' being in the definition environment 
  } ## no transfo defined yet
  
  # This will be called through .calc_optim_args -> .calc_inits -> corr_info$corr_families[[rd]]$calc_inits
  if (is.null(corrfamily$calc_inits)) corrfamily$calc_inits <- function(inits, char_rd, 
                                                                        # optim.scale, # currently ignored (not passed)
                                                                        user.lower,
                                                                        user.upper,
                                                                        ...) { # possibly add moreargs_rd
    # This calls a default generic function that should provide valid inits
    inits <- .calc_inits_corrFamily(corrfamily=corrfamily, init=inits[["init"]], char_rd=char_rd, optim.scale="", 
                                    init.optim=inits$init.optim, init.HLfit=inits$init.HLfit, ranFix=inits$ranFix, 
                                    user.lower=user.lower,user.upper=user.upper)
    return(inits)
  }
  if (is.null(corrfamily$calc_moreargs)) corrfamily$calc_moreargs= function(...) {return(NULL)} ## NULL OK since it returns in [[char_rd]], not [[rd]]
  if (is.null(corrfamily$levels_type)) corrfamily$levels_type <- "data_order" # this is $uGeo_levels_type
  if ( is.function(corrfamily$Af)) assign("levels_type", corrfamily$levels_type, environment(corrfamily$Af)) 
  if (is.null(corrfamily$need_Cnn)) corrfamily$need_Cnn <- TRUE
  if (is.null(corrfamily$possiblyDenseCorr)) corrfamily$possiblyDenseCorr <- TRUE
  if (is.null(corrfamily$sparsePrec)) corrfamily$sparsePrec <- FALSE
  
  if (is.null(corrfamily$calc_corr_from_dist)) corrfamily$calc_corr_from_dist <- 
    .stub("<corrFamily>$calc_corr_from_dist() may be be needed for some correlation models.")
  #

  corrfamily
}

.initialize_corrFamily <- function(corrfamily, Zmatrix) {
  if ( ! is.null(corrfamily$initialize)) corrfamily$initialize(Zmatrix=Zmatrix)
  
  if (is.null(corrfamily$"Cf")) { # not API
    rownames(corrfamily$template) <- colnames(Zmatrix) 
  } else {
    if (is.null(corrfamily$template)) corrfamily$template <- corrfamily$"Cf"(corrfamily$tpar)
    if (!is.null(corrfamily$template)) { # allow NULL result
      dim_template <- dim(corrfamily$template)
      if (dim_template[1L]!=dim_template[2L])  stop("corrMatrix's Cf(tpar) is not square")
    }
  }
  
  if (inherits(corrfamily$template,"CsparseMatrix")) corrfamily$sp_chk <- length(corrfamily$template@i)
  # corrfamily
}

.corrfamily2Matrix <- function(corrfamily, parvec, AUGI0_ZX_envir, rd) {
  if (is.null( Cf <- corrfamily$"Cf")) {
    AUGI0_ZX_envir$updateable[rd] <- (all(parvec!=0)) 
    poslist <- corrfamily$poslist
    template <- corrfamily$template
    for(param in names(parvec)) template@x[poslist[[param]]] <- parvec[[param]]
    template
  } else {
    corrm <- Cf(parvec)
    if ( ! is.null(corrfamily$sp_chk)) { # Was automatically added by .preprocess_corrFamily() for sparse $Cf(tpar).
      if (length(corrm@i)> corrfamily$sp_chk) { # imperfect test, so some user programming bugs might not be caught.  
        warning("Cf(tpar) was not a least sparse matrix in the corrFamily. Check 'tpar'.", immediate. = TRUE)
        corrfamily$sp_chk <- NULL # so that the check no longer performed and is permanently set to FALSE
        AUGI0_ZX_envir$updateable[rd] <- FALSE 
      } else AUGI0_ZX_envir$updateable[rd] <- TRUE 
    } # else do nothing ; Thus:
    # for dense matrices: updateable[rd] was initialized to FALSE for CORR algo, and .init_AUGI0_ZX_envir_spprec_info() set it to TRUE;
    # for sparse matrices: once the check has set it to FALSE, it remains permanently so.
    if (( ! is.null(corrfamily$type)) && corrfamily$type=="precision") {
      # the reason why corrfamily's $f() still returns a matrix rather than list(matrix=<precision mat>) in that case 
      # is the permutation code in .corrfamily_permute(), that changes the definition of f() by only adding the permutation operation.
      # We could redefine .corrfamily_permute(), but the optional presence of the 'type' is 
      # also useful for the test by any(.unlist(lapply(corr_info$corr_families,`[[`, "type"))=="precision") in .determine_spprec().
      structure(list(matrix=corrm),class=c("list","precision")) # so that it can be handled in the same way as a "corrMatrix" precision matrix 
    } else corrm
  }
}

.corrfamily_permute <- function(corrfamily, perm) { # redefines the template/the function once for all; should not be applied repeatedly
  if ( ! is.null(corrfamily$permuted)) stop("Trying to apply .corrfamily_permute() twice")
  perm <- pmatch(perm, rownames(corrfamily$template)) # converts names into indices so that permutations are faster, and that only rownames needed at user level. 
  if (is.null( Cf <- corrfamily$"Cf")) {
    # reconstuct the augmap with NA's is fixed positions to deduce its poslist:
    augmap <- corrfamily$template
    augmap@x <- rep(NA_real_,length(augmap@x))
    levels_ <- seq_along(corrfamily$poslist)
    for (lev in seq_along(levels_)) augmap@x[corrfamily$poslist[[lev]]] <- lev
    augmap <- augmap[perm,perm]
    
    # permuted poslist:
    poslist <- vector("list", length(levels_))
    for (lev in seq_along(levels_)) poslist[[lev]] <- which(augmap@x==levels_[lev])
    names(poslist) <- corrfamily$parnames
    corrfamily$poslist <- poslist
    corrfamily$template <- corrfamily$template[perm,perm] # permuted template with NA's in poslist positions 
  } else {
    corrfamily$"Cf" <- function(parvec) Cf(parvec)[perm,perm]
  }
  corrfamily$permuted <- TRUE
  corrfamily # return value not essential since it's an envir...
}

.subset_corrFamily <- function(corrfamily, ZAnames, corrnames) {
  template <- corrfamily$template
  uZAnames <- unique(ZAnames)
  if ( length(setdiff(corrnames,uZAnames)) || any(corrnames!=uZAnames) ) { # subsetting || reordering necessary
    corrfamily <- .corrfamily_permute(corrfamily, perm=uZAnames)  
  } ## else orders already match
  # and this should be correct in the case of repeated ZAnames (composite ranef) as much as it is with the unpermuted corrmatrix
  return(corrfamily)
}

