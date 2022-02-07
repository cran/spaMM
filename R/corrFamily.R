# build $template$parnames, $diagars, $poslist from user input $map, $fixed ; a corr CHM template may by added by other functions  

.corrfamily2Info <- function(corrfamily) {
  ## Standardize the structure of the object:
  if ( ! inherits(corrfamily,"list")) corrfamily <- list(map=corrfamily)
  corrfamily <- list2env(corrfamily) # so that .corrfamily2Matrix() can write into it...
  #
  if (is.null(corrfamily$"f")) { # private first attempt. 'Real' code is in the alternative
    map <- corrfamily$map
    if ( ! inherits(map,"dgCMatrix")) map <- as(map,"dgCMatrix") # for sp|de corr algo early conversion to dsC leads to 
    #  call of calc_Lunique_for_correl_algos -> mat_sqrt_dsCMatrix on unregularized matrix. Rething later (__F I X M E__)
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
    if (is.null(corrfamily$parnames)) corrfamily$parnames <- paste0("p", levels_)
    names(poslist) <- corrfamily$parnames
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
  } else {
    if (is.null(tpar <- corrfamily$tpar)) stop("'tpar' is required in the covStruct element for a corrFamily term.")
    if (is.null(corrfamily$parnames)) corrfamily$parnames <- names(corrfamily$tpar)
    if (is.null(corrfamily$parnames)) corrfamily$parnames <- paste0("p", seq_along(corrfamily$tpar))
    if (is.null(corrfamily$template)) corrfamily$template <- corrfamily$"f"(tpar)
    if (inherits(corrfamily$template,"CsparseMatrix")) corrfamily$sp_chk <- length(corrfamily$template@i)
    dim_template <- dim(corrfamily$template)
    if (dim_template[1L]!=dim_template[2L])  stop("corrMatrix's f(tpar) is not square")
  }
  corrfamily
}

.corrfamily2Matrix <- function(corrfamily, value) {
  if (is.null( f <- corrfamily$"f")) {
    poslist <- corrfamily$poslist
    template <- corrfamily$template
    for(param in names(value)) template@x[poslist[[param]]] <- value[[param]]
    template
  } else {
    corrm <- f(value)
    if ( ! is.null(corrfamily$sp_chk) && length(corrm@i)> corrfamily$sp_chk) {
      warning("f(tpar) was not a least sparse matrix in the corrFamily. Check 'tpar'.", immediate. = TRUE)
      corrfamily$sp_chk <- NULL # so that the check is performed only once
    }
    corrm
  }
}

.corrfamily_permute <- function(corrfamily, perm) { # redefines the template/the function once for all; should not be applied repeatedly
  if ( ! is.null(corrfamily$permuted)) stop("Trying to apply .corrfamily_permute() twice")
  perm <- pmatch(perm, rownames(corrfamily$template)) # converts names into indices so that permutations are faster, and that only rownames needed at user level. 
  if (is.null( f <- corrfamily$"f")) {
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
    corrfamily$"f" <- function(parvec) f(parvec)[perm,perm]
  }
  corrfamily$permuted <- TRUE
  corrfamily
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

