if (FALSE) {  ## DOC:
  ## Il semble difficile de garder une info sur chaque element de lambda ou ranCoefs, particulierement parce que 
  #  les elements NULL de ranCoefs poseraient probleme pour relist(). Il faut plutôt utiliser les noms.
  essai <- list(a=1,b=NULL,c=2)
  relist(c(3,5),skeleton=essai) ## error mal documentee
  
  ## idiom for merging parameters
  varNames <- setdiff(names(init.HLfit),c("fixef","v_h")) ## maybe also corrPars ? empty list in init.HLfit...
  HLCor.args$fixed <- structure(.modify_list(fixed,init.HLfit[varNames]),
                                  type=.modify_list(.relist_rep("fix",fixed),
                                                    attr(init.HLfit,"type")[varNames]))  
  
  ## idiom for splitting parameters
  rPtype <- attr(ranPars,"type")
  if (is.null(rPtype) && length(ranPars)) { ## direct HLCor call
    HL.info$ranFix <- structure(ranPars,
                                type=.relist_rep("fix",ranPars))
    #HL.info$init.HLfit previously filled by dotlist[good_dotnames]
  } else { ## through corrHLfit or fitme call: ranPars inherits values from <'corrfitme'> (...,init.HLfit(...))
    u_rPtype <- unlist(rPtype)
    varNames <- names(which(u_rPtype=="var"))
    fix_outer_Names <- setdiff(names(u_rPtype),varNames) 
    ## init.HLfit must recover elements from ranPars! (bug detected by AIC( <SEM fit> ) in test-CAR.R where it must get rho...
    if (is.null(fix_outer_Names)) { ## can be NULL for corrMatrix case => not $ranFix    
      HL.info$init.HLfit <- .modify_list(HL.info$init.HLfit,ranPars) 
    } else { # builds a ranFix with types from rPtype (typically "fix" as outer" is set at the end of <'corrfitme'>, but see refit...)
      HL.info$ranFix <- structure(.remove_from_cP(ranPars,u_names=varNames), ## loses attributes
                                  type=.remove_from_cP(rPtype,u_rPtype, u_names=varNames) ) 
      HL.info$init.HLfit <- .modify_list(HL.info$init.HLfit,
                                         .remove_from_cP(ranPars,u_names=fix_outer_Names)) ## loses attributes
    }
  }

}

## derived from utils::modifyList ... works on named vectors! 
.modify_list <- function (x, val, obey_NULLs=TRUE) { # obey_NULLs = FALSE => NULL elements in val are ignored, as if inexistent
  if( is.null(x)) {
    if (is.null(val)) {
      return(NULL)
    } else return(val)
  } else if (is.null(val)) return(x) # but if val is a named list with explicit NULLs, those explicit NULLs will replace the corresponding LHS elements 
  #stopifnot(is.list(x), is.list(val)) # inefficient
  xnames <- names(x)
  vnames <- names(val)
  if ( ! obey_NULLs ) {
    is_null_vec <- sapply(val, is.null)
    vnames <- vnames[which( ! is_null_vec)]
  }
  vnames <- vnames[nzchar(vnames)]
  for (v in vnames) {
    if (v %in% xnames) {
      if ( is.list(x[[v]]) && is.list(val[[v]])) {
        x[[v]] <- .modify_list(x[[v]], val[[v]])
      } else if ( ! is.null(dim(val[[v]]))) { # if val[[v]] is a matrix names(val[[v]]) is not what we need here
        x[[v]] <- val[[v]]
      } else if ( ! is.null(nam <- names(val[[v]]))) { # handles val[[v]] being list, or vector 
        x[[v]][nam] <- val[[v]]
      } else x[[v]] <- val[[v]]
    } else x[[v]] <- val[[v]] 
  }
  x
}

.denullify <- function(x, modifier, vec_nobs=NULL) { # changes NULL's and not to NULLs
  if (is.null(vec_nobs)) {
    if (is.null(x)) x <- modifier
  } else if (.anyNULL(x) ) {
    for (mv_it in seq_along(modifier)) if ( is.null(x[[mv_it]])) x[mv_it] <- list(unlist(modifier[as.character(mv_it)])) # handling missing data properly
  }
  x
}

# getPar extract values from a list of lists, controlling that there is no redundancies between the lists => useful to *merge* lists 
# but in fact I do not seem to use this facility. .getPar() is applied to 'ranFix' (once to 'fixed')
# Argument 'which' can be any way of indexing a list
.getPar <- function(parlist,name,which=NULL, count=FALSE) { ## see .get_cP_stuff() to extract from first level or from an optional corrPars element !
  if ( ! is.null(which)) parlist <- parlist[[which]] 
  val <- parlist[[name]] 
  if (is.null(val)) { ## ie name not found a topmost level; scan sublists: NOT RECURSIVELY
    nmatch <- 0L
    val <- NULL
    for (it in seq_along(parlist)) { ## the sublists are typically lists that we wish to merge (see examples below)
      if (is.list(parlist[[it]]) && length(vv <- parlist[[it]][[name]])) {
        val <- vv
        nmatch <- nmatch+1L
      } 
    }
    if (count) return(nmatch)
    ## ELSE
    if (nmatch>1L) {
      stop(paste0("Found several instances of element '",name,"' in nested list: use 'which' to resolve this."))
    } 
    return(val)
  } else if (count) {return(1L)} else return(val) ## single first-level or [[which]] value
}

# .getPar(list("1"=list(a=1,b=2),"2"=list(a=3,c=4)),"b") ## 2
# .getPar(list("1"=list(a=1,b=2),"2"=list(a=3,c=4)),"c") ## 4
# .getPar(list("1"=list(a=1,b=2),"2"=list(a=3,c=4)),"a") ## error
# .getPar(list("1"=list(a=1,b=2),"2"=list(a=3,c=4)),"a",which=1) ## 1
# .getPar(list("1"=list(a=1,b=2),"2"=list(a=3,c=4)),"d") ## NULL

.get_cP_stuff <- function(typelist,name,which=NULL,count=FALSE) { 
  if (is.null(cP_types <- typelist$corrPars)) {
    .getPar(typelist,name,which=NULL,count=count)
  } else .getPar(cP_types,name,which=which,count=count)
}

.process_HLfit_corrPars <- function(init.HLfit, template) { ## the template should be provided by preprocess
  if (is.null(corrPars <- init.HLfit$corrPars)) {
    if (!is.null(rho <- init.HLfit$rho)) {
      return(relist(rho,template))
    } else return(NULL)
  } else return(corrPars)
}

.set_pars_stuff <- function(lhs_list, value, names_from) {
  u_lhs <- unlist(lhs_list) ## generates automatic names
  u_lhs[names(unlist(names_from))] <- value
  relist(u_lhs,lhs_list)
}

.rmNaN_fn <- function(x) if (is.list(x)) .rmNaN(x) else {if (is.character(x)) x[x!= "NaN"] else {x[!is.nan(x)]}}
## Recursively step down into list, removing all NaN elements from vectors and vectors of NaN from lists
.rmNaN <- function(x) {
  res <- vector("list",length(x))
  for(it in seq_along(x)) res[[it]] <- .rmNaN_fn(x[[it]]) 
  names(res) <- names(x) ## crucial (other attributes are lost !)
  len <- integer(length(res))
  for(it in seq_along(res)) len[it] <- length(res[[it]]) 
  res[len>0L]  
}


.remove_from_cP <- function(parlist, u_list=unlist(parlist), u_names) { ## not simply corrPars...
  if (length(u_names)) { ## if something to subtract
    u_list[u_names] <- rep(NaN,length(u_names))
    u_list <- relist(u_list,parlist)
    return(.rmNaN(u_list)) ## removes attributes
  } else return(parlist) ## DHGLM where all parameters are fixed.
}

remove_from_parlist <- function(parlist, removand=NULL, rm_names=names(unlist(removand))) {
  type <- attr(parlist,"type")
  if ( ! is.null(type)) type <- .remove_from_cP(type, u_names=rm_names)
  structure(.remove_from_cP(parlist,u_names=rm_names),
            type=type )
}

#extract a sublist from a (structured) list x according to a skeleton; used for mv code
.subPars <- function (x, skeleton) { 
  xnames <- names(x)
  sknames <- names(skeleton)
  sknames <- sknames[nzchar(sknames)]
  for (v in sknames) {
    if (v %in% xnames) {
      if (( is.list(x[[v]]) || inherits(x[[v]],"R6")) && is.list(skeleton[[v]])) {
        skeleton[[v]] <- .subPars(x[[v]], skeleton[[v]])
      } else if ( ! is.null(nam <- names(skeleton[[v]]))) { # ideally this test is always TRUE when it is reached 
        if (length(subnames <- intersect(nam, names(x[[v]])))) {
          skeleton[[v]] <- x[[v]][subnames] # sub-vector here
        } else skeleton[v] <- NULL # remove element from list
      } else skeleton[[v]] <- x[[v]]
    } else skeleton[[v]] <- x[[v]]
  }
  skeleton
}
