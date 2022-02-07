setClassUnion("missingOrNULL", c("missing", "NULL"))

# ZAXlist is a representation of ZAL as a list 'LIST' of blocks for each ranef, 
# where each block is a M/matrix, or a ZA_QCHM object which is a list made of ZA and of a chol factor (L=solve(Q_CHM,system="Lt"))
# The aim of this class is to use solve(sparse chol_Q, ...) rather that to use the dense product ([pre-stored] solve(chol_Q)) %*% ... 
setClass("ZAXlist", slots = list( LIST = "list"))

## This is ambiguous, bc result is not necessarily consistent with the generically expected result:
# as.matrix.ZAXlist <- function(x, as_matrix=FALSE, ...) .ad_hoc_cbind(x@LIST, as_matrix=as_matrix )

t.ZAXlist <- function(x) { # for the few uses of t(ZAL) that may occur on a ZAXlist (<-> spprec)
  # such as .calc_lev_from_hat() -> .m_Matrix_times_Dvec(t(ZAL), drop(dh0deta))
  # A test code is wfit <- HLfit(..., resid.model = ~ X3+I(X3^2) , data=wafers) with forced spprec (through test-confint-spprec.R)
  t(.ad_hoc_cbind(x@LIST,as_matrix=FALSE)) # wastes the benefits of ZALlist in spprec __F I X M E__
} 


setMethod("%*%", c(x = "ZAXlist", y= "numeric"),
          definition = function(x, y) {
            nrand <- length(x@LIST)
            res <- vector("list", nrand)
            sum_nc <- 0L
            for (rd in seq_len(nrand)) {
              if (inherits(zax_rd <- x@LIST[[rd]],"ZA_QCHM")) {
                nc <- ncol(zax_rd$Q_CHMfactor)
                lastcol <- sum_nc+nc
                res[[rd]] <- drop(zax_rd$ZA %*% solve(zax_rd$Q_CHMfactor, b=y[(sum_nc+1L):(lastcol)], system="Lt"))
              } else if (inherits(zax_rd,"ZA_Kron")) {
                nc <- ncol(zax_rd$Kronfacto)
                lastcol <- sum_nc+nc
                rhs <- zax_rd$Kronfacto %*% y[(sum_nc+1L):(lastcol)]
                res[[rd]] <- drop(zax_rd$ZA %*% rhs)
              } else {
                nc <- ncol(zax_rd)
                lastcol <- sum_nc+nc
                res[[rd]] <- drop(zax_rd %*% y[(sum_nc+1L):(lastcol)])
              }
              sum_nc <- lastcol
            }
            Reduce("+",res)
          })

for (.inh_y in c("matrix","Matrix")) {
  setMethod("%*%", c(x = "ZAXlist", y= .inh_y),
            definition = function(x, y) {
              nrand <- length(x@LIST)
              res <- vector("list", nrand)
              sum_nc <- 0L
              for (rd in seq_len(nrand)) {
                if (inherits(zax_rd <- x@LIST[[rd]],"ZA_QCHM")) {
                  nc <- ncol(zax_rd$Q_CHMfactor)
                  lastcol <- sum_nc+nc
                  res[[rd]] <- (zax_rd$ZA %*% solve(zax_rd$Q_CHMfactor, b=y[(sum_nc+1L):(lastcol),,drop=FALSE ], system="Lt"))
                } else if (inherits(zax_rd,"ZA_Kron")) {
                  nc <- ncol(zax_rd$Kronfacto)
                  lastcol <- sum_nc+nc
                  stop("code missing in %*%,Kronfacto,<m|M>atrix-method") # 
                  rhs <- zax_rd$Kronfacto %*% y[(sum_nc+1L):(lastcol)]
                  res[[rd]] <- drop(zax_rd$ZA %*% rhs)
                } else {
                  nc <- ncol(zax_rd)
                  lastcol <- sum_nc+nc
                  res[[rd]] <- (zax_rd %*% y[(sum_nc+1L):(lastcol), , drop=FALSE ])
                }
                sum_nc <- lastcol
              }
              Reduce("+",res)
            })
}

setMethod("%*%", c(x = "numeric", y= "ZAXlist"),
          definition = function(x, y) {
            nrand <- length(y@LIST)
            res <- vector("list", nrand)
            for (rd in seq_len(nrand)) {
              if (inherits(zax_rd <- y@LIST[[rd]],"ZA_QCHM")) {
                lhs <- drop(x %*% zax_rd$ZA) # lhs in x . ZA . solve(Q_,"Lt") is rhs in next line:
                res[[rd]] <- drop(solve(zax_rd$Q_CHMfactor, b=lhs, system="L"))
              } else if (inherits(zax_rd,"ZA_Kron")) {
                stop("code missing in %*%,numeric,ZAXlist-method")
              } else {
                res[[rd]] <- drop(x %*% zax_rd)
              }
            }
            unlist(res) # cbind for a matrix x
          })

setMethod("crossprod", c(x = "ZAXlist", y= "numeric"),
          definition = function(x, y) {
            nrand <- length(x@LIST)
            res <- vector("list", nrand)
            for (rd in seq_len(nrand)) {
              if (inherits(zax_rd <- x@LIST[[rd]],"ZA_QCHM")) {
                rhs <- .crossprod(zax_rd$ZA, y)
                res[[rd]] <- drop(solve(zax_rd$Q_CHMfactor, b=rhs, system="L"))
              } else if (inherits(zax_rd,"ZA_Kron")) {
                stop("code missing in crossprod,ZAXlist,numeric-method")
              } else {
                res[[rd]] <- drop(.crossprod(zax_rd, y))
              }
            }
            unlist(res)
          })

if (FALSE) {
  ## This method was not needed until Kronfacto was introduced. But then, using 
  ## ZAL <- get_ZALMatrix(object, force_bind = ! ("AUGI0_ZX_sparsePrecision" %in% object$MME_method) ) appeared better.
  setMethod("crossprod", c(x = "ZAXlist", y= "missingOrNULL"), 
            # Occurs in get_predVar(fit1) -> ... -> .calcD2hDv2 ->  crossprodZAL <- .crossprod(ZAL)
            # There are cross-block products, so... 
            definition = function(x, y) {
              ZAL <- .ad_hoc_cbind(x@LIST, as_matrix=FALSE )
              .crossprod(ZAL)
            })
}


for (.inh_y in c("matrix","Matrix")) {
  setMethod("crossprod", c(x = "ZAXlist", y= .inh_y),
            definition = function(x, y) {
              nrand <- length(x@LIST)
              res <- vector("list", nrand)
              for (rd in seq_len(nrand)) {
                if (inherits(zax_rd <- x@LIST[[rd]],"ZA_QCHM")) {
                  rhs <- .crossprod(zax_rd$ZA, y)
                  res[[rd]] <- (solve(zax_rd$Q_CHMfactor, b=rhs, system="L"))
                } else if (inherits(zax_rd,"ZA_Kron")) {
                  rhs <- .crossprod(zax_rd$ZA, y)
                  res[[rd]] <- as.matrix(crossprod(zax_rd$Kronfacto, rhs))
                } else {
                  res[[rd]] <- (.crossprod(zax_rd, y))
                }
              }
              do.call(rbind,res)
            })
}

setMethod("tcrossprod", c(x = "ZAXlist", y= "missingOrNULL"),
          definition = function(x, y=NULL) {
            nrand <- length(x@LIST)
            tcrossfac <- x@LIST
            for (rd in seq_len(nrand)) {
              if (inherits(zax_rd <- tcrossfac[[rd]],"ZA_QCHM")) {
                tcrossfac[[rd]] <- t(solve(zax_rd$Q_CHMfactor, b=t(zax_rd$ZA), system="L")) # but if $ZA were NULL, solve(,"A)
              } else if (inherits(zax_rd,"ZA_Kron")) {
                stop("code missing in tcrossprod,ZAXlist,missingOrNULL-method")
              }
            } 
            tcrossfac <- do.call(cbind,tcrossfac) # ad_hoc_cbind ??
            .tcrossprod(tcrossfac) # y is NULL
          })

.ZAXlist_times_Dvec <- function(X, Dvec) { # Derived from %*%
  nrand <- length(X@LIST)
  res <- vector("list", nrand)
  sum_nc <- 0L
  for (rd in seq_len(nrand)) {
    if (inherits(zax_rd <- X@LIST[[rd]],"ZA_QCHM")) {
      nc <- ncol(zax_rd$Q_CHMfactor)
      lastcol <- sum_nc+nc
      res[[rd]] <- (zax_rd$ZA %*% solve(zax_rd$Q_CHMfactor, b=Diagonal(x=Dvec[(sum_nc+1L):(lastcol)]), system="Lt"))
    } else {
      nc <- ncol(zax_rd)
      lastcol <- sum_nc+nc
      res[[rd]] <- (zax_rd %*% Diagonal(x=Dvec[(sum_nc+1L):(lastcol)]))
    }
    sum_nc <- lastcol
  }
  do.call(cbind,res)
}

setClass("Kronfacto", slots = list( BLOB = "environment"))

dim.Kronfacto <- function(x) return(x@BLOB$DIM)

# as.matrix.Kronfacto <- function(x, ...) return(as.matrix(x@BLOB$long))

# t.Kronfacto <- function(x) .def_Kronfacto(lhs=t(x@BLOB$lhs),rhs=t(x@BLOB$rhs))

.def_Kronfacto <- function(lhs, rhs) {
  BLOB <- list2env(list(lhs=lhs, rhs=rhs), parent=environment(.AUGI0_ZX_sparsePrecision))
  delayedAssign("long", {.makelong_kronprod(Lcompact=lhs,kron_Y=rhs)}, eval.env = BLOB, assign.env = BLOB)
  delayedAssign("DIM", {
    diml <- dim(lhs)
    dimr <- dim(rhs)
    if (is.null(diml)) diml <- c(1L,length(lhs))
    if (is.null(dimr)) dimr <- c(1L,length(rhs))
    diml*dimr
  }, eval.env = BLOB, assign.env = BLOB)
  #
  new("Kronfacto", BLOB=BLOB)
}

setMethod("%*%", c(x = "Kronfacto", y= "numeric"),
          definition = function(x, y) {
            BLOB <- x@BLOB
            if (.is_evaluated("long",BLOB)) {
              BLOB$long %*% y
            } else {
              yy <- y
              dim(yy) <- c(ncol(BLOB$rhs),ncol(BLOB$lhs)) # nc_r, nc_l
              tmp <- BLOB$rhs %*% yy # nr_r * nc_l
              tmp <- tcrossprod(tmp,BLOB$lhs) # nr_r * nr_l
              dim(tmp) <- c(nrow(x),1L)
              tmp
            }
          })

setMethod("crossprod", c(x = "Kronfacto", y= "numeric"),
          definition = function(x, y) {
            BLOB <- x@BLOB
            if (.is_evaluated("long",BLOB)) {
              .crossprod(BLOB$long, y)
            } else {
              yy <- y
              dim(yy) <- c(nrow(BLOB$rhs),nrow(BLOB$lhs)) # nr_r, nr_l
              tmp <- .crossprod(BLOB$rhs, yy, chk_sparse2mat=FALSE) # nc_r * nr_l
              tmp <- tmp %*% BLOB$lhs # nc_r * nc_l
              dim(tmp) <- c(ncol(x),1L) 
              dimnames(tmp) <- list(NULL,NULL)
              tmp
            }
          })

if (FALSE) { # rethink when the need appears
  setMethod("crossprod", c(x = "Kronfacto", y= "missingOrNULL"),
            definition = function(x, y) {
              BLOB <- x@BLOB
              if (.is_evaluated("long",BLOB)) {
                .crossprod(BLOB$long)
              } else {
                .crossprod(BLOB$long)
                # lhs <- crossprod(x, BLOB$lhs) 
                # .def_Kronfacto(lhs = lhs, rhs=BLOB$rhs)
              }
            })
}


for (.inh_y in c("matrix","Matrix")) {
  setMethod("crossprod", c(x = "Kronfacto", y= .inh_y),
            definition = function(x, y) {
              BLOB <- x@BLOB
              if (.is_evaluated("long",BLOB)) {
                .crossprod(BLOB$long, y)
              } else if (ncol(y)==1L) {
                crossprod(x,y[,1L]) 
              } else {
                tmp <- matrix(NA,nrow=nrow(x),ncol=ncol(y))
                for (jt in seq_len(ncol(y))) tmp[,jt] <- crossprod(x,y[,jt])
                dimnames(tmp) <- list(NULL,NULL)
                tmp
              }
            })
}



if (FALSE) {
  a <- matrix(1:4,ncol=2)
  b <- matrix((1:4)+1,ncol=2)
  kronecker(a,b) %*% c(7,11,13,17) # 340 448 492 648
  Kf <- spaMM:::.def_Kronfacto(a,b)
  Kf %*%  c(7,11,13,17)
}
