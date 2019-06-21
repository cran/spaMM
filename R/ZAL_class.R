setClassUnion("missingOrNULL", c("missing", "NULL"))

# ZAXlist is a representation of ZAL as a list 'LIST' of blocks for each ranef, 
# where each block is a M/matrix, or a ZA_QCHM object which is a list made of ZA and of a chol factor (L=solve(Q_CHM,system="Lt"))
# The aim of this class is to use solve(sparse chol_Q, ...) rather that to use the dense product ([pre-stored] solve(chol_Q)) %*% ... 
setClass("ZAXlist", slots = list( LIST = "list"))

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
              } else {
                res[[rd]] <- drop(.crossprod(zax_rd, y))
              }
            }
            unlist(res)
          })

for (.inh_y in c("matrix","Matrix")) {
  setMethod("crossprod", c(x = "ZAXlist", y= .inh_y),
            definition = function(x, y) {
              nrand <- length(x@LIST)
              res <- vector("list", nrand)
              for (rd in seq_len(nrand)) {
                if (inherits(zax_rd <- x@LIST[[rd]],"ZA_QCHM")) {
                  rhs <- .crossprod(zax_rd$ZA, y)
                  res[[rd]] <- (solve(zax_rd$Q_CHMfactor, b=rhs, system="L"))
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
