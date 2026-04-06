# currently returns a list with single element 'norm'
.qresiduals <- function(object, families=object$families,
                        #
                        fam=family(object),
                        y=response(object),
                        pw=object$prior.weights,
                        pw_resvar= residVar(object, which = "phi"), # includes pw, according to the doc
                        dev_res=residuals(object, type="deviance"),
                        mu_U= if (identical(fam$zero_truncated,TRUE)) {
                          fam$linkinv(predict(object,type="link")) # cf ?Poisson...
                        } else fitted(object),
                        BinomialDen=object$BinomialDen
) {
  if ( ! is.null(families)) { # mv case, list of families
    qres <- vector("list", length(families))
    for (mv_it in seq_along(families)) {
      cum_nobs <- attr(families,"cum_nobs")
      resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
      eta <- predict(object,type="link")
      fam <- families[[mv_it]]
      qres[[mv_it]] <- .qresiduals(fam=fam, families=NULL,
                                   y=y[resp_range], pw=pw[[mv_it]],
                                   pw_resvar=pw_resvar[resp_range],
                                   dev_res=dev_res[resp_range],
                                   mu_U=fam$linkinv(eta[resp_range]),
                                   BinomialDen=BinomialDen[resp_range])
    }
    # prospective generic code allowing additional future elements 
    qres <- do.call(rbind, qres)
    QRES <- list()
    for(st in colnames(qres)) {
      QRES[[st]] <- .unlist(qres[,st])
    }
    QRES
  } else {
    famfam <- fam$family
    qres <- switch(
      famfam,
      "gaussian" = {
        list(norm=dev_res/sqrt(pw_resvar))
      },
      "beta_resp" = {
        prec <- .get_family_par(family=fam)
        precW <- prec*pw
        u.log <- stats::pbeta(y, shape1=mu_U*precW,shape2=(1-mu_U)*precW, log.p = TRUE)
        list(norm=qnorm(u.log, log.p=TRUE)) # qres
      },
      "Gamma" = {
        u.log <- stats::pgamma(y/(mu_U*pw_resvar), 1/pw_resvar, log.p = TRUE)
        list(norm=qnorm(u.log, log.p=TRUE)) # qres
      },
      ## all other cases are presumably count families, for which randomization is used:
      {
        switch(famfam,
               "binomial" = {
                 a <- stats::pbinom(y - 1L, BinomialDen, mu_U)
                 b <- stats::pbinom(y, BinomialDen, mu_U)
               },
               "betabin" = {
                 prec <- .get_family_par(family=fam)
                 precW <- prec*pw
                 ## not the most efficient implem... being reasonably lazy here.
                 dbetabin <- function(y, n, prec, mu) { 
                   lchoose(n,y) + 
                     lbeta(y+mu*prec, n-y+(1-mu)*prec) - 
                     lbeta(mu*prec, (1-mu)*prec)
                 } # as in the family's $logl().
                 a <- rep(0, length(y)) 
                 for (ii in 0L:(max(y)-1L)) {
                   pos <- (ii < y)
                   dbb <- dbetabin(ii, n=BinomialDen[pos], prec=precW[pos], mu=mu_U[pos])
                   a[pos] <- a[pos] + exp(dbb)
                 }
                 dbb <- dbetabin(y, n=BinomialDen, prec=precW, mu=mu_U)
                 b <- a + exp(dbb)
               },
               "poisson" = {
                 a <- stats::ppois(y - 1L, mu_U)
                 b <- stats::ppois(y, mu_U)
               },
               "COMpoisson" = {
                 family_env <- environment(fam$aic)
                 COMP_nu <- family_env$nu # no pw implemented yet for this family
                 lambda <- family_env$mu2lambda(mu_U)
                 nobs <- length(y)
                 a <- b <- numeric(nobs) 
                 for (i in seq_len(nobs)) {
                   aa  <- 0
                   mu_i <- mu_U[i]
                   lambda_i <- lambda[i]
                   maxn_i <- .COMP_maxn(lambda_i,COMP_nu)
                   for (ii in 0L:(y[i]-1L)) { # a sort of pCOMPoisson(), up to y-1
                     dbb <- .dCOMP(ii, mu=mu_i, # family_env,
                                   nu=COMP_nu,
                                   lambda=lambda_i,
                                   log = FALSE, maxn=maxn_i)
                     aa <- aa + dbb
                   }
                   a[i] <- aa 
                   dbb <- .dCOMP(y, mu=mu_i,# family_env,
                                 nu=COMP_nu,
                                 lambda=lambda_i,
                                 log = FALSE, maxn=maxn_i)
                   b[i] <- a[i] + dbb
                 }
               },
               "negbin1" = { 
                 NB_shape <- .get_family_par(family=fam)
                 a <- stats::pnbinom(y - 1L, size=NB_shape*mu_U, mu=mu_U)
                 b <- stats::pnbinom(y, size=NB_shape*mu_U, mu=mu_U)
               },
               "negbin2" = { 
                 NB_shape <- .get_family_par(family=fam)
                 a <- stats::pnbinom(y - 1L, size=NB_shape, mu=mu_U)
                 b <- stats::pnbinom(y, size=NB_shape, mu=mu_U)
               },
               "negbin" = { 
                 NB_shape <- .get_family_par(family=fam)
                 a <- stats::pnbinom(y - 1L, size=NB_shape, mu=mu_U)
                 b <- stats::pnbinom(y, size=NB_shape, mu=mu_U)
               },
               stop("family not handled")
        )
        u <- runif(n = length(y), min = a, max = b)
        if (identical(fam$zero_truncated,TRUE)) {
          p0 <- switch(famfam,
                       "poisson" = {exp(-mu_U)},
                       "negbin1" = {.negbin1_p0(mu_U,NB_shape)},
                       "negbin2" = {.negbin2_p0(mu_U,NB_shape)},
                       "negbin" = {.negbin2_p0(mu_U,NB_shape)},
                       stop("family zero-truncation not handled")
          )
          u <- (u-p0)/(1-p0) 
        }
        # A typical problem is y>> mu_U, a=b=1 and qnorm(u)=Inf for u ~ > 1-5e-17
        qres <- list(norm=qnorm(u))
      }
    )
    qres  
  }
}

gof <- function(object, method="RQR", plot.=FALSE, testfn=stats::shapiro.test, ...) {
  if (method %in% c("RQR","std_RQR")) {
    RQR <- residuals(object, type=method)
    if (any(is.infinite(RQR))) {
      message("Infinite quantile residuals found, indicating large outliers.")
      RQR <- pmax(pmin(RQR,.Machine$double.xmax),-.Machine$double.xmax) # otherwise shapiro.test() returns confusing NaN / NA.
    }
    if (length(RQR)>5000L && missing(testfn)) {
      message("Too many values for shapiro.test(): using ks.test() instead.")
      goftest <- stats::ks.test(RQR, y="pnorm", ...)
    } else goftest <- testfn(RQR, ...)
    goftest$RQR <- RQR
    if (plot.) plot(object, form=RQR, res_type=method, which="mean")
  }
  goftest
}
