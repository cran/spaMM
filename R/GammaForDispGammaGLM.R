GammaForDispGammaGLM <-
function (link = "inverse") {
    linktemp <- substitute(link)
    if (!is.character(linktemp)) 
        linktemp <- deparse(linktemp)
    okLinks <- c("inverse", "log", "identity")
    if (linktemp %in% okLinks) 
        stats <- make.link(linktemp)
    else if (is.character(link)) 
        stats <- make.link(link)
    else {
        if (inherits(link, "link-glm")) {
            stats <- link
            if (!is.null(stats$name)) 
                linktemp <- stats$name
        }
        else {
            stop(gettextf("link \"%s\" not available for gamma family; available links are %s", 
                linktemp, paste(sQuote(okLinks), collapse = ", ")), 
                domain = NA)
        }
    }
    variance <- function(mu) mu^2
    validmu <- function(mu) all(mu > 0)
    dev.resids <- function(y, mu, wt) {
      pmax(-2 * wt * (log(ifelse(y == 0, 1, y/mu)) - (y - mu)/mu),1e-16)  ## FR: added the pmax
    }
    aic <- function(y, n, mu, wt, dev) {
        n <- sum(wt)
        disp <- dev/n
        -2 * sum(dgamma(y, 1/disp, scale = mu * disp, log = TRUE) * 
            wt) + 2
    }
    initialize <- expression({
        if (any(y <= 0)) stop("non-positive values not allowed for the gamma family") 
        n <- rep.int(1, nobs)
        mustart <- y
    })
    simfun <- function(object, nsim) {
        wts <- object$prior.weights
        if (any(wts != 1)) 
            message("using weights as shape parameters")
        ftd <- fitted(object)
        shape <- MASS::gamma.shape(object)$alpha * wts ## MASS:::gamma.shape.glm
        rgamma(nsim * length(ftd), shape = shape, rate = shape/ftd)
    }
    structure(list(family = "Gamma", link = linktemp, linkfun = stats$linkfun, 
        linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, 
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
        validmu = validmu, valideta = stats$valideta, simulate = simfun), 
        class = "family")
}
