bessel_lnKnu <-
function (nu, x, give = FALSE, strict = TRUE) { ## from bessel_lnKnu in gsl package
    jj <- HL_process.args(nu, x)
    nu.vec <- jj$arg1
    x.vec <- jj$arg2
    attr <- jj$attr
    jj <- .C("bessel_lnKnu_e", as.double(nu.vec), as.double(x.vec), 
        as.integer(length(x.vec)), val = as.double(x.vec), err = as.double(x.vec), 
        status = as.integer(0 * x.vec), PACKAGE = "spaMM")
    val <- jj$val
    err <- jj$err
    status <- jj$status
    attributes(val) <- attr
    attributes(err) <- attr
    attributes(status) <- attr
    if (strict) {
        val <- HL_strictify(val, status)
    }
    if (give) {
        return(list(val = val, err = err, status = status))
    }
    else {
        return(val)
    }
}
