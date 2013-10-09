plot.HLfit <-
function(x, which=c(1,2), 
      titles = list(
          meanmodel=list(outer="Mean model",devres="Deviance residuals", absdevres="|Deviance residuals|", resq="Residual quantiles", devreshist="Deviance residuals"), 
          ranef=list(outer="Random effects and leverages",qq="Random effects Q-Q plot", levphi=expression(paste("Leverages for ",phi)), levlambda=expression(paste("Leverages for ",lambda))) 
        ), 
      control = list() , ...) {
	residuals <- x$std_dev_res 
 	fitted.values <- x$fv
	lev_lambda <- x$lev_lambda
	lev_phi <- x$lev_phi
    ranef <- x$ranef ## u
    ##
    pch <- control$pch
    if (is.null(pch)) pch <- "+" 
    pcol <- control$pcol
    if (is.null(pcol)) pcol <- "blue"
    lcol <- control$lcol
    if (is.null(lcol)) lcol <- "red" 
	for (i in which) {
		if (i > 1) dev.new()
		if (i == 1) { ## diagnostic plots for mean model
			par(mfrow = c(2, 2), oma = c( 0, 0, 2, 0 ), pty = "s", ...)
			loess.fit <- loess.smooth(fitted.values, residuals)
			plot(fitted.values, residuals, xlab = "Fitted Values", 
				 ylab = titles$meanmodel$devres, pch = pch, col = pcol, bty = "n", main = titles$meanmodel$devres)
			lines(loess.fit$x, loess.fit$y, col = lcol)
			loess.fit <- loess.smooth(fitted.values, abs(residuals))
			plot(fitted.values, abs(residuals), xlab = "Fitted Values", 
				 ylab = titles$meanmodel$absdevres, pch = pch, col = pcol, bty = "n", main = titles$meanmodel$absdevres)
			lines(loess.fit$x, loess.fit$y, col = lcol)
			qqnorm(residuals, col = pcol, pch = pch, bty = "n", 
				   xlab = "Normal quantiles", ylab = titles$meanmodel$resq, main = titles$meanmodel$resd)
			qqline(residuals, col = lcol)
			hist(residuals, density = 15, xlab = titles$meanmodel$devreshist, main = "", col = pcol)
            title(titles$meanmodel$outer,outer=TRUE)
		}
        nlevplots <- length(c(which(!is.null(lev_lambda)),which(!is.null(lev_phi)),!is.null(ranef)))
		if (i == 2 && nlevplots > 0 ) {
			if (nlevplots<3) {
              par(mfrow = c(1, nlevplots), oma = c( 0, 0, 2, 0 ), pty = "s", ...)
            } else { par(mfrow = c(2, 2), oma = c( 0, 0, 2, 0 ), pty = "s", ...) }
			if (!is.null(ranef)) {
              std_ranef <- ranef/sqrt(1-lev_lambda)
     		  qqnorm(ranef, col = pcol, pch = pch, bty = "n", 
				   xlab = "Normal quantiles", ylab = expression(paste("Standardized ",italic(u)," quantiles")), main = titles$ranef$qq)
	    	  qqline(ranef, col = lcol)
            }
			if (!is.null(lev_phi)) plot(lev_phi, ylab = "", main = titles$ranef$levphi , pch = pch, col = pcol, bty = "n")
			if (!is.null(lev_lambda)) plot(lev_lambda, ylab = "", main= titles$ranef$levlambda , pch = pch, col = pcol, bty = "n")
            title(titles$ranef$outer,outer=TRUE)
		} 
	}
  invisible(x)
}
