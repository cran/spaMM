pdep_effects <- function(object, focal_var, newdata =object$data, length.out=20L, focal_values=NULL, 
                         level=0.95, levels = NULL,
                         intervals = "predVar", indiv=FALSE,...) {
  was_invColdoldList_NULL <- is.null(object$envir$invColdoldList) # to be able to restore initial state 
  if (!focal_var %in% colnames(newdata)) {
    stop("'focal_var' is not found in the data.")
  }
  ori_values <- newdata[, focal_var,drop=TRUE] # drop[] bc tibbles do not automatically drop 
  if (is.character(ori_values)) ori_values <- factor(ori_values)
  if (is.logical(ori_values)) {
    if (is.null(focal.values <- focal_values)) focal.values <- c(FALSE,TRUE) 
  } else if (is.factor(ori_values)) {
    if (is.null(levels)) {
      if (is.null(focal_values)) {
        focal.values <- levels(droplevels(ori_values))
      } else {
        if ( ! inherits(focal_values,"factor")) focal_values <- factor(focal_values)
        focal.values <- levels(focal_values)
        if (!all(focal.values %in% levels(droplevels(ori_values)))) {
          stop("Some of the levels of 'focal_values' are absent from the fitted data.")
        }
      }
    } else {
      if (!all(levels %in% levels(droplevels(ori_values)))) {
        stop("Some of the 'levels' are absent from the fitted data.")
      }
      focal.values <- levels
    }
    focal.values <- factor(focal.values)
    length.out <- length(focal.values)
  } else if (is.numeric(ori_values)) {
    if (is.null(focal.values <- focal_values)) {
      focal.values <- ori_values
    } else if (missing(length.out)) length.out <- 0L # reverse default in that case.
    if (length.out && diff(range(focal.values))>0) {
      focal.values <- seq(min(focal.values),
                          max(focal.values),
                          length.out = length.out)
    } else focal.values <- sort(unique(focal.values))
  } else {
    stop("Unhandled class for 'focal_var'.")
  }
  if (indiv) {
    resu <- lapply(focal.values,function(v)list(focal_var=v))
  } else {
    resu <- data.frame(focal_var = focal.values)
    resu$pointp <- resu$low <- resu$up <- numeric(nrow(resu))
  }
  for (it in seq_along(focal.values)) {
    newdata[,focal_var] <- focal.values[it]
    pred <- predict(object,newdata,intervals = intervals, control=list(fix_predVar=NA), 
                    level=level, ...)
    CIs <- attr(pred,"intervals") ## not intervals <- ... within the loop!...
    if (indiv) {
      resu[[it]]$pointp <- pred[,1]
      resu[[it]]$low <- CIs[,1]
      resu[[it]]$up <- CIs[,2]
    } else {
      resu$pointp[it] <- mean(pred[,1])
      resu$low[it] <- mean(CIs[,1])
      resu$up[it] <- mean(CIs[,2])
    }
  }
  if (was_invColdoldList_NULL) object$envir$invColdoldList <- NULL
  return(resu)
}

#pdep_effects(simple1_ML,"diamZ")

plot_effects <- function(object, focal_var, newdata=object$data, # doc as a data frame, but a matrix may be sufficient
                         focal_values=NULL, effects=NULL, 
                        xlab = focal_var, ylab=NULL, rgb.args=col2rgb("blue"), add=FALSE, ylim=NULL, ...) {
  # If focal_var remains NULL, the idea is probably to run over all predictor variables (not all regressors), 
  #             but this entails other graphic decisions... 
  if (is.null(effects)) effects <- pdep_effects(object, newdata=newdata, focal_var=focal_var, 
                                                indiv=FALSE, focal_values=focal_values, ...) # 'predict on hacked values'
  # : could imagine plotting the results of indiv=TRUE (requires more code)
  if (object$family$family %in% c("binomial","betabin")) {
    resp <- object$y/object$BinomialDen
    if (is.null(ylab)) {
      form <- formula.HLfit(object,which="")
      if (paste(form[[2L]])[[1L]]=="cbind") {
        ylab <- paste("frequency(",form[[2L]][[2L]],")")
      } else ylab <- paste("frequency(",form[[2L]],")")
    }
  } else {
    resp <- object$y
    if (is.null(ylab)) ylab <- paste(formula.HLfit(object,which="")[[2L]])
  }
  if (is.null(ylim)) ylim <- stats::quantile(resp, c(0.025, 0.975))
  rgb.args <- as.list(rgb.args)
  if (is.null(rgb.args$maxColorValue)) rgb.args$maxColorValue <- 255
  colpts <- do.call(rgb,rgb.args)
  rgb.args$alpha <- rgb.args$maxColorValue*0.1 ## should be able to control the alpha factor
  colshd <- do.call(rgb,rgb.args)
  new.x <- effects[,"focal_var"]
  if (inherits(new.x,"factor")) {asnum.x <- seq_along(new.x)} else asnum.x <- new.x
  if (!add) {
    if (inherits(new.x,"factor")) {
      plot(NULL, type = "l", xaxt="n",
           ylab = ylab, xlab = xlab, xlim = range(asnum.x), ylim = ylim)
      axis(1,at=asnum.x, labels=new.x)
    } else plot(NULL, type = "l",
       ylab = ylab, xlab = xlab, xlim = range(new.x), ylim = ylim)
  }
  polygon(c(asnum.x,rev(asnum.x)),c(effects[,"low"],rev(effects[,"up"])),border=NA,col=colshd)
  if (inherits(new.x,"numeric")) {
    lines(new.x,effects[,"pointp"],lwd=2,col=colpts)
  } else points(asnum.x,effects[,"pointp"],col=colpts,pch=19)
  graphics::rug(effects[, "focal_var"], side = 1, col = colpts) ## should be able to control the side
  graphics::rug(resp, side = 2, col = colpts, quiet= ! is.null(focal_values)) ## idem
  invisible(effects)
}

#plot_effects(simple1_ML,focal_var="density")
#plot_effects(simple1_ML,newdata=simple1_ML$data[ ! isMale,],focal_var="density", rgb.args=col2rgb("red"))
#plot_effects(simple1_ML,newdata=simple1_ML$data[ isMale,],focal_var="density", add=TRUE)
