pdep_effects <- function(object, focal_var, newdata =object$data, length.out=20L, focal_values=NULL, 
                         level=0.95, levels = NULL, submodel=NULL,
                         intervals = "predVar", indiv=FALSE, verbose=NULL, ...) {
  verbose <- .modify_list(list(na_once=TRUE), verbose)
  was_invColdoldList_NULL <- is.null(object$envir$invColdoldList) # to be able to restore initial state 
  if (inherits(object,"fitmv") && is.null(submodel)) 
    stop("'submodel' argument required for multivariate-response fits.")
  
  if (!focal_var %in% colnames(newdata)) {
    stop("'focal_var' is not found in the data.")
  }
  ori_values <- newdata[, focal_var,drop=TRUE] # [,drop] bc tibbles do not automatically drop 
  ori_values <- na.omit(ori_values)
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
                    level=level, verbose=verbose, submodel=submodel,
                    ...)
    CIs <- attr(pred,"intervals") ## not intervals <- ... within the loop!... as this would modify the argument of predict()
    if ( ! is.null(submodel)) {
      cumnobs <- cumsum(c(0L,attr(pred,"nobs"))) 
      predrange <- .subrange(cumnobs, submodel)
      minmax <- range(predrange)
      pred <- attr(pred,"mv")[[submodel]]
      dim(pred) <- c(length(pred),1L)
      CIs <- CIs[predrange,,drop=FALSE]
      attr(resu,"range") <- minmax
    }
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
  environment(.warn_NA_in_newdata)$NA_in_newdata_NOT_warned <- TRUE
  return(resu)
}

#pdep_effects(simple1_ML,"diamZ")

plot_effects <- function(object, focal_var, newdata=object$data, # doc as a data frame, but a matrix may be sufficient
                         focal_values=NULL, effects=NULL, submodel=NULL,
                        xlab = focal_var, ylab=NULL, rgb.args=col2rgb("blue"), add=FALSE, ylim=NULL, 
                        ...) {
  # If focal_var remains NULL, the idea is probably to run over all predictor variables (not all regressors), 
  #             but this entails other graphic decisions... 
  if (is.null(effects)) effects <- pdep_effects(object, newdata=newdata, focal_var=focal_var, 
                                                indiv=FALSE, focal_values=focal_values, submodel=submodel,
                                                ...) # 'predict on hacked values'
  # : could imagine plotting the results of indiv=TRUE (requires more code)
  if (family(object, submodel=submodel)$family %in% c("binomial","betabin")) {
    resp <- object$y/object$BinomialDen
    if (is.null(ylab)) {
      form <- formula.HLfit(object,which="")
      if (paste(form[[2L]])[[1L]]=="cbind") {
        ylab <- paste("frequency(",form[[2L]][[2L]],")")
      } else ylab <- paste("frequency(",form[[2L]],")")
    }
  } else {
    # there are predictions even when response values were missing in the data
    # In 'byP3' pois4mlogit 1st submodel for example there are 8 fitted values but  predict can generate 9 ones...
    # This block deals with the 8 fitted values
    resp <- object$y
    if ( ! is.null(submodel)) {
      cum_nobs <- attr(object$families,"cum_nobs")
      yrange <- .subrange(cum_nobs, submodel)
      resp <- resp[yrange] #   8 values
    }
    if (inherits(object,"pois4mlogit")) {
      mnsizes <- object$p4m_info$multinom_info[["mnsizes"]] # 10 including one NA
      mnsizes <- mnsizes[object$p4m_info$multinom_info$mnpos_in_template[,submodel]] # mnsizes of the 8 fitted values
      resp <- resp/mnsizes
    }
    if (is.null(ylab)) {
      if (inherits(object,"pois4mlogit")) {
        form <- formula.HLfit(object,which="")[[submodel]]
        ylab <- paste("frequency(",form[[2L]],")")
      } else if (inherits(object,"fitmv")) {
        form <- formula.HLfit(object,which="")[[submodel]]
        ylab <- paste(form[[2L]])
      } else ylab <- paste(formula.HLfit(object,which="")[[2L]])}
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
