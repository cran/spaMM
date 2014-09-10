validData <-
function(formula,resid.formula=NULL,data) {
  formula <- asStandardFormula(formula) ## removes spatial tags
  frame.form <- subbarsMM(formula) ## this comes from lme4 and converts (...|...) terms to some "+" form
  if (!is.null(resid.formula)) {
    resid.formula <- asStandardFormula(resid.formula) ## removes spatial tags
    frame.resid.form <- subbarsMM(resid.formula) ## this comes from lme4 and converts (...|...) terms to some "+" form
    frame.form <- paste(DEPARSE(frame.form),"+",DEPARSE(frame.resid.form[[2]]))
  }
  frame.form <- as.formula(frame.form)
  environment(frame.form) <- environment(formula)
  mf <- call("model.frame",data=data) ## it adds the formula argument below....
  mf$formula <- frame.form
  mf$drop.unused.levels <- TRUE
  mf <- eval(mf) ## data.frame with many attributes
  return(mf)
}
