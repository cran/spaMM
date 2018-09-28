# Adapted from the code provided by B&G 2009, http://www.math.tau.ac.il/~ybenja/Software/MSFDR.r
MSFDR <- function(  nullfit, fullfit , q = 0.05 , verbose = TRUE) {
  if (is.null(nullfit$REMLformula) || ! identical(attr(nullfit$REMLformula,"isML"),TRUE)) {
    warning("nullfit used REML, hence all further fits will. This is not recommended for inferences about fixed effects.")
  }
  scope <- list(lower=(nullfit$HLframes$fixef_terms), ## only the fixed-effect terms
                upper=(fullfit$HLframes$fixef_terms))
  #trace.stepAIC <- ifelse(print.the.steps , 1,0)
  iter <- 1L
  nslopes_f <- extractAIC(fullfit)[1L] - 1 # B&G comment: 'check if the full model should include the intercept or not !!!!!!'
  current_size <- max(extractAIC(nullfit)[1L]-1, 1)	# B&G comment 'so if the model is with intercept only, the i size won't be 0.'
  # q = .05 # default
  chi2_LR <- qnorm( (1- 0.5 * q* current_size/(nslopes_f+1-current_size*(1-q))  ) )^ 2
  if(verbose ) {print(paste("Starting threshold is: ", chi2_LR)) }
  new_fit <- step(nullfit, direction="forward", scope=scope ,  k = chi2_LR, trace = verbose ) ## stats::step
  new_fit_size <- extractAIC(new_fit)[1] - 1
  
  while(new_fit_size > current_size) {
    iter <- iter + 1L
    if(verbose ) {
      print("=========================================") 
      print(paste0("iteration ", iter , ": new model size is ", new_fit_size)) 
    }
    current_size <-  new_fit_size
    chi2_LR <- qnorm( (1- 0.5 * q* new_fit_size/(nslopes_f+1-new_fit_size*(1-q))  ) )^2
    if(verbose ) {print(paste("new threshold is: ", chi2_LR)) }
    new.fit <- step(new_fit, direction="forward", scope=scope ,  k = chi2_LR, trace = verbose ) ## stats::step
    new_fit_size <- extractAIC(new_fit)[1] - 1
  }	
  return(new_fit)
}
