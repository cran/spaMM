cat("\ntest of bboptim compatibility:")

if (require("blackbox",quietly=TRUE)) {
  example("bboptim") ## detected calcpredVar issue for etaFix   
}

# also feasible but slow:
if (require("Infusion",quietly=TRUE)) {
  #Infusion.options(example_maxtime=20)  ## allows basic example
  #Infusion.options(example_maxtime=120)  ## allows refine()
  example("Infusion") ## detected calc_logdisp_cov issues   
}

if (require("IsoriX",quietly=TRUE)) { ## Checks that exports/imports are OK
  IsoriX.options(spaMM.options("example_maxtime"))
  example("isofit") ## to detect e.g. do_TRACE issues
}