cat("\ntest of bboptim compatibility:")

if (require("blackbox",quietly=TRUE)) {
  example("bboptim") ## detected calcpredVar issue for etaFix   
}

# also feasible but slow:
if ( FALSE && require("Infusion",quietly=TRUE)) {
  Infusion.options(example_maxtime=20)
  example("Infusion") ## detected calc_logdisp_cov issues   
}