cat(crayon::yellow("test aug_ZXy:\n"))

#cat("\ntest no fixef, pw & REMLformula, by augZXy:\n")
if(requireNamespace("lme4", quietly = TRUE)
   && spaMM.getOption("example_maxtime")>0.7) { # but it's faster than that.)
  data("sleepstudy",package = "lme4")
  oldopt <- spaMM.options(allow_augZXy=FALSE)
  (fit <- fitme(Reaction ~ 0 + (1|Subject), data = sleepstudy))
  set.seed(123)
  sleepstudy$pw <- 1 + runif(180)
  (fit <- fitme(Reaction ~ 0 + (1|Subject), data = sleepstudy, prior.weights=pw))
  spaMM.options(allow_augZXy=2)
  (fit <- fitme(Reaction ~ 0 + (1|Subject), data = sleepstudy, prior.weights=pw))
  (fit <- fitme(Reaction ~ 0 + (1|Subject), data = sleepstudy, prior.weights=pw, REMLformula=Reaction ~ 1 + (1|Subject), method="REML"))
  spaMM.options(allow_augZXy=FALSE)
  (fit <- fitme(Reaction ~ 0 + (1|Subject), data = sleepstudy, prior.weights=pw, REMLformula=Reaction ~ 1 + (1|Subject), method="REML"))
  spaMM.options(oldopt) # spaMM.options(allow_augZXy=NULL)
  # works also with .HLfit_body_augZXy_W()
}


if (file.exists((privdata <- "C:/home/francois/travail/stats/spaMMplus/spaMM/package/tests_private/all_fitness.txt"))) {
  my.data <- read.table(privdata, header = TRUE, sep = "\t",dec = ".")
  my.data$line <- factor(as.character(my.data$line))
  my.data <- na.omit(my.data)
  
  # Minimal example (small matrices) for testing ranCoefs code in general: 
  set.seed(666)
  perm <- sample(nrow(my.data))
  if (exists(".HLfit_body_augZXy_invL",envir = asNamespace("spaMM"))) { # actually spprec case...
    oldopt <- spaMM.options(augZXy_fitfn=".HLfit_body_augZXy_invL", check_alt_augZXy=FALSE) # F I X_invL no check_alt_augZXy bc it fails. But the final logLik is checked
    (mini_rC <- fitme(total_red ~ 1 + (sex|env), data = my.data[perm[1:20],], method="ML"))
    spaMM.options(oldopt)
  } else message(".HLfit_body_augZXy_invL() not included in build => test not performed on pckage istalled from build")
  
  # Plus an augZXy test (not ranCoefs):
  (vanilla <- fitme(total_red ~ sex*env + (1|rep) + (1|line), data = my.data, method="ML"))
  if (exists(".HLfit_body_augZXy_W",envir = asNamespace("spaMM"))) {
    oldopt <- spaMM.options(augZXy_fitfn=".HLfit_body_augZXy_W", check_alt_augZXy=TRUE) 
    essainola <- fitme(total_red ~ sex*env + (1|rep) + (1|line), data = my.data, method="ML")
    spaMM.options(oldopt)
    testthat::expect_true((diff(range(logLik(vanilla),logLik(essainola)))<1e-8)) 
  } else message(".HLfit_body_augZXy_W() not included in build => test not performed on pckage istalled from build")
  if (exists(".HLfit_body_augZXy_invL",envir = asNamespace("spaMM"))) {
    oldopt <- spaMM.options(augZXy_fitfn=".HLfit_body_augZXy_invL", check_alt_augZXy=TRUE)
    essainola <- fitme(total_red ~ sex*env + (1|rep) + (1|line), data = my.data, method="ML")
    spaMM.options(oldopt)
    testthat::expect_true((diff(range(logLik(vanilla),logLik(essainola)))<1e-8)) 
  } else message(".HLfit_body_augZXy_invL() not included in build => test not performed on pckage istalled from build")
}

