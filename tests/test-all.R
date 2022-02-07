# tools::check_packages_in_dir(".", check_args = "--no-multiarch")
# devtools::check(".", args ="--no-multiarch", manual=TRUE)
if (Sys.getenv("_LOCAL_TESTS_")=="TRUE") { ## set in <R_HOME>/etc/Renviron.site (cf R Windows FAQ) or by usethis::edit_r_environ() for 'user' one (in C:/Documents) 
  if(requireNamespace("testthat", quietly = TRUE)) {
    pkg   <- "spaMM"
    require(pkg, character.only=TRUE, quietly=TRUE)
    # options(error = quote({dump.frames(to.file = TRUE)})) # useful for bugs in .do_TRACE()
    if (interactive()) {
      # install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
      # install.packages(c("probitgem", "DHARMa", "inlabru"))
      # options(error=recover)
      # spaMM.options(use_ZA_L=NULL)
      # abyss <- matrix(runif(2e7),nrow=1000); gc(reset=TRUE) ## partial control of gc trigger...
      {
        tfun <- function(x) {
          gc()# cf doc of system.time(., gcFirst) => but if gc timings are highly variable, gcFirst=TRUE is pointless (and the whole is misleading). 
          system.time(source(x), gcFirst=FALSE)
        }
        while (dev.cur()>1L) dev.off()
        op <- devAskNewPage(ask=FALSE)
        testfiles <- dir(paste0(projpath(),"/package/tests/testthat/"),pattern="*.R",full.names = TRUE)
        # oldmaxt <- spaMM.options(example_maxtime=60)
        timings <- t(sapply(testfiles, function(fich){tfun(fich)}))
        # spaMM.options(oldmaxt)
        print(sums <- colSums(timings))
      }
      ## testthat::test_package(pkg) ## for an installed package
      if (FALSE) { ## tests not included in package (using unpublished data, etc.)
        #install.packages(c("FactoMineR"))
        # see also includes in tests_private/test-back-compat.R
        priv_testfiles <- dir(paste0(projpath(),"/package/tests_private/"),pattern="*.R",full.names = TRUE)
        priv_testfiles <- setdiff(priv_testfiles,paste0(projpath(),"/package/tests_private/knit_LM2GLMM.R"))
        priv_timings <- t(sapply(priv_testfiles, function(fich){
          cat(crayon::green(paste0("\n",fich)))
          gc()
          tps <- system.time(chk <- try(source(fich)), gcFirst=FALSE)
          if (inherits(chk,"try-error")) warning(paste0(fich," generated an error"))
          tps
        }))
        print(colSums(priv_timings)) # very roughly 1380 s for default maxtime (0.7)
      }
      if (FALSE) { ## kept separate bc obscure interference as if there was a bug in setTimeLimit()*Rstudio ?
        # abyss <- matrix(runif(2e7),nrow=1000); gc(reset=TRUE) ## partial control of gc trigger... but RESTARTING R appears more efficient.
        useR2021_testfiles <- dir(paste0(projpath(),"/package/useR2021/"),pattern="*.R",full.names = TRUE)
        useR2021_timings <- t(sapply(useR2021_testfiles, function(fich){
          cat(crayon::green(paste0("\n",fich)))
          gc()
          tps <- system.time(chk <- try(source(fich)), gcFirst=FALSE)
          if (inherits(chk,"try-error")) warning(paste0(fich," generated an error"))
          tps
        }))
        print(colSums(useR2021_timings))
      }
      devAskNewPage(op)
      if (FALSE) {
        save(timings,file=paste0(projpath(),"/timings_",packageVersion("spaMM"),"_",sums[[1]],"s.rda"))
      }
    } else if (FALSE) { ## for R CMD check (but still assuming _LOCAL_TESTS_), but this does not work on nested files
      library("testthat") # cf ?test_check for using library() here:
      library(pkg, character.only = TRUE)
      oldmaxt <- spaMM.options(example_maxtime=60) ## then slow (Rstudio -> devtools tests) 
      report <- test_check(pkg) 
      spaMM.options(oldmaxt)
      print(warnings()) # TODO? catch most of these by expect_warning(..)
    }
  } else {
    cat( "package 'testthat' not available, cannot run unit tests\n" )
  }
}
