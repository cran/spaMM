if (Sys.getenv("_LOCAL_TESTS_")=="TRUE") { ## set in etc/Renviron.site (cf R Windows FAQ)
  if(requireNamespace("testthat", quietly = TRUE)) {
    pkg   <- "spaMM"
    require(pkg, character.only=TRUE, quietly=TRUE)
    # options(error = quote({dump.frames(to.file = TRUE)})) # useful for bugs in .do_TRACE()
    if (interactive()) {
      # options(error=recover)
      # spaMM.options(use_ZA_L=NULL)
      abyss <- matrix(runif(2e7),nrow=1000); gc(reset=TRUE) ## partial control of gc trigger...
      while (dev.cur()>1L) dev.off()
      op <- devAskNewPage(ask=FALSE)
      testfiles <- dir(paste0(projpath(),"/package/tests/testthat/"),pattern="*.R",full.names = TRUE)
      # oldmaxt <- spaMM.options(example_maxtime=60)
      timings <- t(sapply(testfiles, function(fich){system.time(source(fich))}))
      # spaMM.options(oldmaxt)
      print(sums <- colSums(timings))
      ## testthat::test_package(pkg) ## for an installed package
      if (FALSE) { ## tests not included in package (using unpublished data, etc.)
        priv_testfiles <- dir(paste0(projpath(),"/package/tests_private/"),pattern="*.R",full.names = TRUE)
        priv_testfiles <- setdiff(priv_testfiles,paste0(projpath(),"/package/tests_private/knit_LM2GLMM.R"))
        priv_timings <- t(sapply(priv_testfiles, function(fich){
          cat(crayon::green(paste0("\n",fich)))
          system.time(try(source(fich)))
        }))
        #spaMM.options(oldmaxt)
        print(colSums(priv_timings))
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
