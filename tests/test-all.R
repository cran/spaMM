# tools::check_packages_in_dir(".", check_args = "--no-multiarch")
# devtools::check(".", args ="--no-multiarch", manual=TRUE)
if (Sys.getenv("_LOCAL_TESTS_")=="TRUE") { ## set in <R_HOME>/etc/Renviron.site (cf R Windows FAQ) or by usethis::edit_r_environ() for 'user' one (in C:/Documents)
  # BUT there is a problem... see AAA_Renviron_mess.R 
  if(requireNamespace("testthat", quietly = TRUE)) {
    pkg   <- "spaMM"
    require(pkg, character.only=TRUE, quietly=TRUE)
    # options(error = quote({dump.frames(to.file = TRUE)})) # useful for bugs in .do_TRACE()
    if (interactive()) {
      # **** see # see AAA_install_INLA.memo.txt: to install INLA ****
      # install.packages(c("DHARMa", "inlabru")) # + private pkg "probitgem"
      # + packages in requireNamespace(.)... and install_git[...]
      # + optional but important dependencies of Infusion... aster... hglm...
      # options(error=recover)
      # spaMM.options(use_ZA_L=NULL)
      # abyss <- matrix(runif(2e7),nrow=1000); gc(reset=TRUE) ## partial control of gc trigger...
      {
        # spaMM.options(projpath="D:/home/francois/travail/stats/spaMMplus/spaMM")
        testfiles <- dir(paste0(spaMM::projpath(),"/package/tests/testthat/"),pattern="*.R$",full.names = TRUE)
        #  testfiles <- dir(paste0(spaMM::projpath(),"/package/tests/testthat/"),full.names = TRUE)[-1L] # temporary fix for Rstudio problem
        #  testfiles <- dir(paste0("D:/home/francois/travail/stats/spaMMplus/atuer","/package/tests/testthat/"),pattern="*.R",full.names = TRUE)
        #  testfiles <- dir(paste0("D:/home/francois/travail/stats/spaMMplus/spaMM","/package/tests/testthat/"),pattern="*.R",full.names = TRUE)
        #
        tfun <- function(x) {
          gc()# cf doc of system.time(., gcFirst) => but if gc timings are highly variable, gcFirst=TRUE is pointless (and the whole is misleading). 
          system.time(source(x), gcFirst=FALSE)
        }
        while (dev.cur()>1L) dev.off()
        op <- devAskNewPage(ask=FALSE)
        oldWarnOpt <- options(warnPartialMatchArgs = TRUE, # Hornik, R-devel, 2024/04/23
                              warnPartialMatchAttr = TRUE,
                              warnPartialMatchDollar = TRUE)
        # oldmaxt <- spaMM.options(example_maxtime=60)
        timings <- t(sapply(testfiles, function(fich){tfun(fich)}))
        # spaMM.options(oldmaxt)
        options(oldWarnOpt)
        print(sums <- colSums(timings))
      }
      if (FALSE) { # long mv tests, not really for the timings; important tests, mv_nested notably
        if (TRUE) { # 'pattern' should work, but didn't in Rstudio. Fixed in Rstudio version 2022.07.2+576
          extra_testfiles <- dir(paste0(spaMM::projpath(),"/package/tests/testthat/extralong/"),pattern="*.R$",full.names = TRUE)
        } else {
          extra_testfiles <- dir(paste0(spaMM::projpath(),"/package/tests/testthat/extralong/"),full.names = TRUE)
          extra_testfiles <- extra_testfiles[grep("*.R$",extra_testfiles)]
        }
        # extra_testfiles <- dir(paste0(spaMM::projpath(),"/package/tests/testthat/extralong/"),full.names = TRUE)
        extra_timings <- t(sapply(extra_testfiles, function(fich){
          cat(crayon::green(paste0("\n",fich)))
          gc()
          tps <- system.time(chk <- try(source(fich)), gcFirst=FALSE)
          if (inherits(chk,"try-error")) warning(paste0(fich," generated an error"))
          tps
        }))
        print(colSums(extra_timings)) 
      }
      ## testthat::test_package(pkg) ## for an installed package
      if (FALSE) { ## tests not included in package (using unpublished data, etc.)
        cat(crayon::green("Widen the plot panel!\n"))
        # install.packages("FactoMineR")
        # see also includes in tests_private/test-back-compat.R
        if (TRUE) { # see above comment about Rstudio
          priv_testfiles <- dir(paste0(spaMM::projpath(),"/package/tests_private/"),pattern="*.R$",full.names = TRUE)
        } else {
          priv_testfiles <- dir(paste0(spaMM::projpath(),"/package/tests_private/"),full.names = TRUE)
          priv_testfiles <- priv_testfiles[grep("*.R$",priv_testfiles)]
        }
        priv_testfiles <- setdiff(priv_testfiles,paste0(spaMM::projpath(),"/package/tests_private/knit_LM2GLMM.R"))
        priv_timings <- t(sapply(priv_testfiles, function(fich){
          cat(crayon::green(paste0("\n",fich)))
          gc()
          tps <- system.time(chk <- try(source(fich)), gcFirst=FALSE)
          if (inherits(chk,"try-error")) warning(paste0(fich," generated an error"))
          tps
        }))
        print(colSums(priv_timings)) # very roughly 1205.44 s elapsed for default maxtime (0.7) # much less after improving COMP
      }
      if (FALSE) { ## kept separate bc obscure interference as if there was a bug in setTimeLimit()*Rstudio ?
        # abyss <- matrix(runif(2e7),nrow=1000); gc(reset=TRUE) ## partial control of gc trigger... but RESTARTING R appears more efficient.
        useR2021_testfiles <- dir(paste0(spaMM::projpath(),"/package/useR2021/"),pattern="*.R",full.names = TRUE)
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
        save(timings,file=paste0(spaMM::projpath(),"/timings_",packageVersion("spaMM"),"_",sums[[1]],"s.rda"))
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
