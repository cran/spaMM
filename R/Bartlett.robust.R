Bartlett.robust <-
function(LRTobject,robust=T,verbose=F) {
      bootLRTS <- with(LRTobject,2*(bootreps[,1]-bootreps[,2]))## full -null but can be p_v or p_bv
      # plot(qchisq(ppoints(zut),1),sort(zut)) ## QQplot, MASS p. 108
      filter <- bootLRTS[bootLRTS>-1e-08]
      resu <- list()
      if (robust) {
        ## finds the df of the distribution by robust regression (MM: MASS p. 161) of QQplot, 
        robustMean <- rlm(sort(filter)~qchisq(ppoints(filter),1)-1,method="MM",maxit=200)$coefficients[[1]] 
        ## [[1]] to remove name which otherwise finishes as a rowname in a subsequent dataframe       
        # plot(qchisq(ppoints(filter),1),sort(filter)) ## QQplot, MASS p. 108
        # points(qchisq(ppoints(filter),1),robustMean * qchisq(ppoints(filter),1),pch=".")   
        if (class(robustMean)=="try-error") {
          resu$robustMean <- NA
          mess <- pastefrom("problem in computation of robustMean.")
          warning(mess)
        } else resu$robustMean <- robustMean
      }
      resu$meanPosbootLRT <- mean(filter)
      resu$nPosbootLRT <- length(filter) 
      if (verbose) print(unlist(resu))
      return(resu)
}
