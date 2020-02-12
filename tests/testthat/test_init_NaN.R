cat(crayon::yellow("test NaN's in init\n")) 

if (spaMM.getOption("example_maxtime")>1.6) {
  data("blackcap")
  fitme(migStatus ~ means+ Matern(1|latitude+longitude),data=blackcap) # poor
  #  Compare with the following two ways of avoiding outer-optimization of lambda:
  corrHLfit(migStatus ~ means+ Matern(1|latitude+longitude),data=blackcap,
            HLmethod="ML")
  fitme(migStatus ~ means+ Matern(1|latitude+longitude),data=blackcap, 
        init=list(lambda=NaN))
}
