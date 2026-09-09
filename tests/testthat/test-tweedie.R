cat(cli::col_yellow("\ntest tweedie:\n"))

{ # toy data 
  set.seed(123)
  y <- rgamma(20,shape=rpois(20,lambda=rgamma(20,5)))
  x <- 1:20
  toyTw <- data.frame(x=x,y=y, grp=x)
}

{ # Examples from the doc

  ### Fits
  fitme(y~x+(1|grp),family=spaMM::tweedie(), data=toyTw) 
  
  # tentative estimation of the link exponent:      
  fitme(y~x+(1|grp),family=spaMM::tweedie(link="fit"), data=toyTw, 
        init=list(lambda=NA, phi=NA, Tw_link=0.5),
        lower=list(Tw_link=0.4), upper=list(Tw_link=0.6)) 
  # (this shows that this _can_ work).
  
  ### Interactions with other packages
  tweedie::logLiktweedie(toyglm <- 
                           glm(y~x,family=spaMM::tweedie(index=1.5,link= -0.5), data=toyTw))
  toyglm$fitted.values
  predict(toyglm, type="response")
  tweedie::logLiktweedie(spaMM_glm(y~x,
                                   family=spaMM::tweedie(index=1.5,link= -0.5), data=toyTw))
}