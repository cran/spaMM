checkRandLink <-
function(rand.family) {
  if (class(rand.family)=="family") { ## then check what is feasible
    lcrandfamfam <- tolower(rand.family$family) ## tolower once and for all
    oklink <- F
    ## cases where g(u)=th(u)
    if (lcrandfamfam=="gaussian" && rand.family$link=="identity") oklink <- T          
    if (lcrandfamfam=="gamma" && rand.family$link=="log") oklink <- T          
    if (lcrandfamfam=="inverse.gamma" && rand.family$link=="-1/mu") oklink <- T
    if (lcrandfamfam=="beta" && rand.family$link=="logit") oklink <- T
    ## cases where g(u)!=th(u)
    if (lcrandfamfam=="inverse.gamma" && rand.family$link=="log") oklink <- T 
    if ( ! oklink) {
      allowed <- switch(lcrandfamfam,
                        gaussian= "is 'identity'",
                        gamma= "is 'log'",
                        beta= "is 'logit'",
                        "inverse.gamma" = "are '-1/mu' and 'log'"
                        )
      mess <- paste("(!) rand.family/link combination not handled;\nallowed link(s) for rand.family '",rand.family$family,"' ",allowed,sep="")
      stop(mess)
    }
  } else { ## rand.family is a string ## preprocess -> checkRandLinkS -> here : OK
    lcrandfamfam<-tolower(rand.family) ## tolower once and for all
    rand.family <- switch(lcrandfamfam,
                          gaussian = gaussian(),
                          gamma = Gamma(link=log), ## NOT the default link
                          beta = Beta(), 
                          "inverse.gamma" = inverse.Gamma(),
                          stop("rand.family argument not valid")
    )
  }
  lcrandfamfam
}
