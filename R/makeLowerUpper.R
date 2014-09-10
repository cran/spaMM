makeLowerUpper <-
function(canon.init, ## cf calls: ~ in user scale, must be a full list of relevant params
                           lower,upper, ## ~in transformed scale
                           user.lower=list(),user.upper=list(),
                           corr.model="Matern",nbUnique,
                           ranFix=list(),
                           optim.scale) {
  ## init.optim not further used...
  if (corr.model=="adjacency") { ## adjacency model
    ## no default value, user values are required 
    lower$rho <- user.lower$rho ## no transfo for adjacency model
    upper$rho <- user.upper$rho ## idem
  } else {
    if (corr.model=="AR1") {
      if ( ! is.null(canon.init$ARphi)) {
        ARphi <- user.lower$ARphi
        if (is.null(ARphi)) ARphi <- -0.999999
        lower$ARphi <- ARphi
        ARphi <- user.upper$ARphi
        if (is.null(ARphi)) ARphi <- 0.999999
        upper$ARphi <- ARphi
      }    
    } else { ## then Matern model....
      if (! is.null(canon.init$rho)) {
        rho <- user.lower$rho
        if (is.null(rho)) rho <- canon.init$rho/150
        if (optim.scale=="transformed") {
          lower$trRho <- rhoFn(rho)
        } else lower$rho <- rho
        rho <- user.upper$rho
        if (is.null(rho)) {
          if (inherits(nbUnique,"list")) nbUnique <- mean(unlist(nbUnique))
          rho <- canon.init$rho*2*nbUnique ## The following was a bit too low for experiments with nu=0.5 : 1/(maxrange/(2*nbUnique)) ## nb => unique rows !
          if (optim.scale=="transformed") rho <- rho*.SpaMM$RHOMAX/(1+rho) ## so that it does not exceed RHOMAX
        }
        if (optim.scale=="transformed") {
          upper$trRho <- rhoFn(rho) 
        } else upper$rho <- rho
        rhoForNu <- canon.init$rho
      } else rhoForNu <-ranFix$rho
      if (! is.null(canon.init$nu)) {
        nu <- user.lower$nu
        if (is.null(nu)) nu <- canon.init$nu/100
        if (optim.scale=="transformed") {
          lower$trNu <- nuFn(nu,rhoForNu)
          #print(c(rhoForNu,nu,lower$trNu))
        } else lower$nu <-nu
        nu <- user.upper$nu
        if (is.null(nu)) nu <- canon.init$nu*.SpaMM$NUMAX/(1+canon.init$nu) ## nu should not diverge otherwise it will diverge in Bessel_lnKnu, whatever the transformation used
        if (optim.scale=="transformed") {
          upper$trNu <- nuFn(nu,rhoForNu)
        } else upper$nu <- nu
        #print(c(rhoForNu,nu,upper$trNu))
      }
    } 
    ##### common to the different models except adjacency (because there are several places where NUgget+adjacency is not handled)
    if ( ! is.null(canon.init$Nugget)) {
      lower$Nugget <- 0
      upper$Nugget <- 0.999999
    }
  }
  if (! is.null(canon.init$phi)) {
    phi <- user.lower$phi
    if (is.null(phi)) phi <- canon.init$phi/1000
    lower$trPhi <- dispFn(phi)
    phi <- user.upper$phi
    if (is.null(phi)) phi <- canon.init$phi*1000
    ## if phi is badly initialized then it gets a default which may cause hard to catch problems in the bootstrap...
    upper$trPhi <- dispFn(phi)
  }
  if (! is.null(canon.init$lambda)) {
    lambda <- user.lower$lambda
    if (is.null(lambda)) lambda <- canon.init$lambda/1000
    lower$trLambda <- dispFn(lambda)
    lambda <- user.upper$lambda
    if (is.null(lambda)) lambda <- canon.init$lambda*1000
    upper$trLambda <- dispFn(lambda)
  }
  return(list(lower=lower,upper=upper))
}
