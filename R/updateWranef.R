updateWranef <-
function(rand.family,lambda,u_h,v_h) {
    dudv <- rand.family$mu.eta(v_h) ## general cf InvGamma with log link rand.family$mu.eta(v_h) = exp(v) =u is du/d(log u)   
## compute w.ranef, the scaled d^2f(v)/dv^2. See Appendix 3 of LeeN01 + my notes for what it implies...
    if (rand.family$family=="gaussian") {
       if (rand.family$link=="identity") {
          V_M <- rand.family$variance(u_h) ##rep(1,length(u_h)) ## GLMMs in general
          dlogWran_dv_h <- rep(0L,length(u_h))
       }
     } else if (rand.family$family=="Gamma") { 
       if (rand.family$link=="log") {
          V_M <- u_h ## V(u), canonical conjugate Gamma as in canonical Poisson Gamma HGLM
          dlogWran_dv_h <- rep(1L,length(u_h))
       } 
    } else if (rand.family$family=="inverse.Gamma") { ## for Gamma HGLM 
       ## the canonical form gives the density of theta(u)
       if (rand.family$link=="log") {
          w.ranef <- as.numeric(1/(u_h * lambda)) ## W1/lambda, W1 computation shown in appendix 3 of LeeN01; also in Noh and Lee's code.
          dlogWran_dv_h <- rep(-1L,length(u_h)) ## v=log u, dlogW/dv= dlog(1/u)/dv=-1
          return(list(w.ranef=w.ranef,dlogWran_dv_h=dlogWran_dv_h,dvdu=1/dudv))  ###### return here !
       } else if (rand.family$link=="-1/mu") {
          ## D[\[Nu] (\[Theta][u] - (-Log[-\[Theta][u]])), {\[Theta][u], 2}]
          V_M <- rand.family$variance(u_h) ## u_h^2 ## V(u), canonical conjugate HGLM 
          dlogWran_dv_h <- 2 * u_h ## no independent check 
       }
    } else if (rand.family$family=="Beta") {
       if (rand.family$link=="logit") {
          V_M <- rand.family$variance(u_h) ##  u_h*(1-u_h) ## canonical conjugate HGLM
          dlogWran_dv_h <- 1 - 2 * u_h ## D[Log[u (1 - u)] /. u -> 1/(1 + E^-v), v] /. v -> Log[u/(1 - u)] ; no independent check
       }
    }
    w.ranef <- as.numeric((dudv^2)/(V_M*lambda)) ## semble valide quand v=g(u) = th(u)
return(list(w.ranef=w.ranef,dlogWran_dv_h=dlogWran_dv_h,dvdu=1/dudv))
}
