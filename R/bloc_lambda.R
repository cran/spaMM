.bloc_lambda <- function(models, #init.lambda, 
                         SEMblob=NULL, 
                        calcRanefPars_blob=NULL, ## contains $next_lambda_est (one step ahead of the last updated lambda_est)
                        processed, lambda.Fix, cum_n_u_h, lev_lambda=NULL, next_LMatrices) {
  nrand <- length(processed$ZAlist)
  rand_to_glm_map <- integer(nrand)
  resglm_lambdaS <- list()
  print_lambdas <- as.list(rep(NA,nrand)) ## to be _partially filled_ by this function uing available glm's
  ## je peux avoir SEM sans adjacency (SEM-Matern) et adjacency sans SEM (Poisson-adjacency)
  if (all(models[[2]]=="lamScal")) { 
    ####### includes SEM
    if ( ! is.null(SEMblob)) {
      glm_lambda <- SEMblob$glm_lambda
      attr(glm_lambda,"whichrand") <- done <- 1L
      resglm_lambdaS[["SEM"]] <- glm_lambda
      ## following syntax OK for all adjacency case or for ~1 (but not random slope)
      print_lambdas[attr(glm_lambda,"whichrand")] <- predict(glm_lambda,newdata=data.frame("X.Intercept."=1,adjd=0),type="response")[1L] ## le [1L] en cas d'offset... 
    } else {
      ## checks whether a previous resglm was computed for adjacency model  
      glm_lambda <- calcRanefPars_blob$glm_lambda
      if ( ! is.null(glm_lambda)) { ## includes this adjacency fit
        done <- attr(glm_lambda,"whichrand") ## index according to ordering of attr(ZAlist,"exp_ranef_strings")
        resglm_lambdaS[["adjacency_from_calcRanefPars_blob"]] <- glm_lambda
        print_lambdas[attr(glm_lambda,"whichrand")] <- predict(glm_lambda,newdata=data.frame("X.Intercept."=1,adjd=0),type="response")[1L] ## le [1L] en cas d'offset... 
      } else done <- NULL
    }
    rand_to_glm_map[done] <- length(resglm_lambdaS)
    ## next builds a resglm for all other random effects
    notdone <- setdiff(seq(nrand),done)
    if (length(notdone)>0) {
      cum_Xi_cols <- cumsum(c(0,attr(processed$ZAlist,"Xi_cols")))
      for(it in notdone) { ## CAR here if old corrHLfit method (or rho fixed?)
        if (is.na(lambda.Fix[it])) {
          colrange <- (cum_Xi_cols[it]+1L):cum_Xi_cols[it+1L] 
          u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
          loclamdata <- data.frame(processed$X_lamres[u.range,colrange,drop=FALSE])
          ## part of the problem here is that we need a formula call not an ~X-1 call
          ## part of the problem is that in random-slope model, the Intercept column is not constant  
          colnames(loclamdata)[colnames(loclamdata)=="(Intercept)"] <-"\\(Intercept\\)" 
          ## FR->FR in a later version the formula should be provided by processed$ ?
          loclamformula <- as.formula(paste("~",paste(colnames(loclamdata),collapse="+"),"-1"))
          locetastart <- log(calcRanefPars_blob$next_lambda_est[u.range]) 
          ## !! For random-slope models, this is not the solution, as the next_lambda_est do not store the final variances
          resp_lambda <- calcRanefPars_blob$resp_lambda
          locarglist <- list(formula=loclamformula, dev.res=resp_lambda[u.range],
                             lev=lev_lambda[u.range],
                             data=loclamdata,
                             control=processed[["control.glm"]],
                             etastart=locetastart
          )
          if (NCOL(loclamdata) == 1L && all(loclamdata==1)) { ## 02/2016
            # we practically have the answer and reformat it
            locarglist$method <- ".glm_reformat"
            locarglist$start <- unique(locetastart)
            #!# locarglist$control <- list(maxit=1L)
          } else {
            ## includes random-slope
            #!# locarglist$try <- TRUE 
            ## 'try' removed, 04/2016: spaMM_glm.fit  should always provide a final glm
          }
          glm_lambda <- do.call(".calc_dispGammaGLM",locarglist)
          attr(glm_lambda,"whichrand") <- it
          glmname <- paste(it,"from_post_fit",sep="_")
          resglm_lambdaS[[glmname]] <- glm_lambda
          rand_to_glm_map[it] <- length(resglm_lambdaS)  ## gives the position of just-added glm
          ## following syntax OK for for ~1 and for random-slope (but not adjacency)
          # [single] <- (multiple) is not correctly interpreted
          for (lit in attr(glm_lambda,"whichrand")) print_lambdas[[lit]] <- unique(predict(glm_lambda,type="response"))
        } ## else no glm for outer/fixed lambda
      }
    }
    attr(resglm_lambdaS,"rand.families") <- processed$rand.families ## FR->FR util ?
  } else {
    stop("From HLfit: 'lamHGLM' and 'lamGLM' not fully implemented.")
    ## there is a template code with comments in version 260812 of HLfit
  }
  ## now reformat all resglm's results
  process_resglm_blob <- .process_resglm_list(resglm_lambdaS,nrand=nrand)
  process_resglm_blob$rand_to_glm_map <- rand_to_glm_map
  process_resglm_blob$print_lambdas <- print_lambdas
  return(process_resglm_blob)
} 
