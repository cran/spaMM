.bloc_lambda <- function(HL, SEMblob=NULL, loopout_blob, processed, lam_fix_or_outer_or_NA, cum_n_u_h=processed$cum_n_u_h,
                         models=processed$models, next_LMatrices=loopout_blob$LMatrices,
                         calcRanefPars_blob=loopout_blob$calcRanefPars_blob, ## contains $next_lambda_est (one step ahead of the last updated lambda_est)
                         lev_lambda=loopout_blob$leverages$ranef) {
  nrand <- length(processed$ZAlist)
  rand_to_glm_map <- integer(nrand)
  resglm_lambdaS <- list()
  lambda_pred_list <- as.list(rep(NA,nrand)) ## to be _partially filled_ by this function uing available glm's
  ## je peux avoir SEM sans adjacency (SEM-Matern) et adjacency sans SEM (Poisson-adjacency)
  if (all(models[["lambda"]]=="lamScal")) { # includes ranCoefs
    ####### includes SEM
    if (HL[1]=="SEM") {
      glm_lambda <- SEMblob$glm_lambda
      if ( ! is.null(glm_lambda)) {
        attr(glm_lambda,"whichrand") <- done <- 1L
        resglm_lambdaS[["SEM"]] <- glm_lambda
        ## following syntax OK for all adjacency case or for ~1 (but not random slope)
        lambda_pred_list[attr(glm_lambda,"whichrand")] <- predict(glm_lambda,newdata=data.frame("X.Intercept."=1,adjd=0),type="response")[1L] ## le [1L] en cas d'offset... 
      }
    } else {
      ## checks whether a previous resglm was computed for adjacency model  
      glm_lambda <- calcRanefPars_blob$glm_lambda
      if ( ! is.null(glm_lambda)) { ## includes this adjacency fit
        done <- attr(glm_lambda,"whichrand") ## index according to ordering of attr(ZAlist,"exp_ranef_strings")
        resglm_lambdaS[["adjacency_from_calcRanefPars_blob"]] <- glm_lambda
        lambda_pred_list[attr(glm_lambda,"whichrand")] <- predict(glm_lambda,newdata=data.frame("X.Intercept."=1,adjd=0),type="response")[1L] ## le [1L] en cas d'offset... 
      } else done <- NULL
    }
    rand_to_glm_map[done] <- length(resglm_lambdaS)
    ## next builds a resglm for all other random effects
    notdone <- setdiff(seq_len(nrand),done)
    if (length(notdone)) {
      cum_Xi_cols <- cumsum(c(0,attr(processed$ZAlist,"Xi_cols")))
      for(it in notdone) { ## CAR here if old corrHLfit method (or rho fixed?). But ranCoefs typically...
        if (is.na(lam_fix_or_outer_or_NA[it])) {
          colrange <- (cum_Xi_cols[it]+1L):cum_Xi_cols[it+1L] 
          u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
          loclamdata <- data.frame(processed$X_lamres[u.range,colrange,drop=FALSE]) # virtually these 'datat' are the design matrix of the GLM 
          # which is ans is incidence matrix for each level of the grouping factor of the ranCoef affecting each 'response' (here the leverages...)
          ## part of the problem here is that we need a formula call not an ~X-1 call
          ## part of the problem is that in random-slope model, the Intercept column is not constant  (<= long after: unclear comment)
          colnames(loclamdata)[colnames(loclamdata)=="(Intercept)"] <-"\\(Intercept\\)" 
          ## FR->FR in a later version the formula should be provided by processed$ ?
          loclamformula <- as.formula(paste("~",paste(colnames(loclamdata),collapse="+"),"-1")) # one factor for eahc column of the incidence design matrix
          locetastart <- log(calcRanefPars_blob$next_lambda_est[u.range]) 
          ## !! For random-slope models, locetastart is not the solution, as the next_lambda_est do not store the final variances
          resp_lambda <- calcRanefPars_blob$resp_lambda
          locarglist <- list(formula=loclamformula, dev.res=resp_lambda[u.range],
                             lev=lev_lambda[u.range],
                             data=loclamdata,
                             control=processed[["control.glm"]],
                             etastart=locetastart
          )
          if (NCOL(loclamdata) == 1L && all(loclamdata==1) 
              && length(start <- unique(locetastart))==1L # the last test bc currently .glm_reformat does not handle other cases, 
              #        not all excluded by the two previous tests in the case a prior_lam_fac was used 
              #        (very contrived... see also comments on dev.res argument of .calc_dispGammaGLM).
              # To use .glm_reformat() in the prior_lam_fac case, one would have to get the unique_lambda, which is not in calcRanefPars_blob, as 'start' value.
              # A test of all this is test-negbin1.R, where the final estimate 1.4710451324 is the unique_lambda, 
              #                                             the next_lambda_est are unique_lambda*(prior_lam_fac="wei" variable) 
              #                                             and the locetastart are their log().
              ) {
            # we practically have the answer and reformat it; the .calc_dispGammaGLM() first tries the $method <- ".glm_reformat", and if it fails, 
            # uses spaMM_glm.fit()
            locarglist$method <- ".glm_reformat"
            locarglist$start <- start
            #!# locarglist$control <- list(maxit=1L)
          } else {
            ## includes random-slope
            ## .calc_dispGammaGLM() -> spaMM_glm.fit() should always provide the final glm in these cases.
          }
          glm_lambda <- do.call(".calc_dispGammaGLM",locarglist)
          attr(glm_lambda,"whichrand") <- it
          glmname <- paste(it,"from_post_fit",sep="_")
          resglm_lambdaS[[glmname]] <- glm_lambda
          rand_to_glm_map[it] <- length(resglm_lambdaS)  ## gives the position of just-added glm
          ## following syntax OK for for ~1 and for random-slope (but not adjacency)
          # [single] <- (multiple) is not correctly interpreted
          for (lit in attr(glm_lambda,"whichrand")) lambda_pred_list[[lit]] <- predict(glm_lambda,
                                                                                    newdata=unique(glm_lambda$data),type="response")
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
  process_resglm_blob$lambda_pred_list <- lambda_pred_list
  return(process_resglm_blob)
} 
