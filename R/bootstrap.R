#' Function for one bootstrap sample
#'
#' @param data original data to bootstrap
#' @param Y_name growth outcome variable name
#' @param Z_name vaccination variable name
#' @param X_name covariate name(s)
#' @param S_name infection variable name
#' @param est character vector of names of estimators to use for growth effect
#' @param ml boolean to use SuperLearner models, default FALSE
#' @param Y_X_S1_model optional specify model to be used for fitting growth on covariates in the infected, otherwise growth on all covariates
#' @param Y_X_S0_model optional specify model to be used for fitting growth on covariates in the uninfected, otherwise growth on all covariates
#' @param S_X_model optional specify model to be used for fitting infection on covariates, otherwise infection on all covariates
#' @param Y_X_library optional specify SuperLearner libraries for model fitting growth on covariates, default glm
#' @param S_X_library optional specify SuperLearner libraries for model fitting infection on covariates, default glm
#' @param family family for outcome model, defaults to gaussian for growth
#' @param v_folds number of cross validation folds for SuperLearner, default 3
#' @param effect_dir direction of beneficial effect, defaults to "positive" for beneficial outcome. Used for one-side tests of bounds.  
#' @param epsilon a vector of values for the sensitivity parameter
#' 
#' @returns list containing results for specified estimators on single bootstrap sample
one_boot <- function(
    data,
    Y_name = "Y",
    Z_name = "Z",
    X_name = "X",
    S_name = "S", 
    estimand = c("nat_inf", "doomed", "pop"),
    method = c("gcomp", "ipw", "aipw", "tmle", "bound", "sens"),
    exclusion_restriction = FALSE,
    cross_world = TRUE,
    ml = FALSE, 
    Y_Z_X_model = NULL,
    Y_X_S1_model = NULL,
    Y_X_S0_model = NULL,
    S_X_model = NULL,
    S_Z_X_model = NULL,
    Z_X_model = NULL,
    Y_Z_X_library = c("SL.glm"),
    Y_X_library = c("SL.glm"),
    S_X_library = c("SL.glm"),
    S_Z_X_library = c("SL.glm"),
    Z_X_library = c("SL.mean"),
    family = "gaussian",
    v_folds = 3,
    effect_dir = "positive",
    epsilon = exp(seq(log(0.5), log(2), length = 50)),
    max_resample = 10,
    return_se = TRUE,
    two_part_model = FALSE
){
  
  n <- dim(data)[1]
  boot_row_idx <- sample(1:n, replace=TRUE)
  boot_data <- data[boot_row_idx,]
  
  # If there are more infections in vaccine arm than placebo arm, resample up to max_resample 
  rhobar_0_n <- mean(boot_data[[S_name]][boot_data[[Z_name]] == 0])
  rhobar_1_n <- mean(boot_data[[S_name]][boot_data[[Z_name]] == 1])
  resample <- 0
  
  while(rhobar_0_n <= rhobar_1_n & resample <= max_resample){
    boot_row_idx <- sample(1:n, replace=TRUE)
    boot_data <- data[boot_row_idx,]
    
    rhobar_0_n <- mean(boot_data[[S_name]][boot_data[[Z_name]] == 0])
    rhobar_1_n <- mean(boot_data[[S_name]][boot_data[[Z_name]] == 1])
    resample <- resample + 1
  }
  
  if(resample > max_resample){
    stop(paste0("Exceeded max_resample of ", max_resample, "for given bootstrap replicate"))
  }
  
  # compute estimators using bootstrap data set
  if(ml){
    
    if(any(method %in% c("aipw", "tmle", "sens")) & return_se == FALSE){
      boot_ml_models <- vegrowth::fit_ml_models(data = boot_data, 
                                                 estimand = estimand,
                                                 method = method, 
                                                 exclusion_restriction = exclusion_restriction,
                                                 cross_world = cross_world,
                                                 Y_name = Y_name,
                                                 Z_name = Z_name,
                                                 S_name = S_name,
                                                 X_name = X_name,
                                                 Y_Z_X_library = Y_Z_X_library,
                                                 Y_X_library = Y_X_library,
                                                 S_X_library = S_X_library,
                                                 S_Z_X_library = S_Z_X_library,
                                                 Z_X_library = Z_X_library,
                                                 family = family,
                                                 v_folds = v_folds)
    } 
    
    if(any(method %in% c("gcomp", "ipw"))){
      boot_models <- vegrowth::fit_models(data = boot_data, 
                                          estimand = estimand,
                                          method = method, 
                                          exclusion_restriction = exclusion_restriction,
                                          cross_world = cross_world,
                                          Y_name = Y_name,
                                          Z_name = Z_name,
                                          S_name = S_name,
                                          X_name = X_name,
                                          Y_Z_X_model = Y_Z_X_model,
                                          Y_X_S1_model = Y_X_S1_model,
                                          Y_X_S0_model = Y_X_S0_model,
                                          S_X_model = S_X_model,
                                          S_Z_X_model = S_Z_X_model,
                                          Z_X_model = Z_X_model,
                                          family = family)
    }
    
  } else{
    # GLMS for all
    boot_models <- vegrowth::fit_models(data = boot_data, 
                                        estimand = estimand,
                                        method = method, 
                                        exclusion_restriction = exclusion_restriction,
                                        cross_world = cross_world,
                                        Y_name = Y_name,
                                        Z_name = Z_name,
                                        S_name = S_name,
                                        X_name = X_name,
                                        Y_Z_X_model = Y_Z_X_model,
                                        Y_X_S1_model = Y_X_S1_model,
                                        Y_X_S0_model = Y_X_S0_model,
                                        S_X_model = S_X_model,
                                        S_Z_X_model = S_Z_X_model,
                                        Z_X_model = Z_X_model,
                                        family = family)
  } 
  
  out <- vector("list", length = length(estimand))
  names(out) <- estimand
  
  # Naturally infected --------------------------------------------------------
  
  if("nat_inf" %in% estimand){
    
    for(er in exclusion_restriction){
      
      # Exclustion restriction -- save as <estimator>_ER 
      er_suffix <- if (er) "_ER" else ""
      
      if("ipw" %in% method){
        estimator <- paste0("ipw", er_suffix)
        out$nat_inf[[estimator]] <- do_ipw_nat_inf(data = boot_data, models = boot_models, S_name = S_name, Y_name = Y_name, Z_name = Z_name, exclusion_restriction = er)
      }
      
      # Cross-world assumption can be toggled for AIPW only at this point
      for(cw in cross_world){
        
        cw_suffix <- if(cw) "_CW" else ""
        
        # cannot have scenario where both are false
        if(er == FALSE & cw == FALSE) next
        
        if("gcomp" %in% method){
          estimator <- paste0("gcomp", er_suffix, cw_suffix)
          out$nat_inf[[estimator]] <- do_gcomp_nat_inf(data = boot_data, models = boot_models, Z_name = Z_name, X_name = X_name, exclusion_restriction = er, cross_world = cw, two_part_model = two_part_model)
        }
        
        # if we want bootstrap SE for AIPW (otherwise return closed form SE when we get point estimate)
        if("aipw" %in% method & return_se == FALSE){
          estimator <- paste0("aipw", er_suffix, cw_suffix)
          if(ml){
            out$nat_inf[[estimator]] <- do_aipw_nat_inf(data = boot_data, models = boot_ml_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, X_name = X_name, return_se = return_se, exclusion_restriction = er, cross_world = cw, two_part_model = two_part_model)
          } else{
            out$nat_inf[[estimator]] <- do_aipw_nat_inf(data = boot_data, models = boot_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, X_name = X_name, return_se = return_se, exclusion_restriction = er, cross_world = cw, two_part_model = two_part_model)
          }
        }
      }
      
      
    }
    
    if("tmle" %in% method & return_se == FALSE){
      if(ml){
        out$nat_inf$tmle <- do_tmle_nat_inf(data = boot_data, models = boot_ml_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, return_se = return_se)
      } else{
        out$nat_inf$tmle <- do_tmle_nat_inf(data = boot_data, models = boot_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, return_se = return_se)
      }
    }
    
    if("sens" %in% method & return_se == FALSE){
      if(ml){
        out$nat_inf$sens<- do_sens_aipw_nat_inf(data = boot_data, models = boot_ml_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, epsilon = epsilon, return_se = return_se)
      } else{
        out$nat_inf$sens <- do_sens_aipw_nat_inf(data = boot_data, models = boot_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, epsilon = epsilon, return_se = return_se)
      }
    }
    
    if("bound" %in% method){
      out$nat_inf$bound <- get_bound_nat_inf(data = boot_data, Y_name = Y_name, Z_name = Z_name, S_name = S_name, family = family)
    }
    
  }
  
  # Doomed --------------------------------------------------------------------
  
  if("doomed" %in% estimand){
    
    if("gcomp" %in% method){
      out$doomed$gcomp <- do_gcomp_doomed(data = boot_data, models = boot_models)
    }
    
    if("ipw" %in% method){
      out$doomed$ipw <- do_ipw_doomed(data = boot_data, models = boot_models, S_name = S_name, Y_name = Y_name, Z_name = Z_name)
    }
    
    if("aipw" %in% method & return_se == FALSE){
      if(ml){
        out$doomed$aipw <- do_aipw_doomed(data = boot_data, models = boot_ml_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, return_se = return_se)
      } else{
        out$doomed$aipw <- do_aipw_doomed(data = boot_data, models = boot_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, return_se = return_se)
      }
    }
    
    if("bound" %in% method){
      out$doomed$bound <- get_bound_doomed(data = boot_data, Y_name = Y_name, Z_name = Z_name, S_name = S_name, family = family)
    }
    
  }
  
  # Population ----------------------------------------------------------------
  
  if("pop" %in% estimand){
    
    if("gcomp" %in% method){
      out$pop$gcomp <- do_gcomp_pop(data = boot_data, models = boot_models, Z_name = Z_name, X_name = X_name)
    }
    
    if("ipw" %in% method){
      out$pop$ipw <- do_ipw_pop(data = boot_data, models = boot_models, Y_name = Y_name, Z_name = Z_name)
    }
    
    if("aipw" %in% method & return_se == FALSE){
      if(ml){
        out$pop$aipw <- do_aipw_pop(data = boot_data, models = boot_ml_models, Z_name = Z_name, Y_name = Y_name, X_name = X_name, return_se = return_se, two_part_model = two_part_model)
      } else{
        out$pop$aipw <- do_aipw_pop(data = boot_data, models = boot_models, Z_name = Z_name, Y_name = Y_name, X_name = X_name, return_se = return_se, two_part_model = two_part_model)
      }
    }
    
  }

  return(out)

}

#' Function to replicate n_boot bootstrap samples and get bootstrap standard error
#' 
#' @param data original data to bootstrap
#' @param Y_name growth outcome variable name
#' @param Z_name vaccination variable name
#' @param X_name covariate name(s)
#' @param S_name infection variable name
#' @param n_boot number of bootstrap replicates
#' @param est character vector of names of estimators to use for growth effect
#' @param ml boolean to use SuperLearner models, default FALSE
#' @param Y_Z_X_model optional specify model to be used for fitting growth on vaccine + covariates, otherwise growth on all covariates
#' @param Y_X_S1_model optional specify model to be used for fitting growth on covariates in the infected, otherwise growth on all covariates
#' @param Y_X_S0_model optional specify model to be used for fitting growth on covariates in the uninfected, otherwise growth on all covariates
#' @param S_X_model optional specify model to be used for fitting infection on covariates, otherwise infection on all covariates
#' @param Y_Z_X_library optional specify SuperLearner libraries for model fitting growth on covariates + vaccine, default glm
#' @param Y_X_library optional specify SuperLearner libraries for model fitting growth on covariates, default glm
#' @param S_X_library optional specify SuperLearner libraries for model fitting infection on covariates, default glm
#' @param family family for outcome model, defaults to gaussian for growth
#' @param v_folds number of cross validation folds for SuperLearner, default 3
#' @param effect_dir direction of beneficial effect, defaults to "positive" for beneficial outcome. Used for one-side tests of bounds.  
#' @param epsilon a vector of values for the sensitivity parameter
#' 
#' @returns list containing bootstrap se and 95% CI bounds for estimators specified in est
bootstrap_estimates <- function(
    data, 
    Y_name = "Y",
    Z_name = "Z",
    X_name = "X",
    S_name = "S", 
    n_boot = 1000, 
    estimand = c("nat_inf", "doomed", "pop"),
    method = c("gcomp", "ipw", "aipw", "tmle", "bound", "sens"),
    exclusion_restriction = FALSE,
    cross_world = TRUE,
    ml = ml,
    Y_Z_X_model = NULL,
    Y_X_S1_model = NULL, 
    Y_X_S0_model = NULL, 
    S_X_model = NULL,
    S_Z_X_model = NULL,
    Z_X_model = NULL,
    Y_Z_X_library = c("SL.glm"),
    Y_X_library = c("SL.glm"),
    S_X_library = c("SL.glm"),
    S_Z_X_library = c("SL.glm"),
    Z_X_library = c("SL.mean"),
    family = "gaussian",
    v_folds = 3,
    effect_dir = "positive",
    epsilon = exp(seq(log(0.5), log(2), length = 50)),
    return_se = TRUE,
    two_part_model = FALSE
){
  
  # Initial boot_estimates for all viable estimand & method combinations
  boot_estimates <- replicate(n_boot, one_boot(data, 
                                               Y_name = Y_name,
                                               Z_name = Z_name,
                                               S_name = S_name,
                                               X_name = X_name,
                                               estimand = estimand,
                                               method = method,
                                               exclusion_restriction = exclusion_restriction,
                                               cross_world = cross_world,
                                               ml = ml,
                                               Y_Z_X_model = Y_Z_X_model,
                                               Y_X_S1_model = Y_X_S1_model,
                                               Y_X_S0_model = Y_X_S0_model,
                                               S_X_model = S_X_model,
                                               S_Z_X_model = S_Z_X_model,
                                               Z_X_model = Z_X_model,
                                               Y_Z_X_library = Y_Z_X_library,
                                               Y_X_library = Y_X_library,
                                               S_X_library = S_X_library,
                                               S_Z_X_library = S_Z_X_library,
                                               Z_X_library = Z_X_library,
                                               v_folds = v_folds,
                                               family = family,
                                               epsilon = epsilon,
                                               return_se = return_se,
                                               two_part_model = two_part_model), simplify = FALSE)
    # List to store results
    out <- vector("list", length = length(estimand))
    names(out) <- estimand
  
    # Naturally infected --------------------------------------------------------
    
    if("nat_inf" %in% estimand){
      
      for(er in exclusion_restriction){
        
        # Exclustion restriction -- save as <estimator>_ER 
        er_suffix <- if (er) "_ER" else ""
        
        if("ipw" %in% method){
          estimator <- paste0("ipw", er_suffix)
          out$nat_inf[[estimator]]$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "nat_inf", method = estimator)
        }
        
        for(cw in cross_world){
          
          # can't have both false
          if(er == FALSE & cw == FALSE) next
          
          cw_suffix <- if (er) "_CW" else ""
          
          if("gcomp" %in% method){
            estimator <- paste0("gcomp", er_suffix, cw_suffix)
            out$nat_inf[[estimator]]$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "nat_inf", method = estimator)
          }
          
          # if we want bootstrap SE for AIPW (otherwise return closed form SE when we get point estimate)
          if("aipw" %in% method & return_se == FALSE){
            estimator <- paste0("aipw", er_suffix, cw_suffix)
            out$nat_inf[[estimator]]$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "nat_inf", method = estimator)
          }
          
        }
        
      }

      if("tmle" %in% method & return_se == FALSE){
        out$nat_inf$tmle$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "nat_inf", method = "tmle")
      }
      
      if("sens" %in% method & return_se == FALSE){
        out$nat_inf$sens$boot_se <- get_boot_se_sens(boot_estimates = boot_estimates, estimand = "nat_inf", method = "sens")
      }
      
      if("bound" %in% method){
        out$nat_inf$bound$boot_se <- get_boot_se_bound(boot_estimates = boot_estimates, estimand = "nat_inf", method = "bound")
      }
      
    }
    
    
    # Doomed --------------------------------------------------------------------
    
    if("doomed" %in% estimand){
      
      if("gcomp" %in% method){
        out$doomed$gcomp$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "doomed", method = "gcomp")
      }
      
      if("ipw" %in% method){
        out$doomed$ipw$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "doomed", method = "ipw")
      }
      
      if("aipw" %in% method & return_se == FALSE){
        out$doomed$aipw$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "doomed", method = "aipw")
      }
      
      if("bound" %in% method){
        out$doomed$bound$boot_se <- get_boot_se_bound(boot_estimates = boot_estimates, estimand = "doomed", method = "bound")
      }
      
    }
    
    # Population ----------------------------------------------------------------
    
    if("pop" %in% estimand){
      
      if("gcomp" %in% method){
        out$pop$gcomp$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "pop", method = "gcomp")
      }
      
      if("ipw" %in% method){
        out$pop$ipw$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "pop", method = "ipw")
      }
      
      if("aipw" %in% method & return_se == FALSE){
        out$pop$aipw$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "pop", method = "aipw")
      }
      
    }

  return(out)
    
}