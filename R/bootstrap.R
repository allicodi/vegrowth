#' Function for one bootstrap sample
#'
#' @param data original data to bootstrap
#' @param G_name growth outcome variable name
#' @param V_name vaccination variable name
#' @param X_name covariate name(s)
#' @param Y_name infection variable name
#' @param est character vector of names of estimators to use for growth effect
#' @param ml boolean to use SuperLearner models, default FALSE
#' @param G_X_Y1_model optional specify model to be used for fitting growth on covariates in the infected, otherwise growth on all covariates
#' @param G_X_Y0_model optional specify model to be used for fitting growth on covariates in the uninfected, otherwise growth on all covariates
#' @param Y_X_model optional specify model to be used for fitting infection on covariates, otherwise infection on all covariates
#' @param G_X_library optional specify SuperLearner libraries for model fitting growth on covariates, default glm
#' @param Y_X_library optional specify SuperLearner libraries for model fitting infection on covariates, default glm
#' @param family family for outcome model, defaults to gaussian for growth
#' @param v_folds number of cross validation folds for SuperLearner, default 3
#' @param effect_dir direction of beneficial effect, defaults to "positive" for beneficial outcome. Used for one-side tests of bounds.  
#' 
#' @returns list containing results for specified estimators on single bootstrap sample
one_boot <- function(
    data,
    G_name = "G",
    V_name = "V",
    X_name = "X",
    Y_name = "Y", 
    est = c("gcomp_pop_estimand", "gcomp", "efficient_aipw", "efficient_tmle", "hudgens_adj_upper", "hudgens_adj_lower"),
    ml = FALSE, 
    G_V_X_model = NULL,
    G_X_Y1_model = NULL,
    G_X_Y0_model = NULL,
    Y_X_model = NULL,
    G_V_X_library = c("SL.glm"),
    G_X_library = c("SL.glm"),
    Y_X_library = c("SL.glm"),
    family = "gaussian",
    v_folds = 3,
    effect_dir = "positive"
){
  n <- dim(data)[1]
  boot_row_idx <- sample(1:n, replace=TRUE)
  boot_data <- data[boot_row_idx,]
  
  # compute estimators using bootstrap data set
  if(ml){
    
    if(any(est %in% c("efficient_aipw", "efficient_tmle"))){
      boot_ml_models <- vegrowth::fit_ml_models(data = boot_data, 
                                           est = est, 
                                           G_name = G_name,
                                           V_name = V_name,
                                           Y_name = Y_name,
                                           X_name = X_name,
                                           G_V_X_library = G_V_X_library,
                                           G_X_library = G_X_library,
                                           Y_X_library = Y_X_library,
                                           family = family,
                                           v_folds = v_folds)
    } 
    
    if(any(est %in% c("gcomp_pop_estimand", "gcomp",
                      "hudgens_adj_lower", "hudgens_adj_upper", 
                      "hudgens_lower", "hudgens_upper",
                      "hudgens_upper_doomed", "hudgens_lower_doomed"))){
      boot_models <- vegrowth::fit_models(data = boot_data, 
                                     est = est, 
                                     G_name = G_name,
                                     V_name = V_name,
                                     Y_name = Y_name,
                                     X_name = X_name,
                                     G_V_X_model = G_V_X_model,
                                     G_X_Y1_model = G_X_Y1_model,
                                     G_X_Y0_model = G_X_Y0_model,
                                     Y_X_model = Y_X_model,
                                     family = family)
    }
    
  } else{
    boot_models <- vegrowth::fit_models(data = boot_data, 
                                   est = est, 
                                   G_name = G_name,
                                   V_name = V_name,
                                   Y_name = Y_name,
                                   X_name = X_name,
                                   G_V_X_model = G_V_X_model,
                                   G_X_Y1_model = G_X_Y1_model,
                                   G_X_Y0_model = G_X_Y0_model,
                                   Y_X_model = Y_X_model,
                                   family = family)
  } 
  
  out <- list()
  
  if("gcomp" %in% est){
    out$growth_effect <- do_gcomp(boot_data, boot_models)
  }
  if("gcomp_pop_estimand" %in% est){
    out$growth_effect_pop <- do_gcomp_pop_estimand(boot_data, 
                                                   boot_models, 
                                                   V_name = V_name,
                                                   X_name = X_name)
  }
  if("efficient_aipw" %in% est){
    if(ml){
      out$growth_effect_aipw <- do_efficient_aipw(boot_data, 
                                                  boot_ml_models,
                                                  G_name = G_name,
                                                  V_name = V_name,
                                                  Y_name = Y_name)
    } else{
      out$growth_effect_aipw <- do_efficient_aipw(boot_data, 
                                                  boot_models,
                                                  G_name = G_name,
                                                  V_name = V_name,
                                                  Y_name = Y_name)
    }
  }
  if("efficient_tmle" %in% est){
    if(ml){
      out$growth_effect_tmle <- do_efficient_tmle(boot_data, 
                                                  boot_ml_models,
                                                  G_name = G_name,
                                                  V_name = V_name,
                                                  Y_name = Y_name)
    } else{
      out$growth_effect_tmle <- do_efficient_tmle(boot_data, 
                                                  boot_models,
                                                  G_name = G_name,
                                                  V_name = V_name,
                                                  Y_name = Y_name)
    }
  } 
  if("hudgens_adj_upper" %in% est){
    out$growth_effect_hudgens_adj_upper <- get_adjusted_hudgens_stat(boot_data,
                                                                     boot_models,
                                                                     family = family, 
                                                                     lower_bound = FALSE)
  } 
  if("hudgens_adj_lower" %in% est){
    out$growth_effect_hudgens_adj_lower <- get_adjusted_hudgens_stat(boot_data,
                                                                     boot_models,
                                                                     family = family, 
                                                                     lower_bound = TRUE)
  } 
  if("hudgens_lower" %in% est){
    out$growth_effect_hudgens_lower <- get_hudgens_stat(boot_data, G_name = G_name, V_name = V_name, Y_name = Y_name, lower_bound = TRUE)
  }
  if("hudgens_upper" %in% est){
    out$growth_effect_hudgens_upper <- get_hudgens_stat(boot_data, G_name = G_name, V_name = V_name, Y_name = Y_name, lower_bound = FALSE)
  }
  if("hudgens_lower_doomed" %in% est){
    out$growth_effect_hudgens_lower_doomed <- get_hudgens_stat_doomed(boot_data, G_name = G_name, V_name = V_name, Y_name = Y_name, lower_bound = TRUE)
  }
  if("hudgens_upper_doomed" %in% est){
    out$growth_effect_hudgens_upper_doomed <- get_hudgens_stat_doomed(boot_data, G_name = G_name, V_name = V_name, Y_name = Y_name, lower_bound = FALSE)
  }
  
  return(out)
}

#' Function to replicate n_boot bootstrap samples and get bootstrap standard error
#' 
#' @param data original data to bootstrap
#' @param G_name growth outcome variable name
#' @param V_name vaccination variable name
#' @param X_name covariate name(s)
#' @param Y_name infection variable name
#' @param n_boot number of bootstrap replicates
#' @param est character vector of names of estimators to use for growth effect
#' @param ml boolean to use SuperLearner models, default FALSE
#' @param G_V_X_model optional specify model to be used for fitting growth on vaccine + covariates, otherwise growth on all covariates
#' @param G_X_Y1_model optional specify model to be used for fitting growth on covariates in the infected, otherwise growth on all covariates
#' @param G_X_Y0_model optional specify model to be used for fitting growth on covariates in the uninfected, otherwise growth on all covariates
#' @param Y_X_model optional specify model to be used for fitting infection on covariates, otherwise infection on all covariates
#' @param G_V_X_library optional specify SuperLearner libraries for model fitting growth on covariates + vaccine, default glm
#' @param G_X_library optional specify SuperLearner libraries for model fitting growth on covariates, default glm
#' @param Y_X_library optional specify SuperLearner libraries for model fitting infection on covariates, default glm
#' @param family family for outcome model, defaults to gaussian for growth
#' @param v_folds number of cross validation folds for SuperLearner, default 3
#' @param effect_dir direction of beneficial effect, defaults to "positive" for beneficial outcome. Used for one-side tests of bounds.  
#' 
#' @returns list containing bootstrap se and 95% CI bounds for estimators specified in est
bootstrap_estimates <- function(
    data, 
    G_name = "G",
    V_name = "V",
    X_name = "X",
    Y_name = "Y", 
    n_boot = 1000, 
    est = c("gcomp_pop_estimand", "gcomp", "efficient_aipw", "efficient_tmle", "hudgens_adj_upper", "hudgens_adj_lower"),
    ml = ml,
    G_V_X_model = NULL,
    G_X_Y1_model = NULL, 
    G_X_Y0_model = NULL, 
    Y_X_model = NULL,
    G_V_X_library = c("SL.glm"),
    G_X_library = c("SL.glm"),
    Y_X_library = c("SL.glm"),
    family = "gaussian",
    v_folds = 3,
    effect_dir = "positive"
){
  
  boot_estimates <- replicate(n_boot, one_boot(data, 
                                               G_name = G_name,
                                               V_name = V_name,
                                               Y_name = Y_name,
                                               X_name = X_name,
                                               est = est,
                                               ml = ml,
                                               G_V_X_model = G_V_X_model,
                                               G_X_Y1_model = G_X_Y1_model,
                                               G_X_Y0_model = G_X_Y0_model,
                                               Y_X_model = Y_X_model,
                                               G_V_X_library = G_V_X_library,
                                               G_X_library = G_X_library,
                                               Y_X_library = Y_X_library,
                                               v_folds = v_folds,
                                               family = family))
  
  out <- list()
  
  if("gcomp_pop_estimand" %in% est){
    
    # if there was only one method, unlist boot_est as is
    if(length(boot_estimates) == n_boot){
      #growth_effect_pop <- unlist(boot_estimates)
      growth_effect_pop <- data.frame(do.call(rbind, boot_estimates))
    } else{
      #growth_effect_pop <- unlist(boot_estimates["growth_effect_pop",])
      growth_effect_pop <- data.frame(do.call(rbind, boot_estimates["growth_effect_pop",]))
    }
    
    # Additive
    ci_gcomp_pop_estimand_additive <- quantile(growth_effect_pop$additive_effect, p = c(0.025, 0.975), na.rm=TRUE)
    out$se_gcomp_pop_estimand_additive <- sd(growth_effect_pop$additive_effect, na.rm = TRUE)
    out$lower_ci_gcomp_pop_estimand_additive <- ci_gcomp_pop_estimand_additive[1]
    out$upper_ci_gcomp_pop_estimand_additive <- ci_gcomp_pop_estimand_additive[2]
    
    # Multiplicative 
    ci_gcomp_pop_log_mult <- quantile(growth_effect_pop$log_multiplicative_effect, p = c(0.025, 0.975), na.rm=TRUE)
    out$se_gcomp_pop_estimand_log_mult <- sd(growth_effect_pop$log_multiplicative_effect, na.rm = TRUE)
    out$lower_ci_gcomp_pop_estimand_mult <- exp(ci_gcomp_pop_log_mult[1])
    out$upper_ci_gcomp_pop_estimand_mult <- exp(ci_gcomp_pop_log_mult[2])
    
    # ci_gcomp_pop_estimand <- quantile(growth_effect_pop, p = c(0.025, 0.975), na.rm=TRUE)
    # out$se_gcomp_pop_estimand <- sd(growth_effect_pop, na.rm = TRUE)
    # out$lower_ci_gcomp_pop_estimand <- ci_gcomp_pop_estimand[1]
    # out$upper_ci_gcomp_pop_estimand <- ci_gcomp_pop_estimand[2]
  }
  if("gcomp" %in% est){
    
    if(length(boot_estimates) == n_boot){
      #growth_effect <- unlist(boot_estimates)
      growth_effect <- data.frame(do.call(rbind, boot_estimates))
    } else{
      #growth_effect <- unlist(boot_estimates["growth_effect",])
      growth_effect <- data.frame(do.call(rbind, boot_estimates["growth_effect",]))
    }
    
    # Additive
    ci_gcomp_additive <- quantile(growth_effect$additive_effect, p = c(0.025, 0.975), na.rm=TRUE)
    out$se_gcomp_additive <- sd(growth_effect$additive_effect, na.rm = TRUE)
    out$lower_ci_gcomp_additive <- ci_gcomp_additive[1]
    out$upper_ci_gcomp_additive <- ci_gcomp_additive[2]
    
    # Multiplicative 
    ci_gcomp_log_mult <- quantile(growth_effect$log_multiplicative_effect, p = c(0.025, 0.975), na.rm=TRUE)
    out$se_gcomp_log_mult <- sd(growth_effect$log_multiplicative_effect, na.rm = TRUE)
    out$lower_ci_gcomp_mult <- exp(ci_gcomp_log_mult[1])
    out$upper_ci_gcomp_mult <- exp(ci_gcomp_log_mult[2])
    
  }
  if("efficient_aipw" %in% est){
    
    if(length(boot_estimates) == n_boot){
      #growth_effect_aipw <- unlist(boot_estimates)
      growth_effect_aipw <- data.frame(do.call(rbind, boot_estimates))
    } else{
      #growth_effect_aipw <- unlist(boot_estimates["growth_effect_aipw",])
      growth_effect_aipw <- data.frame(do.call(rbind, boot_estimates["growth_effect_aipw",]))
    }
    
    # Additive
    ci_efficient_aipw_additive <- quantile(growth_effect_aipw$additive_effect, p = c(0.025, 0.975), na.rm=TRUE)
    out$se_efficient_aipw_additive <- sd(growth_effect_aipw$additive_effect, na.rm = TRUE)
    out$lower_ci_efficient_aipw_additive <- ci_efficient_aipw_additive[1]    
    out$upper_ci_efficient_aipw_additive <- ci_efficient_aipw_additive[2]
    
    # Multiplicative
    ci_efficient_aipw_log_mult <- quantile(growth_effect_aipw$log_multiplicative_effect, p = c(0.025, 0.975), na.rm=TRUE)
    out$se_efficient_aipw_log_mult <- sd(growth_effect_aipw$log_multiplicative_effect, na.rm = TRUE)
    out$lower_ci_efficient_aipw_mult <- exp(ci_efficient_aipw_log_mult[1])
    out$upper_ci_efficient_aipw_mult <- exp(ci_efficient_aipw_log_mult[2])
    
  }
  if("efficient_tmle" %in% est){
    
    if(length(boot_estimates) == n_boot){
      #growth_effect_tmle <- unlist(boot_estimates)
      growth_effect_tmle <- data.frame(do.call(rbind, boot_estimates))
    } else{
      #growth_effect_tmle <- unlist(boot_estimates["growth_effect_tmle",])
      growth_effect_tmle <- data.frame(do.call(rbind, boot_estimates["growth_effect_tmle",]))
    }
    
    # Additive
    ci_efficient_tmle_additive <- quantile(growth_effect_tmle$additive_effect, p = c(0.025, 0.975), na.rm=TRUE)
    out$se_efficient_tmle_additive <- sd(growth_effect_tmle$additive_effect, na.rm = TRUE)
    out$lower_ci_efficient_tmle_additive <- ci_efficient_tmle_additive[1]    
    out$upper_ci_efficient_tmle_additive <- ci_efficient_tmle_additive[2]
    
    # Multiplicative
    ci_efficient_tmle_log_mult <- quantile(growth_effect_tmle$log_multiplicative_effect, p = c(0.025, 0.975), na.rm=TRUE)
    out$se_efficient_tmle_log_mult <- sd(growth_effect_tmle$log_multiplicative_effect, na.rm = TRUE)
    out$lower_ci_efficient_tmle <- exp(ci_efficient_tmle_log_mult[1])    
    out$upper_ci_efficient_tmle <- exp(ci_efficient_tmle_log_mult[2])
    
  }
  if("hudgens_adj_upper" %in% est){
    
    # if there was only one method, unlist boot_est as is
    if(length(boot_estimates) == n_boot){
      growth_effect_hudgens_adj_upper <- unlist(boot_estimates)
    } else{
      growth_effect_hudgens_adj_upper <- unlist(boot_estimates["growth_effect_hudgens_adj_upper",])
    }
    
    ci_hudgens_adj_upper <- quantile(growth_effect_hudgens_adj_upper, p = c(0.025, 0.975), na.rm=TRUE)
    out$se_hudgens_adj_upper <- sd(growth_effect_hudgens_adj_upper, na.rm = TRUE)
    out$lower_ci_hudgens_adj_upper <- ci_hudgens_adj_upper[1]
    out$upper_ci_hudgens_adj_upper <- ci_hudgens_adj_upper[2]
    
  } 
  if("hudgens_adj_lower" %in% est){
    
    # if there was only one method, unlist boot_est as is
    if(length(boot_estimates) == n_boot){
      growth_effect_hudgens_adj_lower <- unlist(boot_estimates)
    } else{
      growth_effect_hudgens_adj_lower <- unlist(boot_estimates["growth_effect_hudgens_adj_lower",])
    }
    
    ci_hudgens_adj_lower <- quantile(growth_effect_hudgens_adj_lower, p = c(0.025, 0.975), na.rm=TRUE)
    out$se_hudgens_adj_lower <- sd(growth_effect_hudgens_adj_lower, na.rm = TRUE)
    out$lower_ci_hudgens_adj_lower <- ci_hudgens_adj_lower[1]
    out$upper_ci_hudgens_adj_lower <- ci_hudgens_adj_lower[2]
    
  }
  if("hudgens_lower" %in% est){
    
    # if there was only one method, unlist boot_est as is
    if(length(boot_estimates) == n_boot){
      #growth_effect_hudgens_lower <- unlist(boot_estimates)
      growth_effect_hudgens_lower <- data.frame(do.call(rbind, boot_estimates))
    } else{
      #growth_effect_hudgens_lower <- unlist(boot_estimates["growth_effect_hudgens_lower",])
      growth_effect_hudgens_lower <- do.call(rbind, boot_estimates["growth_effect_hudgens_lower",])
    }
    
    ci_hudgens_lower <- quantile(growth_effect_hudgens_lower, p = c(0.025, 0.975), na.rm=TRUE)
    out$se_hudgens_lower <- sd(growth_effect_hudgens_lower, na.rm = TRUE)
    out$lower_ci_hudgens_lower <- ci_hudgens_lower[1]
    out$upper_ci_hudgens_lower <- ci_hudgens_lower[2]
    
  }
  if("hudgens_upper" %in% est){
    
    # if there was only one method, unlist boot_est as is
    if(length(boot_estimates) == n_boot){
      #growth_effect_hudgens_upper <- unlist(boot_estimates)
      growth_effect_hudgens_upper <- data.frame(do.call(rbind, boot_estimates))
    } else{
      #growth_effect_hudgens_upper <- unlist(boot_estimates["growth_effect_hudgens_upper",])
      growth_effect_hudgens_upper <- do.call(rbind, boot_estimates["growth_effect_hudgens_upper",])
    }
    
    ci_hudgens_upper <- quantile(growth_effect_hudgens_upper, p = c(0.025, 0.975), na.rm=TRUE)
    out$se_hudgens_upper <- sd(growth_effect_hudgens_upper, na.rm = TRUE)
    out$lower_ci_hudgens_upper <- ci_hudgens_upper[1]
    out$upper_ci_hudgens_upper <- ci_hudgens_upper[2]
    
  }
  if("hudgens_lower_doomed" %in% est){
    
    # if there was only one method, unlist boot_est as is
    if(length(boot_estimates) == n_boot){
      #growth_effect_hudgens_lower_doomed <- unlist(boot_estimates)
      growth_effect_hudgens_lower_dooemd <- do.call(rbind, boot_estimates)
    } else{
      #growth_effect_hudgens_lower_doomed <- unlist(boot_estimates["growth_effect_hudgens_lower_doomed",])
      growth_effect_hudgens_lower_doomed <- do.call(rbind, boot_estimates["growth_effect_hudgens_lower_doomed",])
    }
    
    ci_hudgens_lower_doomed <- quantile(growth_effect_hudgens_lower_doomed, p = c(0.025, 0.975), na.rm=TRUE)
    out$se_hudgens_lower_doomed <- sd(growth_effect_hudgens_lower_doomed, na.rm = TRUE)
    out$lower_ci_hudgens_lower_doomed <- ci_hudgens_lower_doomed[1]
    out$upper_ci_hudgens_lower_doomed <- ci_hudgens_lower_doomed[2]
    
  }
  if("hudgens_upper_doomed" %in% est){
    
    # if there was only one method, unlist boot_est as is
    if(length(boot_estimates) == n_boot){
      #growth_effect_hudgens_upper_doomed <- unlist(boot_estimates)
      growth_effect_hudgens_upper_dooemd <- do.call(rbind, boot_estimates)
    } else{
      #growth_effect_hudgens_upper_doomed <- unlist(boot_estimates["growth_effect_hudgens_upper_doomed",])
      growth_effect_hudgens_upper_doomed <- do.call(rbind, boot_estimates["growth_effect_hudgens_upper_doomed",])
    }
    
    ci_hudgens_upper_doomed <- quantile(growth_effect_hudgens_upper_doomed, p = c(0.025, 0.975), na.rm=TRUE)
    out$se_hudgens_upper_doomed <- sd(growth_effect_hudgens_upper_doomed, na.rm = TRUE)
    out$lower_ci_hudgens_upper_doomed <- ci_hudgens_upper_doomed[1]
    out$upper_ci_hudgens_upper_doomed <- ci_hudgens_upper_doomed[2]
    
  }
  
  return(out)
}
