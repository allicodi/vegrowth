#' Function for one bootstrap sample
#'
#' @param data original data to bootstrap
#' @param G_name growth outcome variable name
#' @param V_name vaccination variable name
#' @param X_name covariate name(s)
#' @param Y_name infection variable name
#' @param est character vector of names of estimators to use for growth effect
#' @param G_X_model optional specify model to be used for fitting growth on covariates, otherwise growth on all covariates
#' @param Y_X_model optional specify model to be used for fitting infection on covariates, otherwise infection on all covariates
#' @param family family for outcome model, defaults to gaussian for growth
#' 
#' @returns list containing results for specified estimators on single bootstrap sample
one_boot <- function(
    data,
    G_name = "G",
    V_name = "V",
    X_name = "X",
    Y_name = "Y", 
    est = c("gcomp_pop_estimand", "gcomp", "efficient_aipw", "efficient_tmle", "hudgens_adj_upper", "hudgens_adj_lower"),
    G_X_model = NULL,
    Y_X_model = NULL,
    family = "gaussian"
){
  n <- dim(data)[1]
  boot_row_idx <- sample(1:n, replace=TRUE)
  boot_data <- data[boot_row_idx,]
  
  # compute estimators using bootstrap data set
  boot_models <- vegrowth::fit_models(boot_data, 
                            G_name = G_name,
                            V_name = V_name,
                            X_name = X_name,
                            Y_name = Y_name,
                            est = est, 
                            G_X_model = G_X_model,
                            Y_X_model = Y_X_model,
                            family = family)
  
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
    out$growth_effect_aipw <- do_efficient_aipw(boot_data, 
                                                boot_models,
                                                G_name = G_name,
                                                V_name = V_name,
                                                Y_name = Y_name)
  }
  if("efficient_tmle" %in% est){
    out$growth_effect_tmle <- do_efficient_tmle(boot_data, 
                                                boot_models,
                                                G_name = G_name,
                                                V_name = V_name,
                                                Y_name = Y_name)
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
#' @param G_X_model optional specify model to be used for fitting growth on covariates, otherwise growth on all covariates
#' @param Y_X_model optional specify model to be used for fitting infection on covariates, otherwise infection on all covariates
#' @param family family for outcome model, defaults to gaussian for growth
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
    G_X_model = NULL, 
    Y_X_model = NULL,
    family = "gaussian"
){
  
  boot_estimates <- replicate(n_boot, one_boot(data, 
                                               G_name = G_name,
                                               V_name = V_name,
                                               Y_name = Y_name,
                                               X_name = X_name,
                                               est = est,
                                               G_X_model = G_X_model, 
                                               Y_X_model = Y_X_model,
                                               family = family))
  
  out <- list()
  
  if("gcomp_pop_estimand" %in% est){
    
    # if there was only one method, unlist boot_est as is
    if(length(boot_estimates) == n_boot){
      growth_effect_pop <- unlist(boot_estimates)
    } else{
      growth_effect_pop <- unlist(boot_estimates["growth_effect_pop",])
    }
    
    ci_gcomp_pop_estimand <- quantile(growth_effect_pop, p = c(0.025, 0.975), na.rm=TRUE)
    out$se_gcomp_pop_estimand <- sd(growth_effect_pop, na.rm = TRUE)
    out$lower_ci_gcomp_pop_estimand <- ci_gcomp_pop_estimand[1]
    out$upper_ci_gcomp_pop_estimand <- ci_gcomp_pop_estimand[2]
  }
  if("gcomp" %in% est){
    
    if(length(boot_estimates) == n_boot){
      growth_effect <- unlist(boot_estimates)
    } else{
      growth_effect <- unlist(boot_estimates["growth_effect",])
    }
    
    ci_gcomp <- quantile(growth_effect, p = c(0.025, 0.975), na.rm=TRUE)
    out$se_gcomp <- sd(growth_effect, na.rm = TRUE)
    out$lower_ci_gcomp <- ci_gcomp[1]
    out$upper_ci_gcomp <- ci_gcomp[2]
  }
  if("efficient_aipw" %in% est){
    
    if(length(boot_estimates) == n_boot){
      growth_effect_aipw <- unlist(boot_estimates)
    } else{
      growth_effect_aipw <- unlist(boot_estimates["growth_effect_aipw",])
    }
    
    ci_efficient_aipw <- quantile(growth_effect_aipw, p = c(0.025, 0.975), na.rm=TRUE)
    out$se_efficient_aipw <- sd(growth_effect_aipw, na.rm = TRUE)
    out$lower_ci_efficient_aipw <- ci_efficient_aipw[1]    
    out$upper_ci_efficient_aipw <- ci_efficient_aipw[2]
  }
  if("efficient_tmle" %in% est){
    
    if(length(boot_estimates) == n_boot){
      growth_effect_tmle <- unlist(boot_estimates)
    } else{
      growth_effect_tmle <- unlist(boot_estimates["growth_effect_tmle",])
    }
    
    ci_efficient_tmle <- quantile(growth_effect_tmle, p = c(0.025, 0.975), na.rm=TRUE)
    out$se_efficient_tmle <- sd(growth_effect_tmle, na.rm = TRUE)
    out$lower_ci_efficient_tmle <- ci_efficient_tmle[1]    
    out$upper_ci_efficient_tmle <- ci_efficient_tmle[2]
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
  
  return(out)
}
