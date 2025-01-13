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
#' 
#' @returns list containing results for specified estimators on single bootstrap sample
one_boot <- function(
    data,
    G_name = "G",
    V_name = "V",
    X_name = "X",
    Y_name = "Y", 
    est = c("gcomp_pop_estimand", "gcomp", "efficient_aipw", "efficient_tmle"),
    G_X_model = NULL,
    Y_X_model = NULL
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
                            Y_X_model = Y_X_model)
  
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
#' 
#' @returns list containing bootstrap se and 95% CI bounds for estimators specified in est
bootstrap_estimates <- function(
    data, 
    G_name = "G",
    V_name = "V",
    X_name = "X",
    Y_name = "Y", 
    n_boot = 1000, 
    est = c("gcomp_pop_estimand", "gcomp", "efficient_aipw", "efficient_tmle"),
    G_X_model = NULL, 
    Y_X_model = NULL
){
  
  boot_estimates <- replicate(n_boot, one_boot(data, 
                                               G_name = G_name,
                                               V_name = V_name,
                                               Y_name = Y_name,
                                               X_name = X_name,
                                               est = est,
                                               G_X_model = G_X_model, 
                                               Y_X_model = Y_X_model))
  
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
  
  return(out)
}




# ------------------------------------------------------------

# Not needed?

#' Function to run simulation for single seed 
#' 
#' @param params parameters used for data generation
#' @param n_boot number of bootstrap replicates
#' @param est character vector of names of estimators to use for growth effect
#' @param null_hypothesis_value value of null hypothesis for growth effect hypothesis testing, default 0
#' @param alpha_level one-sided alpha for hypothesis testing, default 0.025
#' 
#' @returns dataframe with parameters, point estimates, confidence intervals for specified estimators
do_one_ci_sim <- function(
    params, n_boot = 1000, 
    est = c("gcomp_pop_estimand", "gcomp",
            "efficient_aipw", "efficient_tmle", "choplump"),
    null_hypothesis_value = 0,
    alpha_level = 0.025
){
  set.seed(params$seed)
  
  # data <- generate_EFGH_complex(n = params$n,
  #                               incidence_shigella = params$incidence_shigella,
  #                               incidence_severe_shigella = params$incidence_severe_shigella,
  #                               VE_mild = params$VE_mild,
  #                               VE_severe = params$VE_severe,
  #                               lfaz_bl_effect_shig = params$lfaz_bl_effect_shig,
  #                               lfaz_bl_effect_severe_shig = params$lfaz_bl_effect_severe_shig,
  #                               shig_coef_mild = params$shig_coef_mild,
  #                               shig_coef_severe = params$shig_coef_severe,
  #                               sd_growth = params$sd,
  #                               zhifei = params$zhifei,
  #                               catch_up = params$catch_up)
  
  data <- generate_EFGH_complex_quad(n = params$n,
                                     incidence_shigella = params$incidence_shigella,
                                     incidence_severe_shigella = params$incidence_severe_shigella,
                                     VE_mild = params$VE_mild,
                                     VE_severe = params$VE_severe,
                                     lfaz_bl_effect_shig = params$lfaz_bl_effect_shig,
                                     lfaz_bl_effect_severe_shig = params$lfaz_bl_effect_severe_shig,
                                     shig_coef_mild = params$shig_coef_mild,
                                     shig_coef_severe = params$shig_coef_severe,
                                     sd_growth = params$sd,
                                     zhifei = params$zhifei,
                                     catch_up = params$catch_up)
  
  models <- vegrowth::fit_models(data = data, 
                       est = est,
                       G_X_model = params$G_X_model,
                       Yinf_X_model = params$Yinf_X_model)
  
  bootstrap_results <- bootstrap_estimates(data, n_boot, est)
  
  # point estimates for effects of interest
  out <- list()
  if("gcomp" %in% est){
    out$pt_est_ge <- do_gcomp(data, models)
    out$bootstrap_se_ge <- bootstrap_results$se_gcomp
    out$lower_ci_ge <- bootstrap_results$lower_ci_gcomp
    out$upper_ci_ge <- bootstrap_results$upper_ci_gcomp
    out$reject_ge <- ((out$pt_est_ge - null_hypothesis_value) / out$bootstrap_se_ge) > qnorm(1 - alpha_level)
  }
  if("efficient_gcomp" %in% est){
    out$pt_est_ge_e <- do_efficient_gcomp(data, models)
    out$bootstrap_se_ge_e <- bootstrap_results$se_efficient_gcomp
    out$lower_ci_ge_e <- bootstrap_results$lower_ci_efficient_gcomp
    out$upper_ci_ge_e <- bootstrap_results$upper_ci_efficient_gcomp
    out$reject_ge_e <- ((out$pt_est_ge_e - null_hypothesis_value) / out$bootstrap_se_ge_e) > qnorm(1 - alpha_level)
  }
  if("efficient_aipw" %in% est){
    out$pt_est_aipw <- do_efficient_aipw(data, models)
    out$bootstrap_se_aipw <- bootstrap_results$se_efficient_aipw
    out$lower_ci_aipw <- bootstrap_results$lower_ci_efficient_aipw
    out$upper_ci_aipw <- bootstrap_results$upper_ci_efficient_aipw
    out$reject_ge_aipw <- ((out$pt_est_aipw - null_hypothesis_value) / out$bootstrap_se_aipw) > qnorm(1 - alpha_level)
  }
  if("efficient_tmle" %in% est){
    out$pt_est_tmle <- do_efficient_tmle(data, models)
    out$bootstrap_se_tmle <- bootstrap_results$se_efficient_tmle
    out$lower_ci_tmle <- bootstrap_results$lower_ci_efficient_tmle
    out$upper_ci_tmle <- bootstrap_results$upper_ci_efficient_tmle
    out$reject_ge_tmle <- ((out$pt_est_tmle - null_hypothesis_value) / out$bootstrap_se_tmle) > qnorm(1 - alpha_level)
  }
  if("gcomp_pop_estimand" %in% est){
    out$pt_est_pop <- do_gcomp_pop_estimand(data, models)
    out$bootstrap_se_pop <- bootstrap_results$se_gcomp_pop_estimand
    out$lower_ci_pop <- bootstrap_results$lower_ci_gcomp_pop_estimand
    out$upper_ci_pop <- bootstrap_results$upper_ci_gcomp_pop_estimand
    out$reject_ge_pop <- ((out$pt_est_pop - null_hypothesis_value) / out$bootstrap_se_pop) > qnorm(1 - alpha_level)
  }
  if("choplump" %in% est){
    choplump_rslt <- do_chop_lump_test(data)
    out$reject_ge_chop_lump <- choplump_rslt$pval < 0.05
  }
  if("hudgens_lower" %in% est){
    hudgens_rslt_lower <- hudgens_test(data, lower_bound = TRUE)#, n_boot = n_boot)
    out$reject_hudgens_lower <- hudgens_rslt_lower$pval < 0.05
  }
  if("hudgens_upper" %in% est){
    hudgens_rslt_upper <- hudgens_test(data, lower_bound = FALSE)#, n_boot = n_boot)
    out$reject_hudgens_upper <- hudgens_rslt_upper$pval < 0.05
  }
  
  return(cbind(params, out))
  
}
