#' Main function to get growth effect point estimate and standard error
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param G_name growth outcome variable name
#' @param V_name vaccination variable name
#' @param X_name covariate name(s)
#' @param Y_name infection variable name
#' @param est character vector of names of estimators to use for growth effect
#' @param n_boot number of bootstrap replicates
#' @param seed seet to set for replicability of bootstrap
#' @param return_se indicator to return closed form standard error for efficient_aipw or efficient_tmle, default FALSE
#' @param ml boolean to use SuperLearner models, default FALSE
#' @param G_V_X_model optional specify model to be used for fitting growth on vaccine + covariates, otherwise growth on all covariates
#' @param G_X_Y1_model optional specify model to be used for fitting growth on covariates in the infected, otherwise growth on all covariates
#' @param G_X_Y0_model optional specify model to be used for fitting growth on covariates in the uninfected, otherwise growth on all covariates
#' @param Y_X_model optional specify model to be used for fitting infection on covariates, otherwise infection on all covariates
#' @param G_V_X_library optional specify SuperLearner libraries for model fitting growth on covariates + vaccine, default glm
#' @param G_X_library optional specify SuperLearner libraries for model fitting growth on covariates, default glm
#' @param Y_X_library optional specify SuperLearner libraries for model fitting infection on covariates, default glm
#' @param null_hypothesis_value null value for hypothesis test effect for VE and population estimand g-comp, default 0
#' @param alpha_level alpha level for hypothesis testing, default 0.025
#' @param return_models boolean return models, default TRUE
#' @param family family for outcome variable 'G', defaults to gaussian for growth
#' @param v_folds number of cross validation folds for SuperLearner, default 3
#'
#' @export
#' 
#' @returns List of class `vegrowth`
vegrowth <- function(data,
                     G_name = "G",
                     V_name = "V",
                     X_name = "X",
                     Y_name = "Y",
                     est = c("gcomp_pop_estimand", "gcomp", "efficient_aipw", "efficient_tmle",
                             "choplump", "hudgens_lower", "hudgens_upper", "hudgens_adj_lower", "hudgens_adj_upper"),
                     n_boot = 1000, 
                     seed = 12345,
                     return_se = FALSE,
                     ml = FALSE,
                     G_V_X_model = NULL,
                     G_X_Y1_model = NULL,
                     G_X_Y0_model = NULL,
                     Y_X_model = NULL,
                     G_V_X_library = c("SL.glm"),
                     G_X_library = c("SL.glm"),
                     Y_X_library = c("SL.glm"),
                     null_hypothesis_value = 0,
                     alpha_level = 0.025,
                     return_models = TRUE,
                     family = "gaussian",
                     v_folds = 3){
  
  set.seed(seed)
  
  # ml_models only for AIPW and TMLE
  model_list <- list(models = NULL, 
                     ml_models = NULL)
  
  # Estimation methods requiring model fitting & bootstrap se
  if(any(est %in% c("gcomp_pop_estimand", "gcomp", "efficient_aipw", "efficient_tmle", "hudgens_adj_lower", "hudgens_adj_upper"))){
    
    # If ML specified and AIPW/TMLE are in est, fit ML models; otherwise only fit GLMs
    if(ml){
      if(any(est %in% c("efficient_aipw", "efficient_tmle"))){
        ml_models <- vegrowth::fit_ml_models(data = data, 
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
        
        model_list$ml_models <- ml_models
      } 
      
      if(any(est %in% c("gcomp_pop_estimand", "gcomp","hudgens_adj_lower", "hudgens_adj_upper"))){
        models <- vegrowth::fit_models(data = data, 
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
        
        model_list$models <- models
      }
    } else{
      models <- vegrowth::fit_models(data = data, 
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
      model_list$models <- models
    } 
      
    # If using return_se is true, do not use bootstrap se for AIPW and TMLE (remove from est for boot, used closed form SE for AIPW/TMLE)
    if(return_se == TRUE){
      bootstrap_results <- bootstrap_estimates(data = data, 
                                               G_name = G_name,
                                               V_name = V_name,
                                               Y_name = Y_name,
                                               X_name = X_name,
                                               n_boot = n_boot, 
                                               family = family,
                                               ml = ml,
                                               G_V_X_model = G_V_X_model,
                                               G_X_Y1_model = G_X_Y1_model,
                                               G_X_Y0_model = G_X_Y0_model,
                                               Y_X_model = Y_X_model,
                                               G_V_X_library = G_V_X_library,
                                               G_X_library = G_X_library,
                                               Y_X_library = Y_X_library,
                                               v_folds = v_folds,
                                               est = setdiff(est, c("efficient_aipw", "efficient_tmle")))
    } else{
      bootstrap_results <- bootstrap_estimates(data = data, 
                                               G_name = G_name,
                                               V_name = V_name,
                                               Y_name = Y_name,
                                               X_name = X_name,
                                               n_boot = n_boot, 
                                               family = family,
                                               ml = ml,
                                               G_V_X_model = G_V_X_model,
                                               G_X_Y1_model = G_X_Y1_model,
                                               G_X_Y0_model = G_X_Y0_model,
                                               Y_X_model = Y_X_model,
                                               G_V_X_library = G_V_X_library,
                                               G_X_library = G_X_library,
                                               Y_X_library = Y_X_library,
                                               v_folds = v_folds,
                                               est = est)
    }
  }
  
  # Point estimates for effects of interest & format results
  out <- list()
  
  if(return_models){
    out$models <- model_list
  }
  
  if("gcomp" %in% est){
    
    gcomp_res <- list()
    
    pt_est <- do_gcomp(data, models)
    gcomp_res$pt_est_additive <- pt_est['additive_effect']
    gcomp_res$pt_est_mult <- exp(pt_est['log_multiplicative_effect'])
    
    gcomp_res$se_additive <- bootstrap_results$se_gcomp_additive
    gcomp_res$se_log_mult <- bootstrap_results$se_gcomp_log_mult
    
    gcomp_res$lower_ci_additive <- bootstrap_results$lower_ci_gcomp_additive
    gcomp_res$upper_ci_additive <- bootstrap_results$upper_ci_gcomp_additive
    
    gcomp_res$lower_ci_mult <- bootstrap_results$lower_ci_gcomp_mult
    gcomp_res$upper_ci_mult <- bootstrap_results$upper_ci_gcomp_mult
    
    # just needed on an additive scale?
    gcomp_res$reject <- ((gcomp_res$pt_est_additive - null_hypothesis_value) / gcomp_res$se_additive) > qnorm(1 - alpha_level)
    
    class(gcomp_res) <- "gcomp_res"
    out$gcomp_res <- gcomp_res
    
  }
  if("gcomp_pop_estimand" %in% est){
    
    pop_gcomp_res <- list() 
    
    pt_est <- do_gcomp_pop_estimand(data, models, V_name = V_name, X_name = X_name)
    pop_gcomp_res$pt_est_additive <- pt_est['additive_effect']
    pop_gcomp_res$pt_est_mult <- exp(pt_est['log_multiplicative_effect'])
    
    pop_gcomp_res$se_additive <- bootstrap_results$se_gcomp_pop_estimand_additive
    pop_gcomp_res$se_log_mult <- bootstrap_results$se_gcomp_pop_estimand_log_mult
    
    pop_gcomp_res$lower_ci_additive <- bootstrap_results$lower_ci_gcomp_pop_estimand_additive
    pop_gcomp_res$upper_ci_additive <- bootstrap_results$upper_ci_gcomp_pop_estimand_additive
    
    pop_gcomp_res$lower_ci_mult <- bootstrap_results$lower_ci_gcomp_pop_estimand_mult
    pop_gcomp_res$upper_ci_mult <- bootstrap_results$upper_ci_gcomp_pop_estimand_mult
    
    pop_gcomp_res$reject <- ((pop_gcomp_res$pt_est_additive - null_hypothesis_value) / pop_gcomp_res$se_additive) > qnorm(1 - alpha_level)
    
    class(pop_gcomp_res) <- "pop_gcomp_res"
    out$pop_gcomp_res <- pop_gcomp_res
    
  }
  if("efficient_aipw" %in% est){
    
    aipw_res <- list()
    
    if(return_se == FALSE){
      # Point est + bootstrap SE
      
      if(ml){
        pt_est <- do_efficient_aipw(data, ml_models, G_name = G_name, V_name = V_name, Y_name = Y_name, return_se = return_se)
      } else{
        pt_est <- do_efficient_aipw(data, models, G_name = G_name, V_name = V_name, Y_name = Y_name, return_se = return_se)
      }
      
      aipw_res$pt_est_additive <- pt_est['additive_effect']
      aipw_res$pt_est_mult <- exp(pt_est['log_multiplicative_effect'])
      
      aipw_res$se_additive <- bootstrap_results$se_efficient_aipw_additive
      aipw_res$se_log_mult <- bootstrap_results$se_efficient_aipw_log_mult
      
      aipw_res$lower_ci_additive <- bootstrap_results$lower_ci_efficient_aipw_additive
      aipw_res$upper_ci_additive <- bootstrap_results$upper_ci_efficient_aipw_additive
      
      aipw_res$lower_ci_mult <- bootstrap_results$lower_ci_efficient_aipw_mult
      aipw_res$upper_ci_mult <- bootstrap_results$upper_ci_efficient_aipw_mult
      
      aipw_res$reject <- ((aipw_res$pt_est_additive - null_hypothesis_value) / aipw_res$se_additive) > qnorm(1 - alpha_level)
      
    } else {
      # Point est + closed form SE
      
      if(ml){
        aipw_result <- do_efficient_aipw(data, ml_models, G_name = G_name, V_name = V_name, Y_name = Y_name, return_se = return_se)
      } else{
        aipw_result <- do_efficient_aipw(data, models, G_name = G_name, V_name = V_name, Y_name = Y_name, return_se = return_se)
      }
      
      aipw_res$pt_est_additive <- aipw_result['additive_effect']
      aipw_res$pt_est_mult <- exp(aipw_result['log_multiplicative_effect'])
      
      aipw_res$se_additive <- aipw_result['additive_se']
      aipw_res$se_log_mult <- aipw_result['log_multiplicative_se']
      
      aipw_res$lower_ci_additive <- aipw_res$pt_est_additive - 1.96*aipw_res$se_additive
      aipw_res$upper_ci_additive <- aipw_res$pt_est_additive + 1.96*aipw_res$se_additive
      
      # CHECK THIS BRAIN NOT WORKING
      # already exponentiated the pt est
      aipw_res$lower_ci_mult <- exp(log(aipw_res$pt_est_mult) - 1.96*aipw_res$se_log_mult)
      aipw_res$upper_ci_mult <- exp(log(aipw_res$pt_est_mult) + 1.96*aipw_res$se_log_mult)
      
      aipw_res$reject <- ((aipw_res$pt_est_additive - null_hypothesis_value) / aipw_res$se_additive) > qnorm(1 - alpha_level)
    }
    
    class(aipw_res) <- "aipw_res"
    out$aipw_res <- aipw_res
    
  }
  if("efficient_tmle" %in% est){
    
    tmle_res <- list()
    
    if(return_se == FALSE){
      # Point est + bootstrap SE
      
      # STOPPED HERE IMPLEMENTING THE MULTIPLICATIVE EFFECT
      
      if(ml){
        pt_est <- do_efficient_tmle(data, ml_models, G_name = G_name, V_name = V_name, Y_name = Y_name, return_se = return_se)
      } else{
        pt_est <- do_efficient_tmle(data, models, G_name = G_name, V_name = V_name, Y_name = Y_name, return_se = return_se)
      }
      
      tmle_res$pt_est_additive <- pt_est['additive_effect']
      tmle_res$pt_est_mult <- exp(pt_est['log_multiplicative_effect'])
      
      tmle_res$se_additive <- bootstrap_results$se_efficient_tmle_additive
      tmle_res$se_log_mult <- bootstrap_results$se_efficient_tmle_log_mult
      
      tmle_res$lower_ci_additive <- bootstrap_results$lower_ci_efficient_tmle_additive
      tmle_res$upper_ci_additive <- bootstrap_results$upper_ci_efficient_tmle_additive
      
      tmle_res$lower_ci_mult <- bootstrap_results$lower_ci_efficient_tmle_mult
      tmle_res$upper_ci_mult <- bootstrap_results$upper_ci_efficient_tmle_mult
      
      tmle_res$reject <- ((tmle_res$pt_est_additive - null_hypothesis_value) / tmle_res$se_additive) > qnorm(1 - alpha_level)
      
    } else {
      # Point est + closed form SE
      
      if(ml){
        tmle_result <- do_efficient_tmle(data, ml_models, G_name = G_name, V_name = V_name, Y_name = Y_name, return_se = return_se)
      } else{
        tmle_result <- do_efficient_tmle(data, models, G_name = G_name, V_name = V_name, Y_name = Y_name, return_se = return_se)
      }
      
      tmle_res$pt_est_additive <- tmle_result['additive_effect']
      tmle_res$pt_est_mult <- exp(tmle_result['log_multiplicative_effect'])
      
      tmle_res$se_additive <- tmle_result['additive_se']
      tmle_res$se_log_mult <- tmle_result['log_multiplicative_se']
      
      tmle_res$lower_ci_additive <- tmle_res$pt_est_additive - 1.96*tmle_res$se_additive
      tmle_res$upper_ci_additive <- tmle_res$pt_est_additive + 1.96*tmle_res$se_additive
      
      # CHECK THIS BRAIN NOT WORKING
      # already exponentiated the pt est
      tmle_res$lower_ci_mult <- exp(log(tmle_res$pt_est_mult) - 1.96*tmle_res$se_log_mult)
      tmle_res$upper_ci_mult <- exp(log(tmle_res$pt_est_mult) + 1.96*tmle_res$se_log_mult)
      
      tmle_res$reject <- ((tmle_res$pt_est_additive - null_hypothesis_value) / aipw_res$se_additive) > qnorm(1 - alpha_level)
    }
    
    class(tmle_res) <- "tmle_res"
    out$tmle_res <- tmle_res
    
  }
  if("choplump" %in% est){
    choplump_rslt <- do_chop_lump_test(data, G_name = G_name, V_name = V_name, Y_name = Y_name)
    choplump_rslt$reject <- choplump_rslt$pval < 0.05
    
    class(choplump_rslt) <- "choplump_res"
    out$chop_lump_res <- choplump_rslt
  }
  if("hudgens_lower" %in% est){
    hudgens_rslt_lower <- hudgens_test(data, G_name = G_name, V_name = V_name, Y_name = Y_name, lower_bound = TRUE)
    hudgens_rslt_lower$reject <- hudgens_rslt_lower$pval < 0.05
    
    hudgens_rslt_lower$se <- bootstrap_results$se_hudgens_lower
    hudgens_rslt_lower$lower_ci <- bootstrap_results$lower_ci_hudgens_lower
    hudgens_rslt_lower$upper_ci <- bootstrap_results$upper_ci_hudgens_lower
    #hudgens_rslt_lower$reject <- ((hudgens_rslt_lower$pt_est - null_hypothesis_value) / hudgens_rslt_lower$se) > qnorm(1 - alpha_level)
    
    class(hudgens_rslt_lower) <- "hudgens_lower_res"
    out$hudgens_rslt_lower <- hudgens_rslt_lower
  }
  if("hudgens_upper" %in% est){
    hudgens_rslt_upper <- hudgens_test(data, G_name = G_name, V_name = V_name, Y_name = Y_name, lower_bound = FALSE)
    hudgens_rslt_upper$reject <- hudgens_rslt_upper$pval < 0.05
    
    hudgens_rslt_upper$se <- bootstrap_results$se_hudgens_lower
    hudgens_rslt_upper$lower_ci <- bootstrap_results$lower_ci_hudgens_upper
    hudgens_rslt_upper$upper_ci <- bootstrap_results$upper_ci_hudgens_upper
    #hudgens_rslt_upper$reject <- ((hudgens_rslt_upper$pt_est - null_hypothesis_value) / hudgens_rslt_upper$se) > qnorm(1 - alpha_level)
    
    class(hudgens_rslt_upper) <- "hudgens_upper_res"
    out$hudgens_rslt_upper <- hudgens_rslt_upper
  }
  if("hudgens_lower_doomed" %in% est){
    hudgens_rslt_lower_doomed <- hudgens_test_doomed(data, G_name = G_name, V_name = V_name, Y_name = Y_name, lower_bound = TRUE)
    hudgens_rslt_lower_doomed$reject <- hudgens_rslt_lower_doomed$pval < 0.05
    
    hudgens_rslt_lower_doomed$se <- bootstrap_results$se_hudgens_lower_doomed
    hudgens_rslt_lower_doomed$lower_ci <- bootstrap_results$lower_ci_hudgens_lower_doomed
    hudgens_rslt_lower_doomed$upper_ci <- bootstrap_results$upper_ci_hudgens_lower_doomed
    
    class(hudgens_rslt_lower_doomed) <- "hudgens_lower_res_doomed"
    out$hudgens_rslt_lower_doomed <- hudgens_rslt_lower_doomed
  }
  if("hudgens_upper_doomed" %in% est){
    hudgens_rslt_upper_doomed <- hudgens_test_doomed(data, G_name = G_name, V_name = V_name, Y_name = Y_name, lower_bound = FALSE)
    hudgens_rslt_upper_doomed$reject <- hudgens_rslt_upper_doomed$pval < 0.05
    
    hudgens_rslt_upper_doomed$se <- bootstrap_results$se_hudgens_upper_doomed
    hudgens_rslt_upper_doomed$lower_ci <- bootstrap_results$lower_ci_hudgens_upper_doomed
    hudgens_rslt_upper_doomed$upper_ci <- bootstrap_results$upper_ci_hudgens_upper_doomed
    
    class(hudgens_rslt_upper_doomed) <- "hudgens_upper_res_doomed"
    out$hudgens_rslt_upper_doomed <- hudgens_rslt_upper_doomed
  }
  if("hudgens_adj_lower" %in% est){
    
    hudgens_adj_rslt_lower <- list()
    
    hudgens_adj_rslt_lower$pt_est <- get_adjusted_hudgens_stat(data, models, family, lower_bound = TRUE)
    hudgens_adj_rslt_lower$se <- bootstrap_results$se_hudgens_adj_lower
    hudgens_adj_rslt_lower$lower_ci <- bootstrap_results$lower_ci_hudgens_adj_lower
    hudgens_adj_rslt_lower$upper_ci <- bootstrap_results$upper_ci_hudgens_adj_lower
    hudgens_adj_rslt_lower$reject <- ((hudgens_adj_rslt_lower$pt_est - null_hypothesis_value) / hudgens_adj_rslt_lower$se) > qnorm(1 - alpha_level)
    
    class(hudgens_adj_rslt_lower) <- "hudgens_adj_lower_res"
    out$hudgens_adj_rslt_lower <- hudgens_adj_rslt_lower
  }
  if("hudgens_adj_lower" %in% est){
    
    hudgens_adj_rslt_upper <- list()
    
    hudgens_adj_rslt_upper$pt_est <- get_adjusted_hudgens_stat(data, models, family, lower_bound = FALSE)
    hudgens_adj_rslt_upper$se <- bootstrap_results$se_hudgens_adj_upper
    hudgens_adj_rslt_upper$lower_ci <- bootstrap_results$lower_ci_hudgens_adj_upper
    hudgens_adj_rslt_upper$upper_ci <- bootstrap_results$upper_ci_hudgens_adj_upper
    hudgens_adj_rslt_upper$reject <- ((hudgens_adj_rslt_upper$pt_est - null_hypothesis_value) / hudgens_adj_rslt_upper$se) > qnorm(1 - alpha_level)
    
    class(hudgens_adj_rslt_upper) <- "hudgens_adj_upper_res"
    out$hudgens_adj_rslt_upper <- hudgens_adj_rslt_upper
  }
  
  class(out) <- "vegrowth"
  
  return(out)

}