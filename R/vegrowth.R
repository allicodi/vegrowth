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
#' @param G_X_model optional specify model to be used for fitting growth on covariates, otherwise growth on all covariates
#' @param Y_X_model optional specify model to be used for fitting infection on covariates, otherwise infection on all covariates
#' @param null_hypothesis_value null value for hypothesis test effect for VE and population estimand g-comp, default 0
#' @param alpha_level alpha level for hypothesis testing, default 0.025
#' @param return_models boolean return models, default TRUE
#' @param family family for outcome variable 'G', defaults to gaussian for growth
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
                             "choplump", "hudgens_lower", "hudgens_upper"),
                     n_boot = 1000, 
                     seed = 12345,
                     return_se = FALSE,
                     G_X_model = NULL,
                     Y_X_model = NULL, 
                     null_hypothesis_value = 0,
                     alpha_level = 0.025,
                     return_models = TRUE,
                     family = "gaussian"){
  
  set.seed(seed)
  
  # Estimation methods requiring model fitting & bootstrap se
  if(any(est %in% c("gcomp_pop_estimand", "gcomp", "efficient_aipw", "efficient_tmle"))){
    models <- vegrowth::fit_models(data = data, 
                                   est = est, 
                                   G_name = G_name,
                                   V_name = V_name,
                                   Y_name = Y_name,
                                   X_name = X_name,
                                   G_X_model = G_X_model,
                                   Y_X_model = Y_X_model,
                                   family = family)
    
    # If using return_se is true, do not use bootstrap se for AIPW and TMLE (remove from est for boot)
    if(return_se == TRUE){
      bootstrap_results <- bootstrap_estimates(data = data, 
                                               G_name = G_name,
                                               V_name = V_name,
                                               Y_name = Y_name,
                                               X_name = X_name,
                                               n_boot = n_boot, 
                                               family = family,
                                               est = setdiff(est, c("efficient_aipw", "efficient_tmle")))
    } else{
      bootstrap_results <- bootstrap_estimates(data = data, 
                                               G_name = G_name,
                                               V_name = V_name,
                                               Y_name = Y_name,
                                               X_name = X_name,
                                               n_boot = n_boot, 
                                               family = family,
                                               est = est)
    }
  }
  
  # Point estimates for effects of interest & format results
  out <- list()
  
  if(return_models){
    out$models <- models
  }
  
  if("gcomp" %in% est){
    
    gcomp_res <- list()
    
    gcomp_res$pt_est <- do_gcomp(data, models)
    gcomp_res$se <- bootstrap_results$se_gcomp
    gcomp_res$lower_ci <- bootstrap_results$lower_ci_gcomp
    gcomp_res$upper_ci <- bootstrap_results$upper_ci_gcomp
    gcomp_res$reject <- ((gcomp_res$pt_est - null_hypothesis_value) / gcomp_res$se) > qnorm(1 - alpha_level)
    
    class(gcomp_res) <- "gcomp_res"
    out$gcomp_res <- gcomp_res
    
  }
  if("gcomp_pop_estimand" %in% est){
    
    pop_gcomp_res <- list() 
    
    pop_gcomp_res$pt_est <- do_gcomp_pop_estimand(data, models, V_name = V_name, X_name = X_name)
    pop_gcomp_res$se <- bootstrap_results$se_gcomp_pop_estimand
    pop_gcomp_res$lower_ci <- bootstrap_results$lower_ci_gcomp_pop_estimand
    pop_gcomp_res$upper_ci <- bootstrap_results$upper_ci_gcomp_pop_estimand
    pop_gcomp_res$reject <- ((pop_gcomp_res$pt_est - null_hypothesis_value) / pop_gcomp_res$se) > qnorm(1 - alpha_level)
    
    class(pop_gcomp_res) <- "pop_gcomp_res"
    out$pop_gcomp_res <- pop_gcomp_res
    
  }
  if("efficient_aipw" %in% est){
    
    aipw_res <- list()
    
    if(return_se == FALSE){
      # Point est + bootstrap SE
      aipw_res$pt_est <- do_efficient_aipw(data, models, G_name = G_name, X_name = X_name, V_name = V_name, Y_name = Y_name, return_se = return_se)
      aipw_res$se <- bootstrap_results$se_efficient_aipw
      aipw_res$lower_ci <- bootstrap_results$lower_ci_efficient_aipw
      aipw_res$upper_ci <- bootstrap_results$upper_ci_efficient_aipw
      aipw_res$reject <- ((aipw_res$pt_est - null_hypothesis_value) / aipw_res$se) > qnorm(1 - alpha_level)
    } else {
      # Point est + closed form SE
      aipw_res <- do_efficient_aipw(data, models, G_name = G_name, X_name = X_name, V_name = V_name, Y_name = Y_name, return_se = return_se)
      aipw_res$pt_est <- aipw_res[1]
      aipw_res$se <- aipw_res[2]
      aipw_res$lower_ci <- aipw_res[1] - 1.96*aipw_res[2]
      aipw_res$upper_ci <- aipw_res[1] + 1.96*aipw_res[2]
      aipw_res$reject <- ((aipw_res$pt_est - null_hypothesis_value) / aipw_res$se) > qnorm(1 - alpha_level)
    }
    
    class(aipw_res) <- "aipw_res"
    out$aipw_res <- aipw_res
    
  }
  if("efficient_tmle" %in% est){
    
    tmle_res <- list()
    
    if(return_se == FALSE){
      # Point est + bootstrap SE
      tmle_res$pt_est <- do_efficient_tmle(data, models, G_name = G_name, V_name = V_name, Y_name = Y_name, return_se = return_se)
      tmle_res$se <- bootstrap_results$se_efficient_tmle
      tmle_res$lower_ci <- bootstrap_results$lower_ci_efficient_tmle
      tmle_res$upper_ci <- bootstrap_results$upper_ci_efficient_tmle
      tmle_res$reject <- ((tmle_res$pt_est - null_hypothesis_value) / tmle_res$se) > qnorm(1 - alpha_level)
    } else {
      # Point est + closed form SE
      tmle_res <- do_efficient_tmle(data, models, G_name = G_name, V_name = V_name, Y_name = Y_name, return_se = return_se)
      tmle_res$pt_est <- tmle_res[1]
      tmle_res$se <- tmle_res[2]
      tmle_res$lower_ci <- tmle_res[1] - 1.96*tmle_res[2]
      tmle_res$upper_ci <- tmle_res[1] + 1.96*tmle_res[2]
      tmle_res$reject <- ((tmle_res$pt_est - null_hypothesis_value) / tmle_res$se) > qnorm(1 - alpha_level)
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
    
    class(hudgens_rslt_lower) <- "hudgens_lower_res"
    out$hudgens_rslt_lower <- hudgens_rslt_lower
  }
  if("hudgens_upper" %in% est){
    hudgens_rslt_upper <- hudgens_test(data, G_name = G_name, V_name = V_name, Y_name = Y_name, lower_bound = FALSE)
    hudgens_rslt_upper$reject <- hudgens_rslt_upper$pval < 0.05
    
    class(hudgens_rslt_upper) <- "hudgens_upper_res"
    out$hudgens_rslt_upper <- hudgens_rslt_upper
  }
  
  class(out) <- "vegrowth"
  
  return(out)

}