#' Main function to get growth effect point estimate and standard error
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param Y_name growth outcome variable name
#' @param Z_name vaccination variable name
#' @param X_name covariate name(s)
#' @param S_name infection variable name
#' @param estimand character vector with name(s) of estimands of interest; "nat_inf" = naturally infected, "doomed" = doomed, "pop" = marginal/population-level
#' @param method character vector with name(s) of methods to use for estimation; "gcomp" = g-computation, "ipw" = inverese probability weighting, "aipw" = augmented inverse probability weighting, "tmle" = targeted maximum likelihood estimation (nat_inf only), "bound" = bounds without cross-world assumptions (nat_inf and doomed only), "sens" = sensitivity analysis (nat_inf only)
#' @param exclusion_restriction boolean or vector of boolean (TRUE,FALSE) indicating version of naturally infected estimators with and/or without exclusion restriction assumptions, default FALSE
#' @param two_part_model If \code{exclusion_restriction} is \code{TRUE} or \code{"pop"} is included in \code{estimand}, should E(Y | Z, X) be estimated using separate models for E(Y | Z, X, S) and P(S | Z, X) (if \code{TRUE}) or using a single model for E(Y | Z, X) (if \code{FALSE}). Currently, this is only implemented for \code{'gcomp'} and \code{'aipw'} naturally infected estimators
#' @param n_boot number of bootstrap replicates
#' @param seed seet to set for replicability of bootstrap
#' @param return_se indicator to return closed form standard error for efficient_aipw or efficient_tmle, default TRUE
#' @param ml boolean to use SuperLearner models for AIPW & TMLE, default TRUE
#' @param Y_Z_X_model optional specify model to be used for fitting growth on vaccine + covariates, otherwise growth on all covariates
#' @param Y_X_S1_model optional specify model to be used for fitting growth on covariates in the infected, otherwise growth on all covariates
#' @param Y_X_S0_model optional specify model to be used for fitting growth on covariates in the uninfected, otherwise growth on all covariates
#' @param S_X_model optional specify model to be used for fitting infection on covariates, otherwise infection on all covariates
#' @param Y_Z_X_library optional specify SuperLearner libraries for model fitting growth on covariates + vaccine, default glm
#' @param Y_X_library optional specify SuperLearner libraries for model fitting growth on covariates, default glm
#' @param S_X_library optional specify SuperLearner libraries for model fitting infection on covariates, default glm
#' @param null_hypothesis_value null value for hypothesis test effect for ZE and population estimand g-comp, default 0
#' @param alpha_level alpha level for hypothesis testing, default 0.025
#' @param return_models boolean return models, default TRUE
#' @param family family for outcome variable 'G', defaults to gaussian for growth
#' @param v_folds number of cross validation folds for SuperLearner, default 3
#' @param effect_dir direction of beneficial effect, defaults to "positive" for beneficial outcome. Used for one-side tests of bounds.  
#' @param epsilon a vector of values for the sensitivity parameter (only applicable for method = "sens")
#'
#' @export
#' 
#' @returns List of class `vegrowth`
vegrowth <- function(data,
                     Y_name = "G",
                     Z_name = "Z",
                     X_name = "X",
                     S_name = "S",
                     estimand = c("nat_inf", "doomed", "pop"),
                     method = c("gcomp", "ipw", "aipw", "tmle", "bound", "sens"),
                     exclusion_restriction = FALSE,
                     two_part_model = FALSE,
                     n_boot = 1000,
                     permutation = FALSE,
                     n_perm = 1000,
                     seed = 12345,
                     return_se = TRUE,
                     ml = TRUE,
                     Y_Z_X_model = NULL,
                     Y_X_S1_model = NULL,
                     Y_X_S0_model = NULL,
                     S_X_model = NULL,
                     S_Z_X_model = NULL,
                     Z_X_model = paste0(Z_name, " ~ 1"),
                     Y_Z_X_library = c("SL.glm"),
                     Y_X_library = c("SL.glm"),
                     S_X_library = c("SL.glm"),
                     S_Z_X_library = c("SL.glm"),
                     Z_X_library = c("SL.mean"),
                     null_hypothesis_value = 0,
                     alpha_level = 0.05,
                     return_models = TRUE,
                     family = "gaussian",
                     v_folds = 3,
                     effect_dir = "positive",
                     epsilon = exp(seq(log(0.5), log(2), length = 49))){
  
  set.seed(seed)
  
  # ----------------------------------------------------------------------------
  # 1. Model Fitting -----------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  # ml_models only for AIPW and TMLE
  model_list <- list(models = NULL, 
                     ml_models = NULL)
  
  # Estimation methods requiring model fitting (everything but bounds)
  if(any(method %in% c("gcomp", "ipw", "aipw", "tmle", "sens"))){
    
    # If ML specified and aipw, tmle, and/or sens are in method, fit ML models; otherwise only fit GLMs
    if(ml){
      if(any(method %in% c("aipw", "tmle", "sens"))){
        ml_models <- vegrowth::fit_ml_models(data = data, 
                                             estimand = estimand,
                                             method = method,
                                             exclusion_restriction = exclusion_restriction,
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
        
        model_list$ml_models <- ml_models
      } 
      
      # still force gcomp & ipw to be glm only
      if(any(method %in% c("gcomp", "ipw"))){
        models <- vegrowth::fit_models(data = data, 
                                       estimand = estimand,
                                       method = method,
                                       exclusion_restriction = exclusion_restriction,
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
        
        model_list$models <- models
      }
    } else{
      # ML not specified; use glms for all 
      
      models <- vegrowth::fit_models(data = data, 
                                     estimand = estimand,
                                     method = method,
                                     exclusion_restriction = exclusion_restriction,
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
      model_list$models <- models
    } 
  }
  
  # ----------------------------------------------------------------------------
  # 2. Bootstrap standard error & confidence intervals -------------------------
  # ----------------------------------------------------------------------------
 out <- bootstrap_estimates(data = data, 
                             Y_name = Y_name,
                             Z_name = Z_name,
                             S_name = S_name,
                             X_name = X_name,
                             n_boot = n_boot, 
                             family = family,
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
                             estimand = estimand, 
                             method = method, 
                             exclusion_restriction = exclusion_restriction,
                             two_part_model = two_part_model,
                             effect_dir = effect_dir,
                             epsilon = epsilon,
                             return_se = return_se)
 
 # ----------------------------------------------------------------------------
 # 3. Point estimates for effects of interest & tests
 # ----------------------------------------------------------------------------
  
  # Naturally infected --------------------------------------------------------
  
  if("nat_inf" %in% estimand){
    
    # For each exclusion restriction option (TRUE, FALSE)
    for(er in exclusion_restriction){
      
      # Exclustion restriction -- save as <estimator>_ER 
      er_suffix <- if (er) "_ER" else ""
      
        if("gcomp" %in% method){
          
          estimator <- paste0("gcomp", er_suffix)
          
          out$nat_inf[[estimator]]$pt_est <- do_gcomp_nat_inf(data = data, models = models, Z_name = Z_name, X_name = X_name, exclusion_restriction = er, two_part_model = two_part_model)
          
          out$nat_inf[[estimator]]$test_stat$additive <- (out$nat_inf[[estimator]]$pt_est['additive_effect'] - null_hypothesis_value) / 
            out$nat_inf[[estimator]]$boot_se$se_additive
          
          out$nat_inf[[estimator]]$test_stat$mult <- (out$nat_inf[[estimator]]$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
            out$nat_inf[[estimator]]$boot_se$se_log_mult
          
          out$nat_inf[[estimator]]$p_val$additive <- 2 * (1 - pnorm(abs(out$nat_inf[[estimator]]$test_stat$additive)))
          out$nat_inf[[estimator]]$p_val$mult <- 2 * (1 - pnorm(abs(out$nat_inf[[estimator]]$test_stat$mult)))
          
          out$nat_inf[[estimator]]$reject$additive <- (abs(out$nat_inf[[estimator]]$pt_est['additive_effect'] - null_hypothesis_value) / 
                                                  out$nat_inf[[estimator]]$boot_se$se_additive) > qnorm(1 - alpha_level/2)
          out$nat_inf[[estimator]]$reject$mult <- (abs(out$nat_inf[[estimator]]$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                              out$nat_inf[[estimator]]$boot_se$se_log_mult) > qnorm(1 - alpha_level/2)
          
          class(out$nat_inf[[estimator]]) <- estimator
        }
        
        if("ipw" %in% method){
          
          estimator <- paste0("ipw", er_suffix)
          
          out$nat_inf[[estimator]]$pt_est <- do_ipw_nat_inf(data = data, models = models, S_name = S_name, Y_name = Y_name, Z_name = Z_name, exclusion_restriction = er)
          
          out$nat_inf[[estimator]]$test_stat$additive <- (out$nat_inf[[estimator]]$pt_est['additive_effect'] - null_hypothesis_value) / 
            out$nat_inf[[estimator]]$boot_se$se_additive
          
          out$nat_inf[[estimator]]$test_stat$mult <- (out$nat_inf[[estimator]]$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
            out$nat_inf[[estimator]]$boot_se$se_log_mult
          
          out$nat_inf[[estimator]]$p_val$additive <- 2 * (1 - pnorm(abs(out$nat_inf[[estimator]]$test_stat$additive)))
          out$nat_inf[[estimator]]$p_val$mult <- 2 * (1 - pnorm(abs(out$nat_inf[[estimator]]$test_stat$mult)))
          
          out$nat_inf[[estimator]]$reject$additive <- (abs(out$nat_inf[[estimator]]$pt_est['additive_effect'] - null_hypothesis_value) / 
                                                out$nat_inf[[estimator]]$boot_se$se_additive) > qnorm(1 - alpha_level/2)
          out$nat_inf[[estimator]]$reject$mult <- (abs(out$nat_inf[[estimator]]$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                            out$nat_inf[[estimator]]$boot_se$se_log_mult) > qnorm(1 - alpha_level/2)
          
          class(out$nat_inf[[estimator]]) <- estimator
        }
        
        if("aipw" %in% method){
          
          estimator <- paste0("aipw", er_suffix)
          
          if(ml){
            out$nat_inf[[estimator]]$pt_est <- do_aipw_nat_inf(data = data, models = ml_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, X_name = X_name, return_se = return_se, exclusion_restriction = er, two_part_model = two_part_model)
          } else{
            out$nat_inf[[estimator]]$pt_est <- do_aipw_nat_inf(data = data, models = models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, X_name = X_name, return_se = return_se, exclusion_restriction = er, two_part_model = two_part_model)
          }
          
          if(is.null(out$nat_inf[[estimator]]$boot_se)){
            # closed form SE
            
            out$nat_inf[[estimator]]$test_stat$additive <- (out$nat_inf[[estimator]]$pt_est['additive_effect'] - null_hypothesis_value) / 
              out$nat_inf[[estimator]]$pt_est['additive_se']
            
            out$nat_inf[[estimator]]$test_stat$mult <- (out$nat_inf[[estimator]]$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
              out$nat_inf[[estimator]]$pt_est['log_multiplicative_se']
            
            out$nat_inf[[estimator]]$p_val$additive <- 2 * (1 - pnorm(abs(out$nat_inf[[estimator]]$test_stat$additive)))
            out$nat_inf[[estimator]]$p_val$mult <- 2 * (1 - pnorm(abs(out$nat_inf[[estimator]]$test_stat$mult)))
            
            out$nat_inf[[estimator]]$reject$additive <- (abs(out$nat_inf[[estimator]]$pt_est['additive_effect'] - null_hypothesis_value) / 
                                                   out$nat_inf[[estimator]]$pt_est['additive_se']) > qnorm(1 - alpha_level / 2)
            out$nat_inf[[estimator]]$reject$mult <- (abs(out$nat_inf[[estimator]]$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                               out$nat_inf[[estimator]]$pt_est['log_multiplicative_se']) > qnorm(1 - alpha_level / 2)
          } else{
            # bootstrap se
            
            out$nat_inf[[estimator]]$test_stat$additive <- (out$nat_inf[[estimator]]$pt_est['additive_effect'] - null_hypothesis_value) / 
              out$nat_inf[[estimator]]$boot_se$se_additive
            
            out$nat_inf[[estimator]]$test_stat$mult <- (out$nat_inf[[estimator]]$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
              out$nat_inf[[estimator]]$boot_se$se_log_mult
            
            out$nat_inf[[estimator]]$p_val$additive <- 2 * (1 - pnorm(abs(out$nat_inf[[estimator]]$test_stat$additive)))
            out$nat_inf[[estimator]]$p_val$mult <- 2 * (1 - pnorm(abs(out$nat_inf[[estimator]]$test_stat$mult)))
            
            out$nat_inf[[estimator]]$reject$additive <- (abs(out$nat_inf[[estimator]]$pt_est['additive_effect'] - null_hypothesis_value) / 
                                                   out$nat_inf[[estimator]]$boot_se$se_additive) > qnorm(1 - alpha_level / 2)
            out$nat_inf[[estimator]]$reject$mult <- (abs(out$nat_inf[[estimator]]$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                               out$nat_inf[[estimator]]$boot_se$se_log_mult) > qnorm(1 - alpha_level / 2)
            
          }
          
          class(out$nat_inf[[estimator]]) <- estimator

        }
      }
      
    if("tmle" %in% method){
      if(ml){
        out$nat_inf$tmle$pt_est <- do_tmle_nat_inf(data = data, models = ml_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, return_se = return_se)
      } else{
        out$nat_inf$tmle$pt_est <- do_tmle_nat_inf(data = data, models = models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, return_se = return_se)
      }
      
      if(is.null(out$nat_inf$tmle$boot_se)){
        # closed form SE
        
        out$nat_inf$tmle$test_stat$additive <- (out$nat_inf$tmle$pt_est['additive_effect'] - null_hypothesis_value) / 
          out$nat_inf$tmle$pt_est['additive_se']
        
        out$nat_inf$tmle$test_stat$mult <- (out$nat_inf$tmle$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
          out$nat_inf$tmle$pt_est['log_multiplicative_se']
        
        out$nat_inf$tmle$p_val$additive <- 2 * (1 - pnorm(abs(out$nat_inf$tmle$test_stat$additive)))
        out$nat_inf$tmle$p_val$mult <- 2 * (1 - pnorm(abs(out$nat_inf$tmle$test_stat$mult)))
        
        out$nat_inf$tmle$reject$additive <- (abs(out$nat_inf$tmle$pt_est['additive_effect'] - null_hypothesis_value) / 
                                               out$nat_inf$tmle$pt_est['additive_se']) > qnorm(1 - alpha_level / 2)
        out$nat_inf$tmle$reject$mult <- (abs(out$nat_inf$tmle$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                           out$nat_inf$tmle$pt_est['log_multiplicative_se']) > qnorm(1 - alpha_level / 2)
      } else{
        # bootstrap se
        
        out$nat_inf$tmle$test_stat$additive <- (out$nat_inf$tmle$pt_est['additive_effect'] - null_hypothesis_value) / 
          out$nat_inf$tmle$boot_se$se_additive
        
        out$nat_inf$tmle$test_stat$mult <- (out$nat_inf$tmle$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
          out$nat_inf$tmle$boot_se$se_log_mult
        
        out$nat_inf$tmle$p_val$additive <- 2 * (1 - pnorm(abs(out$nat_inf$tmle$test_stat$additive)))
        out$nat_inf$tmle$p_val$mult <- 2 * (1 - pnorm(abs(out$nat_inf$tmle$test_stat$mult)))
        
        out$nat_inf$tmle$reject$additive <- (abs(out$nat_inf$tmle$pt_est['additive_effect'] - null_hypothesis_value) / 
                                               out$nat_inf$tmle$boot_se$se_additive) > qnorm(1 - alpha_level / 2)
        out$nat_inf$tmle$reject$mult <- (abs(out$nat_inf$tmle$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                               out$nat_inf$tmle$boot_se$se_log_mult) > qnorm(1 - alpha_level / 2)
      }
      
      class(out$nat_inf$tmle) <- "tmle"

    }
    
    if("sens" %in% method){
      if(ml){
        out$nat_inf$sens$pt_est <- do_sens_aipw_nat_inf(data = data, models = ml_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, epsilon = epsilon, return_se = return_se)
      } else{
        out$nat_inf$sens$pt_est <- do_sens_aipw_nat_inf(data = data, models = models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, epsilon = epsilon, return_se = return_se)
      }
      
      if(is.null(out$nat_inf$sens$boot_se)){
        # closed form SE
        
        # test stats
        test_stat_add <- (out$nat_inf$sens$pt_est$additive_effect - null_hypothesis_value) /
          out$nat_inf$sens$pt_est$additive_se
        
        test_stat_mult <- (out$nat_inf$sens$pt_est$log_multiplicative_effect - null_hypothesis_value) /
          out$nat_inf$sens$pt_est$log_multiplicative_se
        
        # p-values
        p_val_add <- 2 * (1 - pnorm(abs(test_stat_add)))
        p_val_mult <- 2 * (1 - pnorm(abs(test_stat_mult)))
        
        # results data frames
        out$nat_inf$sens$reject$additive <- data.frame(
          epsilon   = out$nat_inf$sens$pt_est$epsilon,
          pt_est    = out$nat_inf$sens$pt_est$additive_effect,
          se        = out$nat_inf$sens$pt_est$additive_se,
          test_stat = test_stat_add,
          p_val      = p_val_add,
          reject    = p_val_add < alpha_level
        )
        
        out$nat_inf$sens$reject$mult <- data.frame(
          epsilon   = out$nat_inf$sens$pt_est$epsilon,
          pt_est    = exp(out$nat_inf$sens$pt_est$log_multiplicative_effect),
          se        = out$nat_inf$sens$pt_est$log_multiplicative_se,
          test_stat = test_stat_mult,
          p_val      = p_val_mult,
          reject    = p_val_mult < alpha_level
        )
        
      } else{
        # bootstrap se
        
        # additive
        test_stat_add <- (out$nat_inf$sens$pt_est$additive_effect - null_hypothesis_value) /
          out$nat_inf$sens$boot_se$se_additive
        p_val_add <- 2 * (1 - pnorm(abs(test_stat_add)))
        
        out$nat_inf$sens$reject$additive <- data.frame(
          epsilon = out$nat_inf$sens$boot_se$epsilon,
          pt_est = out$nat_inf$sens$pt_est$additive_effect,
          se = out$nat_inf$sens$boot_se$se_additive,
          test_stat = test_stat_add,
          p_val = p_val_add,
          reject = p_val_add < alpha_level
        )
        
        # multiplicative
        test_stat_mult <- (out$nat_inf$sens$pt_est$log_multiplicative_effect - null_hypothesis_value) /
          out$nat_inf$sens$boot_se$se_mult
        p_val_mult <- 2 * (1 - pnorm(abs(test_stat$mult)))
        
        out$nat_inf$sens$reject$mult <- data.frame(
          epsilon = out$nat_inf$sens$boot_se$epsilon,
          pt_est = exp(out$nat_inf$sens$pt_est$log_multiplicative_effect),
          se = out$nat_inf$sens$boot_se$se_mult,
          test_stat = test_stat_mult,
          p_val = p_val_mult,
          reject = p_val_mult < alpha_level
        )
        
      }
      
      class(out$nat_inf$sens) <- "sens"

    }
    
    if("bound" %in% method){
      
      out$nat_inf$bound$pt_est <- get_bound_nat_inf(data = data, Y_name = Y_name, Z_name = Z_name, S_name = S_name, family = family)
      
      # Bounds test - one sided
      # If effect direction < 0, test upper bound; Else, test lower bound 
      if(effect_dir == "negative"){
        out$nat_inf$bound$test_stat$additive <- out$nat_inf$bound$pt_est['additive_effect_upper'] / out$nat_inf$bound$boot_se$se_additive_upper
        out$nat_inf$bound$p_val$additive <- pnorm(out$nat_inf$bound$test_stat$additive, lower.tail = TRUE) 
        out$nat_inf$bound$reject$additive <- out$nat_inf$bound$test_stat$additive < qnorm(alpha_level)
      
        out$nat_inf$bound$test_stat$mult <- log(out$nat_inf$bound$pt_est['mult_effect_upper']) / out$nat_inf$bound$boot_se$se_log_mult_upper
        out$nat_inf$bound$p_val$mult <- pnorm(out$nat_inf$bound$test_stat$mult, lower.tail = TRUE) # had note about * 2 but i think we want one-sided right? 
        out$nat_inf$bound$reject$mult <- out$nat_inf$bound$test_stat$mult < qnorm(alpha_level)
      } else{
        out$nat_inf$bound$test_stat$additive <- out$nat_inf$bound$pt_est['additive_effect_lower'] / out$nat_inf$bound$boot_se$se_additive_lower
        out$nat_inf$bound$p_val$additive <- pnorm(out$nat_inf$bound$test_stat$additive, lower.tail = FALSE)
        out$nat_inf$bound$reject$additive <- out$nat_inf$bound$test_stat$additive > qnorm(1-alpha_level)
        
        out$nat_inf$bound$test_stat$mult <- log(out$nat_inf$bound$pt_est['mult_effect_lower']) / out$nat_inf$bound$boot_se$se_log_mult_lower
        out$nat_inf$bound$p_val$mult <- pnorm(out$nat_inf$bound$test_stat$mult, lower.tail = FALSE) # had note about * 2 but i think we want one-sided right? 
        out$nat_inf$bound$reject$mult <- out$nat_inf$bound$test_stat$mult > qnorm(1-alpha_level)
        
      }
      
      # Permutation test
      if(permutation){
        out$nat_inf$bound$permutation <- permutation_bound_nat_inf(data = data, Y_name = Y_name, Z_name = Z_name, S_name = S_name, n_permutations = n_perm, family = family, effect_dir = effect_dir)
      }
      
      class(out$nat_inf$bound) <- "bound"
      
    }
    
    class(out$nat_inf) <- "nat_inf"
    
  }
  
  # Doomed --------------------------------------------------------------------
  
  if("doomed" %in% estimand){
    
    if("gcomp" %in% method){
      out$doomed$gcomp$pt_est <- do_gcomp_doomed(data = data, models = models)
      
      out$doomed$gcomp$test_stat$additive <- (out$doomed$gcomp$pt_est['additive_effect'] - null_hypothesis_value) / 
        out$doomed$gcomp$boot_se$se_additive
      
      out$doomed$gcomp$test_stat$mult <- (out$doomed$gcomp$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
        out$doomed$gcomp$boot_se$se_log_mult
      
      out$doomed$gcomp$p_val$additive <- 2 * (1 - pnorm(abs(out$doomed$gcomp$test_stat$additive)))
      out$doomed$gcomp$p_val$mult <- 2 * (1 - pnorm(abs(out$doomed$gcomp$test_stat$mult)))
      
      out$doomed$gcomp$reject$additive <- (abs(out$doomed$gcomp$pt_est['additive_effect'] - null_hypothesis_value) / 
                                              out$doomed$gcomp$boot_se$se_additive) > qnorm(1 - alpha_level/2)
      out$doomed$gcomp$reject$mult <- (abs(out$doomed$gcomp$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                          out$doomed$gcomp$boot_se$se_log_mult) > qnorm(1 - alpha_level/2)
      
      class(out$doomed$gcomp) <- "gcomp"
    }
    
    if("ipw" %in% method){
      out$doomed$ipw$pt_est <- do_ipw_doomed(data = data, models = models, S_name = S_name, Y_name = Y_name, Z_name = Z_name)
      
      out$doomed$ipw$test_stat$additive <- (out$doomed$ipw$pt_est['additive_effect'] - null_hypothesis_value) / 
        out$doomed$ipw$boot_se$se_additive
      
      out$doomed$ipw$test_stat$mult <- (out$doomed$ipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
        out$doomed$ipw$boot_se$se_log_mult
      
      out$doomed$ipw$p_val$additive <- 2 * (1 - pnorm(abs(out$doomed$ipw$test_stat$additive)))
      out$doomed$ipw$p_val$mult <- 2 * (1 - pnorm(abs(out$doomed$ipw$test_stat$mult)))
      
      out$doomed$ipw$reject$additive <- (abs(out$doomed$ipw$pt_est['additive_effect'] - null_hypothesis_value) / 
                                            out$doomed$ipw$boot_se$se_additive) > qnorm(1 - alpha_level/2)
      out$doomed$ipw$reject$mult <- (abs(out$doomed$ipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                        out$doomed$ipw$boot_se$se_log_mult) > qnorm(1 - alpha_level/2)
      
      class(out$doomed$ipw) <- "ipw"
    }
    
    if("aipw" %in% method){
      if(ml){
        out$doomed$aipw$pt_est <- do_aipw_doomed(data = data, models = ml_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, return_se = return_se)
      } else{
        out$doomed$aipw$pt_est <- do_aipw_doomed(data = data, models = models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, return_se = return_se)
      }
      
      if(is.null(out$doomed$aipw$boot_se)){
        # closed form SE
        
        out$doomed$aipw$test_stat$additive <- (out$doomed$aipw$pt_est['additive_effect'] - null_hypothesis_value) / 
          out$doomed$aipw$pt_est['additive_se']
        
        out$doomed$aipw$test_stat$mult <- (out$doomed$aipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
          out$doomed$aipw$pt_est['log_multiplicative_se']
        
        out$doomed$aipw$p_val$additive <- 2 * (1 - pnorm(abs(out$doomed$aipw$test_stat$additive)))
        out$doomed$aipw$p_val$mult <- 2 * (1 - pnorm(abs(out$doomed$aipw$test_stat$mult)))
        
        out$doomed$aipw$reject$additive <- (abs(out$doomed$aipw$pt_est['additive_effect'] - null_hypothesis_value) / 
                                               out$doomed$aipw$pt_est['additive_se']) > qnorm(1 - alpha_level / 2)
        out$doomed$aipw$reject$mult <- (abs(out$doomed$aipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                           out$doomed$aipw$pt_est['log_multiplicative_se']) > qnorm(1 - alpha_level / 2)
      } else{
        # bootstrap se
        
        out$doomed$aipw$test_stat$additive <- (out$doomed$aipw$pt_est['additive_effect'] - null_hypothesis_value) / 
          out$doomed$aipw$boot_se$se_additive
        
        out$doomed$aipw$test_stat$mult <- (out$doomed$aipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
          out$doomed$aipw$boot_se$se_log_mult
        
        out$doomed$aipw$p_val$additive <- 2 * (1 - pnorm(abs(out$doomed$aipw$test_stat$additive)))
        out$doomed$aipw$p_val$mult <- 2 * (1 - pnorm(abs(out$doomed$aipw$test_stat$mult)))
        
        out$doomed$aipw$reject$additive <- (abs(out$doomed$aipw$pt_est['additive_effect'] - null_hypothesis_value) / 
                                               out$doomed$aipw$boot_se$se_additive) > qnorm(1 - alpha_level / 2)
        out$doomed$aipw$reject$mult <- (abs(out$doomed$aipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                           out$doomed$aipw$boot_se$se_log_mult) > qnorm(1 - alpha_level / 2)
        
      }
      
      class(out$doomed$aipw) <- "aipw"
      
    }
    
    if("bound" %in% method){
      
      out$doomed$bound$pt_est <- get_bound_doomed(data = data, Y_name = Y_name, Z_name = Z_name, S_name = S_name, family = family)
      
      # Bounds test - one sided
      # If effect direction < 0, test upper bound; Else, test lower bound 
      if(effect_dir == "negative"){
        out$doomed$bound$test_stat$additive <- out$doomed$bound$pt_est['additive_effect_upper'] / out$doomed$bound$boot_se$se_additive_upper
        out$doomed$bound$p_val$additive <- pnorm(out$doomed$bound$test_stat$additive, lower.tail = TRUE) 
        out$doomed$bound$reject$additive <- out$doomed$bound$test_stat$additive < qnorm(alpha_level)
        
        out$doomed$bound$test_stat$mult <- log(out$doomed$bound$pt_est['mult_effect_upper']) / out$doomed$bound$boot_se$se_log_mult_upper
        out$doomed$bound$p_val$mult <- pnorm(out$doomed$bound$test_stat$mult, lower.tail = TRUE) 
        out$doomed$bound$reject$mult <- out$doomed$bound$test_stat$mult < qnorm(alpha_level)
        
      } else{
        out$doomed$bound$test_stat$additive <- out$doomed$bound$pt_est['additive_effect_lower'] / out$doomed$bound$boot_se$se_additive_lower
        out$doomed$bound$p_val$additive <- pnorm(out$doomed$bound$test_stat$additive, lower.tail = FALSE)
        out$doomed$bound$reject$additive <- out$doomed$bound$test_stat$additive > qnorm(1-alpha_level)
        
        out$doomed$bound$test_stat$mult <- log(out$doomed$bound$pt_est['mult_effect_lower']) / out$doomed$bound$boot_se$se_log_mult_lower
        out$doomed$bound$p_val$mult <- pnorm(out$doomed$bound$test_stat$mult, lower.tail = FALSE) 
        out$doomed$bound$reject$mult <- out$doomed$bound$test_stat$mult > qnorm(1-alpha_level)
      }
      
      # Permutation test
      if(permutation){
        out$doomed$bound$permutation <- permutation_bound_doomed(data = data, Y_name = Y_name, Z_name = Z_name, S_name = S_name, n_permutations = n_perm, family = family, effect_dir = effect_dir)
      }
      
      class(out$doomed$bound) <- "bound"
      
    }
    
    if(!is.null(out$doomed)){
      class(out$doomed) <- "doomed"
    }

  }
  
  # Population ----------------------------------------------------------------
  
  if("pop" %in% estimand){
    
    if("gcomp" %in% method){
      out$pop$gcomp$pt_est <- do_gcomp_pop(data = data, models = models,  Z_name = Z_name, X_name = X_name, two_part_model = two_part_model)
      
      out$pop$gcomp$test_stat$additive <- (out$pop$gcomp$pt_est['additive_effect'] - null_hypothesis_value) / 
        out$pop$gcomp$boot_se$se_additive
      
      out$pop$gcomp$test_stat$mult <- (out$pop$gcomp$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
        out$pop$gcomp$boot_se$se_log_mult
      
      out$pop$gcomp$p_val$additive <- 2 * (1 - pnorm(abs(out$pop$gcomp$test_stat$additive)))
      out$pop$gcomp$p_val$mult <- 2 * (1 - pnorm(abs(out$pop$gcomp$test_stat$mult)))
      
      out$pop$gcomp$reject$additive <- (abs(out$pop$gcomp$pt_est['additive_effect'] - null_hypothesis_value) / 
                                              out$pop$gcomp$boot_se$se_additive) > qnorm(1 - alpha_level/2)
      out$pop$gcomp$reject$mult <- (abs(out$pop$gcomp$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                          out$pop$gcomp$boot_se$se_log_mult) > qnorm(1 - alpha_level/2)
      
      class(out$pop$gcomp) <- "gcomp"
    }
    
    if("ipw" %in% method){
      out$pop$ipw$pt_est <- do_ipw_pop(data = data, models = models, Y_name = Y_name, Z_name = Z_name)
      
      out$pop$ipw$test_stat$additive <- (out$pop$ipw$pt_est['additive_effect'] - null_hypothesis_value) / 
        out$pop$ipw$boot_se$se_additive
      
      out$pop$ipw$test_stat$mult <- (out$pop$ipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
        out$pop$ipw$boot_se$se_log_mult
      
      out$pop$ipw$p_val$additive <- 2 * (1 - pnorm(abs(out$pop$ipw$test_stat$additive)))
      out$pop$ipw$p_val$mult <- 2 * (1 - pnorm(abs(out$pop$ipw$test_stat$mult)))
      
      out$pop$ipw$reject$additive <- (abs(out$pop$ipw$pt_est['additive_effect'] - null_hypothesis_value) / 
                                            out$pop$ipw$boot_se$se_additive) > qnorm(1 - alpha_level/2)
      out$pop$ipw$reject$mult <- (abs(out$pop$ipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                        out$pop$ipw$boot_se$se_log_mult) > qnorm(1 - alpha_level/2)
      
      class(out$pop$ipw) <- "ipw"
    }
    
    if("aipw" %in% method){
      if(ml){
        out$pop$aipw$pt_est <- do_aipw_pop(data = data, models = ml_models, Y_name = Y_name, Z_name = Z_name, X_name = X_name, return_se = return_se, two_part_model = two_part_model)
      } else{
        out$pop$aipw$pt_est <- do_aipw_pop(data = data, models = models, Y_name = Y_name, Z_name = Z_name, X_name = X_name, return_se = return_se, two_part_model = two_part_model)
      }
      
      if(is.null(out$pop$aipw$boot_se)){
        # closed form SE
        
        out$pop$aipw$test_stat$additive <- (out$pop$aipw$pt_est['additive_effect'] - null_hypothesis_value) / 
          out$pop$aipw$pt_est['additive_se']
        
        out$pop$aipw$test_stat$mult <- (out$pop$aipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
          out$pop$aipw$pt_est['log_multiplicative_se']
        
        out$pop$aipw$p_val$additive <- 2 * (1 - pnorm(abs(out$pop$aipw$test_stat$additive)))
        out$pop$aipw$p_val$mult <- 2 * (1 - pnorm(abs(out$pop$aipw$test_stat$mult)))
        
        out$pop$aipw$reject$additive <- (abs(out$pop$aipw$pt_est['additive_effect'] - null_hypothesis_value) / 
                                               out$pop$aipw$pt_est['additive_se']) > qnorm(1 - alpha_level / 2)
        out$pop$aipw$reject$mult <- (abs(out$pop$aipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                           out$pop$aipw$pt_est['log_multiplicative_se']) > qnorm(1 - alpha_level / 2)
      } else{
        # bootstrap se
        
        out$pop$aipw$test_stat$additive <- (out$pop$aipw$pt_est['additive_effect'] - null_hypothesis_value) / 
          out$pop$aipw$boot_se$se_additive
        
        out$pop$aipw$test_stat$mult <- (out$pop$aipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
          out$pop$aipw$boot_se$se_log_mult
        
        out$pop$aipw$p_val$additive <- 2 * (1 - pnorm(abs(out$pop$aipw$test_stat$additive)))
        out$pop$aipw$p_val$mult <- 2 * (1 - pnorm(abs(out$pop$aipw$test_stat$mult)))
        
        out$pop$aipw$reject$additive <- (abs(out$pop$aipw$pt_est['additive_effect'] - null_hypothesis_value) / 
                                               out$pop$aipw$boot_se$se_additive) > qnorm(1 - alpha_level / 2)
        out$pop$aipw$reject$mult <- (abs(out$pop$aipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                           out$pop$aipw$boot_se$se_log_mult) > qnorm(1 - alpha_level / 2)
        
      }
      
      class(out$pop$aipw) <- "aipw"
      
    }
    
    if(!is.null(out$pop)){
      class(out$pop) <- "pop"
    }
    
  }
  
  if(return_models){
    out$models <- model_list
  }
 
  class(out) <- "vegrowth"
  
  return(out)

}
