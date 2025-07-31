#' Function for g-computation of counterfactual post-infection outcomes in the 
#' doomed principal strata
#' 
#' @param data dataset to predict on
#' @param models list of pre-fit models needed for estimation
#' 
#' @returns g-comp estimate of growth effect in the doomed strata
do_gcomp_doomed <- function(data, models){
  
  # Psi_1 = E[P(Y=1 | V = 0, X) / P(Y = 1 | V = 0) * E[G | V=1, X] ]
  if(inherits(models$fit_G_V1_Y1_X, "SuperLearner")){
    mu_11_X <- predict(models$fit_G_V1_Y1_X, newdata = data)$pred
    mu_10_X <- predict(models$fit_G_V0_Y1_X, newdata = data)$pred
    rho_1_X <- predict(models$fit_Y_V1_X, newdata = data)$pred
  } else{
    mu_11_X <- predict(models$fit_G_V1_Y1_X, newdata = data, type = "response")
    mu_10_X <- predict(models$fit_G_V0_Y1_X, newdata = data, type = "response")
    rho_1_X <- predict(models$fit_Y_V1_X, newdata = data, type = "response")
  }
  rho_bar_1 <- mean(rho_1_X)
  
  psi_1 <- mean(
    rho_1_X / rho_bar_1 * mu_11_X
  )
  psi_0 <- mean(
    rho_1_X / rho_bar_1 * mu_01_X
  )
  
  growth_effect <- psi_1 - psi_0
  growth_effect_log_mult <- log(psi_1 / psi_0)
  
  out <- c(growth_effect, growth_effect_log_mult)
  names(out) <- c("additive_effect","log_multiplicative_effect")
  
  return(out)
}

#' Function for IPW of counterfactual post-infection outcomes in the 
#' doomed principal stratum
#' 
#' @param data dataset to predict on
#' @param models list of pre-fit models needed for estimation
#' 
#' @returns IPW estimate of growth effect in the doomed principal stratum
do_ipw_doomed <- function(
    data, models,
    Y_name, G_name, V_name
){
  
  # Psi_1 = E[P(Y=1 | V = 0, X) / P(Y = 1 | V = 0) * E[G | V=1, X] ]
  if(inherits(models$fit_G_V1_Y1_X, "SuperLearner")){
    rho_1_X <- predict(models$fit_Y_V1_X, newdata = data)$pred
    rho_0_X <- predict(models$fit_Y_V0_X, newdata = data)$pred
    pi_1_X <- predict(models$fit_V_X, newdata = data)$pred
  } else{
    rho_1_X <- predict(models$fit_Y_V1_X, newdata = data, type = "response")
    rho_0_X <- predict(models$fit_Y_V0_X, newdata = data, type = "response")
    pi_1_X <- predict(models$fit_V_X, newdata = data, type = "response")
  }
  pi_0_X <- 1 - pi_1_X
  rho_bar_1 <- mean(rho_1_X)
  
  S <- data[[Y_name]]
  Y <- data[[G_name]]
  Z <- data[[V_name]]
  
  eta_1 <- mean(
    ( Z / pi_1_X ) * ( S / rho_bar_1 ) * Y
  )
  
  eta_0 <- mean(
    ( (1 - Z) / pi_0_X ) * ( S / rho_0_X ) * ( rho_1_X / rho_bar_1 ) * Y
  )
  
  growth_effect <- eta_1 - eta_0
  growth_effect_log_mult <- log(eta_1 / eta_0)
  
  out <- c(growth_effect, growth_effect_log_mult)
  names(out) <- c("additive_effect","log_multiplicative_effect")
  
  return(out)
}

#' Function for efficient AIPW estimator in the doomed stratum
#' 
#' @param data dataset to predict on
#' @param models list of pre-fit models needed for estimation
#' @param G_name name of growth outcome variable, default G
#' @param V_name name of vaccine treatment variable, default V
#' @param Y_name name of infection variable, default Y
#' @param return_se flag to return standard error, defualt FALSE
#' 
#' @returns AIPW estimate of growth effect in doomed stratum (+ standard error if return_se = TRUE)
do_efficient_aipw_doomed <- function(data, 
                                     models,
                                     G_name = "G",
                                     V_name = "V",
                                     Y_name = "Y",
                                     return_se = FALSE){
  
  if(inherits(models$fit_Y_V0_X, "SuperLearner")){
    mu_11_X <- predict(models$fit_G_V1_Y1_X, newdata = data)$pred
    mu_01_X <- predict(models$fit_G_V0_Y1_X, newdata = data)$pred
    
    rho_0_X <- predict(models$fit_Y_V0_X, newdata = data)$pred
    rho_1_X <- predict(models$fit_Y_V1_X, newdata = data)$pred
    
    pi_1_X <- predict(models$fit_V_X, newdata = data)$pred
  } else{
    mu_11_X <- predict(models$fit_G_V1_Y1_X, newdata = data, type = "response")
    mu_01_X <- predict(models$fit_G_V0_Y1_X, newdata = data, type = "response")
    
    rho_0_X <- predict(models$fit_Y_V0_X, newdata = data, type = "response")
    rho_1_X <- predict(models$fit_Y_V1_X, newdata = data, type = "response")
    
    pi_1_X <- models$fit_V_X$fitted.values
  }
  
  pi_0_X <- 1 - pi_1_X
  rho_bar_1 <- mean(rho_0)
  
  eta_tilde_0 <- rho_1_X * mu_11_X / rho_bar_1
  eta_tilde_1 <- rho_1_X * mu_01_X / rho_bar_1
  
  eta_0_plugin <- mean( eta_tilde_0 )
  eta_1_plugin <- mean( eta_tilde_1 )
  
  Z <- data[[V_name]]
  S <- data[[Y_name]]
  Y <- data[[G_name]]
  
  augmentation_0 <- (
    ( (1 - Z) / pi_0_X ) * ( rho_1_X / rho_bar_1 ) * ( S / rho_0_X ) * ( Y - mu_01_X ) + 
      (mu_01_X - eta_0_plugin) / rho_bar_1 * ( Z / pi_1_X ) * ( S - rho_1_X ) -
      eta_0_plugin / rho_bar_1 * ( rho_1_X - rho_bar_1 ) + eta_tilde_0 - eta_0_plugin
  )
  augmentation_1 <- (
    ( Z / pi_1_X ) * ( S / rho_bar_1 ) * ( Y - mu_11_X ) + 
      ( mu_11_X - eta_1_plugin ) / rho_bar_1 * ( Z / pi_1_X ) * ( S - rho_1_X ) - 
      eta_1_plugin / rho_bar_1 * ( rho_1_X - rho_bar_1 ) + eta_tilde_1 - eta_1_plugin
    
  )
  
  eta_0_aipw <- eta_0 + mean(augmentation_0)
  eta_1_aipw <- eta_1 + mean(augmentation_1)
  
  
  # Additive effect
  efficient_growth_effect <- eta_1_aipw - eta_0_aipw
  se <- sqrt(var(augmentation_1 - augmentation_0) / dim(data)[1])
  
  # Multiplicative effect (log scale)
  efficient_growth_effect_log_mult <- log(eta_1_aipw / eta_0_aipw)
  
  # Get SE using IF matrix same way as TMLE
  if_matrix <- cbind(augmentation_1, augmentation_0)
  cov_matrix <- cov(if_matrix) / dim(data)[1]
  
  gradient <- matrix(c(1 / psi_1_aipw, -1 / psi_0_aipw), ncol = 1)
  
  se_log_mult_eff <- sqrt(t(gradient) %*% cov_matrix %*% gradient)
  
  if(return_se){
    out <- c(efficient_growth_effect, se, efficient_growth_effect_log_mult, se_log_mult_eff)
    names(out) <- c("additive_effect", "additive_se", "log_multiplicative_effect", "log_multiplicative_se")
    return(out)
  }else{
    out <- c(efficient_growth_effect, efficient_growth_effect_log_mult)
    names(out) <- c("additive_effect", "log_multiplicative_effect")
    return(out)
  }
}