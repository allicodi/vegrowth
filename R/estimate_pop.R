#' Function for g-computation of average counterfactual post-infection outcome marginally
#' 
#' @param data dataset to predict on
#' @param models list of pre-fit models needed for estimation
#' @param V_name name of vaccine treatment variable, default V
#' @param X_name character vector containing name(s) of covariates, default X
#' 
#' @returns g-comp estimate of growth effect for population estimand
do_gcomp_pop <- function(data, 
                          models,
                          V_name = "V",
                          X_name = c("X")){
  # E[G(1) - G(0)] = E[G(1)] - E[G(0)] = E[E[G | V = 1, X = x]] - E[E[G | V = 0, X = x]]
  
  #df_V1 <- data.frame(V = 1, X = data[[X_name]])
  df_V1 <- data.frame(V = 1, X = data[,colnames(data) %in% X_name, drop = FALSE])
  names(df_V1) <- c(V_name, X_name)
  
  if(inherits(models$fit_G_V_X, "SuperLearner")){
    Qbar_V1 <- predict(models$fit_G_V_X, newdata = df_V1, type = "response")$pred
  } else{
    Qbar_V1 <- predict(models$fit_G_V_X, newdata = df_V1, type = "response")
  }
  
  psi_1 <- mean(Qbar_V1)
  
  df_V0 <- data.frame(V = 0, X = data[,colnames(data) %in% X_name, drop = FALSE])
  names(df_V0) <- c(V_name, X_name)
  
  if(inherits(models$fit_G_V_X, "SuperLearner")){
    Qbar_V0 <- predict(models$fit_G_V_X, newdata = df_V0, type = "response")$pred 
  } else{
    Qbar_V0 <- predict(models$fit_G_V_X, newdata = df_V0, type = "response") 
  }
  
  psi_0 <- mean(Qbar_V0) 
  
  pop_growth_effect <- psi_1 - psi_0
  
  pop_growth_effect_log_mult <- log(psi_1 / psi_0)
  
  out <- c(pop_growth_effect, pop_growth_effect_log_mult)
  names(out) <-  c("additive_effect","log_multiplicative_effect")
  return(out)
}

#' Function for IPW of average counterfactual post-infection outcome marginally
#' 
#' @param data dataset to predict on
#' @param models list of pre-fit models needed for estimation
#' @param V_name name of vaccine treatment variable, default V
#' @param G_name character vector containing name(s) of post-infection outcome, default G
#' 
#' @returns g-comp estimate of growth effect for population estimand
do_ipw_pop <- function(data, models, V_name = "V", X_name = "G"){
  
  
  if(inherits(models$fit_G_V_X, "SuperLearner")){
    pi_1_X <- predict(models$fit_V_X, newdata = data)$pred
  } else{
    pi_1_X <- models$fit_V_X$fitted.values
  }
  pi_0_X <- 1 - pi_1_X
  
  Y <- data[[G_name]]
  Z <- data[[V_name]]
  
  psi_1_ipw <- mean(
    Z / pi_1_X * Y
  ) 
  psi_0_ipw <- mean(
    (1 - Z) / pi_0_X * Y
  ) 
  
  pop_growth_effect <- psi_1_ipw - psi_0_ipw
  pop_growth_effect_log_mult <- log(psi_1_ipw / psi_0_ipw)
  
  out <- c(pop_growth_effect, pop_growth_effect_log_mult)
  names(out) <-  c("additive_effect","log_multiplicative_effect")
  return(out)
}

#' Function for AIPW of average counterfactual post-infection outcome marginally
#' 
#' @param data dataset to predict on
#' @param models list of pre-fit models needed for estimation
#' @param V_name name of vaccine treatment variable, default V
#' @param X_name character vector containing name(s) of covariates, default X
#' @param Y_name character vector containing name(s) of covariates, default Y
#' 
#' @returns g-comp estimate of growth effect for population estimand
do_efficient_aipw_pop <- function(
    data, 
    models,
    V_name = "V",
    X_name = c("X"),
    G_name = "G",
    return_se = TRUE
){
  
  df_V1 <- data.frame(V = 1, X = data[,colnames(data) %in% X_name, drop = FALSE])
  names(df_V1) <- c(V_name, X_name)
  
  if(inherits(models$fit_G_V_X, "SuperLearner")){
    Qbar_V1 <- predict(models$fit_G_V_X, newdata = df_V1)$pred
    pi_1_X <- predict(models$fit_V_X, newdata = data)$pred
  } else{
    Qbar_V1 <- predict(models$fit_G_V_X, newdata = df_V1, type = "response")
    pi_1_X <- models$fit_V_X$fitted.values
  }
  pi_0_X <- 1 - pi_1_X
  
  Y <- data[[G_name]]
  Z <- data[[V_name]]
  
  psi_1_plugin <- mean(Qbar_V1)
  augmentation_1 <- Z / pi_1_X * ( Y - Qbar_V1 ) + Qbar_V1 - psi_1
  psi_1_aipw <- psi_1_plugin + mean(eif_psi_1)
  
  
  df_V0 <- data.frame(V = 0, X = data[,colnames(data) %in% X_name, drop = FALSE])
  names(df_V0) <- c(V_name, X_name)
  
  if(inherits(models$fit_G_V_X, "SuperLearner")){
    Qbar_V0 <- predict(models$fit_G_V_X, newdata = df_V0, type = "response")$pred 
  } else{
    Qbar_V0 <- predict(models$fit_G_V_X, newdata = df_V0, type = "response") 
  }
  
  psi_0_plugin <- mean(Qbar_V0)
  augmentation_0 <- (1 - Z) / pi_0_X * ( Y - Qbar_V0 ) + Qbar_V0 - psi_0
  psi_0_aipw <- psi_0_plugin + mean(eif_psi_0)
  
  # Additive effect
  n <- dim(data)[1]
  efficient_growth_effect <- psi_1_aipw - psi_0_aipw
  se <- sqrt(var(augmentation_1 - augmentation_0) / n)
  
  # Multiplicative effect (log scale)
  efficient_growth_effect_log_mult <- log(psi_1_aipw / psi_0_aipw)
  
  # Get SE using IF matrix same way as TMLE
  if_matrix <- cbind(augmentation_1, augmentation_0)
  cov_matrix <- cov(if_matrix) / n
  
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
