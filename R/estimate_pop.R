#' Function for g-computation of average counterfactual post-infection outcome marginally
#' 
#' @param data dataset to predict on
#' @param models list of pre-fit models needed for estimation
#' @param Z_name name of vaccine treatment variable, default Z
#' @param X_name character vector containing name(s) of covariates, default X
#' 
#' @export
#' 
#' @returns g-comp estimate of growth effect for population estimand
do_gcomp_pop <- function(data, 
                          models,
                          two_part_model = FALSE,
                          Z_name = "Z",
                          X_name = c("X")){
  # E[Y(1) - Y(0)] = E[Y(1)] - E[Y(0)] = E[E[Y | Z = 1, X = x]] - E[E[Y | Z = 0, X = x]]
  
  df_Z1 <- data.frame(Z = 1, X = data[,colnames(data) %in% X_name, drop = FALSE])
  names(df_Z1) <- c(Z_name, X_name)
  df_Z0 <- data.frame(Z = 0, X = data[,colnames(data) %in% X_name, drop = FALSE])
  names(df_Z0) <- c(Z_name, X_name)
  
  if(!two_part_model){
    E_Y_Z1_X <- simple_predict(models$fit_Y_Z_X, newdata = df_Z1)
    E_Y_Z0_X <- simple_predict(models$fit_Y_Z_X, newdata = df_Z0)
  }else{
    E_Y_Z0_S0_X <- simple_predict(models$fit_Y_Z0_S0_X, newdata = data)
    E_Y_Z0_S1_X <- simple_predict(models$fit_Y_Z0_S1_X, newdata = data)

    E_Y_Z1_S0_X <- simple_predict(models$fit_Y_Z1_S0_X, newdata = data)
    E_Y_Z1_S1_X <- simple_predict(models$fit_Y_Z1_S1_X, newdata = data)
    
    rho_0_X <- simple_predict(models$fit_S_Z0_X, newdata = data)
    rho_1_X <- simple_predict(models$fit_S_Z1_X, newdata = data)

    E_Y_Z1_X <- E_Y_Z1_S1_X * rho_1_X + E_Y_Z1_S0_X * (1 - rho_1_X)
    E_Y_Z0_X <- E_Y_Z0_S1_X * rho_0_X + E_Y_Z0_S0_X * (1 - rho_0_X)
  }

  psi_1 <- mean(E_Y_Z1_X)
  psi_0 <- mean(E_Y_Z0_X)

  pop_growth_effect <- psi_1 - psi_0
  
  pop_growth_effect_log_mult <- log(psi_1 / psi_0)
  
  out <- c(pop_growth_effect, pop_growth_effect_log_mult, psi_1, psi_0)
  names(out) <-  c("additive_effect","log_multiplicative_effect", "psi_1", "psi_0")
  return(out)
}

#' Function for IPW of average counterfactual post-infection outcome marginally
#' 
#' @param data dataset to predict on
#' @param models list of pre-fit models needed for estimation
#' @param Z_name name of vaccine treatment variable, default Z
#' @param Y_name character vector containing name(s) of post-infection outcome, default Y
#' 
#' @returns g-comp estimate of growth effect for population estimand
do_ipw_pop <- function(data, models, Z_name = "Z", Y_name = "Y"){
  
  pi_1_X <- simple_predict(models$fit_Z_X, newdata = data)
  pi_0_X <- 1 - pi_1_X
  
  Y <- data[[Y_name]]
  Z <- data[[Z_name]]
  
  psi_1_ipw <- mean(
    Z / pi_1_X * Y
  ) 
  psi_0_ipw <- mean(
    (1 - Z) / pi_0_X * Y
  ) 
  
  pop_growth_effect <- psi_1_ipw - psi_0_ipw
  pop_growth_effect_log_mult <- log(psi_1_ipw / psi_0_ipw)
  
  out <- c(pop_growth_effect, pop_growth_effect_log_mult, psi_1_ipw, psi_0_ipw)
  names(out) <-  c("additive_effect","log_multiplicative_effect", "psi_1", "psi_0")
  return(out)
}

#' Function for AIPW of average counterfactual post-infection outcome marginally
#' 
#' @param data dataset to predict on
#' @param models list of pre-fit models needed for estimation
#' @param Z_name name of vaccine treatment variable, default Z
#' @param X_name character vector containing name(s) of covariates, default X
#' @param S_name character vector containing name(s) of covariates, default Y
#' 
#' @returns g-comp estimate of growth effect for population estimand
do_aipw_pop <- function(
    data, 
    models,
    Z_name = "Z",
    X_name = c("X"),
    Y_name = "Y",
    return_se = TRUE
){
  
  df_Z1 <- data.frame(Z = 1, X = data[,colnames(data) %in% X_name, drop = FALSE])
  names(df_Z1) <- c(Z_name, X_name)
  
  
  Qbar_Z1 <- simple_predict(models$fit_Y_Z_X, newdata = df_Z1)
  pi_1_X <- simple_predict(models$fit_Z_X, newdata = data)
  pi_0_X <- 1 - pi_1_X
  
  Y <- data[[Y_name]]
  Z <- data[[Z_name]]
  
  psi_1_plugin <- mean(Qbar_Z1)
  augmentation_1 <- Z / pi_1_X * ( Y - Qbar_Z1 ) + Qbar_Z1 - psi_1_plugin
  psi_1_aipw <- psi_1_plugin + mean(augmentation_1)
  
  
  df_Z0 <- data.frame(Z = 0, X = data[,colnames(data) %in% X_name, drop = FALSE])
  names(df_Z0) <- c(Z_name, X_name)
  

  Qbar_Z0 <- simple_predict(models$fit_Y_Z_X, newdata = df_Z0)
  
  psi_0_plugin <- mean(Qbar_Z0)
  augmentation_0 <- (1 - Z) / pi_0_X * ( Y - Qbar_Z0 ) + Qbar_Z0 - psi_0_plugin
  psi_0_aipw <- psi_0_plugin + mean(augmentation_0)
  
  # Additive effect
  n <- dim(data)[1]
  growth_effect <- psi_1_aipw - psi_0_aipw
  se <- sqrt(var(augmentation_1 - augmentation_0) / n)
  
  se_psi_1 <- sqrt(var(augmentation_1) / n)
  se_psi_0 <- sqrt(var(augmentation_0) / n)
  
  # Multiplicative effect (log scale)
  growth_effect_log_mult <- log(psi_1_aipw / psi_0_aipw)
  
  # Yet SE using IF matrix same way as TMLE
  if_matrix <- cbind(augmentation_1, augmentation_0)
  cov_matrix <- cov(if_matrix) / n
  
  gradient <- matrix(c(1 / psi_1_aipw, -1 / psi_0_aipw), ncol = 1)
  
  se_log_mult_eff <- sqrt(t(gradient) %*% cov_matrix %*% gradient)
  
  if(return_se){
    out <- c(growth_effect, se, growth_effect_log_mult, se_log_mult_eff, psi_1_aipw, se_psi_1, psi_0_aipw, se_psi_0)
    names(out) <- c("additive_effect", "additive_se", "log_multiplicative_effect", "log_multiplicative_se", "psi_1", "se_psi_1", "psi_0", "se_psi_0")
    return(out)
  }else{
    out <- c(growth_effect, growth_effect_log_mult)
    names(out) <- c("additive_effect", "log_multiplicative_effect")
    return(out)
  }
  
}
