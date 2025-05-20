#' Function for g-computation of VE estimand
#' 
#' @param data dataset to predict on
#' @param models list of pre-fit models needed for estimation
#' 
#' @returns g-comp estimate of growth effect for VE estimand
do_gcomp <- function(data, models){
  
  # Psi_1 = E[P(Y=1 | V = 0, X) / P(Y = 1 | V = 0) * E[G | V=1, X] ]
  if(inherits(models$fit_G_V1_Y1_X, "SuperLearner")){
    E_G_V1_Y1_X <- predict(models$fit_G_V1_Y1_X, newdata = data, type = "response")$pred
    E_G_V1_Y0_X <- predict(models$fit_G_V1_Y0_X, newdata = data, type = "response")$pred
    P_Y1_V1_X <- predict(models$fit_Y_V1_X, newdata = data, type = "response")$pred
    P_Y1_V0_X <- predict(models$fit_Y_V0_X, newdata = data, type = "response")$pred
  } else{
    E_G_V1_Y1_X <- predict(models$fit_G_V1_Y1_X, newdata = data, type = "response")
    E_G_V1_Y0_X <- predict(models$fit_G_V1_Y0_X, newdata = data, type = "response")
    P_Y1_V1_X <- predict(models$fit_Y_V1_X, newdata = data, type = "response")
    P_Y1_V0_X <- predict(models$fit_Y_V0_X, newdata = data, type = "response")
  }
  
  P_Y1_V0 <- mean(P_Y1_V0_X)
  VE_X <- 1 - ( P_Y1_V1_X / P_Y1_V0_X )
  E_G1_Y01_X <- E_G_V1_Y1_X * (1 - VE_X) + E_G_V1_Y0_X * VE_X
  
  psi_1 <- mean(
    ( P_Y1_V0_X / P_Y1_V0 ) * E_G1_Y01_X
  )
  
  # Psi_0 = E[P(Y=1 | V = 0, X) / P(Y = 1 | V = 0) * E[G | V=0, Y = 1, X] ]
  
  # Option 1 for estimation:
  # psi_0 <- mean(sub_V0_Y1$G) 
  
  # Option 2 for estimation:
  if(inherits(models$fit_G_V0_Y1_X, "SuperLearner")){
    E_G_V0_Y1_X <- predict(models$fit_G_V0_Y1_X, newdata = data, type = "response")$pred
  } else {
    E_G_V0_Y1_X <- predict(models$fit_G_V0_Y1_X, newdata = data, type = "response")
  }
  
  psi_0 <- mean(
    ( P_Y1_V0_X / P_Y1_V0 ) * E_G_V0_Y1_X
  )

  growth_effect <- psi_1 - psi_0
  
  return(growth_effect)
}

#' Function for g-computation of traditional population estimand
#' 
#' @param data dataset to predict on
#' @param models list of pre-fit models needed for estimation
#' @param V_name name of vaccine treatment variable, default V
#' @param X_name character vector containing name(s) of covariates, default X
#' 
#' @returns g-comp estimate of growth effect for population estimand
do_gcomp_pop_estimand <- function(data, 
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
  
  return(pop_growth_effect)
}

#' Function for efficient AIPW estimator
#' 
#' @param data dataset to predict on
#' @param models list of pre-fit models needed for estimation
#' @param G_name name of growth outcome variable, default G
#' @param V_name name of vaccine treatment variable, default V
#' @param Y_name name of infection variable, default Y
#' @param return_se flag to return standard error, defualt FALSE
#' 
#' @returns AIPW estimate of growth effect (+ standard error if return_se = TRUE)
do_efficient_aipw <- function(data, 
                              models,
                              G_name = "G",
                              V_name = "V",
                              Y_name = "Y",
                              return_se = FALSE){
  
  # vaccine probabilities
  pi_1 <- mean(data[[V_name]])
  pi_0 <- 1 - pi_1
  
  # Get weight
  sub_V0 <- data[data[[V_name]] == 0,]
  
  rho_bar_0 <- mean(sub_V0[[Y_name]])
  
  if(inherits(models$fit_Y_V0_X, "SuperLearner")){
    rho_0 <- predict(models$fit_Y_V0_X, newdata = data, type = "response")$pred
    mu_01 <- predict(models$fit_G_V0_Y1_X, newdata = data, type = "response")$pred
  } else{
    rho_0 <- predict(models$fit_Y_V0_X, newdata = data, type = "response")
    mu_01 <- predict(models$fit_G_V0_Y1_X, newdata = data, type = "response")
  }
  
  
  # psi_0 = Weight * E[E[G | V = 0, Y = 1, X]]

  psi_tilde_0 <- rho_0 / rho_bar_0 * mu_01
  
  psi_0 <- mean( psi_tilde_0 )
  
  augmentation_0 <- (
    (1 - data[[V_name]]) / pi_0 * ( data[[Y_name]] / rho_bar_0 ) * (data[[G_name]] - mu_01) + 
      (1 - data[[V_name]]) / pi_0 * ( mu_01 - psi_0 ) / rho_bar_0 * ( data[[Y_name]] - rho_0 ) + 
      ( psi_0 / rho_bar_0 ) * ( rho_0 - rho_bar_0 ) + 
      psi_tilde_0 - psi_0
  )
  
  psi_0_aipw <- psi_0 + mean(augmentation_0)
  
  # psi_1 = Weight * E[E[G | V = 1, X]]
  
  if(inherits(models$fit_G_V1_Y1_X, "SuperLearner")){
    mu_11 <- predict(models$fit_G_V1_Y1_X, newdata = data, type = "response")$pred
    mu_10 <- predict(models$fit_G_V1_Y0_X, newdata = data, type = "response")$pred
    rho_1 <- predict(models$fit_Y_V1_X, newdata = data, type = "response")$pred
  } else{
    mu_11 <- predict(models$fit_G_V1_Y1_X, newdata = data, type = "response")
    mu_10 <- predict(models$fit_G_V1_Y0_X, newdata = data, type = "response")
    rho_1 <- predict(models$fit_Y_V1_X, newdata = data, type = "response")
  }
  
  psi_tilde_1 <- rho_1 / rho_bar_0 * mu_11 + ( rho_0 - rho_1 ) / rho_bar_0 * mu_10
  psi_1 <- mean( psi_tilde_1 )
  
  augmentation_1 <- (
    (data[[V_name]] / pi_1) * (data[[Y_name]] / rho_bar_0) * (data[[G_name]] - mu_11) +
      (data[[V_name]] / pi_1) * ((1 - data[[Y_name]]) / (1 - rho_1)) * (rho_0 - rho_1) / rho_bar_0 * (data[[G_name]] - mu_10) + 
      (data[[V_name]] / pi_1) * (mu_11 - mu_10) / rho_bar_0 * (data[[Y_name]] - rho_1) + 
      ((1 - data[[V_name]]) / pi_0) * (mu_10 - psi_1) / rho_bar_0 * (data[[Y_name]] - rho_0) - 
      psi_1 / rho_bar_0 * (rho_0 - rho_bar_0) + psi_tilde_1 - psi_1
  )
  
  psi_1_aipw <- psi_1 + mean(augmentation_1)
  
  efficient_growth_effect <- psi_1_aipw - psi_0_aipw
  # TODO: add in multiplicative effects
  # agumentation_1 = phi_1_data from tmle function
  # similarly for _0
  if(return_se){
    se <- sqrt(var(augmentation_1 - augmentation_0) / dim(data)[1])
    return(c(efficient_growth_effect, se))
  }else{
    return(efficient_growth_effect)
  }
}

#' Function for efficient TMLE estimator
#' 
#' @param data dataset to predict on
#' @param models list of pre-fit models needed for estimation
#' @param G_name name of growth outcome variable, default G
#' @param V_name name of vaccine treatment variable, default V
#' @param Y_name name of infection variable, default Y
#' @param return_se flag to return standard error, defualt FALSE
#' @param max_iter TODO 
#' @param tol TOOD 
#' 
#' @returns TMLE estimate of growth effect (+ standard error if return_se = TRUE)
do_efficient_tmle <- function(
    data, models, G_name = "G", V_name = "V", Y_name = "Y",
    return_se = FALSE, max_iter = 10,
    tol = 1 / (sqrt(dim(data)[1]) * log(dim(data)[1]))
){
  
  idx_V0 <- which(data[[V_name]] == 0)
  idx_V1 <- which(data[[V_name]] == 1)
  idx_V0_Y1 <- which(data[[V_name]] == 0 & data[[Y_name]] == 1)
  l <- min(data[[G_name]])
  u <- max(data[[G_name]])
  
  pi_1 <- rep(mean(data[[V_name]]), dim(data)[1])
  pi_0 <- 1 - pi_1
  
  if(inherits(models$fit_Y_V0_X, "SuperLearner")){
    rho_0 <- predict(models$fit_Y_V0_X, newdata = data, type = "response")$pred
    rho_1 <- predict(models$fit_Y_V1_X, newdata = data, type = "response")$pred
    
    mu_11 <- predict(models$fit_G_V1_Y1_X, newdata = data, type = "response")$pred
    mu_10 <- predict(models$fit_G_V1_Y0_X, newdata = data, type = "response")$pred
    mu_01 <- predict(models$fit_G_V0_Y1_X, newdata = data, type = "response")$pred
  } else{
    rho_0 <- predict(models$fit_Y_V0_X, newdata = data, type = "response")
    rho_1 <- predict(models$fit_Y_V1_X, newdata = data, type = "response")
    
    mu_11 <- predict(models$fit_G_V1_Y1_X, newdata = data, type = "response")
    mu_10 <- predict(models$fit_G_V1_Y0_X, newdata = data, type = "response")
    mu_01 <- predict(models$fit_G_V0_Y1_X, newdata = data, type = "response")
  }
  
  rho_bar_0 <- mean(rho_0)
  
  psi_tilde_0 <- rho_0 / rho_bar_0 * mu_01
  psi_0 <- mean( psi_tilde_0 )
  
  psi_tilde_1 <- rho_1 / rho_bar_0 * mu_11 + ( rho_0 - rho_1 ) / rho_bar_0 * mu_10
  psi_1 <- mean( psi_tilde_1 )
  
  phi_0 <- function(data, V_name, Y_name, G_name, pi_0, rho_0, rho_bar_0, mu_01, psi_tilde_0, psi_0) {
    (
      (1 - data[[V_name]]) / pi_0 * (data[[Y_name]] / rho_bar_0) * (data[[G_name]] - mu_01) +
        (1 - data[[V_name]]) / pi_0 * (mu_01 - psi_0) / rho_bar_0 * (data[[Y_name]] - rho_0) +
        (psi_0 / rho_bar_0) * (rho_0 - rho_bar_0) +
        psi_tilde_0 - psi_0
    )
  }
  
  phi_1 <- function(data, V_name, Y_name, G_name, pi_1, pi_0, rho_0, rho_bar_0, rho_1, mu_11, mu_10, psi_tilde_1, psi_1) {
    (
      (data[[V_name]] / pi_1) * (data[[Y_name]] / rho_bar_0) * (data[[G_name]] - mu_11) +
        (data[[V_name]] / pi_1) * ((1 - data[[Y_name]]) / (1 - rho_1)) * (rho_0 - rho_1) / rho_bar_0 * (data[[G_name]] - mu_10) +
        (data[[V_name]] / pi_1) * (mu_11 - mu_10) / rho_bar_0 * (data[[Y_name]] - rho_1) +
        ((1 - data[[V_name]]) / pi_0) * (mu_10 - psi_1) / rho_bar_0 * (data[[Y_name]] - rho_0) -
        psi_1 / rho_bar_0 * (rho_0 - rho_bar_0) +
        psi_tilde_1 - psi_1
    )
  }
  
  trim_logit <- function(p, tol = 1e-3){ 
    p[p < tol] <- tol
    p[p > 1 - tol] <- 1 - tol
    return(qlogis(p))
  }
  
  scale_01 <- function(x, l, u){
    ( x - l ) / ( u - l )
  }
  
  rescale_01 <- function(x, l, u){
    x * (u - l) + l
  }
  
  phi_0_data <- phi_0(data, V_name, Y_name, G_name, pi_0, rho_0, rho_bar_0, mu_01, psi_tilde_0, psi_0)
  phi_1_data <- phi_1(data, V_name, Y_name, G_name, pi_1, pi_0, rho_0, rho_bar_0, rho_1, mu_11, mu_10, psi_tilde_1, psi_1)
  phi_ge_data <- phi_1_data - phi_0_data
  
  mean_phi_0 <- mean(phi_0_data)
  mean_phi_1 <- mean(phi_1_data)
  mean_phi_ge_data <- mean(phi_ge_data)
  
  iter <- 0
  
  mu_11_star <- mu_11
  mu_10_star <- mu_10
  mu_01_star <- mu_01
  rho_1_star <- rho_1
  rho_0_star <- rho_0
  rho_bar_0_star <- rho_bar_0
  psi_tilde_0_star <- psi_tilde_0
  psi_tilde_1_star <- psi_tilde_1
  psi_0_star <- psi_0
  psi_1_star <- psi_1
  
  while((mean_phi_0^2 + mean_phi_1^2 + mean_phi_ge_data^2) > tol & iter <= max_iter){
    # cat("iter", iter, "\n")
    # cat("mean_eif", mean_phi_ge_data, "\n")
    
    # target mu's
    G_scale <- scale_01(data[[G_name]], l, u)
    
    # target mu_11
    mu_11_star_scale <- scale_01(mu_11_star, l, u)
    logit_mu_11_star_scale <- trim_logit(mu_11_star_scale)
    target_wt <- (
      (data[[V_name]] / pi_1) * (data[[Y_name]] / rho_bar_0_star)
    )
    
    target_data <- data.frame(
      G_scale = G_scale,
      target_wt = target_wt,
      logit_mu_11_star_scale = logit_mu_11_star_scale
    )
    target_fit <- suppressWarnings(glm(
      G_scale ~ offset(logit_mu_11_star_scale), 
      weight = target_wt,
      family = binomial(),
      data = target_data,
      start = c(0)
    ))
    mu_11_star <- rescale_01(target_fit$fitted.values, l, u)
    
    # target mu_01
    mu_01_star_scale <- scale_01(mu_01_star, l, u)
    logit_mu_01_star_scale <- trim_logit(mu_01_star_scale)
    target_wt <- (
      ((1 - data[[V_name]]) / (1 - pi_1)) * (data[[Y_name]] / rho_bar_0_star)
    )
    
    target_data <- data.frame(
      G_scale = G_scale,
      target_wt = target_wt,
      logit_mu_01_star_scale = logit_mu_01_star_scale
    )
    target_fit <- suppressWarnings(glm(
      G_scale ~ offset(logit_mu_01_star_scale), 
      weight = target_wt,
      family = binomial(),
      data = target_data,
      start = c(0)
    ))
    mu_01_star <- rescale_01(target_fit$fitted.values, l, u)
    
    # target mu_10
    mu_10_star_scale <- scale_01(mu_10_star, l, u)
    logit_mu_10_star_scale <- trim_logit(mu_10_star_scale)
    target_wt <- with(data, 
      ( data[[V_name]] / pi_1 ) * ( (1 - data[[Y_name]]) / rho_bar_0_star ) 
    )
    H1 <- ( rho_0_star - rho_1_star ) / ( 1 - rho_1_star )
    target_data <- data.frame(
      G_scale = G_scale,
      target_wt = target_wt,
      H1 = H1,
      logit_mu_10_star_scale = logit_mu_10_star_scale
    )
    target_fit <- suppressWarnings(glm(
      G_scale ~ -1 + offset(logit_mu_10_star_scale) + H1, 
      weight = target_wt,
      family = binomial(),
      data = target_data,
      start = c(0)
    ))
    mu_10_star <- rescale_01(target_fit$fitted.values, l, u)
    
    psi_tilde_1_star <- rho_1_star / rho_bar_0_star * mu_11_star + ( rho_0_star - rho_1_star ) / rho_bar_0_star * mu_10_star
    psi_1_star <- mean( psi_tilde_1_star )
    
    psi_tilde_0_star <- rho_0_star / rho_bar_0_star * mu_01_star
    psi_0_star <- mean( psi_tilde_0_star )
    
    # target rho_0
    H1 <- mu_10_star - psi_1_star
    H0 <- mu_01_star - psi_0_star
    logit_rho_0_star <- trim_logit(rho_0_star)
    target_wt <- (1 - data[[V_name]]) / pi_0
    
    # with linear models, these may be perfectly correlated, but numerically
    # R thinks they are not and tries to fit a glm, which blows up. setting 
    # H0 to a constant in these cases will remove the term from the model because
    # the model also includes an intercept
    if(cor(H1, H0) > 0.99999){
      H0 <- 1
    }
    
    target_data <- data.frame(
      Y_inf = data[[Y_name]],
      target_wt = target_wt,
      H1 = H1,
      H0 = H0,
      logit_rho_0_star = logit_rho_0_star
    )
    target_data <- setNames(target_data, c(Y_name, names(target_data[-1])))
    
    # include intercept so rho_bar_0_star is still mean(Y[V == 0])
    target_fit <- glm(
      as.formula(paste0(Y_name," ~ offset(logit_rho_0_star) + H1 + H0")), 
      family = binomial(),
      weight = target_wt,
      data = target_data,
      start = c(0, 0, 0)
    )
    rho_0_star <- target_fit$fitted.values
    
    # shouldn't change because of intercept, but just in case
    rho_bar_0_star <- mean(rho_0_star)
    
    ## sanity check
    # tmp <- with(data, 
    # ( (1 - V) / pi_0 ) * ( mu_01_star - psi_0_star ) / rho_bar_0_star * ( Y_inf - rho_0_star )
    # )
    # mean(tmp) # should be small
    # tmp <- with(data, 
    #     ( (1 - V) / pi_0 ) * ( (mu_10_star - psi_1_star) / rho_bar_0_star ) * (Y_inf - rho_0_star) 
    # )
    # mean(tmp) # should be small
    
    # target rho_1
    H1 <- mu_11_star - mu_10_star
    logit_rho_1_star <- trim_logit(rho_1_star)
    target_wt <- data[[V_name]] / pi_1
    target_data <- data.frame(
      Y_name = data[[Y_name]],
      target_wt = target_wt,
      H1 = H1,
      logit_rho_1_star = logit_rho_1_star
    )
    target_data <- setNames(target_data, c(Y_name, names(target_data[-1])))
    
    # include intercept so rho_bar_0_star is still mean(Y[V == 0])
    target_fit <- glm(
      as.formula(paste0(Y_name, " ~ -1 + offset(logit_rho_1_star) + H1")), 
      family = binomial(),
      weight = target_wt,
      data = target_data,
      start = c(0)
    )
    rho_1_star <- target_fit$fitted.values
    
    psi_tilde_1_star <- rho_1_star / rho_bar_0_star * mu_11_star + ( rho_0_star - rho_1_star ) / rho_bar_0_star * mu_10_star
    psi_1_star <- mean( psi_tilde_1_star )
    
    psi_tilde_0_star <- rho_0_star / rho_bar_0_star * mu_01_star
    psi_0_star <- mean( psi_tilde_0_star )
    
    phi_0_data <- phi_0(data, V_name, Y_name, G_name, pi_0, rho_0_star, rho_bar_0_star, mu_01_star, psi_tilde_0_star, psi_0_star)
    phi_1_data <- phi_1(data, V_name, Y_name, G_name, pi_1, pi_0, rho_0_star, rho_bar_0_star, rho_1_star, mu_11_star, mu_10_star, psi_tilde_1_star, psi_1_star)
    phi_ge_data <- phi_1_data - phi_0_data
    
    mean_phi_0 <- mean(phi_0_data)
    mean_phi_1 <- mean(phi_1_data)
    mean_phi_ge_data <- mean(phi_ge_data)
    
    iter <- iter + 1
  }
  
  tmle_ge <- psi_1_star - psi_0_star
  tmle_ge_log_mult <- log(psi_1_star / psi_0_star)
  
  if(return_se){
    se <- sqrt(var(phi_ge_data) / dim(data)[1])

    if_matrix <- cbind(phi_0_data, phi_1_data)
    cov_matrix <- cov(if_matrix) / dim(data)[1]
    # 1/psi_1, -1/psi_0
    gradient <- matrix(c(1 / psi_1_star, -1 / psi_0_star), ncol = 1)
    se_log_mult_eff <- sqrt(t(gradient) %*% cov_matrix %*% gradient)

    out <- c(tmle_ge, se, tmle_ge_log_mult, se_log_mult_eff)
    names(out) <- c("additive_effect", "additive_se", "log_multiplicative_effect", "log_multiplicative_se")
    return(out)
  }else{
    return(c(tmle_ge, tmle_ge_log_mult))
  }
  
}

#' Function to fit all models needed by various estimators
#' 
#' @param data dataset used to fit models
#' @param G_name name of growth outcome variable, default G
#' @param V_name name of vaccine treatment variable, default V
#' @param X_name character vector containing name(s) of covariates, default X
#' @param Y_name name of infection variable, default Y
#' @param est character vector of names of estimators to use for growth effect
#' @param G_X_model optional specify model to be used for fitting growth on covariates, otherwise growth on all covariates
#' @param Y_X_model optional specify model to be used for fitting infection on covariates, otherwise infection on all covariates
#' @param family family for outcome model, defaults to gaussian for growth
#' 
#' @export
#' 
#' @returns list with models needed for specified estimators
fit_models <- function(data,
                       G_name = "G",
                       V_name = "V",
                       X_name = c("X"),
                       Y_name = "Y",
                       est = c("gcomp_pop_estimand", 
                               "gcomp", 
                               "efficient_aipw", 
                               "efficient_tmle"),
                       G_V_X_model = NULL,
                       G_X_Y1_model = NULL,
                       G_X_Y0_model = NULL,
                       Y_X_model = NULL,
                       family = "gaussian"){
  
  # Prep model formulas if not pre-specified
  if(is.null(G_V_X_model)){
    G_V_X_model <- as.formula(paste0(G_name, "~",
                                   V_name, "+",
                                   paste0(X_name, collapse = "+")))
  }
  
  if(is.null(G_X_Y1_model)){
    G_X_Y1_model <- as.formula(paste0(G_name, "~",
                                   paste0(X_name, collapse = "+")))
  }
  
  if(is.null(G_X_Y0_model)){
    G_X_Y0_model <- as.formula(paste0(G_name, "~",
                                   paste0(X_name, collapse = "+")))
  }
  
  if(is.null(Y_X_model)){
    Y_X_model <- as.formula(paste0(Y_name, "~",
                                   paste0(X_name, collapse = "+")))
  }
  
  out <- list()
  
  # only needed for population estimator
  if("gcomp_pop_estimand" %in% est){
    out$fit_G_V_X <- glm(G_V_X_model, data = data, family = family)
  } 
  
  # needed for any weight-based estimator
  if(any(c("gcomp", "efficient_aipw", "efficient_tmle") %in% est)){
    sub_V0 <- data[data[[V_name]] == 0,]
    out$fit_Y_V0_X <- glm(Y_X_model, sub_V0, family = "binomial")
    
    sub_V1 <- data[data[[V_name]] == 1,]
    out$fit_Y_V1_X <- glm(Y_X_model, sub_V1, family = "binomial")
    
    sub_V1_Y1 <- data[data[[V_name]] == 1 & data[[Y_name]] == 1,]
    out$fit_G_V1_Y1_X <- glm(G_X_Y1_model, data = sub_V1_Y1, family = family)
    
    sub_V1_Y0 <- data[data[[V_name]] == 1 & data[[Y_name]] == 0,]
    out$fit_G_V1_Y0_X <- glm(G_X_Y0_model, data = sub_V1_Y0, family = family)
    
    sub_V0_Y1 <- data[data[[V_name]] == 0 & data[[Y_name]] == 1,]
    out$fit_G_V0_Y1_X <- glm(G_X_Y1_model, data = sub_V0_Y1, family = family)
    
  }
  
  return(out)
  
}


#' Function to fit all SuperLearner models needed by various estimators
#' 
#' @param data dataset used to fit models
#' @param G_name name of growth outcome variable, default G
#' @param V_name name of vaccine treatment variable, default V
#' @param X_name character vector containing name(s) of covariates, default X
#' @param Y_name name of infection variable, default Y
#' @param est character vector of names of estimators to use for growth effect
#' @param G_X_library specify SuperLearner libraries to be used for fitting growth on covariates. Default GLM 
#' @param Y_X_library specify SuperLearner libraries to be used for fitting infection on covariates. Default GLM 
#' @param family family for outcome model, defaults to gaussian for growth
#' @param v_folds number of cross-validation folds to use in SuperLearner, default 3
#' 
#' @export
#' 
#' @returns list with models needed for specified estimators
fit_ml_models <- function(data,
                       G_name = "G",
                       V_name = "V",
                       X_name = c("X"),
                       Y_name = "Y",
                       est = c("gcomp_pop_estimand", 
                               "gcomp", 
                               "efficient_aipw", 
                               "efficient_tmle"),
                       G_V_X_library = c("SL.glm"),
                       G_X_library = c("SL.glm"),
                       Y_X_library = c("SL.glm"),
                       family = "gaussian",
                       v_folds = 3){
  
  out <- list()
  
  # only needed for population estimator
  if("gcomp_pop_estimand" %in% est){
    out$fit_G_V_X <- SuperLearner::SuperLearner(Y = data[[G_name]],
                                                X = data[, colnames(data) %in% c(V_name, X_name), drop = FALSE],
                                                family = family,
                                                SL.library = G_V_X_library, 
                                                cvControl = list(V = v_folds))
    
  } 
  
  # needed for any weight-based estimator
  if(any(c("gcomp", "efficient_aipw", "efficient_tmle") %in% est)){
    sub_V0 <- data[data[[V_name]] == 0,]
    out$fit_Y_V0_X <- SuperLearner::SuperLearner(Y = sub_V0[[Y_name]],
                                                X = sub_V0[, X_name, drop = FALSE],
                                                family = stats::binomial(),
                                                SL.library = Y_X_library, 
                                                cvControl = list(V = v_folds))
    
    sub_V1 <- data[data[[V_name]] == 1,]
    out$fit_Y_V1_X <- SuperLearner::SuperLearner(Y = sub_V1[[Y_name]],
                                                 X = sub_V1[, X_name, drop = FALSE],
                                                 family = stats::binomial(),
                                                 SL.library = Y_X_library, 
                                                 cvControl = list(V = v_folds))
    
    sub_V1_Y1 <- data[data[[V_name]] == 1 & data[[Y_name]] == 1,]
    out$fit_G_V1_Y1_X <- SuperLearner::SuperLearner(Y = sub_V1_Y1[[G_name]],
                                                 X = sub_V1_Y1[, X_name, drop = FALSE],
                                                 family = family,
                                                 SL.library = G_X_library, 
                                                 cvControl = list(V = v_folds))
    
    sub_V1_Y0 <- data[data[[V_name]] == 1 & data[[Y_name]] == 0,]
    out$fit_G_V1_Y0_X <- SuperLearner::SuperLearner(Y = sub_V1_Y0[[G_name]],
                                                    X = sub_V1_Y0[, X_name, drop = FALSE],
                                                    family = family,
                                                    SL.library = G_X_library, 
                                                    cvControl = list(V = v_folds))
    
    sub_V0_Y1 <- data[data[[V_name]] == 0 & data[[Y_name]] == 1,]
    out$fit_G_V0_Y1_X <- SuperLearner::SuperLearner(Y = sub_V0_Y1[[G_name]],
                                                    X = sub_V0_Y1[, X_name, drop = FALSE],
                                                    family = family,
                                                    SL.library = G_X_library, 
                                                    cvControl = list(V = v_folds))
    
  }
  
  return(out)
  
}
