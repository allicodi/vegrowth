#' Function for g-computation of counterfactual post-infection outcomes in the 
#' doomed principal strata
#' 
#' @param data dataset to predict on
#' @param models list of pre-fit models needed for estimation
#' 
#' @returns g-comp estimate of growth effect in the doomed strata
do_gcomp_doomed <- function(data, models){
  
  # Psi_1 = E[P(S=1 | Z = 0, X) / P(S= 1 | Z = 0) * E[Y | Z=1, X] ]
  if(inherits(models$fit_Y_Z1_S1_X, "SuperLearner")){
    mu_11_X <- predict(models$fit_Y_Z1_S1_X, newdata = data)$pred
    mu_10_X <- predict(models$fit_Y_Z0_S1_X, newdata = data)$pred
    rho_1_X <- predict(models$fit_S_Z1_X, newdata = data)$pred
  } else{
    mu_11_X <- predict(models$fit_Y_Z1_S1_X, newdata = data, type = "response")
    mu_10_X <- predict(models$fit_Y_Z0_S1_X, newdata = data, type = "response")
    rho_1_X <- predict(models$fit_S_Z1_X, newdata = data, type = "response")
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
    S_name, Y_name, Z_name
){
  
  # Psi_1 = E[P(S=1 | Z = 0, X) / P(S= 1 | Z = 0) * E[Y | Z=1, X] ]
  if(inherits(models$fit_Y_Z1_S1_X, "SuperLearner")){
    rho_1_X <- predict(models$fit_S_Z1_X, newdata = data)$pred
    rho_0_X <- predict(models$fit_S_Z0_X, newdata = data)$pred
    pi_1_X <- predict(models$fit_Z_X, newdata = data)$pred
  } else{
    rho_1_X <- predict(models$fit_S_Z1_X, newdata = data, type = "response")
    rho_0_X <- predict(models$fit_S_Z0_X, newdata = data, type = "response")
    pi_1_X <- predict(models$fit_Z_X, newdata = data, type = "response")
  }
  pi_0_X <- 1 - pi_1_X
  rho_bar_1 <- mean(rho_1_X)
  
  S <- data[[S_name]]
  Y <- data[[Y_name]]
  Z <- data[[Z_name]]
  
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
#' @param Y_name name of growth outcome variable, default Y
#' @param Z_name name of vaccine treatment variable, default Z
#' @param S_name name of infection variable, default Y
#' @param return_se flag to return standard error, defualt FALSE
#' 
#' @returns AIPW estimate of growth effect in doomed stratum (+ standard error if return_se = TRUE)
do_aipw_doomed <- function(data, 
                           models,
                           Y_name = "Y",
                           Z_name = "Z",
                           S_name = "S",
                           return_se = FALSE){
  
  if(inherits(models$fit_S_Z0_X, "SuperLearner")){
    mu_11_X <- predict(models$fit_Y_Z1_S1_X, newdata = data)$pred
    mu_01_X <- predict(models$fit_Y_Z0_S1_X, newdata = data)$pred
    
    rho_0_X <- predict(models$fit_S_Z0_X, newdata = data)$pred
    rho_1_X <- predict(models$fit_S_Z1_X, newdata = data)$pred
    
    pi_1_X <- predict(models$fit_Z_X, newdata = data)$pred
  } else{
    mu_11_X <- predict(models$fit_Y_Z1_S1_X, newdata = data, type = "response")
    mu_01_X <- predict(models$fit_Y_Z0_S1_X, newdata = data, type = "response")
    
    rho_0_X <- predict(models$fit_S_Z0_X, newdata = data, type = "response")
    rho_1_X <- predict(models$fit_S_Z1_X, newdata = data, type = "response")
    
    pi_1_X <- models$fit_Z_X$fitted.values
  }
  
  pi_0_X <- 1 - pi_1_X
  rho_bar_1 <- mean(rho_0)
  
  eta_tilde_0 <- rho_1_X * mu_11_X / rho_bar_1
  eta_tilde_1 <- rho_1_X * mu_01_X / rho_bar_1
  
  eta_0_plugin <- mean( eta_tilde_0 )
  eta_1_plugin <- mean( eta_tilde_1 )
  
  Z <- data[[Z_name]]
  S <- data[[S_name]]
  Y <- data[[Y_name]]
  
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
  growth_effect <- eta_1_aipw - eta_0_aipw
  se <- sqrt(var(augmentation_1 - augmentation_0) / dim(data)[1])
  
  # Multiplicative effect (log scale)
  growth_effect_log_mult <- log(eta_1_aipw / eta_0_aipw)
  
  # Yet SE using IF matrix same way as TMLE
  if_matrix <- cbind(augmentation_1, augmentation_0)
  cov_matrix <- cov(if_matrix) / dim(data)[1]
  
  gradient <- matrix(c(1 / psi_1_aipw, -1 / psi_0_aipw), ncol = 1)
  
  se_log_mult_eff <- sqrt(t(gradient) %*% cov_matrix %*% gradient)
  
  if(return_se){
    out <- c(growth_effect, se, growth_effect_log_mult, se_log_mult_eff)
    names(out) <- c("additive_effect", "additive_se", "log_multiplicative_effect", "log_multiplicative_se")
    return(out)
  }else{
    out <- c(growth_effect, growth_effect_log_mult)
    names(out) <- c("additive_effect", "log_multiplicative_effect")
    return(out)
  }
}

#' Function for bounds on doomed estimate without use of cross-world assumption
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param Y_name growth outcome variable name
#' @param Z_name vaccination variable name
#' @param S_name infection variable name
#' @param family gaussian for continuous outcome, binomial for binary outcome
#' 
#' @returns list containing estimate of E[Y(1) | Y(0) = 1], bounds on E[Y(0) | Y(0) = 1], bounds on additive effect, bounds on multiplicative effect
get_bound_doomed <- function(
    data, 
    Y_name = "Y",
    Z_name = "Z",
    S_name = "S",
    family = "gaussian"
){
  
  # Step 1: rhobar_z_n
  
  # 1.1 rhobar_0_n (or mean in subset) -- size of infected in unvaccinated (naturally infected, doomed + protected)
  rhobar_0_n <- mean(data[[S_name]][data[[Z_name]] == 0])
  
  # 1.2 rhobar_1_n -- size of infected in vaccinated (doomed)
  rhobar_1_n <- mean(data[[S_name]][data[[Z_name]] == 1])
  
  if(rhobar_0_n > rhobar_1_n){
    
    # Step 2: q_n (relative size of protected in naturally infected (doomed + protected) in unvax)
    q_n = 1 - rhobar_1_n / rhobar_0_n
    
    # Step 3: q_n^th quintiles of S__Z1_S0 (aka Y__Z1_S0, need to rename everything at some point)
    Y__Z0_S1 <- data[[Y_name]][which(data[[Z_name]] == 0 & data[[S_name]] == 1)]
    q_nth_quintile <- quantile(Y__Z0_S1, probs = q_n)
    one_minus_q_nth_quintile <- quantile(Y__Z0_S1, probs = 1 - q_n)
    
    # Step 4: mubar_10_l,u_n 
    if(family == "gaussian"){
      # get people < quintile lower (protected people)
      mubar_10_u_n <- sum(data[[Y_name]] * as.numeric(data[[S_name]] == 1 & data[[Z_name]] == 0 & data[[Y_name]] > q_nth_quintile )) / 
        sum(as.numeric(data[[S_name]] == 1 & data[[Z_name]] == 0 & data[[Y_name]] > q_nth_quintile ))
      
      # get people > quintile upper (protected people)
      mubar_10_l_n <- sum(data[[Y_name]] * as.numeric(data[[S_name]] == 1 & data[[Z_name]] == 0 & data[[Y_name]] < one_minus_q_nth_quintile )) / 
        sum(as.numeric(data[[S_name]] == 1 & data[[Z_name]] == 0 & data[[Y_name]] < one_minus_q_nth_quintile ))
    } else{
      # Binary outcome
      
      # Unvaccinated, infected (naturally infected in placebo arm)
      data__Z0_S1 <- data[which(data[[Z_name]] == 0 & data[[S_name]] == 1),]
      
      # q_n = proportion of people in naturally infected, placebo arm that are protected
      # 1 - q_n = proportion of people in naturally infected, placebo arm that are doomed
      # target_num = proportion of doomed people * number of people in naturally infected = target number of doomed people
      target_num <- ceiling((1 - q_n)*nrow(data__Z0_S1)) # target number of doomed in placebo arm
      num_0s <- length(which(data__Z0_S1[[Y_name]] == 0))  # observed number of 0s in Doomed + Protected (naturally infected) box in the placebo arm
      num_1s <- length(which(data__Z0_S1[[Y_name]] == 1))  # observed number of 1s in Doomed + Protected (naturally infected) box in the placebo arm
      
      prop_0s__Z0_S1 <- 1 - mean(data__Z0_S1[[Y_name]])
      # check equivalence
      (num_0s >= target_num) == (prop_0s__Z0_S1 > rhobar_1_n / rhobar_0_n)
      ## Lower Bound:
      
      # Check if at least target_num 0s in the vax uninfected 
      if(num_0s >= target_num){
        # If so, muhat_10_l = 0
        mubar_10_l_n <- 0
        mubar_10_l_n_2 <- 0
      } else{
        # Else, determine fraction of 1s that need to be kept
        
        # (number of 0s to be size of doomed - num 0s observed) / number of 0s to be size of doomed = proportion of 1s to be added
        mubar_10_l_n <- ( (target_num - num_0s) / target_num )
      }
      
      ## Upper Bound:
      
      # Check if at least q_n * 100 % 1s in the vax uninfected 
      if(num_1s >= target_num){
        # If so, mubar_10_u = 1
        mubar_10_u_n <- 1
        mubar_10_u_n_2 <- 1
      } else{
        # Else, proportion of 1s in doomed box
        
        mubar_10_u_n <- num_1s / target_num
        
      }
      
    }
    
    l_n <- mubar_10_l_n
    u_n <- mubar_10_u_n
    
    #mean in vaccinated infecteds for comparison
    E_Y1__S0_1 <- mean(data[[Y_name]][data[[S_name]] == 1 & data[[Z_name]] == 1])
    
    out <- list(E_Y1__S0_1 = E_Y1__S0_1,
                E_Y0__S0_1_lower = l_n,
                E_Y0__S0_1_upper = u_n,
                additive_effect_lower = E_Y1__S0_1 - l_n,
                additive_effect_upper = E_Y1__S0_1 - u_n,
                mult_effect_lower = E_Y1__S0_1 / l_n,
                mult_effect_upper = E_Y1__S0_1 / u_n,
                success = 1)
    
  } else{
    print("Method not applicable unless evidence of vaccine protection.")
    
    #mean in vaccinated infecteds for comparison
    E_Y1__S0_1 <- mean(data[[Y_name]][data[[S_name]] == 1 & data[[Z_name]] == 1])
    
    out <- list(E_Y1__S0_1 = NA,
                E_Y0__S0_1_lower = NA,
                E_Y0__S0_1_upper = NA,
                additive_effect_lower = NA,
                additive_effect_upper = NA,
                mult_effect_lower = NA,
                mult_effect_upper = NA,
                success = 0)
  }
  
  class(out) <- "bound_doomed"
  
  return(out)
  
}