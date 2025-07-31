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

#' Function for bounds on doomed estimate without use of cross-world assumption
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param G_name growth outcome variable name
#' @param V_name vaccination variable name
#' @param Y_name infection variable name
#' @param family gaussian for continuous outcome, binomial for binary outcome
#' 
#' @returns list containing estimate of E[G(1) | Y(0) = 1], bounds on E[G(0) | Y(0) = 1], bounds on additive effect, bounds on multiplicative effect
get_bound_doomed <- function(
    data, 
    G_name = "G",
    V_name = "V",
    Y_name = "Y",
    family = "gaussian"
){
  
  # Step 1: rhobar_z_n
  
  # 1.1 rhobar_0_n (or mean in subset) -- size of infected in unvaccinated (naturally infected, doomed + protected)
  rhobar_0_n <- mean(data[[Y_name]][data[[V_name]] == 0])
  
  # 1.2 rhobar_1_n -- size of infected in vaccinated (doomed)
  rhobar_1_n <- mean(data[[Y_name]][data[[V_name]] == 1])
  
  if(rhobar_0_n > rhobar_1_n){
    
    # Step 2: q_n (relative size of protected in naturally infected (doomed + protected) in unvax)
    q_n = 1 - rhobar_1_n / rhobar_0_n
    
    # Step 3: q_n^th quintiles of Y__Z1_S0 (aka G__V1_Y0, need to rename everything at some point)
    G__V0_Y1 <- data[[G_name]][which(data[[V_name]] == 0 & data[[Y_name]] == 1)]
    q_nth_quintile <- quantile(G__V0_Y1, probs = q_n)
    one_minus_q_nth_quintile <- quantile(G__V0_Y1, probs = 1 - q_n)
    
    # Step 4: mubar_10_l,u_n 
    if(family == "gaussian"){
      # get people < quintile lower (protected people)
      mubar_10_u_n <- sum(data[[G_name]] * as.numeric(data[[Y_name]] == 1 & data[[V_name]] == 0 & data[[G_name]] > q_nth_quintile )) / 
        sum(as.numeric(data[[Y_name]] == 1 & data[[V_name]] == 0 & data[[G_name]] > q_nth_quintile ))
      
      # get people > quintile upper (protected people)
      mubar_10_l_n <- sum(data[[G_name]] * as.numeric(data[[Y_name]] == 1 & data[[V_name]] == 0 & data[[G_name]] < one_minus_q_nth_quintile )) / 
        sum(as.numeric(data[[Y_name]] == 1 & data[[V_name]] == 0 & data[[G_name]] < one_minus_q_nth_quintile ))
    } else{
      # Binary outcome
      
      # Unvaccinated, infected (naturally infected in placebo arm)
      data__V0_Y1 <- data[which(data[[V_name]] == 0 & data[[Y_name]] == 1),]
      
      # q_n = proportion of people in naturally infected, placebo arm that are protected
      # 1 - q_n = proportion of people in naturally infected, placebo arm that are doomed
      # target_num = proportion of doomed people * number of people in naturally infected = target number of doomed people
      target_num <- ceiling((1 - q_n)*nrow(data__V0_Y1)) # target number of doomed in placebo arm
      num_0s <- length(which(data__V0_Y1[[G_name]] == 0))  # observed number of 0s in Doomed + Protected (naturally infected) box in the placebo arm
      num_1s <- length(which(data__V0_Y1[[G_name]] == 1))  # observed number of 1s in Doomed + Protected (naturally infected) box in the placebo arm
      
      prop_0s__V0_Y1 <- 1 - mean(data__V0_Y1[[G_name]])
      # check equivalence
      (num_0s >= target_num) == (prop_0s__V0_Y1 > rhobar_1_n / rhobar_0_n)
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
    
  } else{
    stop("Method not applicable unless evidence of vaccine protection.")
  }
  
  #mean in vaccinated infecteds for comparison
  E_G1__Y0_1 <- mean(data[[G_name]][data[[Y_name]] == 1 & data[[V_name]] == 1])
  
  out <- list(E_G1__Y0_1 = E_G1__Y0_1,
              E_G0__Y0_1_lower = l_n,
              E_G0__Y0_1_upper = u_n,
              additive_effect_lower = E_G1__Y0_1 - l_n,
              additive_effect_upper = E_G1__Y0_1 - u_n,
              mult_effect_lower = E_G1__Y0_1 / l_n,
              mult_effect_upper = E_G1__Y0_1 / u_n)
  
  return(out)
  
}