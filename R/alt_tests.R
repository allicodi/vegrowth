
#' Function to get chop-lump style test-statistic
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param G_name growth outcome variable name
#' @param V_name vaccination variable name
#' @param Y_name infection variable name
#' 
#' @returns dataframe with chop-lump style test statistics for mean in vax, mean in placebo
#' 
get_chop_lump_statistic <- function(data,
                                    G_name = "G",
                                    V_name = "V",
                                    Y_name = "Y"){
  
  # 1. Get number of people in relevant groups
  n_no_inf_plc <- sum(data[[Y_name]] == 0 & data[[V_name]] == 0)
  n_no_inf_vax <- sum(data[[Y_name]] == 0 & data[[V_name]] == 1)
  n_plc <- sum(data[[V_name]] == 0)
  n_vax <- sum(data[[V_name]] == 1)
  n_inf_plc <- n_plc - n_no_inf_plc
  n_inf_vax <- n_vax - n_no_inf_vax
  
  # 2. Chop based on group with more infections
  if(n_inf_plc > n_inf_vax){
    
    # placebo - everyone is infected, simple mean
    mean_G_plc <- mean(data[[G_name]][data[[V_name]] == 0 & data[[Y_name]] == 1])
    
    # vaccinated - weighted mean
    mean_G_noinf_vax <- mean(data[[G_name]][data[[V_name]] == 1 & data[[Y_name]] == 0])
    mean_G_inf_vax <- mean(data[[G_name]][data[[V_name]] == 1 & data[[Y_name]] == 1])
    
    mean_G_vax  <- mean_G_noinf_vax * ((n_inf_plc - n_inf_vax) / n_inf_plc) +
      mean_G_inf_vax * (n_inf_vax / n_inf_plc)
    
  } else if (n_inf_plc < n_inf_vax){
    
    # vaccinated - everyone is infected, simple mean
    mean_G_vax <- mean(data[[G_name]][data[[V_name]] == 1 & data[[Y_name]] == 1])
    
    # placebo - weighted mean
    mean_G_noinf_plc <- mean(data[[G_name]][data[[V_name]] == 0 & data[[Y_name]] == 0])
    mean_G_inf_plc <- mean(data[[G_name]][data[[V_name]] == 0 & data[[Y_name]] == 1])
    
    mean_G_plc  <- mean_G_noinf_plc * ((n_inf_vax - n_inf_plc) / n_inf_vax) +
      mean_G_inf_plc * (n_inf_plc / n_inf_vax)
    
  } else{
    
    # vaccinated - everyone is infected, simple mean
    mean_G_vax <- mean(data[[G_name]][data[[V_name]] == 1 & data[[Y_name]] == 1])
    
    # placebo - everyone is infected, simple mean
    mean_G_plc <- mean(data[[G_name]][data[[V_name]] == 0 & data[[Y_name]] == 1])
    
  }
  
  return(data.frame(mean_G_plc = mean_G_plc,
                    mean_G_vax = mean_G_vax))
}

#' Function to do permutation test for chop-lump style test-statistics
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param G_name growth outcome variable name
#' @param V_name vaccination variable name
#' @param Y_name infection variable name
#' @param n_permutations number of permutations to complete
#' 
#' @returns chop lump test statistic
do_chop_lump_test <- function(data, 
                              G_name = "G",
                              V_name = "V",
                              Y_name = "Y",
                              n_permutations = 1e4){
  
  original_means <- get_chop_lump_statistic(data, 
                                            G_name = G_name,
                                            V_name = V_name,
                                            Y_name = Y_name)
  ## Permutation approach
  null_means <- vector("list", length = n_permutations)
  for(i in 1:n_permutations){
    data_shuffle <- data
    data_shuffle[[V_name]] <- sample(data_shuffle[[V_name]])
    
    null_means[[i]] <- get_chop_lump_statistic(data_shuffle,
                                               G_name = G_name,
                                               V_name = V_name,
                                               Y_name = Y_name)
  }
  
  null_df <- do.call(rbind, null_means)
  
  ## Hypothesis test
  observed_diff <- original_means$mean_G_vax - original_means$mean_G_plc
  null_df$mean_diff <- null_df$mean_G_vax - null_df$mean_G_plc
  
  out <- list(
    obs_diff = observed_diff,
    null_diffs = null_df$mean_diff,
    pval = mean(null_df$mean_diff > observed_diff)
  )
  
  class(out) <- "chop_lump_res"
  return(out)
}


#' Helper function to plot chop-lump test results
#' 
#' @param out chop_lump_res object to plot
#' 
#' @returns A histogram plot with the distribution of test statistics under the null hypothesis,
#' and a vertical line showing the observed test statistic.
plot.chop_lump_res <- function(out){
  hist(
    out$null_diffs,
    xlab = "Test statistic",
    xlim = range(c(out$null_diffs, out$observed_diff)),
    main = "Distribution of test statistic under null"
  )
  abline(v = out$observed_diff, col = 2, lwd = 2)
}

# -----------------------------------------------------------------------------

# Hudgens Method --------------------------------------------------------------

#' Function for hudgens-style test statistic
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param G_name growth outcome variable name
#' @param V_name vaccination variable name
#' @param Y_name infection variable name
#' @param lower_bound A boolean. If TRUE, then adds the smallest growth measures 
#'    to the infected vaccinees thereby yielding a lower
#'    bound on the effect of interest. If FALSE, then adds the largest
#'    growth measures to the infected vaccinees thereby yielding an upper
#'    bound on the effect of interest. 
#' 
#' @returns Hudgens-style test statistic
get_hudgens_stat <- function(
    data, 
    G_name = "G",
    V_name = "V",
    Y_name = "Y",
    lower_bound = TRUE
){
  n_no_inf_plc <- sum(data[[Y_name]] == 0 & data[[V_name]] == 0)
  n_no_inf_vax <- sum(data[[Y_name]] == 0 & data[[V_name]] == 1)
  n_plc <- sum(data[[V_name]] == 0)
  n_vax <- sum(data[[V_name]] == 1)
  n_inf_plc <- n_plc - n_no_inf_plc
  n_inf_vax <- n_vax - n_no_inf_vax
  
  n_diff <- n_inf_plc - n_inf_vax
  # rho_0n > rho_1n
  if(n_inf_plc / n_plc > n_inf_vax / n_vax){
    
    # E[G | Y_inf = 1, V = 0]
    mean_G_plc <- mean(data[[G_name]][data[[V_name]] == 0 & data[[Y_name]] == 1])
    
    G_inf_vax <- data[[G_name]][data[[V_name]] == 1 & data[[Y_name]] == 1]
    G_noinf_vax <- data[[G_name]][data[[V_name]] == 1 & data[[Y_name]] == 0]
    
    if(lower_bound){
      # find the n_diff smallest growths in non-infected vaccinees 
      G_aug_vax <- sort(G_noinf_vax, decreasing = FALSE)[1:n_diff]
    }else{
      # find the n_diff largest growths in non-infected vaccinees 
      G_aug_vax <- sort(G_noinf_vax, decreasing = TRUE)[1:n_diff]
    }
    
    mean_G_vax <- mean(c(G_inf_vax, G_aug_vax))
    
  }else{
    
    stop("Method not applicable unless evidence of vaccine protection.")
    
  }
  
  additive <- mean_G_vax - mean_G_plc
  #log_mult <- log(mean_G_vax / mean_G_plc)
  
  #out <- c(additive, multiplicative)
  #names(out) <- c("additive_effect", "log_multiplicative_effect")
  #return(out)
  
  return(additive)
  
}

#' Function for bootstrap test of Hudgens-style test statistic
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param G_name growth outcome variable name
#' @param V_name vaccination variable name
#' @param Y_name infection variable name
#' @param lower_bound A boolean. If TRUE, then adds the smallest growth measures 
#'    to the infected vaccinees thereby yielding a lower
#'    bound on the effect of interest. If FALSE, then adds the largest
#'    growth measures to the infected vacccinees thereby yielding an upper
#'    bound on the effect of interest.
#'  @param n_boot number of bootstrap replicates
#'  @param n_boot_try max number of attempts for bootstrap resampling
#'  @param effect_dir direction of beneficial effect, defaults to "positive" for beneficial outcome. Used for one-side tests of bounds.  
#'  
#'  @returns list with observed difference between plc and vax, Hudgens-style test statistic, and p-value
hudgens_test <- function(
    data, 
    G_name = "G",
    V_name = "V",
    Y_name = "Y",
    lower_bound = TRUE,
    n_boot = 1e3, 
    n_boot_try = n_boot*10,
    effect_dir = "positive"
){
  n_no_inf_plc <- sum(data[[Y_name]] == 0 & data[[V_name]] == 0)
  n_no_inf_vax <- sum(data[[Y_name]] == 0 & data[[V_name]] == 1)
  n_plc <- sum(data[[V_name]] == 0)
  n_vax <- sum(data[[V_name]] == 1)
  n <- n_plc + n_vax
  n_inf_plc <- n_plc - n_no_inf_plc
  n_inf_vax <- n_vax - n_no_inf_vax
  n_diff <- n_inf_plc - n_inf_vax
  
  if(n_diff > 0){
    obs_diff <- get_hudgens_stat(
      data = data, lower_bound = lower_bound, Y_name = Y_name, G_name = G_name, V_name = V_name
    )
    boot_diff <- rep(NA, n_boot)
    success_ct <- 0
    attempt_ct <- 0
    while(success_ct < n_boot & attempt_ct <= n_boot_try){
      attempt_ct <- attempt_ct + 1
      boot_idx <- sample(seq_len(n), replace = TRUE)
      boot_data <- data[boot_idx, , drop = FALSE]
      
      n_no_inf_plc_boot <- sum(boot_data[[Y_name]] == 0 & boot_data[[V_name]] == 0)
      n_no_inf_vax_boot <- sum(boot_data[[Y_name]] == 0 & boot_data[[V_name]] == 1)
      n_plc_boot <- sum(boot_data[[V_name]] == 0)
      n_vax_boot <- sum(boot_data[[V_name]] == 1)
      n_inf_plc_boot <- n_plc_boot - n_no_inf_plc_boot
      n_inf_vax_boot <- n_vax_boot - n_no_inf_vax_boot
      
      if(n_inf_plc_boot > n_inf_vax_boot){
        success_ct <- success_ct + 1
        boot_diff[success_ct] <- get_hudgens_stat(boot_data, lower_bound = lower_bound, Y_name = Y_name, G_name = G_name, V_name = V_name)
      }
    }
    if(sum(!is.na(boot_diff)) != n_boot){
      warning(paste0(
        "Only achieved ", sum(!is.na(boot_diff)), " successful resamples. Try again with larger n_boot_try? \n sample size: ", nrow(data)
      ))
    }
    
    if(sum(!is.na(boot_diff)) > 2){
      sd_boot_diff <- sd(boot_diff, na.rm = TRUE)
      test_stat <- obs_diff / sd_boot_diff
      
      if(effect_dir == "positive"){
        pval <- pnorm(test_stat, lower.tail = FALSE)
      } else{
        pval <- pnorm(test_stat, lower.tail = TRUE)
      }
      
      # one-sided p-value
      # pval <- pnorm(test_stat, lower.tail = FALSE)
      
      out <- list(
        obs_diff = obs_diff,
        test_stat = test_stat,
        pval = pval
      )
      
      class(out) <- "hudgens_res"
      return(out)
    }else{
      stop("Did not achieve 2 successful bootstrap resamples. Try again with larger n_boot_try? Or give up?")
    }
    
  }else{
    stop("Method not applicable unless evidence of vaccine protection.")
  }
}

#' Function for hudgens-style test statistic - match description in the paper
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param G_name growth outcome variable name
#' @param V_name vaccination variable name
#' @param Y_name infection variable name
#' @param family gaussian for continuous outcome, binomial for binary outcome
#' 
#' @returns Hudgens-style test statistic
get_hudgens_stat_new <- function(
    data, 
    G_name = "G",
    V_name = "V",
    Y_name = "Y",
    family = "gaussian"
){
  
  # Step 1: rhobar_z_n
  
  # 1.1 rhobar_0_n (or mean in subset)
  rhobar_0_n <- mean(data[[Y_name]][data[[V_name]] == 0])
  
  # rhobar_0_n <- sum(data[[Y_name]]*as.numeric(data[[V_name]] == 0)) / 
  #   sum(as.numeric(data[[V_name]] == 0))
  
  # 1.2 rhobar_1_n
  rhobar_1_n <- mean(data[[Y_name]][data[[V_name]] == 1])
  
  # rhobar_1_n <- sum(data[[Y_name]]*as.numeric(data[[V_name]] == 1)) / 
  #   sum(as.numeric(data[[V_name]] == 1))
  
  if(rhobar_0_n > rhobar_1_n){
    # Step 2: mubar_11_n 
    mubar_11_n <- sum(data[[G_name]]*data[[Y_name]]*data[[V_name]]) / sum(data[[Y_name]]*data[[V_name]])
    
    # Step 3: q_n (relative size of protected? in (immune + protected) in vax)
    q_n = 1 - (1 - rhobar_0_n) / (1 - rhobar_1_n)
    #one_minus_q_n <- 1 - q_n # relative size of the immune? unclear why added?? 
    
    # Step 4: q_n^th quintiles of Y__Z1_S0 (aka G__V1_Y0, need to rename everything at some point)
    G__V1_Y0 <- data[[G_name]][which(data[[V_name]] == 1 & data[[Y_name]] == 0)]
    q_nth_quintile <- quantile(G__V1_Y0, probs = q_n)
    one_minus_q_nth_quintile <- quantile(G__V1_Y0, probs = 1 - q_n)
    
    # Step 5: mubar_10_l,u_n 
    if(family == "gaussian"){
      mubar_10_l_n <- sum(data[[G_name]] * as.numeric(data[[Y_name]] == 0 & data[[V_name]] == 1 & data[[G_name]] < q_nth_quintile )) / 
        sum(as.numeric(data[[Y_name]] == 0 & data[[V_name]] == 1 & data[[G_name]] < q_nth_quintile ))
      
      mubar_10_u_n <- sum(data[[G_name]] * as.numeric(data[[Y_name]] == 0 & data[[V_name]] == 1 & data[[G_name]] > one_minus_q_nth_quintile )) / 
        sum(as.numeric(data[[Y_name]] == 0 & data[[V_name]] == 1 & data[[G_name]] > one_minus_q_nth_quintile ))
    } else{
      # Binary outcome
      
      # Vaccinated, Uninfected
      data__V1_Y0 <- data[which(data[[V_name]] == 1 & data[[Y_name]] == 0),]
      
      target_num <- ceiling(q_n * nrow(data__V1_Y0))
      num_0s <- length(which(data__V1_Y0[[G_name]] == 0))
      num_1s <- length(which(data__V1_Y0[[G_name]] == 1))
      
      ## Lower Bound:
      
      # Check if at least q_n * 100 % 0s in the vax uninfected 
      if(num_0s >= target_num){
        # If so, muhat_10_l = 0
        mubar_10_l_n <- 0
      } else{
        # Else, mubar_10_l = (q_n - prop zeros in vax uninf) / q_n
        mubar_10_l_n <- (q_n - (num_0s / nrow(data__V1_Y0)) ) /
          q_n
      }
      
      ## Upper Bound:
      
      # Check if at least q_n * 100 % 1s in the vax uninfected 
      if(num_1s >= target_num){
        # If so, mubar_10_u = 1
        mubar_10_u_n <- 1
      } else{
        # Else, mubar_10_u = (q_n - prop ones in vax uninf) / q_n
         
        mubar_10_u_n <- 1 - ((q_n - (num_1s / nrow(data__V1_Y0)) ) /
          q_n)
      }

      # ------ ------ ------ ------ ------ ------ ------ ------ ------ ------
      # ^^^^ NOT WORKING FOR BINARY OUTCOME, Q_NTH_QUINTILE == 0 OR 1 SO NOTHING < OR >
      # randomly sample q%??
      
      # # Get subset uninfected
      # data__Y0 <- data[which(data[[Y_name]] == 0),]
      # data__Y0$idx <- 1:nrow(data__Y0)
      # 
      # # which rows match q_n^th and 1 - q_n^th quintiles (likely 0 and 1)
      # rows_q_nth_quintile <- which(data__Y0[[G_name]] == q_nth_quintile)
      # rows_one_minus_q_nth_quintile <- which(data__Y0[[G_name]] == one_minus_q_nth_quintile)
      # 
      # # number of rows to sample = q_n% of uninfected participants (round up to whole number)
      # num_q_nth <- ceiling(q_n * nrow(data__Y0))
      # 
      # # randomly remove q_n% of the people matching each quintile? not going to be the same estimate every time
      # rm_lower <- sample(rows_q_nth_quintile, size = num_q_nth, replace = FALSE)
      # rm_upper <- sample(rows_one_minus_q_nth_quintile, size = num_q_nth, replace = FALSE)
      # 
      # # add removal flag to data 
      # data__Y0$rm_lower <- ifelse(data__Y0$idx %in% rm_lower, 1, 0)
      # data__Y0$rm_upper <- ifelse(data__Y0$idx %in% rm_upper, 1, 0)
      # 
      # mubar_10_l_n <- sum(data__Y0[[G_name]] * data__Y0$rm_lower) / 
      #   sum(data__Y0$rm_lower)
      # 
      # mubar_10_u_n <- sum(data__Y0[[G_name]] * data__Y0$rm_upper) / 
      #   sum(data__Y0$rm_upper)
      
      # ^^ results from this are equivalent to what i changed it to after discussion??
      
    }
    
    # Step 6: final estimates of the bounds (just doing both each time for now??)
    
    l_n <- mubar_11_n * (rhobar_1_n / rhobar_0_n) + mubar_10_l_n * (1 - (rhobar_1_n / rhobar_0_n))
    u_n <- mubar_11_n * (rhobar_1_n / rhobar_0_n) + mubar_10_u_n * (1 - (rhobar_1_n / rhobar_0_n))
    
  } else{
    stop("Method not applicable unless evidence of vaccine protection.")
  }
  
  #mean in unvaccinated infecteds for comparison
  E_G0__Y0_1 <- mean(data[[G_name]][data[[Y_name]] == 1 & data[[V_name]] == 0])

  out <- list(E_G0__Y0_1 = E_G0__Y0_1,
              E_G1__Y0_1_lower = l_n,
              E_G1__Y0_1_upper = u_n,
              additive_effect_lower = l_n - E_G0__Y0_1,
              additive_effect_upper = u_n - E_G0__Y0_1,
              mult_effect_lower = l_n / E_G0__Y0_1,
              mult_effect_upper = u_n / E_G0__Y0_1)
  
  return(out)

}

#' Function for bootstrap test of Hudgens-style test statistic
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param G_name growth outcome variable name
#' @param V_name vaccination variable name
#' @param Y_name infection variable name
#' @param family gaussian for continuous outcome, binomial for binary
#' @param n_boot number of bootstrap replicates
#' @param n_boot_try max number of attempts for bootstrap resampling
#' @param effect_dir direction of beneficial effect, defaults to "positive" for beneficial outcome. Used for one-side tests of bounds.  
#'  
#' @returns list with observed difference between plc and vax, Hudgens-style test statistic, and p-value
hudgens_test_new <- function(
    data, 
    G_name = "G",
    V_name = "V",
    Y_name = "Y",
    family = "gaussian",
    n_boot = 1e3, 
    n_boot_try = n_boot*10,
    effect_dir = "positive"
){
  
  n <- nrow(data)
  
  # Step 1: rhobar_z_n
  rhobar_0_n <- mean(data[[Y_name]][data[[V_name]] == 0])
  rhobar_1_n <- mean(data[[Y_name]][data[[V_name]] == 1])

  if(rhobar_0_n > rhobar_1_n){
    
    out <- get_hudgens_stat_new(
      data = data, Y_name = Y_name, G_name = G_name, V_name = V_name, family = family
    )
    
    boot_diff <- vector(mode = "list", length = n_boot)
    success_ct <- 0
    attempt_ct <- 0
    while(success_ct < n_boot & attempt_ct <= n_boot_try){
      attempt_ct <- attempt_ct + 1
      boot_idx <- sample(seq_len(n), replace = TRUE)
      boot_data <- data[boot_idx, , drop = FALSE]
      
      rhobar_0_n_boot <- mean(boot_data[[Y_name]][boot_data[[V_name]] == 0])
      rhobar_1_n_boot <- mean(boot_data[[Y_name]][boot_data[[V_name]] == 1])
      
      if(rhobar_0_n_boot > rhobar_1_n_boot){
        success_ct <- success_ct + 1
        boot_diff[[success_ct]] <- get_hudgens_stat_new(
          data = boot_data, Y_name = Y_name, G_name = G_name, V_name = V_name, family = family
        )
      }
    }
    if(sum(!is.na(boot_diff)) != n_boot){
      warning(paste0(
        "Only achieved ", sum(!is.na(boot_diff)), " successful resamples. Try again with larger n_boot_try? \n sample size: ", nrow(data)
      ))
    }
    
    if(sum(!is.na(boot_diff)) > 2){
      
      boot_diff_res <- data.frame(do.call(rbind, boot_diff))
      
      sd_boot_diff_additive_lower <- sd(unlist(boot_diff_res$additive_effect_lower), na.rm = TRUE)
      sd_boot_diff_additive_upper <- sd(unlist(boot_diff_res$additive_effect_upper), na.rm = TRUE)
      sd_boot_diff_mult_lower <- sd(unlist(boot_diff_res$mult_effect_lower), na.rm = TRUE)
      sd_boot_diff_mult_upper <- sd(unlist(boot_diff_res$mult_effect_upper), na.rm = TRUE)

      out <- c(
        out,
        list(
          sd_additive_lower = sd_boot_diff_additive_lower,
          sd_additive_upper = sd_boot_diff_additive_upper,
          sd_mult_lower = sd_boot_diff_mult_lower,
          sd_mult_upper = sd_boot_diff_mult_upper
        )
      )
      
      
      class(out) <- "hudgens_res_new"
      return(out)
    }else{
      stop("Did not achieve 2 successful bootstrap resamples. Try again with larger n_boot_try? Or give up?")
    }
    
  }else{
    stop("Method not applicable unless evidence of vaccine protection.")
  }
}

#' Function for hudgens-style test statistic - match description in the paper
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param G_name growth outcome variable name
#' @param V_name vaccination variable name
#' @param Y_name infection variable name
#' @param family gaussian for continuous outcome, binomial for binary
#' 
#' @returns Hudgens-style test statistic
get_hudgens_stat_doomed_new <- function(
    data, 
    G_name = "G",
    V_name = "V",
    Y_name = "Y",
    family = "gaussian"
){
  
  # Step 1: rhobar_z_n
  
  # 1.1 rhobar_0_n (or mean in subset) -- size of infected in unvaccinated (naturally infected, doomed + protected)
  rhobar_0_n <- mean(data[[Y_name]][data[[V_name]] == 0])
  
  # rhobar_0_n <- sum(data[[Y_name]]*as.numeric(data[[V_name]] == 0)) / 
  #   sum(as.numeric(data[[V_name]] == 0))
  
  # 1.2 rhobar_1_n -- size of infected in vaccinated (doomed)
  rhobar_1_n <- mean(data[[Y_name]][data[[V_name]] == 1])
  
  # rhobar_1_n <- sum(data[[Y_name]]*as.numeric(data[[V_name]] == 1)) / 
  #   sum(as.numeric(data[[V_name]] == 1))
  
  if(rhobar_0_n > rhobar_1_n){
    # Step 2: mubar_11_n - doomed we observe
    # mubar_11_n <- sum(data[[G_name]]*data[[Y_name]]*data[[V_name]]) / sum(data[[Y_name]]*data[[V_name]])
    
    # Step 3: q_n (relative size of protected in naturally infected (doomed + protected) in unvax)
    q_n = 1 - rhobar_1_n / rhobar_0_n
    
    # Step 4: q_n^th quintiles of Y__Z1_S0 (aka G__V1_Y0, need to rename everything at some point)
    G__V0_Y1 <- data[[G_name]][which(data[[V_name]] == 0 & data[[Y_name]] == 1)]
    q_nth_quintile <- quantile(G__V0_Y1, probs = q_n)
    one_minus_q_nth_quintile <- quantile(G__V0_Y1, probs = 1 - q_n)
    
    # Step 5: mubar_10_l,u_n 
    if(family == "gaussian"){
      # get people < quintile lower (protected people)
      mubar_10_u_n <- sum(data[[G_name]] * as.numeric(data[[Y_name]] == 1 & data[[V_name]] == 0 & data[[G_name]] < q_nth_quintile )) / 
        sum(as.numeric(data[[Y_name]] == 1 & data[[V_name]] == 0 & data[[G_name]] < q_nth_quintile ))
      
      # get people > quintile upper (protected people)
      mubar_10_l_n <- sum(data[[G_name]] * as.numeric(data[[Y_name]] == 1 & data[[V_name]] == 0 & data[[G_name]] > one_minus_q_nth_quintile )) / 
        sum(as.numeric(data[[Y_name]] == 1 & data[[V_name]] == 0 & data[[G_name]] > one_minus_q_nth_quintile ))
    } else{
      # Binary outcome
      
      # Unvaccinated, infected
      data__V0_Y1 <- data[which(data[[V_name]] == 0 & data[[Y_name]] == 1),]
      
      # Number of protected
      target_num <- ceiling(q_n * nrow(data__V0_Y1))
      num_0s <- length(which(data__V0_Y1[[G_name]] == 0))
      num_1s <- length(which(data__V0_Y1[[G_name]] == 1))
      
      ## Lower Bound:
      
      # Check if at least q_n * 100 % 0s in the vax uninfected 
      if(num_0s >= target_num){
        # If so, muhat_10_l = 0
        mubar_10_l_n <- 0
      } else{
        # Else, mubar_10_l = (q_n - prop zeros in vax uninf) / q_n
        mubar_10_l_n <- (q_n - (num_0s / nrow(data__V0_Y1)) ) /
          q_n
      }
      
      ## Lower Bound:
      
      # Check if at least q_n * 100 % 1s in the vax uninfected 
      if(num_1s >= target_num){
        # If so, mubar_10_u = 1
        mubar_10_u_n <- 1
      } else{
        # Else, mubar_10_u = (q_n - prop ones in vax uninf) / q_n
        
        mubar_10_u_n <- 1 - ((q_n - (num_1s / nrow(data__V0_Y1)) ) /
                               q_n)
      }

    }
    
    # Step 6: final estimates of the bounds: effect in naturally infected in placebo arm - effect in protected*proportion protected???
    l_n <- mean(data[[G_name]][data[[Y_name]] == 1 & data[[V_name]] == 0]) - mubar_10_l_n * rhobar_1_n / rhobar_0_n
    u_n <- mean(data[[G_name]][data[[Y_name]] == 1 & data[[V_name]] == 0]) - mubar_10_u_n * rhobar_1_n / rhobar_0_n
    
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

#' Function for hudgens-style test statistic - doomed strata
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param G_name growth outcome variable name
#' @param V_name vaccination variable name
#' @param Y_name infection variable name
#' @param lower_bound A boolean. If TRUE, then SUBTRACTS the smallest growth measures 
#'    to the infected vaccinees thereby yielding a lower
#'    bound on the effect of interest. If FALSE, then SUBTRACTS the largest
#'    growth measures to the infected vacccinees thereby yielding an upper
#'    bound on the effect of interest.
#' 
#' @returns Hudgens-style test statistic
get_hudgens_stat_doomed <- function(
    data, 
    G_name = "G",
    V_name = "V",
    Y_name = "Y",
    lower_bound = TRUE
){
  n_no_inf_plc <- sum(data[[Y_name]] == 0 & data[[V_name]] == 0)
  n_no_inf_vax <- sum(data[[Y_name]] == 0 & data[[V_name]] == 1)
  n_plc <- sum(data[[V_name]] == 0)
  n_vax <- sum(data[[V_name]] == 1)
  n_inf_plc <- n_plc - n_no_inf_plc
  n_inf_vax <- n_vax - n_no_inf_vax
  
  n_diff <- n_inf_plc - n_inf_vax
  
  if(n_inf_plc > n_inf_vax){
    
    # E[G | Y_inf = 1, V = 0] (doomed in the vaccinated are observed)
    mean_G_vax <- mean(data[[G_name]][data[[V_name]] == 1 & data[[Y_name]] == 1])
    
    # E[G | Y_inf = 1, V = 0]
    #mean_G_plc <- mean(data[[G_name]][data[[V_name]] == 0 & data[[Y_name]] == 1])
    
    G_inf_plc <- data[[G_name]][data[[V_name]] == 0 & data[[Y_name]] == 1]
    G_noinf_plc <- data[[G_name]][data[[V_name]] == 0 & data[[Y_name]] == 0]
    
    # G_inf_vax <- data[[G_name]][data[[V_name]] == 1 & data[[Y_name]] == 1]
    # G_noinf_vax <- data[[G_name]][data[[V_name]] == 1 & data[[Y_name]] == 0]
    
    if(lower_bound){
      # find the n_diff smallest growths in non-infected vaccinees 
      # G_aug_vax <- sort(G_noinf_vax, decreasing = FALSE)[1:n_diff]
      
      # NEW
      # find the n_diff largest growths in the infected placebos (not necessary, just take rest of vector)
      # G_aug_plc <- sort(G_noinf_plc, decreasing = FALSE)[1:n_diff]
      
      # Get n_diff + 1 to end
      # decreasing = FALSE means keeping bigger people --> lower bound
      G_plc_doomed <- sort(G_noinf_plc, decreasing = FALSE)[(n_diff+1):length(G_noinf_plc)]
      
    }else{
      # find the n_diff largest growths in non-infected vaccinees 
      # G_aug_vax <- sort(G_noinf_vax, decreasing = TRUE)[1:n_diff]
      
      # NEW
      # find the n_diff smallest growths in the infected placebos (not necessary, just take rest of vector)
      #G_aug_plc <- sort(G_noinf_plc, decreasing = TRUE)[1:n_diff]
      
      # Get n_diff + 1 to end
      # decreasing = TRUE means keeping smaller people --> subtracting smaller so upper bound
      G_plc_doomed <- sort(G_noinf_plc, decreasing = TRUE)[(n_diff+1):length(G_noinf_plc)]
    }
    
    mean_G_plc <- mean(G_plc_doomed)
    #mean_G_vax <- mean(c(G_inf_vax, G_aug_vax))
    
  }else{
    
    stop("Method not applicable unless evidence of vaccine protection.")
    
  }
  
  return(mean_G_vax - mean_G_plc)
  
}

#' Function for bootstrap test of Hudgens-style test statistic
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param G_name growth outcome variable name
#' @param V_name vaccination variable name
#' @param Y_name infection variable name
#' @param family gaussian for continuous outcome, binomial for binary
#' @param n_boot number of bootstrap replicates
#' @param n_boot_try max number of attempts for bootstrap resampling
#' @param effect_dir direction of beneficial effect, defaults to "positive" for beneficial outcome. Used for one-side tests of bounds.  
#'  
#' @returns list with observed difference between plc and vax, Hudgens-style test statistic, and p-value
hudgens_test_doomed_new <- function(
    data, 
    G_name = "G",
    V_name = "V",
    Y_name = "Y",
    family = "gaussian",
    n_boot = 1e3, 
    n_boot_try = n_boot*10,
    effect_dir = "positive"
){
  
  n <- nrow(data)
  
  # Step 1: rhobar_z_n
  rhobar_0_n <- mean(data[[Y_name]][data[[V_name]] == 0])
  rhobar_1_n <- mean(data[[Y_name]][data[[V_name]] == 1])
  
  if(rhobar_0_n > rhobar_1_n){
    
    out <- get_hudgens_stat_doomed_new(
      data = data, Y_name = Y_name, G_name = G_name, V_name = V_name, family = family
    )
    
    boot_diff <- vector(mode = "list", length = n_boot)
    success_ct <- 0
    attempt_ct <- 0
    while(success_ct < n_boot & attempt_ct <= n_boot_try){
      attempt_ct <- attempt_ct + 1
      boot_idx <- sample(seq_len(n), replace = TRUE)
      boot_data <- data[boot_idx, , drop = FALSE]
      
      rhobar_0_n_boot <- mean(boot_data[[Y_name]][boot_data[[V_name]] == 0])
      rhobar_1_n_boot <- mean(boot_data[[Y_name]][boot_data[[V_name]] == 1])
      
      if(rhobar_0_n_boot > rhobar_1_n_boot){
        success_ct <- success_ct + 1
        boot_diff[[success_ct]] <- get_hudgens_stat_doomed_new(
          data = boot_data, Y_name = Y_name, G_name = G_name, V_name = V_name, family = family
        )
      }
    }
    if(sum(!is.na(boot_diff)) != n_boot){
      warning(paste0(
        "Only achieved ", sum(!is.na(boot_diff)), " successful resamples. Try again with larger n_boot_try? \n sample size: ", nrow(data)
      ))
    }
    
    if(sum(!is.na(boot_diff)) > 2){
      
      boot_diff_res <- data.frame(do.call(rbind, boot_diff))
      
      sd_boot_diff_additive_lower <- sd(unlist(boot_diff_res$additive_effect_lower), na.rm = TRUE)
      sd_boot_diff_additive_upper <- sd(unlist(boot_diff_res$additive_effect_upper), na.rm = TRUE)
      sd_boot_diff_mult_lower <- sd(unlist(boot_diff_res$mult_effect_lower), na.rm = TRUE)
      sd_boot_diff_mult_upper <- sd(unlist(boot_diff_res$mult_effect_upper), na.rm = TRUE)
      
      out <- c(
        out,
        list(
          sd_additive_lower = sd_boot_diff_additive_lower,
          sd_additive_upper = sd_boot_diff_additive_upper,
          sd_mult_lower = sd_boot_diff_mult_lower,
          sd_mult_upper = sd_boot_diff_mult_upper
        )
      )
      
      
      class(out) <- "hudgens_res_doomed_new"
      return(out)
    }else{
      stop("Did not achieve 2 successful bootstrap resamples. Try again with larger n_boot_try? Or give up?")
    }
    
  }else{
    stop("Method not applicable unless evidence of vaccine protection.")
  }
}

#' Function for bootstrap test of Hudgens-style test statistic - doomed
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param G_name growth outcome variable name
#' @param V_name vaccination variable name
#' @param Y_name infection variable name
#' @param lower_bound A boolean. If TRUE, then adds the smallest growth measures 
#'    to the infected vaccinees thereby yielding a lower
#'    bound on the effect of interest. If FALSE, then adds the largest
#'    growth measures to the infected vacccinees thereby yielding an upper
#'    bound on the effect of interest.
#' @param n_boot number of bootstrap replicates
#' @param n_boot_try max number of attempts for bootstrap resampling
#' @param effect_dir direction of beneficial effect, defaults to "positive" for beneficial outcome. Used for one-side tests of bounds.  
#'  
#'  @returns list with observed difference between plc and vax, Hudgens-style test statistic, and p-value
hudgens_test_doomed <- function(
    data, 
    G_name = "G",
    V_name = "V",
    Y_name = "Y",
    lower_bound = TRUE,
    n_boot = 1e3, 
    n_boot_try = n_boot*10,
    effect_dir = "positive"
){
  n_no_inf_plc <- sum(data[[Y_name]] == 0 & data[[V_name]] == 0)
  n_no_inf_vax <- sum(data[[Y_name]] == 0 & data[[V_name]] == 1)
  n_plc <- sum(data[[V_name]] == 0)
  n_vax <- sum(data[[V_name]] == 1)
  n <- n_plc + n_vax
  n_inf_plc <- n_plc - n_no_inf_plc
  n_inf_vax <- n_vax - n_no_inf_vax
  n_diff <- n_inf_plc - n_inf_vax
  
  if(n_diff > 0){
    obs_diff <- get_hudgens_stat_doomed(
      data = data, lower_bound = lower_bound, Y_name = Y_name, G_name = G_name, V_name = V_name
    )
    boot_diff <- rep(NA, n_boot)
    success_ct <- 0
    attempt_ct <- 0
    while(success_ct < n_boot & attempt_ct <= n_boot_try){
      attempt_ct <- attempt_ct + 1
      boot_idx <- sample(seq_len(n), replace = TRUE)
      boot_data <- data[boot_idx, , drop = FALSE]
      
      n_no_inf_plc_boot <- sum(boot_data[[Y_name]] == 0 & boot_data[[V_name]] == 0)
      n_no_inf_vax_boot <- sum(boot_data[[Y_name]] == 0 & boot_data[[V_name]] == 1)
      n_plc_boot <- sum(boot_data[[V_name]] == 0)
      n_vax_boot <- sum(boot_data[[V_name]] == 1)
      n_inf_plc_boot <- n_plc_boot - n_no_inf_plc_boot
      n_inf_vax_boot <- n_vax_boot - n_no_inf_vax_boot
      
      if(n_inf_plc_boot > n_inf_vax_boot){
        success_ct <- success_ct + 1
        boot_diff[success_ct] <- get_hudgens_stat_doomed(boot_data, lower_bound = lower_bound, Y_name = Y_name, G_name = G_name, V_name = V_name)
      }
    }
    if(sum(!is.na(boot_diff)) != n_boot){
      warning(paste0(
        "Only achieved ", sum(!is.na(boot_diff)), " successful resamples. Try again with larger n_boot_try? \n sample size: ", nrow(data)
      ))
    }
    
    if(sum(!is.na(boot_diff)) > 2){
      sd_boot_diff <- sd(boot_diff, na.rm = TRUE)
      test_stat <- obs_diff / sd_boot_diff
      
      if(effect_dir == "positive"){
        pval <- pnorm(test_stat, lower.tail = FALSE)
      } else{
        pval <- pnorm(test_stat, lower.tail = TRUE)
      }
      
      # one-sided p-value
      #pval <- pnorm(test_stat, lower.tail = FALSE)
      
      out <- list(
        obs_diff = obs_diff,
        test_stat = test_stat,
        pval = pval
      )
      
      class(out) <- "hudgens_res_doomed"
      return(out)
    }else{
      stop("Did not achieve 2 successful bootstrap resamples. Try again with larger n_boot_try? Or give up?")
    }
    
  }else{
    stop("Method not applicable unless evidence of vaccine protection.")
  }
}


# -----------------------------------------------------------------------------

# Covariate-adjusted Hudgens Method -------------------------------------------

#' Function for Hudgens-style bounds on effect that incorporate covariates
#' 
#' Currently assumes that the conditional mean of G follows a linear model
#' with Normal errors.
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param models list of pre-fit models needed for estimation
#' @param family gaussian for continuous outcome G, binomial for binary outcome
#' @param lower_bound A boolean. If TRUE, then adds the smallest growth measures 
#'    to the infected vaccines thereby yielding a lower
#'    bound on the effect of interest. If FALSE, then adds the largest
#'    growth measures to the infected vacccinees thereby yielding an upper
#'    bound on the effect of interest.
#' 
#' @examples
#' 
#' n <- 1e4
#' X <- sample(seq(-1,1), n, replace = TRUE)
#' p_immune <- 0.5 + 0.25 * X
#' p_doomed <- 0.1 + 0.05 * X
#' p_helped <- 1 - (p_immune + p_doomed)
#' ps <- mapply(
#'   p_i = p_immune, p_d = p_doomed, p_h = p_helped, 
#'   FUN = function(p_i, p_d, p_h){
#'     sample(
#'       c("immune", "doomed", "helped"), 
#'       size = 1, prob = c(p_i, p_d, p_h)
#'     )
#'   }
#' )
#' 
#' V <- rbinom(n, 1, 0.5)
#' Y0 <- ifelse(ps == "immune", 0, 1)
#' Y1 <- ifelse(ps == "doomed", 1, 0)
#' Y <- ifelse(V == 1, Y1, Y0)
#' G1 <- 1*X - 0.5 * Y1 + rnorm(n, 0, 0.5)
#' G0 <- 1*X - 0.5 * Y0 + rnorm(n, 0, 0.5)
#' G <- ifelse(V == 1, G1, G0)
#' 
#' marginal_effect <- mean(G1 - G0)
#' ps_effect <- mean(G1[Y0 == 1] - G0[Y0 == 1])
#' 
#' data <- data.frame(X, V, Y, G)
#' models <- fit_models(data)
#' 
#' get_adjusted_hudgens_stat(data, models, lower_bound = TRUE)
#' # compare to unadjusted
#' get_hudgens_stat(data, lower_bound = TRUE)
#' 
#' get_adjusted_hudgens_stat(data, models, lower_bound = FALSE)
#' # compare to unadjusted
#' get_hudgens_stat(data, lower_bound = FALSE)
#' 
#' # binary outcome
#' G_binary <- as.numeric(G > 1)
#' data <- data.frame(X, V, Y, G = G_binary)
#' models <- fit_models(data, family = binomial())
#' get_adjusted_hudgens_stat(data, models, family = "binomial", lower_bound = TRUE)
#' get_adjusted_hudgens_stat(data, models, family = "binomial", lower_bound = FALSE)
#' 
#' 
#' @returns Hudgens-style estimate of bound on effect in naturally infected
get_adjusted_hudgens_stat <- function(
    data, 
    models,
    family = "gaussian",
    lower_bound = TRUE
){
  
  E_G_V0_Y1_X <- predict(models$fit_G_V0_Y1_X, newdata = data, type = "response")

  E_G_V1_Y1_X <- predict(models$fit_G_V1_Y1_X, newdata = data, type = "response")
  E_G_V1_Y0_X <- predict(models$fit_G_V1_Y0_X, newdata = data, type = "response")

  P_Y1_V1_X <- predict(models$fit_Y_V1_X, newdata = data, type = "response")
  P_Y1_V0_X <- predict(models$fit_Y_V0_X, newdata = data, type = "response")
  P_Y0_V1_X <- 1 - P_Y1_V1_X
  P_Y0_V0_X <- 1 - P_Y1_V0_X

  P_Y1_V0 <- mean(P_Y1_V0_X)
  
  VE_X <- 1 - ( P_Y1_V1_X / P_Y1_V0_X )
  if(any(VE_X < 0)){
    warning("Some condtional VE estimates < 0 -- truncating these estimates at 0.")
  }
  VE_X[VE_X < 0] <- 0
  VE_is_zero <- (VE_X == 0)
  VE_is_nonzero <- (VE_X > 0)

  q_X_low <- 1 - P_Y0_V0_X / P_Y0_V1_X
  q_X_high <- 1 - q_X_low
  
  if(family == "gaussian"){
    sd_G <- (mean(models$fit_G_V1_Y0_X$residuals^2))^(1/2)
  }

  E_G_V1_Y0_truncG_X <- E_G_V1_Y0_X

  if(lower_bound){
    
    if(family == "gaussian"){
      # calculate mean of Normal given less than q_X_low
      beta_X <- (q_X_low - E_G_V1_Y0_X) / sd_G
      E_G_V1_Y0_truncG_X[VE_is_nonzero] <- E_G_V1_Y0_X[VE_is_nonzero] - sd_G * dnorm(beta_X[VE_is_nonzero]) / pnorm(beta_X[VE_is_nonzero])
    }else{
      E_G_V1_Y0_truncG_X[VE_is_nonzero] <- as.numeric(E_G_V1_Y0_X[VE_is_nonzero] > q_X_low)
    }
  }else{
    
    if(family == "gaussian"){
      # calculate mean of Normal given greater than q_X_high
      alpha_X <- (q_X_high - E_G_V1_Y0_X) / sd_G
      E_G_V1_Y0_truncG_X[VE_is_nonzero] <- E_G_V1_Y0_X[VE_is_nonzero] + sd_G * dnorm(alpha_X[VE_is_nonzero]) / pnorm(alpha_X[VE_is_nonzero], lower.tail = FALSE)
    }else{
      E_G_V1_Y0_truncG_X[VE_is_nonzero] <- as.numeric(E_G_V1_Y0_X[VE_is_nonzero] > q_X_high)
    }
  }
  
  E_G_V1_Y0_X_bound <- E_G_V1_Y1_X * (1 - VE_X) + E_G_V1_Y0_truncG_X * VE_X

  effect <- mean(
    P_Y1_V0_X / P_Y1_V0 * (E_G_V1_Y0_X_bound - E_G_V0_Y1_X)
  )
  
  return(effect)
  
}