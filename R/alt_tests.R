
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
#'    growth measures to the infected vacccinees thereby yielding an upper
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
  
  if(n_inf_plc > n_inf_vax){
    
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
  
  return(mean_G_vax - mean_G_plc)
  
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
#'  
#'  @returns list with observed difference between plc and vax, Hudgens-style test statistic, and p-value
hudgens_test <- function(
    data, 
    G_name = "G",
    V_name = "V",
    Y_name = "Y",
    lower_bound = TRUE,
    n_boot = 1e3, 
    n_boot_try = n_boot*10
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
      # one-sided p-value
      pval <- pnorm(test_stat, lower.tail = FALSE)
      
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
