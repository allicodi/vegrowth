
#' Helper function to make dataframe from bootstrap results
#'
#' @param boot_estimates list of n_boot bootstrap estimates from one_boot
#' @param estimand estimand of interest
#' @param method method of interest
#' 
#' @returns data frame with bootstrap results for point estimates & bounds
make_boot_df <- function(boot_estimates, estimand = "nat_inf", method = "gcomp") {
  res_list <- lapply(seq_along(boot_estimates), function(i) {
    boot_res <- boot_estimates[[i]]
    
    # Check if the desired estimand/method exists in this iteration
    if (!is.null(boot_res[[estimand]]) && !is.null(boot_res[[estimand]][[method]])) {
      row <- boot_res[[estimand]][[method]]
      
      # If it's just a single number, convert to named data frame
      if (is.atomic(row) && length(row) == 1) {
        row <- data.frame(estimate = row)
      } else if (is.null(dim(row))) {
        row <- as.data.frame(t(row))
      }
      
      row$boot_id <- i
      row$estimand <- estimand
      row$method <- method
      return(row)
    } else {
      return(NULL)  # Skip missing or failed bootstraps
    }
  })
  # Combine all into one data frame
  data.frame(do.call(rbind, res_list))
}

#' Helper function to get SE and 95% CI for additive and multiplicative effects - point estimates
#'
#' @param estimand estimand of interest
#' @param method method of interest
#' 
#' @returns dataframe with bootstrap standard error and confidence intervals on additive and multiplicative scales for point estimates
get_boot_se <- function(boot_estimates, estimand = "nat_inf", method = "gcomp"){
  boot_df <- make_boot_df(boot_estimates = boot_estimates,
                          estimand = estimand,
                          method = method)
  
  data.frame(se_additive = sd(boot_df$additive_effect),
             lower_ci_additive = quantile(boot_df$additive_effect, 0.025),
             upper_ci_additive = quantile(boot_df$additive_effect, 0.975),
             se_log_mult = sd(boot_df$log_multiplicative_effect),
             lower_ci_mult = exp(quantile(boot_df$log_multiplicative_effect, 0.025)),
             upper_ci_mult = exp(quantile(boot_df$log_multiplicative_effect, 0.975)))
}

#' Helper function to get SE and 95% CI for additive and multiplicative effects - bounds
#'
#' @param estimand estimand of interest
#' @param method method of interest
#' 
#' @returns dataframe with bootstrap standard error and confidence intervals on additive and multiplicative scales for bounds
get_boot_se_bound <- function(boot_estimates, estimand = "nat_inf", method = "bound"){
  boot_df <- make_boot_df(boot_estimates = boot_estimates,
                          estimand = estimand,
                          method = method)
  
  data.frame(se_additive_lower = sd(boot_df$additive_effect_lower),
             lower_ci_additive_lower = quantile(boot_df$additive_effect_lower, 0.025),
             upper_ci_additive_lower = quantile(boot_df$additive_effect_lower, 0.975),
             se_additive_upper = sd(boot_df$additive_effect_upper),
             lower_ci_additive_upper = quantile(boot_df$additive_effect_upper, 0.025),
             upper_ci_additive_upper = quantile(boot_df$additive_effect_upper, 0.975),
             # CHECK IF THIS NEEDS TO BE EXPONENTIATED-- I THINK NO?
             se_mult_lower = sd(boot_df$mult_effect_lower),
             lower_ci_mult_lower = quantile(boot_df$mult_effect_lower, 0.025),
             upper_ci_mult_lower = quantile(boot_df$mult_effect_lower, 0.975),
             se_mult_upper = sd(boot_df$mult_effect_upper),
             lower_ci_mult_upper = quantile(boot_df$mult_effect_upper, 0.025),
             upper_ci_mult_upper = quantile(boot_df$mult_effect_upper, 0.975))
}

#' Helper function to get SE and 95% CI for additive and multiplicative effects - bounds
#'
#' @param estimand estimand of interest
#' @param method method of interest
#' 
#' @returns dataframe with bootstrap standard error and confidence intervals on additive and multiplicative scales for each epsilon
get_boot_se_sens <- function(boot_estimates, estimand = "nat_inf", method = "sens"){
  boot_df <- make_boot_df(boot_estimates = boot_estimates,
                          estimand = estimand,
                          method = method)
  
  epsilon <- unique(boot_df$epsilon)
  boot_res_list <- vector(list, length = length(epsilon))
  names(boot_res_list) <- paste0("epsilon_", epsilon)
  
  for(e in 1:length(epsilon)){
    boot_res_list[[e]] <- data.frame(epsilon = epsilon[e],
                                     se_additive = sd(boot_df$additive_effect[boot_df$epsilon == epsilon[e]]),
                                     lower_ci_additive = quantile(boot_df$additive_effect[boot_df$epsilon == epsilon[e]], 0.025),
                                     upper_ci_additive = quantile(boot_df$additive_effect[boot_df$epsilon == epsilon[e]], 0.975),
                                     se_mult = sd(boot_df$log_multiplicative_effect[boot_df$epsilon == epsilon[e]]),
                                     lower_ci_mult = exp(quantile(boot_df$log_multiplicative_effect[boot_df$epsilon == epsilon[e]], 0.025)),
                                     upper_ci_mult = exp(quantile(boot_df$log_multiplicative_effect[boot_df$epsilon == epsilon[e]], 0.975)))
  }
  
  return(do.call(rbind, boot_res_list))
  
}

#' Example simulated dataset based on the PROVIDE study
#' 
#' @format A data frame with X rows and Y variables:
#' \describe{
#'   \item{wk10_haz}{Height-for-age Z-score at 10 weeks of age (numeric)}
#'   \item{gender}{Infant gender, either `"Male"` or `"Female"` (character)}
#'   \item{num_hh_sleep}{Number of people sleeping in the household (integer)}
#'   \item{rotaarm}{Indicator for assignment to rotavirus vaccination arm (0 = control, 1 = vaccine) (integer)}
#'   \item{rotaepi}{Indicator for rotavirus episode during follow-up (0 = no, 1 = yes) (integer)}
#'   \item{any_abx_wk52}{Any antibiotic usage reported through 52 weeks of age (0 = no, 1 = yes) (integer)}
#' }
"provide"

#' Print the output of a \code{"vegrowth"} object
#' 
#' @param x An \code{"vegrowth"} object.
#' @param ... other arguments (not used)
#' 
#' @method print vegrowth
#' @export
print.vegrowth <- function(x, ...) {
  
  cat("                                        Growth Effect Estimation Results: Additive\n")
  cat(paste(rep("-", 130), collapse = ""), "\n")
  
  # Check if any methods with confidence intervals are used
  if (any(sapply(x, function(item) inherits(item, c("gcomp_res", "pop_gcomp_res", "aipw_res", "tmle_res", "hudgens_adj_lower_res", "hudgens_adj_upper_res"))))) {
    col_names <- c("Method", "Point Est.", "SE", "95% CI: Lower", "95% CI: Upper", "Reject (2-sided)")
    
    # Print header with dashed line
    cat(sprintf("%-50s%-15s%-15s%-15s%-15s%-10s\n",
                col_names[1], col_names[2], col_names[3], col_names[4], col_names[5], col_names[6]))
    cat(paste(rep("-", 130), collapse = ""), "\n")
    
    # Iterate through objects and print their results
    lapply(x, function(i) {
      if (inherits(i, c("gcomp_res", "pop_gcomp_res", "aipw_res", "tmle_res", "hudgens_adj_lower_res", "hudgens_adj_upper_res"))) {
        # Determine method name
        method_name <- switch(class(i)[1],
                              "gcomp_res" = "G-Comp VE Estimand",
                              "pop_gcomp_res" = "G-Comp Population Estimand",
                              "aipw_res" = "AIPW Estimand",
                              "tmle_res" = "TMLE Estimand",
                              "hudgens_adj_lower_res" = "Covariate-adjusted Hudgens: Lower Bound",
                              "hudgens_adj_upper_res" = "Covariate-adjusted Hudgens: Upper Bound")
        
        # Print the results
        cat(sprintf("%-50s%-15.4f%-15.4f%-15.4f%-15.4f%-10s\n",
                    method_name,
                    i$pt_est_additive,
                    i$se_additive,
                    i$lower_ci_additive,
                    i$upper_ci_additive,
                    ifelse(i$reject_additive, "Yes", "No")))
      }
    })
    
    
    cat(paste(rep("-", 130), collapse = ""), "\n")
  } 
  
  # If other methods used print separate (chop lump, hudgens)
  if (any(sapply(x, function(item) inherits(item, c("choplump_res", "hudgens_lower_res", "hudgens_upper_res", "hudgens_lower_res_doomed", "hudgens_upper_res_doomed"))))) {
    col_names <- c("Method", "Observed Diff.", "P-Value", "95% CI: Lower", "95% CI: Upper", "Reject (1-sided)")
    
    # Print header with dashed line
    cat(sprintf("%-50s%-15s%-15s%-15s%-15s%-10s\n",
                col_names[1], col_names[2], col_names[3], col_names[4], col_names[5], col_names[6]))
    cat(paste(rep("-", 130), collapse = ""), "\n")
    
    # Iterate through objects and print their results
    lapply(x, function(i) {
      if (inherits(i, c("choplump_res", "hudgens_lower_res", "hudgens_upper_res", "hudgens_lower_res_doomed", "hudgens_upper_res_doomed"))) {
        # Determine method name
        method_name <- switch(class(i)[1],
                              "choplump_res" = "Chop-Lump",
                              "hudgens_lower_res" = "Naturally Infected: Lower Bound",
                              "hudgens_upper_res" = "Naturally Infected: Upper Bound",
                              "hudgens_lower_res_doomed" = "Doomed: Lower Bound",
                              "hudgens_upper_res_doomed" = "Doomed: Upper Bound")
        
        # Print the results
        if( inherits(i, c("choplump_res"))){
          cat(sprintf("%-50s%-15.4f%-15.4f%-15.4f%-15.4f%-10s\n",
                      method_name,
                      i$obs_diff,
                      i$pval,
                      NA, 
                      NA,
                      ifelse(i$reject, "Yes", "No")))
        } else{
          cat(sprintf("%-50s%-15.4f%-15.4f%-15.4f%-15.4f%-10s\n",
                      method_name,
                      i$obs_diff,
                      i$pval,
                      i$lower_ci, 
                      i$upper_ci,
                      ifelse(i$reject, "Yes", "No")))
        }
        
      }
    })
    
  }
  
  # REPEAT FOR MULTIPLICATIVE EFFECTS:
  cat("\n\n")
  cat("                                        Growth Effect Estimation Results: Multiplicative\n")
  cat(paste(rep("-", 130), collapse = ""), "\n")
  
  # Check if any methods with confidence intervals are used
  if (any(sapply(x, function(item) inherits(item, c("gcomp_res", "pop_gcomp_res", "aipw_res", "tmle_res", "hudgens_adj_lower_res", "hudgens_adj_upper_res"))))) {
    col_names <- c("Method", "Point Est.", "SE", "95% CI: Lower", "95% CI: Upper", "Reject (2-sided)")
    
    # Print header with dashed line
    cat(sprintf("%-50s%-15s%-15s%-15s%-15s%-10s\n",
                col_names[1], col_names[2], col_names[3], col_names[4], col_names[5], col_names[6]))
    cat(paste(rep("-", 130), collapse = ""), "\n")
    
    # Iterate through objects and print their results
    lapply(x, function(i) {
      if (inherits(i, c("gcomp_res", "pop_gcomp_res", "aipw_res", "tmle_res", "hudgens_adj_lower_res", "hudgens_adj_upper_res"))) {
        # Determine method name
        method_name <- switch(class(i)[1],
                              "gcomp_res" = "G-Comp VE Estimand",
                              "pop_gcomp_res" = "G-Comp Population Estimand",
                              "aipw_res" = "AIPW Estimand",
                              "tmle_res" = "TMLE Estimand",
                              "hudgens_adj_lower_res" = "Covariate-adjusted Hudgens: Lower Bound",
                              "hudgens_adj_upper_res" = "Covariate-adjusted Hudgens: Upper Bound")
        
        # Print the results
        cat(sprintf("%-50s%-15.4f%-15.4f%-15.4f%-15.4f%-10s\n",
                    method_name,
                    i$pt_est_mult,
                    i$se_log_mult,
                    i$lower_ci_mult,
                    i$upper_ci_mult,
                    ifelse(i$reject_mult, "Yes", "No")))
      }
    })
    
    
    cat(paste(rep("-", 130), collapse = ""), "\n")
  } 
  
}

#' Plot method for sens objects
#'
#' @param object An object of class "sens"
#' @param se Logical; whether to include error bars
#' @param effect_type Character; either "additive" or "multiplicative"
#' @param ... Additional arguments passed to plot (not used here)
#' 
#' @export
plot.sens <- function(
  object, se = TRUE, effect_type = c("additive", "multiplicative"), 
  ...
) {
  effect_type <- match.arg(effect_type)

  if (!inherits(object, "sens")) {
    stop("Object must be of class 'sens'")
  }

  df <- object
  class(df) <- "data.frame"

  if (effect_type == "additive") {
    df$effect <- df$additive_effect
    if (se) {
      df$lower <- df$additive_effect - 1.96 * df$additive_se
      df$upper <- df$additive_effect + 1.96 * df$additive_se
    }
    ylab <- "Additive Effect"
  } else {
    df$effect <- exp(df$log_multiplicative_effect)
    if (se) {
      df$lower <- exp(df$log_multiplicative_effect - 1.96 * df$log_multiplicative_se)
      df$upper <- exp(df$log_multiplicative_effect + 1.96 * df$log_multiplicative_se)
    }
    ylab <- "Multiplicative Effect"
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = epsilon, y = effect)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = expression(epsilon), y = ylab) +
    ggplot2::theme_minimal()

  if (se) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.2)
  }

  print(p)
}
