
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
      } else if(class(row) == "sens"){
        class(row) <- "data.frame"
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
             lower_ci_mult = exp(quantile(boot_df$log_multiplicative_effect, 0.025, na.rm = TRUE)), # added NA rm true for some bootstrap replicates
             upper_ci_mult = exp(quantile(boot_df$log_multiplicative_effect, 0.975, na.rm = TRUE)),
             se_psi_1 = sd(boot_df$psi_1),
             lower_ci_psi_1 = quantile(boot_df$psi_1, 0.025),
             upper_ci_psi_1 = quantile(boot_df$psi_1, 0.975),
             se_psi_0 = sd(boot_df$psi_0),
             lower_ci_psi_0 = quantile(boot_df$psi_0, 0.025),
             upper_ci_psi_0 = quantile(boot_df$psi_0, 0.975))
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
             #original
             #se_mult_lower = sd(boot_df$mult_effect_lower),
             #new
             se_log_mult_lower = sd(log(boot_df$mult_effect_lower)),
             lower_ci_mult_lower = quantile(boot_df$mult_effect_lower, 0.025),
             upper_ci_mult_lower = quantile(boot_df$mult_effect_lower, 0.975),
             #original
             #se_mult_upper = sd(boot_df$mult_effect_upper),
             #new
             se_log_mult_upper = sd(log(boot_df$mult_effect_upper)),
             lower_ci_mult_upper = quantile(boot_df$mult_effect_upper, 0.025, na.rm = TRUE),
             upper_ci_mult_upper = quantile(boot_df$mult_effect_upper, 0.975, na.rm = TRUE))
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
  boot_res_list <- vector("list", length = length(epsilon))
  names(boot_res_list) <- paste0("epsilon_", epsilon)
  
  for(e in 1:length(epsilon)){
    # QUESTION handling NAs 
    boot_res_list[[e]] <- data.frame(epsilon = epsilon[e],
                                     se_additive = sd(boot_df$additive_effect[boot_df$epsilon == epsilon[e]], na.rm = TRUE),
                                     lower_ci_additive = quantile(boot_df$additive_effect[boot_df$epsilon == epsilon[e]], 0.025, na.rm = TRUE),
                                     upper_ci_additive = quantile(boot_df$additive_effect[boot_df$epsilon == epsilon[e]], 0.975, na.rm = TRUE),
                                     se_mult = sd(boot_df$log_multiplicative_effect[boot_df$epsilon == epsilon[e]], na.rm = TRUE),
                                     lower_ci_mult = exp(quantile(boot_df$log_multiplicative_effect[boot_df$epsilon == epsilon[e]], 0.025, na.rm = TRUE)),
                                     upper_ci_mult = exp(quantile(boot_df$log_multiplicative_effect[boot_df$epsilon == epsilon[e]], 0.975, na.rm = TRUE)))
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
#' @param type format to print results in, group by 'estimand' (default) or group by 'effect' (additive and multiplicative) 
#' @param ... other arguments (not used)
#' 
#' @method print vegrowth
#' @export
print.vegrowth <- function(x, scale = "additive", ...) {
  
  # Helper to print one row
  print_row <- function(label, method, est, lower, upper) {
    cat(sprintf("%-25s%-15s%-15.4f%-15.4f%-15.4f\n", label, method, est, lower, upper))
  }
  
  # Helper to get estimates
  extract_estimates <- function(est_obj, method, scale, label_prefix) {
    if (scale == "additive") {
      if (is.null(est_obj$boot_se)) {
        est <- est_obj$pt_est['additive_effect']
        se <- est_obj$pt_est['additive_se']
        print_row(label_prefix, method, est, est - 1.96 * se, est + 1.96 * se)
      } else {
        est <- est_obj$pt_est['additive_effect']
        lower <- est_obj$boot_se$lower_ci_additive
        upper <- est_obj$boot_se$upper_ci_additive
        print_row(label_prefix, method, est, lower, upper)
      }
    } else {
      if (is.null(est_obj$boot_se)) {
        est <- exp(est_obj$pt_est['log_multiplicative_effect'])
        se <- est_obj$pt_est['log_multiplicative_se']
        print_row(label_prefix, method, est,
                  exp(log(est) - 1.96 * se),
                  exp(log(est) + 1.96 * se))
      } else {
        est <- exp(est_obj$pt_est['log_multiplicative_effect'])
        lower <- est_obj$boot_se$lower_ci_mult
        upper <- est_obj$boot_se$upper_ci_mult
        print_row(label_prefix, method, est, lower, upper)
      }
    }
  }
  
  # Helper to print bounds
  print_bounds <- function(est_obj, scale, label_prefix) {
    if (scale == "additive") {
      print_row(label_prefix, "Lower Bound",
                est_obj$pt_est['additive_effect_lower'],
                est_obj$boot_se$lower_ci_additive_lower,
                est_obj$boot_se$upper_ci_additive_lower)
      print_row(label_prefix, "Upper Bound",
                est_obj$pt_est['additive_effect_upper'],
                est_obj$boot_se$lower_ci_additive_upper,
                est_obj$boot_se$upper_ci_additive_upper)
    } else {
      print_row(label_prefix, "Lower Bound",
                est_obj$pt_est['mult_effect_lower'],
                est_obj$boot_se$lower_ci_mult_lower,
                est_obj$boot_se$upper_ci_mult_lower)
      print_row(label_prefix, "Upper Bound",
                est_obj$pt_est['mult_effect_upper'],
                est_obj$boot_se$lower_ci_mult_upper,
                est_obj$boot_se$upper_ci_mult_upper)
    }
  }
  
  # Header
  scale_label <- ifelse(scale == "additive", "Additive", "Multiplicative")
  cat(sprintf("%50s\n", paste("                         Growth Effect Estimation Results:", scale_label)))
  cat(paste(rep("-", 90), collapse = ""), "\n")
  col_names <- c("Estimand", "Method", "Point Est.", "95% CI: Lower", "95% CI: Upper")
  cat(sprintf("%-25s%-15s%-15s%-15s%-15s\n", col_names[1], col_names[2], col_names[3], col_names[4], col_names[5]))
  cat(paste(rep("-", 90), collapse = ""), "\n")
  
  # Loop through estimands
  lapply(x, function(i) {
    if (inherits(i, "nat_inf")) {
      cat(sprintf("%-25s%-15s%-15s%-15s%-15s\n", "Naturally Infected - - - ", "- - - - - - - -", " - - - - - - - ", "- - - - - - - -", " - - - - - - - - - - "))
      lapply(i, function(j) {
        if (inherits(j, "gcomp")) extract_estimates(j, "G-Computation", scale, "")
        if (inherits(j, "ipw")) extract_estimates(j, "IPW", scale, "")
        if (inherits(j, "aipw")) extract_estimates(j, "AIPW", scale, "")
        if (inherits(j, "tmle")) extract_estimates(j, "TMLE", scale, "")
        if (inherits(j, "bound")) print_bounds(j, scale, "")
      })
    }
    
    if (inherits(i, "doomed")) {
      cat(sprintf("%-25s%-15s%-15s%-15s%-15s\n", "Doomed - - - - - - - - -", "- - - - - - - -", " - - - - - - - ", "- - - - - - - -", " - - - - - - - - - - "))
      lapply(i, function(j) {
        if (inherits(j, "gcomp")) extract_estimates(j, "G-Computation", scale, "")
        if (inherits(j, "ipw")) extract_estimates(j, "IPW", scale, "")
        if (inherits(j, "aipw")) extract_estimates(j, "AIPW", scale, "")
        if (inherits(j, "tmle")) extract_estimates(j, "TMLE", scale, "")
        if (inherits(j, "bound")) print_bounds(j, scale, "")
      })
    }
    
    if (inherits(i, "pop")) {
      cat(sprintf("%-25s%-15s%-15s%-15s%-15s\n", "Population - - - - - - - ", "- - - - - - - -", " - - - - - - - ", "- - - - - - - -", " - - - - - - - - - - "))
      lapply(i, function(j) {
        if (inherits(j, "gcomp")) extract_estimates(j, "G-Computation", scale, "")
        if (inherits(j, "ipw")) extract_estimates(j, "IPW", scale, "")
        if (inherits(j, "aipw")) extract_estimates(j, "AIPW", scale, "")
        if (inherits(j, "tmle")) extract_estimates(j, "TMLE", scale, "")
      })
    }
  })
  
  invisible(x)
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
