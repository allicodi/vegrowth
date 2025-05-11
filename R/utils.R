#' Print the output of a \code{"vegrowth"} object
#' 
#' @param x An \code{"vegrowth"} object.
#' @param ... other arguments (not used)
#' 
#' @method print vegrowth
#' @export
print.vegrowth <- function(x, ...) {
  
  cat("                                        Growth Effect Estimation Results\n")
  cat(paste(rep("-", 120), collapse = ""), "\n")
  
  # Check if any methods with confidence intervals are used
  if (any(sapply(x, function(item) inherits(item, c("gcomp_res", "pop_gcomp_res", "aipw_res", "tmle_res", "hudgens_adj_lower_res", "hudgens_adj_upper_res"))))) {
    col_names <- c("Method", "Point Est.", "SE", "95% CI: Lower", "95% CI: Upper", "Reject")
    
    # Print header with dashed line
    cat(sprintf("%-50s%-15s%-15s%-15s%-15s%-10s\n",
                col_names[1], col_names[2], col_names[3], col_names[4], col_names[5], col_names[6]))
    cat(paste(rep("-", 120), collapse = ""), "\n")
    
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
                    i$pt_est,
                    i$se,
                    i$lower_ci,
                    i$upper_ci,
                    ifelse(i$reject, "Yes", "No")))
      }
    })
    
    
    cat(paste(rep("-", 120), collapse = ""), "\n")
  } 
  
  # If other methods used print separate (chop lump, hudgens)
  if (any(sapply(x, function(item) inherits(item, c("choplump_res", "hudgens_lower_res", "hudgens_upper_res", "hudgens_lower_res_doomed", "hudgens_upper_res_doomed"))))) {
    col_names <- c("Method", "Observed Diff.", "P-Value", "Reject")
    
    # Print header with dashed line
    cat(sprintf("%-50s%-30s%-30s%-15s\n",
                col_names[1], col_names[2], col_names[3], col_names[4]))
    cat(paste(rep("-", 120), collapse = ""), "\n")
    
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
        cat(sprintf("%-50s%-30.4f%-30.4f%-10s\n",
                    method_name,
                    i$obs_diff,
                    i$pval,
                    ifelse(i$reject, "Yes", "No")))
      }
    })
    
  }
  
}
