
#' Function for permutation test of bounds in naturally infected
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param Y_name growth outcome variable name
#' @param Z_name vaccination variable name
#' @param S_name infection variable name
#' @param n_permutations number of permutations to run, default 1000
#' @param family defaults to gaussian for continuous outcome
#' @param effect_dir direction of beneficial effect, defaults to "positive" for beneficial outcome. Used for one-side tests of bounds.  
#'  
#' @returns list with original estimate, null estimates for n_permutations, and p_value for one-sided test based on effect_dir
permutation_bound_nat_inf <- function(
    data, 
    Y_name = "Y",
    Z_name = "Z",
    S_name = "S",
    n_permutations = 1e3, 
    family = "gaussian",
    effect_dir = "positive"
){
  
  original_est <- get_bound_nat_inf(
    data = data, S_name = S_name, Y_name = Y_name, Z_name = Z_name, family = family
  )
  
  ## Permutation approach
  null_est <- vector("list", length = n_permutations)
  for(i in 1:n_permutations){
    data_shuffle <- data
    data_shuffle[[Z_name]] <- sample(data_shuffle[[Z_name]])
    
    null_est[[i]] <- get_bound_nat_inf(
      data = data_shuffle, S_name = S_name, Y_name = Y_name, Z_name = Z_name, family = family
    )
  }
  
  null_est <- do.call(rbind, null_est)
  
  out <- list(
    original_est = original_est, 
    null_est = null_est,
    pval_bound = ifelse(effect_dir == "negative",
                        mean(null_est$additive_effect_upper < original_est$additive_effect_upper), # negative effect means we are interested in the upper bound < 0
                        mean(null_est$additive_effect_lower > original_est$additive_effect_lower)  # positive effect means we are interested in the lower bound > 0 
    )
  )
    
    class(out) <- "permutation_bonud_nat_inf"
    return(out)
    
}

#' Function for permutation test of bounds in doomed
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param Y_name growth outcome variable name
#' @param Z_name vaccination variable name
#' @param S_name infection variable name
#' @param n_permutations number of permutations to run, default 1000
#' @param family defaults to gaussian for continuous outcome
#' @param effect_dir direction of beneficial effect, defaults to "positive" for beneficial outcome. Used for one-side tests of bounds.  
#'  
#' @returns list with original estimate, null estimates for n_permutations, and p_value for one-sided test based on effect_dir
permutation_bound_doomed <- function(
    data, 
    Y_name = "Y",
    Z_name = "Z",
    S_name = "S",
    lower_bound = TRUE,
    n_permutations = 1e3, 
    n_boot_try = n_boot*10,
    effect_dir = "positive"
){
  
  original_est <- get_bound_doomed(
    data = data, S_name = S_name, Y_name = Y_name, Z_name = Z_name, family = family
  )
  
  ## Permutation approach
  null_est <- vector("list", length = n_permutations)
  for(i in 1:n_permutations){
    data_shuffle <- data
    data_shuffle[[Z_name]] <- sample(data_shuffle[[Z_name]])
    
    null_est[[i]] <- get_bound_doomed(
      data = data_shuffle, S_name = S_name, Y_name = Y_name, Z_name = Z_name, family = family
    )
  }
  
  null_est <- do.call(rbind, null_est)
  
  out <- list(
    original_est = original_est, 
    null_est = null_est,
    pval_bound = ifelse(effect_dir == "negative",
                        mean(null_est$additive_effect_upper < original_est$additive_effect_upper), # negative effect means we are interested in the upper bound < 0
                        mean(null_est$additive_effect_lower > original_est$additive_effect_lower)  # positive effect means we are interested in the lower bound > 0 
    )
  )
    
    class(out) <- "permutation_bonud_doomed"
    return(out)
    
}