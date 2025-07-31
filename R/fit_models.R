#' Function to fit all models needed by various estimators
#' 
#' @param data dataset used to fit models
#' @param Y_name name of growth outcome variable, default G
#' @param Z_name name of vaccine treatment variable, default Z
#' @param X_name character vector containing name(s) of covariates, default X
#' @param S_name name of infection variable, default Y
#' @param est character vector of names of estimators to use for growth effect
#' @param Y_X_model optional specify model to be used for fitting growth on covariates, otherwise growth on all covariates
#' @param S_X_model optional specify model to be used for fitting infection on covariates, otherwise infection on all covariates
#' @param family family for outcome model, defaults to gaussian for growth
#' 
#' @export
#' 
#' @returns list with models needed for specified estimators
fit_models <- function(data,
                       Y_name = "Y",
                       Z_name = "Z",
                       X_name = c("X"),
                       S_name = "S",
                       est = c("gcomp_pop_estimand", 
                               "gcomp",
                               "ipw",
                               "efficient_aipw",
                               "efficient_aipw_sens",
                               "efficient_tmle"),
                       Y_Z_X_model = NULL,
                       Y_X_S1_model = NULL,
                       Y_X_S0_model = NULL,
                       S_X_model = NULL,
                       Z_X_model = paste0(Z_name, " ~ 1"),
                       family = "gaussian"){
  
  # Prep model formulas if not pre-specified
  if(is.null(Y_Z_X_model)){
    Y_Z_X_model <- as.formula(paste0(Y_name, "~",
                                     Z_name, "+",
                                     paste0(X_name, collapse = "+")))
  }
  
  if(is.null(Y_X_S1_model)){
    Y_X_S1_model <- as.formula(paste0(Y_name, "~",
                                      paste0(X_name, collapse = "+")))
  }
  
  if(is.null(Y_X_S0_model)){
    Y_X_S0_model <- as.formula(paste0(Y_name, "~",
                                      paste0(X_name, collapse = "+")))
  }
  
  if(is.null(S_X_model)){
    S_X_model <- as.formula(paste0(S_name, "~",
                                   paste0(X_name, collapse = "+")))
  }
  
  out <- list()
  
  # only needed for population estimator
  if("gcomp_pop_estimand" %in% est){
    out$fit_Y_Z_X <- glm(Y_Z_X_model, data = data, family = family)
  } 
  
  # needed for any weight-based estimator
  if(any(c("gcomp", "efficient_aipw", "efficient_tmle", "efficient_aipw_sens") %in% est)){
    sub_Z0 <- data[data[[Z_name]] == 0,]
    out$fit_S_Z0_X <- glm(S_X_model, sub_Z0, family = "binomial")
    
    sub_Z1 <- data[data[[Z_name]] == 1,]
    out$fit_S_Z1_X <- glm(S_X_model, sub_Z1, family = "binomial")
    
    sub_Z1_S1 <- data[data[[Z_name]] == 1 & data[[S_name]] == 1,]
    out$fit_Y_Z1_S1_X <- glm(Y_X_S1_model, data = sub_Z1_S1, family = family)
    
    sub_Z1_S0 <- data[data[[Z_name]] == 1 & data[[S_name]] == 0,]
    out$fit_Y_Z1_S0_X <- glm(Y_X_S0_model, data = sub_Z1_S0, family = family)
    
    sub_Z0_S1 <- data[data[[Z_name]] == 0 & data[[S_name]] == 1,]
    out$fit_Y_Z0_S1_X <- glm(Y_X_S1_model, data = sub_Z0_S1, family = family)
    
  }
  
  if(any(c("efficient_aipw_sens") %in% est)){
    out$fit_S_Z_X <- glm(
      paste0(S_name, "~", Z_name, "+", paste0(X_name, collapse = "+")),
      family = "binomial",
      data = data
    )
  }
  
  if(any(c("efficient_aipw_sens", "efficient_aipw", "ipw", "efficient_tmle") %in% est)){
    out$fit_Z_X <- glm(
      Z_X_model, family = "binomial", data = data
    )
  }
  return(out)
  
}


#' Function to fit all SuperLearner models needed by various estimators
#' 
#' @param data dataset used to fit models
#' @param Y_name name of growth outcome variable, default G
#' @param Z_name name of vaccine treatment variable, default Z
#' @param X_name character vector containing name(s) of covariates, default X
#' @param S_name name of infection variable, default Y
#' @param est character vector of names of estimators to use for growth effect
#' @param Y_X_library specify SuperLearner libraries to be used for fitting growth on covariates. Default GLM 
#' @param S_X_library specify SuperLearner libraries to be used for fitting infection on covariates. Default GLM 
#' @param family family for outcome model, defaults to gaussian for growth
#' @param v_folds number of cross-validation folds to use in SuperLearner, default 3
#' 
#' @export
#' 
#' @returns list with models needed for specified estimators
fit_ml_models <- function(data,
                          Y_name = "Y",
                          Z_name = "Z",
                          X_name = c("X"),
                          S_name = "S",
                          est = c("gcomp_pop_estimand", 
                                  "gcomp",
                                  "ipw",
                                  "efficient_aipw", 
                                  "efficient_tmle"),
                          Y_Z_X_library = c("SL.glm"),
                          Y_X_library = c("SL.glm"),
                          S_X_library = c("SL.glm"),
                          Z_X_library = c("SL.mean"),
                          family = "gaussian",
                          v_folds = 3){
  
  out <- list()
  
  # only needed for population estimator
  if("gcomp_pop_estimand" %in% est){
    out$fit_Y_Z_X <- SuperLearner::SuperLearner(Y = data[[Y_name]],
                                                X = data[, colnames(data) %in% c(Z_name, X_name), drop = FALSE],
                                                family = family,
                                                SL.library = Y_Z_X_library, 
                                                cvControl = list(Z = v_folds))
    
  } 
  
  # needed for any weight-based estimator
  if(any(c("gcomp", "efficient_aipw", "efficient_tmle") %in% est)){
    sub_Z0 <- data[data[[Z_name]] == 0,]
    out$fit_S_Z0_X <- SuperLearner::SuperLearner(Y = sub_Z0[[S_name]],
                                                 X = sub_Z0[, X_name, drop = FALSE],
                                                 family = stats::binomial(),
                                                 SL.library = S_X_library, 
                                                 cvControl = list(Z = v_folds))
    
    sub_Z1 <- data[data[[Z_name]] == 1,]
    out$fit_S_Z1_X <- SuperLearner::SuperLearner(Y = sub_Z1[[S_name]],
                                                 X = sub_Z1[, X_name, drop = FALSE],
                                                 family = stats::binomial(),
                                                 SL.library = S_X_library, 
                                                 cvControl = list(Z = v_folds))
    
    sub_Z1_S1 <- data[data[[Z_name]] == 1 & data[[S_name]] == 1,]
    out$fit_Y_Z1_S1_X <- SuperLearner::SuperLearner(Y = sub_Z1_S1[[Y_name]],
                                                    X = sub_Z1_S1[, X_name, drop = FALSE],
                                                    family = family,
                                                    SL.library = Y_X_library, 
                                                    cvControl = list(Z = v_folds))
    
    sub_Z1_S0 <- data[data[[Z_name]] == 1 & data[[S_name]] == 0,]
    out$fit_Y_Z1_S0_X <- SuperLearner::SuperLearner(Y = sub_Z1_S0[[Y_name]],
                                                    X = sub_Z1_S0[, X_name, drop = FALSE],
                                                    family = family,
                                                    SL.library = Y_X_library, 
                                                    cvControl = list(Z = v_folds))
    
    sub_Z0_S1 <- data[data[[Z_name]] == 0 & data[[S_name]] == 1,]
    out$fit_Y_Z0_S1_X <- SuperLearner::SuperLearner(Y = sub_Z0_S1[[Y_name]],
                                                    X = sub_Z0_S1[, X_name, drop = FALSE],
                                                    family = family,
                                                    SL.library = Y_X_library, 
                                                    cvControl = list(Z = v_folds))
    
  }
  
  if(any(c("ipw", "efficient_aipw", "efficient_tmle") %in% est)){
    out$fit_Z_X <- SuperLearner::SuperLearner(Y = data[[Z_name]],
                                              X = data[, X_name, drop = FALSE],
                                              family = binomial(),
                                              SL.library = Z_X_library, 
                                              cvControl = list(Z = v_folds))
  }
  
  return(out)
  
}
