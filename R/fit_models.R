#' Function to fit all models needed by various estimators
#' 
#' @param data dataset used to fit models
#' @param Y_name name of growth outcome variable, default G
#' @param Z_name name of vaccine treatment variable, default Z
#' @param X_name character vector containing name(s) of covariates, default X
#' @param S_name name of infection variable, default Y
#' @param estimand character vector with name(s) of estimands of interest; "nat_inf" = naturally infected, "doomed" = doomed, "pop" = marginal/population-level
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
                       estimand = c("nat_inf", "doomed", "pop"),
                       method = c("gcomp", "ipw", "aipw", "tmle", "bound", "sens"),
                       exclusion_restriction = c(TRUE, FALSE),
                       Y_Z_X_model = NULL,
                       Y_X_S1_model = NULL,
                       Y_X_S0_model = NULL,
                       S_X_model = NULL,
                       S_Z_X_model = NULL,
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
  
  if(
    "pop" %in% estimand | 
    ("nat_inf" %in% estimand & TRUE %in% exclusion_restriction & any(c("gcomp", "aipw", "tmle") %in% method))
  ){
    out$fit_Y_Z_X <- glm(Y_Z_X_model, data = data, family = family)
  } 
  
  if(any(c("nat_inf", "doomed") %in% estimand)){
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
  
  # only needed for AIPW sensitivity analysis
  if(any(c("sens") %in% method)){
    if(is.null(S_Z_X_model)){
      S_Z_X_model <- paste0(S_name, "~", Z_name, "+", paste0(X_name, collapse = "+"))
    }
    
    out$fit_S_Z_X <- glm(
      S_Z_X_model,
      family = "binomial",
      data = data
    )
  }
  
  # not needed for gcomp
  if(any(c("aipw", "sens", "ipw", "tmle") %in% method)){
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
                          estimand = c("nat_inf", "doomed", "pop"),
                          method = c("gcomp", "ipw", "aipw", "tmle", "bound", "sens"),
                          Y_Z_X_library = c("SL.glm"),
                          Y_X_library = c("SL.glm"),
                          S_X_library = c("SL.glm"),
                          S_Z_X_library = c("SL.glm"),
                          Z_X_library = c("SL.mean"),
                          family = "gaussian",
                          v_folds = 3){
  
  out <- list()
  
  # only needed for population estimator
  if("pop" %in% estimand){
    out$fit_Y_Z_X <- SuperLearner::SuperLearner(Y = data[[Y_name]],
                                                X = data[, colnames(data) %in% c(Z_name, X_name), drop = FALSE],
                                                family = family,
                                                SL.library = Y_Z_X_library, 
                                                cvControl = list(V = v_folds))
    
  } 
  
  # needed for any weight-based estimator
  if(any(c("nat_inf", "doomed") %in% estimand)){
    sub_Z0 <- data[data[[Z_name]] == 0,]
    out$fit_S_Z0_X <- SuperLearner::SuperLearner(Y = sub_Z0[[S_name]],
                                                 X = sub_Z0[, X_name, drop = FALSE],
                                                 family = stats::binomial(),
                                                 SL.library = S_X_library, 
                                                 cvControl = list(V = v_folds))
    
    sub_Z1 <- data[data[[Z_name]] == 1,]
    out$fit_S_Z1_X <- SuperLearner::SuperLearner(Y = sub_Z1[[S_name]],
                                                 X = sub_Z1[, X_name, drop = FALSE],
                                                 family = stats::binomial(),
                                                 SL.library = S_X_library, 
                                                 cvControl = list(V = v_folds))
    
    sub_Z1_S1 <- data[data[[Z_name]] == 1 & data[[S_name]] == 1,]
    out$fit_Y_Z1_S1_X <- SuperLearner::SuperLearner(Y = sub_Z1_S1[[Y_name]],
                                                    X = sub_Z1_S1[, X_name, drop = FALSE],
                                                    family = family,
                                                    SL.library = Y_X_library, 
                                                    cvControl = list(V = v_folds))
    
    sub_Z1_S0 <- data[data[[Z_name]] == 1 & data[[S_name]] == 0,]
    out$fit_Y_Z1_S0_X <- SuperLearner::SuperLearner(Y = sub_Z1_S0[[Y_name]],
                                                    X = sub_Z1_S0[, X_name, drop = FALSE],
                                                    family = family,
                                                    SL.library = Y_X_library, 
                                                    cvControl = list(V = v_folds))
    
    sub_Z0_S1 <- data[data[[Z_name]] == 0 & data[[S_name]] == 1,]
    out$fit_Y_Z0_S1_X <- SuperLearner::SuperLearner(Y = sub_Z0_S1[[Y_name]],
                                                    X = sub_Z0_S1[, X_name, drop = FALSE],
                                                    family = family,
                                                    SL.library = Y_X_library, 
                                                    cvControl = list(V = v_folds))
    
  }
  
  # only needed for AIPW sensitivity analysis
  if(any(c("sens") %in% method)){
    out$fit_S_Z_X <- SuperLearner::SuperLearner(
      Y = data[[S_name]],
      X = data[, c(Z_name, X_name), drop = FALSE],
      family = "binomial",
      SL.library = S_Z_X_library,
      cvControl = list(V = v_folds))
  }
  
  # needed for all but gcomp
  if(any(c("aipw", "sens", "tmle") %in% method)){
    out$fit_Z_X <- SuperLearner::SuperLearner(Y = data[[Z_name]],
                                              X = data[, X_name, drop = FALSE],
                                              family = binomial(),
                                              SL.library = Z_X_library, 
                                              cvControl = list(V = v_folds))
  }
  
  return(out)
  
}
