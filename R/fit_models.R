#' Function to fit all models needed by various estimators
#' 
#' @param data dataset used to fit models
#' @param G_name name of growth outcome variable, default G
#' @param V_name name of vaccine treatment variable, default V
#' @param X_name character vector containing name(s) of covariates, default X
#' @param Y_name name of infection variable, default Y
#' @param est character vector of names of estimators to use for growth effect
#' @param G_X_model optional specify model to be used for fitting growth on covariates, otherwise growth on all covariates
#' @param Y_X_model optional specify model to be used for fitting infection on covariates, otherwise infection on all covariates
#' @param family family for outcome model, defaults to gaussian for growth
#' 
#' @export
#' 
#' @returns list with models needed for specified estimators
fit_models <- function(data,
                       G_name = "G",
                       V_name = "V",
                       X_name = c("X"),
                       Y_name = "Y",
                       est = c("gcomp_pop_estimand", 
                               "gcomp",
                               "ipw",
                               "efficient_aipw",
                               "efficient_aipw_sens",
                               "efficient_tmle"),
                       G_V_X_model = NULL,
                       G_X_Y1_model = NULL,
                       G_X_Y0_model = NULL,
                       Y_X_model = NULL,
                       V_X_model = paste0(V_name, " ~ 1"),
                       family = "gaussian"){
  
  # Prep model formulas if not pre-specified
  if(is.null(G_V_X_model)){
    G_V_X_model <- as.formula(paste0(G_name, "~",
                                     V_name, "+",
                                     paste0(X_name, collapse = "+")))
  }
  
  if(is.null(G_X_Y1_model)){
    G_X_Y1_model <- as.formula(paste0(G_name, "~",
                                      paste0(X_name, collapse = "+")))
  }
  
  if(is.null(G_X_Y0_model)){
    G_X_Y0_model <- as.formula(paste0(G_name, "~",
                                      paste0(X_name, collapse = "+")))
  }
  
  if(is.null(Y_X_model)){
    Y_X_model <- as.formula(paste0(Y_name, "~",
                                   paste0(X_name, collapse = "+")))
  }
  
  out <- list()
  
  # only needed for population estimator
  if("gcomp_pop_estimand" %in% est){
    out$fit_G_V_X <- glm(G_V_X_model, data = data, family = family)
  } 
  
  # needed for any weight-based estimator
  if(any(c("gcomp", "efficient_aipw", "efficient_tmle", "efficient_aipw_sens") %in% est)){
    sub_V0 <- data[data[[V_name]] == 0,]
    out$fit_Y_V0_X <- glm(Y_X_model, sub_V0, family = "binomial")
    
    sub_V1 <- data[data[[V_name]] == 1,]
    out$fit_Y_V1_X <- glm(Y_X_model, sub_V1, family = "binomial")
    
    sub_V1_Y1 <- data[data[[V_name]] == 1 & data[[Y_name]] == 1,]
    out$fit_G_V1_Y1_X <- glm(G_X_Y1_model, data = sub_V1_Y1, family = family)
    
    sub_V1_Y0 <- data[data[[V_name]] == 1 & data[[Y_name]] == 0,]
    out$fit_G_V1_Y0_X <- glm(G_X_Y0_model, data = sub_V1_Y0, family = family)
    
    sub_V0_Y1 <- data[data[[V_name]] == 0 & data[[Y_name]] == 1,]
    out$fit_G_V0_Y1_X <- glm(G_X_Y1_model, data = sub_V0_Y1, family = family)
    
  }
  
  if(any(c("efficient_aipw_sens") %in% est)){
    out$fit_Y_V_X <- glm(
      paste0(Y_name, "~", V_name, "+", paste0(X_name, collapse = "+")),
      family = "binomial",
      data = data
    )
  }
  
  if(any(c("efficient_aipw_sens", "efficient_aipw", "ipw", "efficient_tmle") %in% est)){
    out$fit_V_X <- glm(
      V_X_model, family = "binomial", data = data
    )
  }
  return(out)
  
}


#' Function to fit all SuperLearner models needed by various estimators
#' 
#' @param data dataset used to fit models
#' @param G_name name of growth outcome variable, default G
#' @param V_name name of vaccine treatment variable, default V
#' @param X_name character vector containing name(s) of covariates, default X
#' @param Y_name name of infection variable, default Y
#' @param est character vector of names of estimators to use for growth effect
#' @param G_X_library specify SuperLearner libraries to be used for fitting growth on covariates. Default GLM 
#' @param Y_X_library specify SuperLearner libraries to be used for fitting infection on covariates. Default GLM 
#' @param family family for outcome model, defaults to gaussian for growth
#' @param v_folds number of cross-validation folds to use in SuperLearner, default 3
#' 
#' @export
#' 
#' @returns list with models needed for specified estimators
fit_ml_models <- function(data,
                          G_name = "G",
                          V_name = "V",
                          X_name = c("X"),
                          Y_name = "Y",
                          est = c("gcomp_pop_estimand", 
                                  "gcomp",
                                  "ipw",
                                  "efficient_aipw", 
                                  "efficient_tmle"),
                          G_V_X_library = c("SL.glm"),
                          G_X_library = c("SL.glm"),
                          Y_X_library = c("SL.glm"),
                          Z_X_library = c("SL.mean"),
                          family = "gaussian",
                          v_folds = 3){
  
  out <- list()
  
  # only needed for population estimator
  if("gcomp_pop_estimand" %in% est){
    out$fit_G_V_X <- SuperLearner::SuperLearner(Y = data[[G_name]],
                                                X = data[, colnames(data) %in% c(V_name, X_name), drop = FALSE],
                                                family = family,
                                                SL.library = G_V_X_library, 
                                                cvControl = list(V = v_folds))
    
  } 
  
  # needed for any weight-based estimator
  if(any(c("gcomp", "efficient_aipw", "efficient_tmle") %in% est)){
    sub_V0 <- data[data[[V_name]] == 0,]
    out$fit_Y_V0_X <- SuperLearner::SuperLearner(Y = sub_V0[[Y_name]],
                                                 X = sub_V0[, X_name, drop = FALSE],
                                                 family = stats::binomial(),
                                                 SL.library = Y_X_library, 
                                                 cvControl = list(V = v_folds))
    
    sub_V1 <- data[data[[V_name]] == 1,]
    out$fit_Y_V1_X <- SuperLearner::SuperLearner(Y = sub_V1[[Y_name]],
                                                 X = sub_V1[, X_name, drop = FALSE],
                                                 family = stats::binomial(),
                                                 SL.library = Y_X_library, 
                                                 cvControl = list(V = v_folds))
    
    sub_V1_Y1 <- data[data[[V_name]] == 1 & data[[Y_name]] == 1,]
    out$fit_G_V1_Y1_X <- SuperLearner::SuperLearner(Y = sub_V1_Y1[[G_name]],
                                                    X = sub_V1_Y1[, X_name, drop = FALSE],
                                                    family = family,
                                                    SL.library = G_X_library, 
                                                    cvControl = list(V = v_folds))
    
    sub_V1_Y0 <- data[data[[V_name]] == 1 & data[[Y_name]] == 0,]
    out$fit_G_V1_Y0_X <- SuperLearner::SuperLearner(Y = sub_V1_Y0[[G_name]],
                                                    X = sub_V1_Y0[, X_name, drop = FALSE],
                                                    family = family,
                                                    SL.library = G_X_library, 
                                                    cvControl = list(V = v_folds))
    
    sub_V0_Y1 <- data[data[[V_name]] == 0 & data[[Y_name]] == 1,]
    out$fit_G_V0_Y1_X <- SuperLearner::SuperLearner(Y = sub_V0_Y1[[G_name]],
                                                    X = sub_V0_Y1[, X_name, drop = FALSE],
                                                    family = family,
                                                    SL.library = G_X_library, 
                                                    cvControl = list(V = v_folds))
    
  }
  
  if(any(c("ipw", "efficient_aipw", "efficient_tmle") %in% est)){
    out$fit_V_X <- SuperLearner::SuperLearner(Y = data[[V_name]],
                                              X = data[, X_name, drop = FALSE],
                                              family = binomial(),
                                              SL.library = V_X_library, 
                                              cvControl = list(V = v_folds))
  }
  
  return(out)
  
}
