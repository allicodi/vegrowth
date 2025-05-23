% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap.R
\name{bootstrap_estimates}
\alias{bootstrap_estimates}
\title{Function to replicate n_boot bootstrap samples and get bootstrap standard error}
\usage{
bootstrap_estimates(
  data,
  G_name = "G",
  V_name = "V",
  X_name = "X",
  Y_name = "Y",
  n_boot = 1000,
  est = c("gcomp_pop_estimand", "gcomp", "efficient_aipw", "efficient_tmle",
    "hudgens_adj_upper", "hudgens_adj_lower"),
  ml = ml,
  G_V_X_model = NULL,
  G_X_Y1_model = NULL,
  G_X_Y0_model = NULL,
  Y_X_model = NULL,
  G_V_X_library = c("SL.glm"),
  G_X_library = c("SL.glm"),
  Y_X_library = c("SL.glm"),
  family = "gaussian",
  v_folds = 3,
  effect_dir = "positive"
)
}
\arguments{
\item{data}{original data to bootstrap}

\item{G_name}{growth outcome variable name}

\item{V_name}{vaccination variable name}

\item{X_name}{covariate name(s)}

\item{Y_name}{infection variable name}

\item{n_boot}{number of bootstrap replicates}

\item{est}{character vector of names of estimators to use for growth effect}

\item{ml}{boolean to use SuperLearner models, default FALSE}

\item{G_V_X_model}{optional specify model to be used for fitting growth on vaccine + covariates, otherwise growth on all covariates}

\item{G_X_Y1_model}{optional specify model to be used for fitting growth on covariates in the infected, otherwise growth on all covariates}

\item{G_X_Y0_model}{optional specify model to be used for fitting growth on covariates in the uninfected, otherwise growth on all covariates}

\item{Y_X_model}{optional specify model to be used for fitting infection on covariates, otherwise infection on all covariates}

\item{G_V_X_library}{optional specify SuperLearner libraries for model fitting growth on covariates + vaccine, default glm}

\item{G_X_library}{optional specify SuperLearner libraries for model fitting growth on covariates, default glm}

\item{Y_X_library}{optional specify SuperLearner libraries for model fitting infection on covariates, default glm}

\item{family}{family for outcome model, defaults to gaussian for growth}

\item{v_folds}{number of cross validation folds for SuperLearner, default 3}

\item{effect_dir}{direction of beneficial effect, defaults to "positive" for beneficial outcome. Used for one-side tests of bounds.}
}
\value{
list containing bootstrap se and 95\% CI bounds for estimators specified in est
}
\description{
Function to replicate n_boot bootstrap samples and get bootstrap standard error
}
