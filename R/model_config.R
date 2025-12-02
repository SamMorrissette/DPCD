#' Configure Nimble models
#'
#' Helper functions to configure nimble models with model-specific hyperparameters and starting values.
#'
#' @noRd

get_config <- function(model_name, d_obs, p, trunc_value, user_hyper_params, user_inits) {
  switch(
    model_name,
    "UU" = configure_uu(d_obs, p, trunc_value, user_hyper_params, user_inits),
    "EU" = configure_eu(d_obs, p, trunc_value, user_hyper_params, user_inits),
    "UD" = configure_ud(d_obs, p, trunc_value, user_hyper_params, user_inits),
    "ED" = configure_ed(d_obs, p, trunc_value, user_hyper_params, user_inits),
    "US" = configure_us(d_obs, p, trunc_value, user_hyper_params, user_inits),
    "ES" = configure_es(d_obs, p, trunc_value, user_hyper_params, user_inits)
  )
}

configure_model_common <- function(d_obs,
                                   p,
                                   trunc_value,
                                   user_hyper_params,
                                   user_inits,
                                   model_hyper_params,
                                   model_inits,
                                   model_inits_override) {
  # Hyperparameters
  common_hyper_params <- list(alpha_0 = 1,
                       a_0 = 1,
                       b_0 = 1,
                       lambda = 1,
                       mu_0 = rep(0, p))

  default_hyper_params <- modifyList(common_hyper_params, model_hyper_params)

  hyper_params <- NULL
  if (is.null(user_hyper_params) || length(user_hyper_params) == 0) {
    hyper_params <- default_hyper_params
  } else {
    if (!is.list(user_hyper_params)) {
      stop("`hyper_params` must be a list")
    }

    if (!all(names(user_hyper_params) %in% names(default_hyper_params))) {
      stop("`hyper_params` must be a named list containing values for one or more hyperparameters for the selected model.")
    }
    hyper_params <- modifyList(default_hyper_params, user_hyper_params)
  }

  # Constants
  n <- nrow(d_obs)
  common_constants <- list(p = p,
                           n = n,
                           N = trunc_value)
  constants <- c(common_constants, hyper_params)

  # Initial values
  mds_init <- cmdscale(d_obs, k = p)
  common_inits <- list(
    x = mds_init,
    delta = as.matrix(dist(mds_init)),
    z = sample(1:trunc_value, n, replace = TRUE),
    beta = rep(0.5, trunc_value - 1),
    mu = matrix(0, nrow = trunc_value, ncol = p),
    sigma_sq = 1
  )

  default_inits <- modifyList(common_inits, model_inits)

  if (is.null(user_inits) || length(user_inits) == 0) {
    init_params <- default_inits
  } else {
    if (!is.list(user_inits)) {
      stop("`init_params` must be a list")
    }

    if (!all(names(user_inits) %in% names(default_inits))) {
      stop("`init_params` must be a named list containing initial values for one or more model parameters for the selected model.")
    }
    init_params <- modifyList(default_inits, user_inits)
  }

  list(constants = constants,
       init_params = init_params)
}


configure_uu <- function(d_obs, p, trunc_value, user_hyper_params, user_inits) {
  model_hyper_params <- list("Psi_0" = diag(p),
                             "nu_0" = p + 2)

  model_inits <-  list(Sigma = array(diag(p), dim = c(p,p,trunc_value)))

  configure_model_common(d_obs, p, trunc_value,
                       user_hyper_params, user_inits,
                       model_hyper_params, model_inits)
}

configure_eu <- function(d_obs, p, trunc_value, user_hyper_params, user_inits) {
  model_hyper_params <- list("Psi_0" = diag(p),
                             "nu_0" = p + 2)

  model_inits <-  list(Sigma = array(diag(p), dim = c(p,p)))

  configure_model_common(d_obs, p, trunc_value,
                       user_hyper_params, user_inits,
                       model_hyper_params, model_inits)
}

configure_ud <- function(d_obs, p, trunc_value, user_hyper_params, user_inits) {
  model_hyper_params <- list("alpha_d" = 1,
                             "beta_d" = 1)

  model_inits <-  list(tau_vec = matrix(1, nrow = p, ncol = trunc_value))

  configure_model_common(d_obs, p, trunc_value,
                       user_hyper_params, user_inits,
                       model_hyper_params, model_inits)
}

configure_ed <- function(d_obs, p, trunc_value, user_hyper_params, user_inits) {
  model_hyper_params <- list("alpha_d" = 1,
                             "beta_d" = 1)

  model_inits <-  list(tau_vec = rep(1, p))

  configure_model_common(d_obs, p, trunc_value,
                       user_hyper_params, user_inits,
                       model_hyper_params, model_inits)
}

configure_us <- function(d_obs, p, trunc_value, user_hyper_params, user_inits) {
  model_hyper_params <- list("alpha_tau" = 1,
                             "beta_tau" = 1)

  model_inits <-  list(tau_sq = rep(1, trunc_value))

  configure_model_common(d_obs, p, trunc_value,
                       user_hyper_params, user_inits,
                       model_hyper_params, model_inits)
}
configure_es <- function(d_obs, p, trunc_value, user_hyper_params, user_inits) {
  model_hyper_params <- list("alpha_tau" = 1,
                             "beta_tau" = 1)

  model_inits <-  list(tau_sq = 1)

  configure_model_common(d_obs, p, trunc_value,
                       user_hyper_params, user_inits,
                       model_hyper_params, model_inits)
}
