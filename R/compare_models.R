compare_models <- function(...) {

  models <- list(...)
  if (!length(models) >= 2) {
    stop("Please provide at least two models to compare.")
  }

  waic_exists <- sapply(models, function(model) {
    "WAIC" %in% names(model)
  })
  if (!all(waic_exists)) {
    stop("Not all models have WAIC calculated. Ensure that WAIC = TRUE if using run_dpcd.")
  }

  model_names <- sapply(substitute(list(...))[-1], deparse)
  waic_values <- sapply(models, function(x) x$WAIC$WAIC)
  waic_df <- data.frame(
    "Name" = model_names,
    "WAIC" = waic_values
  )

  waic_df[order(waic_df$WAIC), ]
}
