# AIC/BIC select the model with the smaller criterion.
benchmark_aic_bic <- function(models) {
    c(AIC = names(models)[which.min(vapply(models, AIC, numeric(1)))],
      BIC = names(models)[which.min(vapply(models, BIC, numeric(1)))])
}
