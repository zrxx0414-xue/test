# Linear regression simulation and model selection.

structure_levels <- c("Additional main effects", "Additional quadratic term", "Additional interaction term")

truth_levels <- c("M1", "M2")

standardised_random <- function(n, distribution = "Normal") {
    values <- switch(distribution,
        "Normal" = rnorm(n), "Uniform" = runif(n, -sqrt(3), sqrt(3)),
        "Heavy-tailed t(3)" = rt(n, df = 3) / sqrt(3),
        "Bimodal" = (sample(c(-1, 1), n, replace = TRUE) + rnorm(n, sd = 0.5)) / sqrt(1.25),
        "Right-skewed" = rexp(n) - 1, stop("Unknown distribution: ", distribution))
    as.numeric(scale(values))
}

generate_predictors <- function(n, distribution = "Normal") {
    data.frame(x1 = standardised_random(n, distribution), x2 = standardised_random(n, distribution),
               x3 = standardised_random(n, distribution))
}

shared_signal <- function(data, structure) {
    switch(structure,
        "Additional main effects" = 1 + data$x1,
        "Additional quadratic term" = 1 + data$x1,
        "Additional interaction term" = 1 + data$x1 + data$x2)
}

additional_component <- function(data, structure) {
    switch(structure,
        "Additional main effects" = (data$x2 + data$x3) / sqrt(2),
        "Additional quadratic term" = (data$x1^2 - 1) / sqrt(2),
        "Additional interaction term" = data$x1 * data$x2)
}

# coefficient scales the complete non-intercept signal; M1 omits the additional term.
generate_signal <- function(data, structure, truth, coefficient = 1) {
    extra <- if (truth == "M2") additional_component(data, structure) else 0
    1 + coefficient * (shared_signal(data, structure) - 1 + extra)
}

# Calibrate on the training design; apply the SAME parameters to the test design.
generate_experiment <- function(x_train, x_test, structure, truth, snr, snr_by,
                                coefficient = 1, predictor_variance = 1, error_variance = 1,
                                error_distribution = "Normal") {
    signal_at <- function(v) generate_signal(x_train * sqrt(v), structure, truth, coefficient)
    if (snr_by == "error_variance") {
        error_variance <- var(signal_at(predictor_variance)) / snr
    } else if (snr_by == "coefficient") {
        unit_signal <- generate_signal(x_train * sqrt(predictor_variance), structure, truth, 1)
        coefficient <- sqrt(snr * error_variance / var(unit_signal))
    } else {
        # Nonlinear terms mean signal variance is not simply proportional to Var(X).
        target <- snr * error_variance
        upper <- 1
        while (var(signal_at(upper)) < target) upper <- 2 * upper
        predictor_variance <- uniroot(function(v) var(signal_at(v)) - target,
                                     interval = c(0, upper), tol = 1e-10)$root
    }
    train <- x_train * sqrt(predictor_variance)
    test <- x_test * sqrt(predictor_variance)
    mu <- generate_signal(train, structure, truth, coefficient)
    train$y <- mu + sqrt(error_variance) * standardised_random(nrow(train), error_distribution)
    test$y <- generate_signal(test, structure, truth, coefficient)
    parameters <- data.frame(coefficient = coefficient, predictor_variance = predictor_variance,
                             error_variance = error_variance, train_signal_variance = var(mu),
                             achieved_snr = var(mu) / error_variance)
    list(train = train, test = test, parameters = parameters)
}

fit_candidate_models <- function(data, structure) {
    switch(structure,
        "Additional main effects" = list(M1 = lm(y ~ x1, data), M2 = lm(y ~ x1 + x2 + x3, data)),
        "Additional quadratic term" = list(M1 = lm(y ~ x1, data), M2 = lm(y ~ x1 + I(x1^2), data)),
        "Additional interaction term" = list(M1 = lm(y ~ x1 + x2, data), M2 = lm(y ~ x1 * x2, data)))
}

candidate_parameter_count <- function(structure) {
    switch(structure,
        "Additional main effects" = c(M1 = 2L, M2 = 4L),
        "Additional quadratic term" = c(M1 = 2L, M2 = 3L),
        "Additional interaction term" = c(M1 = 3L, M2 = 4L))
}
