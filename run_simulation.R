source("drm_selection.R")
source("data_generating.R")
source("benchmark_aic_bic.R")
source("main.R")
library(dplyr)

run_linear_simulation <- function(structure = "Additional main effects", truth = "M2", snr = 1,
    n_train = 200L, n_test = 200L, repetitions = 500L, basis = "adaptive",
    kernel_bw = "adaptive", dim_basis = "adaptive", predictor_distribution = "Normal",
    error_distribution = "Normal", test_seed = 4100L, training_seed = 42L, alpha = 0.05, maxit = 2000L,
    snr_by = c("error_variance", "coefficient", "predictor_variance"),
    coefficient = 1, predictor_variance = 1, error_variance = 1) {
    snr_by <- match.arg(snr_by)
    truth <- match.arg(truth, truth_levels)
    structure <- match.arg(structure, structure_levels)
    controls <- c(snr, coefficient, predictor_variance, error_variance)
    if (length(controls) != 4L || any(!is.finite(controls) | controls <= 0))
        stop("snr, coefficient, predictor_variance and error_variance must be positive finite numbers.")
    simulation_start <- proc.time()[["elapsed"]]
    complexity <- candidate_parameter_count(structure)
    set.seed(test_seed)
    x_test <- generate_predictors(n_test, predictor_distribution)
    set.seed(training_seed)
    results <- parameters <- vector("list", repetitions)

    for (i in seq_len(repetitions)) {
        x_train <- generate_predictors(n_train, predictor_distribution)
        data <- generate_experiment(x_train, x_test, structure, truth, snr, snr_by,
                                    coefficient, predictor_variance, error_variance, error_distribution)
        parameters[[i]] <- data$parameters
        models <- fit_candidate_models(data$train, structure)
        prediction_1 <- as.numeric(predict(models$M1, newdata = data$test))
        prediction_2 <- as.numeric(predict(models$M2, newdata = data$test))
        drm_start <- proc.time()[["elapsed"]]
        drm <- safe_drm_model_selection(data$test$y, prediction_1, prediction_2, basis,
            complexity = complexity, alpha = alpha, maxit = maxit, kernel_bw = kernel_bw, dim_basis = dim_basis)
        drm$runtime_seconds <- proc.time()[["elapsed"]] - drm_start
        drm$method <- "DRM"
        benchmarks <- benchmark_aic_bic(models)
        reference <- data.frame(method = names(benchmarks), decision = unname(benchmarks), valid = TRUE)
        results[[i]] <- mutate(bind_rows(drm, reference), replicate = i)
    }

    parameters <- bind_rows(parameters) |> mutate(replicate = row_number())
    raw <- bind_rows(results) |> left_join(parameters, by = "replicate") |>
        mutate(structure = structure, truth = truth, snr = snr, snr_by = snr_by,
               n_train = n_train, n_test = n_test, error_sd = sqrt(error_variance), correct = decision == truth)
    summary <- raw |> group_by(method) |>
        summarise(valid_runs = sum(valid), invalid_runs = sum(!valid),
                  correct_selections = sum(correct, na.rm = TRUE),
                  correct_rate = if (any(valid)) mean(correct[valid]) else NA_real_, .groups = "drop")
    failures <- raw |> filter(method == "DRM") |>
        mutate(outcome = case_when(!valid ~ "Invalid DRM fit", correct ~ "Correct",
                                  !rejected ~ "Failed to reject equality", rejected ~ "Wrong Wald direction")) |>
        count(outcome, name = "runs") |> mutate(proportion = runs / sum(runs))
    error_summary <- raw |> filter(method == "DRM", !valid) |>
        count(failure_stage, error_type, error_message, name = "runs")
    settings <- list(structure = structure, truth = truth, snr = snr, snr_by = snr_by,
        coefficient = coefficient, predictor_variance = predictor_variance, error_variance = error_variance,
        snr_definition = "Full signal sample variance / error sample variance, conditional on training design",
        test_response = "Noiseless true signal; same DGP parameters as the corresponding training sample",
        n_train = n_train, n_test = n_test, repetitions = repetitions, basis = basis,
        kernel_bw = kernel_bw, dim_basis = dim_basis, predictor_distribution = predictor_distribution,
        error_distribution = error_distribution, test_seed = test_seed, training_seed = training_seed,
        alpha = alpha, maxit = maxit)
    runtime <- list(total_seconds = proc.time()[["elapsed"]] - simulation_start,
                    drm_total_seconds = sum(raw$runtime_seconds[raw$method == "DRM"]),
                    drm_mean_seconds = mean(raw$runtime_seconds[raw$method == "DRM"]))
    list(settings = settings, parameters = parameters, raw = raw, summary = summary,
         drm_outcomes = failures, error_summary = error_summary, runtime = runtime)
}
