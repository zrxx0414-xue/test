# Linear regression simulation and model selection.

library(drmdel)

default_basis_functions <- function() {
    list(`q(y) = (y)` = list(d = 1L, fun = function(y) c(y = y)), `q(y) = (y^2)` = list(d = 1L, fun = function(y) c(y2 = y^2)), 
        `q(y) = (y, y^2)` = list(d = 2L, fun = function(y) c(y = y, y2 = y^2)), `q(y) = (y, y^2, y^3)` = list(d = 3L, 
            fun = function(y) c(y = y, y2 = y^2, y3 = y^3)))
}

quadratic_form <- function(vector, covariance, scale) {
    scale * drop(t(vector) %*% solve(covariance, vector))
}

diagnostic_drmdel <- function(x, n_samples, basis_func, control, diagnostics) {
    fitter <- drmdel::drmdel
    local_environment <- new.env(parent = environment(fitter))
    local_environment$optim <- function(...) {
        diagnostics$stage <- "drm_optimization"
        opt <- stats::optim(...)
        diagnostics$convergence <- opt$convergence
        diagnostics$optimizer_message <- if (is.null(opt$message)) NA_character_ else opt$message
        if (opt$convergence != 0L) {
            stop("DRM optimizer did not converge; optim code = ", opt$convergence, call. = FALSE)
        }
        diagnostics$stage <- "drm_postprocessing"
        opt
    }
    environment(fitter) <- local_environment
    diagnostics$stage <- "drm_setup"
    fitter(x = x, n_samples = n_samples, basis_func = basis_func, control = control)
}

drm_model_selection <- function(baseline, prediction_1, prediction_2, basis_name, complexity = c(M1 = 1L, M2 = 2L), 
    alpha = 0.05, maxit = 2000L, basis_functions = default_basis_functions(), kernel_bw = "adaptive", dim_basis = "adaptive", 
    .diagnostics = new.env(parent = emptyenv())) {
    .diagnostics$stage <- "basis_construction"
    .diagnostics$basis <- basis_name
    .diagnostics$basis_ready <- FALSE
    .diagnostics$convergence <- NA_integer_
    .diagnostics$optimizer_message <- NA_character_
    .diagnostics$dimension <- NA_integer_
    basis <- if (identical(basis_name, "adaptive")) {
        fit_adaptive_basis(data = list(baseline, prediction_1, prediction_2), kernel_bw = kernel_bw, dim_basis = dim_basis)
    } else {
        basis_functions[[basis_name]]
    }
    d <- basis$d
    .diagnostics$dimension <- d
    .diagnostics$basis_ready <- TRUE
    pooled_values <- c(baseline, prediction_1, prediction_2)
    total_n <- length(pooled_values)
    sample_sizes <- c(length(baseline), length(prediction_1), length(prediction_2))
    fit <- diagnostic_drmdel(x = pooled_values, n_samples = sample_sizes, basis_func = basis$fun, control = list(maxit = maxit), 
        diagnostics = .diagnostics)
    estimate <- fit$mele
    .diagnostics$stage <- "covariance"
    covariance <- as.matrix(meleCov(fit))
    .diagnostics$stage <- "wald_contrast"
    beta_1_indices <- 2:(d + 1L)
    beta_2_indices <- (d + 3L):(2L * d + 2L)
    contrast <- matrix(0, nrow = d, ncol = length(estimate))
    for (component in seq_len(d)) {
        contrast[component, beta_1_indices[component]] <- 1
        contrast[component, beta_2_indices[component]] <- -1
    }
    beta_difference <- drop(contrast %*% estimate)
    difference_covariance <- contrast %*% covariance %*% t(contrast)
    .diagnostics$stage <- "wald_equality"
    equality_wald <- quadratic_form(beta_difference, difference_covariance, total_n)
    equality_p_value <- pchisq(equality_wald, df = d, lower.tail = FALSE)
    beta_1 <- estimate[beta_1_indices]
    beta_2 <- estimate[beta_2_indices]
    .diagnostics$stage <- "wald_M1"
    wald_1 <- quadratic_form(beta_1, covariance[beta_1_indices, beta_1_indices, drop = FALSE], total_n)
    .diagnostics$stage <- "wald_M2"
    wald_2 <- quadratic_form(beta_2, covariance[beta_2_indices, beta_2_indices, drop = FALSE], total_n)
    .diagnostics$stage <- "model_selection"
    rejected <- equality_p_value < alpha
    if (!rejected) {
        minimum_complexity <- min(complexity[c("M1", "M2")])
        candidates <- names(complexity)[complexity == minimum_complexity]
        decision <- if (length(candidates) == 1L) candidates else "Tie"
        reason <- "Failed to reject equality; complexity tie-break"
    } else if (wald_1 < wald_2) {
        decision <- "M1"
        reason <- "Rejected equality; M1 has smaller individual Wald"
    } else if (wald_2 < wald_1) {
        decision <- "M2"
        reason <- "Rejected equality; M2 has smaller individual Wald"
    } else {
        decision <- "Tie"
        reason <- "Rejected equality; individual Wald statistics are equal"
    }
    data.frame(decision = decision, reason = reason, rejected = rejected, equality_wald = equality_wald, equality_p = equality_p_value, 
        wald_1 = wald_1, wald_2 = wald_2, basis = basis_name, basis_dimension = d, valid = TRUE, failure_stage = NA_character_, 
        basis_ready = TRUE, drm_convergence = .diagnostics$convergence, optimizer_message = .diagnostics$optimizer_message, 
        error_message = NA_character_, stringsAsFactors = FALSE)
}

safe_drm_model_selection <- function(...) {
    diagnostics <- new.env(parent = emptyenv())
    diagnostics$stage <- "setup"
    diagnostics$basis <- NA_character_
    diagnostics$dimension <- NA_integer_
    diagnostics$basis_ready <- FALSE
    diagnostics$convergence <- NA_integer_
    diagnostics$optimizer_message <- NA_character_
    error_call <- NA_character_
    warnings <- character()
    result <- tryCatch(withCallingHandlers(drm_model_selection(..., .diagnostics = diagnostics), error = function(e) {
        calls <- vapply(sys.calls(), function(x) paste(deparse(x[[1]]), collapse = ""), character(1))
        if ("adaptive_bw" %in% calls) 
            diagnostics$stage <- "adaptive_bandwidth"
        else if ("adaptive_dim" %in% calls) 
            diagnostics$stage <- "adaptive_dimension"
        else if ("fit_adaptive_basis" %in% calls && "eigen" %in% calls) 
            diagnostics$stage <- "basis_eigensystem"
        else if ("fit_adaptive_basis" %in% calls) 
            diagnostics$stage <- "basis_construction"
        error_call <<- if (is.null(conditionCall(e))) 
            diagnostics$stage
        else substr(paste(deparse(conditionCall(e)), collapse = " "), 1L, 400L)
    }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = function(error) data.frame(decision = NA_character_, reason = "Invalid DRM fit", rejected = NA, 
        equality_wald = NA_real_, equality_p = NA_real_, wald_1 = NA_real_, wald_2 = NA_real_, basis = diagnostics$basis, 
        basis_dimension = diagnostics$dimension, valid = FALSE, failure_stage = diagnostics$stage, basis_ready = diagnostics$basis_ready, 
        drm_convergence = diagnostics$convergence, optimizer_message = diagnostics$optimizer_message, error_message = conditionMessage(error), 
        stringsAsFactors = FALSE))
    msg <- result$error_message
    result$error_type <- if (result$valid) NA_character_ else if (!is.na(result$drm_convergence) && result$drm_convergence != 0) "nonconvergence" else "other"
    if (!result$valid && grepl("infinite|missing values|non-finite|NA/NaN/Inf", msg)) result$error_type <- "nonfinite_values"
    if (!result$valid && grepl("singular|reciprocal condition", msg)) result$error_type <- "singular_matrix"
    if (!result$valid && grepl("dimension of the basis", msg)) result$error_type <- "invalid_dimension"
    result$error_call <- error_call
    result$warning_message <- if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_
    result
}

drm_selection_rule <- drm_model_selection

safe_drm_selection_rule <- safe_drm_model_selection

