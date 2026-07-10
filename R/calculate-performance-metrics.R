# ─────────────────────────────────────────────────────────────────────────────
# Sheet-free scoring
# ─────────────────────────────────────────────────────────────────────────────

#' Predict dementia risk from exported Fine-Gray model artefacts.
#'
#' Sheet-free scoring function. Takes raw data and the exported model parameters
#' (coefficients, baseline hazard, transformation specs with knots/centres) and
#' produces predicted cumulative incidence at the specified horizon.
#'
#' @param data          Data frame with raw predictor columns (pre-transformation).
#' @param coefficients  Named numeric vector of log-SHR betas.
#' @param baseline_hazard Data frame with columns \code{time} and \code{hazard}
#'   (cumulative baseline subdistribution hazard, centred on training means).
#' @param specs         Transformation specs from extract_transformation_specs()
#'   or manually constructed. See apply_transformations() for format.
#' @param transformation_info Named list with knots and center_values for each
#'   transformed variable (as exported in transformation_info.csv).
#' @param horizon       Prediction horizon in years (default 5).
#' @return Numeric vector of predicted cumulative incidence at horizon.
#' @export
predict_fg_risk <- function(data, coefficients, baseline_hazard, specs,
                            transformation_info, horizon = 5) {
  result <- apply_transformations(data, specs, transformation_info)
  model_data <- result$data

  predictor_names <- names(coefficients)
  missing_preds <- setdiff(predictor_names, colnames(model_data))
  if (length(missing_preds) > 0) {
    stop(paste("Missing predictor columns:", paste(missing_preds, collapse = ", ")))
  }

  x <- as.matrix(model_data[, predictor_names, drop = FALSE])
  lp <- as.numeric(x %*% coefficients)

  idx_h <- findInterval(horizon, baseline_hazard$time)
  H0_h <- baseline_hazard$hazard[idx_h]

  # CIF = 1 - exp(-H0(t))^exp(lp)
  # baseline_hazard was computed with centered = TRUE, so lp must be centred
  # on the same training means (which it is, because the transformations use
  # the same centre values from transformation_info).
  as.numeric(1 - exp(-H0_h)^exp(lp))
}

# ─────────────────────────────────────────────────────────────────────────────
# Individual metric functions
# Each takes pre-computed inputs so they can be called independently or shared
# across apparent and bootstrap validation without redundant computation.
# ─────────────────────────────────────────────────────────────────────────────

#' Predicted CIF at horizon from a Fine-Gray model.
#' Baseline subdistribution survival S0(t) is computed once; subject-specific
#' CIF = 1 - S0(t)^exp(lp_i).
compute_fg_pred_risk_old <- function(model, data, horizon) {
  sf0   <- survival::survfit(model)
  idx_h <- findInterval(horizon, sf0$time)
  lp    <- predict(model, newdata = data, type = "lp")
  as.numeric(1 - sf0$surv[idx_h]^exp(lp))
}

compute_fg_pred_risk <- function(model, data, horizon) {
  # JL
  sf0   <- survival::survfit(model, data)
  as.numeric(1 - summary(sf0, times = horizon)$surv)
}

#' IPCW-weighted C-statistic for a Fine-Gray model.
#' Recreates the Fine-Gray expansion of `data` so that concordance is computed
#' on the same weighted pseudo-observation scale used during model fitting.
#' This matches `model$concordance` for the apparent case and correctly handles
#' competing events without treating them as simple censoring.
#' @return list(c_stat, c_lo, c_hi)
compute_c_statistic_fg <- function(model, data, time_to_event = "time_to_event", event_type = "event_type") {
  # Some modification by JL
  names(data)[which(names(data) == time_to_event)] <- "time_to_event"
  names(data)[which(names(data) == event_type)] <- "event_type"
  fg_data <- survival::finegray(
    survival::Surv(data$time_to_event, factor(data$event_type)) ~ .,
    data = data, etype = "1"
  )
  lp   <- predict(model, newdata = fg_data, type = "lp")
  conc <- survival::concordance(
    survival::Surv(fg_data$fgstart, fg_data$fgstop, fg_data$fgstatus) ~ lp,
    weights = fg_data$fgwt,
    reverse = TRUE
  )
  c_stat <- as.numeric(conc$concordance)
  c_se   <- as.numeric(sqrt(conc$var))
  list(c_stat = c_stat, c_lo = c_stat - 1.96 * c_se, c_hi = c_stat + 1.96 * c_se)
}

#' Ratio of predicted risks at 95th vs 5th percentile.
compute_risk_ratio_fg <- function(pred_risk) {
  as.numeric(quantile(pred_risk, 0.95, na.rm = TRUE) / quantile(pred_risk, 0.05, na.rm = TRUE))
}

#' Observed CIF at horizon via Aalen-Johansen estimator.
#' @return list(obs_cif, obs_ci_lo, obs_ci_hi)
compute_obs_cif_fg <- function(data, horizon, time_to_event = "time_to_event", event_type = "event_type") {
  # Some modification by JL
  names(data)[which(names(data) == time_to_event)] <- "time_to_event"
  names(data)[which(names(data) == event_type)] <- "event_type"
  if (!any(data$event_type == 1, na.rm = TRUE)) {
    return(list(obs_cif = 0, obs_ci_lo = 0, obs_ci_hi = 0))
  }
  ci      <- cmprsk::cuminc(ftime = data$time_to_event, fstatus = data$event_type, cencode = 0)
  tp      <- cmprsk::timepoints(ci, horizon)
  if (!"1 1" %in% rownames(tp$est)) {
    return(list(obs_cif = 0, obs_ci_lo = 0, obs_ci_hi = 0))
  }
  obs_cif <- tp$est["1 1", ]
  obs_var <- tp$var["1 1", ]
  list(
    obs_cif   = as.numeric(obs_cif),
    obs_ci_lo = as.numeric(obs_cif - 1.96 * sqrt(obs_var)),
    obs_ci_hi = as.numeric(obs_cif + 1.96 * sqrt(obs_var))
  )
}

#' Traditional calibration slope from a secondary Fine-Gray model.
#' Fits a Fine-Gray model with the linear predictor (LP) from the primary model
#' as the sole covariate. Perfect calibration yields a coefficient of 1 on the LP.
#' The calibration intercept is not estimated here; the baseline subdistribution
#' hazard (saved in the model export) serves that role.
#' @param pred_risk numeric vector of predicted CIFs (unused, kept for API compatibility).
#' @param data      data frame with time_to_event and event_type columns.
#' @param horizon   prediction horizon in years (unused, kept for API compatibility).
#' @param model_for_slope the primary coxph model, used to extract the LP.
#' @return list(slope)
compute_calibration_fg <- function(pred_risk, data, horizon, model_for_slope = NULL, time_to_event = "time_to_event", event_type = "event_type") {
  # Some modification by JL
  names(data)[which(names(data) == time_to_event)] <- "time_to_event"
  names(data)[which(names(data) == event_type)] <- "event_type"

  keep      <- !is.na(pred_risk)
  data      <- data[keep, ]

  lp <- predict(model_for_slope, newdata = data, type = "lp")
  data$lp_primary <- lp

  fg_data <- survival::finegray(
    survival::Surv(time_to_event, factor(event_type)) ~ lp_primary,
    data = data, etype = "1"
  )
  slope_model <- survival::coxph(
    survival::Surv(fgstart, fgstop, fgstatus) ~ lp_primary,
    data = fg_data, weights = fg_data$fgwt
  )

  list(slope = as.numeric(coef(slope_model)[["lp_primary"]]))
}

#' Annotation string of calibration + discrimination metrics for a calibration plot.
#'
#' Reports the OPTIMISM-CORRECTED calibration slope, Integrated Calibration Index
#' (ICI), and Harrell's c-statistic, each with 95% CIs, from a
#' bootstrap_fg_validation() object, formatted as a two-block label for placement
#' on the calibration plot. Optimism-corrected values are used because the
#' apparent calibration slope is ~1 by construction on the training data and
#' therefore uninformative. A calibration intercept is omitted: in the Fine-Gray
#' framework it is not separately estimated (the baseline subdistribution hazard
#' plays that role; see compute_calibration_fg()).
#' @param bv bootstrap_fg_validation() output (uses its $optimism_corrected table).
#' @return A single character string with embedded newlines.
calibration_plot_annotation <- function(bv) {
  oc <- stats::setNames(bv$optimism_corrected$optimism_corrected,
                        bv$optimism_corrected$metric)
  f <- function(x, d = 4) formatC(as.numeric(x), digits = d, format = "f")
  paste0(
    "Calibration (optimism-corrected)\n",
    "  slope: ", f(oc["calibration_slope"]),
    " (", f(oc["calibration_slope_ci_lo"]), "–", f(oc["calibration_slope_ci_hi"]), ")\n",
    "  ICI: ", f(oc["ici"]),
    " (", f(oc["ici_ci_lo"]), "–", f(oc["ici_ci_hi"]), ")\n",
    "Discrimination (optimism-corrected)\n",
    "  c-statistic: ", f(oc["c_statistic"]),
    " (", f(oc["c_statistic_ci_lo"]), "–", f(oc["c_statistic_ci_hi"]), ")"
  )
}

#' IPCW Brier score properly accounting for competing risks.
#' The censoring distribution G(t) is estimated treating both dementia and
#' death as observed events, so that only true censoring (event_type == 0)
#' contributes to G(t). This is the correct IPCW formulation for competing risks.
#' @param pred_risk numeric vector of predicted CIFs at horizon.
#' @param data data frame with time_to_event and event_type columns.
#' @param horizon prediction horizon in years.
#' @return list(brier, brier_scaled)
compute_brier_fg <- function(pred_risk, data, horizon, time_to_event = "time_to_event") {
  # Some modification by JL
  names(data)[which(names(data) == time_to_event)] <- "time_to_event"

  followup <- pmin(data$time_to_event, horizon)
  # Observed status at horizon: 1 if dementia before horizon, 0 otherwise
  event_horizon <- as.integer(data$time_to_event <= horizon & data$event_type == 1)

  # Censoring indicator: 1 = censored (event_type == 0 only).
  # Competing events (event_type == 2) are NOT treated as censoring for IPCW.
  censor_ind <- as.integer(data$event_type == 0 & data$time_to_event <= horizon)

  # KM estimate of censoring distribution G(t), treating both dementia and
  # death as "observed" (not censored) for the purposes of IPCW.
  observed_ind <- as.integer(data$event_type != 0)
  censor_km_fit <- survival::survfit(
    survival::Surv(followup, 1L - observed_ind) ~ 1
  )
  g_fun <- function(t) {
    idx <- max(c(0L, which(censor_km_fit$time < t)))
    if (idx == 0L) 1.0 else censor_km_fit$surv[idx]
  }
  g_at_followup <- vapply(followup, g_fun, numeric(1))
  g_at_horizon  <- g_fun(horizon)

  weight_i <- ifelse(
    data$time_to_event <= horizon & data$event_type != 0,
    1.0 / pmax(g_at_followup, 1e-8),
    ifelse(
      data$time_to_event > horizon,
      1.0 / pmax(g_at_horizon, 1e-8),
      0.0
    )
  )

  # Observed CIF at horizon for null model
  obs_cif <- compute_obs_cif_fg(data, horizon)$obs_cif

  brier      <- mean(weight_i * (pred_risk - event_horizon)^2)
  brier_null <- mean(weight_i * (obs_cif - event_horizon)^2)
  list(brier = brier, brier_scaled = 1 - brier / brier_null)
}

#' Integrated Calibration Index (ICI) for a Fine-Gray model.
#' ICI is the weighted mean absolute difference between observed (from the
#' Austin et al. flexible calibration curve) and predicted risks across the
#' distribution of predictions.
#' @param pred_risk numeric vector of predicted CIFs.
#' @param data      data frame with time_to_event and event_type columns.
#' @param horizon   prediction horizon in years.
#' @param n_knots   number of RCS knots for the flexible model (default 4).
#' @param n_grid    number of grid points for the curve (default 200).
#' @return numeric scalar: the ICI value.
compute_ici_fg <- function(pred_risk, data, horizon, n_knots = 4, n_grid = 200, time_to_event = "time_to_event") {
  # Some modification by JL
  names(data)[which(names(data) == time_to_event)] <- "time_to_event"

  keep      <- !is.na(pred_risk)
  pred_risk <- pred_risk[keep]
  data      <- data[keep, ]
  data$p_hat <- pred_risk

  fg_data <- survival::finegray(
    survival::Surv(time_to_event, factor(event_type)) ~ p_hat,
    data = data, etype = "1"
  )
  flex_model <- survival::coxph(
    survival::Surv(fgstart, fgstop, fgstatus) ~ rms::rcs(p_hat, n_knots),
    data = fg_data, weights = fg_data$fgwt
  )

  pred_grid  <- seq(min(pred_risk), max(pred_risk), length.out = n_grid)
  obs_smooth <- compute_fg_pred_risk(flex_model, data.frame(p_hat = pred_grid), horizon)

  mean(abs(obs_smooth - pred_grid), na.rm = TRUE)
}

#' Nagelkerke R² from Fine-Gray partial log-likelihoods.
#' Uses nrow(data) (original subjects) not model$n (expanded finegray rows).
compute_nagelkerke_r2_fg <- function(model, data) {
  ll_null <- model$loglik[[1]]
  ll_full <- model$loglik[[2]]
  n       <- nrow(data)
  r2_max  <- 1 - exp(2 * ll_null / n)
  (1 - exp(-2 * (ll_full - ll_null) / n)) / r2_max
}

# ─────────────────────────────────────────────────────────────────────────────
# Internal aggregator
# Collects all scalar metrics given a pre-computed pred_risk vector.
# Avoids recomputing pred_risk when sharing it across boot_app / boot_test.
# ─────────────────────────────────────────────────────────────────────────────

.compute_fg_metrics_vec <- function(model, data, pred_risk, horizon) {
  c_res     <- compute_c_statistic_fg(model, data)
  rr        <- compute_risk_ratio_fg(pred_risk)
  obs_res   <- compute_obs_cif_fg(data, horizon)
  cal_res   <- compute_calibration_fg(pred_risk, data, horizon, model_for_slope = model)
  brier_res <- compute_brier_fg(pred_risk, data, horizon)
  ici       <- compute_ici_fg(pred_risk, data, horizon)
  r2        <- compute_nagelkerke_r2_fg(model, data)

  c(
    c_statistic           = c_res$c_stat,
    c_stat_ci_lo          = c_res$c_lo,
    c_stat_ci_hi          = c_res$c_hi,
    risk_ratio_95_5       = rr,
    obs_cif               = obs_res$obs_cif,
    obs_cif_ci_lo         = obs_res$obs_ci_lo,
    obs_cif_ci_hi         = obs_res$obs_ci_hi,
    mean_predicted_risk   = mean(pred_risk),
    op_relative_diff_pct  = (obs_res$obs_cif - mean(pred_risk)) / obs_res$obs_cif * 100,
    calibration_slope     = cal_res$slope,
    ici                   = ici,
    brier_score           = brier_res$brier,
    brier_scaled          = brier_res$brier_scaled,
    nagelkerke_r2         = r2
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# Public wrappers
# ─────────────────────────────────────────────────────────────────────────────


#' Optimism-corrected performance via the bootstrap (Efron's method).
#'
#' For each of n_boot iterations:
#'   1. Draw bootstrap sample (with replacement).
#'   2. Refit Fine-Gray model on bootstrap sample.
#'   3. Compute baseline survfit ONCE for the bootstrap model.
#'   4. Derive pred_risk for bootstrap sample and original data by applying
#'      the bootstrap model's linear predictor to the shared baseline.
#'   5. Compute all metrics for both datasets (boot_apparent, boot_test).
#'   6. optimism_b = boot_apparent_b - boot_test_b.
#' Optimism-corrected = apparent - mean(optimism_b).
#' Bootstrap percentile 95% CIs are computed for all correctable metrics
#' (both apparent and optimism-corrected). Apparent CIs use the distribution
#' of boot_apparent values; optimism-corrected CIs use the distribution of
#' (apparent - optimism_b) across iterations.
#'
#' @param model   coxph object from fit_fine_and_gray_regression_model().
#' @param data    Original model dataset.
#' @param n_boot  Bootstrap iterations. Default 200.
#' @param horizon Prediction horizon in years. Default 5.
#' @param seed    Optional integer seed.
#'
#' @return list(apparent, optimism_corrected, n_boot_used)
#' @export
bootstrap_fg_validation <- function(model, data, n_boot = 200, horizon = 5,
                                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  correctable <- c(
    "c_statistic",
    "risk_ratio_95_5",
    "calibration_slope",
    "ici",
    "brier_score", "brier_scaled", "nagelkerke_r2"
  )

  predictors   <- labels(terms(model))
  apparent_vec <- .compute_fg_metrics_vec(
    model, data, compute_fg_pred_risk(model, data, horizon), horizon
  )

  optimism_mat <- matrix(
    NA_real_,
    nrow     = n_boot,
    ncol     = length(correctable),
    dimnames = list(NULL, correctable)
  )
  # Store boot_apparent values to compute apparent CIs
  boot_app_mat <- matrix(
    NA_real_,
    nrow     = n_boot,
    ncol     = length(correctable),
    dimnames = list(NULL, correctable)
  )
  n_failed <- 0L

  for (b in seq_len(n_boot)) {
    boot_idx  <- sample(nrow(data), replace = TRUE)
    boot_data <- data[boot_idx, ]

    result <- tryCatch({
      boot_model <- fit_fine_and_gray_regression_model(
        time_var   = "time_to_event",
        event_var  = "event_type",
        predictors = predictors,
        data       = boot_data
      )

      # Baseline computed once and shared for both pred_risk calculations
      sf0   <- survival::survfit(boot_model)
      idx_h <- findInterval(horizon, sf0$time)
      S0_h  <- sf0$surv[idx_h]

      pr_boot <- as.numeric(1 - S0_h^exp(predict(boot_model, newdata = boot_data, type = "lp")))
      pr_orig <- as.numeric(1 - S0_h^exp(predict(boot_model, newdata = data,      type = "lp")))

      boot_app_vec  <- .compute_fg_metrics_vec(boot_model, boot_data, pr_boot, horizon)
      boot_test_vec <- .compute_fg_metrics_vec(boot_model, data,      pr_orig, horizon)

      list(
        optimism     = boot_app_vec[correctable] - boot_test_vec[correctable],
        boot_apparent = boot_app_vec[correctable]
      )
    }, error = function(e) NULL)

    if (is.null(result)) {
      n_failed <- n_failed + 1L
    } else {
      optimism_mat[b, ]  <- result$optimism
      boot_app_mat[b, ]  <- result$boot_apparent
    }
  }

  if (n_failed > 0L) {
    warning(n_failed, " bootstrap iteration(s) failed and were skipped.")
  }

  mean_optimism <- colMeans(optimism_mat, na.rm = TRUE)
  corrected_vec <- apparent_vec[correctable] - mean_optimism

  # Bootstrap percentile CIs for apparent metrics (from boot_apparent distribution)
  app_ci_lo <- numeric(length(correctable))
  app_ci_hi <- numeric(length(correctable))
  names(app_ci_lo) <- names(app_ci_hi) <- correctable
  for (m in correctable) {
    app_ci_lo[m] <- as.numeric(quantile(boot_app_mat[, m], 0.025, na.rm = TRUE))
    app_ci_hi[m] <- as.numeric(quantile(boot_app_mat[, m], 0.975, na.rm = TRUE))
  }

  # Bootstrap percentile CIs for optimism-corrected metrics
  oc_ci_lo <- numeric(length(correctable))
  oc_ci_hi <- numeric(length(correctable))
  names(oc_ci_lo) <- names(oc_ci_hi) <- correctable
  for (m in correctable) {
    boot_corrected <- apparent_vec[m] - optimism_mat[, m]
    oc_ci_lo[m] <- as.numeric(quantile(boot_corrected, 0.025, na.rm = TRUE))
    oc_ci_hi[m] <- as.numeric(quantile(boot_corrected, 0.975, na.rm = TRUE))
  }

  # Build apparent data frame with CIs
  app_metrics <- character(0)
  app_values  <- numeric(0)
  for (m in names(apparent_vec)) {
    app_metrics <- c(app_metrics, m)
    app_values  <- c(app_values, apparent_vec[m])
    if (m %in% correctable) {
      app_metrics <- c(app_metrics, paste0(m, "_ci_lo"), paste0(m, "_ci_hi"))
      app_values  <- c(app_values, app_ci_lo[m], app_ci_hi[m])
    }
  }

  # Build optimism-corrected data frame with CIs
  oc_metrics <- character(0)
  oc_apparent <- numeric(0)
  oc_mean_opt <- numeric(0)
  oc_corrected <- numeric(0)
  for (m in correctable) {
    oc_metrics   <- c(oc_metrics,   m, paste0(m, "_ci_lo"), paste0(m, "_ci_hi"))
    oc_apparent  <- c(oc_apparent,  apparent_vec[m], NA_real_, NA_real_)
    oc_mean_opt  <- c(oc_mean_opt,  mean_optimism[m], NA_real_, NA_real_)
    oc_corrected <- c(oc_corrected, corrected_vec[m], oc_ci_lo[m], oc_ci_hi[m])
  }

  list(
    apparent = data.frame(
      metric = app_metrics,
      value  = as.numeric(app_values)
    ),
    optimism_corrected = data.frame(
      metric             = oc_metrics,
      apparent           = oc_apparent,
      mean_optimism      = oc_mean_opt,
      optimism_corrected = oc_corrected,
      row.names          = NULL
    ),
    n_boot_used = n_boot - n_failed
  )
}

#' Austin et al. smooth calibration curve for a Fine-Gray model.
#'
#' Fits a flexible Fine-Gray model with restricted cubic splines (RCS) of the
#' primary model's predicted CIF as the sole predictor, following Austin et al.
#' (2020). The flexible model's predicted CIF on a fine grid gives a smooth
#' "observed vs predicted" calibration curve that properly accounts for
#' competing risks — unlike naive 1-KM or decile groupings.
#'
#' @param model    coxph object from fit_fine_and_gray_regression_model().
#' @param data     Data frame with time_to_event, event_type, and predictors.
#' @param horizon  Prediction horizon in years (default 5).
#' @param n_knots  Number of RCS knots for the flexible model (default 4).
#' @param n_grid   Number of points on the calibration curve grid (default 200).
#'
#' @return list with:
#'   * curve: data.frame(predicted, observed) — the smooth calibration curve
#'   * p_hat: numeric vector of per-subject predicted CIFs (for rug plot)
#' @export
compute_fg_calibration_curve <- function(model, data, horizon = 5,
                                         n_knots = 4, n_grid = 200,
                                        time_to_event = "time_to_event") {
  # Some modification by JL
  names(data)[which(names(data) == time_to_event)] <- "time_to_event"

  # Step 1: predicted CIF from the primary model
  p_hat    <- compute_fg_pred_risk(model, data, horizon)
  keep     <- !is.na(p_hat)
  p_hat    <- p_hat[keep]
  data     <- data[keep, ]
  data$p_hat <- p_hat

  # Step 2: flexible Fine-Gray model — single predictor is RCS of p_hat
  fg_data <- survival::finegray(
    survival::Surv(time_to_event, factor(event_type)) ~ p_hat,
    data  = data,
    etype = "1"
  )
  flex_model <- survival::coxph(
    survival::Surv(fgstart, fgstop, fgstatus) ~ rms::rcs(p_hat, n_knots),
    data    = fg_data,
    weights = fg_data$fgwt
  )

  # Step 3: evaluate calibration curve on a grid
  pred_grid  <- seq(min(p_hat), max(p_hat), length.out = n_grid)
  obs_smooth <- compute_fg_pred_risk(flex_model, data.frame(p_hat = pred_grid), horizon)

  curve_data <- data.frame(predicted = pred_grid, observed = obs_smooth)
  curve_data <- curve_data[complete.cases(curve_data), ]

  list(
    curve = curve_data,
    p_hat = p_hat
  )
}

#' OvsP calibration by subgroup for a Fine-Gray model (sheet-free core).
#'
#' For each subgroup variable and each of its categories, computes the
#' Aalen-Johansen observed CIF and the mean model-predicted risk.
#'
#' @param model          A Fine-Gray model object (coxph from finegray data).
#' @param model_data     Data frame with predictor columns and subgroup variables.
#' @param sub_group_vars Character vector of categorical subgroup variable names.
#' @param category_labels Named character vector: names are "variable|recEnd",
#'   values are human-readable labels. If NULL, category codes are used as labels.
#' @param horizon        Prediction horizon in years (default 5).
#' @param min_events_pct Minimum observed CIF to include a subgroup (default 0.05).
#' @return A data frame with columns: variable, category, cat_label, n, observed, predicted, pct_diff.
#' @export
compute_fg_ovsp_subgroup_core <- function(model, model_data, sub_group_vars,
                                          category_labels = NULL, horizon = 5,
                                          min_events_pct = 0.05) {
  rows <- list()
  for (sg_var in sub_group_vars) {
    cats <- unique(as.character(model_data[[sg_var]]))
    cats <- cats[!is.na(cats)]
    for (cat_val in cats) {
      sub_data <- model_data[model_data[[sg_var]] == cat_val, ]
      if (nrow(sub_data) == 0) next
      obs <- compute_obs_cif_fg(sub_data, horizon)$obs_cif
      if (obs < min_events_pct) next
      pred <- mean(compute_fg_pred_risk(model, sub_data, horizon), na.rm = TRUE)
      lookup_key <- paste0(sg_var, "|", cat_val)
      cat_label <- if (!is.null(category_labels) && lookup_key %in% names(category_labels)) {
        category_labels[[lookup_key]]
      } else {
        cat_val
      }
      rows[[length(rows) + 1]] <- data.frame(
        variable         = sg_var,
        category         = as.character(cat_val),
        cat_label        = as.character(cat_label),
        n                = nrow(sub_data),
        observed         = obs,
        predicted        = pred,
        pct_diff         = round((obs - pred) / obs * 100, 1),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, Filter(Negate(is.null), rows))
}

#' OvsP calibration by subgroup for a Fine-Gray model (sheet wrapper).
#'
#' Thin wrapper around compute_fg_ovsp_subgroup_core() that reads subgroup
#' variables and category labels from the project worksheets.
#'
#' @param model                A Fine-Gray model object (coxph from finegray data).
#' @param model_data           Data frame used to fit the model.
#' @param variables_sheet      The project variables sheet (data.frame).
#' @param variable_details_sheet The project variable details sheet (data.frame).
#' @param horizon              Prediction horizon in years (default 5).
#' @param min_events_pct       Minimum observed CIF to include a subgroup (default 0.05).
#' @return A data frame with columns: variable, category, cat_label, n, observed, predicted, pct_diff.
#' @export
compute_fg_ovsp_subgroup <- function(model, model_data, variables_sheet,
                                     variable_details_sheet, horizon = 5,
                                     min_events_pct = 0.05) {
  sub_group_vars <- recodeflow:::select_vars_by_role("sub-group", variables_sheet)

  # Build category label lookup from variable_details_sheet
  category_labels <- setNames(
    as.character(variable_details_sheet$catLabelLong),
    paste0(variable_details_sheet$variable, "|", variable_details_sheet$recEnd)
  )

  compute_fg_ovsp_subgroup_core(
    model, model_data, sub_group_vars,
    category_labels, horizon, min_events_pct
  )
}

#' Decision Curve Analysis for a Fine-Gray model.
#'
#' Computes net benefit across a range of threshold probabilities using the
#' competing risks framework via dcurves::dca(). The "Treat All" and
#' "Treat None" reference strategies are included for comparison.
#'
#' @param model      coxph object (Fine-Gray).
#' @param data       Data frame with time_to_event, event_type, and predictors.
#' @param horizon    Prediction horizon in years (default 5).
#' @param thresholds Numeric vector of threshold probabilities (default seq(0, 0.5, by = 0.01)).
#' @return A data frame with columns: variable, label, threshold, net_benefit (and others from dcurves).
compute_fg_dca <- function(model, data, horizon = 5,
                           thresholds = seq(0, 0.5, by = 0.01), 
                           time_to_event = "time_to_event", event_type = "event_type") {
  # Some modification by JL
  names(data)[which(names(data) == time_to_event)] <- "time_to_event"
  names(data)[which(names(data) == event_type)] <- "event_type"
  
  data$.pred_risk <- compute_fg_pred_risk(model, data, horizon)

  dca_result <- dcurves::dca(
    survival::Surv(time_to_event, factor(event_type)) ~ .pred_risk,
    data       = data,
    time       = horizon,
    thresholds = thresholds
  )

  as.data.frame(dca_result$dca)
}

#' Plot decision curve analysis results.
#'
#' @param dca_data Data frame from compute_fg_dca().
#' @param title    Plot title.
#' @return ggplot object.
plot_dca <- function(dca_data, title = "Decision Curve Analysis") {
  # dcurves labels the model series by its column name (".pred_risk");
  # rename it to a reader-friendly label for the legend.
  dca_data$label <- ifelse(
    dca_data$label %in% c("Treat All", "Treat None"),
    as.character(dca_data$label),
    "DemPoRTv2"
  )
  dca_data$label <- factor(
    dca_data$label,
    levels = c("DemPoRTv2", "Treat All", "Treat None")
  )

  ggplot2::ggplot(dca_data,
                  ggplot2::aes(x = threshold, y = net_benefit,
                               color = label, linetype = label)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    ggplot2::scale_color_manual(
      values = c("DemPoRTv2" = "darkred", "Treat All" = "steelblue",
                 "Treat None" = "grey40")
    ) +
    ggplot2::scale_linetype_manual(
      values = c("DemPoRTv2" = "solid", "Treat All" = "dashed",
                 "Treat None" = "dotted")
    ) +
    ggplot2::coord_cartesian(ylim = c(-0.05, NA)) +
    ggplot2::labs(title = title, x = "Threshold Probability",
                  y = "Net Benefit", color = NULL, linetype = NULL) +
    ggplot2::theme_minimal()
}
