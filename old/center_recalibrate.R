require("stringr")
require("dplyr")
require("labelled")
require("rms")
require("survival")
require("riskRegression")
require("prodlim")

# ─────────────────────────────────────────────────────────────────────────────
#        The proper "transform_then_centering" pipeline
# ─────────────────────────────────────────────────────────────────────────────

#' Main function to fit a centered model based on the original model before centering
#' Currently consider the following types:
#' Categorical (with dummy)
#' Continuous, linear
#' Continuous, polynomial (raw and not orthogonal polynomials)
#' Continuous, rcs
#' Interaction terms
#'
#' @param fit  Fitted original model before centering
#' @param data The development data used to fit the original model.
#' @param time Only for survival model, time points for cumulative risk
#' @param fine_gray Flag if the model is fine gray
#' @return     A list containing the following components:
#' * `fit.c`: Fitted model with centering.
#' * `vars_con`: Vector of continuous variables.
#' * `vars_cat`: Vector of categorical variables.
#' * `labList`: List of levels of each categorical variable.
#' * `vars_rcs`: Vector of rcs variables.
#' * `predictors_rcs`: Calls of rcs terms.
#' * `rcsList`: List of rcs terms, including variable name and knot locations.
#' * `predictors_poly`: Calls of raw polinomial terms.
#' * `interaction_list`: List of interaction terms.
#' * `vec_mean`: Vector of mean values of each term.
#' * `outcome_mean`: Mean value of the outcome.
#' * `timeHorizon`: timeHorizon for survival models
#' * `maxtime`: Maximum time for survival models
#' @export
fit_center <- function(fit, data, timeHorizon = NULL, fine_gray = FALSE) {
  # extract the formula
  formula <- formula(fit)
  if (is.null(formula)) {
    stop("Cannot find the formula.")
  }

  # parse the formula
  outcome <- as.character(formula)[2] # outcome
  surv <- str_detect(outcome, "Surv") # check if it is a survival model

  # prepare the data
  res_prepared <- prepare_data_0(formula, data, fit = fit)
  data_processed <- res_prepared$data

  # ------- refit the centerred model -------
  if (surv) {
    fit.c <- update(fit, nocenter = NULL)
  } else {
    formula_C <- res_prepared$formula_prepared
    fit.c <- update(fit, formula = formula_C, data = data_processed)
  }

  # ------- outcome mean -------
  outcome_mean <- NULL
  maxtime <- NULL
  if (surv) {
    outcome_vec <- str_split_1(outcome, "[[:punct:]]+")[2:3] %>% str_trim()
    # estimated risk at a given time period using Kaplan-Meier survival function
    if (is.null(timeHorizon)) {
      timeHorizon <- median(data_processed[, outcome_vec[1]], na.rm = T)
    }
    formula_surv <- as.formula(paste0(outcome, " ~ 1"))
    surv_fit <- survfit(formula_surv, data = data_processed)
    S_KM <- summary(surv_fit, times = timeHorizon)$surv
    outcome_mean <- 1 - S_KM
    maxtime <- max(data_processed[, outcome_vec[1]], na.rm = T)
  } else {
    outcome_mean <- mean(data[, outcome], na.rm = T)
  }

  # ------- return the results-------
  return(list(
    fit.c = fit.c,
    vars_con = res_prepared$vars_con,
    vars_cat = res_prepared$vars_cat,
    labList = res_prepared$labList,
    vars_rcs = res_prepared$vars_rcs,
    predictors_rcs = res_prepared$predictors_rcs,
    rcsList = res_prepared$rcsList,
    predictors_poly = res_prepared$predictors_poly,
    interaction_list = res_prepared$interaction_list,
    vec_mean = res_prepared$vec_mean,
    outcome_mean = outcome_mean,
    timeHorizon = timeHorizon,
    maxtime = maxtime
  ))
}

# ----------------------------------------------------------
#' The function transports the centered model to `newdata`,
#' and computes predicted values.
#'
#' @param centerList The returned list from `fit_center`
#' @param newdata The dataset to apply the centered model
#' @param delta_recal The workaround in coxph to reflect the recalibration change
#' @param times The specific times for predict survival
#' @param fine_gray Flag if the model is fine gray
#' @return      A vector of predicted values
#' @export
predict_center <- function(
  centerList,
  newdata,
  delta_recal = NULL,
  times = NULL,
  fine_gray = FALSE,
  ...
) {
  # -------- determine model type --------
  fit.c <- centerList$fit.c
  outcome <- as.character(formula(fit.c))[2]
  surv <- str_detect(outcome, "Surv") # check if it is a survival model

  pred <- NULL
  if (!surv) {
    # Not a survival (coxph) model
    # -------- process newdata for centered model --------
    newdata_processed <- prepare_data(centerList, newdata, add.outcome = F)
    # -------- get predicted values --------
    pred <- predict(fit.c, newdata = newdata_processed, ...)
  } else {
    # Survival models did internal centering in the original model,
    # no additional data preparation required

    # get the "..." parameters
    captured_call <- match.call(expand.dots = FALSE)
    dots <- captured_call$...

    if (dots$type == "survival" | dots$type == "probability") {
      if (is.null(times)) {
        stop("To predict 'survival`, specific `times` are required.")
      } else {
        # use below method instead of `riskRegression::predictRisk` 
        # because this does not require `x=TRUE, y=TRUE` when fitting the model

        # First, calculate survival probability at specific times
        surv.fit <- survfit(fit.c, newdata = newdata)
        surv <- t(summary(surv.fit, times = times)$surv)
        colnames(surv) <- times
        if (!is.null(delta_recal)) {
          surv <- surv * exp(delta_recal)
        }
        
        if (dots$type == "probability") {
          # event probability at specific times
          pred <- 1 - surv
        } else if (dots$type == "survival") {
          # survival probability at specific times
          pred <- surv
        }
      }
    } else {
      pred <- predict(fit.c, newdata = newdata, ...)
      if (!is.null(delta_recal)) {
        # recalibrated survival model
        # update predicted values by incorporating the recalibrated shift
        if (dots$type == "lp") {
          pred <- pred + delta_recal
        } else if (dots$type == "risk") {
          pred <- pred * exp(delta_recal)
        }
      }
    }
  }

  # ------- return the results-------
  return(pred)
}

# ----------------------------------------------------------
#' The function recalibrates the centered model using `newdata`.
#'
#' @param centerList The returned list from `fit_center`
#' @param newdata The dataset to apply the centered model
#' @return      A list containing the following components, will be used in `recalibrate_center`:
#' * `centerList`: The returned list from `fit_center`.
#' * `delta_recal`: recalibration offset
#' @export
recalibrate_center <- function(centerList, newdata) {
  # model object of the centered model
  fit.c <- centerList$fit.c
  outcome <- as.character(formula(fit.c))[2]
  surv <- str_detect(outcome, "Surv") # check if it is a survival model

  delta_recal <- NULL
  if (all(class(fit.c) == "lm")) {
    # linear regression
    # ------- get processed data without centering -------
    newdata_processed <- prepare_data(
      centerList,
      newdata,
      center = F,
      add.outcome = F
    )
    # ------- center each term using mean values from the new dataset-------
    vars_mean <- names(newdata_processed)
    # get mean values
    vec_mean_newdata <- get_mean(newdata_processed, vars_mean)
    # centering
    newdata_processed <- step_center(
      newdata_processed,
      vars_mean,
      vec_mean_newdata
    )
    newdata_processed <- newdata_processed %>% select(-all_of(vars_mean))
    # -------- manually add the outcome variable --------
    newdata_processed <- bind_cols(
      newdata %>% select(all_of(outcome)),
      newdata_processed
    )
    # -------- outcome mean --------
    # in the original dataset
    outcome_mean_train <- centerList$outcome_mean
    # in the new dataset
    outcome_mean_test <- mean(newdata_processed[, outcome], na.rm = T)
    # -------- update `fit.c` in `centerList` --------
    fit.c.recalibrated <- fit.c
    fit.c.recalibrated$coefficients[1] <- outcome_mean_test
    centerList$fit.c <- fit.c.recalibrated
    # -------- calculate the recalibration offset --------
    vec_mean_train <- centerList$vec_mean
    delta_recal <- outcome_mean_test -
      vec_mean_newdata %*%
        coef(fit.c.recalibrated)[2:length(coef(fit.c.recalibrated))] -
      (outcome_mean_train -
        vec_mean_train %*% coef(fit.c)[2:length(coef(fit.c))])
    # -------- update `vec_mean` in `centerList` --------
    centerList$vec_mean <- vec_mean_newdata
  } else if (length(class(fit.o)) == 2) {
    # logistic regression
    # -------- get processed data with centering --------
    newdata_processed <- prepare_data(centerList, newdata)
    # -------- outcome mean --------
    # in the original dataset
    outcome_mean_train <- centerList$outcome_mean
    # in the new dataset
    outcome_mean_test <- mean(newdata_processed[, outcome], na.rm = T)
    # -------- update `fit.c` in `centerList` --------
    outcome_mean_train_logit <- log(
      outcome_mean_train / (1 - outcome_mean_train)
    )
    outcome_mean_test_logit <- log(
      outcome_mean_test / (1 - outcome_mean_test)
    )
    delta_recal <- outcome_mean_test_logit - outcome_mean_train_logit
    # update `fit.c` in `centerList` by updating the intercept
    fit.c.recalibrated <- fit.c
    fit.c.recalibrated$coefficients[1] <- fit.c$coefficients[1] + delta_recal
    centerList$fit.c <- fit.c.recalibrated
  } else if (surv) {
    # survival models
    # -------- outcome mean --------
    # in the original dataset
    outcome_mean_train <- centerList$outcome_mean
    # in the new dataset: estimated risk at a given time period using Kaplan-Meier survival function
    time <- centerList$time
    formula_surv <- as.formula(paste0(outcome, " ~ 1"))
    surv_fit <- survfit(formula_surv, data = newdata)
    S_KM <- summary(surv_fit, times = time)$surv
    outcome_mean_test <- 1 - S_KM

    # -------- calculate the shift: delta_recal --------
    # the model object cannot directly updated to reflect this change
    # Cox proportional hazards model, clog-log transforms
    # see doi:10.1016/j.jclinepi.2025.111895.
    outcome_mean_train_clog_log <- log(-log(1 - outcome_mean_train))
    outcome_mean_test_clog_log <- log(-log(1 - outcome_mean_test))
    delta_recal <- outcome_mean_test_clog_log - outcome_mean_train_clog_log
  }

  # ------- return the results-------
  return(list(
    centerList = centerList,
    delta_recal = delta_recal
  ))
}

# ----------------------------------------------------------
#' The function process the data for the centered model
#' base on known `centerList`, see `fit_center()`
#'
#' @param centerList The returned list from `fit_center`
#' @param data The dataset to apply the centered model
#' @param center Flag of whether doing centering, default is TRUE
#' @param add.outcome  Flag of whether including the outcome variable(s), default is TRUE
#' @return      The processed dataset
#' @export
prepare_data <- function(
  centerList,
  data,
  center = TRUE,
  add.outcome = TRUE
) {
  # ------- prepare data -------
  # categorical predictors, dummy
  vars_cat <- centerList$vars_cat
  if (length(vars_cat) > 0) {
    data <- step_dummy(data, vars_cat, labList = centerList$labList)
    data <- data %>% select(-all_of(vars_cat))
  }

  # continuous predictors, rcs terms
  rcsList <- centerList$rcsList
  vars_rcs <- centerList$vars_rcs
  predictors_rcs <- centerList$predictors_rcs
  if (length(vars_rcs) > 0) {
    data <- step_rcs(data, vars_rcs, rcsList, predictors_rcs) %>%
      select(-all_of(vars_rcs))
  }

  # continuous predictors, raw polinomial terms
  predictors_poly <- centerList$predictors_poly
  if (length(predictors_poly) > 0) {
    res_poly <- step_poly(data, predictors_poly)
    var_poly <- res_poly$var_poly
    data <- res_poly$newdata %>% select(-all_of(var_poly))
  }

  # add interaction terms
  interaction_list <- centerList$interaction_list
  data <- step_interaction(data, interaction_list)

  # select and reorder predictors
  data_processed <- data[, names(centerList$vec_mean)]

  if (center) {
    # ------- center each term using mean values from the original dataset-------
    # get mean values
    vec_mean <- centerList$vec_mean
    vars_mean <- names(vec_mean)
    # centering
    data_processed <- step_center(data_processed, vars_mean, vec_mean)
    data_processed <- data_processed %>% select(-all_of(vars_mean))
  }

  # model object of the centered model
  fit.c <- centerList$fit.c
  outcome <- as.character(formula(fit.c))[2]
  surv <- str_detect(outcome, "Surv") # check if it is a survival model

  if (add.outcome) {
    if (surv) {
      temp <- str_split_1(outcome, "[[:punct:]]+")
      outcome_vec <- temp[2:(length(temp) - 1)] %>% str_trim()
      data_processed <- bind_cols(
        data %>% select(all_of(outcome_vec)),
        data_processed
      )
    } else {
      data_processed <- bind_cols(
        data %>% select(all_of(outcome)),
        data_processed
      )
    }
  }

  return(data_processed)
}

# ----------------------------------------------------------
#' The function process the data using `formula`,
#' used within `fit_center()`,
#' and prepare data for `riskRegression::FGR()`
#'
#' @param formula The `formula`
#' @param data The dataset
#' @param center Flag of whether doing centering, default is TRUE
#' @param add.outcome  Flag of whether including the outcome variable(s), default is TRUE
#' @param fit Fit of the original model, if known, default is NULL
#' @return      The processed dataset
#' @export
prepare_data_0 <- function(
  formula,
  data,
  center = TRUE,
  add.outcome = TRUE,
  fit = NULL
) {
  # parse the formula
  outcome <- as.character(formula)[2] # outcome
  surv <- str_detect(outcome, "Surv") # check if it is a survival model

  predictors <- str_split(
    as.character(formula)[3],
    fixed("+"),
    simplify = TRUE
  ) %>%
    str_trim()

  # get all the terms and categorical variable labels in the original model
  terms <- NULL
  labList <- NULL
  if (!is.null(fit)) {
    if (surv) {
      terms <- names(coef(fit))
    } else {
      terms <- names(coef(fit))[-1]
    }
    labList <- fit$xlevels
  }

  # ------- extract predictors and their terms based on their types -------
  # interaction terms
  interactions <- predictors[str_detect(predictors, "\\*")]
  interaction_list <- list()
  if (length(interactions) > 0) {
    for (i_interaction in seq_len(length(interactions))) {
      interaction_list[[i_interaction]] <- str_split_1(
        interactions[[i_interaction]],
        "\\*"
      ) %>%
        str_trim()
    }
  }
  # update predictors
  predictors <- unique(c(
    predictors[!(predictors %in% interactions)],
    unlist(interaction_list)
  ))

  # continuous predictors, rcs terms
  rcsList <- list()
  predictors_rcs <- predictors[str_detect(predictors, "rcs")]
  vars_rcs <- c()
  if (length(predictors_rcs) > 0) {
    for (temp in predictors_rcs) {
      var <- str_split_1(temp, "[[:punct:]]+")[2]
      k <- as.numeric(str_split_1(temp, "[[:punct:]]+")[3])
      vars_rcs <- c(vars_rcs, var)

      # fit rcs
      rcs.fit <- rcs(data[, var], k)

      # save knot locations
      rcsList[[var]] <- attributes(rcs.fit)$parms

      # extract rcs terms
      res <- labelled::remove_attributes(
        rcs.fit,
        names(attributes(rcs.fit))[3:length(names(attributes(rcs.fit)))]
      ) # don't print the attributes
      colnames(res) <- str_replace_all(colnames(res), "data", paste0(temp, var))

      data <- bind_cols(data, res) %>% select(-all_of(var))
    }
  }

  # continuous predictors, raw polinomial terms
  predictors_poly <- predictors[str_detect(predictors, "poly")]
  if (length(predictors_poly) > 0) {
    for (temp in predictors_poly) {
      var <- str_split_1(temp, "[[:punct:]]+")[2]
      degree <- as.numeric(str_split_1(temp, "[[:punct:]]+")[3])

      # fit poly
      res <- poly(data[, var], degree, raw = TRUE)
      colnames(res) <- paste0(temp, colnames(res))

      data <- bind_cols(data, res) %>% select(-all_of(var))
    }
  }

  # continuous predictors, linear term
  if (!is.null(terms)) {
    vars_con <- predictors[predictors %in% terms]
  } else {
    vars_con <- names(data)[sapply(data, is.numeric)]
    vars_con <- vars_con[which(vars_con %in% predictors)]
  }

  # categorical predictors, dummy
  vars_cat <- predictors[which(
    !(predictors %in% c(vars_con, predictors_rcs, predictors_poly))
  )]
  if (length(vars_cat) > 0) {
    data <- step_dummy(data, vars_cat, labList)
    data <- data %>% select(-all_of(vars_cat))
  }

  # add interaction terms
  data <- step_interaction(data, interaction_list)

  # select predictors, and reorder them to match the original model if needed
  data_processed <- data.frame(temp = rep(1, nrow(data)))
  if (!is.null(terms)) {
    data_processed <- data[, terms]
  } else {
    for (var in predictors) {
      df_temp <- data %>% select(all_of(str_subset(names(data), fixed(var))))
      vars <- names(df_temp)[which(
        !(names(df_temp) %in% names(data_processed))
      )]
      data_processed <- bind_cols(
        data_processed,
        df_temp %>% select(all_of(vars))
      )
    }
    data_processed <- data_processed %>% select(-temp)
  }

  # center each variable
  vec_mean <- NULL
  if (center) {
    vars_mean <- names(data_processed)
    # get mean values
    vec_mean <- get_mean(data_processed, vars_mean)
    # centering
    data_processed <- step_center(data_processed, vars_mean, vec_mean)
    data_processed <- data_processed %>% select(-all_of(vars_mean))
  }

  # variable names of all the predictors
  name_predictor <- names(data_processed)

  # update formula
  formula_prepared <- as.formula(paste0(
    outcome,
    " ~ ",
    paste(paste0(paste0("`", name_predictor), "`"), collapse = " + ")
  ))

  # model object of the centered model
  outcome <- as.character(formula)[2]
  surv <- str_detect(outcome, "Surv") # check if it is a survival model

  if (add.outcome) {
    if (surv) {
      temp <- str_split_1(outcome, "[[:punct:]]+")
      outcome_vec <- temp[2:(length(temp) - 1)] %>% str_trim()
      data_processed <- bind_cols(
        data %>% select(all_of(outcome_vec)),
        data_processed
      )
    } else {
      data_processed <- bind_cols(
        data %>% select(all_of(outcome)),
        data_processed
      )
    }
  }

  return(list(
    data = data_processed,
    formula_prepared = formula_prepared,
    vars_con = vars_con,
    vars_cat = vars_cat,
    labList = labList,
    vars_rcs = vars_rcs,
    predictors_rcs = predictors_rcs,
    rcsList = rcsList,
    predictors_poly = predictors_poly,
    interaction_list = interaction_list,
    vec_mean = vec_mean
  ))
}

# ----------------------------------------------------------
# dummy variables
step_dummy <- function(df, vars_cat, labList = NULL) {
  if (!is.null(labList)) {
    if (length(labList) != length(vars_cat)) {
      stop("labList must have the same length as vars_cat")
    }
  } else {
    labList <- list()
    for (var in vars_cat) {
      if (!is.null(levels(df[, var]))) {
        cats <- levels(df[, var])
      } else {
        cats <- sort(unique(df[, var]))
      }
      labList[[var]] <- cats
    }
  }

  for (var in vars_cat) {
    cats <- labList[[var]]
    for (cat in cats[-1]) {
      df$new <- ifelse(df[, var] == cat, 1, 0)
      names(df)[ncol(df)] <- paste0(var, cat)
    }
  }
  return(df)
}

# ----------------------------------------------------------
# rcs
step_rcs <- function(df, vars_rcs, rcsList, predictors_rcs) {
  for (i_var in seq_len(length(vars_rcs))) {
    var <- vars_rcs[i_var]
    knots <- rcsList[[var]]
    k <- length(knots)
    rcs.fit <- rcs(df[, var], knots)

    # extract rcs terms
    res <- labelled::remove_attributes(
      rcs.fit,
      names(attributes(rcs.fit))[3:length(names(attributes(rcs.fit)))]
    ) # don't print the attributes
    colnames(res) <- str_replace_all(
      colnames(res),
      "df",
      paste0(predictors_rcs[i_var], var)
    )

    df <- bind_cols(df, res)
  }
  return(df)
}

# ----------------------------------------------------------
# polynominal
step_poly <- function(df, predictors_poly) {
  var_poly <- c()
  for (temp in predictors_poly) {
    var <- str_split_1(temp, "[[:punct:]]+")[2]
    degree <- as.numeric(str_split_1(temp, "[[:punct:]]+")[3])

    # fit poly
    res <- poly(df[, var], degree, raw = TRUE)
    colnames(res) <- paste0(temp, colnames(res))

    df <- bind_cols(df, res)
    var_poly <- c(var_poly, var)
  }
  return(list(newdata = df, var_poly = var_poly))
}

# ----------------------------------------------------------
# add interaction terms
step_interaction <- function(df, interaction_list) {
  df_old <- df
  for (interaction in interaction_list) {
    var1 <- interaction[1]
    var2 <- interaction[2]

    # in case of categorical variable, or rcs, list all the possible combinations
    vars1 <- names(df_old)[stringr::str_detect(names(df_old), fixed(var1))]
    if (length(vars1) > 1) {
      vars1 <- vars1[which(vars1 != var1)]
    }
    vars2 <- names(df_old)[str_detect(names(df_old), fixed(var2))]
    if (length(vars2) > 1) {
      vars2 <- vars2[which(vars2 != var2)]
    }

    for (term2 in vars2) {
      for (term1 in vars1) {
        df$new <- df[, term1] * df[, term2]
        names(df)[ncol(df)] <- paste0(term1, ":", term2)
      }
    }
  }
  return(df)
}

# ----------------------------------------------------------
# get mean values from the original dataset
get_mean <- function(df, vars_mean) {
  means <- list()
  for (var in vars_mean) {
    if (is.numeric(df[, var])) {
      means[[var]] <- mean(df[, var])
    }
  }
  return(unlist(means))
}

# ----------------------------------------------------------
# centering all variables, including the interaction terms
step_center <- function(df, vars_center, means) {
  for (var in vars_center) {
    df$new <- df[, var] - means[var]
    names(df)[ncol(df)] <- paste0(var, "_C")
  }
  return(df)
}

# ─────────────────────────────────────────────────────────────────────────────
# Claude: riskRegression integration for Fine-Gray models fitted with
#         survival::finegray + survival::coxph.
#
# Defines an S3 class `fg_coxph` with a `predictRisk` method so that
# riskRegression::predictRisk() and riskRegression::plotCalibration() work
# on these models (including centered and recalibrated variants).
#
# Requires: riskRegression and its dependency prodlim to be loaded by the caller.
# ─────────────────────────────────────────────────────────────────────────────

# ----------------------------------------------------------
#' Wrap a Fine-Gray coxph model for use with riskRegression
#'
#' Creates an `fg_coxph` object that provides the `predictRisk` interface
#' expected by `riskRegression::predictRisk()` and `riskRegression::plotCalibration()`.
#'
#' @param fit        A `coxph` object fitted on `finegray()`-weighted data.
#' @param time_var   Name of the time variable in the **original** (non-fg) data.
#' @param event_var  Name of the multi-level event variable in the original data
#'                   (0 = censored, 1 = event of interest, 2+ = competing events).
#' @param centerList Optional `centerList` from `fit_center()`. When supplied,
#'                   `newdata` is preprocessed through the centering pipeline
#'                   before prediction.
#' @param delta_recal Optional scalar recalibration offset from `recalibrate_center()`.
#'                    Applied as CIF * exp(delta_recal), clipped to [0, 1].
#' @param cause      Integer code for the event of interest (default 1).
#' @return An object of class `"fg_coxph"`.
#' @export
# Claude
as_fg_coxph <- function(
  fit,
  time_var,
  event_var,
  centerList = NULL,
  delta_recal = NULL,
  cause = 1
) {
  structure(
    list(
      fit = fit,
      centerList = centerList,
      delta_recal = delta_recal,
      cause = cause,
      # Hist() formula is what riskRegression::Score / plotCalibration uses to
      # locate the competing-risk outcome columns in newdata
      formula = as.formula(
        paste0("Hist(", time_var, ", ", event_var, ") ~ 1")
      )
    ),
    class = "fg_coxph"
  )
}

# ----------------------------------------------------------
# Claude: formula method — riskRegression calls formula(object) to extract the
# outcome formula when building calibration assessments via Score()
#' @export
formula.fg_coxph <- function(x, ...) x$formula

# ----------------------------------------------------------
# Claude: print method for quick inspection
#' @export
print.fg_coxph <- function(x, ...) {
  cat("Fine-Gray coxph wrapper (fg_coxph)\n")
  cat("Outcome :", deparse(x$formula[[2]]), "\n")
  cat("Cause   :", x$cause, "\n")
  cat(
    "Centered:",
    !is.null(x$centerList),
    " | Recalibrated:",
    !is.null(x$delta_recal),
    "\n"
  )
  invisible(x)
}

# ----------------------------------------------------------
#' predictRisk S3 method for fg_coxph
#'
#' Called automatically by `riskRegression::predictRisk()` and internally by
#' `riskRegression::plotCalibration()`. Returns predicted CIFs at the requested
#' time points.
#'
#' @param object  An `fg_coxph` object created by `as_fg_coxph()`.
#' @param newdata Data frame of new subjects using the **original** data layout
#'                (not the finegray-expanded data).
#' @param times   Numeric vector of time points at which to evaluate the CIF.
#' @param cause   Ignored — the cause is fixed at object creation; kept for
#'                interface compatibility with riskRegression.
#' @return An n_subjects × n_times numeric matrix of predicted CIFs.
#' @export
# Claude
predictRisk.fg_coxph <- function(object, newdata, times, cause, ...) {
  # preprocess newdata through the centering pipeline for centered models
  if (!is.null(object$centerList)) {
    newdata_processed <- prepare_data(
      object$centerList,
      newdata,
      add.outcome = FALSE
    )
  } else {
    newdata_processed <- newdata
  }

  surv_pred <- survfit(object$fit, newdata = newdata_processed)

  # summary() with a scalar `times` returns a vector of length n_subjects;
  # loop over time points to build an n_subjects × n_times matrix
  n <- nrow(newdata)
  cif_mat <- matrix(NA_real_, nrow = n, ncol = length(times))
  for (j in seq_along(times)) {
    s <- summary(surv_pred, times = times[j], extend = TRUE)$surv
    cif_mat[, j] <- 1 - s
  }

  # apply recalibration shift and clip to valid probability range
  if (!is.null(object$delta_recal)) {
    cif_mat <- pmin(pmax(cif_mat * exp(object$delta_recal), 0), 1)
  }

  cif_mat
}

# ─────────────────────────────────────────────────────────────────────────────
# Claude: Score()-based bridge for Fine-Gray models and riskRegression.
#
# `as_fg_pred` wraps a fitted model (bare coxph or centerList) into a minimal
# S3 object. `predictRisk.fg_pred` is then dispatched by riskRegression::Score()
# to supply predicted CIFs, so `plotCalibration()` can be called on the Score
# result rather than on the model object directly.
#
# Typical workflow:
#   test$etime <- etime          # outcome columns must exist in the data frame
#   test$event <- event          # passed to Score()
#
#   s <- Score(
#     list("Original"    = as_fg_pred(fit.o,       time = horizon),
#          "Centered"    = as_fg_pred(centerList,   time = horizon),
#          "Recalibrated"= as_fg_pred(recalibrateList$centerList, time = horizon,
#                                     delta_recal = recalibrateList$delta_recal)),
#     formula = Hist(etime, event) ~ 1,
#     data    = test, times = horizon, cause = 1, plots = "calibration"
#   )
#   plotCalibration(s)
# ─────────────────────────────────────────────────────────────────────────────

# ----------------------------------------------------------
#' Wrap a Fine-Gray model for use with riskRegression::Score()
#'
#' Creates an `fg_pred` object whose `predictRisk` method is dispatched by
#' `riskRegression::Score()`. Use the `Score()` result with `plotCalibration()`.
#'
#' @param fit_or_centerList Either a `coxph` fit (original, uncentered model) or
#'   a `centerList` returned by `fit_center()` (centered / recalibrated model).
#'   When a `centerList` is supplied, `newdata` is preprocessed through the
#'   centering pipeline before each prediction call.
#' @param time  Scalar time horizon at which the CIF is evaluated.
#' @param delta_recal  Optional recalibration offset from `recalibrate_center()`.
#'   Applied as CIF * exp(delta_recal), clipped to [0, 1].
#' @return An object of class `"fg_pred"`.
#' @export
# Claude
as_fg_pred <- function(fit_or_centerList, time, delta_recal = NULL) {
  if (inherits(fit_or_centerList, "coxph")) {
    centerList <- NULL
    fit <- fit_or_centerList
  } else if (is.list(fit_or_centerList) && !is.null(fit_or_centerList$fit.c)) {
    centerList <- fit_or_centerList
    fit <- centerList$fit.c
  } else {
    stop(
      "fit_or_centerList must be a coxph fit or a centerList from fit_center()"
    )
  }
  structure(
    list(
      fit = fit,
      centerList = centerList,
      time = time,
      delta_recal = delta_recal
    ),
    class = "fg_pred"
  )
}

# ----------------------------------------------------------
#' predictRisk S3 method for fg_pred
#'
#' Called automatically by `riskRegression::Score()`. Returns an
#' n_subjects x n_times matrix of predicted CIFs.
#'
#' @param object  An `fg_pred` object from `as_fg_pred()`.
#' @param newdata Data frame of new subjects (original layout, not finegray-expanded).
#'   Must contain the outcome columns (`time_var`, `event_var`) specified in the
#'   `Score()` formula so that `Score()` can assess observed outcomes.
#' @param times   Time points requested by `Score()`.
#' @param cause   Ignored; kept for interface compatibility.
#' @return An n_subjects x n_times numeric matrix of predicted CIFs.
#' @export
# Claude
predictRisk.fg_pred <- function(object, newdata, times, cause, ...) {
  # preprocess through centering pipeline for centered models
  if (!is.null(object$centerList)) {
    newdata_processed <- prepare_data(
      object$centerList,
      newdata,
      add.outcome = FALSE
    )
  } else {
    newdata_processed <- newdata
  }

  surv_pred <- survfit(object$fit, newdata = newdata_processed)

  n <- nrow(newdata)
  cif_mat <- matrix(NA_real_, nrow = n, ncol = length(times))
  for (j in seq_along(times)) {
    s <- summary(surv_pred, times = times[j], extend = TRUE)$surv
    cif_mat[, j] <- 1 - s
  }

  if (!is.null(object$delta_recal)) {
    cif_mat <- pmin(pmax(cif_mat * exp(object$delta_recal), 0), 1)
  }

  cif_mat
}