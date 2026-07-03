require(dplyr)
require(stringr)
require(survival)

source("R/center_prepare.R", echo = FALSE) 

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
#' @param times Only for cox and fine-gray, the time point(s) at which to calculate the `outcome_mean`
#' @param data_fg Only for fine-gray, the extended dataset by `survival::finegray`
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
#' * `outcome_mean`: Mean value of the outcome in `data` (scalar for linear and logistic regression, vector for cox and fine-gray).
#' * `maxtime`: Maximum time for survival models
#' @export
center_fit <- function(fit, data, times = NULL, data_fg = NULL) {
  # extract the formula
  formula <- formula(fit)
  if (is.null(formula)) {
    stop("Cannot find the formula.")
  }

  # parse the formula
  outcome <- as.character(formula)[2] # outcome
  surv <- str_detect(outcome, "Surv") # check if it is a survival model
  outcome_vec <- str_split_1(outcome, "[[:punct:]]+") %>% str_trim()
  outcome_vec <- outcome_vec[-c(1, length(outcome_vec))]

  # prepare the data
  if (length(outcome_vec) == 3) {
    # fine-gray: fgstart, fgstop, fgstatus
    if (is.null(data_fg)) {
      stop("The fine-gray model requires both data and data_fg.")
    }
    # get mean values of `data`
    res_prepared_0 <- center_prepare_0(fit, data)
    # recentering `data_fg` using the mean values above
    res_prepared <- center_prepare_0(fit, data_fg, means_ext = res_prepared_0$vec_mean)
  } else {
    res_prepared <- center_prepare_0(fit, data)
  }
  data_processed <- res_prepared$data

  # ------- refit the centerred model -------
  if (surv) {
    if (length(outcome_vec) == 2) {
      # cox: time, tatus 
      # only need to make sure `nocenter = NULL`
      fit.c <- update(fit, nocenter = NULL)
    } else {
      # fine-gray: fgstart, fgstop, fgstatus
      # Need to manually center data_fg with mean values in data
      fit.c <- update(fit, formula = res_prepared$formula_prepared, data = data_processed, nocenter = NULL)
    }
  } else {
    # linear, logistic
    # Need to manually center all terms
    fit.c <- update(fit, formula = res_prepared$formula_prepared, data = data_processed)
  }

  # ------- outcome mean -------
  outcome_mean <- NULL
  maxtime <- NULL
  if (surv) {
    formula_surv <- as.formula(paste0(outcome, " ~ 1"))
    # NOTE: for fine-gray, this is a Kaplan-Meier (KM) estimator on the prepared `data_fg`
    # May not be correct.
    # Recalibration of fine-gray is still under development!!
    surv_fit <- survfit(formula_surv, data = data_processed)
    # probability of event at specific `times`
    pred.prob <- 1 - t(summary(surv_fit, times = times)$surv)
    colnames(pred.prob) <- times
    outcome_mean <- pred.prob
    # maxtime
    if (length(outcome_vec) == 2) {
      # cox: time, tatus
      maxtime <- max(data_processed[, outcome_vec[1]], na.rm = T)
    } else {
      # fine-gray: fgstart, fgstop, fgstatus
      maxtime <- max(data_processed[, outcome_vec[2]], na.rm = T)
    }
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
    maxtime = maxtime
  ))
}