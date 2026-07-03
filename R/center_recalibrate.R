require(dplyr)
require(stringr)
require(survival)

source("R/center_prepare.R", echo = FALSE) 

#' The function recalibrates the centered model using `newdata`.
#'
#' @param centerList The returned list from `center_fit`
#' @param newdata The dataset to apply the centered model
#' @param timeHorizon Only for survival models, the time point of interest
#' @return      A list containing the following components, will be used in `center_recalibrate`:
#' * `centerList`: The returned list from `center_fit`.
#' * `delta_recal`: recalibration offset
#' @export
center_recalibrate <- function(centerList, newdata, timeHorizon = NULL) {
  # model object of the centered model
  fit.c <- centerList$fit.c
  outcome <- as.character(formula(fit.c))[2]
  surv <- str_detect(outcome, "Surv") # check if it is a survival model

  # check timeHorizon
  if (surv) {
    if (is.null(timeHorizon)) {
      stop("Recalibrating the survival models requires a specific `timeHorizon`.")
    } else {
      times <- as.numeric(colnames(centerList$outcome_mean))
      if (!(timeHorizon %in% times)) {
        stop(
          "`timeHorizon` is not within the pre-specify `times` used in `center_fit`."
        )
      }
    }
  }

  delta_recal <- NULL
  if (all(class(fit.c) == "lm")) {
    # linear regression
    # ------- get processed data without centering -------
    newdata_processed <- center_prepare(
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
    outcome_mean_deriv <- centerList$outcome_mean
    # in the new dataset
    outcome_mean_target <- mean(newdata_processed[, outcome], na.rm = T)
    # -------- update `fit.c` in `centerList` --------
    fit.c.recalibrated <- fit.c
    fit.c.recalibrated$coefficients[1] <- outcome_mean_target
    centerList$fit.c <- fit.c.recalibrated
    # -------- calculate the recalibration offset --------
    vec_mean_deriv <- centerList$vec_mean
    delta_recal <- outcome_mean_target -
      vec_mean_newdata %*%
        coef(fit.c.recalibrated)[2:length(coef(fit.c.recalibrated))] -
      (outcome_mean_deriv -
        vec_mean_deriv %*% coef(fit.c)[2:length(coef(fit.c))])
    # -------- update `vec_mean` in `centerList` --------
    centerList$vec_mean <- vec_mean_newdata
  } else if (length(class(fit.o)) == 2) {
    # logistic regression
    # -------- get processed data with centering --------
    newdata_processed <- center_prepare(centerList, newdata)
    # -------- outcome mean --------
    # in the original dataset
    outcome_mean_deriv <- centerList$outcome_mean
    # in the new dataset
    outcome_mean_target <- mean(newdata_processed[, outcome], na.rm = T)
    # -------- update `fit.c` in `centerList` --------
    outcome_mean_deriv_logit <- log(
      outcome_mean_deriv / (1 - outcome_mean_deriv)
    )
    outcome_mean_target_logit <- log(
      outcome_mean_target / (1 - outcome_mean_target)
    )
    delta_recal <- outcome_mean_target_logit - outcome_mean_deriv_logit
    # update `fit.c` in `centerList` by updating the intercept
    fit.c.recalibrated <- fit.c
    fit.c.recalibrated$coefficients[1] <- fit.c$coefficients[1] + delta_recal
    centerList$fit.c <- fit.c.recalibrated
  } else if (surv) {
    # survival models
    # -------- outcome mean --------
    # in the original dataset
    outcome_mean_vec <- centerList$outcome_mean
    outcome_mean_deriv <- as.numeric(outcome_mean_vec[,which(colnames(outcome_mean_vec) == timeHorizon)])
    # in the new dataset
    formula_surv <- as.formula(paste0(outcome, " ~ 1"))
    surv_fit_target <- survfit(formula_surv, data = newdata)
    S_KM_target <- summary(surv_fit_target, times = timeHorizon)$surv
    outcome_mean_target <- 1 - S_KM_target

    # -------- calculate the shift: delta_recal --------
    # the model object cannot directly updated to reflect this change
    # Cox proportional hazards model, clog-log transforms
    # see doi:10.1016/j.jclinepi.2025.111895.
    outcome_mean_deriv_clog_log <- log(-log(1 - outcome_mean_deriv))
    outcome_mean_target_clog_log <- log(-log(1 - outcome_mean_target))
    delta_recal <- outcome_mean_target_clog_log - outcome_mean_deriv_clog_log
  }

  # ------- return the results-------
  return(list(
    centerList = centerList,
    delta_recal = delta_recal
  ))
}

