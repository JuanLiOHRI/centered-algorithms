#' The function transports the centered model to `newdata`,
#' and computes predicted values.
#'
#' @param centerList The returned list from `center_fit`
#' @param newdata The dataset to apply the centered model
#' @param delta_recal The workaround in coxph to reflect the recalibration change
#' @param times The specific times for predict survival
#' @return      A vector of predicted values
#' @export
center_predict <- function(
  centerList,
  newdata,
  delta_recal = NULL,
  times = NULL,
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
    newdata_processed <- center_prepare(centerList, newdata, add.outcome = F)
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
        stop("To predict 'survival` or `probability`, specific `times` are required.")
      } else {
        # First, calculate survival probability at specific times
        surv.fit <- survfit(fit.c, newdata = newdata)
        surv <- t(summary(surv.fit, times = times)$surv)
        colnames(surv) <- times
        if (!is.null(delta_recal)) {
          surv <- surv ^ exp(delta_recal)
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