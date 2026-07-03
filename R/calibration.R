source("R/calibration.con.R", echo = FALSE) 

#' A wrapper function to call functions for calibration curves and metrics
#'
#' @param prediction  The vector of predicted values.
#' @param outcome The vector of true outcome.
#' @param package The package for calibration curves for logistic regression.
#' @param ... Optional arguments for insider functions
#' @return     A named vector of calibration stats
#' @export
calibration <- function(
  prediction,
  outcome,
  package = "CalibrationCurves",
  ...
) {
  if (length(outcome) != length(prediction)) {
    stop("lengths of outcome and prediction do not agree")
  }

  if (length(unique(outcome)) == 2) {
    # Logistic regression: based on rms::val.prob
    if (min(prediction, na.rm = TRUE) < 0 | max(prediction, na.rm = TRUE) > 1) {
      stop(
        "Logistic regression: predicted probability should be between 0 and 1."
      )
    }

    # remove probability 0 or 1
    outcome <- outcome[prediction > 0 & prediction < 1]
    prediction <- prediction[prediction > 0 & prediction < 1]
    if (package == "CalibrationCurves") {
      stats <- CalibrationCurves::val.prob.ci.2(prediction, outcome, ...)
    } else {
      stats <- rms::val.prob(prediction, outcome, ...)
    }
  } else {
    stats <- calibration.con(prediction, outcome, ...)
  }

  return(stats)
}

