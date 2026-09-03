#' The outter layer function to call `calibration.surv.core` then `valProbSurvival.2` below for Cox PH model
#'
#' @param centerList The returned list from `fit_center`
#' @param newdata The dataset to apply the centered model
#' @param timeHorizon The specific times for predict survival
#' @param delta_recal The workaround in coxph to reflect the recalibration change
#' @param group The vector of group labels if want to assess calibration within subgroups
#' @param group_name A string vector, levels of group
#' @param g.group Number of groups when the subgroup variable is continuous
#' @param ... Optional arguments for insider functions
#' @return     A named vector of calibration stats
#' @export
calibration.surv <- function(
  centerList,
  newdata,
  timeHorizon = NULL,
  delta_recal = NULL,
  group = NULL,
  group_name = NULL,
  g.group = 4,
  ...
) {
  if (is.null(timeHorizon)) {
    stop("Specifit time horizon is required.")
  }
  # Save metrics in each subgroup and overall
  metrics <- c(
    "group",
    "Obs",
    "Pavg",
    "OE",
    "Slope",
    "ICI",
    "E50",
    "E90",
    "Emax"
  )
  res <- data.frame(matrix(ncol = length(metrics), nrow = 0))
  colnames(res) <- metrics

  if (!is.null(group)) {
    if (nrow(newdata) != length(group)) {
      stop("nrow(newdata) != length(group)")
    }
    if (is.factor(group) | is.logical(group) | is.character(group)) {
      if (!is.null(group_name)) {
        if (length(unique(group)) == length(group_name)) {
          newdata$group_fct <- as.character(group)
          newdata$group_fct <- factor(newdata$group_fct, levels = group_name)
        } else {
          stop(
            "categorical labels and group do not match, please check group_name"
          )
        }
      } else {
        newdata$group_fct <- factor(group)
        group_name <- levels(newdata$group_fct)
      }
    } else {
      newdata$group_fct <- cut2(group, g = g.group)
      group_name <- levels(newdata$group_fct)
    }

    # each group
    for (i in 1:length(group_name)) {
      group_i <- group_name[i]
      df_i <- newdata %>% filter(group_fct == group_i)
      res_i <- calibration.surv.core(
        centerList,
        newdata = df_i,
        timeHorizon = timeHorizon,
        delta_recal = delta_recal,
        print.plot = FALSE
      )
      vec <- unname(c(
        res_i$stats$Calibration$InTheLarge["Obs"],
        res_i$stats$Calibration$InTheLarge["Pavg"],
        res_i$stats$Calibration$InTheLarge["OE"],
        res_i$stats$Calibration$Slope["calibration slope"],
        res_i$stats$Calibration$Statistics["ICI"],
        res_i$stats$Calibration$Statistics["E50"],
        res_i$stats$Calibration$Statistics["E90"],
        res_i$stats$Calibration$Statistics["Emax"]
      ))
      res[i, 1] <- group_i # to make sure the rest values are numeric
      res[i, 2:ncol(res)] <- vec
    }
  }

  # Overall
  res_i <- calibration.surv.core(
    centerList,
    newdata = newdata,
    timeHorizon = timeHorizon,
    delta_recal = delta_recal,
    ...
  )
  vec <- unname(c(
    res_i$stats$Calibration$InTheLarge["Obs"],
    res_i$stats$Calibration$InTheLarge["Pavg"],
    res_i$stats$Calibration$InTheLarge["OE"],
    res_i$stats$Calibration$Slope["calibration slope"],
    res_i$stats$Calibration$Statistics["ICI"],
    res_i$stats$Calibration$Statistics["E50"],
    res_i$stats$Calibration$Statistics["E90"],
    res_i$stats$Calibration$Statistics["Emax"]
  ))
  i <- nrow(res) + 1
  res[i, 1] <- "Overall" # to make sure the rest values are numeric
  res[i, 2:ncol(res)] <- vec

  if (is.null(group)) {
    return(res_i)
  } else {
    return(res)
  }
}

# ========================================================================
#' The middle layer function to call `valProbSurvival.2` below for Cox PH model
#'
#' @param centerList The returned list from `fit_center`
#' @param newdata The dataset to apply the centered model.
#'                For fine-gray, it is the weighted data create by `survival::finegray`
#' @param timeHorizon The specific times for predict survival
#' @param delta_recal The workaround in coxph to reflect the recalibration change
#' @param print.plot Control if plot the calibration plot
#' @param ... Optional arguments for insider functions
#' @return     A named vector of calibration stats
#' @export
calibration.surv.core <- function(
  centerList,
  newdata,
  fg_time = NULL,
  fg_event = NULL,
  timeHorizon = NULL,
  delta_recal = NULL,
  print.plot = TRUE,
  ...
) {
  pred.lp <- center_predict(
    centerList,
    newdata = newdata,
    delta_recal = delta_recal,
    type = "lp"
  )
  pred.prob <- center_predict(
    centerList,
    newdata = newdata,
    delta_recal = delta_recal,
    type = "probability",
    times = timeHorizon
  )
  pred.prob.brier <- center_predict(
    centerList,
    newdata = newdata,
    delta_recal = delta_recal,
    type = "probability",
    times = timeHorizon - 0.01
  )

  # outcome
  formula <- centerList$fit.c$formula
  outcome <- as.character(formula)[2] # # parse the formula
  outcome_vec <- str_split_1(outcome, "[[:punct:]]+") %>% str_trim()
  outcome_vec <- outcome_vec[-c(1,length(outcome_vec))]

  if (length(outcome_vec) == 2) {
    # cox
    outcome.time <- as.numeric(unname(unlist(newdata[outcome_vec[1]])))
    outcome.event <- as.numeric(unname(unlist(newdata[outcome_vec[2]])))
    valProbSurvival.2(
      pred.lp,
      pred.prob,
      pred.prob.brier,
      outcome.time,
      outcome.event,
      print.plot = print.plot,
      timeHorizon = timeHorizon,
      ...
    )
  } else {
    # fine gray
    outcome.fgstart <- as.numeric(unname(unlist(newdata[outcome_vec[1]])))
    outcome.fgstop <- as.numeric(unname(unlist(newdata[outcome_vec[2]])))
    outcome.fgstatus <- as.numeric(unname(unlist(newdata[outcome_vec[3]])))

    valProbSurvival.fg(
      pred.lp,
      pred.prob,
      pred.prob.brier,
      outcome.fgstart,
      outcome.fgstop,
      outcome.fgstatus,
      print.plot = print.plot,
      timeHorizon = timeHorizon,
      ...
    )
  }
}
