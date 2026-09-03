#' The function process the data for the centered model base on known `fit` object
#' Normally used inside `center_fit()`
#'
#' @param fit Fit of the original model
#' @param data The dataset
#' @param center Flag of whether doing centering, default is TRUE
#' @param add.outcome  Flag of whether including the outcome variable(s), default is TRUE
#' @param means_ext Vector of mean values of an external dataset
#' @return      The processed dataset
#' @export
center_prepare_0 <- function(
  fit,
  data,
  center = TRUE,
  add.outcome = TRUE,
  means_ext = NULL
) {
  # extract the formula
  formula <- formula(fit)
  if (is.null(formula)) {
    stop("Cannot find the formula.")
  }

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
  if (surv) {
    terms <- names(coef(fit))
  } else {
    terms <- names(coef(fit))[-1]
  }
  labList <- fit$xlevels

  # check means_ext
  if (!is.null(means_ext)) {
    if (length(means_ext) != length(terms)) {
      stop("length of means_ext does not match with number of terms.")
    }
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

  # continuous predictors, raw polynomial terms
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
    if (is.null(names(data_processed))) {
      names(data_processed) <- terms
    }
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
    if (is.null(means_ext)) {
      # get mean values
      vec_mean <- get_mean(data_processed, vars_mean)
    } else {
      vec_mean <- means_ext
    }
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
#' The function process the data for the centered model
#' base on known `centerList`, see `center_fit()`
#' Normally used inside `center_predict()` and `center_recalibrate()`
#'
#' @param centerList The returned list from `center_fit`
#' @param data The dataset to apply the centered model
#' @param center Flag of whether doing centering, default is TRUE
#' @param add.outcome  Flag of whether including the outcome variable(s), default is TRUE
#' @param means_ext Vector of mean values of an external dataset
#' @return      The processed dataset
#' @export
center_prepare <- function(
  centerList,
  data,
  center = TRUE,
  add.outcome = TRUE,
  means_ext = NULL
) {
  # check means_ext
  if (!is.null(means_ext)) {
    if (length(means_ext) != length(centerList$vec_mean)) {
      stop("length of means_ext does not match centerList$vec_mean.")
    }
  }

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

  # continuous predictors, raw polynomial terms
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
    # get mean values
    vec_mean <- centerList$vec_mean
    vars_mean <- names(vec_mean)
    if (!is.null(means_ext)) {
      vec_mean <- means_ext
    }
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