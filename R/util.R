require("dplyr")
# implement the formula of rcs components
get_rcs <- function(x, knots) {
  k <- length(knots)
  res <- data.frame(rcs.1 = x)
  for (j in 1:(k - 2)) {
    kd <- (knots[k] - knots[1])^(2 / 3) # (knot_k-knot_1)^(2/3)

    vec1 <- x - knots[j] # X-knot_j

    val2 <- knots[k - 1] - knots[j] # knot_{k-1}-knot_j
    vec2 <- x - knots[k] # X-knot_k

    val3 <- knots[k] - knots[j] # knot_k-knot_j
    vec3 <- x - knots[k - 1] # X-knot_{k-1}
    val4 <- knots[k] - knots[k - 1] # knot_k-knot_{k-1}

    vec_j <- pmax(vec1 / kd, 0)^3 +
      (val2 * pmax(vec2 / kd, 0)^3 - val3 * pmax(vec3 / kd, 0)^3) / val4

    # View(Hmisc::rcspline.eval): line 111-114
    # vec_j <- pmax((x - knots[j])/kd, 0)^3 +
    #   ((knots[k-1] - knots[j]) * pmax((x - knots[k])/kd, 0)^3 -
    #      (knots[k] - knots[j]) * (pmax((x - knots[k-1])/kd, 0)^power))/(knots[k] - knots[k-1])

    res$new <- vec_j
    names(res)[ncol(res)] <- paste0("rcs.", j + 1)
  }
  res <- as.matrix(res)
  return(res)
}

# --------------------------
# define function for data processing

# -------- dummy variables --------
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
      names(df)[ncol(df)] <- paste0(var, "_", cat)
    }
  }
  return(df)
}

# -------- rcs --------
step_rcs <- function(df, vars_rcs, knots_list) {
  for (var in vars_rcs) {
    knots <- knots_list[[var]]
    k <- length(knots)
    res <- get_rcs(df[, var], knots)
    colnames(res) <- paste0(var, "_rcs_", 1:(k - 1))
    df <- bind_cols(df, res)
  }
  return(df)
}

# -------- add interaction terms --------
step_interaction <- function(df, interaction_list) {
  df_old <- df
  for (interaction in interaction_list) {
    var1 <- interaction[1]
    var2 <- interaction[2]

    # in case of categorical variable, or rcs, list all the possible combinations
    vars1 <- names(df_old)[stringr::str_detect(names(df_old), var1)]
    if (length(vars1) > 1) {
      vars1 <- vars1[which(vars1 != var1)]
    }
    vars2 <- names(df_old)[stringr::str_detect(names(df_old), var2)]
    if (length(vars2) > 1) {
      vars2 <- vars2[which(vars2 != var2)]
    }

    for (term2 in vars2) {
      for (term1 in vars1) {
        df$new <- df[, term1] * df[, term2]
        names(df)[ncol(df)] <- paste0(term1, "_by_", term2)
      }
    }
  }
  return(df)
}

# -------- get mean values from the original dataset --------
get_mean <- function(df, vars_mean) {
  means <- list()
  for (var in vars_mean) {
    if (is.numeric(df[, var])) {
      means[[var]] <- mean(df[, var])
    }
  }
  return(unlist(means))
}

# -------- centering all variables, including the interaction terms --------
step_center <- function(df, vars_center, means) {
  for (var in vars_center) {
    df$new <- df[, var] - means[var]
    names(df)[ncol(df)] <- paste0(var, "_C")
  }
  return(df)
}

# --------------------------
# a function to manually calculate log likehood when some of the coefficients are manually changed
# modified from https://stackoverflow.com/questions/39627615/calculating-loglik-by-hand-from-a-logistic-regression
get_logLik <- function(outcome, prediction, coefficients, family = binomial(logit)) {
    n <- length(outcome)
    wt <- rep(1, n) # or s_model$prior_weights if field exists
    deviance <- sum(family$dev.resids(outcome, prediction, wt))
    mod_rank <- sum(!is.na(coefficients)) # or s_model$rank if field exists

    aic <- family$aic(outcome, rep(1, n), prediction, wt, deviance) + 2 * mod_rank
    log_lik <- mod_rank - aic/2
    return(log_lik)
}

# --------------------------
# Root searching
root.search <- function(
  model,
  data,
  rcs_vars,
  k_vec,
  y_mean,
  rngList = NULL,
  length.out = 101,
  err = 1e-4
) {
  n_rcs <- length(rcs_vars)
  # check if rcs_vars and k_vec has the same length
  if (length(k_vec) != n_rcs) {
    stop("rcs_vars and k_vec must have the same length")
  }

  # create the grid
  df_grid <- data.frame(ind = 1:length.out)
  for (i in 1:n_rcs) {
    var <- rcs_vars[i]
    if (is.null(rngList)) {
      # without range
      vec <- data[[var]]
      df_grid$new <- seq(
        min(vec, na.rm = TRUE),
        max(vec, na.rm = TRUE),
        length.out = length.out
      )
    } else {
      # range determined
      rng <- rngList[[i]]
      df_grid$new <- seq(rng[1], rng[2], length.out = length.out)
    }
    names(df_grid)[which(names(df_grid) == "new")] <- var
  }
  df_grid <- df_grid %>% select(-ind)

  # check if there are other variables in the model other than those in rcs_vars
  vars <- names(coef(model))
  test <- sum(stringr::str_detect(
    vars,
    paste(c("(Intercept)", rcs_vars), collapse = "|")
  )) # add intercept
  if (test != length(vars)) {
    # other variables in the model
    vars_other <- vars[
      !stringr::str_detect(
        vars,
        paste(c("(Intercept)", rcs_vars), collapse = "|")
      )
    ]

    for (var in vars_other) {
      df_grid$new <- rep(0, length.out)
      names(df_grid)[which(names(df_grid) == "new")] <- var
    }
  }

  df_grid$y <- predict(model, newdata = df_grid)
  # check if roots are found, otherwise, reduce the ranges and search again.
  df_grid$diff <- df_grid$y - y_mean
  start <- df_grid %>% filter(diff >= 0)
  start <- start %>% filter(diff == min(diff))
  end <- df_grid %>% filter(diff < 0)
  end <- end %>% filter(diff == max(diff))

  if (abs(start$diff) <= err) {
    return(start)
  } else if (abs(end$diff) <= err) {
    return(end)
  } else {
    rngList <- list()
    for (i in 1:n_rcs) {
      rngList[[i]] <- sort(c(start[, rcs_vars[i]], end[, rcs_vars[i]]))
    }

    root.search(
      model,
      data,
      rcs_vars,
      k_vec,
      y_mean,
      rngList,
      length.out,
      err
    )
  }
}