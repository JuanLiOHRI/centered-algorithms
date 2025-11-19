require("dplyr")

# Create a wrapper function to compare calibration under different noise to data ratio
compare <- function(data, sd_noise, seed) {
  # ========== Set seed ==========
  set.seed(seed)

  # ========== Gaussian (Normal) Noise ==========
  noise <- rnorm(n = nrow(data), mean = 0, sd = sd_noise)

  # ========== Add noise ==========
  data$y <- data$y + noise

  # ========== Split data ==========
  data.dev <- data %>% filter(type == "development")
  data.ext <- data %>% filter(type == "external")

  # ========== Centering ==========
  # -------- DEVELOPMENT --------
  vars_mean <- names(data.dev)
  # get mean values
  means.dev <- get_mean(data.dev, vars_mean)
  # centering vars in the data.dev dataset
  vars_center <- names(means.dev)
  data.dev <- step_center(data.dev, vars_center, means.dev)

  names(data.dev)

  # -------- EXTERNAL: `data.ext.d` --------
  data.ext.d <- step_center(data.ext, vars_center, means.dev) # means from the DEVELOPMENT dataset
  names(data.ext.d)

  # -------- EXTERNAL: `data.ext.t` --------
  vars_mean <- names(data.ext)
  # get mean values
  means.ext <- get_mean(data.ext, vars_mean)
  # centering vars in the data.ext dataset
  vars_center <- names(means.ext)
  data.ext.t <- step_center(data.ext, vars_center, means.ext) # means from the EXTERNAL dataset
  names(data.ext.t)

  # ========== The original model ==========
  fit.o <- lm(
    y ~ x1 +
      x2 +
      x3 +
      rcs(x4, k_x4) +
      x5 +
      rcs(x6, k_x6) +
      x1 * x2 +
      x2 * x3 +
      x3 * x5 +
      x2 * rcs(x4, k_x4) +
      x3 * rcs(x4, k_x4) +
      rcs(x4, k_x4) * rcs(x6, k_x6),
    data = data.dev
  )

  # ========== The centered model on the development dataset ==========
  fit.c <- lm(
    y ~ x1_cat2_C +
      x2_cat2_C +
      x2_cat3_C +
      x3_C +
      x4_rcs_1_C +
      x4_rcs_2_C +
      x4_rcs_3_C +
      x5_C +
      x6_rcs_1_C +
      x6_rcs_2_C +
      # interaction terms of x1 * x2
      x1_cat2_by_x2_cat2_C +
      x1_cat2_by_x2_cat3_C +
      # interaction terms of x2 * x3
      x2_cat2_by_x3_C +
      x2_cat3_by_x3_C +
      # interaction terms of x3 * x5
      x3_by_x5_C +
      # interaction terms of x2 * rcs(x4, k_x4)
      x2_cat2_by_x4_rcs_1_C +
      x2_cat3_by_x4_rcs_1_C +
      x2_cat2_by_x4_rcs_2_C +
      x2_cat3_by_x4_rcs_2_C +
      x2_cat2_by_x4_rcs_3_C +
      x2_cat3_by_x4_rcs_3_C +
      # interaction terms of x3 * rcs(x4, k_x4)
      x3_by_x4_rcs_1_C +
      x3_by_x4_rcs_2_C +
      x3_by_x4_rcs_3_C +
      # interaction terms of x4 * rcs(x6, k_x6)
      x4_rcs_1_by_x6_rcs_1_C +
      x4_rcs_2_by_x6_rcs_1_C +
      x4_rcs_3_by_x6_rcs_1_C +
      x4_rcs_1_by_x6_rcs_2_C +
      x4_rcs_2_by_x6_rcs_2_C +
      x4_rcs_3_by_x6_rcs_2_C,
    data = data.dev
  )

  # ========== Model Transport ==========
  # Model 1: The original model
  pred.1 <- predict(fit.o, newdata = data.ext)

  # Model 2: The centered model: y = mean(y.dev) + sum\[beta_i\*(x_i - mean(x_i.dev))\]
  pred.2 <- predict(fit.c, newdata = data.ext.d)

  # Update model intercept to the new outconme mean
  fit.c.2 <- fit.c
  fit.c.2$coefficients[1] <- mean(data.ext$y)

  # Model 3: The centered model: y = mean(y.ext) + sum\[beta_i\*(x_i - mean(x_i.ext))\]
  pred.3 <- predict(fit.c.2, newdata = data.ext.t)

  # ========== Calibration ==========
  res1 <- calibration(pred.1, data.ext$y, labelPos.x = 2500, labelPos.y = 10000, marginPlt = TRUE)
  res2 <- calibration(pred.2, data.ext$y, labelPos.x = 2500, labelPos.y = 10000, marginPlt = TRUE)
  res3 <- calibration(pred.3, data.ext$y, labelPos.x = 2500, labelPos.y = 10000, marginPlt = TRUE)

  #  ========== Output ==========
  res <- bind_rows(res1, res2, res3) %>% select(-group)
  res$model <- c("Model 1", "Model 2", "Model 3")
  res$sd_noise <- rep(sd_noise,3)
  res$y_mean_dev <- rep(mean(data.dev$y),3)
  res$y_mean_ext <- rep(mean(data.ext$y),3)

  return(res)
}