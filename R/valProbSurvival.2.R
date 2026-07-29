require(dplyr)
require(survival)
require(timeROC)
require(rms) # `cph`, `survest`
require(ggplot2)
require(riskRegression)

#' Modify `CalibrationCurves::valProbSurvival` to work with recalibrated cox model.
#'
#' @param pred.lp Vector of the predicted linear predictors
#' @param pred.prob  Vector of the predicted event probabilities
#' @param pred.prob.brier  Vector of the predicted event probabilities estimeated at timeHorizon - 0.01 for brier score
#' @param outcome.time Vector of the time variable
#' @param outcome.event Vector of the event variable
#' @param print.plot Control if plot the calibration plot
#' For other parameters, see `CalibrationCurves::valProbSurvival`
#' @return     Calibration curve and performance metrics
#' @export
valProbSurvival.2 <- function(
  pred.lp,
  pred.prob,
  pred.prob.brier,
  outcome.time,
  outcome.event,
  print.plot = TRUE,
  alpha = 0.05,
  timeHorizon = NULL,
  nk = 3,
  plotCal = c("none", "base", "ggplot"),
  addCox = FALSE,
  addRCS = TRUE,
  CL.cox = c("fill", "line"),
  CL.rcs = c("fill", "line"),
  xlab = "Predicted probability",
  ylab = "Observed proportion",
  xlim = c(-0.02, 1),
  ylim = c(-0.15, 1),
  lty.ideal = 1,
  col.ideal = "red",
  lwd.ideal = 1,
  lty.cox = 1,
  col.cox = "grey",
  lwd.cox = 1,
  fill.cox = "lightgrey",
  lty.rcs = 1,
  col.rcs = "black",
  lwd.rcs = 1,
  fill.rcs = rgb(177, 177, 177, 177, maxColorValue = 255),
  riskdist = "predicted",
  d0lab = "0",
  d1lab = "1",
  size.d01 = 5,
  dist.label = 0.01,
  line.bins = -0.05,
  dist.label2 = 0.04,
  length.seg = 0.85,
  legendloc = c(0.5, 0.27)
) {
  callFn <- match.call()
  pred <- obs <- lower <- upper <- xend <- yend <- NULL
  # JL: set default value of timeHorizon = NULL, force it to be checked
  if (is.null(timeHorizon)) {
    stop("Specifit time horizon is required.")
  }
  # JL: using specific prediciton and outcome variables instead of the whole test set
  valdata <- data.frame(
    time = outcome.time,
    event = outcome.event,
    LP = pred.lp,
    pred = pred.prob
  )
  # make sure the variables names are correct
  names(valdata) <- c("time", "event", "LP", "pred")
  # if (!inherits(fit, "coxph"))
  #     stop("Only model fits of class coxph are allowed")
  plotCal <- match.arg(plotCal)
  CL.cox <- match.arg(CL.cox)
  CL.rcs <- match.arg(CL.rcs)
  stats <- list()
  #valdata$LP = predict(fit, newdata = valdata, type = "lp")
  argzConc <- alist(data = valdata, reverse = TRUE)
  argzConc$obj <- as.formula("Surv(time, event) ~ LP") #argzConc$obj = update(fit$formula, "~ - . + LP")
  HarrellC <- do.call("concordance", argzConc)
  argzConc$timewt <- "n/G2"
  UnoC <- do.call("concordance", argzConc)
  res_C <- matrix(
    c(
      with(HarrellC, {
        c(
          concordance,
          concordance - qnorm(1 - alpha / 2) * sqrt(var),
          concordance + qnorm(1 - alpha / 2) * sqrt(var)
        )
      }),
      with(UnoC, {
        c(
          concordance,
          concordance - qnorm(1 - alpha / 2) * sqrt(var),
          concordance + qnorm(1 - alpha / 2) * sqrt(var)
        )
      })
    ),
    nrow = 2,
    ncol = 3,
    byrow = T,
    dimnames = list(c("Harrell C", "Uno C"), c("Estimate", "2.5 %", "97.5 %"))
  )
  stats$Concordance <- res_C
  res <- timeROC(
    T = valdata$time,
    delta = valdata$event,
    marker = valdata$LP,
    cause = 1,
    weighting = "marginal",
    times = c(timeHorizon, max(valdata$time) - 0.01), # max(as.numeric(fit$y)) - 0.01,
    iid = TRUE
  )
  UnoTDAUC <- with(
    res,
    matrix(
      c(
        c(
          AUC[[1]],
          AUC[[1]] - qnorm(1 - alpha / 2) * inference$vect_sd_1[[1]],
          AUC[[1]] + qnorm(1 - alpha / 2) * inference$vect_sd_1[[1]]
        ),
        c(
          AUC[[2]],
          AUC[[2]] - qnorm(1 - alpha / 2) * inference$vect_sd_1[[2]],
          AUC[[2]] + qnorm(1 - alpha / 2) * inference$vect_sd_1[[2]]
        )
      ),
      nrow = 2,
      ncol = 3,
      byrow = T,
      dimnames = list(names(res$AUC), c("Uno AUC", "2.5 %", "97. 5 %"))
    )
  )
  stats$TimeDependentAUC <- UnoTDAUC
  adjFormula <- as.formula("Surv(time, event) ~ 1") # adjFormula <- update(fit$formula, ~ -. + 1)
  obj <- summary(survfit(adjFormula, data = valdata), times = timeHorizon)
  obs_t <- 1 - obj$surv
  # valdata$pred <- predictRisk(fit, newdata = valdata, times = timeHorizon)
  exp_t <- mean(valdata$pred)
  OE_t <- obs_t / exp_t
  OE_summary <- c(
    OE = OE_t,
    `2.5 %` = OE_t *
      exp(
        -qnorm(
          1 -
            alpha / 2
        ) *
          sqrt(1 / obj$n.event)
      ),
    `97.5 %` = OE_t *
      exp(
        +qnorm(
          1 -
            alpha / 2
        ) *
          sqrt(1 / obj$n.event)
      ),
      Obs = obs_t, # JL OvsP
    Pavg = exp_t # JL OvsP
  )
  stats$Calibration$InTheLarge <- OE_summary
  calFormula <- as.formula("Surv(time, event) ~ LP") # calFormula <- update(fit$formula, ~ -. + LP)
  calCox <- cph(calFormula, x = TRUE, y = TRUE, surv = TRUE, data = valdata)
  # Claude: Use paste0() instead of substitute() to avoid scoping issues with nk parameter
  calRCSFormula <- as.formula(paste0("Surv(time, event) ~ rcs(LP, ", nk, ")"))
  vcal <- cph(calRCSFormula, x = TRUE, y = TRUE, surv = TRUE, data = valdata)
  datCox <- cbind.data.frame(
    obs = 1 - survest(calCox, times = timeHorizon, newdata = valdata)$surv,
    lower = 1 - survest(calCox, times = timeHorizon, newdata = valdata)$upper,
    upper = 1 -
      survest(calCox, times = timeHorizon, newdata = valdata)$lower,
    pred = valdata$pred
  )
  datCox <- datCox[order(datCox$pred), ]
  dat_cal <- cbind.data.frame(
    obs = 1 - survest(vcal, times = timeHorizon, newdata = valdata)$surv,
    lower = 1 - survest(vcal, times = timeHorizon, newdata = valdata)$upper,
    upper = 1 - survest(vcal, times = timeHorizon, newdata = valdata)$lower,
    pred = valdata$pred
  )
  dat_cal <- dat_cal[order(dat_cal$pred), ]
  absdiff_cph <- abs(dat_cal$pred - dat_cal$obs)
  stats$Calibration$Statistics <- c(
    ICI = mean(absdiff_cph),
    setNames(quantile(absdiff_cph, c(0.5, 0.9)), c("E50", "E90")),
    Emax = max(absdiff_cph)
  )
  gval <- coxph(calFormula, data = valdata)
  stats$Calibration$Slope <- c(
    `calibration slope` = unname(gval$coef),
    `2.5 %` = gval$coef - qnorm(1 - alpha / 2) * sqrt(gval$var),
    `97.5 %` = gval$coef + qnorm(1 - alpha / 2) * sqrt(gval$var)
  )
  # stats$Calibration$BrierScore <- Score(
  #   list(cox = fit), 
  #   formula = adjFormula,
  #   data = valdata, 
  #   conf.int = TRUE, 
  #   times = timeHorizon - 0.01, 
  #   cens.model = "km", 
  #   metrics = "brier", 
  #   summary = "ipa")$Brier$score
  stats$Calibration$BrierScore <- Score(
    list(cox = pred.prob.brier), 
    formula = Hist(time, event) ~ 1,
    data = valdata, 
    conf.int = TRUE, 
    times = timeHorizon - 0.01, 
    cens.model = "km", 
    metrics = "brier", 
    summary = "ipa")$Brier$score

  calCurves <- list(CoxCalibration = datCox, RCS = dat_cal)
  if (is.character(riskdist)) {
    if (riskdist == "calibrated") {
      x <- datCox$obs
      x[datCox$pred == 0] <- 0
      x[datCox$pred == 1] <- 1
    } else {
      x <- datCox$pred
    }
    y <- calCox$y[, 2]
    bins <- seq(0, min(1, max(xlim)), length = 101)
    x <- x[x >= 0 & x <= 1]
    f0 <- table(cut(x[y == 0], bins))
    f1 <- table(cut(x[y == 1], bins))
    j0 <- f0 > 0
    j1 <- f1 > 0
    bins0 <- (bins[-101])[j0]
    bins1 <- (bins[-101])[j1]
    f0 <- f0[j0]
    f1 <- f1[j1]
    maxf <- max(f0, f1)
    f0 <- (0.1 * f0) / maxf
    f1 <- (0.1 * f1) / maxf
    if (plotCal == "base") {
      par(xaxs = "i", yaxs = "i", las = 1)
      plot(
        dat_cal$pred,
        dat_cal$obs,
        type = "l",
        col = "white",
        lty = 1,
        xlim = xlim,
        ylim = ylim,
        lwd = 2,
        xlab = xlab,
        ylab = ylab
      )
      abline(0, 1, lty = lty.ideal, col = col.ideal, lwd = lwd.ideal)
      legCol <- c(Ideal = col.ideal)
      lt <- lty.ideal
      lw.d <- lwd.ideal
      marks <- NA
      if (addCox) {
        legCol <- c(legCol, `Cox calibration` = col.cox)
        lt <- c(lt, lty.cox)
        lw.d <- c(lw.d, lwd.cox)
        marks <- c(marks, NA)
        if (CL.cox == "line") {
          legCol <- c(legCol, `CL Cox` = col.cox)
          lt <- c(lt, 2)
          lw.d <- c(lw.d, 1)
          marks <- c(marks, NA)
          lines(datCox$pred, datCox$lower, type = "l", lty = 2, lwd = 2)
          lines(datCox$pred, datCox$upper, type = "l", lty = 2, lwd = 2)
        } else {
          polygon(
            x = with(datCox, c(pred, rev(pred))),
            y = with(datCox, c(upper, rev(lower))),
            col = fill.cox,
            border = NA
          )
        }
        lines(datCox$pred, datCox$obs, type = "l", lwd = lwd.cox, lty = lty.cox)
      }
      if (addRCS) {
        legCol <- c(legCol, `Flexible calibration (rcs)` = col.rcs)
        lt <- c(lt, lty.rcs)
        lw.d <- c(lw.d, lwd.rcs)
        marks <- c(marks, NA)
        if (CL.rcs == "line") {
          legCol <- c(legCol, `CL flexible (rcs)` = col.rcs)
          lt <- c(lt, 2)
          lw.d <- c(lw.d, 1)
          marks <- c(marks, NA)
          lines(dat_cal$pred, dat_cal$lower, type = "l", lty = 2, lwd = 2)
          lines(dat_cal$pred, dat_cal$upper, type = "l", lty = 2, lwd = 2)
        } else {
          polygon(
            x = with(dat_cal, c(pred, rev(pred))),
            y = with(dat_cal, c(upper, rev(lower))),
            col = fill.cox,
            border = NA
          )
        }
        lines(
          dat_cal$pred,
          dat_cal$obs,
          type = "l",
          lwd = lwd.rcs,
          lty = lty.rcs
        )
      }
      ll <- legendloc
      if (!is.logical(ll)) {
        if (!is.list(ll)) {
          ll <- list(x = ll[1], y = ll[2])
        }
      }
      legend(
        ll,
        c("Ideal calibration", "Cox calibration", "95% confidence interval"),
        col = c(2, 1, 1),
        lty = c(2, 1, 2),
        lwd = c(2, 2, 2),
        bty = "n",
        cex = 0.85
      )
      segments(bins1, line.bins, bins1, length.seg * f1 + line.bins)
      segments(bins0, line.bins, bins0, length.seg * -f0 + line.bins)
      lines(
        c(
          min(bins0, bins1) - 0.01,
          max(bins0, bins1) +
            0.01
        ),
        c(line.bins, line.bins)
      )
      text(
        max(bins0, bins1) + dist.label,
        line.bins +
          dist.label2,
        d1lab,
        cex = size.d01 * 0.25
      )
      text(
        max(bins0, bins1) + dist.label,
        line.bins -
          dist.label2,
        d0lab,
        cex = size.d01 * 0.25
      )
    } else if (plotCal == "ggplot") {
      gg <- ggplot(data.frame()) +
        geom_line(
          data = data.frame(x = 0:1, y = 0:1),
          aes(x = x, y = y, colour = "Ideal"),
          linewidth = lwd.ideal,
          show.legend = TRUE
        ) +
        labs(x = xlab, y = ylab)
      legCol <- c(Ideal = col.ideal)
      lt <- lty.ideal
      lw.d <- lwd.ideal
      marks <- NA
      if (addCox) {
        gg <- gg +
          geom_line(
            data = datCox,
            aes(x = pred, y = obs, color = "Cox calibration"),
            linetype = lty.cox,
            linewidth = lwd.cox
          )
        legCol <- c(legCol, `Cox calibration` = col.cox)
        lt <- c(lt, lty.cox)
        lw.d <- c(lw.d, lwd.cox)
        marks <- c(marks, NA)
        gg <- if (CL.cox == "line") {
          legCol <- c(legCol, `CL Cox` = col.cox)
          lt <- c(lt, 2)
          lw.d <- c(lw.d, 1)
          marks <- c(marks, NA)
          gg +
            geom_line(
              data = datCox,
              aes(x = pred, y = lower, color = "CL Cox"),
              linetype = 2,
              linewidth = 1
            ) +
            geom_line(
              data = datCox,
              aes(x = pred, y = upper, color = "CL Cox"),
              linetype = 2,
              linewidth = 1
            )
        } else {
          gg +
            geom_ribbon(
              data = datCox,
              aes(x = pred, ymin = lower, ymax = upper),
              fill = fill.cox
            )
        }
      }
      if (addRCS) {
        gg <- gg +
          geom_line(
            data = dat_cal,
            aes(x = pred, y = obs, color = "Flexible calibration (rcs)"),
            linetype = lty.cox,
            linewidth = lwd.cox
          )
        legCol <- c(legCol, `Flexible calibration (rcs)` = col.rcs)
        lt <- c(lt, lty.rcs)
        lw.d <- c(lw.d, lwd.rcs)
        marks <- c(marks, NA)
        gg <- if (CL.rcs == "line") {
          legCol <- c(legCol, `CL flexible (rcs)` = col.rcs)
          lt <- c(lt, 2)
          lw.d <- c(lw.d, 1)
          marks <- c(marks, NA)
          gg +
            geom_line(
              data = dat_cal,
              aes(x = pred, y = lower, color = "CL flexible (rcs)"),
              linetype = 2,
              linewidth = 1
            ) +
            geom_line(
              data = dat_cal,
              aes(x = pred, y = upper, color = "CL flexible (rcs)"),
              linetype = 2,
              linewidth = 1
            )
        } else {
          gg +
            geom_ribbon(
              data = dat_cal,
              aes(x = pred, ymin = lower, ymax = upper),
              fill = fill.rcs
            )
        }
      }
      gg <- gg +
        geom_segment(
          data = data.frame(
            x = bins1,
            xend = bins1,
            y = rep(line.bins, length(bins1)),
            yend = c(length.seg * f1 + line.bins)
          ),
          aes(x = x, y = y, xend = xend, yend = yend)
        ) +
        geom_segment(
          data = data.frame(
            x = bins0,
            xend = bins0,
            y = rep(line.bins, length(bins0)),
            yend = c(length.seg * -f0 + line.bins)
          ),
          aes(x = x, y = y, xend = xend, yend = yend)
        ) +
        geom_line(
          data = data.frame(
            x = c(min(bins0, bins1) - 0.01, max(bins0, bins1) + 0.01),
            y = c(line.bins, line.bins)
          ),
          aes(x = x, y = y)
        ) +
        annotate(
          geom = "text",
          x = max(bins0, bins1) + dist.label,
          y = line.bins +
            dist.label2,
          label = d1lab,
          size = size.d01
        ) +
        annotate(
          geom = "text",
          x = max(bins0, bins1) +
            dist.label,
          y = line.bins - dist.label2,
          label = d0lab,
          size = size.d01
        )
      gg <- gg +
        scale_color_manual("", values = legCol, breaks = names(legCol)) +
        guides(
          colour = guide_legend(
            override.aes = list(
              linetype = lt,
              shape = marks,
              linewidth = lw.d * 0.5,
              size = 7
            )
          )
        ) +
        theme_bw() +
        theme(
          plot.background = element_blank(),
          panel.border = element_rect(
            colour = "black",
            fill = NA,
            linewidth = 1
          ),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.margin = margin(11, 11, 5.5, 5.5, "points"),
          legend.position = "bottom"
        )
      gg <- gg + coord_cartesian(xlim = xlim, ylim = ylim)
      if (print.plot) {
        print(gg)
      }
    }
  }
  Results <- structure(
    list(
      call = callFn,
      stats = stats,
      Calibration = list(
        Slope = c(
          `Point estimate` = unname(stats$Calibration$Slope[1]),
          `Lower confidence limit` = unname(stats$Calibration$Slope[2]),
          `Upper confidence limit` = unname(stats$Calibration$Slope[3])
        )
      ),
      alpha = alpha,
      CalibrationCurves = calCurves
    ),
    class = "SurvivalCalibrationCurve"
  )
  if (plotCal == "ggplot") {
    Results$ggPlot <- gg
  }
  return(Results)
}