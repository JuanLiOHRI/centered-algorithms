require(survival)

#' Modify `survival::brier` to work with prediction and outcomes instead
#' of fitted models
#' This only works for Cox model, not fine-gray model
#'
#' @param pred.prob  Vector of the predicted event probabilities
#' @param outcome.time Vector of the time variable
#' @param outcome.event Vector of the event variable
#' For other parameters, see `survival::brier`
#' @return     Dataframe of metrics
#' @export
brier.cox <- function(
  pred.prob,
  outcome.time,
  outcome.event,
  times,
  ties = TRUE,
  detail = FALSE,
  timefix = TRUE
) {
  Call <- match.call()
  # if (!inherits(fit, "coxph")) {
  #   stop("fit must be a coxph object")
  # }
  # if (missing(newdata)) {
  #   mf <- stats::model.frame(fit)
  # } else {
  #   mf <- stats::model.frame(fit, data = newdata)
  # }
  Y <- Surv(outcome.time, outcome.event) # Y <- mf[[1]]
  # if (!is.Surv(Y)) {
  #   stop("response must be a Surv object")
  # }
  # type <- attr(Y, "type")
  # if (!(type %in% c("right", "mright", "counting", "mcounting"))) {
  #   stop("response must be right censored")
  # }
  n <- nrow(Y)
  # ny <- ncol(Y)
  if (!is.logical(timefix) || length(timefix) > 1) {
    stop("invalid value for timefix option")
  }
  if (timefix) {
    Y <- aeqSurv(Y)
  }
  casewt <- rep(1, n)
  # casewt <- model.weights(mf)
  # if (is.null(casewt)) {
  #   casewt <- rep(1, n)
  # } else {
  #   if (!is.numeric(casewt)) {
  #     stop("weights must be numeric")
  #   }
  #   if (any(!is.finite(casewt))) {
  #     stop("weights must be finite")
  #   }
  #   if (any(casewt < 0)) {
  #     stop("weights must be non-negative")
  #   }
  #   casewt <- as.numeric(casewt)
  # }
  # id <- model.extract(mf, "id")
  # if (ny == 3 && (is.null(id))) {
  #   stop("id is required for start-stop data")
  # }
  # if (!is.null(id)) {
  #   if (is.null(attr(Y, "states"))) {
  #     ytemp <- Y
  #     attr(ytemp, "states") <- "event"
  #     check <- survcheck2(ytemp, id)
  #   } else {
  #     check <- survcheck2(Y, id)
  #   }
  #   if (any(check$flag > 0)) {
  #     stop("one or more flags are >0 in survcheck")
  #   }
  #   n.startstate <- sum(check$transitions[, 1] > 1)
  #   if (ny == 2) {
  #     samestart <- TRUE
  #   } else {
  #     etemp <- tapply(Y[, 1], id, min)
  #     samestart <- all(etemp == etemp[1])
  #   }
  # } else {
  #   check <- NULL
  # }
  # if (length(id) == 0) {
  #   id <- seq.int(n)
  # }
  # if (is.null(check) || (n.startstate == 1 & samestart)) {
  #   simple <- TRUE
  # } else {
  #   simple <- FALSE
  # }
  # if (!simple) {
  #   stop("delayed entry is not yet implemented")
  # }
  # if (efron && fit$method == "efron") {
  #   s0 <- survfit(Y ~ 1, weights = casewt, se.fit = FALSE, ctype = 2, stype = 2)
  # } else {
  #   s0 <- survfit(Y ~ 1, weights = casewt, se.fit = FALSE, stype = 1)
  # }
  # if (missing(times)) {
  #   times <- s0$time[s0$n.event > 0]
  # }
  s0 <- s0 <- survfit(Y ~ 1, weights = casewt, se.fit = FALSE, stype = 1)
  p0 <- 1 - summary(s0, times, extend = TRUE)$surv
  # if (missing(newdata)) {
  #   s1 <- survfit(fit, newdata = fit$call$data, se.fit = FALSE)
  # } else {
  #   s1 <- survfit(fit, newdata = newdata, se.fit = FALSE)
  # }
  p1 <- t(pred.prob) # p1 <- 1 - summary(s1, times, extend = TRUE)$surv
  dtime <- outcome.time # dtime <- Y[, ny - 1]
  dstat <- outcome.event # dstat <- Y[, ny]
  # The current version only works with one time point
  ntime <- length(times)
  if (ties) {
    mindiff <- min(diff(sort(unique(dtime))))
    dtime <- dtime + ifelse(dstat == 0, mindiff / 2, 0)
  }
  c0 <- survfit0(survfit(Surv(dtime, 1 - dstat) ~ 1, weights = casewt))
  b0 <- b1 <- matrix(0, nrow = n, ncol = ntime)
  casewt <- casewt / sum(casewt)
  brier <- matrix(0, ntime, 2)
  eff.n <- double(ntime)
  for (i in 1:ntime) {
    indx <- findInterval(pmin(dtime, times[i]), c0$time)
    wt <- ifelse(dtime < times[i] & dstat == 0, 0, casewt / c0$surv[indx])
    eff.n[i] <- 1 / sum(wt^2)
    b0[, i] <- ifelse(
      dtime > times[i],
      p0[i]^2,
      (dstat -
        p0[i])^2
    )
    b1[, i] <- ifelse(
      dtime > times[i],
      p1[i, ]^2,
      (dstat -
        p1[i, ])^2
    )
    brier[i, ] <- c(sum(wt * b0[, i]), sum(wt * b1[, i])) / sum(wt)
  }
  # ret <- list(
  #   rsquared = 1 - brier[, 2] / brier[, 1],
  #   brier = brier[,
  #     2
  #   ],
  #   times = times
  # )
  ret <- data.frame(
    model = c("Null model", "cox"),
    times = rep(times, 2),
    IPA = c(0, 1 - brier[, 2] / brier[, 1]),
    Brier = c(brier[, 1], brier[, 2])
  )

  if (detail) {
    ret$p0 <- p0
    ret$phat <- p1
    ret$eff.n <- eff.n
  }
  ret
}
