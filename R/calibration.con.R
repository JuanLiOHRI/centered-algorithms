require(dplyr)
require(ggplot2)
require(ggpubr) # `ggarrange`
require(Hmisc) # `cut2`

#' Calibration curves and metrics for continuous outcome
#'
#' @param prediction  The vector of predicted values.
#' @param outcome The vector of true outcome.
#' @param group The vector of group labels if want to assess calibration within subgroups
#' @param group_name A string vector, levels of group
#' @param g.group Number of groups when the subgroup variable is continuous
#' @param g Number of bins if want to bin the data instead of showing all data points
#' @param violin Set to TRUE if want to using violin plot instead of boxplot when binning the data
#' @param marginPlt Set to TRUE if want to show distributions at the margin
#' @param alpha Alpha value of data points
#' @param nrow Number of rows in the facet
#' @param labelPos.x x position of the label
#' @param labelPos.y y position of the label
#' @return     A named vector of calibration stats
#' @export
calibration.con <- function(
  prediction,
  outcome,
  group = NULL,
  group_name = NULL,
  g.group = 4,
  g = 1,
  violin = FALSE,
  marginPlt = FALSE,
  alpha = 0.5,
  nrow = NULL,
  labelPos.x = NULL,
  labelPos.y = NULL
) {
  df <- data.frame(prediction = prediction, outcome = outcome)
  group_con <- FALSE # flag of continuous group variable
  if (!is.null(group)) {
    if (length(outcome) != length(group)) {
      stop("lengths of outcome and group do not match")
    }
    if (is.factor(group) | is.logical(group) | is.character(group)) {
      if (!is.null(group_name)) {
        if (length(unique(group)) == length(group_name)) {
          df$group_fct <- as.character(group)
          df$group_fct <- factor(df$group_fct, levels = group_name)
        } else {
          stop(
            "categorical labels and group do not match, please check group_name"
          )
        }
      } else {
        df$group_fct <- factor(group)
        group_name <- levels(df$group_fct)
      }
    } else {
      group_con <- TRUE
      df$group_fct <- cut2(group, g = g.group)
      group_name <- levels(df$group_fct)
    }
  } else {
    df <- df %>% mutate(group_fct = "Overall")
    df$group_fct <- factor(df$group_fct)
    group_name <- levels(df$group_fct)
  }
  df <- df %>% mutate(er = abs(outcome - prediction))

  min <- min(min(outcome, na.rm = TRUE), min(prediction, na.rm = TRUE))
  max <- max(max(outcome, na.rm = TRUE), max(prediction, na.rm = TRUE))

  if (is.null(labelPos.x)) {
    labelPos.x <- min * 1.1
  }
  if (is.null(labelPos.y)) {
    labelPos.y <- max * 0.9
  }

  # ------- stats -----------
  # whole data
  rmse <- sqrt(mean((df$outcome - df$prediction)^2))
  lr <- lm(outcome ~ prediction, data = df)
  eavg <- mean(df$er, na.rm = TRUE)
  e50 <- median(df$er, na.rm = TRUE)
  e90 <- unname(quantile(df$er, 0.9, na.rm = TRUE))
  emax <- max(df$er, na.rm = TRUE)

  stats <- data.frame(
    group = "Overall",
    Obs = mean(df$outcome, na.rm = T),
    Pavg = mean(df$prediction, na.rm = T),
    RMSE = rmse,
    Intercept = unname(lr$coefficients[1]),
    Slope = unname(lr$coefficients[2]),
    Eavg = eavg,
    E50 = e50,
    E90 = e90,
    Emax = emax
  )
  stats_vec <- c("RMSE", "Intercept", "Slope", "Eavg", "E50")

  # each group
  if (!is.null(group)) {
    for (i in 1:length(group_name)) {
      temp <- df %>% filter(group_fct == group_name[i])
      lr <- lm(outcome ~ prediction, data = temp)
      rmse <- sqrt(mean((temp$outcome - temp$prediction)^2))
      eavg <- mean(temp$er, na.rm = TRUE)
      emax <- max(temp$er, na.rm = TRUE)
      e90 <- unname(quantile(temp$er, 0.9, na.rm = TRUE))
      stats <- rbind(
        stats,
        data.frame(
          group = group_name[i],
          Obs = mean(temp$outcome, na.rm = T),
          Pavg = mean(temp$prediction, na.rm = T),
          RMSE = rmse,
          Intercept = unname(lr$coefficients[1]),
          Slope = unname(lr$coefficients[2]),
          Eavg = eavg,
          E50 = e50,
          E90 = e90,
          Emax = emax
        )
      )
    }
  }

  # ------- calibration plot -------

  if (g == 1) {
    # single group, scatter plot
    pList <- list()

    for (i in seq_len(length(levels(df$group_fct)))) {
      res <- unlist(
        stats %>% filter(group == group_name[i]) %>% select(all_of(stats_vec))
      )
      res <- as.character(paste(
        paste(names(res), round(as.numeric(res), 2), sep = ": "),
        collapse = "\n"
      ))
      p <- ggplot(
        df %>% filter(group_fct == group_name[i]) %>% droplevels(),
        aes(prediction, outcome)
      ) +
        geom_point(alpha = alpha) +
        geom_smooth(method = 'loess', formula = 'y ~ x', color = "red") +
        labs(
          title = group_name[i],
          x = "Predicted value",
          y = "Observed value",
          color = NULL
        ) +
        xlim(min, max) +
        ylim(min, max) +
        annotate(
          "segment",
          x = min,
          y = min,
          xend = max,
          yend = max,
          linetype = "dotted",
          color = "black"
        ) +
        annotate(
          "text",
          x = labelPos.x,
          y = labelPos.y,
          label = res,
          hjust = "left"
        ) +
        coord_fixed() +
        theme_bw()

      if (marginPlt) {
        p <- ggExtra::ggMarginal(p, type = "histogram")
      }

      pList[[i]] <- p
    }

    if (is.null(group)) {
      p <- pList[[1]]
    } else {
      if (!is.null(nrow)) {
        p <- ggarrange(plotlist = pList, nrow = nrow)
      } else {
        p <- ggarrange(plotlist = pList)
      }
    }

    print(p)
  } else {
    # bin the data to box plot or violin plot
    pList <- list()

    for (i in seq_len(length(levels(df$group_fct)))) {
      res <- unlist(
        stats %>% filter(group == group_name[i]) %>% select(all_of(stats_vec))
      )
      res <- as.character(paste(
        paste(names(res), round(as.numeric(res), 2), sep = ": "),
        collapse = "\n"
      ))

      dfi <- df %>%
        filter(group_fct == group_name[i]) %>%
        droplevels() %>%
        mutate(prediction_bin = factor(ntile(prediction, n = g)))

      data.bin <- plyr::ddply(dfi, "prediction_bin", function(DF) {
        data.frame(median = plyr::numcolwise(median)(DF))
      })
      levels(dfi$prediction_bin) <- data.bin$median.prediction

      p <- ggplot(dfi, aes(x = prediction, y = outcome, group = prediction_bin))
      if (violin) {
        p <- p +
          geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))
      } else {
        p <- p +
          geom_boxplot(orientation = "x")
      }

      p <- p +
        geom_smooth(
          method = "loess",
          formula = 'y ~ x',
          color = "red",
          group = 1
        ) +
        labs(
          title = group_name[i],
          x = "Predicted value",
          y = "Observed value",
          color = NULL
        ) +
        xlim(min, max) +
        ylim(min, max) +
        annotate(
          "segment",
          x = min,
          y = min,
          xend = max,
          yend = max,
          linetype = "dotted",
          color = "black"
        ) +
        annotate(
          "text",
          x = labelPos.x,
          y = labelPos.y,
          label = res,
          hjust = "left"
        ) +
        coord_fixed() +
        theme_bw()

      if (marginPlt) {
        p <- p + geom_point(alpha = 0)
        p <- ggExtra::ggMarginal(p, type = "histogram")
      }

      pList[[i]] <- p
    }

    if (is.null(group)) {
      p <- pList[[1]]
    } else {
      p <- ggarrange(plotlist = pList, nrow = nrow)
    }

    print(p)
  }

  return(stats)
}
